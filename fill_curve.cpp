// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file fill_curve.cpp
 * @brief Fill the interior of a closed curve in an image
 * 
 * (C) 2011-2014, 2019, Pascal Monasse <pascal.monasse@enpc.fr>
 */

#ifdef FILL_CURVE_H

#include <algorithm>
#include <cassert>

/// Sign of f2-f1
inline signed char sign(float f1, float f2) {
    return ((f1<f2)? +1: -1);
}

/// Is f an integer?
inline bool is_integer(float f) {
    return (f==(float)(int)f);
} 

/// Used to iterate over polyline vertex by vertex
struct PolyIterator {
    Point p; ///< Current vertex
    bool bHorizontal; ///< Along horizontal edgel?
    signed char dir; ///< right(+1)/left(-1) if horizontal, else down(+1)/up(-1)
    PolyIterator(const std::vector<Point>& curve);
    void add_point(const Point& pi, std::vector< std::vector<float> >& inter);
};

/// Find index of last point of the curve
static int last_point(const std::vector<Point>& curve) {
    Point p0 = curve.front();
    std::vector<Point>::const_reverse_iterator it=curve.rbegin(), end=curve.rend();
    size_t i = curve.size()-1;
    for(; it!=end; ++it, --i)
        if(*it!=p0)
            return i;
    return 0; // Single vertex
}

/// Constructor
PolyIterator::PolyIterator(const std::vector<Point>& curve) {
    size_t i = last_point(curve);
    Point q = curve[i]; // Previous vertex
    p = curve[0];
    if(q.y==p.y) {
        bHorizontal = is_integer(p.y);
        dir = sign(q.x,p.x);
    } else {
        bHorizontal = false;
        dir = sign(q.y,p.y);
    }
}

/// Allocate inter
static void init_inter(std::vector< std::vector<float> >& inter, size_t size)
{
    inter.resize(size);
    std::vector< std::vector<float> >::iterator it=inter.begin();
    for(; it!=inter.end(); ++it)
        it->clear();
}

/// Add bound of interval on line iy at position x
static void bound(std::vector< std::vector<float> >& inter, float x, int iy) {
    if(0<=iy && iy<(int)inter.size())
        inter[iy].push_back(x);
}

/// Add segment to point i to current polyline: see [2]Figure 4 for the rules.
void PolyIterator::add_point(const Point& pi,
                             std::vector< std::vector<float> >& inter) {
    Point q = p;
    p = pi;
    signed char dirP = dir; // Previous direction

    if(q.y==p.y) { // Horizontal segment
        if(q.x!=p.x && is_integer(q.y)) {
            dir = sign(q.x,p.x);
            if(bHorizontal) { // Half-turn, rule (f)
                if(dirP!=dir)
                    bound(inter, q.x, (int)q.y); 
            } else { // Rules (b), (c)
                bHorizontal = true; // First among horizontal edgels
                if(dirP==dir) // Rule (b)
                    bound(inter, q.x, (int)q.y);
            }
        }
        return;
    }

    dir = sign(q.y,p.y);
    int iy1 = (int)q.y;
    int iy2 = (int)p.y + dir;
    float a = (q.x-p.x) / (q.y-p.y); // Slope

    if(bHorizontal) { // Away from horizontal edgel, rules (d), (e)
        bHorizontal = false;
        if(dirP!=dir) // Rule (d)
            bound(inter, q.x, iy1);
        iy1 += dir;
    } else if(dir!=dirP && q.y==(float)iy1) { // Local peak, rule (g)
        bound(inter, q.x, iy1); // Single point interval
        bound(inter, q.x, iy1);
        iy1 += dir;
    } else if(dir>0 && (float)iy1<q.y)
        iy1 += dir;

    for(int j=iy1; j!=iy2; j+=dir) { // Rule (a)
        if(dir > 0) {
            if(p.y<=(float)j) continue; // Out of bounds
        } else
            if((float)j<=p.y) continue; // Out of bounds
        float xj = q.x + a*((float)j-q.y);
        assert((q.x<=xj && xj<=p.x) || (p.x<=xj && xj<=q.x));
        bound(inter, xj, j);
    }
}

/// Fill curve with a single vertex
template <class T>
void fill_point(Point p, T value, T* out, size_t w) {
    if(is_integer(p.x) && is_integer(p.y))
        out[(int)p.y*w+(int)p.x] = value;
}

/// Fill line of image
template <typename T>
void fill_line(T value, T* im, T* end, std::vector<float>& inter) {
    std::sort(inter.begin(), inter.end());
    bool bIn = false; // inside/outside
    std::vector<float>::const_iterator it=inter.begin();
    for(;it!=inter.end() && *it<0; ++it) // Handle curve outside left boundary
        bIn = !bIn;
    if(it==inter.end()) return;
    if(bIn)
        std::fill(im, std::min(end,im+(int)*it), value);
    float i = (float)(int)*it; // Current pixel number (int stored as float)
    im += (int)i;
    for(; im<end; i++, im++) {
        while(*it<i) {
            bIn = !bIn;
            if(++it == inter.end()) {
                assert(!bIn);
                return;
            }
        }
        if(bIn || *it==i)
            *im = value;
    }
}

/// Fill in intervals defined by inter
template <typename T>
void fill_inter(T value, T* im, size_t w, size_t h,
                std::vector< std::vector<float> >& inter) {
    for(size_t i=0; i<h; i++)
        if(! inter[i].empty())
            fill_line(value, im+i*w, im+(i+1)*w, inter[i]);
}

/// Fill interior region of curve.
template <typename T>
void fill_curve(const std::vector<Point>& line, T value,
                T* out, size_t w, size_t h,
                std::vector< std::vector<float> >* inter) {
    if(line.empty())
        return;
    PolyIterator p(line);
    if(p.dir==0) { // Single vertex
        fill_point(line.front(), value, out,w);
        return;
    }

    bool bAllocInter=(inter==0);
    if(bAllocInter)
        inter = new std::vector< std::vector<float> >(h);
    else
        init_inter(*inter, h);

    std::vector<Point>::const_iterator it=line.begin()+1;
    for(; it!=line.end(); ++it)
        p.add_point(*it, *inter);
    p.add_point(line.front(), *inter); // Close polygon

    fill_inter(value, out, w, h, *inter);
    if(bAllocInter)
        delete inter;
}

#endif
