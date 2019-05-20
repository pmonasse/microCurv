/**
 * @file levelLine.cpp
 * @brief Extraction of level lines from an image
 * @author Pascal Monasse <monasse@imagine.enpc.fr>
 * 
 * Copyright (c) 2011-2016, 2019 Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "levelLine.h"
#include <cmath>
#include <limits>
#include <cassert>

/// South, East, North, West: directions of entry/exit in a dual pixel.
/// Left turn=+1. Right turn=-1. Opposite=+2.
typedef signed char Dir;
static const Dir S=0, /*E=1,*/ N=2 /*, W=3*/;

/// Vectors associated to the 4 directions
static const Point delta[] = {Point(0,+1), //S
                              Point(+1,0), //E
                              Point(0,-1), //N
                              Point(-1,0), //W
                              Point(0,+1)};//S again, to avoid modulo

/// Vector subtraction
inline Point operator-(Point p1, Point p2) {
    return Point(p1.x-p2.x, p1.y-p2.y);
}

/// Vector multiplication
inline Point operator*(float f, Point p) {
    return Point(f*p.x, f*p.y);
}

/// Print list of points of the line (but not the level): x0 y0 x1 y1...
std::ostream& operator<<(std::ostream& str, const LevelLine& l) {
    std::vector<Point>::const_iterator it;
    for(it=l.line.begin(); it!=l.line.end(); ++it)
        str << it->x << " " << it->y << " ";
    return str;
}

/// Apply zoom factor to points of \a line.
void zoom_line(std::vector<Point>& line, float zoom) {
    std::vector<Point>::iterator it=line.begin(), end=line.end();
    for(; it!=end; ++it) {
        it->x *= zoom;
        it->y *= zoom;
    }
}

/// A mobile dual pixel, square whose vertices are 4 data points.
/// This is the main structure to extract a level line, moving from dual pixel
/// to an adjacent one until coming back at starting point. The entry direction
/// of the level line is recorded: south means it enters from the top horizontal
/// edgel, east from the right vertical, north from the bottom horizontal and
/// west from the right.
/// The object stores the levels at its 4 vertices (data points of the image),
/// in clockwise order starting from the top left vertex.
class DualPixel {
public:
    DualPixel(Point& p, float l, const unsigned char* im, size_t w);
    void follow(Point& p, float l, int ptsPixel, std::vector<Point>& line);
    bool mark_visit(std::vector<bool>& visit,
                    std::vector< std::vector<Inter> >* inter, size_t idx,
                    const Point& p) const;
private:
    const unsigned char* _im; ///< The image stored as 1D array.
    const size_t _w; ///< Number of columns of image.
    unsigned char _level[4]; /// The levels at the 4 data points.
    Point _pos; /// The position of the top-left vertex of the dual pixel.
    Dir _d; /// Direction of entry into dual pixel.

    void move(bool left, bool right);
    float numSaddle() const; ///< Numerator of saddle value
    float denomSaddle() const; ///< Denominator of saddle value
    bool find_corner(Point& p, float l, Point& s, float& delta) const;
    void sample(const Point& p1, const Point& p2, int ptsPixel,
                const Point& s, float delta, std::vector<Point>& line);
};

/// Return x for y=v on line joining (0,v0) and (1,v1).
inline float linear(float v0, float v, float v1) {
    return (v-v0)/(v1-v0);
}

/// Constructor.
/// \param[in,out] p the edgel endpoint at the right of incoming direction.
/// \param l the level of the level line.
/// \param im the values of pixels in a 1D array.
/// \param w the number of pixel columns in \a im.
/// The incoming direction is always supposed to be south, so the level line is
/// crossing the edgel from \a p to \a p+(1,0). It means the starting point of
/// the level line is at \a p+(x,0), with 0<x<1. As output, \a p is moved to
/// this position.
DualPixel::DualPixel(Point& p, float l, const unsigned char* im, size_t w)
: _im(im), _w(w), _pos(p), _d(S) {
    Dir d=_d;
    for(int i=0; i<4; i++) {
        _level[i] = im[(size_t)p.y*w+(size_t)p.x];
        p += delta[d];
        if(++d==4) d=0; // Not useful if S=0, but better be safe
    }
    p += linear(_level[0],l,_level[3])*delta[_d+1]; // Safe: delta[4]==delta[0]
}

/// The saddle value can be expressed as a fraction. This is its numerator.
inline float DualPixel::numSaddle() const {
    return (_level[0]*(float)_level[2]-_level[1]*(float)_level[3]);
}

/// The saddle value can be expressed as a fraction. This is its denominator.
inline float DualPixel::denomSaddle() const {
    return ((float)_level[0]+_level[2]-_level[1]-_level[3]);
}

inline float sign(float f) {
    return (f>0)? +1: -1;
}

/// Decompose hyperbola branch.
/// Inside the dual pixel, the level set has implicit equation
/// \f[ D*(x-xs)(y-ys)+N/D = l. \f]
/// This is true only if \f$D\neq0\f$ (otherwise we have just a line segment).
/// The center of the hyperbola (xs,ys) is a saddle point, its level is N/D.
/// Of interest, the vertex of the hyperbola is a point of maximal curvature.
/// It is located at
/// \f[(xs,ys)+(\pm\sqrt{|(Dl-N)/D^2|},\pm\sqrt{|(Dl-N)/D^2|}).\f]
/// The signs are determined by the fact that the vertex is in the same quadrant
/// as the input point p with respect to (xs,ys).
/// The equation of hyperbola is written
/// \f[ (x-xs)(y-ys) = \delta. \f]
/// \param[in] p a point on an edgel of the dual pixel and on the hyperbola.
/// \param[out] p vertex of hyperbola (point of extremal curvature)
/// \param l level.
/// \param[out] s center of the hyperbola, saddle point of the image.
/// \param[out] delta parameter of hyperbola (sqrt(2*delta) is semi major axis)
/// \return Do we really have a hyperbola? It may be a segment.
bool DualPixel::find_corner(Point& p, float l, Point& s, float& delta) const {
    float D=denomSaddle();
    if(D*D<1e-4) return false; // Degenerate hyperbola
    delta = (D*l-numSaddle())/(D*D);
    s.x = _pos.x + (_level[0]-_level[1])/D;
    s.y = _pos.y + (_level[0]-_level[3])/D;
    float d = sqrt(std::abs(delta));
    p.x = s.x + sign(p.x-s.x)*d;
    p.y = s.y + sign(p.y-s.y)*d;
    return (_pos.x<p.x && p.x<_pos.x+1 && _pos.y<p.y && p.y<_pos.y+1);
}

/// Sample branch of hyperbola from p1 to p2 of equation (x-xs)(y-ys)=delta.
/// \param p1 start point.
/// \param p2 end point.
/// \param ptsPixel number of points of discretization per pixel.
/// \param s center of hyperbola.
/// \param delta parameter of hyperbola (sqrt(2) times semi major axis).
/// \param[out] line where the sampled points are stored.
void DualPixel::sample(const Point& p1, const Point& p2, int ptsPixel,
                       const Point& s, float delta, std::vector<Point>& line) {
    if(ptsPixel<2) return;
    Point p = p2-p1;
    if(p.x<0) p.x=-p.x;
    if(p.y<0) p.y=-p.y;
    if(p.x>p.y) { // Uniform sample along x
        int n = int(p.x*ptsPixel+0.99f);
        float dx = (p2.x-p1.x)/n;
        p = p1;
        for(int i=1; i<n; i++) {
            p.x += dx;
            p.y = s.y + delta/(p.x-s.x);
            line.push_back(p);
        }
    } else { // Uniform sample along y
        int n = int(p.y*ptsPixel+0.99f);
        float dy = (p2.y-p1.y)/n;
        p = p1;
        for(int i=1; i<n; i++) {
            p.y += dy;
            p.x = s.x + delta/(p.y-s.y);
            line.push_back(p);
        }
    }
}

/// Move to next adjacent dual pixel, knowing whether we turn.
/// \param left make a left turn?
/// \param right make a right turn?
/// Parameters \a left and \a right are exclusive. They can be both \c false,
/// meaning no turn occurs.
void DualPixel::move(bool left, bool right) {
    // update direction
    if(left  && ++_d>3) _d=0;
    if(right && --_d<0) _d=3;
    // update top-left vertex
    _pos += delta[_d];
    // update levels at vertices
    size_t ind = (size_t)_pos.y*_w+(size_t)_pos.x;
    _level[0] = _im[ind];    _level[3] = _im[ind+1];
    _level[1] = _im[ind+_w]; _level[2] = _im[ind+_w+1];
}

/// The dual pixel is moved to the adjacent one. Find exit point of level line
/// entering at \a p in the dual pixel. The level line is sampled up to there. 
/// \param[in,out] p entry point into the dual pixel
/// \param l level of the level line
/// \param ptsPixel number of points of discretization per pixel.
/// \param[out] line intermediate samples stored here.
void DualPixel::follow(Point& p, float l, int ptsPixel,
                       std::vector<Point>& line) {
    assert(_level[_d]<l && l<_level[(_d+3)%4]);
    Point b(p), m(p), s(p); // b=init, m=corner, s=saddle
    float xy = std::numeric_limits<float>::quiet_NaN();
    bool corner = (ptsPixel>0) && find_corner(m,l,s,xy);
    bool left  = (l>_level[(_d+2)%4]); // Is there an exit at the left?
    bool right = (l<_level[(_d+1)%4]); // Is there an exit at the right?
    if(left && right) { // saddle point
        float num=numSaddle(), denom=denomSaddle();
        right = (denom>0)? (l*denom<num): (l*denom>num);// l<num/denom? w/o div
        left = !right;
    }
    move(left,right);
    float coord = linear(_level[_d], l, _level[(_d+3)%4]);
    p = _pos;
    for(Dir d=0; d<_d; d++) p += delta[d];
    p += coord*delta[_d+1]; // Safe: delta[4]==delta[0]
    if(xy == xy) { // Hyperbola: do not sample otherwise (ie, straight)
        if(corner) { // Sample until corner point m
            sample(b,m, ptsPixel, s,xy, line);
            line.push_back(b=m);
        }
        sample(b,p, ptsPixel, s,xy, line); // Sample until end point
    }
}

/// Mark the edge as "visited", return \c false if already visited.
/// \param visit stores the edgels traversed from the south at current level.
/// \param inter (optional) rows of image traversed are marked with \a idx.
/// \param idx a unique identifier for the level line.
/// \param p the current position in the tracking of the level line.
/// \return whether the tracking must continue (loop not closed yet).
/// When we go through a horizontal data row and going south, we store the
/// visit. If the edgel was already visited at current level, we came back
/// at starting point and must stop.
bool DualPixel::mark_visit(std::vector<bool>& visit,
                           std::vector< std::vector<Inter> >* inter,
                           size_t idx, const Point& p) const {
    bool cont=true;
    if(_d==S) {
        size_t i = (size_t)_pos.y*_w+(size_t)_pos.x;
        cont = !visit[i];
        visit[i] = true;
    }
    if(inter && cont && (_d==S||_d==N))
        (*inter)[(size_t)p.y].push_back( Inter(p.x,idx) );
    return cont;
}

/// Extract level line passing through a given starting point. 
/// \param data the values of pixels in a 1D array.
/// \param w the number of pixel columns in \a data.
/// \param visit array to store the visited explored horizontal edgels.
/// \param ptsPixel number of points of discretization per pixel.
/// \param p the starting point.
/// \param[in,out] ll the level line: level already stored, find its line.
/// \param idx a unique identifier for the level line.
/// \param inter[out] (optional) rows of image traversed are marked with \a idx.
/// \a inter is used to recover the tree hierarchy at the end, could be
/// omitted if the tree is not required, in which case \a idx is unused.
static void extract(const unsigned char* data, size_t w,
                    std::vector<bool>& visit, int ptsPixel,
                    Point p, LevelLine& ll, size_t idx,
                    std::vector< std::vector<Inter> >* inter) {
    DualPixel dual(p, ll.level, data, w);
    while(true) {
        ll.line.push_back(p);
        if(! dual.mark_visit(visit,inter,idx,p))
            break;
        dual.follow(p,ll.level,ptsPixel,ll.line);
    }
    const Point delta(.5f, .5f);
    for(std::vector<Point>::iterator it=ll.line.begin(); it!=ll.line.end();++it)
        *it += delta;
}

/// Level lines extraction algorithm.
/// \param data the values of pixels in a 1D array.
/// \param w the number of pixel columns in \a data.
/// \param h the number of pixel lines in \a data.
/// \param offset the minimal level to be extracted.
/// \param step the level step for extraction.
/// \param ptsPixel number of points of discretization per pixel.
/// \param[out] ll storage for the extracted level lines.
/// \param inter[out] (optional) rows of image traversed by ll are marked.
void extract(const unsigned char* data, size_t w, size_t h,
             float offset, float step, int ptsPixel,
             std::vector<LevelLine*>& ll,
             std::vector< std::vector<Inter> >* inter) {
    std::vector<bool> visit(w*h, false);
    if(inter) {
        assert(inter->empty());
        inter->resize(h);
    }
    for(float l=offset; l<255.0f; l+=step) {
        bool found=false;
        std::vector<bool>::const_iterator it=visit.begin();
        for(size_t i=0; i+1<h; i++) {
            for(size_t j=0; j+1<w; j++, ++it)
                if(data[i*w+j]<l && l<data[i*w+j+1] && !*it) {
                    found = true;
                    LevelLine* line = new LevelLine(l);
                    Point p((float)j,(float)i);
                    extract(data,w, visit, ptsPixel, p, *line, ll.size(),inter);
                    ll.push_back(line);
                }
            ++it;
        }
        if(found) // Reinit for next level
            std::fill(visit.begin(), visit.end(), false);
    }
}
