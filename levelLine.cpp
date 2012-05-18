#include "levelLine.h"
#include <algorithm>
#include <cmath>
#include <cassert>

/// North, West, South, East: direction of entry/exit in a dual pixel
typedef signed char Dir;
static const Dir S=0, E=1, N=2, W=3;

/// Vectors associated to the 4 directions
static const Point delta[] = {Point(0,+1), //S
                              Point(+1,0), //E
                              Point(0,-1), //N
                              Point(-1,0), //W
                              Point(0,+1)};//S again, to avoid modulo

/// Vector addition
inline Point operator+(Point p1, Point p2) {
    return Point(p1.x+p2.x, p1.y+p2.y);
}

inline Point& operator+=(Point& p1, Point p2) {
    p1.x += p2.x;
    p1.y += p2.y;
    return p1;
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

/// A dual pixel, square whose vertices are 4 data points
struct DualPixel {
    DualPixel(Point& p, float l, const unsigned char* im, size_t w);
    void follow(Point& p, float l, std::vector<Point>& line);
    bool mark_visit(std::vector<bool>& visit) const;
    unsigned char level[4];
    Point vertex[4];
private:
    const unsigned char* _im;
    const size_t _w;
    Dir _d;
    void update(bool left, bool right);
    float numSaddle() const; ///< Numerator of saddle value
    float denomSaddle() const; ///< Denominator of saddle value
    bool find_corner(Point& p, float l) const;
};

/// Return x for y=v on line joining (0,v0) and (1,v1).
inline float linear(float v0, float v, float v1) {
    return (v-v0)/(v1-v0);
}

/// Constructor
DualPixel::DualPixel(Point& p, float l, const unsigned char* im, size_t w)
: _im(im), _w(w), _d(S) {
    Dir d=_d;
    for(int i=0; i<4; i++) {
        level[i] = im[(size_t)p.y*w+(size_t)p.x];
        vertex[i] = p;
        p += delta[d];
        if(++d==4) d=0;
    }
    p += linear(level[0],l,level[3])*delta[_d+1]; // Safe: delta[4]==delta[0]
}

/// The saddle value can be expressed as a fraction. This is its numerator.
inline float DualPixel::numSaddle() const {
    return (level[0]*(float)level[2]-level[1]*(float)level[3]);
}

/// The saddle value can be expressed as a fraction. This is its denominator.
inline float DualPixel::denomSaddle() const {
    return ((float)level[0]+level[2]-level[1]-level[3]);
}

inline float sign(float f) {
    return (f>0)? +1: -1;
}

/// Find extremal curvature point on hyperbola branch
bool DualPixel::find_corner(Point& p, float l) const {
    float D=denomSaddle();
    if(D*D<1e-4) return false; // Degenerate hyperbola
    float delta = (D*l-numSaddle())/(D*D);
    bool vert = (_d==S||_d==N); // Vertical direction?
    if(! vert) { D=-D; delta=-delta; }
#define MOD4(d) ((d)&0x3) // Modulo 4
    int base=MOD4(S-_d+4); // S->0, E->3, N->2, W->1
    Point s((level[base]-level[MOD4(base+1)])/D,
            (level[base]-level[MOD4(base+3)])/D);
#undef MOD4
    float s1 = vert?
        sign((p.x-vertex[base].x)-s.x):
        sign((p.y-vertex[base].y)-s.y);
    float s2 = sign(delta);
    delta = s1*sqrt(s2*delta);
    s1 = 1;
    if(! vert) std::swap(s1,s2);
    s.x += s1*delta;
    s.y += s2*delta;
    if(! (0<s.x && s.x<1 &&
          0<s.y && s.y<1))
        return false;
    p = vertex[base] + s;
    return true;
}

// Update fields vertex and level, knowing whether we turn.
void DualPixel::update(bool left, bool right) {
    int base= (right? 0: (left? 2:1));
    if(! right) {
        vertex[0] = vertex[base];
        level[0] = level[base];
    }
    if(! left) {
        vertex[3] = vertex[++base];
        level[3] = level[base];
    }
    Dir d=_d;
    Point p=vertex[0];
    for(int i=1; i<=2; i++) {
        p += delta[d];
        if(++d==4) d=0;
        level[i] = _im[(size_t)p.y*_w+(size_t)p.x];
        vertex[i] = p;
    }
}

/// Find exit point of level line in dual pixel entering at \a p.
/// Return false if it closes the level line.
void DualPixel::follow( Point& p, float l, std::vector<Point>& line) {
    assert(level[0]<l && l<level[3]);
    Point c(p);
    if( find_corner(c,l) )
        line.push_back(c);
    bool right= (l<level[1]);
    bool left = (level[2]<l);
    if(left && right) { // saddle point: test l<saddle_level without division
        float num=numSaddle(), denom=denomSaddle();
        right = (denom>0)? (l*denom<num): (l*denom>num);
        left = !right;
    }
    if(left  && ++_d>3) _d=0;
    if(right && --_d<0) _d=3;
    update(left,right);
    float coord = linear(level[0], l, level[3]);
    p = vertex[0] + coord*delta[_d+1]; // Safe: delta[4]==delta[0]
}

/// Mark the edge as "visited", return \c false if already visited.
bool DualPixel::mark_visit(std::vector<bool>& visit) const {
    bool cont=true;
    if(_d==S) {
        size_t i = (size_t)vertex[0].y*_w+(size_t)vertex[0].x;
        cont = !visit[i];
        visit[i] = true;
    }
    return cont;
}

/// Extract line at level l
static void extract(const unsigned char* data, size_t w,
                    std::vector<bool>& visit, float l,
                    Point p, std::vector<Point>& line) {
    DualPixel dual(p, l, data, w);
    for(bool cont=true; cont;) {
        line.push_back(p);
        cont = dual.mark_visit(visit);
        dual.follow(p,l,line);
    }
    const Point delta(.5f, .5f);
    for(std::vector<Point>::iterator it=line.begin(); it!=line.end(); ++it)
        *it += delta;
}

/// Level lines extraction algorithm
void extract(const unsigned char* data, size_t w, size_t h,
             float offset, float step,
             std::list<LevelLine>& ll) {
    std::vector<bool> visit(w*h, false);
    for(float l=offset; l<255.0f; l+=step) {
        bool found=false;
        std::vector<bool>::const_iterator it=visit.begin();
        for(size_t i=0; i+1<h; i++) {
            for(size_t j=0; j+1<w; j++, ++it)
                if(data[i*w+j]<l && l<data[i*w+j+1] && !*it) {
                    found = true;
                    LevelLine line;
                    line.level=l;
                    ll.push_back(line);
                    Point p((float)j,(float)i);
                    extract(data,w, visit, l, p, ll.back().line);
                }
            ++it;
        }
        if(found) // Reinit for next level
            std::fill(visit.begin(), visit.end(), false);
    }
}
