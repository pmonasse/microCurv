#include "levelLine.h"
#include <algorithm>
#include <cmath>
#include <cassert>

/// North, West, South, East: direction of entry/exit in a dual pixel
typedef short Dir;
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
    bool follow(std::vector<bool>& visit, float l, Point& p);
    unsigned char level[4];
    Point vertex[4];
private:
    const unsigned char* _im;
    const size_t _w;
    Dir _d;
    void complete();
    float numSaddle() const; ///< Numerator of saddle value
    float denomSaddle() const; ///< Denominator of saddle value
    void find_corner(const Point& p, float l) const;
    bool mark_visit(std::vector<bool>& visit) const;
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

// Complete init by filling points 1 and 2, assuming 0 and 3 are already set.
void DualPixel::complete() {
    Dir d=_d, base=0;
    Point p=vertex[base];
    for(int i=0, base=0; i<2; i++) {
        p += delta[d];
        if(++d==4) d=0;
        level[++base] = _im[(size_t)p.y*_w+(size_t)p.x];
        vertex[base] = p;
    }
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
void DualPixel::find_corner(const Point& p, float l) const {
    float N=numSaddle(), D=denomSaddle();
    if(_d != S) return; // For the moment, only handle south orientation
    if(D*D<1e-4) return; // Degenerate hyperbola
    float delta = (D*l-N)/(D*D);
    float xs = (level[0]-level[1])/D;
    float ys = (level[0]-level[3])/D;
    float s1 = sign((p.x-vertex[0].x)-xs);
    float s2 = sign(delta);
    delta = s1*sqrt(s2*delta);
    xs += delta;
    ys += s2*delta;
    if(0<xs && xs<1 &&
       0<ys && ys<1) {
        xs += vertex[0].x;
        ys += vertex[0].y;
        return;
    }
}

/// Mark the edge as "visited", return \c false if already visited.
bool DualPixel::mark_visit(std::vector<bool>& visit) const {
    bool cont=true;
    if(_d==N || _d==S) {
        int base = (_d==S)? 0: 3;
        size_t i = (size_t)vertex[base].y*_w+(size_t)vertex[base].x;
        if(visit[i])
            cont = false;
        visit[i] = true;
    }
    return cont;
}

/// Find exit point of level line in dual pixel entering at \a p.
/// Return false if it closes the level line.
bool DualPixel::follow(std::vector<bool>& visit, float l, Point& p) {
    complete();
    bool cont=mark_visit(visit);
    find_corner(p, l);
    assert(level[0]<l && l<level[3]);
    bool right= (l<level[1]);
    bool left = (level[2]<l);
    if(left && right) { // saddle point: test l<saddle_level without division
        float num=numSaddle(), denom=denomSaddle();
        right = (denom>0)? (l*denom<num): (l*denom>num);
        left = !right;
    }
    int base=1;
    if(left) {
        if(++_d==4) _d=0;
        base = 2;
    } else if(right) {
        if(--_d<0) _d=3;
        base = 0;
    }
    float coord = linear(level[base], l, level[base+1]);
    p = vertex[base] + coord*delta[_d+1]; // Safe: delta[4]==delta[0]
    if(! right) {
        vertex[0] = vertex[base];
        level[0] = level[base];
    }
    if(! left) {
        vertex[3] = vertex[++base];
        level[3] = level[base];
    }
    return cont;
}

/// Extract line at level l
static void extract(const unsigned char* data, size_t w,
                    std::vector<bool>& visit, float l,
                    Point p, std::vector<Point>& line) {
    const Point delta(.5f, .5f);
    DualPixel dual(p, l, data, w);
    do
        line.push_back(p+delta);
    while(dual.follow(visit, l, p));
    line.push_back(p+delta); // close loop
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
