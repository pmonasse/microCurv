#include "levelLine.h"
#include <algorithm>
#include <cassert>

/// North, West, South, East: direction of entry/exit in a dual pixel
typedef short Dir;
static const Dir N=0, W=1, S=2, E=3;

/// Vectors associated to the 4 directions
static const Point delta[] = {Point(0,-1), //N
                              Point(-1,0), //W
                              Point(0,+1), //S
                              Point(+1,0)};//E

/// Vector addition
inline Point operator+(Point p1, Point p2) {
    return Point(p1.x+p2.x, p1.y+p2.y);
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
    DualPixel(Point p, Dir d, const float* im, size_t w);
    float numSaddle() const; ///< Numerator of saddle value
    float denomSaddle() const; ///< Denominator of saddle value
    float level[4];
    Point vertex[4];
};

/// Constructor
DualPixel::DualPixel(Point p, Dir d, const float* im, size_t w) {
    switch(d) {
    case N:
        p.x = static_cast<float>((int)p.x+1);
        break;
    case W:
        p.y = static_cast<float>((int)p.y);
        break;
    case S:
        p.x = static_cast<float>((int)p.x);
        break;
    case E:
        p.y = static_cast<float>((int)p.y+1);
        break;
    }
    for(int i=0; i<4; i++) {
        level[i] = im[(size_t)p.y*w+(size_t)p.x];
        vertex[i] = p;
        p = p+delta[d];
        if(++d==4) d=0;
    }
}

/// The saddle value can be expressed as a fraction. This is its numerator.
inline float DualPixel::numSaddle() const {
    return (level[0]*level[2]-level[1]*level[3]);
}

/// The saddle value can be expressed as a fraction. This is its denominator.
inline float DualPixel::denomSaddle() const {
    return (level[0]+level[2]-level[1]-level[3]);
}

/// Return x for y=v on line joining (0,v0) and (1,v1).
inline float linear(float v0, float v, float v1) {
    return (v-v0)/(v1-v0);
}

/// Mark the edge as "visited", return \c false if already visited.
static bool mark_visit(const DualPixel& dual, std::vector<bool>& visit,size_t w,
                       Dir d) {
    bool cont=true;
    if(d==N || d==S) {
        int base = (d==S)? 0: 3;
        size_t i = (size_t)dual.vertex[base].y*w+(size_t)dual.vertex[base].x;
        if(visit[i])
            cont = false;
        visit[i] = true;
    }
    return cont;
}

/// Find exit point of level line in dual pixel entering at \a p.
/// Return false if it closes the level line.
static bool follow(const float* data, size_t w, size_t h,
                   std::vector<bool>& visit, float l,
                   Point& p, Dir& d) {
    DualPixel dual(p, d, data, w);
    bool cont=mark_visit(dual, visit, w, d);
    assert(dual.level[0]<l && l<dual.level[3]);
    bool right= (l<dual.level[1]);
    bool left = (dual.level[2]<l);
    if(left && right) { // saddle point: test l<saddle_level without division
        float num=dual.numSaddle(), denom=dual.denomSaddle();
        right = ((denom>0 && l*denom<num) || (denom<0 && l*denom>num));
        left = !right;
    }
    int base=1;
    if(left) {
        if(++d==4) d=0;
        base = 2;
    }
    if(right) {
        if(--d<0) d=3;
        base = 0;
    }
    float coord = linear(dual.level[base], l, dual.level[base+1]);
    p = dual.vertex[base] + coord*delta[(d+1)%4];
    return cont;
}

/// Extract line at level l
static void extract(const float* data, size_t w, size_t h,
                    std::vector<bool>& visit, float l,
                    Point p, std::vector<Point>& line) {
    const Point delta(.5f, .5f);
    Dir d=S;
    do
        line.push_back(p+delta);
    while(follow(data, w, h, visit, l, p, d));
    line.push_back(p+delta); // close loop
}

/// Level lines extraction algorithm
void extract(const float* data, size_t w, size_t h, float offset, float step,
             std::list<LevelLine>& ll) {
    std::vector<bool> visit(w*h, false);
    for(float l=offset; l<255.0f; l+=step) {
        bool found=false;
        for(size_t i=0; i+1<h; i++)
            for(size_t j=0; j+1<w; j++)
                if(! visit[i*w+j] &&
                   data[i*w+j]<l && l<data[i*w+j+1]) {
                    found = true;
                    LevelLine line;
                    line.level=l;
                    ll.push_back(line);
                    Point p(j+linear(data[i*w+j],l,data[i*w+j+1]), (float)i);
                    extract(data,w,h, visit, l, p, ll.back().line);
                    
                }
        if(found) // Reinit for next level
            std::fill(visit.begin(), visit.end(), false);
    }
}
