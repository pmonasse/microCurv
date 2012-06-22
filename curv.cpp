#include "curv.h"
#include "levelLine.h"
#include <algorithm>
#include <cmath>

#define SIGN(x) (((x)>0.0)?1: -1)

/// Distance between points
static float dist(Point p, const Point& q) {
    p.x -= q.x; p.y -= q.y;
    return sqrt(p.x*p.x+p.y*p.y);
}

inline float det(const Point& a, const Point& b, const Point& c) {
    return -((b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x));
}

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

static void curv(const std::vector<Point>& curve, signed char sign, int w,
                 std::vector< std::vector<float> >& inter)
{
    int n = last_point(curve);
    const Point* q = &curve[n]; 
    const Point* p = &curve[0];
    const Point* r = &curve[1];
    float u = dist(*p,*q);

    for(int k=0; k<n; k++,r++){
        if(k+1==n)
            r = &curve[0];
        float v = dist(*p,*r);
        float d = u * v * dist(*q,*r); 
        float K = ((d==0)? 0: 2*sign*det(*p,*q,*r)/d);
        int i = (int)(curve[k].x-0.5f);
        int j = (int)(curve[k].y-0.5f);
        inter[j*w+i].push_back(K);
        q = p;
        p = r;
        u = v;
    }
}

float median(std::vector<float>& list) {
    size_t size = list.size();
    std::vector<float>::iterator m=list.begin()+(size-1)/2;
    std::nth_element(list.begin(), m, list.end());
    float v = *m;
    if(size%2 == 0)
        v = (v+*std::min_element(m+1, list.end()))*0.5f;
    return v;
}

void curv(const std::vector<LevelLine*>& ll, const std::vector<bool>& positive,
          float* out, int w, int h) {
    std::vector< std::vector<float> > inter(w*h);
 
    // Compute curvature
    std::vector<LevelLine*>::const_iterator it=ll.begin();
    for(int i=0; it!=ll.end(); ++it, ++i)
        if((*it)->line.size()>10)
            curv((*it)->line, positive[i]? +1: -1, w, inter);
       
    // Median curvature
    for(int i=0; i<w*h; i++, out++)
        if(! inter[i].empty()) {
            float m =  median(inter[i]);
            *out = - SIGN(m)*std::min(std::abs(m),1.0f);
        }
 }
