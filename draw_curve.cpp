#include "draw_curve.h"

/// max(0,min(abs(v),m-1))
static int saturate(float v, int m) {
    if(v<0)
        return 0;
    if(v>=m)
        return m-1;
    return (int)v;
}

/// Draw line in image
void draw_line(const Point& p, const Point& q, unsigned char v,
               unsigned char* im, int w, int h) {
    int x0=saturate(p.x,w), x1=saturate(q.x,w);
    int y0=saturate(p.y,h), y1=saturate(q.y,h);
    int sx = (x0<x1)? +1: -1;
    int sy = (y0<y1)? +1: -1;
    int dx=x1-x0, dy=y1-y0;
    int adx=sx*dx, ady=sy*dy;
    int x=0,y=0;
    if(adx>=ady) {
        int z=-adx/2;
        while(x!=dx) {
            im[(y+y0)*w+(x+x0)] = v;
            x += sx;
            z += ady;
            if(z>0) {
                y += sy;
                z -= adx;
            }
        }
    } else {
        int z=-ady/2;
        while(y!=dy) {
            im[(y+y0)*w+(x+x0)] = v;
            y += sy;
            z += adx;
            if(z>0) {
                x += sx;
                z -= ady;
            }
        }
    }
}

/// Draw curve in image
void draw_curve(const std::vector<Point>& curve, unsigned char v,
                unsigned char* im, int w, int h) {
    if(curve.empty())
        return;
    std::vector<Point>::const_iterator it=curve.begin();
    Point o = *it++;
    while(it != curve.end()) {
        draw_line(o, *it, v, im,w,h);
        o = *it++;
    }
}
