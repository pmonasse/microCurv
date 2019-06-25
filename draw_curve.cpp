/**
 * @file draw_curve.cpp
 * @brief Draw a curve in an image
 * @author Pascal Monasse <monasse@imagine.enpc.fr>
 * 
 * Copyright (c) 2011-2014, Pascal Monasse
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

#include "draw_curve.h"

/// max(0,min(abs(v),m-1))
static int clip(float v, int m) {
    if(v<0)
        return 0;
    if(v>=m)
        return m-1;
    return (int)v;
}

/// Draw line in image
void draw_line(const Point& p, const Point& q, unsigned char v,
               unsigned char* im, int w, int h) {
    int x0=clip(p.x,w), x1=clip(q.x,w);
    int y0=clip(p.y,h), y1=clip(q.y,h);
    if(x0==x1 && y0==y1) {
        im[y0*w+x0] = v;
        return;
    }
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
    Point delta(.5f, .5f);
    std::vector<Point>::const_iterator it=curve.begin();
    Point o = *it++;
    while(it != curve.end()) {
        draw_line(o+delta, *it+delta, v, im,w,h);
        o = *it++;
    }
}
