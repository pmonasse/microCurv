// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file image.cpp
 * @brief Image expansion and crop
 * @author Pascal Monasse <monasse@imagine.enpc.fr>
 * 
 * Copyright (c) 2011-2016, Pascal Monasse
 * All rights reserved.
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Pulic License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "image.h"

/// Input \a R. Format: wxh+x0+y0, eg 100x100+0+0
std::istream& operator>>(std::istream& str, Rect& R) {
    char c;
    str >> R.w >> c;
    if(str.fail() || c!='x') return str;
    str >> R.h >> c;
    if(str.fail() || c!='+') return str;
    str >> R.x >> c;
    if(str.fail() || c!='+') return str;
    str >> R.y;
    return str;
}

/// Intersection of rectangles \a R1 and \a R2.
Rect operator*(Rect R1, Rect R2) {
    int x = std::min(R1.x+R1.w, R2.x+R2.w);
    int y = std::min(R1.y+R1.h, R2.y+R2.h);
    if(R1.x < R2.x)
        R1.x = R2.x;
    if(R1.y < R2.y)
        R1.y = R2.y;
    R1.w = x-R1.x;
    R1.h = y-R1.y;
    return R1;
}

/// Extend line horizontally by mirror effect of \a m pixels in the direction
/// opposite of \a dir. \a in points to the the first pixel to be duplicated.
/// The first pixel written is the neighbor of \a in at position \a in-dir.
///
/// If the line is xxxxab, running with \a in pointing over value a, \a m=4 and
/// \a dir=+1 will result in abbaab.
static void mirrorH(unsigned char* in, size_t w, int m, int dir) {
    const int dirOut = -dir;
    unsigned char* out = in+dirOut;
    size_t n = w;
    for(int i=m-1; i>=0; i--) {
        *out = *in;
        out += dirOut;
        if(--n == 0) { // Change read direction when w values are copied
            n = w;
            dir = -dir;
        } else
            in += dir;
    }
}

/// Extend image vertically by mirror effect of \a m pixels in the direction
/// opposite of \a dir. \a in points to the first pixel to be duplicated.
/// The first line written is the one next to \a in in direction \a -dir.
///
/// If the image is
///
/// xxxx
/// xxxx
/// xxxx
/// abcd
/// efgh
///
/// \a in points over value a, \a w=4, \a h=2, \m=3, \dir=+1 will result in
///
/// efgh
/// efgh
/// abcd
/// abcd
/// efgh
static void mirrorV(unsigned char* in, size_t w, size_t h, int m, int dir) {
    dir *= w; // Operate line by line
    const int dirOut = -dir;
    unsigned char* out = in+dirOut;
    size_t n = h;
    for(int i=m-1; i>=0; i--) {
        std::copy(in, in+w, out);
        out += dirOut;
        if(--n == 0) { // Change read direction when h lines are copied
            n = h;
            dir = -dir;
        } else
            in += dir;
    }
}

/// Extract a (possibly too large) rectangle \a R from an image.
/// Missing pixels are filled in by mirror and periodicity.
/// Return a newly allocated array, remember to \c delete[] it.
unsigned char* extract(unsigned char* inImage, size_t w, size_t h, Rect R) {
    unsigned char* outImage = new unsigned char[R.w*R.h];
    std::fill(outImage, outImage+(R.w*R.h), 0);
    Rect R0 = {0,0,(int)w,(int)h};
    Rect R1 = R*R0;

    // Copy original image in center
    for(int i=R1.y; i<R1.y+R1.h; i++)
        std::copy(inImage+i*w+R1.x, inImage+i*w+R1.x+R1.w,
                  outImage+(i-R.y)*R.w+R1.x-R.x);

    // Horizontal mirror
    for(int i=R1.y-R.y; i<R1.y-R.y+R1.h; i++) {
        mirrorH(outImage + R.w*i+R1.x-R.x,        w, R1.x-R.x,          +1);
        mirrorH(outImage + R.w*i+R1.x-R.x+R1.w-1, w, R.x+R.w-R1.x-R1.w, -1);
    }

    // Vertical mirror
    mirrorV(outImage + R.w*(R1.y-R.y),        R.w, h, R1.y-R.y,          +1);
    mirrorV(outImage + R.w*(R1.y-R.y+R1.h-1), R.w, h, R.y+R.h-R1.y-R1.h, -1);

    return outImage;
}

/// Fill border of image \a im with constant value \a v.
void fill_border(unsigned char* im, size_t w, size_t h, unsigned char v) {
    for(size_t i=0; i<h; i++) // Left and right
        im[i*w]=im[i*w+w-1]=v;
    std::fill(im, im+w, v); // Top
    std::fill(im+(h-1)*w, im+h*w, v); // Bottom
}
