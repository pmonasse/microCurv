/**
 * @file image.cpp
 * @brief Image expansion and crop
 * @author Pascal Monasse <monasse@imagine.enpc.fr>
 * 
 * Copyright (c) 2011-2014, Pascal Monasse
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

/// Extend line horizontally by mirror effect of \a m pixels in the direction
/// opposite of \a dir. \a in points to the the first pixel to be duplicated.
///
/// If the line is xxxxab, running with \a in pointing over value a, \a m=2 and
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

/// Symmetrize an image by mirror effect
///
/// Return a newly allocated array, remember to \c delete[] it.
unsigned char* mirror(unsigned char* inImage, size_t w, size_t h, int m)
{
    // Prepare output image
    size_t ncol=w+2*m;
    size_t nrow=h+2*m;
    unsigned char* outImage = new unsigned char[ncol*nrow];
    std::fill(outImage, outImage+(ncol*nrow), 0);

    // Copy original image in center
    for(size_t i=0; i<h; i++)
        std::copy(inImage+i*w, inImage+(i+1)*w, outImage+ncol*(i+m)+m);

    // Horizontal mirror
    for(int i=0; i<h; i++) {
        mirrorH(outImage + ncol*(i+m)+m,     w, m, +1); // Left
        mirrorH(outImage + ncol*(i+m)+m+w-1, w, m, -1); // Right
    }

    // Vertical mirror
    mirrorV(outImage + ncol*m,       ncol, h, m, +1); // Top
    mirrorV(outImage + ncol*(m+h-1), ncol, h, m, -1); // Bottom

    return outImage;
}

/// Symmetrize by \a m pixels and paint border with value 0
///
/// Return a newly allocated array, remember to \c delete[] it.
unsigned char* sym_enlarge(unsigned char* in, size_t w, size_t h, int m) {
    unsigned char* sym = mirror(in, w, h, m);
    w += 2*m;
    h += 2*m;
    // Paint image border
    for(size_t i=0; i<h; i++) // Left and right
        sym[i*w]=sym[i*w+w-1]=0;
    std::fill(sym, sym+w, 0); // Top
    std::fill(sym+(h-1)*w, sym+h*w, 0); // Bottom
    return sym;
}
