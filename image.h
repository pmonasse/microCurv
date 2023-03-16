// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file image.h
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

#ifndef IMAGE_H
#define IMAGE_H

#include <algorithm>
#include <istream>
#include <cstddef>

/// A rectangle.
struct Rect {
    int x, y, w, h;
};
std::istream& operator>>(std::istream& str, Rect& R);
inline std::ostream& operator<<(std::ostream& str, const Rect&) { return str; }

unsigned char* extract(unsigned char* inImage, size_t w, size_t h, Rect R);
void fill_border(unsigned char* im, size_t w, size_t h, unsigned char v=0);

/// Crop rectangle \a r from image \a in having \a ch channels.
///
/// Return a newly allocated array, remember to \c delete[] it.
template <typename T>
static T* crop(const T* in, size_t w, size_t h, Rect r, int ch=1) {
    const int npix = r.w*r.h;
    T* out = new T[ch*npix];
    for(int j=0; j<ch; j++, in+=w*h)
        for(int i=0; i<r.h; i++) {
            const T* p = in+w*(r.y+i)+r.x;
            std::copy(p, p+r.w, out+j*npix+i*r.w);
        }
    return out;
}

#endif
