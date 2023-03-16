// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file gass.h
 * @brief Geometric Affine Scale Space
 * @author Lionel Moisan <Lionel.Moisan@parisdescartes.fr>
 *         Pascal Monasse <monasse@imagine.enpc.fr>
 * 
 * Copyright (c) 2002, 2012-2014, Lionel Moisan, Pascal Monasse
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

#ifndef GASS_H
#define GASS_H

#include <vector>

struct DPoint {
    double x,y;
    DPoint(double x0, double y0): x(x0), y(y0) {}
    bool operator==(const DPoint& p) const { return (x==p.x&&y==p.y); }
    bool operator!=(const DPoint& p) const { return !operator==(p); }
    DPoint operator+(const DPoint& p) const { return DPoint(x+p.x,y+p.y); }
    DPoint operator-(const DPoint& p) const { return DPoint(x-p.x,y-p.y); }
    DPoint& operator/=(double d) { x/=d; y/=d; return *this; }
};

inline DPoint operator*(double d, const DPoint& p)
{ return DPoint(d*p.x,d*p.y); }

void gass(std::vector<DPoint>& curve,
          double first=0, double last=1, double prec=5.0, double maxStep=1.0);

#endif
