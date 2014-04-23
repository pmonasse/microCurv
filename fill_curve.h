/**
 * @file fill_curve.h
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

#ifndef FILL_CURVE_H
#define FILL_CURVE_H

#include "levelLine.h"

template <typename T>
void fill_curve(const std::vector<Point>& line, T value,
                T* data, size_t w, size_t h,
                std::vector< std::vector<float> >* inter=0);

// Templates must have their implementation nearby
#include "fill_curve.cpp"

#endif
