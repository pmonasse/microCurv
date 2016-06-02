/**
 * @file levelLine.h
 * @brief Extraction of level lines from an image
 * @author Pascal Monasse <monasse@imagine.enpc.fr>
 * 
 * Copyright (c) 2011-2016, Pascal Monasse
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

#ifndef LEVELLINE_H
#define LEVELLINE_H

#include <vector>
#include <iostream>

struct Point {
    float x, y;
    Point() {}
    Point(float x0, float y0): x(x0), y(y0) {}
    bool operator==(const Point& p) const;
    bool operator!=(const Point& p) const;
};

inline bool Point::operator==(const Point& p) const {
    return (x==p.x && y==p.y); }
inline bool Point::operator!=(const Point& p) const {
    return !operator==(p); }

/// Vector addition
inline Point operator+(Point p1, Point p2) {
    return Point(p1.x+p2.x, p1.y+p2.y); }
inline Point& operator+=(Point& p1, Point p2) {
    p1.x += p2.x;
    p1.y += p2.y;
    return p1;
}

void zoom_line(std::vector<Point>& line, float zoom);

/// Level line: a level and a polygonal line
struct LevelLine {
    float level;
    std::vector<Point> line;
    LevelLine(float l): level(l) {}
    void fill(unsigned char* data, size_t w, size_t h,
              std::vector< std::vector<float> >* inter=0) const;
};

std::ostream& operator<<(std::ostream& str, const LevelLine& line);

/// Abscissa (Inter.first) of intersection of level line of index (Inter.second)
typedef std::pair<float,size_t> Inter;

void extract(const unsigned char* data, size_t w, size_t h,
             float offset, float step, int ptsPixel,
             std::vector<LevelLine*>& ll,
             std::vector< std::vector<Inter> >* inter=0);

#endif
