// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file levelLine.h
 * @brief Extraction of level lines from an image
 * 
 * (C) 2011-2016, 2019, Pascal Monasse <pascal.monasse@enpc.fr>
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
