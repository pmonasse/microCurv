#ifndef LEVELLINE_H
#define LEVELLINE_H

#include <vector>
#include <list>
#include <iostream>

struct Point {
    float x, y;
};

struct LevelLine {
    float level;
    std::vector<Point> line;
};

std::ostream& operator<<(std::ostream& str, const LevelLine& line);

void extract(const float* data, size_t w, size_t h, float offset, float step,
             std::list<LevelLine>& ll);

#endif
