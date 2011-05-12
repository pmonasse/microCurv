#include "levelLine.h"

/// Print list of points of the line (but not the level): x0 y0 x1 y1...
std::ostream& operator<<(std::ostream& str, const LevelLine& l) {
    std::vector<Point>::const_iterator it;
    for(it=l.line.begin(); it!=l.line.end(); ++it)
        str << it->x << " " << it->y << " ";
    return str;
}

/// Level lines extraction algorithm
void extract(const float* data, size_t w, size_t h, float offset, float step,
             std::list<LevelLine>& ll) {
}
