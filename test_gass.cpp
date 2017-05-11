#include "gass.h"

// This program tests a rare situation in split_convex that resulted in a bug
// (inifinite loop).
// There is an inflection point between 1 and 2, but 0, 1 and 3 are aligned.
// When it compares orientations of 2-3-0 and 3-0-1, since 3-0-1 is flat, it
// tries 3-0-2, which has the same orientation as 2-3-0. Therefore, the algo
// can begin with an inflection point and not find it as an inflection point
// when it loops.

int main() {
    std::vector<DPoint> line;
    line.push_back( DPoint(0,0) );
    line.push_back( DPoint(1,0) );
    line.push_back( DPoint(2,1) );
    line.push_back( DPoint(3,0) );
    line.push_back( line.front());
    gass(line, 0.0, 2.0);
    return 0;
}
