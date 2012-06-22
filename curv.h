#ifndef CURVE_H
#define CURVE_H

#include <vector>
class LevelLine;

void curv(const std::vector<LevelLine*>& ll, const std::vector<bool>& positive,
          float* out, int w, int h);

#endif
