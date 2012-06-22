#ifndef DRAWCURVE_H
#define DRAWCURVE_H

#include "levelLine.h"

void draw_curve(const std::vector<Point>& curve, unsigned char v,
                unsigned char* im, int w, int h);

#endif
