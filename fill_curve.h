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
