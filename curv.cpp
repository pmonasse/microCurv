/**
 * @file curv.cpp
 * @brief Compute mean curvatures
 * @author Adina Ciomaga <adina@math.uchicago.edu>
 *         Pascal Monasse <monasse@imagine.enpc.fr>
 * 
 * Copyright (c) 2011-2016, Adina Ciomaga, Pascal Monasse
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

#include "curv.h"
#include "levelLine.h"
#include <algorithm>
#include <cmath>
#include <cassert>

/// Sign (+/-1) of \a x
#define SIGN(x) (((x)>0.0)?1: -1)

/// Minimum number of points in a curve to compute curvatures. The bare minimum
/// is 3, but small curves result in large curvatures, so we ignore them.
static const std::vector<float>::size_type MIN_PTS_CURV=10;

/// Distance between points.
static float dist(Point p, const Point& q) {
    p.x -= q.x; p.y -= q.y;
    return sqrt(p.x*p.x+p.y*p.y);
}

/// Compute 2x2 determinant: det(b-a c-a).
/// The minus sign of the result comes from the fact that the pixel coordinate
/// system is left-handed.
inline float det(const Point& a, const Point& b, const Point& c) {
    return -((b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x));
}

/// Find index of last point of the curve.
///
/// Take into account mulitple points: return the max index of point that is
/// different from point of index 0, or 0 if no such point exists.
static int last_point(const std::vector<Point>& curve) {
    assert(! curve.empty());
    Point p0 = curve.front();
    std::vector<Point>::const_reverse_iterator it =curve.rbegin(),
                                               end=curve.rend();
    size_t i = curve.size()-1;
    for(; it!=end; ++it, --i)
        if(*it!=p0)
            return i;
    return 0; // Single vertex
}

/// Record curvatures inside the pixels the \a curve goes through.
///
/// The orientation is said positive (wrt trigonometric orientation) if the
/// gradient of the image is to the left when following the level line.
/// \param curve is the polygonal curve
/// \param curvatures is the array to store in
/// \param w the width of the array \a curvatures.
static void curv(const std::vector<Point>& curve, int w,
                 std::vector<float>* curvatures)
{
    assert(curve.size()>=3); // Need at least 3 points to compute a curvature
    int n = last_point(curve);
    const Point* q = &curve[n]; 
    const Point* p = &curve[0];
    const Point* r = &curve[1];
    float u = dist(*p,*q);

    for(int k=0; k<n; k++,r++){
        if(k+1==n)
            r = &curve[0];
        float v = dist(*p,*r);
        float d = u * v * dist(*q,*r); 
        float K = ((d==0)? 0: 2*det(*p,*q,*r)/d);
        int i = (int)(curve[k].x);
        int j = (int)(curve[k].y);
        curvatures[j*w+i].push_back(K);
        q = p;
        p = r;
        u = v;
    }
}

/// Median value of array.
///
/// If there is an even number of elements, the average of the two middle ones
/// is returned.
float median(std::vector<float>& list) {
    assert(! list.empty());
    size_t size = list.size();
    std::vector<float>::iterator m=list.begin()+(size-1)/2;
    std::nth_element(list.begin(), m, list.end());
    float v = *m;
    if(size%2 == 0)
        v = (v+*std::min_element(m+1, list.end()))*0.5f;
    return v;
}

/// Compute the curvature inside each pixel.
///
/// This is the median of the curvatures of level lines at points inside the
/// pixel. Only level lines having more than \c MIN_PTS_CURV are considered,
/// and after taking median, each value is truncated in interval [-1,+1].
void curv(const std::vector<LevelLine*>& ll,
          float zoom, float* out, int w, int h) {
    // List of curvatures inside each pixel
    std::vector<float>* curvatures = new std::vector<float>[w*h];

    // Compute curvature
    std::vector<LevelLine*>::const_iterator it=ll.begin();
    for(int i=0; it!=ll.end(); ++it, ++i)
        if((*it)->line.size()>=MIN_PTS_CURV) {
            if(zoom==1.0f)
                curv((*it)->line, w, curvatures);
            else {
                std::vector<Point> line = (*it)->line;
                zoom_line(line, zoom);
                curv(line, w, curvatures);
            }
        }

    // Median curvature inside each pixel
    for(int i=0; i<w*h; i++, out++)
        if(! curvatures[i].empty()) {
            float m =  median(curvatures[i]);
            *out = SIGN(m)*std::min(std::abs(m),1.0f);
        }

    delete [] curvatures;
}
