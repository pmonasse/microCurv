/**
 * @file microCurv.cpp
 * @brief Compute mean curvatures
 * @author Adina Ciomaga <adina@math.uchicago.edu>
 *         Pascal Monasse <monasse@imagine.enpc.fr>
 * 
 * Copyright (c) 2011-2017, Adina Ciomaga, Pascal Monasse
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

#include "microCurv.h"
#include "fill_curve.h"
#include "gass.h"
#include "image.h"
#include <algorithm>
#include <cmath>

/// Shift the level of each line: up for positive line, down for negative
static void fix_level(LLTree& tree, float offset) {
    for(LLTree::iterator it=tree.begin(PostOrder); it!=tree.end(); ++it) {
        float level=it->ll->level;
        float d=offset;
        if(it->parent && level<it->parent->ll->level) // Negative line
            d=-d;
        it->ll->level += d;
        assert(0<=it->ll->level && it->ll->level<=255);
    }
}

/// Extract all level lines from image with quantization 1.
LLTree* extract_tree(const unsigned char* im, int w, int h){
    const float offset=0.5f;
    LLTree* tree = new LLTree(im, w, h, offset, 1.0f, 5);
    fix_level(*tree, offset);
    return tree;
}

/// Extract level lines from \a tree at levels multiple of \a qstep.
///
/// Output is in \a qll.
void quantize(LLTree& tree, int qstep, std::vector<LevelLine*>& qll) {
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
		if((int)it->ll->level%qstep == 0)
            qll.push_back(it->ll);
}

/// Reconstruct (zoomed) image from level sets stored in \a tree.
///
/// The zoom factor is applied to \a w, \a h and \a R.
unsigned char* reconstruct(LLTree& tree, int& w, int& h, Rect& R, float zoom) {
    w = static_cast<int>( std::ceil(zoom*w) );
    h = static_cast<int>( std::ceil(zoom*h) );
    unsigned char* outImage = new unsigned char[w*h];
    std::fill_n(outImage, w*h, 0);
    std::vector< std::vector<float> > inter;
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
        if(zoom==1.0f)
            fill_curve(it->ll->line,(unsigned char)it->ll->level,
                       outImage,w,h, &inter);
        else {
            std::vector<Point> line = it->ll->line;
            zoom_line(line, zoom);
            fill_curve(line,(unsigned char)it->ll->level,
                       outImage,w,h, &inter);
        }
    R.x = static_cast<int>( std::floor(zoom*R.x) );
    R.y = static_cast<int>( std::floor(zoom*R.y) );
    R.w = static_cast<int>( std::floor(zoom*R.w) );
    R.h = static_cast<int>( std::floor(zoom*R.h) );
    unsigned char* out = crop(outImage,w,h, R);
    delete [] outImage;
    return out;
}

/// Smooth level line
static void smooth(std::vector<Point>& line, double lastScale) {
    std::vector<DPoint> dline;
    for(std::vector<Point>::iterator it=line.begin(); it!=line.end(); ++it)
        dline.push_back( DPoint((double)it->x,(double)it->y) );
    assert(dline.front()==dline.back());
    gass(dline, 0.0, lastScale);
    line.clear();
    for(std::vector<DPoint>::iterator it=dline.begin(); it!=dline.end(); ++it)
        line.push_back( Point((float)it->x,(float)it->y) );
}

/// Smooth each level line in \a tree by affine shortening, up to \a scale.
void smooth_ll(LLTree& tree, float scale) {
    if(scale<=0) return;
    const int size=(int)tree.nodes().size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<size; i++)
        smooth(tree.nodes()[i].ll->line, scale);
}
