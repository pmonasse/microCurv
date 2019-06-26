/**
 * @file bilines.cpp
 * @brief Display bilinear level lines.
 * 
 * Copyright (c) 2019, Pascal Monasse
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

#include "lltree.h"
#include "draw_curve.h"
#include "fill_curve.h"
#include "cmdLine.h"
#include "xmtime.h"
#include "io_png.h"
#include <map>
#include <cmath>

/// Timer class to measure real time (not CPU time)
class Timer {
    unsigned long t; ///< Current time in milliseconds
public:
    Timer() { tick(); } ///< Constructor
    unsigned long tick() { return t=xmtime(); } ///< Reset time
    void time() { ///< Display elapsed time and reset current time
        unsigned long told = t;
        std::cout << "Time = " << (tick()-told)/1000.0f << "s" << std::endl;
    }
};

/// Compute histogram of level at pixels at the border of the image.
static void histogram(unsigned char* im, size_t w, size_t h, size_t histo[256]){
    size_t j;
    for(j=0; j<w; j++) // First line
        ++histo[im[j]];
    for(size_t i=1; i+1<h; i++) { // All lines except first and last
        ++histo[im[j]];  // First pixel of line
        j+= w-1;
        ++histo[im[j++]]; // Last pixel of line
    }
    for(; j<w*h; j++) // Last line
        ++histo[im[j]];    
}

/// Put pixels at border of image to value \a v.
static void put_border(unsigned char* im, size_t w, size_t h, unsigned char v) {
    size_t j;
    for(j=0; j<w; j++)
        im[j] = v;
    for(size_t i=1; i+1<h; i++) {
        im[j] = v;
        j+= w-1;
        im[j++] = v;
    }
    for(; j<w*h; j++)
        im[j] = v;
}

/// Set all pixels at border of image to their median level.
static unsigned char fill_border(unsigned char* im, size_t w, size_t h) {
    size_t histo[256] = {0}; // This puts all values to zero
    histogram(im, w, h, histo);
    size_t limit=w+h-2; // Half number of pixels at border
    size_t sum=0;
    int i=-1;
    while((sum+=histo[++i]) < limit);
    put_border(im,w,h, (unsigned char)i);
    return (unsigned char)i;
}

/// Find upper/lower level sets and shift level accodingly.
void fix_levels(LLTree& tree, unsigned char bg, int qStep) {
    std::map<LLTree::Node*,bool> upper;
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it) {
        float parentLevel = it->parent? it->parent->ll->level: bg;
        bool up = it->ll->level > parentLevel;
        if(it->ll->level == parentLevel)
            up = !upper[it->parent];
        upper[&*it] = up;
    }
    float delta = 0.5f*(float)qStep;
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it) {
        if(upper[&*it])
            it->ll->level = std::min(it->ll->level+delta,255.0f);
        else
            it->ll->level = std::max(it->ll->level-delta,0.0f);
    }
}

/// Return depth of node in tree. Roots are at depth 0.
static int depth(const LLTree::Node& node) {
    const LLTree::Node* n=&node;
    int d=0;
    while(n->parent) {
        n = n->parent;
        ++d;
    }
    return d;
}

/// Palette 'rainbow' of gnuplot, from purple to red through blue and yellow.
static void palette(float x,
                    unsigned char& r, unsigned char& g, unsigned char& b) {
    r = (unsigned char)(255*std::min(1.0f, std::abs(2*x-.5f)));
    g = (unsigned char)(255*std::sin(M_PI*x));
    b = (unsigned char)(255*std::cos(M_PI_2*x));
}

/// Main procedure for curvature microscope.
int main(int argc, char** argv) {
    int qstep=32;
    CmdLine cmd; cmd.prefixDoc = "\t";
    cmd.add( make_option('q',qstep).doc("Quantization step (integer)") );
    cmd.process(argc, argv);
    if(argc!=3) {
        std::cerr << "Usage: " << argv[0]
                  << " [options] in.png out.png" << std::endl;
        std::cerr << "Option:\n" << cmd;
        return 1;
    }

    size_t w, h;
    unsigned char* in = io_png_read_u8_gray(argv[1], &w, &h);
    if(! in) {
        std::cerr << "Error reading as PNG image: " << argv[1] << std::endl;
        return 1;
    }
    unsigned char bg = fill_border(in, w, h); // Background gray of output

    Timer timer;

    std::cout << " 1. Extract level lines: " << std::flush;
    timer.tick();
    const float offset=0.5f;
    LLTree tree(in, (int)w, (int)h, offset, qstep, 0);
    free(in);
    std::cout << tree.nodes().size() << " level lines. ";
    timer.time();

    fix_levels(tree, bg, qstep);

    std::cout << " 2. Reconstruct from level lines. " << std::flush;
    unsigned char* out = new unsigned char[3*w*h];
    timer.tick();
    std::fill(out, out+3*w*h, bg);
    std::vector< std::vector<float> > inter;
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
        fill_curve(it->ll->line,(unsigned char)it->ll->level,
                   out,(int)w,(int)h, &inter);
    std::copy(out, out+w*h, out+1*w*h); // Copy to green channel
    std::copy(out, out+w*h, out+2*w*h); // Copy to red channel
    timer.time();

    int depthMax=0;
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it) {
        int d = depth(*it);
        if(depthMax<d)
            depthMax = d;
    }
    std::cout << "Max depth: " << depthMax << std::endl;

    std::cout << " 3. Draw curves. " << std::flush;
    timer.tick();
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it) {
        int d = depth(*it);
        unsigned char r,g,b;
        palette(d/(float)depthMax, r,g,b);
        draw_curve(it->ll->line,r, out+0*w*h,(int)w,(int)h);
        draw_curve(it->ll->line,g, out+1*w*h,(int)w,(int)h);
        draw_curve(it->ll->line,b, out+2*w*h,(int)w,(int)h);
    }
    timer.time();

    if(io_png_write_u8(argv[2], out, (int)w, (int)h, 3)!=0) {
        std::cerr << "Error writing image file " << argv[2] << std::endl;
        return 1;
    }
    delete [] out;

    return 0;
}
