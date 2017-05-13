/**
 * @file main.cpp
 * @brief Handle i/o of program microCurv.
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
#include "lltree.h"
#include "draw_curve.h"
#include "curv.h"
#include "image.h"
#include "cmdLine.h"
#include "xmtime.h"
#include "io_png.h"
#include "io_tiff.h"
#include <algorithm>
#include <fstream>
#include <cmath>
#include <limits>

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

/// The image is enlarged by this number of pixels on each side before
/// extraction of level lines, in order to reduce border effects.
static const int MARGIN=20;

/// Write in SVG file the level lines of \a tree at level multiple of \a qstep.
static bool output_svg(LLTree& tree, int qstep, int w, int h, Rect R,
                       const std::string& fileName) {
    std::ofstream file(fileName.c_str());
    file << "<svg xmlns=\"http://www.w3.org/2000/svg\" ";
    file << "width=\"" << w-R.x << "\" " << "height=\"" << h-R.y << "\" ";
    file << "viewBox=\"0 0 " << R.w << ' ' << R.h << "\">" << std::endl;
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
        if((int)it->ll->level%qstep == 0) {
            file << "<polygon stroke=\"blue\" stroke-width=\".1\" ";
            file << "fill=\"white\" points=\"";
            std::vector<Point>::const_iterator it2, end=it->ll->line.end();
            for(it2=it->ll->line.begin(); it2!=end; ++it2) {
                file << it2->x-R.x << ' ' << it2->y-R.y << ' ';
            }
            file << "\" />";
        }
    file << "</svg>" << std::endl;
    return file.good();
}

/// Write in PNG file the level lines of \a tree at level multiple of \a qstep.
/// \a z is a zoom factor.
static bool output_png(LLTree& tree, int qstep, int w, int h, Rect R, float z,
                       const std::string& fileName) {
    w = static_cast<int>( ceil(z*w) );
    h = static_cast<int>( ceil(z*h) );
    unsigned char* imgLL = new unsigned char[3*w*h];
    std::fill(imgLL, imgLL+3*w*h, 255); // White image
    // Set R and G channels of pixels on curves at 0, so that they become blue
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
        if((int)it->ll->level%qstep == 0)
            if(z==1.0f)
                draw_curve(it->ll->line, 0, imgLL, w, h);
            else {
                std::vector<Point> line = it->ll->line;
                zoom_line(line, z);
                draw_curve(line, 0, imgLL, w, h);                
            }
    std::copy(imgLL, imgLL+w*h, imgLL+w*h); // No need to redraw in G, copy R

    R.x = static_cast<int>( floor(z*R.x) );
    R.y = static_cast<int>( floor(z*R.y) );
    R.w = static_cast<int>( floor(z*R.w) );
    R.h = static_cast<int>( floor(z*R.h) );
    unsigned char* ll = crop(imgLL,w,h,R,3);
    bool ok = (io_png_write_u8(fileName.c_str(),ll,R.w,R.h,3)==0);
    delete [] ll;
    delete [] imgLL;
    return ok;
}

/// Color in blue the level lines of \a tree at a level multiple of \a qstep.
///
/// The dimensions of the image are \a w and \a h, but only the crop region \a R
/// is used to output to file \a fileName. If file name is empty, do nothing.
static bool output_curves(LLTree& tree, int qstep, int w, int h, Rect R,
                          float zoom,
                          const std::string& fileName) {
    if(fileName.empty()) return true;
    size_t end = fileName.rfind(".svg");
    bool ok = (end!=std::string::npos && fileName.size()==end+4)?
        output_svg(tree, qstep, w, h, R, fileName):
        output_png(tree, qstep, w, h, R, zoom, fileName);
    if(! ok)
        std::cerr << "Error writing image file " << fileName << std::endl;
    return ok;
}

/// Fonction mapping curvature value to intensity
inline unsigned char fct(float f) {
    return (unsigned char)(255*(1.0f-log(1000*f+1)/log(1000*1.0f+1)));
}

/// Display curvature map as color superimposed on faded image.
/// Red=negative, green=positive, yellow=zero.
bool colorCurv(const unsigned char* in, const float* map, int w, int h,
               const char* sFileName) {
    unsigned char *r=new unsigned char[3*w*h], *g=r+w*h, *b=g+w*h;
    std::fill(r, r+3*w*h, 0);
    for(int i=0; i<w*h; i++, in++, map++) {
        if(*map==255.0f)
            r[i]=g[i]=b[i] = (unsigned char)(128+0.5*(*in-128));
        else {
            r[i]=g[i]= fct(0);
            if(*map>=0)
                r[i] = fct(+*map);
            else
                g[i] = fct(-*map);
        }
    }
    bool ok = (io_png_write_u8(sFileName,r,w,h,3)==0);
    delete [] r;
    return ok;
}

/// Main procedure for curvature microscope.
int main(int argc, char** argv) {
    int qstep=8;
    float scale=2.0f;
    float zoom=1.0f;
    Rect rectSelect;
    std::string inLL, outLL, sOutImage;
    CmdLine cmd;
    cmd.add( make_option('q', qstep) );
    cmd.add( make_option('s', scale) );
    cmd.add( make_option('z', zoom) );
    cmd.add( make_option('I', inLL) );
    cmd.add( make_option('O', outLL) );
    cmd.add( make_option('o', sOutImage) );
    cmd.add( make_option('r', rectSelect) );
    cmd.process(argc, argv);
    if(argc!=3 && argc!=4) {
        std::cerr << "Usage: " << argv[0] << ' '
                  << "[-q qstep] [-s scale] [-z zoom] [-r rect] "
                  << "[-I inLL.<png|svg>] [-O outLL.<png|svg>] "
                  << "[-o outImage.png] "
                  << "in.png outCurv.png [outCurv.tif]" << std::endl
                  << "\trect: wxh+x+y" << std::endl;
        return 1;
    }

    size_t w, h;
    unsigned char* inIm = io_png_read_u8_gray(argv[1], &w, &h);
    if(! inIm) {
        std::cerr << "Error reading as PNG image: " << argv[1] << std::endl;
        return 1;
    }

    Timer timer;

    if(! cmd.used('r')) {
        rectSelect.x = rectSelect.y = 0;
        rectSelect.w = (int)w;
        rectSelect.h = (int)h;
    }
    const Rect R={MARGIN, MARGIN, rectSelect.w, rectSelect.h}; // Image ROI
    rectSelect.x -= MARGIN;
    rectSelect.y -= MARGIN;
    rectSelect.w += 2*MARGIN;
    rectSelect.h += 2*MARGIN;
    const int ncol=rectSelect.w, nrow=rectSelect.h; // Dim of image with margins

    unsigned char* inImage = extract(inIm,w,h, rectSelect);
    free(inIm);
    fill_border(inImage,ncol,nrow);

    std::cout << " 1. Extract level lines: " << std::flush;
    timer.tick();
    LLTree* tree = extract_tree(inImage, ncol, nrow);
    std::cout << tree->nodes().size() << " level lines. ";
    timer.time();

    delete [] inImage;
    if(! output_curves(*tree,qstep,ncol,nrow,R,zoom,inLL))
        return 1;

    std::cout << " 2. Smooth level lines by affine shortening. " << std::flush;
    timer.tick();
    smooth_ll(*tree, scale);
    timer.time();

    if(! output_curves(*tree,qstep,ncol,nrow,R,zoom,outLL))
        return 1;

    std::cout << " 3. Reconstruct smoothed image. " << std::flush;
    timer.tick();
    int zncol=ncol, znrow=nrow; Rect zR=R;
    unsigned char* llas = reconstruct(*tree, zncol, znrow, zR, zoom);
    timer.time();

    if(!sOutImage.empty() &&
       io_png_write_u8(sOutImage.c_str(), llas, zR.w, zR.h, 1)!=0) {
        std::cerr << "Error writing image file " << sOutImage << std::endl;
        return 1;
    }

    std::cout << " 4. Compute the curvature map: " << std::flush;
    timer.tick();
    std::vector<LevelLine*> qll;
    quantize(*tree, qstep, qll);
    std::cout << " based on " << qll.size() << " level lines. " << std::flush;
    float* outCurv = new float[znrow*zncol];
    std::fill(outCurv, outCurv+znrow*zncol, 255.0f);
    curv(qll, zoom, outCurv,zncol,znrow);
    timer.time();

    float* cmap = crop(outCurv,zncol,znrow, zR); // Curvature map
    delete [] outCurv;

    if(! colorCurv(llas, cmap,zR.w,zR.h, argv[2])) {
        std::cerr << "Error writing image file " << argv[2] << std::endl;
        return 1;
    }
    if(argc>3) {
        for(int i=zR.w*zR.h-1; i>=0; i--)
            if(cmap[i]==255.0f)
                cmap[i] = std::numeric_limits<float>::quiet_NaN();
        if(io_tiff_write_f32(argv[3], cmap,zR.w,zR.h,1)!=0) {
            std::cerr << "Error writing image file " << argv[3] << std::endl;
            return 1;
        }
    }

    delete tree;
    delete [] llas;
    delete [] cmap;

    return 0;
}
