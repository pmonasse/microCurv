/**
 * @file microCurv.cpp
 * @brief Compute mean curvatures
 * @author Adina Ciomaga <adina@math.uchicago.edu>
 *         Pascal Monasse <monasse@imagine.enpc.fr>
 * 
 * Copyright (c) 2011-2014, Adina Ciomaga, Pascal Monasse
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
#include "gass.h"
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
    void tick() { t=xmtime(); } ///< Reset time
    void time() { ///< Display elapsed time and reset current time
        unsigned long told = t;
        tick();
        std::cout << "Time = " << (t-told)/1000.0f << "s" << std::endl;
    }
};

/// The image is enlarged by this number of pixels on each side before
/// extraction of level lines, in order to reduce border effects.
static const int MARGIN=20;

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

/// Extract level lines from \a tree at levels multiple of \a qstep.
///
/// Output is in \a qll. In addition, orientations are stored in \a positive.
static void quantize(LLTree& tree, int qstep,
                     std::vector<LevelLine*>& qll, std::vector<bool>& positive){
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
		if((int)it->ll->level%qstep == 0) {
            qll.push_back(it->ll);
            bool up = (it->parent==0 || it->parent->ll->level<it->ll->level);
            positive.push_back(up);
        }
}

/// Write in SVG file the level lines of \a tree at level multiple of \a qstep.
static bool output_svg(LLTree& tree, int qstep, int w, int h, Rect R,
                       const std::string& fileName) {
    std::ofstream file(fileName.c_str());
    file << "<svg xmlns=\"http://www.w3.org/2000/svg\" ";
    file << "width=\"" << w-R.x << "\" " << "height=\"" << h-R.y << "\" ";
    file << "viewBox=\"0 0 " << R.w << ' ' << R.h << "\">" << std::endl;
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
        if((int)it->ll->level%qstep == 0) {
            file << "<polygon stroke=\"blue\" stroke-width=\".5\" ";
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
                std::vector<Point> line;
                zoom_line(it->ll->line, line, z);
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
/// is used to output to file \a fileName.
static bool output_curves(LLTree& tree, int qstep, int w, int h, Rect R,
                          float zoom,
                          const std::string& fileName) {
    size_t end = fileName.rfind(".svg");
    bool ok = (end!=std::string::npos && fileName.size()==end+4)?
        output_svg(tree, qstep, w, h, R, fileName):
        output_png(tree, qstep, w, h, R, zoom, fileName);
    if(! ok)
        std::cerr << "Error writing image file " << fileName << std::endl;
    return ok;
}

/// Reconstruct image from level sets stored in \a tree.
static unsigned char* reconstruct(LLTree& tree, int& w, int& h, Rect& R,
                                  float zoom) {
    w = static_cast<int>( ceil(zoom*w) );
    h = static_cast<int>( ceil(zoom*h) );
    unsigned char* outImage = new unsigned char[w*h];
    std::fill_n(outImage, w*h, 0);
    std::vector< std::vector<float> > inter;
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
        if(zoom==1.0f)
            fill_curve(it->ll->line,(unsigned char)it->ll->level,
                       outImage,w,h, &inter);
        else {
            std::vector<Point> line;
            zoom_line(it->ll->line, line, zoom);
            fill_curve(line,(unsigned char)it->ll->level,
                       outImage,w,h, &inter);
        }
    R.x = static_cast<int>( floor(zoom*R.x) );
    R.y = static_cast<int>( floor(zoom*R.y) );
    R.w = static_cast<int>( floor(zoom*R.w) );
    R.h = static_cast<int>( floor(zoom*R.h) );
    unsigned char* out = crop(outImage,w,h, R);
    delete [] outImage;
    return out;
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
    std::string inLL, outLL, sOutImage;
    CmdLine cmd;
    cmd.add( make_option('q', qstep) );
    cmd.add( make_option('s', scale) );
    cmd.add( make_option('z', zoom) );
    cmd.add( make_option('I', inLL) );
    cmd.add( make_option('O', outLL) );
    cmd.add( make_option('o', sOutImage) );
    cmd.process(argc, argv);
    if(argc!=3 && argc!=4) {
        std::cerr << "Usage: " << argv[0] << ' '
                  << "[-q qstep] [-s scale] [-z zoom] "
                  << "[-I inLL.<png|svg>] [-O outLL.<png|svg>] "
                  << "[-o outImage.png] "
                  << "in.png outCurv.png [outCurv.tif]" << std::endl;
        return 1;
    }

    size_t w, h;
    unsigned char* inIm = io_png_read_u8_gray(argv[1], &w, &h);
    if(! inIm) {
        std::cerr << "Error reading as PNG image: " << argv[1] << std::endl;
        return 1;
    }

    Timer timer;

    unsigned char* inImage = sym_enlarge(inIm,w,h, MARGIN);
    free(inIm);
    const Rect R={MARGIN, MARGIN, w, h};
    const int ncol=w+2*MARGIN, nrow=h+2*MARGIN;

    std::cout << " 1. Extract level lines. " << std::flush;
    const float offset=0.5f;
    LLTree tree(inImage, ncol, nrow, offset, 1.0f, 5);
    fix_level(tree, offset);
    delete [] inImage;
    if(!inLL.empty() && !output_curves(tree,qstep,ncol,nrow,R,zoom,inLL))
        return 1;

    std::vector<LevelLine*> qll;
    std::vector<bool> positive;
    quantize(tree, qstep, qll, positive);
    std::cout << qll.size() << "/" <<tree.nodes().size() << " level lines. ";
    timer.time();

    std::cout << " 2. Smooth level lines by affine shortening. " <<std::flush;
    if(scale>0) {
        const int size=(int)tree.nodes().size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(int i=0; i<size; i++)
            smooth(tree.nodes()[i].ll->line, scale);
    }
    if(!outLL.empty() && !output_curves(tree,qstep,ncol,nrow,R,zoom,outLL))
        return 1;
    timer.time();

    std::cout << " 3. Reconstruct LLAS evolution. " << std::flush;
    int zncol=ncol, znrow=nrow;
    Rect zR=R;
    unsigned char* llas = reconstruct(tree, zncol, znrow, zR, zoom);
    if(!sOutImage.empty() && io_png_write_u8(sOutImage.c_str(),
                                             llas, zR.w, zR.h, 1)!=0){
        std::cerr << "Error writing image file " << sOutImage << std::endl;
        return 1;
    }
    timer.time();

    std::cout << " 4. Construct the curvature map. " << std::flush;
    float* outCurv = new float[znrow*zncol];
    std::fill(outCurv, outCurv+znrow*zncol, 255.0f);
    curv(qll, positive, zoom, outCurv,zncol,znrow);
    float* cmap = crop(outCurv,zncol,znrow, zR);
    delete [] outCurv;

    if(! colorCurv(llas, cmap,zR.w,zR.h, argv[2])) {
        std::cerr << "Error writing image file " << argv[2] << std::endl;
        return 1;
    }
    for(int i=znrow*zncol-1; i>=0; i--)
        if(cmap[i]==255.0f)
            cmap[i] = std::numeric_limits<float>::quiet_NaN();
    if(argc>3 && io_tiff_write_f32(argv[3], cmap,zR.w,zR.h,1)!=0) {
        std::cerr << "Error writing image file " << argv[3] << std::endl;
        return 1;
    }
    timer.time();

    delete [] llas;
    delete [] cmap;

    return 0;
}
