#include "lltree.h"
#include "draw_curve.h"
#include "fill_curve.h"
#include "gass.h"
#include "curv.h"
#include "mirror.h"
#include "cmdLine.h"
#include "io_png.h"
#include "io_tiff.h"
#include <algorithm>
#include <cmath>
#include <ctime>

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

/// Crop rectangle wxh+x+y from image
template <typename T>
static T* crop(const T* in, size_t wIn, size_t hIn, int x, int y, int w, int h,
               int ch=1) {
    T* out = new T[ch*w*h];
    for(int j=0; j<ch; j++, in+=wIn*hIn)
        for(int i=0; i<h; i++) {
            const T* p = in+wIn*(y+i)+x;
            std::copy(p, p+w, out+j*w*h+i*w);
        }
    return out;
}

/// Fonction mapping curvature value to intensity
inline unsigned char fct(float f) {
    return (unsigned char)(255*(1.0f-log(1000*f+1)/log(1000*1+1)));
}

/// Display curvature map as color superimposed on faded image
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
    bool ok = (write_png_u8(sFileName,r,w,h,3)==0);
    delete [] r;
    return ok;
}

int main(int argc, char** argv) {
    int qstep=8;
    float last=2.0f;
    std::string inLL, outLL, sOutImage;
    CmdLine cmd;
    cmd.add( make_option('q', qstep) );
    cmd.add( make_option('l', last) );
    cmd.add( make_option('I', inLL) );
    cmd.add( make_option('O', outLL) );
    cmd.add( make_option('o', sOutImage) );
    cmd.process(argc, argv);
    if(argc!=3 && argc!=4) {
        std::cerr << "Usage: " << argv[0] << ' '
                  << "[-q qstep] [-l last] [-I inLL] [-O outLL] [-o outImage] "
                  << "in.png outCurv.png [outCurv.tif]" << std::endl;
        return 1;
    }

    size_t w, h;
    unsigned char* inIm = read_png_u8_gray(argv[1], &w, &h);
    if(! inIm) {
        std::cerr << "Unable to open image " << argv[1] << std::endl;
        return 1;
    }

    std::clock_t t = std::clock();

    int margin=20;
    unsigned char* inImage = sym_enlarge(inIm,w,h, margin);
    free(inIm);
    const int nrow=w+2*margin, ncol=h+2*margin;

    unsigned char* imgLL=0;
    if(! inLL.empty() || ! outLL.empty())
        imgLL = new unsigned char[3*nrow*ncol];

    std::cout << " 1. Extract level lines. " << std::flush;
    float offset=0.5f;
    LLTree tree(inImage, nrow, ncol, offset, 1.0f, 5);
    fix_level(tree, offset);
    delete [] inImage;
    if(! inLL.empty()) {
        std::fill(imgLL, imgLL+3*nrow*ncol, 255);
        for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
            if((int)it->ll->level%qstep == 0)
                draw_curve(it->ll->line, 0, imgLL, ncol, nrow);
        std::copy(imgLL, imgLL+ncol*nrow, imgLL+ncol*nrow);
        unsigned char* ll = crop(imgLL,ncol,nrow, margin,margin,w,h,3);
        if(! write_png_u8(inLL.c_str(),ll,w,h,3)==0) {
            std::cerr << "Error writing image file " << inLL << std::endl;
            return 1;
        }
        delete [] ll;
    }
    // Quantize level lines
    std::vector<LevelLine*> qll;
    std::vector<bool> positive;
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
		if((int)it->ll->level%qstep == 0) {
            qll.push_back(it->ll);
            bool up = (it->parent==0 || it->parent->ll->level<it->ll->level);
            positive.push_back(up);
        }
    t = std::clock()-t;
    std::cout << qll.size() << "/" <<tree.nodes().size() << " level lines. ";
    std::cout << "Time = " << t/(float)CLOCKS_PER_SEC << std::endl;

    std::cout << " 2. Evolve level lines by affine shortening. " <<std::flush;
    if(last>0) {
        const size_t size=tree.nodes().size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for(size_t i=0; i<size; i++)
            smooth(tree.nodes()[i].ll->line, last);
        //        for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
        //            smooth(it->ll->line, last);
    }
    if(! outLL.empty()) {
        std::fill(imgLL, imgLL+3*nrow*ncol, 255);
        for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
            if((int)it->ll->level%qstep == 0)
                draw_curve(it->ll->line, 0, imgLL, ncol, nrow);
        std::copy(imgLL, imgLL+ncol*nrow, imgLL+ncol*nrow);
        unsigned char* ll = crop(imgLL,ncol,nrow, margin,margin,w,h,3);
        if(! write_png_u8(outLL.c_str(),ll,w,h,3)==0) {
            std::cerr << "Error writing image file " << outLL << std::endl;
            return 1;
        }
        delete [] ll;
    }
    delete [] imgLL;
    t = std::clock()-t;
    std::cout << "Time = " << t/(float)CLOCKS_PER_SEC << std::endl;

    std::cout << " 3. Reconstruct LLAS evolution. " << std::flush;
    unsigned char* outImage = new unsigned char[nrow*ncol];
    std::fill_n(outImage, nrow*ncol, 0);
    std::vector< std::vector<float> > inter;
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
        fill_curve(it->ll->line,(unsigned char)it->ll->level,
                   outImage,ncol,nrow, &inter);
    unsigned char* llas = crop(outImage,ncol,nrow, margin,margin,w,h);
    if(! sOutImage.empty() && write_png_u8(sOutImage.c_str(), llas,w,h,1)!=0) {
        std::cerr << "Error writing image file " << sOutImage << std::endl;
        return 1;
    }
    delete [] outImage;
    t = std::clock()-t;
    std::cout << "Time = " << t/(float)CLOCKS_PER_SEC << std::endl;

    std::cout << " 4. Construct the curvature map. " << std::flush;
    float* outCurv = new float[nrow*ncol];
    std::fill(outCurv, outCurv+nrow*ncol, 255.0f);
    curv(qll, positive, outCurv,ncol,nrow);
    float* cmap = crop(outCurv,ncol,nrow, margin, margin, w, h);
    if(! colorCurv(llas, cmap,w,h, argv[2])) {
        std::cerr << "Error writing image file " << argv[2] << std::endl;
        return 1;
    }
    if(argc>3 && write_tiff_f32(argv[3], cmap,w,h,1)!=0) {
        std::cerr << "Error writing image file " << argv[3] << std::endl;
        return 1;
    }
    t = std::clock()-t;
    std::cout << "Time = " << t/(float)CLOCKS_PER_SEC << std::endl;

    delete [] llas;
    delete [] cmap;
    delete [] outCurv;

    return 0;
}
