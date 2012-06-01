#include "lltree.h"
#include "fill_curve.h"
#include "gass.h"
#include "cmdLine.h"
#include "io_png.h"
#include <fstream>

/// Put one pixel wide blank strips at border of image
static void blank_border(unsigned char* data, size_t w, size_t h) {
    for(int i=1; i<h; i++) // Vertical strips
        data[i*w-1] = data[i*w] = 0;
    for(int i=0; i<w; i++) // First line
        data[i] = 0;
    data += (h-1)*w;
    for(int i=0; i<w; i++) // Last line
        data[i] = 0;
}

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

int main(int argc, char** argv) {
    CmdLine cmd;
    int ptsPixel=1;
    float offset=0.5f, step=10.0f;
    double lastScale=0;
    std::string imgOut;
    cmd.add( make_option('p', ptsPixel, "precision") );
    cmd.add( make_option('o', offset, "offset") );
    cmd.add( make_option('s', step, "step") );
    cmd.add( make_option('r', imgOut, "reconstruct") );
    cmd.add( make_option('l', lastScale) );
    try {
        cmd.process(argc, argv);
    } catch(std::string s) { std::cout << s << std::endl; }
    if(argc!=3) {
        std::cout << "Usage: " << argv[0] << ' '
                  << "[-p|--precision prec] "
                  << "[-o|--offset o] "
                  << "[-s|--step s] "
                  << "[-l lastScale]"
                  << "im.png lines.txt" <<std::endl;
        return 1;
    }

    // Input
    size_t w, h;
    unsigned char* data = read_png_u8_gray(argv[1], &w, &h);
    if(! data) {
        std::cout << "Impossible to read PNG image " << argv[1] <<std::endl;
        return 1;
    }

    // Work
    blank_border(data, w, h);
    LLTree tree(data,w,h, offset,step,ptsPixel);

    //Smooth
    if(lastScale>0)
        for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
            smooth(it->ll->line, lastScale);


    // Output
    std::ofstream file(argv[2]);
    for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
        file << *it->ll << "e" <<std::endl; // Required by megwave2's flreadasc
    file << "q" <<std::endl; // Required by megwave2's flreadasc

    if(!imgOut.empty()) {
        std::fill(data, data+w*h, 0);
        std::vector< std::vector<float> > inter;
        for(LLTree::iterator it=tree.begin(); it!=tree.end(); ++it)
            fill_curve(it->ll->line, (unsigned char)it->ll->level,
                       data,w,h, &inter);
        write_png_u8(imgOut.c_str(), data, w, h, 1);
    }

    free(data);
    return 0;
}
