#include "levelLine.h"
#include "io_png.h"
#include <sstream>
#include <fstream>
#include <cstdlib>

/// Put one pixel wide blank strips at border of image
static void blank_border(float* data, size_t w, size_t h) {
    for(int i=1; i<h; i++) // Vertical strips
        data[i*w-1] = data[i*w] = 0;
    for(int i=0; i<w; i++) // First line
        data[i] = 0;
    data += (h-1)*w;
    for(int i=0; i<w; i++) // Last line
        data[i] = 0;
}

int main(int argc, char** argv) {
    if(argc != 5) {
        std::cout << "Usage: " << argv[0]
                  << " im.png offset step lines.txt"<<std::endl;
        return 1;
    }

    // Input
    size_t w, h;
    float* data = read_png_f32_gray(argv[1], &w, &h);
    if(! data) {
        std::cout << "Impossible to read PNG image " << argv[1] <<std::endl;
        return 1;
    }

    float offset, step;
    if((std::istringstream(argv[2])>>offset).fail()) {
        std::cout <<"Unable to interpret "<<argv[2]<< " as number" <<std::endl;
        return 1;
    }
    if((std::istringstream(argv[3])>>step).fail()) {
        std::cout <<"Unable to interpret "<<argv[3]<< " as number" <<std::endl;
        return 1;
    }

    // Work
    blank_border(data, w, h);
    std::list<LevelLine> ll;
    extract(data, w, h, offset, step, ll);

    // Output
    std::ofstream file(argv[4]);
    std::list<LevelLine>::const_iterator it;
    for(it=ll.begin(); it!=ll.end(); ++it)
        file << *it << "e" <<std::endl; //as required by megwave2's flreadasc
    file << "q" <<std::endl; //as required by megwave2's flreadasc


    free(data);
    return 0;
}
