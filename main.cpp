#include "levelLine.h"
#include "io_png.h"
#include <sstream>
#include <fstream>

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
    std::list<LevelLine> ll;
    extract(data, w, h, offset, step, ll);

    // Output
    std::ofstream file(argv[4]);
    std::list<LevelLine>::const_iterator it;
    for(it=ll.begin(); it!=ll.end(); ++it)
        file << *it << "e" <<std::endl; //as required by megwave2's flreadasc
    file << "q" <<std::endl; //as required by megwave2's flreadasc

    return 0;
}
