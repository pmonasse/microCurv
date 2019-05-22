#include "levelLine.h"
#include <sstream>

static const unsigned char image1[4*4]= {
    0, 0, 0, 0,
    0, 0,12, 0,
    0, 5, 2, 0,
    0, 0, 0, 0
};
// Saddle value: 4

static const unsigned char image2[4*4]= {
   20,20,20,20,
   20, 0,12,20,
   20, 5, 2,20,
   20,20,20,20
};
// Saddle value: 4

static const unsigned char image3[4*4]= {
    0, 0, 0, 0,
    0, 0,12, 0,
    0, 2, 8, 0,
    0, 0, 0, 0
};
// Saddle value: 4, saddle point outside

void test_image(const unsigned char* im, size_t w, size_t h, float level) {
    for(int i=0; i<h; i++) {
        for(int j=0; j<w; j++)
            std::cout << "\t" << (int)im[i*w+j];
        std::cout << std::endl;
    }
    std::vector<LevelLine*> ll;
    extract(im,w,h, level, 256.0f, 5, ll);
    for(std::vector<LevelLine*>::iterator it=ll.begin();it!=ll.end();++it) {
        std::cout << **it << std::endl;
        delete *it;
    }
}

int main(int argc, char** argv) {
    if(argc>2) {
        std::cerr << "Usage: " << argv[0] << " level" << std::endl;
        return 1;
    }
    float level=4.0f;
    if(argc>1) {
        std::istringstream str(argv[1]);
        str >> level;
        if(str.fail()) {
            std::cerr<< "Parameter must be a floating point number" <<std::endl;
            return 1;
        }
    } else std::cout << "Default level: " << level << std::endl;

    std::cout << "image 1: saddle level at 4, two lines" << std::endl;
    test_image(image1,4,4, level);

    std::cout << "image 2: saddle level at 4, single line" << std::endl;
    test_image(image2,4,4, level);

    std::cout << "image 3: saddle level at 4, saddle pt outside" << std::endl;
    test_image(image3,4,4, level);

    return 0;
}
