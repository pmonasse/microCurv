#include "levelLine.h"
#include <sstream>

static const unsigned char image[4*4]= {
    0, 0, 0, 0,
    0, 0,12, 0,
    0, 5, 2, 0,
    0, 0, 0, 0
};
// Saddle value: 4

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

    std::vector<LevelLine*> ll;
    extract(image,4,4, level, 256.0f, 5, ll);

    for(std::vector<LevelLine*>::const_iterator it=ll.begin();it!=ll.end();++it)
        std::cout << **it << std::endl;
    return 0;
}
