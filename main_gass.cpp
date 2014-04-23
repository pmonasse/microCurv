/**
 * @file main_gass.cpp
 * @brief Geometric affine smoothing of a curve
 * @author Pascal Monasse <monasse@imagine.enpc.fr>
 * 
 * Copyright (c) 2011-2014, Pascal Monasse
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

#include "gass.h"
#include "cmdLine.h"
#include <fstream>
#include <sstream>

int main(int argc, char** argv) {
    CmdLine cmd;
    double lastScale=0;
    cmd.add( make_option('l', lastScale) );
    cmd.process(argc, argv);
    if(argc!=3) {
        std::cout << argv[0] << ' ' << " [-l lastScale] flistIn flistOut" <<std::endl;
        return 1;
    }

    // Input
    std::ifstream in(argv[1]);
    std::vector< std::vector<DPoint> > ll;
    std::vector<DPoint> line;
    while(true) {
        DPoint d(0,0);
        std::streampos where = in.tellg();
        in >> d.x >> d.y;
        if(in.fail()) {
            in.clear();
            in.seekg(where);
            char c=0;
            in >> c;
            if(c=='e') {
                ll.push_back(line);
                std::cout << line.size() << ' ' << std::flush;
                line.clear();
            } else if(c=='q')
                break;
            else {
                std::cerr << "Problem reading at pos " << (int)where<<std::endl;
                return 1;
            }
        } else
            line.push_back(d);
    }
    std::cout << std::endl;

    // Process
    for(std::vector< std::vector<DPoint> >::iterator it=ll.begin();
        it!=ll.end(); ++it)
        gass(*it, 0, lastScale);

    // Output
    std::ofstream file(argv[2]);
    for(std::vector< std::vector<DPoint> >::iterator it=ll.begin();
        it!=ll.end(); ++it) {
        for(std::vector<DPoint>::iterator it2=it->begin(); it2!=it->end(); ++it2)
            file << it2->x << ' ' << it2->y << ' ';
        file << "e" <<std::endl; // Required by megwave2's flreadasc
    }
    file << "q" <<std::endl; // Required by megwave2's flreadasc

    return 0;
}
