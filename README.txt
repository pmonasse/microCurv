# Extraction of the Level Lines of a Bilinear Image #

## Summary ##
This software computes the level lines of a bilinear interpolated image. It is
linked to an IPOL article [1]. The level lines are organized in a tree structure
showing how they enclose one another. This is at the basis of the image
curvature microscope [2].

[1] Extraction of the Level Lines of a Bilinear Image, IPOL [URL]

[2] IPOL article, <https://doi.org/10.5201/ipol.2017.212>

## Author ##
Pascal Monasse <pascal.monasse@enpc.fr>

Laboratoire d'Informatique Gaspard Monge (LIGM)/
Ecole des Ponts ParisTech

## Version ##
Version 1.0 released on 2019/06/26.

Future releases and updates:
<https://github.com/pmonasse/microCurv.git>

## Build ##
Prerequisites: CMake version 2.6 or later

- Unix, MacOS:
    $ cd /pathToThisFile/
    $ mkdir Build && cd Build && cmake -DCMAKE_BUILD_TYPE:bool=Release ..
    $ make

- Windows with MinGW:
    $ cd /path_to_this_file/
    $ mkdir Build 
    $ cd Build
    $ cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE:bool=Release ..
    $ mingw32-make

CMake tries to find libPNG on your system, otherwise it uses the embedded version. The library needs to come with header files, which are provided by "lipng-dev" package under Linux Debian or Ubuntu.

- Test:
    $ ./bilines ../data/flower.png out.png
Compare out.png with ../data/flowerOut.png

## Usage ##
bilines [options] in.png out.png
  - in.png: PNG input image
  - out.png: PNG output image
Option:
  -q <int>: Quantization step (default 32)

## Files (Only those with &ast; are reviewed through IPOL) ##

    bilines.cpp (*): main function, calls Algo. 1 and implements Algo. 6.
    fill_curve.cpp (*) fill_curve.h (*): implements Algo. 5.
    levelLine.cpp (*)  levelLine.h (*): implements Algo. 4.
    lltree.cpp (*)     lltree.h (*): implements Algos. 1-3.

    CMakeLists.txt
    cmdLine.h
    draw_curve.cpp draw_curve.h
    io_png.c       io_png.h
    LICENSE.txt
    README.txt
    test_extract.cpp
    data/flower.png
    data/flowerOut.png
    third_party/...

## Limitations ##
The program extracts bilinear level lines at half-integer values (0.5, 1.5, etc). This can be modified (constant 'offset'), but it is important that the quantization avoids initial image levels (since iso-level sets can be quite complex at these levels).
