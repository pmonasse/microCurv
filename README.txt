# Extraction of the Level Lines of a Bilinear Image #

## Summary ##
This software computes the mean curvature map of an image. It is linked to an IPOL article [1]. This is based on level line tree extraction of the bilinear interpolation of a digital image [2] and affine invariant smoothing. These functionalities were inspired from the equivalent MegaWave2 [3] functions 'flst_bilinear' and 'gass', though the tree extraction is based on a much simpler algorithm than 'flst_bilinear'.

[1] The Image Curvature Microscope, IPOL, https://doi.org/10.5201/ipol.2017.212
[2] Level Lines of a Bilinear Image, IPOL, https://doi.org/10.5201/ipol.2019.269
[3] MegaWave2, http://megawave.cmla.ens-cachan.fr/


## Authors ##
Pascal Monasse <pascal.monasse@enpc.fr>
Adina Ciomaga <adina@math.univ-paris-diderot.fr>
Lionel Moisan <Lionel.Moisan@parisdescartes.fr>

Laboratoire d'Informatique Gaspard Monge (LIGM)/
Ecole des Ponts ParisTech

## Version ##
Version 2.0, released on 2023/03/21

Old versions:
[1] microCurv Version 1.0 released on 2017/06/29
[2] bilines Version 1.0 released on 2019/06/26

Future releases and updates:
<https://github.com/pmonasse/microCurv.git>

## Build ##
Prerequisites: CMake version 2.8 or later

- Unix, MacOS:
  $ cd /path_to_this_file/
  $ mkdir Build && cd Build && cmake -DCMAKE_BUILD_TYPE:bool=Release ..
  $ cmake --build .
- Windows with MinGW:
  $ cd /path_to_this_file/
  $ mkdir Build 
  $ cd Build
  $ cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE:bool=Release ..
  $ cmake --build .

CMake tries to find libPNG and libTIFF on your system, otherwise it uses the embedded version. The library needs to come with header files, which are provided by "...-dev" packages under Linux Debian or Ubuntu.

- Test:
    $ ./bilines ../data/flower.png out.png
Compare out.png with ../data/flowerLL.png
    $ ./microCurv ../data/flower.png out.png
Compare out.png with ../data/flowerCurv.png

## Usage ##
bilines [options] in.png out.png
  - in.png: PNG input image
  - out.png: PNG output image
Option:
  -q <int>: Quantization step (default 32)

microCurv [options] in.png outCurv.png [outCurv.tif]
  - in.png: PNG input image
  - outCurv.png: PNG output image with superimposed curvatures
  - outCurv.tif: float TIFF output image of mean curvature map
Options:
  -q <int>: quantization level step of level lines (default 8)
  -s <float>: scale of smoothing by affine scale-space (default 2)
  -z <float>: zoom factor for output bitmap images (default 1)
  -o <fileName.png>: output image after level line smoothing
  -I <fileName>: output file for level lines before smoothing
  -O <fileName>: output file for level lines after smoothing
  -r <ROI>: region of interest in input image. Format: wxh+x0+y0
The extension of <fileName> in -I and -O options determines the file format: SVG (Scalable Vector Graphics, extension .svg) or PNG (any other extension). SVG output keeps the original size of the image. All other output images are bitmaps and scaled by the zoom factor. In particular, notice that the optional TIFF output image records the original image curvatures but in a zoomed bitmap image.

## Files (Only those with &ast; are reviewed through IPOL) ##

    bilines.cpp (*): main function, calls [2] Algo. 1 and implements Algo. 6.
    fill_curve.cpp (*) fill_curve.h (*): implements [2] Algo. 5.
    levelLine.cpp (*)  levelLine.h (*): implements [2] Algo. 4.
    lltree.cpp (*)     lltree.h (*): implements [2] Algos. 1-3.
    curv.cpp (*)       curv.h (*): implements [1] Algo. 1 14-16 and (10),(11).
    mainMicroCurv.cpp (*): main function, [1] Algo. 1
    microCurv.cpp (*)  microCurv.h (*)

    CMakeLists.txt
    cmdLine.h
    draw_curve.cpp draw_curve.h
    fill_curve.cpp fill_curve.h
    gass.cpp       gass.h
    image.cpp      image.h
    io_png.c       io_png.h
    io_tiff.c      io_tiff.h
    LICENSE.txt
    README.txt
    test_extract.cpp
    data/flower.png data/flowerCurv.png data/flowerLL.png
    third_party/...

## Limitations ##
The program extracts bilinear level lines at half-integer values (0.5, 1.5, etc). This can be modified (constant 'offset'), but it is important that the quantization avoids initial image levels (since iso-level sets can be quite complex at these levels).
