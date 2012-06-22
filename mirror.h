#ifndef MIRROR_H
#define MIRROR_H

#include <cstddef>

unsigned char* extract_sub(unsigned char*in,size_t w,size_t h,
                           unsigned char* bg,size_t bg_w,size_t bg_h,
                           int X1,int Y1,int X2,int Y2,int Xc,int Yc);
unsigned char* mirror(unsigned char* inImage,size_t w,size_t h, int m);
unsigned char* sym_enlarge(unsigned char* in, size_t w,size_t h, int m);

#endif
