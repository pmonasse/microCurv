#include "mirror.h"
#include <algorithm>
#include <iostream>
#include <cstring>

unsigned char* extract_sub(unsigned char*in,size_t w,size_t h,unsigned char* bg,size_t bg_w,size_t bg_h, int X1,int Y1,int X2,int Y2,int Xc,int Yc)
{

size_t out_w,out_h;
unsigned char* out;
  /* test relative coordinates */
  if (X2<0) X2=w+X2-1;
  if (Y2<0) Y2=h+Y2-1;

  if (X2<X1 || Y2<Y1) {
      std::cerr<< "empty region to extract" << std::endl;
      return 0;
  }

  if (bg) {
      out_w=std::max((int)bg_w,Xc+X2-X1+1);
      out_h=std::max((int)bg_h,Yc+Y2-Y1+1);
    out =new unsigned char[out_w*out_h];
    std::fill(out,out+out_w*out_h,0);
    for (int x=0;x<bg_w;x++)
      for (int y=0;y<bg_h;y++)
	out[y*out_w+x] = bg[y*bg_w+x];
  } else {
      out_w=X2-X1+1;
      out_h=Y2-Y1+1;

    out =new unsigned char[out_h*out_w];
    std::fill(out,out+out_h*out_w,0);
  }


    for (int y=Y1;y<=Y2;y++)
      for (int x=X1;x<=X2;x++) {
     int pos1 = y*w+x;
     int pos2 = (Yc+y-Y1) * out_w+(Xc+x-X1);
      if ((Yc+y-Y1>=0) && (Yc+y-Y1<out_h) &&
	  (Xc+x-X1>=0) && (Xc+x-X1<out_w) &&
	  (x>=0) && (x<w) && (y>=0) && (y<h)) {
	out[pos2] = in[pos1];
      }
    }

  return out;
}

/* Symetrize an image*/
unsigned char* mirror(unsigned char* inImage,size_t w,size_t h, int m)
{
    int dir,n;
    unsigned char *in, *out;
    int j=0;

 /* Dimensions of the output image */
    size_t Ncol=w+2*m;
    size_t Nrow=h+2*m;

    unsigned char* bg = new unsigned char[Ncol*Nrow];
    std::fill(bg,bg+(Ncol*Nrow),0);
    unsigned char* outImage=extract_sub(inImage, w, h, bg, Ncol, Nrow, 0, 0, w, h, m, m);


    /* Left */
    for(int i = 0; i < h; i++) {
        in = outImage + Ncol*(i+m)+m;
        out = in-1;
        dir = +1;
        n = m;
        while(n != 0) {
            j = std::min(n,(int)w);
            n -= j;
            while(j--) {

                *out-- = *in;
                in  += dir;
            }

            dir = -dir;
            in += dir;
        }
    }

    /* Right */
    for(int i = 0; i < h; i++) {
        in = outImage + Ncol*(i+m)+m+w-1;
        out = in+1;
        dir = -1;
        n = m;
        while(n != 0) {
            j = std::min(n,(int)w);
            n -= j;
            while(j--) {
                *out++ = *in;
                in  += dir;
            }
            dir = -dir;
            in += dir;
        }
    }

    /* Top */
    n = h;
    dir = +1;
    in = outImage + Ncol*m;
    out = in - Ncol;
    for(int i = m-1; i >= 0; i--) {
        memcpy(out, in, Ncol*sizeof(unsigned char));
        in += dir*Ncol;
        out -= Ncol;
        if(--n == 0) {
            n = h;
            dir = -dir;
            in += dir*Ncol;
        }
    }

    /* Down */
    n = h;
    dir = -1;
    in = outImage + Ncol*(m+h-1);
    out = in + Ncol;
    for(int i = m-1; i >= 0; i--) {
        memcpy(out, in, Ncol*sizeof(unsigned char));
        in += dir*Ncol;
        out += Ncol;
        if(--n == 0) {
            n =h;
            dir = -dir;
            in += dir*Ncol;
        }
    }

    delete[] bg;
    return outImage;

}

/// Symmetrize by m pixels and paint border with constant gray
unsigned char* sym_enlarge(unsigned char* in, size_t w,size_t h, int m) {
    unsigned char* sym = mirror(in, w, h, m);
    w += 2*m;
    h += 2*m;
    for(int i=0; i<h; i++)
        sym[i*w]=sym[i*w+w-1]=0;
    std::fill(sym, sym+w, 0);
    std::fill(sym+(h-1)*w, sym+h*w, 0); 
    return sym;
}
