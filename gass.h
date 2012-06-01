#ifndef GASS_H
#define GASS_H

#include <vector>

struct DPoint {
    double x,y;
    DPoint(double x0, double y0): x(x0), y(y0) {}
    bool operator==(const DPoint& p) const { return (x==p.x&&y==p.y); }
    bool operator!=(const DPoint& p) const { return !operator==(p); }
    DPoint operator+(const DPoint& p) const { return DPoint(x+p.x,y+p.y); }
    DPoint operator-(const DPoint& p) const { return DPoint(x-p.x,y-p.y); }
    DPoint& operator/=(double d) { x/=d; y/=d; return *this; }
};

inline DPoint operator*(double d, const DPoint& p)
{ return DPoint(d*p.x,d*p.y); }

void gass(std::vector<DPoint>& curve,
          double first=0, double last=1, double prec=5.0, double maxStep=1.0);

#endif
