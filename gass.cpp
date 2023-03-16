// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file gass.cpp
 * @brief Geometric Affine Scale Space
 * @author Lionel Moisan <Lionel.Moisan@parisdescartes.fr>
 *         Pascal Monasse <monasse@imagine.enpc.fr>
 * 
 * Copyright (c) 2002, 2012-2017, Lionel Moisan, Pascal Monasse
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
#include <cmath>
#include <cstddef>
#include <cassert>

#define ABS(x)       ( (x)>0?(x):-(x) )
#define SGN(x)       (((x)==0.0)? 0: (((x)>0.0)?1: -1))

inline double det3(const DPoint& a, const DPoint& b, const DPoint& c) {
    return ((b.x-a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x));
}

inline double area3(const DPoint& a, const DPoint& b, const DPoint& c) {
    double tmp = det3(a,b,c)/2.0;
    return ABS(tmp);
}

#define EPSILON  1e-15          // relative precision for a double

/// return +1, 0 or -1, the sign of det(b-a,c-b) modulo double precision
static int dir(const DPoint& a, const DPoint& b, const DPoint& c)
{
    double det = det3(a,b,c);
    double prec = EPSILON * (ABS(b.x - a.x) * (ABS(c.y) + ABS(b.y)) +
                             ABS(c.y - b.y) * (ABS(b.x) + ABS(a.x)) +
                             ABS(b.y - a.y) * (ABS(c.x) + ABS(b.x)) +
                             ABS(c.x - b.x) * (ABS(b.y) + ABS(a.y)));
    if(ABS(det) <= prec)
        det = 0.0;
    return SGN(det);
}

/*----------------- Split a curve into convex components -----------------*/

static void split_convex(const std::vector<DPoint>& in,
                         std::vector<DPoint>& out, std::vector<size_t>& cvx)
{
    typedef std::vector<DPoint>::const_iterator DPointPtr;
    // initialization
    DPointPtr p=in.begin(), first=p, pmax=in.end();
    bool closed = (in.front()==in.back());

    if(in.size()<4) {
        out = in;
        cvx.push_back( out.size() );
        return;
    }

    int ni = 0; // Number of inflection points
    DPointPtr p1(p++),p2(p1),p3(p++);
    out.push_back(*p1);

    // Find first angle at p2
    int d1=0;
    while(d1==0 && p!=pmax) {
        p2 = p3;
        p3 = p++;
        d1 = dir(*p1, *p2, *p3);
    }

    while(p!=pmax) {
        DPointPtr p4=p++;
        int d2 = dir(*p2, *p3, *p4);
        if(d2!=0) {
            out.push_back(*p2);
            if(d1==d2) {
                if(p == first) { // Rare, but can happen
                    assert(closed);
                    out.push_back(out.front());
                    cvx.push_back(out.size());
                    break;
                }
            } else {
                d1 = d2; // For next iteration
                DPoint m = .5 * (*p2 + *p3);
                out.push_back(m);
                cvx.push_back( out.size() );
                if(p == first) // End of loop
                    break;
                // First inflection point of closed curve? Start over here
                else if(closed && ni==0) {
                    out.clear();
                    cvx.clear();
                    cvx.push_back(0);
                    first = p;
                }
                // Begin new convex part
                out.push_back(m);
                ++ ni;
            }
            p1=p2; p2=p3;
        }
        p3=p4; // Move triangle p1-p2-p3 to p2-p3-p4
        
        if(p==pmax && closed) // Loop for closed curve
            p = in.begin()+1;
        if(ni==0 && p==in.begin()+3) { // Convex closed curve
            out.push_back( out.front() );
            cvx.push_back( out.size() );
            break;
        }
    }

    if (! closed) {
        out.push_back(*p2);
        out.push_back(*p3);
        cvx.push_back( out.size() );
    }
}

/*------------------------------- SAMPLING  -------------------------------*/

/// sample a curve
static void sample(DPoint *in, int size, std::vector<DPoint>& out, double eps2)
{
    // return if the curve has less than 3 points
    if(size<3) {
        for(int i=0; i<size; i++)
            out.push_back(*in++);
        return;
    }

    DPoint* pmax = in+size;
    size = out.size(); // keep start of out vector for second pass
    DPoint p=*in++;
    out.push_back(p);

    //--- first pass: insert points
    while(in!=pmax) {
        DPoint op = p;
        p = *in++;
        DPoint d = p-op;
        double d2 = d.x*d.x + d.y*d.y;
        if(d2>=2.*eps2) { // insert n points
            int n = (int)floor( sqrt(d2/eps2) );
            d /= (double)n;
            for(int k=1; k<n; k++)
                out.push_back(op + (double)k*d);
        }
        out.push_back(p);
    }

    //--- second pass: remove points
    pmax = &out.back();
    std::vector<DPoint>::iterator it=out.begin()+size;
    in = &(*it++);
    DPoint op = *in++;

    while(in!=pmax) {
        p = *in++;
        DPoint d = p-op;
        if(d.x*d.x + d.y*d.y >= eps2)
            *it++ = op = p;
    }
    *it++ = *in;
    out.erase(it, out.end());
}

/// signed area of a polygonal sector p-q1-q2-p
static double area_pol(DPoint& p, DPoint* q1, DPoint* q2) {
    double area = 0.;
    for(DPoint* q=q1; q!=q2; ++q)
        area += det3(q[0], q[1], p);
    return (area/2.0);
}

/*------------------------- AFFINE CONVEX EROSION -------------------------*/
static void aceros(DPoint* in, int size,
                   std::vector<DPoint>& out, double area_sz) {
    // deal with singular cases (less than 2 points)
    if(size<2)
        return;

    // test if the curve is closed
    DPoint* pmax = in+size;
    bool closed = (*in==pmax[-1]);
    if(closed)
        --pmax;

    // deal with singular cases (2 or 3 points, closed)
    if(size<4 && closed)
        return;

    // return input if segment or area_sz=0
    if(size==2 || area_sz==0.) {
        while(in!=pmax)
            out.push_back(*in++);
        return;
    }

    // compute total area
    double tot_area = area_pol(*in, in+1, pmax-1);
    tot_area = ABS(tot_area);

    // check extinction
    if (closed) {
        if(area_sz>=tot_area/2.1) // theoretically: 2.0
            return;
    } else if(area_sz>=tot_area) {
        out.push_back(*in);
        out.push_back(pmax[-1]);
        return;
    }

    if(!closed)
        out.push_back(*in);
    DPoint *p=in, *p0=p, *q0=in+1, *p1=in+1, *q1=in+2;
    double cur_area=0.0;
    bool okp=false, okq=false;

    do { // MAIN LOOP: compute the middle points of significative chords
        if(cur_area<=area_sz) {
            double inc_area = area3(*q0, *q1, *p0);

            if(cur_area+inc_area>=area_sz) { // compute middle point
                double lambda = double((area_sz-cur_area)/inc_area);
                out.push_back(.5*(*p0+(1-lambda)**q0+lambda**q1));
            }
            if(cur_area+inc_area-area3(*p0,*p1,*q1)>area_sz) {
                cur_area -= area3(*p0, *p1, *q0);
                p0 = p1++;
                if(closed && p1==pmax)
                    p1 -= size-1;
                if(p0==p)
                    okp = true;
            }
            else {
                cur_area += inc_area;
                q0 = q1++;
                if(closed && q1==pmax)
                    q1 -= size-1;
                if(q0==p+1)
                    okq = true;
            }
        } else {
            double inc_area = area3(*p0, *p1, *q0);
            if(cur_area-inc_area<=area_sz) { // compute middle point
                double lambda = double((cur_area-area_sz)/inc_area);
                out.push_back(.5*(*q0+(1-lambda)**p0+lambda**p1));
            }
            if(!closed && q1!=pmax &&
               cur_area-inc_area+area3(*p1, *q0, *q1)<area_sz) {
                cur_area += area3(*p0, *q0, *q1);
                q0 = q1++;
                if(closed && q1==pmax)
                    q1 -= size-1;
                if(q0==p+1)
                    okq = true;
            } else {
                cur_area -= inc_area;
                p0 = p1++;
                if(closed && p1==pmax)
                    p1 -= size-1;
                if(p0==p)
                    okp = true;
            }
        }

        // more precise computation of cur_area if needed
        if(p1==q0)
            cur_area = 0.0;
        DPoint* p2 = p1+1;
        if(closed && p2==pmax)
            p2 -= size-1;
        if(p2==q0)
            cur_area = area3(*p0, *p1, *q0);
    } while(!(closed? (okp&&okq): ((q0+1==pmax&&cur_area<=area_sz))));

    // add last point to output
    out.push_back(closed? out.front(): pmax[-1]);
}

/*----------------------- DISCRETE AFFINE EROSION  -----------------------*/

/// input/output curve
/// desired absolute area step (real one is returned)
/// absolute precision squared
static void dafferos(std::vector<DPoint>& l, double& area_sz, double eps2,
                     std::vector<DPoint>& pts1, std::vector<DPoint>& pts2) {
    double min_area = area_sz/8.0;// critical area for effective erosion

    // compute convex components
    pts1.clear();
    pts2.clear();
    std::vector<size_t> cvx;
    cvx.push_back(0);
    split_convex(l, pts1, cvx);

    // compute minimal area
    double narea = area_sz;
    std::vector<size_t> pos2;
    pos2.push_back(0);
    for(std::vector<size_t>::iterator it=cvx.begin(),itn=it+1;
        itn!=cvx.end(); it=itn++) { // resample curve
        size_t front=pts2.size();
        sample(&pts1[*it], *itn-*it, pts2, eps2);
        pos2.push_back(pts2.size());
        DPoint* first = &pts2[front];
        DPoint* last  = &pts2.back();
        double a = (*first==*last)? narea: fabs(area_pol(*first,first,last));
        if(min_area<=a && a<narea)
            narea = a;
    }
    if(area_sz>narea)
        area_sz = narea;

    // apply aceros to convex components and link result
    l.clear();
    for(std::vector<size_t>::iterator it=pos2.begin(),itn=it+1;
        itn!=pos2.end(); it=itn++)
        aceros(&pts2[*it], *itn-*it, l, narea);
}

/*------------------------------ MAIN MODULE  ------------------------------*/

void gass(std::vector<DPoint>& curve,
          double first, double last, double prec, double maxStep) {
    double eps2 = 1.0/(prec*prec); // absolute precision (squared)
    static const double omega = 0.5*pow(1.5,2.0/3.0); // normalization constant
    double step_a = maxStep*maxStep*pow(0.75/omega,1.5); // scale (area) step

    // remaining scale (additive normalization)
    double remaining_h = (pow(last,4.0/3.0)-pow(first, 4.0/3.0)) * 0.75/omega;

    std::vector<DPoint> pts1, pts2;

    while(!curve.empty() && remaining_h>0.0) {
        double remaining_a = pow(remaining_h, 1.5);
        double a = std::min(remaining_a,step_a);
        dafferos(curve, a, eps2, pts1, pts2);
        remaining_h -= pow(a, 2.0 / 3.0);
    }

    if (!curve.empty()) {
        pts1.clear();
        for(std::vector<DPoint>::iterator it=curve.begin();it!=curve.end();++it)
            pts1.push_back(*it);
        curve.clear();
        sample(&(*pts1.begin()), pts1.size(), curve, eps2);
    }
}
