/*
    geneticterrain.c:

    (C) 2020 Christian Bacher

    This file is part of Csound.

    The Csound Library is free software; you can redistribute it
    and/or modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    Csound is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with Csound; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
    02110-1301 USA
*/

#include <csdl.h>
#include <math.h>

/*  Wave-terrain synthesis opcode with genetic functions
 *
 */
typedef struct {

  OPDS h;

  MYFLT *aout;
  MYFLT *kamp;
  MYFLT *kpch;
  MYFLT *kx;
  MYFLT *ky;
  MYFLT *krx;
  MYFLT *kry;
  MYFLT *krot; // rotation of the curve
  MYFLT *kfunc; // the curve index
  MYFLT *kparam;
  MYFLT *ksize;
  MYFLT *kgen[VARGMAX-1];
  double theta;
  int len;
  AUXCH genParam;

} GENETICTR;
static double sqr(double v) {
  if(v<0) return -sqrt(-v);
  return sqrt(v);
}

static double n(double v) {
   return v==0?0.000001:v;
}

  static double f00(double col, double x, double y) { return col+sin(x/n(y))*cos(x/n(y));} //D
  static double f01(double col, double x, double y) { return col+cos(x/n(y));} //D
  static double f02(double col, double x, double y) { return col+sin(y/n(x));} //D
  static double f03(double col, double x, double y) { return col+sin(x*x/n(y)-y*y/n(x));} //D
  static double f04(double col, double x, double y) { return col+cos(x*x/n(y))+sin(y*y/n(x));} //D
  static double f05(double col, double x, double y) { return col*sin(x*y*x)+cos(y*x*y);} //D
  static double f06(double col, double x, double y) { return col+((x+y*x)+sin(x*y)+cos(y/x));} //D
  static double f07(double col, double x, double y) { return col+sqr(sin(sqr(x)/n(sqr(y))));} //D

  static double f08(double col,double x, double y) { return col+sin(x*x+y*y);} //C
  static double f09(double col,double x, double y) { return col+sin(x*x)*cos(y*y);} //C
  static double f10(double col, double x, double y) { return col+fabs(x)*fabs(y);} //C
  static double f11(double col, double x, double y) { return col+sin(x*y)*cos(x*y);} //C
  static double f12(double col, double x, double y) { return col+sin(x*x-y*y);}  //C
  static double f13(double col, double x, double y) { return col+sin(cos(x)*fabs(y)*fabs(y));} //C
  static double f14(double col, double x, double y) { return col+sin(x*x*x-y*y*y);} //C
  static double f15(double col, double x, double y) { return col+sin(y*y*y)+sin(x*x*x);} //C
  static double f16(double col, double x, double y) { return col+cos(y*y*y+x*x*x);} //C
  static double f17(double col, double x, double y) { return col+cos(y*y*y)+cos(x*x*x);} //C
  static double f18(double col, double x, double y) { return col-tan(cos(sqrt(x*y*x*y)));} //C
  static double f19(double col, double x, double y) { return col+sin(x*x);} //C
  static double f20(double col, double x, double y) { return col+sin(x+y*x*y+x*x);} //C
  static double f21(double col, double x, double y) { return col+sin(y+x*y*x+y*y);} //C
  static double f22(double col, double x, double y) { return col+fabs(x*y+x*x+y*y);} //C
  static double f23(double col, double x, double y) { return col+((x+y)*y*x*sin(x)*cos(y));} //C
  static double f24(double col, double x, double y) { return col+sin(fabs(cos(x+y))+fabs(cos(y*x*y)));} //C

  static double f25(double col, double x, double y) { return col+sin(x)*cos(y);} //B
  static double f26(double col, double x, double y) { return col+cos(x)*sin(y)*cos(x*y);} //B
  static double f27(double col, double x, double y) { return col+sin(x)+sin(x)+cos(y)+cos(y);} //B
  static double f28(double col, double x, double y) { return col+cos(x)+cos(x)+sin(y)+sin(y);} //B
  static double f29(double col, double x, double y) { return col+sin(x)+cos(x)+sin(y)+cos(y);} //B
  static double f30(double col, double x, double y) { return col*cos(y)+sin(y)+cos(x)+sin(x);} //B
  static double f31(double col, double x, double y) { return col*sqrt(fabs(x)+fabs(y));} //B
  static double f32(double col, double x, double y) { return col+sqr(cos(x)+sqr(x)*sin(y)+sqr(y));} //B
  static double f33(double col, double x, double y) { return col*sin(col)*cos(x)*sin(x*y);} //B
  static double f34(double col, double x, double y) { return col*sin(col)*cos(y)*sin(x*y);} //B
  static double f35(double col, double x, double y) { return col+sin(x*y+x)+cos(y*x+y);} //B
  static double f36(double col, double x, double y) { return col+cos(sqr(x+y))*y+sqr(cos(y)*sin(x));} //B
  static double f37(double col, double x, double y) { return col+sin(sqr(y+x))*x+sqr(sin(x)*cos(y));} //B
  static double f38(double col, double x, double y) { return col+cos(x)*sin(x)+cos(y)*sin(y);} //B
  static double f39(double col, double x, double y) { return col+sin(fabs(cos(x))+fabs(sin(y)));} //B
  static double f40(double col, double x, double y) { return col+sin(sqrt(fabs(x)))-cos(sqrt(fabs(y)));} //B
  
  static double f41(double col, double x, double y) { return col+fabs(y)-x;} //A
  static double f42(double col, double x, double y) { return col+x+fabs(y);} //A
  static double f43(double col, double x, double y) { return col+fabs(x);} //A
  static double f44(double col, double x, double y) { return col+fabs(y);} //A
  static double f45(double col, double x, double y) { return col+y-fabs(x);} //A
  static double f46(double col, double x, double y) { return col+fabs(x)+y;} //A
  static double f47(double col, double x, double y) { return col+fabs(y*3);} //A
  static double f48(double col, double x, double y) { return col+fabs(x*3);} //A
  static double f49(double col, double x, double y) { return col*y-sin(x);} //A
  static double f50(double col, double x, double y) { return col*x-cos(y);} //A
  static double f51(double col, double x, double y) { return col*cos(x+y)*sin(x+y)/2;} //A
  static double f52(double col, double x, double y) { return atan(((y)+tan((x+y)-sin(x+M_PI)-sin(x*y/M_PI)*sin(((y*x+M_PI))))));} //A


static double (*tfuncs[53])(double,double,double) = {f00,f01,f02,f03,f04,f05,f06,f07,f08,f09,f10,f11,f12,f13,f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25,f26,f27,f28,f29,f30,f31,f32,f33,f34,f35,f36,f37,f38,f39,f40,f41,f42,f43,f44,f45,f46,f47,f48,f49,f50,f51,f52};
static double genomFunc(double *gen, int len, double x, double y) {
   double v = 1;
   for(int k=0;k<len;k++) {
     int index = -1;
     if(gen[k]>=0 && gen[k]<53) {
        index = (int)floor(gen[k]);
     } 
     if(index >=0)
       v = tfuncs[index](v,x,y);
   }
   return v;
}


static void rotate_point(MYFLT  cx, MYFLT  cy, MYFLT  angle, MYFLT *x, MYFLT *y)
{
  if(angle == 0) return;
  MYFLT s = SIN(angle);
  MYFLT c = COS(angle);

  *x -= cx;
  *y -= cy;

  float xnew = *x * c - *y * s;
  float ynew = *x * s + *y * c;

  *x = xnew + cx;
  *y = ynew + cy;
}

/* the normal eclipse function with center kx,ky and radius krx and kry */

static void ellipse(MYFLT t, MYFLT kx, MYFLT ky, MYFLT krx, MYFLT kry, MYFLT kparam, MYFLT *outX, MYFLT *outY ) {
    double x = t+kparam*SIN(t);
    *outX = kx + krx * SIN(x);
    *outY = ky + kry * COS(x);
}

/* the limacon curve parametrized by the kparam value 
   see e.g. http://www.2dcurves.com/roulette/roulettel.html
   for kparam = 1 we have a cardioid
*/

static void limacon(MYFLT t, MYFLT kx, MYFLT ky, MYFLT krx, MYFLT kry, MYFLT kparam, MYFLT *outX, MYFLT *outY ) {
    *outX = kx + krx * SIN(t) * (COS(t) + kparam); 
    *outY = ky + kry * COS(t) * (COS(t) + kparam);
}

/* a simple 8 */

static void lemniskateG(MYFLT t, MYFLT kx, MYFLT ky, MYFLT krx, MYFLT kry, MYFLT kparam, MYFLT *outX, MYFLT *outY ) {
    double x = t+kparam*SIN(t);
    *outX = kx + krx * COS(x);
    *outY = ky + kry * SIN(x)*COS(x);
}

/* the cornoid curve
   see e.g. http://www.2dcurves.com/sextic/sexticco.html
*/
static void cornoid(MYFLT t, MYFLT kx, MYFLT ky, MYFLT krx, MYFLT kry, MYFLT kparam, MYFLT *outX, MYFLT *outY ) {
    *outX = kx + krx * COS(t) * COS(2*t);
    *outY = ky + kry * SIN(t) * (kparam + COS(2*t));
}

/* Chevas trisextix
   see e.g. http://www.2dcurves.com/sextic/sextict.html
*/
static void trisec(MYFLT t, MYFLT kx, MYFLT ky, MYFLT krx, MYFLT kry, MYFLT kparam, MYFLT *outX, MYFLT *outY ) {
    *outX = kx + krx * COS(t) * (1+kparam*SIN(2*t));
    *outY = ky + kry * SIN(t) * (1+kparam*SIN(2*t));
}

/* Scarabeus curve see e.g http://www.2dcurves.com/sextic/sexticsc.html
*/

static void scarabeus(MYFLT t, MYFLT kx, MYFLT ky, MYFLT krx, MYFLT kry, MYFLT kparam, MYFLT *outX, MYFLT *outY ) {
    *outX = kx + krx * COS(t) * (kparam*SIN(2*t)+SIN(t));
    *outY = ky + kry * SIN(t) * (kparam*SIN(2*t)+SIN(t));
}
/* folium see http://www.2dcurves.com/quartic/quarticfo.html */
static void folium(MYFLT t, MYFLT kx, MYFLT ky, MYFLT krx, MYFLT kry, MYFLT kparam, MYFLT *outX, MYFLT *outY ) {
    double sint = SIN(t);
    double cost = COS(t);
    *outX = kx + krx * cost * cost * (sint*sint - kparam);
    *outY = ky + kry * sint * cost * (sint*sint - kparam);
}

/* talbot see http://www.2dcurves.com/trig/trigta.html */
static void talbot(MYFLT t, MYFLT kx, MYFLT ky, MYFLT krx, MYFLT kry, MYFLT kparam, MYFLT *outX, MYFLT *outY ) {
    double sint = SIN(t);
    double cost = COS(t);
    *outX = kx + krx * cost * (1 + kparam * sint*sint);
    *outY = ky + kry * sint * (1 - kparam - kparam*cost*cost);
}

static void (*ifuncs[8])(MYFLT,MYFLT,MYFLT,MYFLT,MYFLT,MYFLT,MYFLT*,MYFLT*) = { ellipse, lemniskateG, limacon, cornoid, trisec, scarabeus, folium, talbot }; 

static int32_t wtinit(CSOUND *csound, GENETICTR *p)
{
    p->theta = 0.0;
    p->len = csound->GetInputArgCnt(p) - 11;
    csound->AuxAlloc(csound, (p->len + 1)*sizeof(MYFLT), &(p->genParam));
    return OK;
}

static int32_t wtPerf(CSOUND *csound, GENETICTR *p)
{
    uint32_t offset = p->h.insdshead->ksmps_offset;
    uint32_t early  = p->h.insdshead->ksmps_no_end;
    uint32_t i, nsmps = CS_KSMPS;
    int32_t xloc, yloc;
    MYFLT xc, yc;
    MYFLT amp = *p->kamp;
    MYFLT pch = *p->kpch;
    MYFLT size = *p->ksize;

    uint32_t kfunc = (uint32_t)*p->kfunc;
    if(kfunc>7) kfunc = 7; 
    if(kfunc<0) kfunc = 0;

    MYFLT **genP = p->kgen;
    MYFLT *genParam = (MYFLT*)p->genParam.auxp; 
    for(int k=0;k<p->len;k++) {
      genParam[k] = genP[k][0];
    }

    MYFLT period = 1;
    MYFLT theta = p->theta;
    MYFLT *aout = p->aout;
    if (UNLIKELY(offset)) memset(aout, '\0', offset*sizeof(MYFLT));
    if (UNLIKELY(early)) {
      nsmps -= early;
      memset(&aout[nsmps], '\0', early*sizeof(MYFLT));
    }
    for (i=offset; i<nsmps; i++) {

      /* COMPUTE LOCATION OF SCANNING POINT */
      ifuncs[kfunc](theta,*p->kx*size,*p->ky*size,*p->krx,*p->kry,*p->kparam,&xc,&yc);
      rotate_point(*p->kx*size,*p->ky*size,*p->krot,&xc,&yc);
      /* OUTPUT AM */
      aout[i] = genomFunc(genParam,p->len,xc,yc) * amp;

      /* MOVE SCANNING POINT ROUND THE CURVE */
      theta += pch*((period*TWOPI_F) / csound->GetSr(csound));
    }
    
    p->theta = theta;
    return OK;
}

#define S(x)    sizeof(x)

static OENTRY localops[] = {
  { "geneticterrain", S(GENETICTR), TR, 3,  "a", "kkkkkkkkkkz",
    (SUBR)wtinit, (SUBR)wtPerf },
};

LINKAGE
