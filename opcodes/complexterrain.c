/*
    complexterrain.c:

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
#include <complex.h>
#define CA(z_,r_,i_) (z_) = (r_) + I * (i_)

/*  Wave-terrain synthesis opcode with complex functions
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
  MYFLT *cfunc;
  MYFLT *ktamp;
  MYFLT *ksize;
  MYFLT *kmode;

  double theta;

} COMPLEXTR;
#define FUNC(f_) double complex (*f_)(double complex, double, double)
static double complex poly4(double complex z, double ksize, double kamp) {
    return z*(z-(ksize+I*ksize))*(z-(-ksize+I*ksize))*(z-(ksize-I*ksize))*(z-(-ksize-I*ksize))/(cabs(z-(ksize+I*ksize))*cabs(z-(-ksize+I*ksize))*cabs(z-(ksize-I*ksize))*cabs(z-(-ksize-I*ksize))) + kamp;
}
static double complex poly2(double complex z, double ksize, double kamp) {
    return z*z/cabs(z) + kamp;
}

static double complex id(double complex z, double ksize, double kamp) {
    return kamp+sin(ksize*cabs(z))*z;
}

static double complex id2(double complex z, double ksize, double kamp) {
    return kamp+sin(ksize*creal(z))*cos(ksize*cimag(z))*z;
}

static double complex cos1(double complex z, double ksize, double kamp) {
    return kamp+ccos((1+ksize/100)*cabs(z)/z);
}

static double complex sin1(double complex z, double ksize, double kamp) {
    return kamp+csin((1+ksize/100)*cabs(z)/z);
}

static double complex fexp(double complex z, double ksize, double kamp ) {
    double complex ex;
    CA(ex,0,(-7.3+ksize)/(1/cpow(cabs(z),2)));
    return kamp + z * cexp(ex);
}
static double fradius(FUNC(f),double x, double y, double ksize, double kamp) {
    return cabs(f(x+I*y,ksize,kamp));
}
static double fphi(FUNC(f),double x, double y, double ksize, double kamp) {
    return carg(f(x+I*y,ksize,kamp))/M_PI;
}
static double fsinphi(FUNC(f),double x, double y, double ksize, double kamp) {
    return sin(carg(f(x+I*y,ksize,kamp)));
}
static double freal(FUNC(f),double x, double y, double ksize, double kamp) {
    return creal(f(x+I*y,ksize,kamp));
}
static double fimag(FUNC(f),double x, double y, double ksize, double kamp) {
    return cimag(f(x+I*y,ksize,kamp));
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
static double (*tfuncs[5])(double complex (*f)(double complex,double,double),double,double,double,double) = {fradius,fphi,fsinphi,freal,fimag};
static double complex (*cfuncs[7])(double complex,double,double) = {fexp,poly4,poly2,id,id2,cos1,sin1};
static int32_t wtinit(CSOUND *csound, COMPLEXTR *p)
{
    p->theta = 0.0;
    return OK;
}

static int32_t wtPerf(CSOUND *csound, COMPLEXTR *p)
{
    uint32_t offset = p->h.insdshead->ksmps_offset;
    uint32_t early  = p->h.insdshead->ksmps_no_end;
    uint32_t i, nsmps = CS_KSMPS;
    int32_t xloc, yloc;
    MYFLT xc, yc;
    MYFLT amp = *p->kamp;
    MYFLT pch = *p->kpch;
    MYFLT size = *p->ksize;
    MYFLT tamp = *p->ktamp;

    uint32_t kfunc = (uint32_t)*p->kfunc;
    if(kfunc>7) kfunc = 7; 
    if(kfunc<0) kfunc = 0; 
    uint32_t tfunc = (uint32_t)*p->kmode;
    if(tfunc<0) tfunc = 0;
    if(tfunc>4) tfunc = 4;
    uint32_t cfunc = (uint32_t)*p->cfunc;
    if(cfunc<0) cfunc = 0;
    if(cfunc>6) cfunc = 6;
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
      ifuncs[kfunc](theta,*p->kx,*p->ky,*p->krx,*p->kry,*p->kparam,&xc,&yc);
      rotate_point(*p->kx,*p->ky,*p->krot,&xc,&yc);
      /* OUTPUT AM */
      aout[i] = tfuncs[tfunc](cfuncs[cfunc],xc,yc,size,tamp) * amp;

      /* MOVE SCANNING POINT ROUND THE CURVE */
      theta += pch*((period*TWOPI_F) / csound->GetSr(csound));
    }
    
    p->theta = theta;
    return OK;
}

#define S(x)    sizeof(x)

static OENTRY localops[] = {
  { "complexterrain", S(COMPLEXTR), TR, 3,  "a", "kkkkkkkkkkkkk",
    (SUBR)wtinit, (SUBR)wtPerf },
};

LINKAGE
