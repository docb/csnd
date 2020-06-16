/*
    wterr.c:

    Copyright (C) 2002 Matt Gilliard, John ffitch
    for the original file wave-terrain.c from the csound distribution

    Modifications and enhancements by (C) 2020 Christian Bacher

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

/*  Wave-terrain synthesis opcode
 *
 *  author: m gilliard
 *          en6mjg@bath.ac.uk
 *
 *  enhancements and modifications
 *  Christian Bacher docb22@googlemail.com
 *  Changes to the original:
 *  - Added curves: limacon with parameter, lemniskate (G), lissajous (4 variants), 
 *  -               rhodonea (5 variants), cornoid with parameter, trisec (Ceva) with parameter
 *  - tables are krate
 *  - added k parameter for rotating the curve arround the current x,y
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
  MYFLT *kparam; // curve parameter

  MYFLT *ktabx, *ktaby;       /* Table numbers */

/* Internals */
  MYFLT oldfnx;  // storage of the current table for k-rate table change
  MYFLT oldfny;  // storage of the current table for k-rate table change

  MYFLT *xarr, *yarr;           /* Actual tables */

  MYFLT sizx, sizy;
  double theta;

} WAVETER;

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
    *outX = kx + krx * SIN(t);
    *outY = ky + kry * COS(t);
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
    *outX = kx + krx * COS(t);
    *outY = ky + kry * SIN(t)*COS(t);
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

/*
    some lissajous curves, see parameters in LISSPARAMS
*/

typedef struct lissparams {
  double n;
  double m;
  double phi;
} LISSPARAMS;

static LISSPARAMS lp[] = { {1,2,M_PI/2}, {3,2,0}, {3,2,M_PI/2}, {3,2,M_PI/4},{3,4,M_PI/4}};


static void lissajous(MYFLT t,  MYFLT kx, MYFLT ky, MYFLT krx, MYFLT kry, MYFLT kparam, MYFLT *outX, MYFLT *outY) {
    int index = (int)kparam;
    if(index > 4) index = 4;
    *outX = kx + krx * SIN(t);
    *outY = ky + kry * SIN((lp[index].n/lp[index].m)*t + lp[index].phi);
}

/*
    some rhodonea curves, see parameters in RHODOPARAMS
*/
typedef struct rhodoparams {
  double n;
  double m;
} RHODOPARAMS;

static RHODOPARAMS rp[] = { {1,2}, {2,1}, {3,1}, {3,2}, {4,1}};

static void rhodonea(MYFLT t,  MYFLT kx, MYFLT ky, MYFLT krx, MYFLT kry, MYFLT kparam , MYFLT *outX, MYFLT *outY) {
    int index = (int)kparam;
    if(index > 4) index = 4;
    *outX = kx + krx*SIN((rp[index].n/rp[index].m)*t)*COS(t);
    *outY = ky + kry*SIN((rp[index].n/rp[index].m)*t)*SIN(t);
}



static MYFLT lissp(MYFLT kp) {
    int index = (int)kp;
    if(index > 4) index = 4;
    return (MYFLT)lp[index].m;
}

static MYFLT rhodop(MYFLT kp) {
    int index = (int)kp;
    if(index > 4) index = 4;
    return (MYFLT)rp[index].m;
}

static MYFLT one(MYFLT kp) {
   return 1;
}


static void (*ifuncs[7])(MYFLT,MYFLT,MYFLT,MYFLT,MYFLT,MYFLT,MYFLT*,MYFLT*) = { ellipse, limacon, lemniskateG, lissajous, rhodonea, cornoid, trisec }; 
static MYFLT (*pfuncs[7])(MYFLT) = { one, one, one, lissp, rhodop, one,one }; 

static int32_t wtinit(CSOUND *csound, WAVETER *p)
{
    p->xarr = NULL;
    p->yarr = NULL;

    p->oldfnx = -1;
    p->oldfny = -1;
    p->sizx = 0;
    p->sizy = 0;
    p->theta = 0.0;
    return OK;
}

static int32_t wtPerf(CSOUND *csound, WAVETER *p)
{
    uint32_t offset = p->h.insdshead->ksmps_offset;
    uint32_t early  = p->h.insdshead->ksmps_no_end;
    uint32_t i, nsmps = CS_KSMPS;
    int32_t xloc, yloc;
    MYFLT xc, yc;
    MYFLT amp = *p->kamp;
    MYFLT pch = *p->kpch;

    if (*(p->ktabx) != p->oldfnx || p->xarr == NULL) {
      p->oldfnx = *(p->ktabx);
      FUNC *ftp = csound->FTFindP(csound, p->ktabx);    /* new table parameters */
      if (UNLIKELY((ftp == NULL) || ((p->xarr = ftp->ftable) == NULL))) return NOTOK;
      p->sizx = (MYFLT)ftp->flen;
    }
    if (*(p->ktaby) != p->oldfny || p->yarr == NULL) {
      p->oldfny = *(p->ktaby);
      FUNC *ftp = csound->FTFindP(csound, p->ktaby);    /* new table parameters */
      if (UNLIKELY((ftp == NULL) || ((p->yarr = ftp->ftable) == NULL))) return NOTOK;
      p->sizy = (MYFLT)ftp->flen;
    }


    uint32_t kfunc = (uint32_t)*p->kfunc;

    MYFLT period = pfuncs[kfunc](*p->kparam);

    MYFLT sizx = p->sizx, sizy = p->sizy;
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
      /* MAP SCANNING POINT TO BE IN UNIT SQUARE */
      xc = xc-FLOOR(xc);
      yc = yc-FLOOR(yc);

      /* SCALE TO TABLE-SIZE SQUARE */
      xloc = (int32_t)(xc * sizx);
      yloc = (int32_t)(yc * sizy);

      /* OUTPUT AM OF TABLE VALUES * KAMP */
      aout[i] = p->xarr[xloc] * p->yarr[yloc] * amp;

      /* MOVE SCANNING POINT ROUND THE ELLIPSE */
      theta += pch*((period*TWOPI_F) / csound->GetSr(csound));
    }
    
    p->theta = theta;
    return OK;
}

#define S(x)    sizeof(x)

static OENTRY localops[] = {
  { "wterr", S(WAVETER), TR, 3,  "a", "kkkkkkkkkkk",
    (SUBR)wtinit, (SUBR)wtPerf },
};

LINKAGE
