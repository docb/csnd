/*
    lissterrain.c:

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
 *  - uses lissajous curves 
 *  
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
  MYFLT *ktabx, *ktaby;       /* Table numbers */
  MYFLT *lm1;
  MYFLT *lm2;
  MYFLT *la;
  MYFLT *lc;
  MYFLT *ln1;
  MYFLT *ln2;
  MYFLT *lb;
  MYFLT *speriod;


/* Internals */
  MYFLT oldfnx;  // storage of the current table for k-rate table change
  MYFLT oldfny;  // storage of the current table for k-rate table change

  MYFLT *xarr, *yarr;           /* Actual tables */

  MYFLT sizx, sizy;
  double theta;

} LISSTERR;

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

typedef struct lissparams {
  double m1;
  double m2;
  double n1;
  double n2;
  double a;
  double b;
  double c;
} LISSPARAMS;

static void lissformula(MYFLT t, MYFLT kx, MYFLT ky, MYFLT krx, MYFLT kry, LISSPARAMS *sp, MYFLT *outX, MYFLT *outY) {
    if(sp->m2 == 0) return;
    if((sp->c!=0) && (sp->n2!=0)) { 
      *outX = kx + krx*SIN(t);
      *outY = ky + kry*(SIN(t*sp->m1/sp->m2 + sp->a) + sp->c*SIN(t*sp->n1/sp->n2 + sp->b));
    } else {
      *outX = kx + krx*SIN(t);
      *outY = ky + kry*(SIN(t*sp->m1/sp->m2 + sp->a));
    }
}

static int32_t wtinit(CSOUND *csound, LISSTERR *p)
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

static int32_t wtPerf(CSOUND *csound, LISSTERR *p)
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

    LISSPARAMS s;
    s.m1 = FLOOR(*p->lm1);
    s.m2 = FLOOR(*p->lm2);
    s.n1 = FLOOR(*p->ln1);
    s.n2 = FLOOR(*p->ln2);
    s.a = *p->la;
    s.b = *p->lb;
    s.c = *p->lc;
    MYFLT period = 1;
    if(*p->speriod != 0) period = 1/(*p->speriod);

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
      lissformula(theta,*p->kx,*p->ky,*p->krx,*p->kry,&s,&xc,&yc);
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
  { "lterrain", S(LISSTERR), TR, 3,  "a", "kkkkkkkkkkkkkkkkk",
    (SUBR)wtinit, (SUBR)wtPerf },
};

LINKAGE
