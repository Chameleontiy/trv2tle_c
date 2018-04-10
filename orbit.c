#include "orbit.h"
//#include "constants.h"
#include "math.h"

SatElem OE;
double rr[3], vv[3];

double acose(double x)
{
    double rval;

    if( x >= 1.)
        rval = 0.;
    else if( x <= -1.)
        rval = pi;
    else rval = acos( x);
    return( rval);
}
//   v1 . v2
double dot(double* v1, double* v2)
{
    double sum = 0;
    int i;
    for (i = 0; i < 3; i++)
        sum += v1[i] * v2[i];
    return sum;
}
//   ||v||
double norm(double* v)
{
    return sqrt(dot(v, v));
}
//  a * v = av
void smult(double a, double* v, double* av)
{
    int i;
    for (i = 0; i < 3; i++)
        av[i] = a * v[i];
}
//   v1 + v2 = s
void vadd(double* v1, double* v2, double* res)
{
    int i;
    for (i = 0; i < 3; i++)
        res[i] = v1[i] + v2[i];
}
//   s = v1 + a * v2   (used for subtraction)
void vmadd(double* v1, double* v2, double* res, double a)
{
    int i;
    for (i = 0; i < 3; i++)
        res[i] = v1[i] + a * v2[i];
}
//   v1 x v2 = b
void cross(double v1[3], double v2[3], double b[3])
{
    b[0] = v1[1] * v2[2] - v1[2] * v2[1];
    b[1] = v1[2] * v2[0] - v1[0] * v2[2];
    b[2] = v1[0] * v2[1] - v1[1] * v2[0];
}
//   u = v / ||v||
void unitv(double v[], double u[])
{
    double no = norm(v);
    int i;
    for (i = 0; i < 3; i++)
        u[i] = v[i] / no;
}
// unit vector projection of v1 onto v2
void proj(double v2[3], double v1[3], double vp[3])
{
    double b, temp[3];
    b = dot(v2, v1)/dot(v2, v2);
    smult(b, v2, temp);
    unitv(temp, vp);
}
// floating point mod function
double mod(double var, double div)
{
    double val = fmod(var, div);
    return((val < 0) ? div + val : val);
}
double fmod2ps( const double x)
{
    double rval = fmod( x, twopi);
    if(rval < 0.)
        rval += twopi;
    return(rval);
}
// round value to decimal place, p
//double round(double value, double p)
//{
//    double place = pow(10, p);
//	return ( floor(value * place + 0.5) / place );
//}

void sgp4(double tsince)
{
   // adapted from projectpluto.com

   // global tle variables
   // jd, xndt2o, xndd6o, bstar, c2;
   // xincl, xnodeo, eo, omegao, xmo, xno;

   double c1, c4, xnodcf, t2cof, /*aodp,*/ cosio, sinio,
          omgdot, xmdot, xnodot, xnodp, c5, d2, d3,
          d4, delmo, eta, omgcof, sinmo, t3cof, t4cof,
          t5cof, xmcof, eosq, betao, betao2;
   int simple_flag;
   double MINIMAL_E = 1.e-4;
   double ECC_EPS = 1.e-6;     /* Too low for computing further drops. */

   double
         coef, coef1, eeta, etasq, perige, pinv, pinvsq,
         psisq, qoms24, s4, temp1, temp2, temp3,
         theta4, tsi, tsi_squared, x1mth2, x1m5th, xhdot1;

   const double a1 = pow(xke/ OE.xno,two_thirds);
   double del1, ao, delo, x3thm1, tval, theta2;

   /* Recover original mean motion (xnodp) and   */
   /* semimajor axis (aodp) from input elements. */
   cosio = cos( OE.xincl);
   theta2 =  cosio* cosio;
   x3thm1 = 3. *  theta2 - 1.;
   eosq =  OE.eo* OE.eo;
   betao2 = 1- eosq;
   betao = sqrt( betao2);
   tval = 1.5 * ck2 * x3thm1 / ( betao *  betao2);
   del1 = tval / (a1 * a1);
   ao = a1*(1-del1*(0.5*two_thirds+del1*(1.+134./81.*del1)));
   delo = tval / (ao * ao);
   xnodp =  OE.xno/(1+delo);
   OE.aodp = ao/(1-delo);

   x3thm1 = 3* theta2-1;
   /* For perigee below 156 km, the values */
   /* of s and qoms2t are altered.         */
   s4 = s;
   qoms24 = qoms2t;
   perige = ( OE.aodp*(1- OE.eo)-1.)*xkmper;
   if(perige < 156)
   {
      double temp_val, temp_val_squared;

      if(perige <= 98)
          s4 = 20;
      else
          s4 = perige-78.;
      temp_val = (120. -  s4) / xkmper;
      temp_val_squared = temp_val * temp_val;
      qoms24 = temp_val_squared * temp_val_squared;
       s4 =  s4/xkmper+1.;
   }  /* End of if(perige <= 156) */

   pinv = 1. / ( OE.aodp* betao2);
   pinvsq = pinv * pinv;
   tsi = 1. / ( OE.aodp -  s4);
   eta =  OE.aodp* OE.eo* tsi;
   etasq =  eta* eta;
   eeta =  OE.eo* eta;
   psisq = fabs(1-etasq);
   tsi_squared =  tsi *  tsi;
   coef = qoms24 * tsi_squared * tsi_squared;
   coef1 =  coef / pow(psisq,3.5);
   OE.c2 =  coef1 *  xnodp * ( OE.aodp*(1+1.5*etasq+eeta*
   (4+etasq))+0.75*ck2* tsi/psisq*x3thm1*(8+3*etasq*(8+etasq)));
   c1 =  OE.bstar*OE.c2;
   sinio = sin( OE.xincl);
   x1mth2 = 1- theta2;
   c4 = 2* xnodp* coef1* OE.aodp* betao2*
        ( eta*(2+0.5*etasq)+ OE.eo*(0.5+2*etasq)-2*ck2* tsi/
        ( OE.aodp*psisq)*(-3*x3thm1*(1-2*eeta+etasq*
        (1.5-0.5*eeta))+0.75*x1mth2*(2*etasq-eeta*(1+etasq))*
        cos(2* OE.omegao)));
   theta4 =  theta2* theta2;
   temp1 = 3*ck2*pinvsq* xnodp;
   temp2 = temp1*ck2*pinvsq;
   temp3 = 1.25*ck4*pinvsq*pinvsq* xnodp;
   xmdot =  xnodp+0.5*temp1* betao*
               x3thm1+0.0625*temp2* betao*
                    (13-78* theta2+137*theta4);
   x1m5th = 1-5* theta2;
   omgdot = -0.5*temp1*x1m5th+0.0625*temp2*
                     (7-114* theta2+395*theta4)+
                temp3*(3-36* theta2+49*theta4);
   xhdot1 = -temp1* cosio;
   xnodot = xhdot1+(0.5*temp2*(4-19* theta2)+
            2*temp3*(3-7* theta2))* cosio;
   xnodcf = 3.5* betao2*xhdot1*c1;
   t2cof = 1.5*c1;

   eeta = OE.eo*eta;

   /* For perigee less than 220 kilometers, the "simple" flag is set */
   /* and the equations are truncated to linear variation in sqrt a  */
   /* and quadratic variation in mean anomaly.  Also, the c3 term,   */
   /* the delta omega term, and the delta m term are dropped.        */
   simple_flag = (OE.aodp*(1- OE.eo) < (220./xkmper+1.));
   if( !simple_flag)
      {
      const double c1sq = c1*c1;
      double temp;

      simple_flag = 0;
      delmo = 1. + eta * cos( OE.xmo);
      delmo *= delmo * delmo;
      d2 = 4*OE.aodp*tsi*c1sq;
      temp = d2* tsi*c1/3;
      d3 = (17*OE.aodp+ s4)*temp;
      d4 = 0.5*temp* OE.aodp* tsi*(221*OE.aodp+31* s4)*c1;
      t3cof = d2+2*c1sq;
      t4cof = 0.25*(3*d3+c1*(12*d2+10*c1sq));
      t5cof = 0.2*(3*d4+12*c1*d3+6*d2*d2+15*c1sq*(2*d2+c1sq));
      sinmo = sin( OE.xmo);
      if(  OE.eo < MINIMAL_E)
         omgcof = xmcof = 0.;
      else
      {
         const double c3 =
              coef * tsi * a3ovk2 * xnodp * sinio /  OE.eo;

         xmcof = -two_thirds * coef *  OE.bstar / eeta;
         omgcof =  OE.bstar*c3*cos( OE.omegao);
      }
   } /* End of if (isFlagClear(SIMPLE_FLAG)) */
   etasq =  eta *  eta;
   c5 = 2* coef1* OE.aodp * betao2*(1+2.75*(etasq+eeta)+eeta*etasq);

  double
        a, e, omega, omgadf,
        temp, tempa, tempe, templ, tsq,
        xl, xmdf, xmp, xnoddf, xnode;

  /* Update for secular gravity and atmospheric drag. */
  xmdf =  OE.xmo+ xmdot*tsince;
  omgadf =  OE.omegao+ omgdot*tsince;
  xnoddf =  OE.xnodeo+ xnodot*tsince;
  omega = omgadf;
  xmp = xmdf;
  tsq = tsince*tsince;
  xnode = xnoddf+xnodcf*tsq;
  tempa = 1-c1*tsince;
  tempe =  OE.bstar*c4*tsince;
  templ = t2cof*tsq;
  if( !simple_flag)
  {
      const double delomg = omgcof*tsince;
      double delm = 1. +  eta * cos(xmdf);
      double tcube, tfour;

      delm = xmcof * (delm * delm * delm - delmo);
      temp = delomg+delm;
      xmp = xmdf+temp;
      omega = omgadf-temp;
      tcube = tsq*tsince;
      tfour = tsince*tcube;
      tempa = tempa-d2*tsq-d3*tcube-d4*tfour;
      tempe = tempe+ OE.bstar*c5*(sin(xmp)-sinmo);
      templ = templ+t3cof*tcube+tfour*(t4cof+tsince*t5cof);
  }; /* End of if (isFlagClear(SIMPLE_FLAG)) */

  a =  OE.aodp*tempa*tempa;
  e =  OE.eo-tempe;

  /* A highly arbitrary lower limit on e,  of 1e-6: */
  if( e < ECC_EPS)
     e = ECC_EPS;
  xl = xmp+omega+xnode+xnodp*templ;

  /* Long period periodics */
  double axn = e*cos(omega);
  temp = 1/(a*(1.-e*e));
  double xlcof = .125 * a3ovk2 * sinio * (3+5*cosio)/ (1. + cosio);
  double aycof = 0.25 * a3ovk2 * sinio;
  double xll = temp*xlcof*axn;
  double aynl = temp*aycof;
  double xlt = xl+xll;
  double ayn = e*sin(omega)+aynl;
  double elsq = axn*axn+ayn*ayn;
  double capu = fmod( xlt - xnode, twopi);
  double chicken_factor_on_eccentricity = 1.e-6;
  double epw = capu;
  double ecosE, esinE, pl, r;
  double betal;
  double u, sinu, cosu, sin2u, cos2u;
  double rk, uk, xnodek, xinck;
  double sinuk, cosuk, sinik, cosik, sinnok, cosnok, xmx, xmy;
  double sinEPW, cosEPW;
  double ux, uy, uz;
  double tempv[3];
  int i, rval = 0;

  /* Dundee changes:  items dependent on cosio get recomputed: */
  double cosio_squared = cosio * cosio;
  x3thm1 = 3.0 * cosio_squared - 1.0;
  x1mth2 = 1.0 - cosio_squared;
  double x7thm1 = 7.0 * cosio_squared - 1.0;

         /* Added 29 Mar 2003,  modified 26 Sep 2006:  extremely    */
         /* decayed satellites can end up "orbiting" within the     */
         /* earth.  Eventually,  the semimajor axis becomes zero,   */
         /* then negative.  In that case,  or if the orbit is near  */
         /* to parabolic,  we zero the posn/vel and quit.  If the   */
         /* object has a perigee or apogee indicating a crash,  we  */
         /* just flag it.  Revised 28 Oct 2006.                     */

  if( elsq > 1. - chicken_factor_on_eccentricity)
     rval = -1;

//  if( rval)
//  {
     /*  this is Bill's fix
     for( i = 0; i < 3; i++)
     {
        rr[i] = 0.;
        vv[i] = 0.;
     }
     exit( 0);
     */
     /* my fix is to circularize the orbit on the surface of the earth */
     /* to prevent user programs from crashing.  User programs must    */
     /* recognize that pgee = 0                                        */
     /*
     unitv(rr, rr);
     unitv(vv, vv);
     cross(vv, rr, tempv);
     cross(rr, tempv, vv);
     smult(xke, vv, vv);
     exit( 0);
     */
//  }
  if( a * (1. - e) < 1. && a * (1. + e) < 1.)   /* entirely within earth */
     rval = -3;                      /* remember, e can be negative */
  if( a * (1. - e) < 1. || a * (1. + e) < 1.)   /* perigee within earth */
     rval = -4;
  /* Solve Kepler's' Equation */
  for( i = 0; i < 10; i++)
  {
     double newton_raphson_epsilon = 1e-12;
     double f, fdot, delta_epw;
     int do_second_order_newton_raphson = 1;

     sinEPW = sin(epw);
     cosEPW = cos(epw);
     ecosE = axn * cosEPW + ayn * sinEPW;
     esinE = axn * sinEPW - ayn * cosEPW;
     f = capu - epw + esinE;
     if ( fabs(f) < newton_raphson_epsilon) break;
     fdot = 1. - ecosE;
     delta_epw = f / fdot;
     if( !i)
     {
        double max_newton_raphson = 1.25 * fabs(e);

        do_second_order_newton_raphson = 0;
        if( delta_epw > max_newton_raphson)
            delta_epw = max_newton_raphson;
        else if( delta_epw < -max_newton_raphson)
                 delta_epw = -max_newton_raphson;
        else
            do_second_order_newton_raphson = 1;
     }
     if( do_second_order_newton_raphson)
         delta_epw = f / (fdot + 0.5*esinE*delta_epw);
                            /* f/(fdot - 0.5*fdotdot * f / fdot) */
     epw += delta_epw;
  }

  /* Short period preliminary quantities */
  temp = 1-elsq;
  pl = a*temp;
  r = a*(1-ecosE);
  temp2 = a / r;
  betal = sqrt(temp);
  temp = esinE/(1+betal);
  cosu = temp2 * (cosEPW - axn + ayn * temp);
  sinu = temp2 * (sinEPW - ayn - axn * temp);
  u = atan2( sinu, cosu);
  sin2u = 2*sinu*cosu;
  cos2u = 2*cosu*cosu-1;
  temp1 = ck2 / pl;
  temp2 = temp1 / pl;

  /* Update for short periodics */
  rk = r*(1-1.5*temp2*betal*x3thm1)+0.5*temp1*x1mth2*cos2u;
  uk = u-0.25*temp2*x7thm1*sin2u;
  xnodek = xnode+1.5*temp2*cosio*sin2u;
  xinck = OE.xincl+1.5*temp2*cosio*sinio*cos2u;

  /* Orientation vectors */
  sinuk = sin(uk);
  cosuk = cos(uk);
  sinik = sin(xinck);
  cosik = cos(xinck);
  sinnok = sin(xnodek);
  cosnok = cos(xnodek);
  xmx = -sinnok*cosik;
  xmy = cosnok*cosik;
  ux = xmx*sinuk+cosnok*cosuk;
  uy = xmy*sinuk+sinnok*cosuk;
  uz = sinik*sinuk;

  /* Global Position, er, and velocity, er/min */
  rr[0] = rk*ux;
  rr[1] = rk*uy;
  rr[2] = rk*uz;

  double rdot = xke*sqrt(a)*esinE/r;
  double rfdot = xke*sqrt(pl)/r;
  double xn = xke/(a * sqrt(a));
  double rdotk = rdot-xn*temp1*x1mth2*sin2u;
  double rfdotk = rfdot+xn*temp1*(x1mth2*cos2u+1.5*x3thm1);
  double vx = xmx*cosuk-cosnok*sinuk;
  double vy = xmy*cosuk-sinnok*sinuk;
  double vz = sinik*cosuk;

  vv[0] = rdotk*ux+rfdotk*vx;
  vv[1] = rdotk*uy+rfdotk*vy;
  vv[2] = rdotk*uz+rfdotk*vz;

  if(rval)
  {
	 unitv(rr, rr);
	 unitv(vv, vv);
	 cross(vv, rr, tempv);
	 cross(rr, tempv, vv);
	 smult(xke, vv, vv);    // circular velocity at surface
  }

}

void rvel(double* rr2, double* vv2)
{
    //double twopi = 2 * 3.14159265358979323846;
    int i;

    /* classical osculating orbit elements calculated from vectors rr2, vv2  */
    double xinck, xnodek, ek, mk, wk, xn, rk, uk, /*aodp,*/ pl, rdotk, rfdotk, temp;

    double h[3], n[3], vec[3], vk[3];
    smult(1. / xke, vv2, vk);
    cross(rr2, vk, h);
    pl = dot(h, h);
    double vz[] = {0, 0, 1};
    double vy[3], t;
    cross(vz, h, n);
    if(n[0] == 0. && n[1] == 0.) n[0] = 1.;
    unitv(n, n);
    rk = norm(rr2);
    rdotk = dot(rr2, vv2) / rk;
    rfdotk = norm(h) * xke / rk;
    temp = dot(rr2, n) / rk;
    uk = acose(temp);
    if(rr2[2] < 0.) uk = twopi - uk;
    cross(vk, h, vz);
    smult(-1. / rk, rr2, vy);
    vadd(vz, vy, vec);
    ek = norm(vec);
    if( ek >= 1.) return;         // open orbit
    xnodek = atan2( n[1], n[0]);
    if( xnodek < 0.)
       xnodek += twopi;
    temp = sqrt( h[0] * h[0] + h[1] * h[1]);
    xinck =  atan2( temp, h[2]);
    temp = dot(vec, n) / ek;
    wk = acose( temp);
    if(vec[2] < 0.)
      wk = fmod2ps(twopi - wk);
    OE.aodp = pl / (1. - ek*ek);
    xn = xke * pow(OE.aodp, -1.5);

    double cosio, sinio, sin2u, cos2u, temp1, temp2,
     rdot, rfdot, theta2, betal, x3thm1, x1mth2, x7thm1,
     esine, ecose, elsq, cosepw, sinepw, axn, ayn,
     cosu, sinu, capu, /*a3ovk2,*/ xlcof, aycof, aynl, xll,
     xl, a0, a1, a2, d0, d1, beta, beta2, r, u;

/*
In the first loop the osculating elements rk, uk, xnodek, xinck, rdotk,
and rfdotk are used as anchors to find the corresponding final SGP4
mean elements r, u, xnodeo, xincl, rdot, and rfdot.  Several other final
mean values based on these are also found: betal, cosio, sinio, theta2,
cos2u, sin2u, x3thm1, x7thm1, x1mth2.  In addition, the osculating values
initially held by aodp, pl, and xn are replaced by intermediate
(not osculating and not mean) values used by SGP4.  The loop converges
on the value of pl in about four iterations.
*/

    /*  seed value for first loop */
    OE.xincl = xinck;
    u = uk;

    for (i = 0; i < 99; ++i)
    {
      a2 = pl;
      betal = sqrt(pl / OE.aodp);
      temp1 = ck2  / pl;
      temp2 = temp1 / pl;
      cosio = cos(OE.xincl);
	  sinio = sin(OE.xincl);
      sin2u = sin(2.*u);
      cos2u = cos(2.*u);
      theta2 = cosio * cosio;
      x3thm1 = 3. * theta2 - 1.;
      x1mth2 = 1. - theta2;
      x7thm1 = 7. * theta2 - 1.;
      r = (rk - .5 * temp1 * x1mth2 * cos2u)
         / (1. - 1.5 * temp2 * betal * x3thm1);
      u = uk + .25 * temp2 * x7thm1 * sin2u;
      OE.xnodeo = xnodek - 1.5 * temp2 * cosio * sin2u;
      OE.xincl = xinck - 1.5 * temp2 * cosio * sinio * cos2u;
      rdot = rdotk + xn * temp1 * x1mth2 * sin2u;
      rfdot = rfdotk - xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1);
      temp = r * rfdot / xke;
      pl = temp * temp;

      // vis-viva equation
      temp = 2. / r - (rdot*rdot + rfdot*rfdot) / (xke*xke);
      OE.aodp = 1. / temp;

      xn = xke * pow(OE.aodp, -1.5);
      if(fabs(a2 - pl) < 1.e-13) break;
    }

/*
The next values are calculated from constants and a combination of mean
and intermediate quantities from the first loop.  These values all remain
fixed and are used in the second loop.
*/


    // preliminary values for the second loop
    ecose = 1. - r / OE.aodp;
    esine = r * rdot / (xke * sqrt(OE.aodp));   /* needed for Kepler's eqn.  */
    elsq = 1. - pl / OE.aodp;               /* intermediate eccentricity squared */
    //a3ovk2 = -xj3 / ck2;
    xlcof = .125 * a3ovk2 * sinio * (3. + 5. * cosio)
          / (1. + cosio);
    aycof = .25 * a3ovk2 * sinio;
    temp1 = esine / (1. + sqrt(1. - elsq));
    cosu = cos(u);
    sinu = sin(u);

/*
The second loop normally converges in about six iterations to the final
mean value for the eccentricity, eo.  The mean perigee, omegao, is also
determined.  Cosepw and sinepw are found to high accuracy and
are used to calculate an intermediate value for the eccentric anomaly,
temp2.  Temp2 is then used in Kepler's equation to find an intermediate
value for the true longitude, capu.
*/
    /*  seed values for loop  */
    OE.eo = sqrt(elsq);
    OE.omegao = wk;
    axn = OE.eo * cos(OE.omegao);

    for (i = 0; i < 99; ++i)
    {
       a2 = OE.eo;
       beta = 1. - OE.eo*OE.eo;
       temp = 1. / (OE.aodp * beta);
       aynl = temp * aycof;
       ayn = OE.eo * sin(OE.omegao) + aynl;
       cosepw = r * cosu / OE.aodp + axn - ayn * temp1;
       sinepw = r * sinu / OE.aodp + ayn + axn * temp1;
       axn = cosepw * ecose + sinepw * esine;
       ayn = sinepw * ecose - cosepw * esine;
       OE.omegao = fmod2ps(atan2(ayn - aynl, axn));
       //eo = axn / cos(omegao);
       // use weighted average to tame instability at high eccentricities
       OE.eo = .9*OE.eo + .1*(axn / cos(OE.omegao));
       if(OE.eo > .999) OE.eo = .999;
       if(fabs(a2 - OE.eo) < 1.e-13) break;
    }

    temp2 = atan2(sinepw, cosepw);
    capu = temp2 - esine;             /* Kepler's equation */
    xll = temp * xlcof * axn;

    /* xll adjusts the intermediate true longitude,  */
    /* capu, to the mean true longitude, xl          */
    xl = capu - xll;

    OE.xmo = fmod2ps(xl - OE.omegao);        /* mean anomaly */

/*
The third loop usually converges after three iterations to the
mean semi-major axis, a1, which is then used to find the mean motion, xno.
*/

    a0 = OE.aodp;
    a1 = a0;
    beta2 = sqrt(beta);
    temp = 1.5 * ck2 * x3thm1 / (beta * beta2);
    for (i = 0; i < 99; ++i)
    {
       a2 = a1;
       d0 = temp / (a0*a0);
       a0 = OE.aodp * (1. - d0);
       d1 = temp / (a1*a1);
       a1 = a0 / (1. - d1 / 3. - d1*d1 - 134. * d1*d1*d1 / 81.);
       if(fabs(a2 - a1) < 1.e-13) break;
    }
    OE.xno = xke * pow(a1 , -1.5);

} /* end rvel  */

void Satellite()
{

}

// vectors to SGP4 mean elements
void rv2el(double* rr, double* vv)
{
   double ik, ok, ek, wk, mk, nk;
   double iz, oz, ez, wz, mz, nz;
   double rr1[3], vv1[3];

   rvel(rr, vv);           // SGP4 x-elements from state vectors

   // These elements are pretty close.  Save the
   // x-elements as the k-element reference elements
   ik = OE.xincl;
   ok = OE.xnodeo;
   ek = OE.eo;
   wk = OE.omegao;
   mk = OE.xmo;
   nk = OE.xno;

   // state vectors are global so they must be stored
   // before they are changed by sgp4(0.0)
   rr1[0] = rr[0];
   rr1[1] = rr[1];
   rr1[2] = rr[2];
   vv1[0] = vv[0];
   vv1[1] = vv[1];
   vv1[2] = vv[2];

   sgp4(0.0);              // SGP4 propagation of k-elements to rr', vv'
   rvel(rr, vv);           // SGP4 of rr', vv' to x-elements'

   // first correction to k-elements is (k-elements - x-elements')
   OE.xincl  = ik + ik - OE.xincl;
   OE.xnodeo = ok + ok - OE.xnodeo;
   OE.eo     = ek + ek - OE.eo;
   OE.omegao = wk + wk - OE.omegao;
   OE.xmo    = mk + mk - OE.xmo;
   OE.xno    = nk + nk - OE.xno;

   // save z-elements. These elements are very close.
   iz = OE.xincl;
   oz = OE.xnodeo;
   ez = OE.eo;
   wz = OE.omegao;
   mz = OE.xmo;
   nz = OE.xno;

   sgp4(0.0);             // SGP4 propagation of z-elements to rr", vv"
   rvel(rr, vv);          // SGP4 of rr", vv" to x-elements"

   // second correction is small adjustment to z-elements
   // final elements are corrected by (k-elements - x-elements")
   OE.xincl  = iz + ik - OE.xincl;
   OE.xnodeo = oz + ok - OE.xnodeo;
   OE.eo     = ez + ek - OE.eo;
   OE.omegao = wz + wk - OE.omegao;
   OE.xmo    = mz + mk - OE.xmo;
   OE.xno    = nz + nk - OE.xno;

   // ensure 0 <= angle < twopi
   OE.xincl  = fabs(OE.xincl);
   OE.xnodeo = fmod2ps(OE.xnodeo);
   OE.omegao = fmod2ps(OE.omegao);
   OE.xmo    = fmod2ps(OE.xmo);

   // restore state vectors
   rr[0] = rr1[0];
   rr[1] = rr1[1];
   rr[2] = rr1[2];
   vv[0] = vv1[0];
   vv[1] = vv1[1];
   vv[2] = vv1[2];

} /* end trv2tle  */

