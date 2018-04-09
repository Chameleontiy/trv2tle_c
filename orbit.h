/* 
 * File:   orbit.h
 * Author: User
 *
 * Created on 28 Март 2018 г., 12:14
 */

#ifndef ORBIT_H
#define	ORBIT_H

//Математические константы
#define pi          ( 3.141592653589793238462643383279502884197 )
#define twopi       ( 2. * pi )
#define de2ra       ( pi / 180 )
#define nocon       ( twopi / 1440 )
#define two_thirds    ( 2. / 3. )

//Орбитальные константы
#define xj3     ( -2.53881E-6  )
#define ck2     ( 5.413079E-4  )
#define ck4     ( 6.2098875E-7 )
#define xke     ( 0.0743669161331734132 )/* = (G*M)^(1/2)*(er/min)^(3/2) where G =
                               Newton's grav const, M = earth mass */
#define xkmper  ( 6378.135 )
#define qoms2t  ( 1.880279159015270643865E-9 )
#define a3ovk2  ( -1.*xj3/ck2 )
#define s       ( 1.0122292801892716 )

#ifdef	__cplusplus
extern "C" {
#endif
    typedef struct OrbitElements {
        double jd;
        double mjd;
        double Epoch;
        double xincl;       // inclination
        double xnodeo;      // right ascension of ascending node
        double eo;          // eccentricity
        double omegao;      // argument of the perigee
        double xmo;         // mean anomaly
        double xno;         // mean motion [rad/min]
        double aodp;        // semi-major axis
        double c2;          // internal drag term for bstar conversion
        double bstar;       // bstar drag term
    } OE_t;

    typedef struct rvVec {
        double rr[3];
        double vv[3];
    } RV_t;

	double *rr2,
			*vv2;

    void sgp4(double tsince);
    void vadd(double* v1, double* v2, double* res);
    void vmadd(double* v1, double* v2, double* res, double a);
	void rvel(double* rr2, double* vv2);
	void rv2el(double* rr2, double* vv2);

#ifdef	__cplusplus
}
#endif

#endif	/* ORBIT_H */

