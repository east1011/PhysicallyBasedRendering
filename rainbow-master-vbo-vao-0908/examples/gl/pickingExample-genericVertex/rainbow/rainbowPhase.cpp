//
//  main.cpp
//  GOMSphereCPP
//
//  Created by jd on 13. 7. 11..
//  Copyright (c) 2013년 jd. All rights reserved.
//

#include "StdAfx.h"

#include <iostream>
#include <fstream>
#include <complex> // for complex numbers
#include <cmath> // cmath is the C++ wrapper of C math.h

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <iomanip>


#include "rainbowPhase.h"  // which is equivalent to "rainbowPhase.h" because
                                   // rainbowPhase.h is in the same folder as rainbowPhase.cpp
// rainbowPhase.h includes spectrum.h 

extern FILE *debugFile;

using namespace std; // complex belongs to the std name space; by Moon Jung 2013/7/17
//using namespace tr1;

double lambda_phasefunction(double lambda, double radius, double theta, double FCscat);

typedef complex<double> dcomplex;

const double lambdaStart = 400e-9;   // 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
const double lambdaEnd = 700e-9;
    
const double stepLambda = ( lambdaEnd  - lambdaStart ) / double( nSpectralSamples) ;

Spectrum particle_phasefunction( double radius, double theta, const Spectrum FCscat) {
	Spectrum phaseSpectrum; // phaseSpectrum instance is created by the default constructor of Spectrum class
	                        // The default constructor sets the wavelengths which are sampled
	                        // to represent the continuous spectrum.

	float lambdas[ nSpectralSamples ];
	float phases[ nSpectralSamples ];
	float lambda = lambdaStart;

	for (int i = 0; i < nSpectralSamples; i++ ) {
  
       
	   fprintf(debugFile, "\n lambda[%d] = %f", i, lambda);

       phases[i] = lambda_phasefunction( lambda * 1.0e-9, radius, theta, FCscat.c[i] ); 
	   lambdas[i] = lambda;

	  
	   fprintf(debugFile, "\n phases[%d] = %f, with total scat crossSection=%f", i, phases[i],
		                   FCscat.c[i] );
	   lambda += stepLambda;

	}

	return phaseSpectrum.FromSampled(lambdas, phases, nSpectralSamples); // sample from [lambdas, phases]
 	
}

/* FORTRAN subroutines work in a way of Call-by-Reference. Therefore C functions arguments should be pointers.   - Minoong Ji */

/*
 included <complex> to deal with complex numbers.
 Module is cabs(a)/cabsl(c)/cabsf(b);
 Real part is creal(a),
 Imaginary is cimag(a).
 carg(a) is for complex argument.
 - Minoong Ji */

#define NUMGAUSS 20 // used for GaussLegendre integration
#define MAXDIV 10 //  the initial guess for the number of subdivision of the interval for root finding
#define NMAXJ 100  // const used for iteration over the interface j's: 
#define NMAXROOT    1000  //  the maximum number of roots. 


/*
 PURPOSE: To calculate the various transmission and reflection coefficients
 INPUT:
 in_angle: incident angle
 ref_index: complex refractive index
 OUTPUT:
 Trans: transmission for natural light
 Reflec: reflectance for natural light
 Trp: transmission for parallel-polarized light
 Trv: transmission for vertical-polarized light
 Rep: reflection for parallel-polarized light
 Rev: reflection for vertical-polarized light
 tp: real, amplitude coefficient of transmission for parallel-polarization
 tv: real, amplitude coefficient of transmission for vertical-polarization
 rp: real, amplitude coefficient of reflection for parallel-polarization
 rv: real, amplitude coefficient of reflection for vertical-polarization
 t: transmission coefficient for natural light t=sqrt(0.5(tp*2+tv**2))
 r: reflection coefficient for natural light r=sqrt(0.5(rp*2+rv**2))
 u: mcos(theta)=u+iv
 v:
 LOCAL:
 thetat: real, refraction angle in redians */

void TransReflec (double in_angle, dcomplex ref_index, double &Trans, double &Reflec,
	double &Trp, double &Trv, dcomplex &Rep, dcomplex &Rev, dcomplex &tp, dcomplex &tv,
	 dcomplex &rp, dcomplex &rv, dcomplex &t, dcomplex &r, double &u, double &v )
{
    double re, im, pi;
    double delta, gama;
    //complex uniti;
    
    //uniti = complex(0.0e0, 1.0e0);
    
    pi = acos(-1.0e0);
    //acos : cos(theta)의역함수
    re = ref_index.real();
    im = ref_index.imag();
    delta = re*re - im*im - pow((sin(in_angle)), 2);
    
    gama = re*im;
    
    u = sqrt(2.0e0) / 2.0e0 * sqrt( sqrt(pow(delta, 2) + 4.0e0 * pow(gama, 2)) + delta);
    v = sqrt(2.0e0) / 2.0e0 * sqrt( sqrt(pow(delta, 2) + 4.0e0 * pow(gama, 2)) - delta);
    //tp = 2. * ( re + uniti * im) * cos(in_angle) / (((re*re - im*im) * cos(in_angle) + u) + uniti * (2. * re * im * cos(in_angle) + v));
    
	tp = 2. * dcomplex(re, im) * cos(in_angle) /
		     dcomplex( (re*re  - im*im) * cos(in_angle) + u,  2. * re * im * cos(in_angle) + v  );
    
    
    //tv = 2. * cos(in_angle) / (cos(in_angle) + u + uniti * v);

    tv = 2. * cos(in_angle) / dcomplex(cos(in_angle)+u, v);
    
    /*rp = (((re*re - im*im) * cos(in_angle) - u) +
          uniti * (2. * re * im * cos(in_angle) - v)) /
    (((re*re - im*im) * cos(in_angle) + u) +
     uniti * (2. * re * im * cos(in_angle) + v)); */
    
    rp = dcomplex( (re*re - im*im) * cos(in_angle) - u, 2. * re * im * cos(in_angle) - v) / 
		 dcomplex ( (re*re - im*im) * cos(in_angle) + u, 2. * re * im * cos(in_angle) + v);

    
    
    //rv = (cos(in_angle) - u - uniti * v) / (cos(in_angle) + u + uniti * v);

    rv =  dcomplex( cos(in_angle) - u, - v) / dcomplex ( cos(in_angle) + u, v);

    t = sqrt( 0.5 * (tp * tp + tv * tv) ); // cdsqrt
    r = sqrt( 0.5 * (rp * rp + rv * rv) ); // cdsqrt
    
    Trp = re * sqrt( pow( ( u * re + v * im), 2 ) + pow ( (v * re - u * im), 2 ) ) /
    ((re*re + im*im) * cos(in_angle)) * abs(tp) * abs(tp); // in the Fortran code, the function was cdabs(), which
	                                                       // computes the magnitude of a dcomplex number.
    
    Trv = re * sqrt( pow( ( u * re + v * im), 2 ) + pow ( (v * re - u * im), 2 ) ) /
    ((re*re + im*im) * cos(in_angle)) * abs(tv) * abs(tv);
    
    Rep = rp * rp;
    Rev = rv * rv;
    Trans = re * sqrt( pow( (u * re + v * im), 2 ) + pow( (v * re - u * im), 2 ) ) /
    ((re*re + im*im) * cos(in_angle)) * abs(t) * abs(t);
    
    Reflec = abs(r) * abs(r);
    
    return;
} // TransReflec

/* To calculate path length between two consecutive incident interfaces
 INPUT:
 radius: radius of the sphere
 in_angle: incident angle in radians
 ref_index: dcomplex refractive index
 OUTPUT:
 xi: path length bwtween two consecutive incident points in m
 LOCAL:
 re: real part of ref_index
 im: imaginary part of ref_index */

void path_length (double radius, double in_angle, dcomplex ref_index, double &xi)
{
    double re, im;
    double Trans, Reflec, Trp, Trv;
    dcomplex Rep, Rev, tp, tv, rp, rv, t, r;
    double u, v;
    
    re = ref_index.real();
    im = ref_index.imag();
    
    TransReflec(in_angle,ref_index,Trans,Reflec,Trp,Trv,Rep,Rev,tp,tv,rp,rv,t,r,u,v);
    
    xi = 2. * radius / (re * re + im * im) * sqrt( pow((u * re + v * im), 2) + pow ((v * re - u * im), 2) );
    
    return;
} // path_length




void errMsg ( char *MESSAG, bool FATAL )
{
    bool MsgLim, Cray;
    int MaxMsg, NumMsg;
    //       SAVE          MaxMsg, NumMsg, MsgLim ?
    
    NumMsg = 0;
    MaxMsg = 100;
    MsgLim = false;
    Cray = false;
    
    if (FATAL)
        std::cout << " ******* ERROR >>>>>>  " << MESSAG << std::endl;
    
    NumMsg = NumMsg + 1;
    
    if (MsgLim)
        return;
    
    if (NumMsg <= MaxMsg) {
        std::cout << " ******* WARNING >>>>>> " << MESSAG << std::endl;
    }
    else {
        MsgLim = true;
        std::cout << " >>>>>> TOO MANY WARNING MESSAGES -- They will no longer printed <<<<<<<" << std::endl;
    }//endif
    
    return;
}




const int   NWL = 57;
//const int   NWLT = 62;


//const int   NWL = 574;
//const int   NWLT = 62;
//double  TABIM[NWL], TABIMT[NWLT, 4], TABRE[ NWL ],  TABRET[NWLT, 4], TEMREF[4], WL[NWL], WLT[NWLT];
//double  TABIM[NWL], TABIMT[NWLT][4], TABRE[NWL],  TABRET[NWLT][4], WL[NWL], WLT[NWLT];
struct RefractiveIndex
{
    double WL;
    double TABRE;
    double TABIM;
};

//double TEMREF[5]  = { 0.0, 272.16, 268.16, 253.16, 213.16 }; // 0.0 added as the first element by Moon Jung 2013/7/17

// from 213 to 272: for the case of the drop whose size is above

/*
 struct TempRefractiveIndex
 {
 double WLT;
 struct {
 double TABRET;
 double  TABIMT;
 }  T[5];
 
 };
 struct TempRefractiveIndex1
 {
 double WLT;
 struct {
 double TABRET;
 double  TABIMT;
 }  T[2];
 
 };
 struct TempRefractiveIndex2
 {
 //double WLT;
 struct {
 double TABRET;
 double  TABIMT;
 }  T[2];
 };
 
 TempRefractiveIndex tempRefraIndex[NWLT+1];
 
 */

RefractiveIndex refraIndex [NWL+1] = {
    
    { 0, 0, 0},
    //파장, 굴절게수 실수부, 굴절계수 허수부 (0)
    
    /* Refractive Indices of distilled water from the near-infrared region to the ultraviolet region , in 24.0 celcius.
     
     from 'Measurement of the refractive index of distilled water from the near-infrared region to the ultraviolet region', Masahiko Daimon and Akira Masumura(masumura@ohara-inc.co.jp), 2006. */
    
    
    { 181.78736, 1.468492, 0.0},
    { 182.61377, 1.465508, 0.0},
    { 183.65075, 1.461964, 0.0},
    { 184.57517, 1.458976, 0.0},
    { 184.96831, 1.457748, 0.0},
    { 185.34523, 1.456594, 0.0},
    { 186.71302, 1.452589, 0.0},
    { 188.30587, 1.448250, 0.0},
    { 188.95226, 1.446576, 0.0},
    { 189.90445, 1.444204, 0.0},
    
    //DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 11, 20 ) /
    { 190.74940, 1.442176, 0.0},
    { 191.60818, 1.440187, 0.0},
    { 192.92449, 1.437273, 0.0},
    { 193.43690, 1.436177, 0.0},
    { 193.74245, 1.435535, 0.0},
    { 194.44617, 1.434084, 0.0},
    { 195.47439, 1.432034, 0.0},
    { 196.31429, 1.430413, 0.0},
    { 197.15374, 1.428844, 0.0},
    { 197.97647, 1.427348, 0.0},
    
    //DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 21, 30 ) /
    { 198.91056, 1.425704, 0.0},
    { 200.84066, 1.422463, 0.0},
    { 204.22306, 1.417245, 0.0},
    { 207.60622, 1.412533, 0.0},
    { 210.40121, 1.408969, 0.0},
    { 213.13806, 1.405734, 0.0},
    { 217.53669, 1.400987, 0.0},
    { 222.57015, 1.396146, 0.0},
    { 226.34187, 1.392869, 0.0},
    { 231.16733, 1.389045, 0.0},
    
    //DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 31, 40 ) /
    { 239.02641, 1.383575, 0.0},
    { 244.08006, 1.380471, 0.0},
    { 248.79191, 1.377822, 0.0},
    { 253.99695, 1.375140, 0.0},
    { 262.88107, 1.371063, 0.0},
    { 273.47659, 1.366882, 0.0},
    { 283.11246, 1.363603, 0.0},
    { 289.44400, 1.361673, 0.0},
    { 296.814, 1.359619, 0.0},
    { 312.657, 1.355795, 0.0},
    
    //DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 41, 50 ) /
    { 334.244, 1.351588, 0.0},
    { 365.119, 1.346994, 0.0},
    { 404.770, 1.342710, 0.0},
    { 435.957, 1.340179, 0.0},
    { 480.126, 1.337418, 0.0},
    { 486.269, 1.337091, 0.0},
    { 546.227, 1.334435, 0.0},
    { 587.725, 1.333012, 0.0},
    { 644.025, 1.331436, 0.0},
    { 656.454, 1.331126, 0.0},
    
    //DATA  ( WL( I ), TABRE( I ), TABIM( I ) , I = 51, 60 ) /
    { 706.714, 1.329996, 0.0},
    { 780.237, 1.328583, 0.0},
    { 852.344, 1.327365, 0.0},
    { 894.596, 1.326702, 0.0},
    { 1014.26, 1.324921, 0.0},
    { 1083.33, 1.323902, 0.0},
    { 1128.95, 1.323216, 0.0}
};




void  refice( double WAVMET, double TEMP, dcomplex  &ref_ice )  {
    
    /*
     c        Calculates complex refractive index of Ice 1H for wavelengths
     c        between 45 nm and 8.6 m.  ****For wavelengths above 167 microns,
     c        temperature dependence is included for temperatures between
     c        213 and 272K. ***  Mainly intended for applications in Earth ice
     c        clouds and snow, not other planets or interstellar space;
     c        the temperature dependence or crystalline form of ice may be
     c        incorrect for these latter applications.
     c      I N P U T :
     c                   WAVMET = wavelength (meters)
     c                   WAVLEN = wavelength (microns)
     c                            (range:  0.0443 to 8.600E+06)
     c                   TEMP   = temperature (K) ( for WAVLEN.GT.167 only )
     c                            (range:  213.16 to 272.16)
     c      O U T P U T :  ref_ice = complex refractive index
     c                              ( with positive imaginary part )
     c      (WARNING:  input out of range will print a warning message and
     c                 return ref_ice=(0.,0.) in order not to unnecessarily
     c                 halt the calling program;  the calling program should
     c                 test the real part of ref_ice to catch these errors)
     c      METHOD :  Tabular interpolation, assuming
     c                (1) real index is linear in log(wavelength)
     c                    and linear in temperature
     c                (2) log(imag. index) is linear in log(wavelength)
     c                    and linear in temperature
     c     AUTHORS OF subroutine refice( WAVMET, TEMP, ref_ice ):
     c             Stephen Warren, Univ. of Washington (1983)
     c               (sgw@cloudy.atmos.washington.edu)
     c             Bo-Cai Gao, JCESS, Univ. of Maryland (1995)
     c               (gao@imagecube.gsfc.nasa.gov)
     c             Warren Wiscombe, NASA Goddard (1995)
     c               (wiscombe@climate.gsfc.nasa.gov)
     c     MODIFICATIONS IN 1995 :
     c       Gao, Warren, and (to a small extent) Wiscombe modified the
     c       original Warren refice program from 1984 to change values of
     c       imaginary refractive index in the 0.161-0.410 and 1.445-2.50
     c       micron regions.  The values in 0.161-0.410 were incorrect and
     c       the values in 1.445-2.50 were among the most uncertain in 1984.
     c       New measurements have made it possible to improve both regions.
     
     c       No changes were made to real refractive indices (by re-doing a
     c       Kramers-Kronig analysis), because the values of imaginary
     c       index MIM involved are so small (below 0.001) that the
     c       resulting changes in real index MRE would be in the third
     c       decimal place at most.  (MIM has negligible effect on MRE
     c       when MIM << MRE.)
     
     c       The 0.161-0.410 micron region was changed using //DATA provided
     c       by Warren, which correct his misinterpretation of Minton's
     c       measurements for 0.181-0.185 micron, and incorporate new
     c       measurements of Perovich and Govoni (1991) for 0.250-0.400
     c       micron.  Warren (1984) correctly represented UV measurements
     c       of Seki et al. and visible measurements of Grenfell/Perovich,
     c       but he plotted Minton's measurements a factor of 2.3 too low
     c       because he misinterpreted base-10 as base-e.  (The UV
     c       measurements of Dressler/Schnepp and Shibaguchi et al are also
     c       probably expressed as absorption coefficients on base-10;
     c       therefore those values also were probably plotted a factor of
     c       2.3 too low in Warren's (1984) Figure 2.)
     
     c       The details of how the present imaginary index //DATA for
     c       0.161-0.410 micron is obtained are as follows.  Point A in
     c       Warren's Figure 2 at 161 nm is joined with a straight line to
     c       Minton's corrected point B at 181 nm.  Minton's reported
     c       values for 181-185 nm have been smoothed within his stated
     c       uncertainty.  Now a smooth curve is drawn to join Minton at
     c       185 nm to Perovich/Govoni (PG) at 250 nm.  PG's values from
     c       their Table 1 show some unrealistic wiggles that are smaller
     c       than their error bars, so a smooth curve was fitted through
     c       them and values were taken from the smoothed curve at 10-nm
     c       intervals.  PG ends at 400 nm, where Grenfell/Perovich (GP)
     c       starts.  At 400 nm we take imaginary index=2.82E-9, the
     c       average of PG (2.93E-9) and GP (2.71E-9).
     
     c       The Warren (1984) values of imaginary index in the 1.445-2.50
     c       micron region were replaced by those of Kou et al.(1994).  In
     c       order to remove the resulting discontinuities near 1.445 and
     c       2.5 micron, the Warren values at 1.43 and 1.44 micron were
     c       changed to 0.9E-04 and 1.3E-04 respectively, and his values at
     c       2.52, 2.55, and 2.565 micron were changed to 8.255E-04,
     c       8.578E-04}, {and 8.739E-04}, {respectively. The latter change
     c       eliminated a small local maximum at 2.5 micron which was not
     c       realistic and has never been seen in spectra of snow bracketing
     c       that wavelength.
     
     
     c     REFERENCES :
     
     c       Warren, S., 1984: Optical Constants of Ice from the Ultraviolet
     c          to the Microwave, Appl. Opt. 23, 1206-1225
     
     c       Kou, L., D. Labrie, and P. Chylek, 1994: Refractive indices
     c          of water and ice in the 0.65- to 2.5-micron spectral range,
     c          Appl. Opt. 32, 3531-3540
     
     c       Perovich, D., and J. Govoni, 1991: Absorption Coefficients
     c          of Ice from 250 to 400 nm, Geophys. Res. Lett. 18, 1233-1235
     c ======================================================================
     
     */
    //      IMPLICIT NONE
    
    
    //      PARAMETER ( NWL = 574, NWLT = 62)
    //c     ..
    //      const int   NWL = 574;
    //      const int   NWLT = 62;
    /* The above has moved outside of the function, because they are needed for the common block declaration*/
    
    
    //c     .. Scalar Arguments ..
    
    //      REAL*8 TEMP, WAVLEN,WAVMET
    double WAVLEN;
    
    char MESSAG[40];
    bool PASS1;
    int i;
    double    FRAC, MIM, MRE, YHI, YLO;
    //c     ..
    //c     .. Intrinsic Functions ..
    
    //      INTRINSIC LOG, CMPLX, EXP
    //c     ..
    
    //c** Refractive index table
    
    //      double  TABIM[NWL], TABIMT[NWLT, 4], TABRE[ NWL ],  TABRET[NWLT, 4], \
    //	          TEMREF[4], WL[NWL], WLT[NWLT];
    
    //      COMMON / ICEREF / WL, WLT, TABRE, TABRET, TABIM, TABIMT, TEMREF
    /* COMMON declaration is the same as global declaration, so we have moved the above declaration outside of the function */
    
    PASS1 = true;
    
    WAVLEN = WAVMET * 1.0e+9;  // meter to nano meter
/*
    if ( PASS1 )  {
        
        PASS1 =  false;
        
         sprintf(buff, "%s", "Hello");
 
        
        //c    ** Superficially test if main table messed up
        
        if ( NWL < 100 )
            errMsg( "refice--NWL value bad", true);
        if ( refraIndex[1].WL > 0.045 )
            errMsg( "refice--WL(1) bad", true);
        if ( refraIndex[NWL].WL < 166.)
            errMsg("refice--WL[NWL] bad", true);
 
        
        for (i = 1 ; i <= NWL ; i++ ) {
            
            if ( refraIndex[i].WL <= 0.0 || refraIndex[i].TABRE <= 0.5 || refraIndex[i].TABRE > 2.0
                ||  refraIndex[i].TABIM < 0.0 || refraIndex[i].TABIM  > 10.0 )  {
                sprintf( MESSAG, "%s %d %s",  "refice--table value  ", i,"  out of bounds " );
                
                errMsg( MESSAG, true);
                
            }
            
            if ( i > 1 &&  refraIndex[i].WL <= refraIndex[i-1].WL  )  {
                sprintf( MESSAG, "%s %d %s", "refice--table WL[", i, "] not increasing  ");
                errMsg( MESSAG, true);
            }
            
        }// endfor
         

    }// endif
*/
    //    std::cout << "WAVELEN :" << WAVLEN << "  refraIndex[1].WL :" << refraIndex[1].WL << std::endl;
    //    std::cout << "NWLT :" << NWLT << "  refraIndex[NWLT].WL :" << refraIndex[NWLT].WL << std::endl;
    
   // if ( WAVLEN < refraIndex[1].WL || WAVLEN > tempRefraIndex1[NWLT].WLT  ) {

	if ( WAVLEN < refraIndex[1].WL  ) {
        
        errMsg("refice--wavelength outside table boundaries", false);
        ref_ice = dcomplex(0.,0.);
        return;
    }
    
	// temperature dependency ignored

    //if ( WAVLEN <= 167.) {
        //c                                  ** Wavelength between 0.045 and 167
        //c                                  ** microns. No temperature dependence
     for (i = 2 ; i <= NWL ; i++ ) {    // DO 10 I = 2, NWL
            if ( WAVLEN <= refraIndex[i].WL)
                break;   
        }
        
		// interpolate the refractive index of WAVLEN by using the i-1 th and ith refractivce index
		// These two refractive indices are those whose wavelengths include the given WAVLEN

      FRAC   = log ( WAVLEN / refraIndex[i-1].WL ) /  log ( refraIndex[i].WL  / refraIndex[i-1].WL  );
      MRE    = refraIndex[i-1].TABRE + FRAC * ( refraIndex[i].TABRE - refraIndex[i-1].TABRE );
    
    MIM = 0.0e0;
    //MIM    = refraIndex[i-1].TABIM * pow( ( refraIndex[i].TABIM / refraIndex[i-1].TABIM ) , FRAC );
        
    //}
    //else {
        //c               ** Wavelength greater than 167 microns
        //c               ** (temperature-dependent case)
        //c        write (*,*)'temp=',TEMP
        
    //    std::cout << "temp=" << TEMP << std::endl;
        
    //    if ( TEMP < TEMREF[4] || TEMP > TEMREF[1] ) { // the index 3, 0 were incremented by Moon Jung , 2013/7/17
    //
            
     //       errMsg("refice--temperature outside table boundaries", false);
            
     //       ref_ice = dcomplex(0.,0.);
     //       return;
            
     //   }
        //c                         ** Find position in temperature array
       
	  //int t;
      //for (t = 2; t <= 4; t++ ) {  // DO 30 L = 2, 4  //?question - why t = 1? // updated by Moon Jung 2013/7/17
      //      if ( TEMP >= TEMREF[t] )
      //          break;
      //
        
      //  //c                         ** Find position in wavelength array
      //  for (i = 2 ; i <= NWLT; i++ ) {  //DO 50 I = 2, NWLT  //?question - why i = 1?
       //     if ( WAVLEN <= tempRefraIndex[i].WLT )
       //         break;
       // }
        // t may refer to 4 or less
        // i may refer to NWLT or less, because of the break from the loop
        
		// interpolate the time dependent refractive index of WAVLEN

       // FRAC   = log( WAVLEN / tempRefraIndex[i-1].WLT ) / log( tempRefraIndex[i].WLT / tempRefraIndex[i-1].WLT );
        
       // YLO    = tempRefraIndex[i-1].T[t].TABRET + FRAC * ( tempRefraIndex[i].T[t].TABRET - tempRefraIndex[i-1].T[t].TABRET );
        
       // YHI    = tempRefraIndex[i-1].T[t-1].TABRET +  FRAC * ( tempRefraIndex[i].T[t-1].TABRET - tempRefraIndex[i-1].T[t-1].TABRET );
        
       // MRE    = YLO + ( YHI - YLO) * ( TEMP - TEMREF[t] ) 
	//	/   ( TEMREF[t-1] - TEMREF[t] );
        
      //  YLO    = log( tempRefraIndex[i-1].T[t].TABIMT) +    FRAC * log( tempRefraIndex[i].T[t].TABIMT / tempRefraIndex[i-1].T[t].TABIMT );
        
        //YHI    = log( tempRefraIndex[i-1].T[t-1].TABIMT ) +    FRAC * log( tempRefraIndex[i].T[t-1].TABIMT / tempRefraIndex[i-1].T[t-1].TABIMT );
        
        //MIM    = exp( YLO + (YHI - YLO) * (TEMP - TEMREF[t]) /  (TEMREF[t-1] - TEMREF[t]) );
        
    //} // endif
    
    ref_ice = dcomplex(MRE, MIM); // complex() is a constructor of a complex class
    
} // refice


/* calculate absorption coefficient
 INPUT:
 lambda:wavelength (m)
 ref_index: complex refractive index
 OUTPUT:
 beta: absorption coefficient (1/m) */
void abs_coef(double lambda, dcomplex ref_index, double &beta)
{
    double pi, im;
    pi = acos(-1.0e0);
    
    im = ref_index.imag();
    
    beta = 4.0 * pi * im / lambda;
    
    return;
} // abs_coef

double absFF(double x, double lambda, double radius)
{
    dcomplex m;
    double beta;
    double xi, in_angle, TT, RR, Trp, Trv;
    
    /* 원래 u, v도 complex지만 TransReflec의 interface와 맞추기 위해 double로 수정 */
    dcomplex Rep, Rev, tp, tv, rp, rv, t, r;
    double u, v;
    
    double absff;

	const double  pi = acos( -1.0);

	const int CC = 30;

    double tempKK = CC + 273; // CC to Kelvin temperature. 
    
    //To obtain refractive index m
    refice( lambda, tempKK, m );
    
    //To obtain absorption coefficient beta
    abs_coef(lambda, m, beta);
    
    //To obtain a path length xi, get the incidence angle from relation x = u_i = cos (in_angle)
    in_angle = acos(x);
    
    path_length(radius, in_angle, m, xi);

    TransReflec(in_angle,m,TT,RR,Trp,Trv,Rep,Rev,tp,tv,rp,rv,t,r,u,v); // compute TT and RR for absff
    
    absff = 2.0 *  TT / (1.0 - RR * exp(-beta * xi) ) * (1.0e0 - exp(-beta * xi) ) * x;   // eq (21a) in GOMsphere document
    
    return absff;
    
} // absFF




void GaussLeg (int n, double x[], double w[])
{
    double x1, x2;
    int m,j,i;
    double eps,z1,z,xm,xk,pp,p1,p2,p3,pi;
    eps = 1.0e-10;

    m  = (n+1)/2;

    x1 = -1.0e0;
    x2 = 1.0e0;
    xm = 0.5e0 * (x2 + x1);
    xk = 0.5e0 * (x2 - x1);
    pi = acos(-1.0e0);
    z1 = 0.0e0;
    for (i = 1 ; i <= m ; i ++)
    {
        z = cos( pi * ( (double)i - 0.25e0)/( (double)n + 0.5e0));
        while ( abs(z-z1) > eps)
        {
            p1 = 1.0e0;
            p2 = 0.0e0;
            
            for (j = 1; j <= n; j++){
                p3 = p2;
                p2 = p1;
                p1 = ((2.0e0 * (double)(j) - 1.0e0) * z * p2 - ((double)(j) - 1.0e0) * p3) / (double)(j);
            }
            pp = (double)(n) * ( z * p1 - p2 ) / ( z * z - 1.0e0);
            z1 = z;
            z = z1 - p1 / pp;
            
        }
        /*
         x(i)=xm-xk*z
         x(n+1-i)=xm+xk*z
         w(i)=2.0*xk/((1.0d0-z*z)*pp*pp)
         w(n+1-i)=w(i)
         */
        x[i] = xm - xk * z;
        x[n+1-i] = xm + xk * z;
        w[i] = 2.0e0 * xk / ( (1.0e0 - z * z) * pp * pp );
        w[n+1-i] = w[i];
    }
    
    return;
} // GaussLeg



void sumNumber1 ( double beta, double xi, double r, double eps, int &sumNumber)
{
    double tem;
    tem = 0.0e0;
    tem = beta * xi * 0.5 - 2. * log(r) - log(eps);
    tem = tem / (beta * xi * 0.5 - log(r));
    sumNumber = (int)tem;
    
    if (sumNumber < 2)
        sumNumber = 2;
    
    return;
}


void sumNumber2 ( double beta, double xi, double r, double t, double eps, int &sumNumber)
{
    double tem;
    tem = 0.0e0;
    tem = beta * xi * 0.5 + log(t) - 2. * log(r) - log(eps);
    tem = tem / (beta * xi * 0.5 - log(r));
    sumNumber = (int)tem;
    
    if (sumNumber < 2)
        sumNumber = 2;
    
    return;
}

void sumNumber3 ( double beta, double xi, dcomplex t, double Reflec, double Trans, double eps, int &sumNumber)
{
    double tem;
    tem = 0.0e0;
    tem = beta * xi + 2.0 * log(abs(t)) - 2. * log(Reflec) - log(eps);
    tem = tem / (beta * xi - log(Reflec) );
    sumNumber = (int)tem;
    
    if (sumNumber < 2)
        sumNumber = 2;
    
    return;
}


// scatFNF: scattering function near field given cos theta_i

double scatFNF(double x, double lambda, double radius) //not used
{
    int NN;
    dcomplex m;
    double  beta, thetai, xi, eps;
    dcomplex Rep, Rev, tp, tv, rp, rv, t, r;
    double Trans, Reflec, Trp, Trv, u, v;
    
    double scatfnf = 0.0e0;
    eps = 1.0e-15;
 
	const double CC = 30.0; 

	double tempKK = CC + 273.0;
    
    refice (lambda, tempKK, m);
    abs_coef(lambda, m, beta);

    thetai = acos(x);
    path_length(radius, thetai, m, xi);
    
    //cout << "radius: "<< radius << " Reflec: " << Reflec << " NN:" << NN << endl;
    TransReflec(thetai, m, Trans, Reflec, Trp, Trv, Rep, Rev, tp, tv, rp, rv, t, r, u, v);

    
    sumNumber3(beta, xi, t, Reflec, Trans, eps, NN); // NN = NN(theta_i): Theta_i corresponds to the points on the sphere which the incident light intersects
	                                                 //    All these points might contribute to the scattering of the incident light.                
    
    //cout << "xi: "<< xi << " Trans: " << Trans << " Reflec: " << Reflec << " NN:" << NN << endl;
    
    for (int j = 2 ; j <= NN ; j++)
        scatfnf = scatfnf +  (Trp * Trp * pow( (abs(Rep)), j-2) + Trv*Trv* pow( (abs(Rev)), j-2) ) * exp(-beta*xi*(double)(j-1)); // eq. (23b)
    
    scatfnf = scatfnf +  2. * Reflec; // eq. (43), cf. (23a) + (23b); Reflec = the term when j=1.
    scatfnf = scatfnf * x; // x = u_i = cos theta_i
    
    return scatfnf;
} // not used



/* scat_efficiency
 c PURPOSE:
 c     To calculate scattering efficiency for both near filed (QNscat) and far filed
 c     (QFscat), given wavelength, radius of sphere, number of roots for Gaussian
 c     quadratures. eq.(44) { eq.(48)
 c REFERENCE:
 c     Davis P J, P Rainowitz, "Methods of numerical integration",Second edition,
 c     Orlando, FL: Academic Press,1984. p.481-483.
 c INPUT:
 c     lambda: wavelength
 c     radius: radius of sphere
 c     m: refractive index
 c OUTPUT:
 c     QFscat: scattering efficiency for far field
 c     QNscat: scattering efficiency for near field
 */
void scat_efficiency(double lambda, double radius, dcomplex m, double &QNscat, double &QFscat){
    
    double a, b;
    //const int n = NUMGAUSS;
    double pi;
    double x[ NUMGAUSS +1], w[NUMGAUSS +1], xm, xk, dx;
    //double xm_minus_dx, xm_plus_dx;
    //double (*_scatFNF)(double, double &, double &) = scatFNF;
    int mm, k, i;
    //n = 20;
    
    pi = acos(-1.0e0);
    a = 0.0e0;
    b = 1.0e0;
    
    GaussLeg(NUMGAUSS, x, w);
    xm = 0.5e0 * (b+a);
    xk = 0.5e0 * (b-a);
    mm =  ( NUMGAUSS  + 1 )/2;
    k = 2 * mm - NUMGAUSS + 1;
    
    QNscat = 0.0e0;
    QFscat = 0.0e0;
    
    for (i = mm + 1; i <= NUMGAUSS ; i++)
    {
        dx = xk * x[i];
        QNscat = QNscat + w[i] * (scatFNF(xm+dx, lambda, radius) + scatFNF(xm-dx, lambda, radius)); // compute the integrand of eq (43)
        // used to compute Q_{sca}^{N}
        
        //cout << dx << "     " << xk << "     " << x[i] << "     " << w[i] << endl;
        
    }
    
    
    
    QNscat = QNscat * xk;  // QNscat = Q_{sca}^N
    QFscat = QNscat + 1.0e0;  // cf. eq (49) Q_{sca}^F = Q_{sca}^F + 1
    
    //cout << QNscat << "     " << QFscat << endl;
    
    return;
}



/* ScatFieldAmp
 
 c PURPOSE:
 c     To evaluate amplitude of scattered electric field in term of incident amplitude.
 c  REFERENCE:
 c   Zhou, X., S. Li, and K. Stamnes, A new geometrical optics code for computing optical properties of large
 c   dielectric spheres, Applied Optics, 42 (21), 4295-4306, 2003.
 c  INPUT VARIABLES:
 c      thetai: real, the incident angle (in radians)
 c      m: complex, refractive index
 c      tp: real, amplitude coefficient of transmission for parallel-polarization
 c      tv: real, amplitude coefficient of transmission for vertical-polarization
 c      rp: real, amplitude coefficient of reflection for parallel-polarization
 c      rv: real, amplitude coefficient of reflection for vertical-polarization
 c      beta: absorption coefficient
 c      xi: path length between two reflection events inside the sphere
 c      j: integer, jth interface of a ray inteacting with the sphere interface
 c  OUTPUT VARIABLES:
 c      epspjp: real, amplitude of scattered parallel-polarized electric field
 c      epspjv: real, amplitude of scattered vertically-polarized electric field
 c  LOCAL VARIABLES:
 c       thetat: real, refraction angle (in radians)
 
 */

void ScatFieldAmp (double thetai, dcomplex tp, dcomplex tv, dcomplex rp, dcomplex rv, 
	double beta, double xi, int j, double u, double v, dcomplex &epspjp, dcomplex &epspjv)
{
    dcomplex a;
    a = dcomplex(u,v);
    if (j == 1) {
        epspjp = rp;
        epspjv = rv;
    }
    else {
        epspjp = a / cos(thetai) * tp * tp * pow( (-rp), j-2) * exp( -beta * xi * (j-1) / 2.0e0 );  //notsure
        epspjv = a / cos(thetai) * tv * tv * pow( (-rv), j-2) * exp( -beta * xi * (j-1) / 2.0e0 );  //notsure
    }
    
    return;
} // ScatFieldAmp


/* BessJ1
 c Purpose:
 c    Returns the value of the first order Bessel function J1(x) for any real x
 c Reference:
 c    Press W H , S A Teukolsky, W T Vetterling and B P Flannery,"Numerical
 c    recipes in C: the art of scientific computing", Second edition, Cambridge
 c    University Press,1992. p.233.
 */

double BessJ1 ( double x )
{
    double ax, z;
    double y, ans1, ans2, xx;
    ax = abs(x);
    double bessj1;
    if (ax < 8.0) {
        y = x * x;
        ans1 = x*(72362614232.0+y*(-7895059235.0+y*(242396853.1+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
        ans2 = 144725228442.0+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y*1.0))));
        bessj1 = ans1 / ans2;
    }
    else {
        z = 8.0 / ax;
        y = z*z;
        xx = ax-2.356194491;
        ans1 = 1.0+y*(0.183105e-2+y*(-0.3516396496e-4+y*(0.2457520174e-5+y*(-0.240337019e-6))));
        ans2 = 0.04687499995+y*(-0.2002690873e-3+y*(0.8449199096e-5+y*(-0.88228987e-6+y*0.105787412e-6)));
        bessj1=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
        if(x < 0.0e0)
            bessj1 = -bessj1; //notsure
        //bessj1 =- bessj1;
    }
    
    return bessj1;
    
} // BessJ1

void rtfunction (int j, dcomplex m, double theta, double T, double &f, double &df)
{
    
    double re;
    double thetat, thetatd;
    
    re = m.real();
    thetat = asin(2.0 * T / re / (1.0 + T * T));
    
    thetatd = sqrt( re * re * pow( (1.0 + T * T), 2) - 4.0 * T * T );
    
    f = pow ( (-1.0e0), j-2 ) * ((1. - 6. * T * T + pow(T, 4)) * cos(2. * (double)(j-1) * thetat) + 4. * T * (1. - T * T) * sin(2. * (double)(j-1) * thetat)) / pow( (1. + T * T ), 2) - cos(theta);
    
    df = 4. * pow( (-1.e0), j-2) / pow( (1. + T * T), 3) * (1. + (double)(j-1) * (T * T - 1.) / thetatd ) * (4. * T * (T * T - 1.) * cos(2. * (double)(j-1) * thetat) + (1. - 6. * T * T + pow( T, 4 )) * sin(2.0 * (double)(j-1) * thetat));
    
    return;
} // rtfunction

void brackets (int j, dcomplex m, double theta, double x1, double x2, int n, double xb1[], double xb2[], int &nb)
{
    int nbb, i;
    double x, fp, fc, dx, f, df;
    
    nbb = 0;
    
    dx = (x2-x1)/(float)n;
    
	x = x1;
    rtfunction(j,m,theta,x,f,df);
    fp = f;
    
	for ( i = 1 ; i <= n; i++ ) {
        
        x = x + dx;  // dx is changed as n is changed => the number of roots is potentially n = NMAXROOT
        rtfunction(j,m,theta,x,f,df);
        fc = f;
    
		if ((fc*fp) <= (-1.0e-15)) {
            nbb = nbb + 1;
            xb1[nbb] = x - dx;
            xb2[nbb] = x;
            // maximum number of roots is reached
            if (nbb == nb) // the actual number of roots could reach nb = NMAXROOT, because n will be
				           // incremented if roots are not found in the caller of brackets
                return;
        }
        fp = fc;
    }
    
	nb = nbb; // nb, which was given as NMAXROOT (the bound of the array xb1) called brackets called and is set to the actual number of
	           // roots found. The array access does not exceed the bound of the array, because n will never be
	            // greater than NMAXROOT in the caller of brackets.
    return;
} // brackets

void rootfind (int j, dcomplex m, double theta, double x1, double x2, double xacc, double &rt, int &flag)
{
    
    int i, maxit;
    double df, dx, dxold, f, fl, fh, temp, xl, xh;
    
    maxit = 10000;

    flag = 0;
    
    rtfunction(j, m, theta, x1, fl, df);
    rtfunction(j, m, theta, x2, fh, df);
    
    if (fl * fh > 0.0e0)
        errMsg("rootfind -- root must be bracketed", true);
    
    if (fl == 0.0e0) {
        rt = x1;
        return;
    }
    
    if (fh == 0.0e0) {
        rt = x2;
        return;
    }
    
    if (fl < 0.0e0) {
        xl = x1;
        xh = x2;
    } else {
        xl = x2;
        xh = x1;
    }
    
    rt = 0.5 * (x1 + x2);
    dxold = abs(x2 - x1);
    dx = dxold;
    
    rtfunction(j, m, theta, rt, f, df);
    
    for ( i = 1 ; i <= maxit ; i ++) {
        if ((((rt-xh)*df-f) * ((rt-xl)*df-f) > 0.0e0) || (abs(2.0*f) > abs(dxold * df))) {
            dxold = dx;
            dx = 0.5 * (xh - xl);
            rt = xl + dx;
            if (xl == rt)
                return;
        } else {
            dxold = dx;
            dx = f/df;
            temp = rt;
            rt = rt - dx;
            if (temp == rt)
                return;
        }//endif
        
        if (abs(dx) < xacc)
            return;
        rtfunction(j, m, theta, rt, f, df);
        if (f < 0.0e0)
            xl = rt;
        else
            xh = rt;
    }//endfor
    
    flag = 1;
    return;
}

void roots(int j, dcomplex m, double theta, double x1, double x2, double xacc, int nsegments, int &nb, double rts[])
{
    double xb1[ NMAXROOT + 1], xb2[ NMAXROOT + 1 ], rt;

    int i, ll;
    //int nmax, flag;
    int flag;

    //nmax = 100; changed to Macro constant NMAXROOT; by Moon Jung, 2013/7/17

    int n, nm;
    
    nb = NMAXROOT;  // by Moon Jung, 2013/7/17
    n = nsegments;

    brackets(j, m, theta, x1, x2, n, xb1, xb2, nb); // find nb brackets { [ xb1[nb], xb2[nb] ] }
	                                                // within each of which there exists a root
    
   
    
    while (nb == 0 && n <= NMAXROOT) { //by Moon Jung, 2013/7/17
        n = n + 1;                                       // n in incremented from nsegments
        nb = NMAXROOT; // by Moon Jung, 2013/7/17
        brackets(j, m, theta, x1, x2, n, xb1, xb2, nb);  
    }
    
    if (nb == 0 || n > NMAXROOT) {
        nb = 0;
        return;
    }

    ll = 0;
    if (nb >= 1) {
        
		for (i = 1 ; i <= nb ; i++) {
            rootfind (j, m, theta, xb1[i], xb2[i], xacc, rt, flag);
            if (flag == 1)
                continue;
            
            if (nb >= 2) {
                for (nm = 1 ; nm <= i-1 ; nm++) {
                    if (rt == rts[nm]) {
                        break;
                    }
                }//endfor
            }//endif
        
			ll = ll + 1;
            rts[ll] = rt;
        
		}//endfor
    }//endif
    //LINE_TEN:    //notsure
    nb = ll;
    
    return;
} // roots()

/*  phasefunction
 c  Purpose:
 c     To obtain phase function for given wavelength and radius of a sphere
 c  Input:
 c     lambda: wavelength, in meters
 c     radius: radius of the scatter, in meters
 c     theta: scattering angle
 c     QFscat: scattering efficiency for far field used in scat_efficiency(...)in Q_scat2.f
 c     QNscat: scattering efficiency for near field used in scat_efficiency(...)in Q_scat2.f
 c  Output:
 c     Nphase: phase function for near field scattering
 c     Fphase: phase function for far field scattering
 */
double lambda_phasefunction(double lambda, double radius, double theta, double FCscat)
{
    double PQN, PQF, a, b, TEMP, xacc, epsilon, x1, x2;
    double thetai, Trans, Reflec, Trp, Trv;
    double beta, xi, pi, x, u, v;// BessJ1;
    dcomplex Rep, Rev, tp, tv, rp, rv, t, r, thetat;
    dcomplex m, epspjp, epspjv;
    int N, j, nb, nsegments;
    //const int Nmax = 100;
    
    double re, D;
    int nbb, i, mm, N2;
    double rts[MAXDIV+1];
    
    nsegments = MAXDIV;
    
    //nb = NMAXROOT; // nb is set to NMAXROOT within the function roots()
    
    xacc = 1.0e-15;
    x1 = 0.0e0;
    x2 = 1.0e0;
    pi = acos(-1.0e0);
    
    
    // To obtain dcomplex refractive index (m). subroutine refice( ... ) is in ref_ice.f
    TEMP = 270.0;
    refice( lambda, TEMP, m );
    
    re = m.real();

    PQN = 0.0e0;
    PQF = 0.0e0;
    mm = 0;
    j = 0;
    
	// compute the numerator of eq (40)

    while (j <= NMAXJ ){
        j = j + 1;
        roots(j, m, theta, x1, x2, xacc, nsegments, nb, rts);  // to solve equation eq (45) in GOMSphere doc. solution is within [x1, x2] = [0,1]
        
        if (nb == 0)
            continue;
        
        
        for (i = 1 ; i <= nb ; i ++) {
            thetai = 2.0 * atan(rts[i]);
            abs_coef(lambda, m, beta);
            path_length(radius, thetai, m, xi);

            TransReflec(thetai, m, Trans, Reflec, Trp, Trv, Rep, Rev, tp, tv, rp, rv, t, r, u, v);  // compute things needed to compute epspjp and epspjv

            ScatFieldAmp(thetai, tp, tv, rp, rv, beta, xi, j, u, v, epspjp, epspjv);
            
            D = sin(2.*thetai) / (4.0 * sin(theta)) / abs((double)(j-1) * cos(thetai) / sqrt(u*u + v*v) -1.0);
            
			PQN = PQN + 2.0 *  (pow( abs(epspjp), 2) + pow(abs(epspjv), 2)) * D; // eq (40)
            mm = mm + 1;
        }//9 endfor
        
    }//10 endwhile
    
		
    x = 2. * pi * radius / lambda;
    
	if (theta < pi/2.) {
        PQF = PQN + 4.0 * BessJ1( x * sin(theta) ) * BessJ1( x * sin(theta) ) / ( sin(theta) * sin(theta) ); // eq (40) note that the Bessel function was integrated from 0 to pi/2
    } else { // the diffusion component is assumed to be zero

        PQF = PQN;
	}

	double QFscat = FCscat / ( pi * (radius * radius) ); // Q_sca = C_sca / pi * r^2

	double phase = PQF / QFscat;

	double standardPhase = phase / ( 4.0 * pi ); // phase computed by this code integrated over sphere to 4 * pi, the standard version is integrated to 1.

    return standardPhase;

} // phasefunction

void scattering_function (double theta, dcomplex m, double lambda, double radius, double &scatfN ) {
    
    double beta, xi;
    dcomplex thetat;
    int NN;
    int nb, nsegments;
    double thetai, Trans, Reflec, Trp, Trv;
    dcomplex Rep,Rev,tp,tv,rp,rv,t,r;
    dcomplex epspjp,epspjv;
    
   // const int Nmax = 1000; // commented out by Moon Jung, 2013/7/17; moved to the global Macro definition NMAXJ.
	                         // As the interface number, 1000 is too much; it is set to 100.
    
    double re, D;
    int nbb, i, j, mm, N2;
    double xb1[MAXDIV+1], xb2[MAXDIV+1], rts[MAXDIV+1];
    double pi, u, v, a, b, xacc, x1, x2, epsilon, x;
    
	nsegments = MAXDIV;

	//nb = MAXDIV;

	//nb = NMAXROOT; nb is set to NMAXROOT within function roots()
            
	xacc = 1.0e-15;

    x1 = 0.0e0;
    x2 = 1.0e0;

    pi = acos(-1.0e0);
    	
    j = 0;
    mm = 0;
    
    double abs_epspjp, abs_epspjv;
    
    scatfN = 0.;
    
    while (j <= NMAXJ) {  // changed to macro constant NMAX // by Moon Jung, 2013/7/17
        
		j = j+1;
        roots(j, m, theta, x1, x2, xacc, nsegments, nb, rts);
        
		if (nb == 0)
            continue; // try the next interface j

        for (i = 1 ; i <= nb ; i ++)
        {
            thetai = 2.0 * atan( rts[i] );

            TransReflec(thetai, m, Trans, Reflec, Trp, Trv, Rep, Rev, tp, tv, rp, rv, t, r, u, v);
            path_length(radius, thetai, m, xi);
            abs_coef(lambda, m, beta);

            ScatFieldAmp(thetai, tp, tv, rp, rv, beta, xi, j, u, v, epspjp, epspjv);

            D = sin(2. * thetai) / ( 4.0 * sin(theta)) / abs( (float)(j-1) * cos(thetai) / sqrt(u*u+v*v) - 1.0);
            
            abs_epspjp =  abs(epspjp) * abs(epspjp);
            abs_epspjv =  abs(epspjv) * abs(epspjv);
            
            scatfN =  scatfN + (abs_epspjp + abs_epspjv) * D;   // eq (26b) in GOMsphere doc.
            
            mm = mm + 1;
            
        }//  endfor
        
    }//endwhile
    
    scatfN = pi * (radius * radius) * scatfN * sin(theta);  // eq (26b) in GOMsphere doc.  // scat function for near field
	        
    return;
    
} // scattering_function

void diffusion_function (double theta, dcomplex m, double lambda, double radius, double &scatfFB ) {
       
            
    double pi = acos(-1.0e0);
          
	
    double x = 2.0*pi*radius/lambda;  // size parameter of the sphere

    scatfFB = 2.0 * pi * (radius * radius) * BessJ1( x * sin(theta)) * BessJ1( x * sin(theta)) / sin(theta);  //eq (35) in GOMSphere; not use (35) itself
	                                                                                                          // because the integrand of (35) is integrated
	                                                                                                          // over [0, pi/2] rather than [0, pi]
	// scatfFB = scat function for far field diffusion (Bessel function)
    
   // asymNF = scatfN * cos(theta);   //  eq (37) in GOMsphere doc.
   // asymFFB = scatfFB * cos(theta);
    
    return;
    
} // diffusion_function


/* QscatAndG
 c     ************************************************************************************
 c     *********To calculate the scattering efficiency Qscat and asymmetry factor**********
 c     ************************************************************************************
 c PURPOSE:
 c     To calculate scattering efficiency and asymmetry factor for both near filed (QNscat)
 c     and far filed (QFscat), given wavelength, radius of sphere, number of roots for
 c     Gaussian quadratures. Eqs.(43)-(44) and (47-49).
 c REFERENCE:
 c     Davis P J, P Rainowitz, "Methods of numerical integration",Second edition,
 c     Orlando, FL: Academic Press,1984. p.481-483.
 c INPUT:
 c     lambda: wavelength
 c     radius: radius of sphere
 c     n: number of roots for Gaussian quadratures.
 c OUTPUT:
 c     QFscat: scattering efficiency for far field
 c     QNscat: scattering efficiency for near field
 c     gN: asymmetry factor for near field
 c     gF: asymmetry factor for far field
 */
double Fscat_crossSection(double lambda, double radius)
{
    //double gNQscat, gFQscat;
    double Fscat, Nscat, Dscat;

    double scatfN1, scatfN2, scatfFB1, scatfFB2;
    //double asymNF1, asymNF2, asymFFB1, asymFFB2;
    double temp, theta, in_angle, Trans, Reflec, Trp, Trv;
    dcomplex Rep, Rev, tp, tv, rp, rv, t, r, thetat;
    
	double beta, eps, re, im, pi, xi, u, v;
    dcomplex m;
    double scatF, xm, xk, dx, xmp, xkp, dxp;
    int mm, k, i, NN;
    double a, b, bp, xacc, x1, x2, epsilon;
    
    double D;
    int nbb, N2;
    
	double xb1[MAXDIV+1], xb2[MAXDIV+1], rts[MAXDIV+1]; // by Moon Jung, 2013.7.17
    
    int nsegments;
    
    int nb;
    
        
    double x[NUMGAUSS + 1], w[NUMGAUSS +1]; // by Moon Jung, 2013.7.17
    
    //nsegments = MAXDIV;
    //nb = MAXDIV;
    pi = acos(-1.0e0);
    
    a = 0.0e0;
    b = pi;

    bp = pi / 2.0;  //  int_{0}^{pi} P^{dif} approx = int_{0}^{pi/2} P^{diff} ?? 
        
	    
	const int CC = 30;

    double tempKK = CC + 273.0e0;
        
    refice (lambda, tempKK, m);

    GaussLeg(NUMGAUSS, x, w);
    
    
    xm = 0.5e0 * (b+a);
    xk = 0.5e0 * (b-a);

    xmp = 0.5e0 * (bp+a);
    xkp = 0.5e0 * (bp-a);
    mm = ( NUMGAUSS + 1) / 2;
    
	
    double FCscat;
	double NCscat =0.0;
	double DCscat = 0.0;
    double QNscat = 0.0;
    double QFscat = 0.0;
    
    
    scat_efficiency (lambda, radius, m, QNscat, QFscat);  // eq (43) and (47) in GOMsphere doc.
    
	/*
    for (int i = mm + 1; i <= NUMGAUSS; i++){  // by Moon Jung, 2013/7/17
        
        dx = xk * x[i];
        dxp = xkp * x[i];
      
		scattering_function (xm+dx, m, lambda, radius, scatfN1); // for the scattering function
        scattering_function (xm-dx, m, lambda, radius, scatfN2); 
        
        NCscat = NCscat + w[i] * (scatfN1 + scatfN2);
      
	//	gNQscat = gNQscat + w[i] * (asymNF1 + asymNF2);
        
        diffusion_function (xmp+dxp, m, lambda, radius, scatfFB1 );  // for the diffusion function
        diffusion_function (xmp-dxp, m, lambda, radius, scatfFB2 );
        
        DCscat = DCscat + w[i] * (scatfFB1 + scatfFB2);
        
	//	gFQscat = gFQscat + w[i] * (asymFFB1 + asymFFB2);
        
    }
    
    NCscat = NCscat * xk;
    //gNQscat = gNQscat * xk;
    
	//gN = gNQscat / Nscat;
    
    DCscat = DCscat * xkp;
    //gFQscat = gFQscat * xkp;
    
    //gF = gFQscat / Fscat;
    //gF = ((QNscat * gN) + gF) / (1. + QNscat);  // QFscat is not used
     
	FCscat = NCscat  +  DCscat;  // eq (22) in GOMsphere; Nscat = C_ref + C_tra, Dscat = C_dif
    
	*/

	FCscat = pi * (radius * radius) * QFscat;

	return FCscat;

    
} // Fscat_crossSection


double abs_crossSection (double lambda, double radius )
{
            
    double x[NUMGAUSS +1 ];
    double w[NUMGAUSS +1 ];

        
    double a, b; 
    int mm;
    double xm, xk, dx;
    
    const double pi = acos(-1.0e0);
    
    a = 0.0e0;
    b = 1.0e0;
    
    xm = .5e0 * (b + a);
    xk = .5e0 * (b - a);
    mm = (NUMGAUSS + 1) / 2;
    
    
	double Qabs = 0.0e0;

    GaussLeg(NUMGAUSS,x,w);
    for (int i = mm + 1 ; i <= NUMGAUSS ; i++){
        dx = xk * x[i];
        Qabs = Qabs + w[i] * ( absFF( xm+dx, lambda, radius ) + absFF( xm-dx, lambda, radius ) );
    }
    
	Qabs = Qabs * xk;
    
	double Cabs = pi * ( radius * radius) * Qabs;

    return Cabs;

} // abs_crossSection



/*    


// main program
int main(int argc, char *argv[]) {
    
 
    
//    for (int i = 1; i <= NWLT; i++ ) {
//        
//        tempRefraIndex[i].WLT = tempRefraIndex1[i].WLT;
//        tempRefraIndex[i].T[1].TABRET = tempRefraIndex1[i].T[0].TABRET;
//        tempRefraIndex[i].T[1].TABIMT = tempRefraIndex1[i].T[0].TABIMT;
//        tempRefraIndex[i].T[2].TABRET = tempRefraIndex1[i].T[1].TABRET;
//        tempRefraIndex[i].T[2].TABIMT = tempRefraIndex1[i].T[1].TABIMT;
//        tempRefraIndex[i].T[3].TABRET = tempRefraIndex2[i].T[0].TABRET;
//        tempRefraIndex[i].T[3].TABIMT = tempRefraIndex2[i].T[0].TABIMT;
//        tempRefraIndex[i].T[4].TABRET = tempRefraIndex2[i].T[1].TABRET;
//        tempRefraIndex[i].T[4].TABIMT = tempRefraIndex2[i].T[1].TABIMT;
//    }
//    
    // basic file operations
    
	// Compute the absorption cross section, scattering cross section of the spheres for plane waves with given wavelengths

	// sphere: radius = 0.05 mm ~ 2 mm   // sample 100 radii between the interval 
	// wavelength: 400 nm ( 400 * 1e-9 m ) ~ 700 nm ( 700 * 1e-9 m ); sample 30 wavelengths between the interval
    
	// Compute the phase function P (theta) for theta =0 ~ 180: sample 90 angles between the interval ( stepsize = 2 degree)


    std::ofstream task1_file, task2_file, task3_file;
    
    task1_file.open ("AbsCrossSection.dat");
    
    task2_file.open ("ScatCrossSection.dat");
    task3_file.open ("PhaseFunction.dat");
    
        
    task1_file << std::setw(15) << "radius[mm]" <<  std::setw(20) << "wavelength[nm]" << std::setw(15) << "Cabs" << std::endl;
    
    task2_file << std::setw(15) << "radius[mm]" <<  std::setw(20) << "wavelength[nm]" << std::setw(15) << "Csca" << std::endl;
    
    task3_file << std::setw(15) << "radisu[mm]" << std::setw(20) <<  "wavelength[nm]" <<  std::setw(15) << "scat angle" << std::setw(15) << "phase" << std::endl;
    
    
    double lambdaStart = 400e-9;   // 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
    double lambdaEnd = 700e-9;
    

//	const int nSpectralSamples = 30;
//    const int nRadii = 50;

	const int nSpectralSamples = 10;
    const int nRadii = 10;
    
    double radiusStart =  0.05e-3;  // 0.05 mm = 0.05 * e-3 m;
    double radiusEnd = 2e-3;           // 2 mm = 2 * e-3 m

    double stepLambda = ( lambdaEnd  - lambdaStart ) / double( nSpectralSamples) ;

	double stepRadius = (radiusEnd - radiusStart) /  double( nRadii );


    
	 // task 1 absorption cross section of each sphere for each wavelength

        
	task1_file << std::setiosflags (std::ios::scientific );

	double Cabs[ nRadii ] [nSpectralSamples];
	double radius, lambda;
    	double FCscat[nRadii][nSpectralSamples];

	radius = radiusStart;

    for ( int i =0; i < nRadii; i++ ) {
        
        lambda = lambdaStart;
		
        for (int j=0; j < nSpectralSamples; j ++ ) {
            
			Cabs[ i ][ j ] = abs_crossSection(lambda, radius);
            
            task1_file << std::setw(15) << std::setprecision(3) << radius * 1.0e+3 << std::setw(20)
            << std::setprecision(4) << lambda * 1.0e+9  << std::setw(15) << Cabs[ i ][ j ] << std::endl;
            
            lambda += stepLambda;
        }
        radius += stepRadius;
    }
    

      // task 2: compute the scattering cross section of each sphere for each wavelength
    
	task2_file << std::setiosflags (std::ios::scientific );


    
	radius = radiusStart;
    
    
    for (int i=0; i< nRadii; i++ ) {
        
        lambda = lambdaStart;
		for (int j =0; j < nSpectralSamples; j++ )  {
            
			FCscat[ i ][ j ] = Fscat_crossSection(lambda, radius );
            
            
            task2_file << std::setw(15) << std::setprecision(3) << radius * 1.0e+3 << std::setw(20)
            << std::setprecision(4) << lambda * 1.0e+9  << std::setw(15) << FCscat[ i ][ j ] << std::endl;
            
            lambda += stepLambda;
            
		}
		
        radius += stepRadius;
    }
    
    
	task3_file << std::setiosflags (std::ios::scientific );


//	const int nThetas = 60 ; // 3 degree in each step between theta 0 and theta 180;
    const int nThetas = 5 ; // 3 degree in each step between theta 0 and theta 180;

	double theta;
	 
	const  double pi = acos( -1.0e0);  // -1.0 const double in default; -1.0f = float

	double thetaStart = 3. * pi / 180.;
	double thetaEnd = pi;

	double stepTheta = ( thetaEnd - thetaStart ) / double (nThetas);

	double phase [nSpectralSamples] [nRadii] [nThetas];


	radius = radiusStart;


    
    for ( int i=0; i  < nRadii; i++ ) {
        
        lambda = lambdaStart;
        for (int j=0; j < nSpectralSamples ; j++ ) {
            
            theta = thetaStart;
            for (int k =0; k < nThetas; k++ ) {
                
		        phase[ i ][ j ][ k ] = phasefunction (lambda, radius, theta, FCscat[ i ][ j ] );
                
                task3_file << std::setw(15) << std::setprecision(3) << radius * 1.0e+3 << std::setw(20)
                << std::setprecision(4) << lambda * 1.0e+9  << std::setw(15) << theta * 180./pi <<
                std::setw(15) << phase[ i ][ j ] [ k ]  << std::endl;
                
			    theta += stepTheta;
            }

            lambda += stepLambda;

            
        }
        radius += stepRadius;

	}
    
} // main
*/