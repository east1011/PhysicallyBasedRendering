

#ifndef dropvolume_h
#define drovvolume_h


#include "cvec.h"

#include "matrix4.h"

#include "sphere.h"
#include "spa.h"  //include the SPA header file
//#include "optical_depth.h"
#include "irradiance_sun.h"
#include "phaseFunction.h"

#include "spectrum_to_rgb.h"

#include <ctype.h>  // from pbrt for readFloatFile
#include <stdlib.h>  // from pbrt for readFloatFile


using namespace std; // string type is defined in namespace 

#ifdef __linux__
// "Compiled for Linux
#else


// Windows doesn't define these values by default, Linux does
#define M_PI 3.141592653589793
#define INFINITY 1e8
#endif

#define ZERO_TOLERANCE = 1e-08



#define Assert(expr) \
    ((expr) ? (void)0 : \
        fprintf(stderr, "Assertion \"%s\" failed in %s, line %d", \
               #expr, __FILE__, __LINE__))



/* set bit masks for the possible character types */

#define _UPPER          0x1     /* upper case letter */ 
#define _LOWER          0x2     /* lower case letter */
#define _DIGIT          0x4     /* digit[0-9] */
#define _SPACE          0x8     /* tab, carriage return, newline, */
                                /* vertical tab or form feed */
#define _PUNCT          0x10    /* punctuation character */
#define _CONTROL        0x20    /* control character */
#define _BLANK          0x40    /* space char = 32 */
#define _HEX            0x80    /* hexadecimal digit */

#define _LEADBYTE       0x8000                  /* multibyte leadbyte */
#define _ALPHA          (0x0100|_UPPER|_LOWER)  /* alphabetic character */


// externals whose contents are defined in dropvolume.cpp
// extern declares a variable without defining it; In the case of a function, extern is not 
// needed because by default functions are external. Do prevent this default behavior, functions
// should be declared as static. Variables can be externally declared and defined in the same
// file.

extern double Cabs[ nRadii ] [nSpectralSamples];
extern double FCscat[nRadii][nSpectralSamples];


extern double particlePhase [nRadii] [nSpectralSamples] [nThetas];
extern double  radii[ nRadii ];
extern double  spectralSamples [ nSpectralSamples ];
extern double  thetas[ nThetas ];


extern Cvec3  g_sunRay;
extern Cvec3  g_light_xyzColor, g_light_rgbColor;
  
void setupVolume();


#define Max(a, b)   (((a) > (b)) ? (a) : (b))


void setSunLightRGBColor();

Cvec3 calculate_rainbowColor (double radius, double uDropDensity, double phi_sun, double phi_obs, 
							   double psi, double tau_N);
Cvec3 rayIntersectVolume(   double uDropDensity, double radius,
								  Cvec3 rayorig, Cvec3 raydir, Cvec3 sunRay);

void  readScatCrossFile(const char *filename, int nRadii, int nSpectralSamples );
void  readPhaseFile(const char *filename, int nRadii, int nSpectralSamples, int nThetas );

void writeScatCrossFile( );
void writePhaseFile( );
void writeRadiusFile();
void writeSpectrumFile();
void writeThetaFile();




#endif
