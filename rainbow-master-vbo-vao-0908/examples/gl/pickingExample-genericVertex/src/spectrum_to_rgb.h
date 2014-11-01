#ifndef spectrum_to_rgb_h
#define spectrum_to_rgb_h

#include <stdio.h>
#include <math.h>
#include "cvec.h"

// CODE FROM https://www.fourmilab.ch/documents/specrend/
/* A colour system is defined by the CIE x and y coordinates of
   its three primary illuminants and the x and y coordinates of
   the white point. */

struct colourSystem {
    char *name;     	    	    /* Colour system name */
    double xRed, yRed,	    	    /* Red x, y */
           xGreen, yGreen,  	    /* Green x, y */
           xBlue, yBlue,    	    /* Blue x, y */
           xWhite, yWhite,  	    /* White point x, y */
	   gamma;   	    	    /* Gamma correction for system */
};

/* White point chromaticities. */


#define IlluminantC     0.3101, 0.3162	    	/* For NTSC television */
#define IlluminantD65   0.3127, 0.3291	    	/* For EBU and SMPTE */
#define IlluminantE 	0.33333333, 0.33333333  /* CIE equal-energy illuminant */

/*  Gamma of nonlinear correction.

    See Charles Poynton's ColorFAQ Item 45 and GammaFAQ Item 6 at:
    
       http://www.poynton.com/ColorFAQ.html
       http://www.poynton.com/GammaFAQ.html
 
*/

#define GAMMA_REC709	0		/* Rec. 709 */

static struct colourSystem
                  /* Name                  xRed    yRed    xGreen  yGreen  xBlue  yBlue    White point        Gamma   */
    NTSCsystem  =  { "NTSC",               0.67,   0.33,   0.21,   0.71,   0.14,   0.08,   IlluminantC,    GAMMA_REC709 },
    EBUsystem   =  { "EBU (PAL/SECAM)",    0.64,   0.33,   0.29,   0.60,   0.15,   0.06,   IlluminantD65,  GAMMA_REC709 },
    SMPTEsystem =  { "SMPTE",              0.630,  0.340,  0.310,  0.595,  0.155,  0.070,  IlluminantD65,  GAMMA_REC709 },
    HDTVsystem  =  { "HDTV",               0.670,  0.330,  0.210,  0.710,  0.150,  0.060,  IlluminantD65,  GAMMA_REC709 },
    CIEsystem   =  { "CIE",                0.7355, 0.2645, 0.2658, 0.7243, 0.1669, 0.0085, IlluminantE,    GAMMA_REC709 },
    Rec709system = { "CIE REC 709",        0.64,   0.33,   0.30,   0.60,   0.15,   0.06,   IlluminantD65,  GAMMA_REC709 };



/*  	    	    	    UPVP_TO_XY

    Given 1976 coordinates u', v', determine 1931 chromaticities x, y
    
*/

void upvp_to_xy(double up, double vp, double *xc, double *yc);

/*  	    	    	    XY_TO_UPVP

    Given 1931 chromaticities x, y, determine 1976 coordinates u', v'
    
*/


void xy_to_upvp(double xc, double yc, double *up, double *vp);

/*                            INSIDE_GAMUT

     Test whether a requested colour is within the gamut
     achievable with the primaries of the current colour
     system.  This amounts simply to testing whether all the
     primary weights are non-negative. */


int inside_gamut(double r, double g, double b);

/*                          CONSTRAIN_RGB

    If the requested RGB shade contains a negative weight for
    one of the primaries, it lies outside the colour gamut 
    accessible from the given triple of primaries.  Desaturate
    it by adding white, equal quantities of R, G, and B, enough
    to make RGB all positive.  The function returns 1 if the
    components were modified, zero otherwise.
    
*/


int constrain_rgb(double *r, double *g, double *b);

/*                          GAMMA_CORRECT_RGB

    Transform linear RGB values to nonlinear RGB values. Rec.
    709 is ITU-R Recommendation BT. 709 (1990) ``Basic
    Parameter Values for the HDTV Standard for the Studio and
    for International Programme Exchange'', formerly CCIR Rec.
    709. For details see
    
       http://www.poynton.com/ColorFAQ.html
       http://www.poynton.com/GammaFAQ.html
*/


void gamma_correct(const struct colourSystem *cs, double *c);
void gamma_correct_rgb(const struct colourSystem *cs, double *r, double *g, double *b);

/*  	    	    	    NORM_RGB

    Normalise RGB components so the most intense (unless all
    are zero) has a value of 1.
    
*/


void norm_rgb(double *r, double *g, double *b);

/*                          SPECTRUM_TO_XYZ

    Calculate the CIE X, Y, and Z coordinates corresponding to
    a light source with spectral distribution given by  the
    function SPEC_INTENS, which is called with a series of
    wavelengths between 380 and 780 nm (the argument is 
    expressed in meters), which returns emittance at  that
    wavelength in arbitrary units.  The chromaticity
    coordinates of the spectrum are returned in the x, y, and z
    arguments which respect the identity:

            x + y + z = 1.
*/

void norm_rgb(Cvec3 &rgb);
int constrain_rgb(Cvec3 &rgb);

void spectrum_to_xyz(double (*spec_intens)(double wavelength),
                     double *x, double *y, double *z);

/*                            BB_SPECTRUM

    Calculate, by Planck's radiation law, the emittance of a black body
    of temperature bbTemp at the given wavelength (in metres).  */


extern double bbTemp;                 /* Hidden temperature argument
                                         to BB_SPECTRUM. */


/*                             XYZ_TO_RGB

    Given an additive tricolour system CS, defined by the CIE x
    and y chromaticities of its three primaries (z is derived
    trivially as 1-(x+y)), and a desired chromaticity (XC, YC,
    ZC) in CIE space, determine the contribution of each
    primary in a linear combination which sums to the desired
    chromaticity.  If the  requested chromaticity falls outside
    the Maxwell  triangle (colour gamut) formed by the three
    primaries, one of the r, g, or b weights will be negative. 

    Caller can use constrain_rgb() to desaturate an
    outside-gamut colour to the closest representation within
    the available gamut and/or norm_rgb to normalise the RGB
    components so the largest nonzero component has value 1.
    
*/


Cvec3  xyz_to_rgb(struct colourSystem *cs, const Cvec3 xyz);

//static struct colourSystem *cs = &SMPTEsystem;
static struct colourSystem *cs = &Rec709system;

// if you want to define a variable in .h file, you had better define it as static,
//  so that even if the header file is included into several files, there would be 
//  no multiple definition error.

/* CIE colour matching functions xBar, yBar, and zBar for
wavelengths from 380 through 780 nanometers, every 5
nanometers.  For a wavelength lambda in this range:

cie_colour_match[(lambda - 380) / 5][0] = xBar
cie_colour_match[(lambda - 380) / 5][1] = yBar
cie_colour_match[(lambda - 380) / 5][2] = zBar

To save memory, this table can be declared as floats
rather than doubles; (IEEE) float has enough 
significant bits to represent the values. It's declared
as a double here to avoid warnings about "conversion
between floating-point types" from certain persnickety
compilers. */

static double cie_colour_match[81][3] = {
	{0.0014,0.0000,0.0065}, {0.0022,0.0001,0.0105}, {0.0042,0.0001,0.0201}, {0.0076,0.0002,0.0362}, 
	// 400
	{0.0143,0.0004,0.0679}, {0.0232,0.0006,0.1102}, {0.0435,0.0012,0.2074}, {0.0776,0.0022,0.3713}, {0.1344,0.0040,0.6456}, 
	// 425
	{0.2148,0.0073,1.0391}, {0.2839,0.0116,1.3856}, {0.3285,0.0168,1.6230},	{0.3483,0.0230,1.7471}, {0.3481,0.0298,1.7826},
	// 450
	{0.3362,0.0380,1.7721}, {0.3187,0.0480,1.7441}, {0.2908,0.0600,1.6692}, {0.2511,0.0739,1.5281}, {0.1954,0.0910,1.2876},
	// 475
	{0.1421,0.1126,1.0419}, {0.0956,0.1390,0.8130}, {0.0580,0.1693,0.6162}, {0.0320,0.2080,0.4652}, {0.0147,0.2586,0.3533},
	// 500
	{0.0049,0.3230,0.2720}, {0.0024,0.4073,0.2123}, {0.0093,0.5030,0.1582}, {0.0291,0.6082,0.1117}, {0.0633,0.7100,0.0782},
	// 525
	{0.1096,0.7932,0.0573}, {0.1655,0.8620,0.0422}, {0.2257,0.9149,0.0298}, {0.2904,0.9540,0.0203}, {0.3597,0.9803,0.0134}, 
	// 550
	{0.4334,0.9950,0.0087}, {0.5121,1.0000,0.0057}, {0.5945,0.9950,0.0039}, {0.6784,0.9786,0.0027}, {0.7621,0.9520,0.0021},
	// 575
	{0.8425,0.9154,0.0018}, {0.9163,0.8700,0.0017}, {0.9786,0.8163,0.0014}, {1.0263,0.7570,0.0011}, {1.0567,0.6949,0.0010},
	// 600
	{1.0622,0.6310,0.0008}, {1.0456,0.5668,0.0006}, {1.0026,0.5030,0.0003}, {0.9384,0.4412,0.0002}, {0.8544,0.3810,0.0002},
	// 625
	{0.7514,0.3210,0.0001}, {0.6424,0.2650,0.0000}, {0.5419,0.2170,0.0000}, {0.4479,0.1750,0.0000}, {0.3608,0.1382,0.0000},
	// 650
	{0.2835,0.1070,0.0000}, {0.2187,0.0816,0.0000}, {0.1649,0.0610,0.0000}, {0.1212,0.0446,0.0000}, {0.0874,0.0320,0.0000}, 
	// 675
	{0.0636,0.0232,0.0000}, {0.0468,0.0170,0.0000}, {0.0329,0.0119,0.0000}, {0.0227,0.0082,0.0000}, {0.0158,0.0057,0.0000}, 
	// 700
	{0.0114,0.0041,0.0000}, {0.0081,0.0029,0.0000}, {0.0058,0.0021,0.0000}, {0.0041,0.0015,0.0000}, {0.0029,0.0010,0.0000},
	// 725
	{0.0020,0.0007,0.0000}, {0.0014,0.0005,0.0000}, {0.0010,0.0004,0.0000}, {0.0007,0.0002,0.0000}, {0.0005,0.0002,0.0000}, 
	// 750
	{0.0003,0.0001,0.0000}, {0.0002,0.0001,0.0000}, {0.0002,0.0001,0.0000}, {0.0001,0.0000,0.0000}, {0.0001,0.0000,0.0000},
	// 775
	{0.0001,0.0000,0.0000}, {0.0000,0.0000,0.0000}
};

// Raymond Lee's sun radiance data



float ISpectrumSun  [61] = {
//wavelength (nm)	normalized irradiance (arbitrary units)
// 700 - 400 = 300 / 5 = 60 + 1 = nSpectralSamples= 61 
// The original wavelength starts from 380 to 700 with step 5 nm
//380
//	0.4651826, 	0.3830807, 	0.4409676, 	0.4415811, 	
//400
	0.5531539, 	0.7498946, 0.746196, 0.7945927,
//420
	0.8004899, 	0.7938542, 	0.7373845, 	0.7613856,   
	0.8335401,  0.9062264,   0.9455949, 0.9749124,
//460
	0.9885056, 	0.9996937, 	0.9703307, 	0.9874594, 	1, 	
	0.9923742, 	0.9242011, 	0.9640443, 	0.9617534,
//505
	0.9306258, 	0.9463357, 	0.9275853, 	0.8597925, 	
	0.9149771, 	0.9065616, 	0.9197563, 	0.9100043,
//545
	0.8943484, 0.8951805, 0.8915439, 0.8640333,  
	0.8611063,   0.8411431,  	0.8347283,  0.82272,
// 585	
   0.8433171, 	0.7949293,  	0.7905713,  0.7947434, 
   0.7932764,  0.7953064, 	0.769975, 	0.7614744,
//625
	0.7494112,  0.7195148,  	0.7222514,  0.7319365,
	0.7231518,  0.6976414,  0.6917001,  0.6692451,
//665
	0.7032121,  0.6983164,  	0.6871358,  0.6830449, 
	0.6650319, 	0.5508773, 	0.5994822, 	0.6185095
};


#endif


