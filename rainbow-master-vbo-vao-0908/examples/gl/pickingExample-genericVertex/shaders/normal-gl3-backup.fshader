#version 130

#define MAX_ITERATIONS 100

out vec4 fragColor;

// opengl 3.0, 2008 => GLSL 1.3
// Opengl 3.1, 2009 => GLSL 1.3 (ATI and NVIDiA)
// Opengl 3.2, 2009 => GLSL 1.5 ( //): Compatibility Profile: Openframework uses this version

// Opengl 3.3, 2010 => GLSL 3.3 (NVIDIA)

 // As of March, 2010: Opengl 4.0, GLSL 4.0
 // My computer: GL 4.0 support

uniform  float uTextureFlag;

// uniform variables can be changed between draw commands

uniform mat4 uEyeMatrix;
uniform mat4 uInvEyeMatrix;

uniform mat4 uModelMatrixAABB;
uniform mat4 uModelViewMatrix;

uniform sampler2D uTexUnit0, vTexCoord0;
uniform sampler2D uTexNormal;

uniform sampler3D uDistTex;
uniform sampler3D uPhaseTex;
uniform sampler2D uScatTex;

// lights in eye space
uniform vec3 uLight1Pos;
uniform vec3 uLight2Pos;


// eyePos and aabb box for rainbow
//uniform vec3 uEyePos; // in eye space, so it should be (0,0,0): see below

uniform vec3 uLight1Pos, uLight2Pos;

uniform vec4 uLight1Color, uLight2Color;

uniform vec4 uMaterialColor;
uniform vec3 uLightRGBColor;

in vec2 vTexCoord0, vTexCoord1;
in vec3 vNormal;


uniform float uRadius;
uniform vec3  uSunRayDir;
uniform float uDropDensity;

// AABB bounding box
uniform vec3[2] uAABB; 
// equivalent to uniform vec3 uAABB[2]

/* The ways arrays are used: 
float[5] foo() { }
as a constructor of an array
float[5](3.4, 4.2, 5.0, 5.2, 1.1)

void foo(float[5])
float[5] a;

Arrays can have initializers formed from array constructors:
float a[5] = float[5](3.4, 4.2, 5.0, 5.2, 1.1);
float a[5] = float[](3.4, 4.2, 5.0, 5.2, 1.1); 
*/ 

in vec2 vTexCoord;
in mat3 vNTMat;

in vec3 vPosInEyeSpace;
 // varying geometry position in eye space



vec3 uEyePos = vec3(0,0,0);
float uTolerance = 1.0;

// parameters for rainbow computation

const  double pi = acos( -1.0e0);  // -1.0 const double in default; -1.0f = float

const int nSpectralSamples = 120;
const int nRadii = 2;
const int nThetas = 100; // 1 degree in each step between theta 130 deg and theta 142;

const double lambdaStart = 400;   // 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
const double lambdaEnd = 700;
     
const  double radiusStart =  0.05e-3;  // 0.05 mm = 0.05 * e-3 m;
const  double radiusEnd = 2e-3;           // 2 mm = 2 * e-3 m

const	double thetaStart = 130.0 / 180.0 * pi;
const	double thetaEnd = 142.0 / 180.0 *  pi;

const double lambdaStep = ( lambdaEnd - lambdaStart) / double( nSpectralSamples );
const double thetaStep = ( thetaEnd - thetaStart) / double (nThetas);
const double radiusStep = (radiusEnd - radiusStart) / double (nRadii);

const integer nLambdasForColorMatch = 81;
const integer nRGB = 3;


/* A colour system is defined by the CIE x and y coordinates of
   its three primary illuminants and the x and y coordinates of
   the white point. */

struct colourSystem {
    uint name[15];     	    	    // Colour system name 
    double xRed, yRed,	    	    // Red x, y
           xGreen, yGreen,  	    // Green x, y 
           xBlue, yBlue,    	    // Blue x, y 
           xWhite, yWhite,  	    // White point x, y 
	       gamma;       	        // Gamma correction for system
};


/* White point chromaticities. */


#define IlluminantC     0.3101, 0.3162	    	// For NTSC television 
#define IlluminantD65   0.3127, 0.3291	    	// For EBU and SMPTE
#define IlluminantE 	0.33333333, 0.33333333  // CIE equal-energy illuminant
#define GAMMA_REC709	0		                // Rec. 709


colourSystem
                  // Name                  xRed    yRed    xGreen  yGreen  xBlue  yBlue    White point        Gamma
    NTSCsystem  =  { "NTSC",               0.67,   0.33,   0.21,   0.71,   0.14,   0.08,   IlluminantC,    GAMMA_REC709 },
    EBUsystem   =  { "EBU (PAL/SECAM)",    0.64,   0.33,   0.29,   0.60,   0.15,   0.06,   IlluminantD65,  GAMMA_REC709 },
    SMPTEsystem =  { "SMPTE",              0.630,  0.340,  0.310,  0.595,  0.155,  0.070,  IlluminantD65,  GAMMA_REC709 },
    HDTVsystem  =  { "HDTV",               0.670,  0.330,  0.210,  0.710,  0.150,  0.060,  IlluminantD65,  GAMMA_REC709 },
    CIEsystem   =  { "CIE",                0.7355, 0.2645, 0.2658, 0.7243, 0.1669, 0.0085, IlluminantE,    GAMMA_REC709 },
    Rec709system = { "CIE REC 709",        0.64,   0.33,   0.30,   0.60,   0.15,   0.06,   IlluminantD65,  GAMMA_REC709 };

//colourSystem *cs = &SMPTEsystem;

colourSystem g_cs = SMPTEsystem;

double cie_colour_match[ nLambdasForColorMatch * nRGB ] = {
        {0.0014,0.0000,0.0065, 0.0022,0.0001,0.0105, 0.0042,0.0001,0.0201,
        0.0076,0.0002,0.0362, 0.0143,0.0004,0.0679, 0.0232,0.0006,0.1102,
        0.0435,0.0012,0.2074, 0.0776,0.0022,0.3713, 0.1344,0.0040,0.6456,
        0.2148,0.0073,1.0391, 0.2839,0.0116,1.3856, 0.3285,0.0168,1.6230,
        0.3483,0.0230,1.7471, 0.3481,0.0298,1.7826, 0.3362,0.0380,1.7721,
        0.3187,0.0480,1.7441, 0.2908,0.0600,1.6692, 0.2511,0.0739,1.5281,
        0.1954,0.0910,1.2876, 0.1421,0.1126,1.0419, 0.0956,0.1390,0.8130,
        0.0580,0.1693,0.6162, 0.0320,0.2080,0.4652, 0.0147,0.2586,0.3533,
        0.0049,0.3230,0.2720, 0.0024,0.4073,0.2123, 0.0093,0.5030,0.1582,
        0.0291,0.6082,0.1117, 0.0633,0.7100,0.0782, 0.1096,0.7932,0.0573,
        0.1655,0.8620,0.0422, 0.2257,0.9149,0.0298, 0.2904,0.9540,0.0203,
        0.3597,0.9803,0.0134, 0.4334,0.9950,0.0087, 0.5121,1.0000,0.0057,
        0.5945,0.9950,0.0039, 0.6784,0.9786,0.0027, 0.7621,0.9520,0.0021,
        0.8425,0.9154,0.0018, 0.9163,0.8700,0.0017, 0.9786,0.8163,0.0014,
        1.0263,0.7570,0.0011, 1.0567,0.6949,0.0010, 1.0622,0.6310,0.0008,
        1.0456,0.5668,0.0006, 1.0026,0.5030,0.0003, 0.9384,0.4412,0.0002,
        0.8544,0.3810,0.0002, 0.7514,0.3210,0.0001, 0.6424,0.2650,0.0000,
        0.5419,0.2170,0.0000, 0.4479,0.1750,0.0000, 0.3608,0.1382,0.0000,
        0.2835,0.1070,0.0000, 0.2187,0.0816,0.0000, 0.1649,0.0610,0.0000,
        0.1212,0.0446,0.0000, 0.0874,0.0320,0.0000, 0.0636,0.0232,0.0000,
        0.0468,0.0170,0.0000, 0.0329,0.0119,0.0000, 0.0227,0.0082,0.0000,
        0.0158,0.0057,0.0000, 0.0114,0.0041,0.0000, 0.0081,0.0029,0.0000,
        0.0058,0.0021,0.0000, 0.0041,0.0015,0.0000, 0.0029,0.0010,0.0000,
        0.0020,0.0007,0.0000, 0.0014,0.0005,0.0000, 0.0010,0.0004,0.0000,
        0.0007,0.0002,0.0000, 0.0005,0.0002,0.0000, 0.0003,0.0001,0.0000,
        0.0002,0.0001,0.0000, 0.0002,0.0001,0.0000, 0.0001,0.0000,0.0000,
        0.0001,0.0000,0.0000, 0.0001,0.0000,0.0000, 0.0000,0.0000,0.0000};
		

struct Ray {
  vec3 origin;
  vec3 direction;
  vec3 inv_direction;
  int sign[3];
};

