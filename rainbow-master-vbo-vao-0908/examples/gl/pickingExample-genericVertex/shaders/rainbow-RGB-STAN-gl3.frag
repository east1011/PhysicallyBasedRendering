#version 130

#define MAX_ITERATIONS 100

uniform int uWindowWidth; 
uniform int uWindowHeight; 

// opengl 3.0, 2008 => GLSL 1.3
// Opengl 3.1, 2009 => GLSL 1.3 (ATI and NVIDiA)
// ***Opengl 3.2, 2009 => GLSL 1.5 ( //): Compatibility nSpectralSamples: Openframework uses this version

// Opengl 3.3, 2010 => GLSL 3.3 (NVIDIA)


uniform mat4 uEyeMatrix;
uniform mat4 uInvEyeMatrix;

uniform mat4 uModelMatrixAABB;
uniform mat4 uInvModelMatrixAABB;
uniform mat4 uModelViewMatrix;
uniform mat4 uProjectionMatrix;

uniform sampler2D uColorTex;
uniform sampler2D uDepthTex;

uniform sampler3D uDistTex;
uniform sampler3D uPhaseTex;
uniform sampler2D uScatTex;

uniform sampler1D uPrimarySpectrums;
uniform sampler1D uXYZSpectrums;

// lights in eye space
uniform vec3 uLight1Pos;
uniform vec3 uLight2Pos;

uniform vec4 uLight1Color, uLight2Color;

uniform vec4 uMaterialColor;
uniform vec3 uLightRGBColor;

// eyePos and aabb box for rainbow
uniform vec3 uEyeOriginInBoxFrame;
 

uniform float uRadius;
uniform vec3  uSunRayDir;
uniform float uDropDensity;

// AABB bounding box
uniform vec3 uAABBmin; 
uniform vec3 uAABBmax;

// the normal to the front face of AABB box in eye space

uniform vec3 uFrontNormalToAABB;

in vec2 vTexCoord;
in vec3 vNormal;
in vec3 vPosition;
in mat3 vNTMat;

out vec4 fragColor;
out vec4 spect1;

// parameters for rainbow computation

float pi = acos(-1.0);  // -1.0 const float in default; -1.0f = float

//const int nSpectralSampleSteps = 120;
const int nThetaSteps = 100; // 100 steps between theta 130 deg and theta 142;

const int nRadiusSteps  = 2;

//const int nSpectralSamples = nSpectralSampleSteps +1;


const int nSpectralSamples = 61; // from 400 nm to 700nm
const int nSamples = nSpectralSamples; // from 400 nm to 700nm

const int nRadii = nRadiusSteps +1;
const int nThetas = nThetaSteps +1; // 100 steps between theta 130 deg and theta 142;

//uniform float uPhaseFunction[nRadii * nSpectralSamples * nThetas];

//float irainbow[ nSpectralSamples ];
float particlePhase[nSpectralSamples ];

// uniform block for primary spectrums

const float CIE_Y_integral = 106.856895;
const float lambdaStart = 400;		// 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
const float lambdaEnd = 700;
     
const float radiusStart =  1.0e-3;  
const float radiusEnd = 2.0e-3;		// 2 mm = 2 * e-3 m

const float thetaStart = 130.0;
const float thetaEnd = 142.0;

const float lambdaStep = ( lambdaEnd - lambdaStart) / float( nSpectralSamples -1); //lambdaStep = (700-400)/(61-1)
const float thetaStep = ( thetaEnd - thetaStart) / float (nThetas -1);
const float radiusStep = (radiusEnd - radiusStart) / float (nRadii -1);

float h0 = 2;

struct Ray {
	vec3 origin;
	vec3 direction;
	vec3 inv_direction;
	int sign[3];
};




Ray makeRay( vec3 origin, vec3 direction) {

	vec3 inv_direction = vec3(1.0) / direction;

	Ray r;
	r.origin = origin;
	r.direction = direction;
	r.inv_direction = inv_direction;
	r.sign[0] = (inv_direction.x < 0) ? 1 : 0; 
	r.sign[1] =  (inv_direction.y < 0)  ?  1:  0;
	r.sign[2] =  (inv_direction.z < 0)  ?  1:  0;
 
	return r;

}


// OLD RAINBOW
const int nLambdasForColorMatch = 81;
const int nRGB = 3;


/* A colour system is defined by the CIE x and y coordinates of
   its three primary illuminants and the x and y coordinates of
   the white point. */

struct colourSystem {
	//uint name[15];     	    /* Colour system name */
	float xRed, yRed,	    	/* Red x, y */
		xGreen, yGreen,  	    /* Green x, y */
		xBlue, yBlue,    	    /* Blue x, y */
		xWhite, yWhite,  	    /* White point x, y */
		gamma;   	    		/* Gamma correction for system */
};


/* White point chromaticities. */


#define IlluminantCxWhite   0.3101	    /* For NTSC television */
#define IlluminantCyWhite	0.3162
#define IlluminantD65xWhite 0.3127		/* For EBU and SMPTE */
#define IlluminantD65yWhite 0.3291	
#define IlluminantExWhite 	0.33333333	/* CIE equal-energy illuminant */
#define IlluminantEyWhite	0.33333333 
#define GAMMA_REC709		0	        /* Rec. 709 */

/* Name         xRed    yRed    xGreen  yGreen  xBlue  yBlue    White point        Gamma   */
/*
colourSystem                   
    NTSCsystem  =  { "NTSC",               0.67,   0.33,   0.21,   0.71,   0.14,   0.08,   IlluminantC,    GAMMA_REC709 },
    EBUsystem   =  { "EBU (PAL/SECAM)",    0.64,   0.33,   0.29,   0.60,   0.15,   0.06,   IlluminantD65,  GAMMA_REC709 },
    SMPTEsystem =  { "SMPTE",              0.630,  0.340,  0.310,  0.595,  0.155,  0.070,  IlluminantD65,  GAMMA_REC709 },
    HDTVsystem  =  { "HDTV",               0.670,  0.330,  0.210,  0.710,  0.150,  0.060,  IlluminantD65,  GAMMA_REC709 },
    CIEsystem   =  { "CIE",                0.7355, 0.2645, 0.2658, 0.7243, 0.1669, 0.0085, IlluminantE,    GAMMA_REC709 },
    Rec709system = { "CIE REC 709",        0.64,   0.33,   0.30,   0.60,   0.15,   0.06,   IlluminantD65,  GAMMA_REC709 };
*/

//colourSystem *cs = &SMPTEsystem;

//colourSystem g_cs = colourSystem(0.67, 0.33, 0.21, 0.71, 0.14, 0.08, IlluminantCxWhite, IlluminantCyWhite, GAMMA_REC709);
colourSystem g_cs = colourSystem(0.64, 0.33, 0.30, 0.60, 0.15, 0.06, IlluminantD65xWhite, IlluminantD65yWhite, GAMMA_REC709);

float cie_colour_match[ nLambdasForColorMatch * nRGB ] = 
	float[] ( 
	    0.0014,0.0000,0.0065, 0.0022,0.0001,0.0105, 0.0042,0.0001,0.0201,
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
        0.0001,0.0000,0.0000, 0.0001,0.0000,0.0000, 0.0000,0.0000,0.0000);
		
float [nSpectralSamples] ISpectrumSun = float[] (
//wavelength (nm)	normalized irradiance (arbitrary units)
// 700 - 400 = 300 / 5 = 60 + 1 = nSpectralSamples= 61 
// The original wavelength starts from 380 to 700 with step 5 nm
//380
//	0.4651826, 	0.3830807, 	0.4409676, 	0.4415811, 	
//400
	0.5531539, 	0.7498946, 0.746196, 0.7945927,
//420
	0.8004899, 	0.7938542, 	0.7373845, 	0.7613856,   0.8335401,  0.9062264,   0.9455949, 0.9749124,
//460
	0.9885056, 	0.9996937, 	0.9703307, 	0.9874594, 	1, 	0.9923742, 	0.9242011, 	0.9640443, 	0.9617534,
//505
	0.9306258, 	0.9463357, 	0.9275853, 	0.8597925, 	0.9149771, 	0.9065616, 	0.9197563, 	0.9100043,
//545
	0.8943484, 0.8951805, 0.8915439, 0.8640333,  	0.8611063,   0.8411431,  	0.8347283,  0.82272,
// 585	
   0.8433171, 	0.7949293,  	0.7905713,  0.7947434,  0.7932764,  0.7953064, 	0.769975, 	0.7614744,
//625
	0.7494112,  0.7195148,  	0.7222514,  0.7319365,  0.7231518,  0.6976414,  0.6917001,  0.6692451,
//665
	0.7032121,  0.6983164,  	0.6871358,  0.6830449, 	0.6650319, 	0.5508773, 	0.5994822, 	0.6185095
	);


bool approximatelyEqual(float a, float b, float epsilon)
{  // compare the difference to the larger of the two numbers
// If the diff is smaller than n% of the larger, both can be considered equal 

	return abs(a - b) <= ( ( abs(a) < abs(b) ? abs(b) : abs(a)  ) * epsilon);
}

bool essentiallyEqual(float a, float b, float epsilon)
{   // check if diff is smaller than n% of the smaller 
    return abs(a - b) <= ( (abs(a) > abs(b) ? abs(b) : abs(a) ) * epsilon );
}

 bool definitelyGreaterThan(float a, float b, float epsilon)
{
    return (a - b) > ( (abs(a) < abs(b) ? abs(b) : abs(a) ) * epsilon );
}

 bool definitelyLessThan(float a, float b, float epsilon)
{
    return (b - a) > ( (abs(a) < abs(b) ? abs(b) : abs(a)) * epsilon);
}


vec3 posInAABB (vec3 posInEyeSpace) {

// aabb = AABB in eye space

//posInEyeSpace = uAABBmin + posInAABB
 
	vec3 widths = uAABBmax - uAABBmin; 

	vec3 posInAABB  = posInEyeSpace - uAABBmin;
	
	vec3 texCoord = posInAABB / widths;

	return texCoord;
	
}


    
int intersection_distance_no_if ( 
   in Ray ray,   in vec3[2] aabb,
   out float tmin, out float tmax,  out vec3 normalAtFront ) {
    
		tmin = ( aabb[ ray.sign[0] ].x - ray.origin.x ) * ray.inv_direction.x;
		tmax = ( aabb[1- ray.sign[0] ].x - ray.origin.x ) * ray.inv_direction.x;
		float tymin = ( aabb[ ray.sign[1]  ].y - ray.origin.y ) * ray.inv_direction.y;
		float tymax = ( aabb[ 1 - ray.sign[1]  ].y - ray.origin.y ) * ray.inv_direction.y;
		float tzmin = ( aabb[ ray.sign[2] ].z - ray.origin.z ) * ray.inv_direction.z;
		float tzmax = ( aabb[ 1- ray.sign[2] ].z - ray.origin.z ) * ray.inv_direction.z;
		tmin = max( max( tmin, tymin), tzmin );
		tmax = min( min( tmax, tymax), tzmax );

   
		// check if the ray intersects the AABB

		int hitType;

		if ( tmin >  tmax ) {
				hitType = 0;  // ray does not hit AABB box
	 
				return hitType;

		}
		// compute the normal at the front face relative to the eye, where
		// the hit position is at tmin along the ray direction.
		// Assume that the front face is hit for now  	  
    
		normalAtFront =  vec3( 0,0,1); 

		hitType =1;
		return hitType;


	
}

    
vec3 XYZToRGB(vec3 xyz);


//Cvec3  xyz_to_rgb(struct colourSystem *cs,   const vec3 xyz) {

vec3  xyz_to_rgb( colourSystem cs, const vec3 xyz ) {

    float xr, yr, zr, xg, yg, zg, xb, yb, zb;
    float xw, yw, zw;
    float rx, ry, rz, gx, gy, gz, bx, by, bz;
    float rw, gw, bw;

	xr = cs.xRed;    yr = cs.yRed;    zr = 1 - (xr + yr);
    xg = cs.xGreen;  yg = cs.yGreen;  zg = 1 - (xg + yg);
    xb = cs.xBlue;   yb = cs.yBlue;   zb = 1 - (xb + yb);

    xw = cs.xWhite;  yw = cs.yWhite;  zw = 1 - (xw + yw);

    // xyz -> rgb matrix, before scaling to white. 
    
    rx = (yg * zb) - (yb * zg);  ry = (xb * zg) - (xg * zb);  rz = (xg * yb) - (xb * yg);
    gx = (yb * zr) - (yr * zb);  gy = (xr * zb) - (xb * zr);  gz = (xb * yr) - (xr * yb);
    bx = (yr * zg) - (yg * zr);  by = (xg * zr) - (xr * zg);  bz = (xr * yg) - (xg * yr);

    // White scaling factors.
    //  Dividing by yw scales the white luminance to unity, as conventional. 
       
    rw = ((rx * xw) + (ry * yw) + (rz * zw)) / yw;
    gw = ((gx * xw) + (gy * yw) + (gz * zw)) / yw;
    bw = ((bx * xw) + (by * yw) + (bz * zw)) / yw;

    // xyz -> rgb matrix, correctly scaled to white. 
    
    rx = rx / rw;  ry = ry / rw;  rz = rz / rw;
    gx = gx / gw;  gy = gy / gw;  gz = gz / gw;
    bx = bx / bw;  by = by / bw;  bz = bz / bw;

    // rgb of the desired point 
    
	vec3 rgbColor = vec3((rx * xyz[0]) + (ry * xyz[1]) + (rz * xyz[2]), 
					(gx * xyz[0]) + (gy * xyz[1]) + (gz * xyz[2]), 
					(bx * xyz[0]) + (by * xyz[1]) + (bz * xyz[2]));

	return rgbColor;

	
} // xyz_to_rgb



bool constrain_rgb(inout vec3 rgb) {
    float w;

    // Amount of white needed is w = - min(0, *r, *g, *b) 
	    
	if (0.0 < rgb[0]) w = 0.0; else w = rgb[0];
	if  (w < rgb[1])  w =w; else w =  rgb[1];
	if  (w < rgb[2])  w =w; else w = rgb[2];
    w = -w;

    // Add just enough white to make r, g, b all positive. 
    
    if (w > 0.0) {
		rgb += w;
        return true;                     // Colour modified to fit RGB gamut 
    }

    return false;                         // Colour within RGB gamut 
}

void normalize_positive_vector( inout float[nSpectralSamples] vector, in  int size) {

    float max = -1.0;

	for ( int i = 0; i < size ; i++ ) {
		if ( vector[i] > max ) max = vector[i];
	}

	if ( max > 0) {   
	  
		for ( int i = 0; i < size; i++ ) {
			vector[i] /= max;
		}
		  
	}
}

void norm_rgb( inout vec3 rgb) {
	#define Max(a, b)   (((a) > (b)) ? (a) : (b))
	float greatest = Max(rgb[0], Max(rgb[1], rgb[2]));
    
    if (greatest > 0) {
		rgb[0] /= greatest;
		rgb[1] /= greatest;
		rgb[2] /= greatest;
    }
	#undef Max
}

vec3 oldToRGB( float [nSpectralSamples] irainbow ) {

	// convert spectral color to rgb color
	
	float X = 0;
	float Y = 0;
	float Z = 0;
	float XYZ;

	float xBar, yBar, zBar;

	// cie_colour_match[][0] = xBar[]
    // cie_colour_match[][1] = yBar[]
    // cie_colour_match[][2] = zBar[]
	//  ==> cie_coulor_match [] from 380 nanometer to 780 nanometer
    
	
	// lambda: 400 ~ 700 into 60 steps (61 samples) with lamdbaStep - 5nm
	// cie_colour_match: from 380 to 780 nm with 5 nm step
	// Sun spectrum: from 380 to 700 ( 65 samples) with step = 5nm
	float yInt = 0.0;

	//float lambdaStep = (lambdaEnd - lambdaStart)/(float)(nSpectralSamples - 1);

	for (int j= 0; j < nSpectralSamples; j++) {
		
		int k = j + 4;

		xBar =  cie_colour_match[ k * 3 + 0 ]; 
		yBar =  cie_colour_match[ k * 3 + 1 ];
		zBar =  cie_colour_match[ k * 3 + 2 ];

		yInt += yBar;

		X += irainbow[j] * xBar * lambdaStep; 
		Y += irainbow[j] * yBar * lambdaStep;
		Z += irainbow[j] * zBar * lambdaStep;
		
    }
	
	XYZ = (X+Y+Z);
	//vec3 xyzColor = vec3( X/yInt, Y/yInt, Z/yInt );
	vec3 xyzColor = vec3( X/XYZ, Y/XYZ, Z/XYZ );
	vec3 rgbColor = xyz_to_rgb(g_cs, xyzColor);

	//vec3 rgbColor = XYZToRGB(xyzColor);
	/*
	if (constrain_rgb(rgbColor)) {  // Colour modified to fit RGB gamut 
		norm_rgb(rgbColor);
		//output << rgbColor << "(Approximation)" << endl;
	} else {
		norm_rgb(rgbColor);
		//output << rgbColor << endl;
	}
	*/

	//constrain_rgb(rgbColor);

	return rgbColor;

} // oldToRGB()




// NEW RAINBOW

float[nSpectralSamples] X;
float[nSpectralSamples] Y;
float[nSpectralSamples] Z;

float[nSpectralSamples] rgbRefl2SpectWhite;
float[nSpectralSamples] rgbRefl2SpectCyan;
float[nSpectralSamples] rgbRefl2SpectMagenta;
float[nSpectralSamples] rgbRefl2SpectYellow;
float[nSpectralSamples] rgbRefl2SpectRed;
float[nSpectralSamples] rgbRefl2SpectGreen;
float[nSpectralSamples] rgbRefl2SpectBlue;


float[nSpectralSamples] XSpectrum() {
	float[nSpectralSamples]  r;
	for (int j=0; j < nSpectralSamples; j++ ) {
		r[j] = texelFetch(uXYZSpectrums, 0 * nSpectralSamples + j, 0).r; 
	}
	return r;
}

float[nSpectralSamples] YSpectrum() {
	float[nSpectralSamples]  r;
	for (int j=0; j < nSpectralSamples; j++ ) {
		r[j] = texelFetch(uXYZSpectrums, 1 * nSpectralSamples + j, 0).r; 
	}
	return r;
}

float[nSpectralSamples] ZSpectrum() {
	float[nSpectralSamples]  r;
	for (int j=0; j < nSpectralSamples; j++ ) {
		r[j] = texelFetch(uXYZSpectrums, 2 * nSpectralSamples + j, 0).r; 
	}
	return r;
}
float[nSpectralSamples] whiteSpectrum() {
	float[nSpectralSamples]  r;
	for (int j=0; j < nSpectralSamples; j++ ) {
		r[j] = texelFetch(uPrimarySpectrums, 0 * nSpectralSamples + j, 0).r; 
	}
	return r;
}
float[nSpectralSamples] cyanSpectrum() {
	float[nSpectralSamples]  r;
	for (int j=0; j < nSpectralSamples; j++ ) {
		r[j] = texelFetch(uPrimarySpectrums, 1 * nSpectralSamples + j, 0).r; 
	}
	return r;
}
float[nSpectralSamples] magentaSpectrum() {
	float[nSpectralSamples]  r;
	for (int j=0; j < nSpectralSamples; j++ ) {
		r[j] = texelFetch(uPrimarySpectrums, 2 * nSpectralSamples + j, 0).r; 
	}
	return r;
}
float[nSpectralSamples] yellowSpectrum() {
	float[nSpectralSamples]  r;
	for (int j=0; j < nSpectralSamples; j++ ) {
		r[j] = texelFetch(uPrimarySpectrums, 3 * nSpectralSamples + j, 0).r; 
	}
	return r;
}

float[nSpectralSamples] redSpectrum() {
	float[nSpectralSamples]  r;
	for (int j=0; j < nSpectralSamples; j++ ) {
		r[j] = texelFetch(uPrimarySpectrums, 4 * nSpectralSamples + j, 0).r; 
	}
	return r;
}
float[nSpectralSamples] greenSpectrum() {
	float[nSpectralSamples]  r;
	for (int j=0; j < nSpectralSamples; j++ ) {
		r[j] = texelFetch(uPrimarySpectrums, 5 * nSpectralSamples + j, 0).r; 
	}
	return r;
}

float[nSpectralSamples] blueSpectrum() {
	float[nSpectralSamples]  r;
	for (int j=0; j < nSpectralSamples; j++ ) {
		r[j] = texelFetch(uPrimarySpectrums, 6 * nSpectralSamples + j, 0).r; 
	}
	return r;
}


float[nSpectralSamples] clamp (float[nSpectralSamples] r) {

 for ( int j=0; j < nSpectralSamples; j++ ) {
   if (r[j] < 0 ) r[j] = 0.0;
 }
 return r;
}  // clamp()

float[nSpectralSamples] spectralColor(vec3 rgb) {
// primarySpectrums
// As background -- GLSL looks a lot like C, but compiles a bit different. 
// Things are very unrolled, and conditionals may be executed in parallel and switched at the end,
// that sort of thing. Depends on the hardware...

// You can use **loop indices or constants** [== contant-index-expression] to index into arrays. 
// The assignment in your loop is ok, but the access by tileID isn't.

	float[nSpectralSamples]  r;
	for (int j=0; j < nSpectralSamples; j++ ) {
		r[j] = 0.0;
	}

 
	// Convert reflectance spectrum to RGB
	if (rgb[0] <= rgb[1] && rgb[0] <= rgb[2]  && rgb[1] <= rgb[2] ) {
							// rgb[0] <= rgb[1] << rgb[2], case (1.1)	  

		// Compute reflectance _SampledSpectrum_ with _rgb[0]_ as minimum

		for (int j =0; j < nSpectralSamples; j++) {
			r[j] += rgb[0] * rgbRefl2SpectWhite[j];

		} // for

	 
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] += (rgb[1] - rgb[0]) * rgbRefl2SpectCyan[j];		
			  
			r[j] += (rgb[2] - rgb[1]) * rgbRefl2SpectBlue[j];
	
		} // for
		 
		 
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] *= .94;
		}

		return clamp(r); 

	} // if case (1.1)
  
	// return r;

	if ( rgb[0] <= rgb[1] && rgb[0] <= rgb[2]  && rgb[2] <= rgb[1]   ) { 
						// rgb[0] << rgb[2] << rgb[1]: case (1.2)

		// Compute reflectance _SampledSpectrum_ with _rgb[0]_ as minimum

		for (int j =0; j < nSpectralSamples; j++) {
			r[j] += rgb[0] * rgbRefl2SpectWhite[j];

		} // for

		for (int j =0; j < nSpectralSamples; j++) {
			r[j] += (rgb[2] - rgb[0]) * rgbRefl2SpectCyan[j];
		
			r[j] += (rgb[1] - rgb[2]) * rgbRefl2SpectGreen[j];
			
		} // for

        
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] *= .94;
		}

		return clamp(r); 
		  
	} //  case (1.2)




	if (rgb[1] <= rgb[0] && rgb[1] <= rgb[2]  &&  rgb[0] <= rgb[2] ) { // case (2.1)
						// rgb[1] <= rgb[0] <= rgb[2]: case (2.1)
	   
		// Compute reflectance _SampledSpectrum_ with _rgb[1]_ as minimum
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] += rgb[1] * rgbRefl2SpectWhite[j];

		} // for

   
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] += (rgb[0] - rgb[1]) *rgbRefl2SpectMagenta[j];
			
			r[j] += (rgb[2] - rgb[0]) *rgbRefl2SpectBlue[j];
		
		} // for
		
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] *= .94;
		}

		return clamp(r); 
			
	}// if case (2.1)




	if (  rgb[1] <= rgb[0] && rgb[1] <= rgb[2]  &&  rgb[2] <= rgb[0]  )  { 
					// rgb[1] <= rfb[2] <= rgb[0]: case (2.2)
	      
		// Compute reflectance _SampledSpectrum_ with _rgb[1]_ as minimum
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] += rgb[1] * rgbRefl2SpectWhite[j];

		} // for

		for (int j =0; j < nSpectralSamples; j++) { 
			r[j] += (rgb[2] - rgb[1]) * rgbRefl2SpectMagenta[j];

			r[j] += (rgb[0] - rgb[2]) * rgbRefl2SpectRed[j];
			
		} // for

		//return rgbRefl2SpectWhite;
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] *= .94;
		}

		return clamp(r); 
	} // if case (2.2)



	if ( rgb[2] <= rgb[0] &&  rgb[2] <= rgb[1] && rgb[0] <= rgb[1] )  { 
						// rgb[2] <= rgb[0] <= rgb[1]: case (3.1) 

		// Compute reflectance _SampledSpectrum_ with _rgb[2]_ as minimum
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] += rgb[2] * rgbRefl2SpectWhite[j];
		} // for
       
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] += (rgb[0] - rgb[2]) * rgbRefl2SpectYellow[j];
	       
			r[j] += (rgb[1] - rgb[0]) * rgbRefl2SpectGreen[j]; 

		} // for
		 
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] *= .94;
		}

		return clamp(r); 
	} // if case (3.1) 


  



	if (  rgb[2] <= rgb[0] &&  rgb[2] <= rgb[1] && rgb[1] <= rgb[0]  )   { 
					// rgb[2] <= rgb[1] <= rgb[0]: case (3.2)
			
		// Compute reflectance _SampledSpectrum_ with _rgb[2]_ as minimum
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] += rgb[2] * rgbRefl2SpectWhite[j];
		} // for
		   
		for (int j =0; j < nSpectralSamples; j++) {

			r[j] += (rgb[1] - rgb[2]) * rgbRefl2SpectYellow[j];
			r[j] += (rgb[0] - rgb[1]) * rgbRefl2SpectRed[j];
			   
		} // for

		// return rgbRefl2SpectWhite;
		for (int j =0; j < nSpectralSamples; j++) {
			r[j] *= .94;
		}

		return clamp(r); 

	} // if case (3.2)
 


} // spectralColor

vec3 XYZToRGB(vec3 xyz);

vec3 toXYZ(float[nSamples] spect);
 
vec3 toRGB(float[nSamples] spect) {
    vec3 xyz;
    xyz = toXYZ(spect);
    return XYZToRGB(xyz);
	
}

vec3 toXYZ(float[nSamples] spect) {
    vec3 xyz;
    xyz[0] = xyz[1] = xyz[2] = 0.f;
	
	  
    for (int i = 0; i < nSamples; ++i) {
		
        xyz[0] += X[i] * spect[i];
        xyz[1] += Y[i] * spect[i];
        xyz[2] += Z[i] * spect[i];
    }

    float scale = float( lambdaEnd - lambdaStart) /
                        float(CIE_Y_integral * nSamples);


    xyz[0] *= scale;
    xyz[1] *= scale;
    xyz[2] *= scale;
		
	return xyz;
}//toXYZ

/*

inline void XYZToRGB(const float xyz[3], float rgb[3]) {
    rgb[0] =  3.240479f*xyz[0] - 1.537150f*xyz[1] - 0.498535f*xyz[2];
    rgb[1] = -0.969256f*xyz[0] + 1.875991f*xyz[1] + 0.041556f*xyz[2];
    rgb[2] =  0.055648f*xyz[0] - 0.204043f*xyz[1] + 1.057311f*xyz[2];
}


inline void RGBToXYZ(const float rgb[3], float xyz[3]) {
    xyz[0] = 0.412453f*rgb[0] + 0.357580f*rgb[1] + 0.180423f*rgb[2];
    xyz[1] = 0.212671f*rgb[0] + 0.715160f*rgb[1] + 0.072169f*rgb[2];
    xyz[2] = 0.019334f*rgb[0] + 0.119193f*rgb[1] + 0.950227f*rgb[2];
}*/

vec3 XYZToRGB(vec3 xyz) {
    vec3 rgb;
    rgb[0] =  3.240479f*xyz[0] - 1.537150f*xyz[1] - 0.498535f*xyz[2];
    rgb[1] = -0.969256f*xyz[0] + 1.875991f*xyz[1] + 0.041556f*xyz[2];
    rgb[2] =  0.055648f*xyz[0] - 0.204043f*xyz[1] + 1.057311f*xyz[2];

	return rgb;
}//XYZToRGB





bool rayAABBIntersect ( Ray ray, vec3[2] aabb,
                        out float tmin, out float tmax );


//www.gamedev.net/topic/429443-obb-ray-and-obb-plane-intersection
//www.gamedev.net/topic/346956-adding-vector-from-another-moving-object
//github.com/hpicgs/cgsee/wiki/Ray-Box-Intersection-on-the-GPU  
//http://people.csail.mit.edu/amy/papers/box-jgt.pdf


bool rayOBBIntersect ( vec3 uEyeOriginInBoxFrame, vec3 rayDirInBoxFrame, vec3 minPosInBoxFrame, 
                        vec3 maxPosInBoxFrame, out float tmin, out float tmax) {

	// current position along the ray

	Ray ray = makeRay( uEyeOriginInBoxFrame, rayDirInBoxFrame); 
	vec3[2] AABB = vec3[2]( minPosInBoxFrame, maxPosInBoxFrame); 
 
	return rayAABBIntersect ( ray, AABB, tmin, tmax );

}


 
bool rayAABBIntersect ( Ray ray,  vec3[2] aabb,
                        out float tmin, out float tmax ) {
   // it is assumed that the interval [tmin, tmax] statisy tmin << tmax throughout the code; 
   // To ensure it, we check the sign of each component direction; helps to check whether the intervals overlap
   // Some components of ray vector may be zero or near zero when the ray is parallel or near parallel to axes
   // Improved method for x component
/*  divx = 1 / r.direction.x();
    if (divx >= 0) {
      tmin = (bounds[0].x() - r.origin.x()) * divx; // index 0 [xmin line] is the entering point
      tmax = (bounds[1].x() - r.origin.x()) * divx;
    }
    else {

      tmin = (bounds[1].x() - r.origin.x()) * divx;  // index 1 [xmax line] is the entering point of the ray
      tmax = (bounds[0].x() - r.origin.x()) * divx;
    } 
*/
	tmin = ( aabb[ ray.sign[0] ].x - ray.origin.x ) * ray.inv_direction.x; // ray.sign[0] = 0 when the ray direction is positive

	tmax = ( aabb[1- ray.sign[0] ].x - ray.origin.x ) * ray.inv_direction.x;
	float tymin = ( aabb[ ray.sign[1]  ].y - ray.origin.y ) * ray.inv_direction.y;
	float tymax = ( aabb[ 1 - ray.sign[1]  ].y - ray.origin.y ) * ray.inv_direction.y;
	float tzmin = ( aabb[ ray.sign[2] ].z - ray.origin.z ) * ray.inv_direction.z;
	float tzmax = ( aabb[ 1- ray.sign[2] ].z - ray.origin.z ) * ray.inv_direction.z;
	tmin = max( max( tmin, tymin), tzmin );
	tmax = min( min( tmax, tymax), tzmax );

   
	// check if the ray intersects the AABB

	if ( tmin >  tmax ) {
		// ray does not hit AABB box
	 
		return false;

	}
 
	return true;
		
 }

    


float computeScatCrossSection( float radius, float lambda) {

	// get the index to the texture for (radius, lambda)

	vec2 index = vec2((lambda - lambdaStart) / (lambdaEnd - lambdaStart),
					  (radius - radiusStart) / (radiusEnd - radiusStart));
	
	float scatCrossSection = texture( uScatTex, index).r;
		
    return scatCrossSection;
	
}


float computePhaseFunction( float radius, float lambda, float theta) {
  	
	
	
    float z0 = (radius - radiusStart) / (radiusEnd - radiusStart );
	float y0 = (lambda - lambdaStart) / (lambdaEnd - lambdaStart );
	float x0 = (theta - thetaStart ) / ( thetaEnd - thetaStart );
	


	float  phase = texture(uPhaseTex, vec3(x0,y0,z0) ).r;
	
    return phase;

}

float attenFactor(float height, float h0) {

	if ( height <= h0 ) {
		return 1;
	}
	else {
		return exp(-(height-h0));
	}
}

float avgPhaseFunction(float radius, float lambda, float psi, float height) {

    float particlePhase = computePhaseFunction( radius, lambda, psi ); 

	//float particlePhase = 1.230634e-002;

    float avgPhasePerVolume = uDropDensity * attenFactor(height, h0) * particlePhase;
    return avgPhasePerVolume;

}

float calculate_radiance_of_sun(float lambda) {
    // http://www.oceanopticsbook.info/view/light_and_radiometry/level_2/blackbody_radiation
	// Radiation in thermodynamic equilibrium is isotropic and unpolarized ( eq. (4) )
    float R_Earth = 1.496e+8;
    float R_Sun = 6.95e+5;
    float h = 6.6261e-34;      // Planck's constant
    float c = 2.9979e+8;       // speed of light in vacuo
    float k = 1.3806e-23;      // Boltzmann's constant
    float T = 5782;            // Sun's temperature in Kelvin
    float Isun;

			
	float lambda_m = lambda * 1.0e-9;				// convert lambda to meters from nanometers
	// blackbody radiation spectrum [W/m^2/nm]
    Isun = pow( (R_Sun/R_Earth),2.0 )  * 2*pi*h*pow(c,2.0)  / (   pow(lambda_m, 5.0) * (  exp(  h*c/ ( lambda_m*k*T ) ) -1 )  ) * 1.0e-9  ;
	    //Moon:  1.0e-9 is multiplied to convert the irradiance per meter to the irradiance per nanometer
	    // Original had 1.0e-10.
	float Lsun = Isun / pi; //convert the irradiance to radiance 

	return Lsun;
		
} //calculate_radiance_of_sun()


float[nSpectralSamples] calculate_radiance_of_sun_Lee() {
    
    float R_Earth = 1.496e+8;
    float R_Sun = 6.95e+5;
    float h = 6.6261e-34;      // Planck's constant
    float c = 2.9979e+8;       // speed of light in vacuo
    float k = 1.3806e-23;      // Boltzmann's constant
    float T = 5782;            // Sun's temperature in Kelvin
    float Lsun;

	float[nSpectralSamples] LSpectrumSun;
	
	float lambda = lambdaStart;
	
    for ( int j = 0; j < nSpectralSamples;  j++) {
		
		//Lsun = calculate_radiance_of_sun( lambda );
		//LSpectrumSun[ j ] = Lsun;		
		LSpectrumSun[ j ] = ISpectrumSun[j] / pi; ;		
		//lambda += lambdaStep;
		
	}  // for 
    
    return LSpectrumSun;

} //calculate_radiance_of_sun_Lee()

float[nSpectralSamples] calculate_radiance_of_sun_Plank() {
    
    float R_Earth = 1.496e+8;
    float R_Sun = 6.95e+5;
    float h = 6.6261e-34;      // Planck's constant
    float c = 2.9979e+8;       // speed of light in vacuo
    float k = 1.3806e-23;      // Boltzmann's constant
    float T = 5782;            // Sun's temperature in Kelvin
    float Lsun;

	float[nSpectralSamples] LSpectrumSun;
	
	float lambda = lambdaStart;
	
    for ( int j = 0; j < nSpectralSamples;  j++) {
		
		Lsun = calculate_radiance_of_sun(lambda);
		//Lsun = calculate_radiance_of_sun( lambda );
		//LSpectrumSun[ j ] = Lsun;		
		LSpectrumSun[ j ] = Lsun;
		lambda += lambdaStep;
		
	}  // for 
    
    return LSpectrumSun;

} //calculate_radiance_of_sun_Plank()



float[nSpectralSamples] spectralColor(vec3  rgb );
vec3 toRGB( float[nSpectralSamples] spect );

float[nSpectralSamples] LSpectrumIn( float[nSpectralSamples] LSpectrumSun, 
                                      vec3 currPos, vec3 dirToLight, float sigma_e) {

	// attenuate the sun light Lsun0 from the entering point to AABB to Pos

	float tmin, tmax;
	float [nSpectralSamples] LSpectIn;

	vec3 lightDirInBoxFrame = vec3( uInvModelMatrixAABB * uEyeMatrix * vec4( dirToLight, 0) );
	vec3 lightOriginInBoxFrame = vec3( uInvModelMatrixAABB * uEyeMatrix * vec4( currPos,1) );

	tmin = 0.0; tmax = 0.0;

	bool isHit = rayOBBIntersect ( lightOriginInBoxFrame, lightDirInBoxFrame, 
									uAABBmin, uAABBmax, tmin, tmax);


	if ( isHit ) {
		//spect1 = vec4(tmin, tmax, 0.0, 1.0);
		return LSpectrumSun;

		// the light origin is at currPos [with t_i = 0, within the water view volume] goes to tmax 
		// at which it exits AABB; tmin < 0 is the opposite direction. 

		//   L(t_i, w_fromLight(t_i) ) = LSpectrumIn is given by
		//   tau(tmax_sun, t_i) Lsun( tmax_sun, dirFromLight(tmax_sun) ) 
		//   = exp( -sigma_e ( tmax_sun - t_i) )  Lsun( tmax_sun, dirFromLight(tmax_sun) ) 

		//for (int j=0; j < nSpectralSamples; j++ ) {
		//	LSpectIn[j] = exp( -sigma_e * tmax ) * LSpectrumSun[j];
		//}
	
		//return LSpectIn; //for debugging

		//return LSpectrumSun;

	}
	else {
		//spect1 = vec4(-1,-1, 0.0, 0.0);
		return LSpectrumSun;
	}
  
} // LSpectrumIn()

float sumOfExtinctionCoefficients(float tmax, float tmin, vec3 rayDir, float sigma_e) {
	int N = 20;
	float deltaT = (tmax - tmin) / float(N);

	float sumExtinction = 0.0;

	//for (float t = tmax; t > tmin; t-= deltaT  ) { 
	for ( int i = N-1; i >= 0; i-- ) {
		
		float t = deltaT * i;

		vec3 currPos = vec3(0,0,0) + rayDir * t;

		// apply the exponential density of water drops depending on the height from ground

		vec4 currPosGlobal = uEyeMatrix * vec4(currPos, 1.0); // the global reference frame is ON the ground
		float heightFromGround = currPosGlobal[1]; // y component

		float currPosSigma_e = sigma_e * exp(-attenFactor(heightFromGround, h0));	// assume the original sigma_e is the value on the ground where
																					// height is zero 
		sumExtinction += currPosSigma_e * deltaT;

		//return sumExtinction; 				   
	} // for t: all positions on the current ray	

	return sumExtinction;
}

vec3 singleScatteringAndAttenuation(float radius, vec3 surfaceColor, float zEye, 
                              bool isPointLight, vec3 lightPosOrDir, float[nSpectralSamples] LSpectrumSun,
	                          vec3 rayDir, float sigma_s, float sigma_e )  {
							  

	// (TV_{0}L) (x,w) = int^{x}_{xdV} tau(x',x) sigma_s(x') 
	//                         int_{S^2} f_p (w, x', w') L(x',w') d(w' - dirFromLight(x') ) dw'dx' 
	//                    +  tau(xdV, x) int_{S^2} f_s(w, xdV, w') L(xdV,w') d(w' - dirFromLight(xdV)) dot(n(xdV), w') ] dw'
	//                  = int^{x}_{xdV} tau(x',x) sigma_s(x') [  f_p (w, x', dirFromLight(x') ) L(x', dirFromLight(x') )dx' 
	//                      +  tau(xdV, x)  f_s(w, xdV, dirFromLight(xdV) ) L(xdV,dirFromLight(xdV) dot( n(xdV),  dirFromLight(xdV) )  
	//xdV = tmax, x = tmin:

	// = int^{tmin}_{tmax} tau(t, tmin) sigma_s(t) f_p( dot(-rayDir, dirFromLight(t) ) ) L(t, dirFromLight(t) ) dt 
	// = SUM_{i=0, i=N-1} exp( - sigma_e ( t_i - tmin) ) sigma_s(t_i) f_p( dot(-rayDir, dirFromLight(t_i) ) ) L(t_i, dirFromLight(t_i) ) 
	// = sigma_s SUM_{i=0, i=N-1} exp( - sigma_e ( t_i - tmin) )  f_p( dot(-rayDir, dirFromLight(t_i) ) ) L(t_i, dirFromLight(t_i) )
	//  where   dot(-rayDir, dirFromLight(t_i) ) is the variable for scattering angle.
  
	// To solve it, we need to compute:  L(t_i, w_fromLight(t_i) ), which is given by
	//   L(t_i, w_fromLight(t_i) ) = tau(tmin_sun, t_i) Lsun( tmin_sun, dirFromLight(tmin_sun) ) 
	//   = exp( -sigma_e ( t_i - tmin_sun) )  Lsun( tmin_sun, dirFromLight(tmin_sun) ) 
  

	// check if the current ray through the fullscreen quad intersects the water volume AABB.
	// otherwise, the background color is used as the color at the pixel color of the ray.

	float tmin, tmax;
	float[nSpectralSamples] LSpectrumOut, LSpectrumOut0;
	
	
	vec3 rayDirInBoxFrame = vec3( uInvModelMatrixAABB * uEyeMatrix * vec4(rayDir, 0) );

	bool isHit = rayOBBIntersect ( uEyeOriginInBoxFrame, rayDirInBoxFrame, 
									uAABBmin, uAABBmax, tmin, tmax);
  
  
	//Loop iterations that only some fragments execute are non-uniform control flow.
	// In that case, if the accessed texture uses **mipmapping or anisotropic** filtering of any kind, 
	// then any texture function that is not "Lod" or "Grad" will retrieve undefined results

	if ( !isHit ) { // the AABB is not hit? 	 getBackgroundColor:
  		 
		return surfaceColor;
	  	  	              
	} // if (!isHit) // NOT HIT

	else { //  hit the front at tmin and the back face of AABB box at tmax;
  

		// if the background scene lies within the AABB box, the attenuation of the background color
		// begins at the intersection point between the ray and the background scene. The 
		// attenuation of the volume color begins also at the intersection point. 
		// Otherwise, the attenuation of background and volume color begins at the back face of the box.

	  	  	  
		// compute the intensity of the light with wavelenth lambda scattered and attenuated 
		// in the direction of -rayDir
     
	  
		float tmaxAtten = tmax; // the t value at which attenuation begins

		float distFromEye = -zEye;

		if ( distFromEye < tmin ) { // the surface is before the AABB box => return the surfaceColor
			
			return surfaceColor;
		
		}

		if ( distFromEye >= tmin && distFromEye <= tmax ) {
			tmaxAtten = distFromEye;
		}
		else { // the surface is behind the AABB box
			tmaxAtten = tmax; 
		}
	  	    
		int N =30;
		float deltaT = (tmaxAtten - tmin) / float(N);
   
   
		// get the light which has been scattered in the eye direction by each region in the volume 
		// starting from tmaxAtten.

		float [nSpectralSamples] LSpectIn;

		for (int j = 0; j < nSpectralSamples; j++ ) {
			LSpectrumOut[j] = 0.0;
		}

		//float [nSpectralSamples] spectralSurfaceColor = spectralColor (surfaceColor);

		// LIGHT DIRECTIONS
		if ( isPointLight ) { // point light
      
			vec3 lightPos = lightPosOrDir;

			for (float t = tmaxAtten; t >= tmin; t -= deltaT) {  // for each position on the ray
      
				vec3 currPos = vec3(0,0,0) + rayDir * t;

				vec3 dirFromPointLight = currPos - lightPos;

				dirFromPointLight = normalize( dirFromPointLight );

				LSpectIn = LSpectrumIn( LSpectrumSun, currPos, -dirFromPointLight, sigma_e);
		
				// The light with given lambda is scattered in the direction of -rayDir with varying amount
				// depending on the scattering angle between the light direction and -rayDir, which varies
				// with each t. Consider only the scattering angles  within the range [ thetaStart, thetaEnd]. 
				// The other scattering angles constributs less to the scattered intensity of lambda
				// They are ignored for computational reasons. Some scattering spectrum is computed for
				// every ray direction (-rayDir)?? 
	  
	            // apply the exponential density of water drops depending on the height from ground

				vec4 currPosGlobal  = uEyeMatrix * vec4( currPos, 1.0); // the global reference frame is ON the ground
				float heightFromGround = currPosGlobal[1]; // y component

				sigma_e = sigma_e * exp( -heightFromGround); // assume the original sigma_e is the value on the ground where
				                                             // height is zero 
			    sigma_s = sigma_s * exp( -heightFromGround);


				// compute the scattered light scattered  along the ray
				float psi = acos ( dot ( -rayDir, dirFromPointLight) ) * 180.0 / pi;
		
				if (  psi >=  thetaStart  &&  psi <= thetaEnd  ) {	// only the light [at the current position on the ray]
																	// whose scattering angle is  within these angles
																	// contribute to the rainbow color along the current ray.
																	// Otherwise, no color is contributed.  
                 

					float lambda = lambdaStart;
					int   lambdaIndex = 0;	
					for ( int j = 0; j < nSpectralSamples;  j++) {

		                float avgPhaseExp =  avgPhaseFunction(radius, lambda, psi, heightFromGround) * exp( - heightFromGround);

						LSpectrumOut[j] += exp( - sigma_e * ( t - tmin) )
										* sigma_s * avgPhaseExp * LSpectIn[j];
						lambda += lambdaStep;
			   
			   
					} // for
		   	   
				}   //if 

				else continue; // go to the next position

			} // for all positions on the current ray

			// for debugging	
			//for ( int j =0; j < nSpectralSamples; j++) {	
			//		LSpectrumOut[j] += spectralSurfaceColor [j]* exp( - sigma_e * ( tmaxAtten - tmin) );
			//}

			return oldToRGB( LSpectrumOut  ) + surfaceColor  * exp( - sigma_e * ( tmaxAtten - tmin) ); 
			//return toRGB( LSpectrumOut );

		} // if (point light)
		
		else { // directional light

			vec3 dirFromSun = lightPosOrDir;
			float psi = acos ( dot ( -rayDir, dirFromSun ) ) * 180.0 / pi;

			// for debugging
			//spect1 = vec4(dirFromSun, 1.0);

			//spect1 = vec4(dot ( -rayDir, dirFromSun ), acos(dot ( -rayDir, dirFromSun ) ), psi, 1.0 );
			//return LSpectrumOut;

			// consider only the cases where the light is scattered in the direction of the eye
			// with scattering angle within [ thetaStart, thetaEnd]. These scattering angles contricute
			// to the rainbow and its neighborhood. The cases with other scattering
			// angles contribute to the further region outside of the rainbow and its neighbor.
			// These regions are not considered and the background color are assumed to dominate. 


			// SUM_{i=0, i=N-1} exp( - sigma_e ( t_i - tmin) ) sigma_s  f_p( dot(-rayDir, dirFromLight(t_i) ) ) L(t_i, dirFromLight(t_i) )
			// L(t_i, dirFromLight(t_i) ) = LSpectrumIn is given by
			// tau(tmin_sun, t_i) Lsun( tmin_sun, dirFromLight(tmin_sun) ) 
			// = exp( -sigma_e ( t_i - tmin_sun) )  Lsun( tmin_sun, dirFromLight(tmin_sun) ) 
			// where tmin_sun is the point at which the sun light enters the volume to reach t_i.
	  

			if ( !( psi >= thetaStart  &&  psi <= thetaEnd ) ) {
				
				return surfaceColor; // just return the background color
		 
			}  

			else {	// psi is a rainbow scattering angle
				// LSpectrumOut = LSpectrumOut + SUM_{t = tmin, tmaxAtten} exp( - sigma_e * ( t - tmin) )  * sigma_s 
				//                                          * avgPhaseFunction( dot ( -rayDir, dirFromSun(t) ) ) 
				//                                          * LSpectIn(t, dirFromLight(t) )

				//return LSpectrumOut;
				//float [nSpectralSamples] LSpectIn;

				float [nSpectralSamples] incrementalSpect;
				
				//for (float t = tmaxAtten; t > tmin; t-= deltaT  ) { 
				for (int i = N-1; i >= 0; i--) { 
					
					float t = deltaT * i;
					//LSpectrumOut = LSpectrumOut +  exp( - sigma_e * ( t - tmin) )  * sigma_s 
					//									 * avgPhaseFunction( dot ( -rayDir, dirFromSun(t) ) * LSpectrumIn
				  
					vec3 currPos = vec3(0,0,0) + rayDir * t;
		 		     
					LSpectIn = LSpectrumIn(LSpectrumSun, currPos, -dirFromSun, sigma_e);

					// apply the exponential density of water drops depending on the height from ground

					vec4 currPosGlobal  = uEyeMatrix * vec4(currPos, 1.0); // the global reference frame is ON the ground
					float heightFromGround = currPosGlobal[1]; // y component

					float currPosSigma_e = sigma_e * exp(-attenFactor(heightFromGround, h0)); // assume the original sigma_e is the value on the ground where
																// height is zero 
					float currPosSigma_s = sigma_s * exp(-attenFactor(heightFromGround, h0));
		
		         
					float lambda = lambdaStart;
		    		
					float sumOfExtinctionToFrontSide = sumOfExtinctionCoefficients(t, tmin, rayDir, sigma_e);

					for ( int j =0; j < nSpectralSamples; j++) {	
							  
				        float avgPhaseExp = avgPhaseFunction(radius, lambda, psi, heightFromGround);
						 //* exp(-attenFactor(heightFromGround, h0));

						//incrementalSpect[j] = exp( - sigma_e * ( t - tmin) ) * sigma_s  
						//							* avgPhaseExp * LSpectIn[j];

						incrementalSpect[j] = exp(-sumOfExtinctionToFrontSide) * currPosSigma_s 
													* avgPhaseExp * LSpectIn[j] * deltaT;

						//incrementalSpect[j] = exp( - sigma_e * ( t - tmin) ) * sigma_s  
						//                          * avgPhaseFunction(radius, lambda, psi) * LSpectrumSun[j];
						//incrementalSpect[j] = exp( - sigma_e * ( t - tmin) ) * sigma_s  
						//                          * avgPhaseFunction(radius, lambda, psi); 
						//LSpectrumOut[j] += exp( - sigma_e * ( t - tmin) ) * sigma_s  
						//							* avgPhaseFunction(radius, lambda, psi) * LSpectIn[j];
			   
						lambda += lambdaStep;			  
					} // for each lambda	

			
					for ( int j =0; j < nSpectralSamples; j++) {			  
						LSpectrumOut[j] += incrementalSpect[j];						  			  
					} // for each lambda	

				} // for t: all positions on the current ray	


				//float weight1 = (psi - 130) / (137.5 - 130);
				//float weight2 = (142 - psi) / (142 - 140.5);

				//vec3 rgb = oldToRGB( LSpectrumOut );

				//if ( psi < 137.5 ) {
				//	return rgb*(weight1) + surfaceColor * exp(-sumOfExtinctionCoefficients( tmaxAtten, tmin, rayDir, sigma_e ) );
				//}
				//else if ( psi > 140.5 ) {
				//	return rgb*(weight2) + surfaceColor * exp(-sumOfExtinctionCoefficients( tmaxAtten, tmin, rayDir, sigma_e ) );
				//}
				//else {
				//	return rgb + surfaceColor * surfaceColor * exp(-sumOfExtinctionCoefficients( tmaxAtten, tmin, rayDir, sigma_e ) );
				//}

				return oldToRGB( LSpectrumOut ) 
						+ surfaceColor * exp(-sumOfExtinctionCoefficients( tmaxAtten, tmin, rayDir, sigma_e ));


				//return oldToRGB( LSpectrumOut ) * 0.7 + surfaceColor * 0.3; // the light accumulated along the current ray
				//return toRGB( LSpectrumOut );
				
			} // else ( psi is rainbow angle)

	    
		} // else (DIRECTIONal light) 	  


	} // else (the ray HITs the volume )

} // singleScatteringAndAttenuation


vec3 singleScatteringAndAttenuation_stan(float radius, vec3 surfaceColor, float zEye,
                              vec3  dirFromSun,  float [nSpectralSamples] LSpectrumSun, 
							  vec3 rayDir, float sigma_s, float sigma_e ) {

	// check if the current ray through the fullscreen quad intersects the water volume AABB.
	// otherwise, the background color is used as the color at the pixel color of the ray.

	float tmin, tmax;
	float[nSpectralSamples] LSpectrumOut, LSpectrumOut0;
	
	
	vec3 rayDirInBoxFrame = vec3( uInvModelMatrixAABB * uEyeMatrix * vec4(rayDir, 0) );

	bool isHit =  rayOBBIntersect ( uEyeOriginInBoxFrame, rayDirInBoxFrame, 
									uAABBmin, uAABBmax, tmin, tmax);
  
	if ( !isHit ) { // the AABB  is not hit? 	 getBackgroundColor:
  		//ivec2 pixelPos = ivec2( gl_FragCoord.xy);
		//vec3 surfaceColor = vec3( texelFetch(uColorTex, pixelPos, 0) ); // color
		//return vec3(0,0,0);      	     	 
		return  surfaceColor;
	  	              
	} // if (!isHit)

	else { //  hit the front at tmin and the back face of AABB box at tmax;
  

		// if the background scene lies within the AABB box, the attenuation of the background color
		// begins at the intersection point between the ray and the background scene. The 
		// attenuation of the volume color begins also at the intersection point. 
		// Otherwise, the attenuation of background and volume color begins at the back face of the box.
	 
	     
		// compute and attenuate the intensity of the light scattered by the water volume
		// in the direction of -rayDir

		float scatCrossSection = computeScatCrossSection( radius, lambdaStart );

		float psi = acos ( dot ( -rayDir, dirFromSun ) ) * 180.0 / pi;

		if ( !( psi >=  thetaStart  &&  psi <= thetaEnd ) ) {
			//ivec2 pixelPos = ivec2( gl_FragCoord.xy);
			//vec3  surfaceColor = vec3( texelFetch(uColorTex, pixelPos, 0) ); // color
   			//return vec3(0,0,0);

			return  surfaceColor; // just return the background color
		}  

		// within the rainbow angles
   
		float phi_obs = acos( dot( uFrontNormalToAABB, -rayDir) );       //  rayDir in eyeSpace
		float phi_sun = acos( dot( uFrontNormalToAABB,  uSunRayDir )  );

		// compute the optical depth normal to the AABB

		float rayPathNormal =  (tmax - tmin) * cos (phi_obs);

		float tau_N = uDropDensity * scatCrossSection  * rayPathNormal;	
	  	
		// compute the surface color at vPosition and attenuated the color through
		// volume space between tmin and tmax.

		//  tau_N = 1.0;



		float lambda = lambdaStart;

		for (int j = 0; j < nSpectralSamples; lambda += lambdaStep, j++) {
		
			//float Lsun = calculate_radiance_of_sun(lambda); // Lsun: radiance per nanometer
			float Lsun = LSpectrumSun[j];

			// ignore the attenuation of the sun light while it passes through the drop volume
				
			float particlePhase= computePhaseFunction( radius,  lambda, psi); 
		
			//float particlePhase = 1.230634e-002;

			float avgPhasePerVolume  =  uDropDensity * particlePhase;
		        
			LSpectrumOut[j] = ( ( Lsun* avgPhasePerVolume  / cos(phi_obs) ) / ( 1/cos(phi_sun) + 1/cos(phi_obs) ) )
							* ( 1 - exp(-tau_N * ( 1/cos(phi_sun) + 1/cos(phi_obs) ) ) );
         
        
		} // for j
		//return vec3(1,0,0);

		return oldToRGB( LSpectrumOut );
		//return toRGB( LSpectrumOut );

	} // else (the ray hits the volume )

} // singleScatteringAndAttenuation_stan


vec3 calculate_rainbowColor (float radius, vec3 surfaceColor, float zEye,  bool isPointLight,  vec3 lightPosOrDir, vec3 rayDir, 
                             float sigma_s, float sigma_e) {
	// compute the scattered light, which is attenuated along the direction 
	// -rayDir 
	
	vec3 rainbowRGB;
		
	//float[nSpectralSamples] LSpectrumSun = calculate_radiance_of_sun_Lee(); 
	float[nSpectralSamples] LSpectrumSun = calculate_radiance_of_sun_Plank(); 



	// LSpectrumSun: radiances for wavelengths [nanometer]

	rainbowRGB = singleScatteringAndAttenuation(radius,  surfaceColor, zEye,
	                          isPointLight, lightPosOrDir, LSpectrumSun,
	                          rayDir, sigma_s, sigma_e );
    
	return rainbowRGB;
	  
}//calculate_rainbowColor

vec3 calculate_rainbowColor_stan (float radius,  vec3 surfaceColor, float zEye, vec3 dirFromSun, vec3 rayDir, 
                             float sigma_s, float sigma_e) {
	// compute the scattered light, which is attenuated along the direction 
	// -rayDir 
	
	vec3 rainbowRGB;
		
	float [nSpectralSamples] LSpectrumSun = calculate_radiance_of_sun_Lee(); 



	// LSpectrumSun: radiances for wavelengths [nanometer]

	rainbowRGB = singleScatteringAndAttenuation_stan(radius,  surfaceColor, zEye, 
	                          dirFromSun, LSpectrumSun, rayDir, sigma_s, sigma_e );
    return rainbowRGB;
	 
}//calculate_rainbowColor_stan



void main() {
	// access textures in uniform contrl flow

	// access textures in uniform contrl flow

	X = XSpectrum();
	Y = YSpectrum();
	Z = ZSpectrum();


	rgbRefl2SpectWhite = whiteSpectrum();
	rgbRefl2SpectCyan = cyanSpectrum();
	rgbRefl2SpectMagenta = magentaSpectrum();
	rgbRefl2SpectYellow = yellowSpectrum();
	rgbRefl2SpectRed = redSpectrum();
	rgbRefl2SpectGreen = greenSpectrum();
	rgbRefl2SpectBlue = blueSpectrum();


	// get the pixel coordinates of the current fragment
  
	ivec2 pixelPos = ivec2( gl_FragCoord.xy );

	vec3 surfaceColor = vec3( texelFetch(uColorTex, pixelPos, 0) ); // color
  
	float zWin  = texelFetch(uDepthTex, pixelPos, 0).r;  // zWin in non linear range [0,1]
	// conversion range ([0,1] into NDC range  [-1,1]
      
	float zNDC = zWin * 2.0 - 1.0; // ( zWin = 1/2 zNDC + 1/2 by Viewport transformation: zWin in [0,1] )

	// This uProjectionMatrix should be the same as the projectionMatrix when the background scene was rendered
	float alpha = uProjectionMatrix[2][2];
	float beta = uProjectionMatrix[3][2]; // glsl matrix is accessed with column first order
	float zEye = - beta / ( zNDC + alpha ); 

	vec3 volumeColor;
  
   
	float scatCrossSection = computeScatCrossSection( uRadius,  lambdaStart);	// scatCrossSection is independent of
																				// lambda, so use lambdaStart
  																			    
	vec3 rayDir = normalize( vPosition ); // the direction from pixel to fragment position vPosition
  
     
	float sigma_s = uDropDensity * scatCrossSection;	// scattering coefficient [1/m]
	float sigma_e = sigma_s; // assume that the water drop not absorb light, the extinction coeff = scat coeff. 

	// compute the surface color at vPosition and attenuate the color through
	// volume space between tmin and tmax.

	bool isPointLight = false;
	vec3 lightPosOrDir;

	// compute the spectrum scattered in the direction of -rayDir 
	// and get the RGB color from the spectrum

	if ( isPointLight ) {
        

		vec3 lightPos = uLight1Pos;  
		volumeColor = calculate_rainbowColor (uRadius, surfaceColor, zEye, isPointLight, lightPos, rayDir, sigma_s, sigma_e); 
		
		//lightPos = uLight2Pos;
		//volumeColor += calculate_rainbowColor (uRadius, surfaceColor, zEye, isPointLight, lightPos, rayDir,sigma_s, sigma_e); 
	
    }
	else {

		// uSunRayDir is  already conveted into eye space in the main program.
	    vec3 dirFromSun  = -uSunRayDir;
	    
		//volumeColor = calculate_rainbowColor_stan (uRadius,  surfaceColor, zEye, dirFromSun, rayDir, sigma_s, sigma_e);  
		
		volumeColor = calculate_rainbowColor (uRadius, surfaceColor, zEye, isPointLight, dirFromSun, rayDir, sigma_s, sigma_e);  

	}
	
	fragColor = vec4(volumeColor, 1);
	   
	// no volume color for rainbow:  The red emerges at a smaller angle psi than the violet
	// scattering is concentrated at theta_0 = 137.97 but also occurs for angles greater than that. 
	// This is what causes the white appearance `insidea the primary bow (and `outsidea the secondary bow).
	// Note the gap between 129 and 138. This is a region of negligible scattering (from higher orders) and appears 
	// as a dark space between the primary and secondary rainbows
   
	   
} // main




