
#version 140
// OpenGL defines a pixel's center, whose location is 0.5 pixels up and to the right of the lower left 
//corner of the pixel (which defines the pixel). In other words, pixel 0,0 goes from 0,0 to 1,1 
//and its sample-center is 0.5, 0.5. To take a simple polygon example: if we draw a rectangle 
//from 0,0 to 8,8 then the first 8x8 pixels will be included. Why? Well, the sampling centers of
// that first row are: 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5. 
// Those 8 pixels are all clearly inside the range (0-8).


// The GL will sample the texture for "each pixel filled in by our polygon". 
// So what texture coordinates are used per pixel?
// The texture coordinates (and all fragment information) is interpolated per fragment (pixel) 
// at the pixel center. we are sampling at 0.5, 1.5, 2.5, etc. So our texture coordinates are
//  going to be: 0.0625, 0.1875, 0.3125, ... 0.9375.
// if you use integer input coordinates, the corners of a rectangle do not receive 
// the unmodified input coordinates from GL vertex! There will be a very slight modification
//  to the values due to the need to sample the interior.
// Multiplying out that horrible sequence for our 4x4 texture on an 8x8 quad we get: 
// 0.25, 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75.

// Nearest neighbor texture filtering is a simple floor function on the UV coordinates. 
// If there's a rule of thumb, it is: given integral alignment of a texture's texels over screen pixels,
//  nearest neighbor will give you a clean copy.


// 
// Shaders have access to all the texture units, they don't care about the active texture.



#define MAX_ITERATIONS 100

/* http://http.developer.nvidia.com/GPUGems2/gpugems2_chapter27.html
So the first step in writing a filter shader is to understand how to compute the integer sample coordinate 
corresponding to the current texture coordinate and how to convert sample coordinates back into conventional
 texture coordinates. The following equations comply with Direct3D's texture-mapping specifications:

sampleCoord = floor(textureCoord * texSize);

textureCoord = (sampleCoord + 0.5) / texSize;

One additional tip is that it is not necessary to snap texture coordinates to the center of the nearest texel 
as long as texture lookups do not perform any type of filtering. This means setting filtering to GL_NEAREST in OpenGL.
 This way, coordinate snapping is performed automatically during texture lookups
*/

//our screen resolution, set whenever the display is resized

//uniform float uResolution[2];


uniform int uWindowWidth; 
uniform int uWindowHeight; 

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
uniform mat4 uInvModelMatrixAABB;
uniform mat4 uModelViewMatrix;
uniform mat4 uProjectionMatrix;

uniform sampler2D uColorTex;
uniform sampler2D uDepthTex;

uniform sampler3D uDistTex;
uniform sampler3D uPhaseTex;
uniform sampler2D uScatTex;

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


in vec2 vTexCoord;
in vec3 vNormal;
in vec3 vPosition;
in   mat3 vNTMat;

out vec4 fragColor;
out vec4 phaseSpectrum;
out vec4 rainbowColorSpectrum;


// parameters for rainbow computation

float pi = acos( -1.0);  // -1.0 const float in default; -1.0f = float

const int nSpectralSampleSteps = 120;
const int nThetaSteps = 100; // 100 steps between theta 130 deg and theta 142;

const int nRadiusSteps  = 2;

const int nSpectralSamples = nSpectralSampleSteps +1;
const int nRadii = nRadiusSteps +1;
const int nThetas = nThetaSteps +1; // 100 steps between theta 130 deg and theta 142;

//uniform float uPhaseFunction[nRadii * nSpectralSamples * nThetas];

float irainbow[ nSpectralSamples ];
float particlePhase[nSpectralSamples ];


const float  lambdaStart = 400;   // 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
const float lambdaEnd = 700;
     
const  float radiusStart =  1.0e-3;  // 
const  float radiusEnd = 2.0e-3;           // 2 mm = 2 * e-3 m

const	float thetaStart = 130.0;
const	float thetaEnd = 142.0;

const float lambdaStep = ( lambdaEnd - lambdaStart) / float( nSpectralSamples -1);
const float thetaStep = ( thetaEnd - thetaStart) / float (nThetas -1);
const float radiusStep = (radiusEnd - radiusStart) / float (nRadii -1);

const int nLambdasForColorMatch = 81;
const int nRGB = 3;


/* A colour system is defined by the CIE x and y coordinates of
   its three primary illuminants and the x and y coordinates of
   the white point. */

struct colourSystem {
 //   uint name[15];     	    	    /* Colour system name */
    float xRed, yRed,	    	    /* Red x, y */
           xGreen, yGreen,  	    /* Green x, y */
           xBlue, yBlue,    	    /* Blue x, y */
           xWhite, yWhite,  	    /* White point x, y */
	   gamma;   	    	    /* Gamma correction for system */
};


/* White point chromaticities. */


#define IlluminantCxWhite     0.3101	    	/* For NTSC television */
#define IlluminantCyWhite  0.3162
#define IlluminantD65xWhite   0.3127    	/* For EBU and SMPTE */
#define IlluminantD65yWhite 0.3291	
#define IlluminantExWhite 	0.33333333 /* CIE equal-energy illuminant */
#define IlluminantEyWhite  0.33333333 
#define GAMMA_REC709	0	                 	/* Rec. 709 */

               /* Name                  xRed    yRed    xGreen  yGreen  xBlue  yBlue    White point        Gamma   */
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

colourSystem g_cs = colourSystem(  0.67,  0.33,  0.21,  0.71,  0.14, 0.08, IlluminantCxWhite, IlluminantCyWhite, GAMMA_REC709);


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
		



// uniform block for primary spectrums

uniform PrimarySpectrums {

float [nSpectralSamples] rgbRefl2SpectWhite;
float [nSpectralSamples] rgbRefl2SpectCyan;
float [nSpectralSamples] rgbRefl2SpectMagenta;
float [nSpectralSamples] rgbRefl2SpectYellow;
float [nSpectralSamples] rgbRefl2SpectRed;
float [nSpectralSamples] rgbRefl2SpectGreen;
float [nSpectralSamples] rgbRefl2SpectBlue;
} primarySpectrums;




struct Ray {
  vec3 origin;
  vec3 direction;
  vec3 inv_direction;
  int sign[3];
};




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

vec3 posInAABB (vec3 posInEyeSpace) {

// aabb = AABB in eye space

//  posInEyeSpace = uAABBmin + posInAABB
 
	vec3 widths = uAABBmax - uAABBmin; 

	vec3 posInAABB  = posInEyeSpace - uAABBmin;
	
	vec3 texCoord = posInAABB  / widths;

	return texCoord;
	
}


bool rayAABBIntersect ( Ray ray,  vec3[2] aabb,
                        out float tmin, out float tmax );


//www.gamedev.net/topic/429443-obb-ray-and-obb-plane-intersection
//www.gamedev.net/topic/346956-adding-vector-from-another-moving-object
//github.com/hpicgs/cgsee/wiki/Ray-Box-Intersection-on-the-GPU  
//http://people.csail.mit.edu/amy/papers/box-jgt.pdf


bool rayOBBIntersect ( vec3 uEyeOriginInBoxFrame, vec3 rayDirInBoxFrame, vec3 minPosInBoxFrame, 
                        vec3 maxPosInBoxFrame, out float tmin, out float tmax) {

 // current position along the ray


  
 Ray ray =  makeRay( uEyeOriginInBoxFrame,  rayDirInBoxFrame); 
 vec3[2]  AABB = vec3[2]( minPosInBoxFrame, maxPosInBoxFrame); 
 
 return rayAABBIntersect (ray, AABB,  tmin,  tmax );

}


 
bool rayAABBIntersect ( Ray ray,  vec3[2] aabb,
                        out float tmin, out float tmax )
   
 {
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

    
 


//Cvec3  xyz_to_rgb(struct colourSystem *cs,   const vec3 xyz) {

vec3  xyz_to_rgb( colourSystem cs,   const vec3 xyz) {

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

	
}



bool constrain_rgb(inout vec3 rgb) {
    float w;

    // Amount of white needed is w = - min(0, *r, *g, *b) 
	    
	if (0.0 < rgb[0]) w = 0.0; else w = rgb[0];
	if  (w < rgb[1])  w = 2; else w =  rgb[1];
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

void norm_rgb( inout vec3 rgb)
{
#define Max(a, b)   (((a) > (b)) ? (a) : (b))
	float greatest = Max(rgb[0], Max(rgb[1], rgb[2]));
    
    if (greatest > 0) {
		rgb[0] /= greatest;
		rgb[1] /= greatest;
		rgb[2] /= greatest;
    }
#undef Max
}


// To be modified

float computeScatCrossSection( float radius, float lambda) {

	// get the index to the texture for (radius, lambda)

	vec2 index = vec2(  
	               (lambda - lambdaStart) / ( lambdaEnd - lambdaStart),
				   (radius  - radiusStart) / (radiusEnd -  radiusStart)
				    );
	
	float scatCrossSection = texture( uScatTex, index).r;

    return scatCrossSection;
	
}


float computePhaseFunction( float radius, float lambda, float theta) {
  	
	
	
    float z0 =   (radius  - radiusStart) / (radiusEnd - radiusStart );
	float y0 =  (lambda - lambdaStart) / (lambdaEnd - lambdaStart  );
	float x0 =  (theta - thetaStart ) / ( thetaEnd - thetaStart );
	


	float  phase = texture(uPhaseTex, vec3(x0,y0,z0) ).r;
	
    return phase;

}

float avgPhaseFunction(float lambda, float psi) {

    float particlePhase  = computePhaseFunction( uRadius,  lambda, psi); 

	//float particlePhase = 1.230634e-002;

    float avgPhasePerVolume  =  uDropDensity * particlePhase;
    return avgPhasePerVolume;

}

float[nSpectralSamples] calculate_irradiance_of_sun() {
    
    float R_Earth = 1.496e+8;
    float R_Sun = 6.95e+5;
    float h = 6.6261e-34;      // Planck's constant
    float c = 2.9979e+8;       // speed of light in vacuo
    float k = 1.3806e-23;      // Boltzmann's constant
    float T = 5782;            // Sun's temperature in Kelvin
    float Isun;

	float[nSpectralSamples] ISpectrumSun;

    for (float lambda = lambdaStart, j = 0; j < nSpectralSamples; lambda += lambdaStep, j++) {
		
		lambda_m = lambda * 1.0e-9;				// convert lambda to meters from nanometers

        Isun = pow( (R_Sun/R_Earth),2.0 )  * 2*pi*h*pow(c,2.0)  / ( pow(lambda, 5.0) * ( exp( h*c/ (lambda*k*T) ) -1 ) ) * 1.0e-9  ;
	      //Moon:  1.0e-9 is multiplied to convert the irradiance per meter to the irradiance per nanometer
	      // Original had 1.0e-10.
		ISpectrumSun[j] = Isun;
	}   
    
    return ISpectrumSun;

} //calculate_irradiance_of_sun()


float[nSpectralSamples] spectralColor(vec3  rgb);

float[nSpectralSamples]  LSpectrumIn(float[nSpectralSamples] LSpectrumSun,  vec3 currPos, vec3 dirToLight, float sigma_e) {

  // attenuate the sun light Lsun0 from the entering point to AABB to Pos

  float tmin, tmax;


  vec3 lightDirInBoxFrame = vec3( uInvModelMatrixAABB * uEyeMatrix * vec4( dirToLight, 0) );
  vec3 lightOriginInBoxFrame = vec3( uInvModelMatrixAABB * uEyeMatrix * vec4( currPos,1) );

  bool isHit =   rayOBBIntersect ( lightOriginInBoxFrame, lightDirInBoxFrame, 
                                  uAABBmin, uAABBmax, tmin, tmax);

   //bool isHit = false;

  if ( isHit ) {
      // 
	 // the light origin is at currPos [with t = 0, within the water view volume] goes to tmax 
	 // at which it exits AABB; tmin < 0 is the opposite direction. 
	  
	  return exp( - sigma_e * tmax ) * LSpectrumSun;
  }
  else {
    return LSpectrumSun;
	}
  
} // LSpectrumIn()

//singleScatteringAndAttenuation(isPointLight, lambda, Lsun, sigma_s, sigma_e, lightPosOrDir, rayDir, tmin, tmax);

float[nSpectralSamples] singleScatteringAndAttenuation(bool isPointLight, 
                                float[nSpectralSamples] LSpectrumSun, 
                               float  sigma_s, float sigma_e, vec3 lightPosOrDir, vec3 rayDir, float tmin, float tmax) {

  // (TV_{0}L) (x,w) = int^{x}_{x'dV} tau(x",x) sigma_s(x") [ int_{S^2} f_p (w, x", w') L(x",w') d(w' - w_fromLight(x") ) dw' 
  //   +  f_s(w, x', w') L_surface(x'dV) ] dx"
  //                  = int^{x}_{x'dV} tau(x",x) sigma_s(x") [  f_p (w, x", w_fromLight(x") ) L(x", w_fromLight(x") ) 
  //     +    f_s(w, x', w') L_surface(x'dV)  ] dx'
  //x'dV = tmax, x = tmin:

  // = int^{tmin}_{tmax} tau(t, tmin) sigma_s(t) f_p( dot(-rayDir, w_fromLight(t) ) ) L(t, w_fromLight(t) ) dt 
  // = SUM_{i=0, i=N-1} exp( - sigma_e ( t_i - tmin) ) sigma_s(t_i) f_p( dot(-rayDir, w_fromLight(t_i) ) ) L(t_i, w_fromLight(t_i) ) 
  // = sigma_s SUM_{i=0, i=N-1} exp( - sigma_e ( t_i - tmin) )  f_p( dot(-rayDir, w_fromLight(t_i) ) ) L(t_i, w_fromLight(t_i) )
  //  where   dot(-rayDir, w_fromLight(t_i) ) is the variable for scattering angle.
  
  // To solve it, we need to compute:  L(t_i, w_fromLight(t_i) ), which is given by
  //   L(t_i, w_fromLight(t_i) ) = tau(tmin_sun, t_i) Lsun( tmin_sun, w_fromLight(tmin_sun) ) 
  //   = exp( -sigma_e ( t_i - tmin_sun) )  Lsun( tmin_sun, w_fromLight(tmin_sun) ) 
  // = exp( -sigma_e ( t_i - tmin_sun) ) Lsun0 ( simplification of L)

  // check if the current ray through the fullscreen quad intersects the water volume AABB.
  // otherwise, the background color is used as the color at the pixel color of the ray.

  float[nSpectralSamples] LSpectrumOut;
  
  for (int i=0; i < nSpectralSamples; i++ ) {
    LSpectrumOut[i] =0.0;
  }


  
  vec3 rayDirInBoxFrame = vec3( uInvModelMatrixAABB * uEyeMatrix * vec4(rayDir, 0) );

  bool isHit =  rayOBBIntersect ( uEyeOriginInBoxFrame, rayDirInBoxFrame, 
                                  uAABBmin, uAABBmax, tmin, tmax);
   
  if ( !isHit ) { // the AABB  is not hit? 	 getBackgroundColor:

      // get the pixel coordinates of the current fragment
	  ivec2 pixelPos = ivec2( gl_FragCoord.xy);
	  vec3  surfaceColor = vec3( texelFetch(uColorTex, pixelPos, 0) ); // color
	  LSpectrumOut = spectralColor(surfaceColor);
	  return LSpectrumOut; // this surface color is NOT attenuated through the water volume
	                       // and reaches the eye without attenuation
	              
  } // if (!isHit)

  else { //  hit the front at tmin and the back face of AABB box at tmax; 	 
	
	 // if the background scene lies within the AABB box, the attenuation of the background color
	 //  begins at the intersection point between the ray and the background scene. The 
	 //  attenuation of the volume color begins also at the intersection point. 
	 //  Otherwise, the attenuation of background and volume color begins at the back face of the box.

	 // get the pixel coordinates of the current fragment

	  ivec2 pixelPos = ivec2( gl_FragCoord.xy);

	  // get the surface color and depth of the current fragment
	 
	  vec3  surfaceColor = vec3( texelFetch(uColorTex, pixelPos, 0) ); // color

	  float[nSpectralSamples] LSpectrumOut0 = spectralColor(surfaceColor);


	  
  // compute the intensity of the light with wavelenth lambda scattered and attenuated 
  // in the direction of -rayDir

	  float zpixel  = texelFetch(uDepthTex, pixelPos, 0);  // z in non linear range [0,1]
     
      // conversion into NDC [-1,1]
      float zndc = zpixel * 2.0 - 1.0;

	  float zeye = uProjectionMatrix[3][2]/ ( zndc - uProjectionMatrix[2][2] ); 
	  float tmaxAttenuation; // the t value at which attenuation begins

	  if ( zeye >= tmin && zeye <= tmax ) {
	     tmaxAttenuation = zeye;
	  }
	  else {
	   tmaxAttenuation = tmax; 
	  }



  
     int N = 30;
     float deltaT  =  (tmax - tmin) / N;
   

  // get the surface reflection color which has been attenuated through the volume
  // assume that the surface color attenuated only through water volume and not through atmosphere

     LSpectrumOut = exp( - sigma_e * ( tmaxAttenuation - tmin) )  * LSpectrumOut0;

  // get the light which has been scattered in the eye direction by the volume 

     if ( isPointLight) { // point light
       vec3 lightPos = lightPosOrDir;

       for (float t = tmin; t <= tmaxAttenuation; t+= deltaT  ) { 
  
    
	    vec3 currPos =  vec3(0,0,0) +  rayDir * t;

		vec3 dirFromPointLight =  currPos - lightPos;

	    dirFromPointLight = normalize( dirFromPointLight );

		
		// The light with given lambda is scattered in the direction of -rayDir with varying amount
		// depending on the scattering angle between the light direction and -rayDir, which varies
		// with each t. Consider only the scattering angles  within the range [ thetaStart, thetaEnd]. 
		//  The other scattering angles constributs less to the scattered intensity of lambda
		//   They are ignored for computational reasons. Some scattering spectrum is computed for
		// every ray direction (-rayDir)?? 
	  
	    float psi = acos ( dot ( -rayDir, dirFromPointLight) ) * 180.0 / pi;
			

		if (  psi >=  thetaStart  &&  psi <= thetaEnd    ) { 
		   
		   LSpectrumOut += exp( - sigma_e * ( t - tmin) )
		         * sigma_s * avgPhaseFunction(lambda, psi)
	               *  LSpectrumIn( LSpectrumSun, currPos, -dirFromPointLight, sigma_e) * deltaT;

		   
        }   //if 
       
       } // for 

       return LSpectrumOut; 

     } // if (point light)
	
     else { // directional light
      vec3 dirFromSun  = lightPosOrDir;
	  float psi = acos ( dot ( -rayDir, dirFromSun ) ) * 180.0 / pi;

	 // consider only the cases where the light is scattered in the direction of the eye
	 // with scattering angle within [ thetaStart, thetaEnd]. These scattering angles contricute
	 // to the rainbow and its neighborhood. The cases with other scattering
	 // angles contribute to the further region outside of the rainbow and its neighbor.
	 // These regions are not considered and the background color are assumed to dominate. 

      if (  psi >=  thetaStart  &&  psi <= thetaEnd    ) {	

        for (float t = tmaxAttenuation; t >= tmin; t-= deltaT  ) { 
      
	       vec3 currPos = vec3(0,0,0) +  rayDir * t;
	    		 
           LSpectrumOut += exp( - sigma_e * ( t - tmin) ) * sigma_s  
		           * avgPhaseFunction(lambda, psi)
	               *  LSpectrumIn(LSpectrumSun, currPos, -dirFromSun, sigma_e) * deltaT ;
        } // for
			  
      } // if

	  return LSpectrumOut; 
  
    } // else (parallel light) 	  

  } // else (isHit)

} // singleScatteringAndAttenuation

vec3 calculate_rainbowColor (bool isPointLight,  vec3 lightPosOrDir, vec3 rayDir, 
                             float sigma_s, float sigma_e, float tmin, float tmax) {
							

    float Isun,  lambda_m, lambda;
	int j;

	// For each wavelength lambda, compute the scattered light, which is attenuated along the direction 
	// -rayDir 

	/*
    for (lambda = lambdaStart, j = 0; j < nSpectralSamples; lambda += lambdaStep, j++) {
		
		lambda_m = lambda * 1.0e-9;				// convert lambda to meters from nanometers

        Isun = calculate_irradiance_of_sun(lambda_m); // isun: irradiance per nanometer

		
		float Lsun = Isun /  pi; // convert the irradiance to radiance 

        		
		//float avgPhasePerVolume  =  uDropDensity * particlePhase[j];
		        
       // irainbow[j] = ( ( Lsun* avgPhasePerVolume  / cos(phi_obs) ) / ( 1/cos(phi_sun) + 1/cos(phi_obs) ) )
	   //	                * ( 1 - exp(-tau_N * ( 1/cos(phi_sun) + 1/cos(phi_obs) ) ) );
         
        // the attenuation of the sun light while it passes through the drop volume IS considered, as well as
		// the attenuation of the scattered light

		
		*/

		//irainbow[j] = singleScatteringAndAttenuation(isPointLight, lambda, Lsun, sigma_s, sigma_e, lightPosOrDir, rayDir, tmin, tmax);
		
		ISpectrumSun = calculate_irradiance_of_sun(); // isun: irradiance per nanometer

		
		float[nSpectalSamples] LSpectrumSun = ISpectrumSun /  pi; // convert the irradiance to radiance 


		irainbow = singleScatteringAndAttenuation(isPointLight, LSpectrumSun, sigma_s, sigma_e, lightPosOrDir, rayDir, tmin, tmax);

		
	// convert spectral color to rgb color
	
	float X = 0, Y = 0, Z = 0;
	float XYZ;

	float xBar, yBar, zBar;

	// cie_colour_match[(lambda - 380) / 5][0] = xBar
    //  cie_colour_match[(lambda - 380) / 5][1] = yBar
    //  cie_colour_match[(lambda - 380) / 5][2] = zBar
	//  ==> cie_coulor_match [] from 380 nanometer to 780 nanometer
    
	
	// lambda: 400 ~ 700 into 120 steps with lamdbaStep - 2.5nm
	// cie_colour_match: from 380 to 780 nm


	for (lambda = lambdaStart, j =0; j < nSpectralSamples; lambda += lambdaStep, j++) {
		
		
		// cie_color_match[ (int (lambda) - 380) /  int( lambdaStep*2) ][k] 
		// => cie_color_match[ (int(lambaa) - 380) / int(lambdaStep*2) * 3 + k ] 
		// =	xBar( cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep*2)] [0]
		//         + cie_colour_match[ ( int(lambda +lambdaStep*2) - 380 ) / int( lambdaStep*2)] [0] ) / 2.0;

     /*
		int i1 =  ( int(lambda) - 380 ) / int( lambdaStep*2 );
		int i2 =  ( int(lambda + lambdaStep*2) - 380 ) / int( lambdaStep*2 );
		
		xBar = ( cie_colour_match[ i1 * 3  + 0] + cie_colour_match[ i2 * 3  + 0] ) / 2.0; 
		yBar = ( cie_colour_match[ i1 * 3 + 1 ] + cie_colour_match[ i2 * 3 + 1 ] ) /2.0;
		zBar = ( cie_colour_match[ i1 * 3 + 2 ] + cie_colour_match[ i2 * 3 + 2 ] ) /2.0;

		X += irainbow[j] * xBar * lambdaStep; 
		Y += irainbow[j] * yBar * lambdaStep;
		Z += irainbow[j] * zBar * lambdaStep;
	 */
	 
	 	
	   int i1 = int ( (lambda - 380 ) / ( lambdaStep*2 ) );
	   int i2 = int ( ( lambda + lambdaStep - 380) / (lambdaStep*2 ) );

	  
	   if ( i1 == i2) { 
		   xBar =  cie_colour_match[ i1 * 3 +  0]; 
		   yBar =  cie_colour_match[ i1 *3  + 1 ];
		   zBar =  cie_colour_match[ i1  *3 +  2 ];
		
		}  
	   else { // i2 = i1 + 1 => lambda is between two sample points 
		   
		  xBar = ( cie_colour_match[ i1 * 3 + 0] + cie_colour_match[ i2 *3 + 0] ) / 2.0; 
		  yBar = ( cie_colour_match[ i1 * 3 + 1 ] + cie_colour_match[ i2 *3 + 1 ] ) /2.0;
		  zBar = ( cie_colour_match[ i1 *3 +  2  ] + cie_colour_match[ i2 *3 + 2 ] ) /2.0;
        }

		X += irainbow[j] * xBar * lambdaStep; 
		Y += irainbow[j] * yBar * lambdaStep;
		Z += irainbow[j] * zBar * lambdaStep;
		

		
    }
	
	
	XYZ = (X+Y+Z);
	vec3 xyzColor = vec3( X/XYZ, Y/XYZ, Z/XYZ );

	
	vec3 rgbColor = xyz_to_rgb(g_cs, xyzColor);

	
/*
	if (constrain_rgb(rgbColor)) {
		norm_rgb(rgbColor);
		//output << rgbColor << "(Approximation)" << endl;
	} else {
		norm_rgb(rgbColor);
		//output << rgbColor << endl;
	}

 */
  constrain_rgb(rgbColor);

  return rgbColor;
  
}



// from pbrt book
float[nSpectralSamples] spectralColor(vec3 rgb) {
// primarySpectrums
    float[nSpectralSamples]  r;

        // Convert reflectance spectrum to RGB
        if (rgb[0] <= rgb[1] && rgb[0] <= rgb[2]) {
            // Compute reflectance _SampledSpectrum_ with _rgb[0]_ as minimum
            r += rgb[0] * rgbRefl2SpectWhite;
            if (rgb[1] <= rgb[2]) {
                r += (rgb[1] - rgb[0]) * rgbRefl2SpectCyan;
                r += (rgb[2] - rgb[1]) * rgbRefl2SpectBlue;
            }
            else {
                r += (rgb[2] - rgb[0]) * rgbRefl2SpectCyan;
                r += (rgb[1] - rgb[2]) * rgbRefl2SpectGreen;
            }
        }
        else if (rgb[1] <= rgb[0] && rgb[1] <= rgb[2]) {
            // Compute reflectance _SampledSpectrum_ with _rgb[1]_ as minimum
            r += rgb[1] * rgbRefl2SpectWhite;
            if (rgb[0] <= rgb[2]) {
                r += (rgb[0] - rgb[1]) * rgbRefl2SpectMagenta;
                r += (rgb[2] - rgb[0]) * rgbRefl2SpectBlue;
            }
            else {
                r += (rgb[2] - rgb[1]) * rgbRefl2SpectMagenta;
                r += (rgb[0] - rgb[2]) * rgbRefl2SpectRed;
            }
        }
        else {
            // Compute reflectance _SampledSpectrum_ with _rgb[2]_ as minimum
            r += rgb[2] * rgbRefl2SpectWhite;
            if (rgb[0] <= rgb[1]) {
                r += (rgb[0] - rgb[2]) * rgbRefl2SpectYellow;
                r += (rgb[1] - rgb[0]) * rgbRefl2SpectGreen;
            }
            else {
                r += (rgb[1] - rgb[2]) * rgbRefl2SpectYellow;
                r += (rgb[0] - rgb[1]) * rgbRefl2SpectRed;
            }
        }

        r *= .94;
		return r;

} // spectralColor

	
void main() {
  
  // get the pixel coordinates of the current fragment

 //	  ivec2 pixelPos = ivec2( gl_FragCoord.xy);

	  // get the surface color and depth of the current fragment
	 
//	  vec3  surfaceColor = vec3( texelFetch(uColorTex, pixelPos, 0) ); // color

//	  float zpixel  = texelFetch(uDepthTex, pixelPos, 0);  // z in non linear range [0,1]
     
      // conversion into NDC [-1,1]
 //     float zndc = zpixel * 2.0 - 1.0;

//	  float zeye = uProjectionMatrix[3][2]/ ( zndc - uProjectionMatrix[2][2] ); 

	  //fragColor = vec4(surfaceColor, 1);
	  //fragColor = vec4(1,0,0,1);

	 //return;


  vec3 diffuse, specular;
  vec3 volumeColor;

  float tmin, tmax;

      
  float scatCrossSection = computeScatCrossSection( uRadius,  lambdaStart);  // scatCrossSection is independent of
                                                                               // lambda, so use lambdaStart
    
  // uSunRayDir is  already conveted into eye space in the main program.
  
  vec3 rayDir = normalize( vPosition ); // the direction from pixel to fragment position vPosition

  
	     float sigma_s = uDropDensity * scatCrossSection;	// scattering coefficient [1/m]
	  	 float sigma_e = sigma_s; // assume that the water drop not absorb light, the extinction coeff = scat coeff. 

    // compute the surface color at vPosition and attenuated the color through
    // volume space between tmin and tmax.

	     bool isPointLight = false;
		 vec3 lightPosOrDir;

		 // compute the spectrum scattered in the direction of -rayDir and get the RGB color from the 
		 // spectrum

		 if ( isPointLight ) {
		      lightPosOrDir = uLight1Pos;  
		      volumeColor = calculate_rainbowColor ( isPointLight, lightPosOrDir,  rayDir, sigma_s, sigma_e, tmin, tmax); 
			  
			  //lightPosOrDir = uLight2Pos;
			  //volumeColor += calculate_rainbowColor ( isPointLight, lightPosOrDir, rayDir,sigma_s, sigma_e, tmin, tmax); 

		 }
		 else {
		      lightPosOrDir = -uSunRayDir;
		      volumeColor = calculate_rainbowColor ( isPointLight, lightPosOrDir, rayDir, sigma_s, sigma_e, tmin, tmax);  
			  		 		   
		 }

		fragColor  = vec4( volumeColor, 1);	
       
					
  

   // no volume color for rainbow:  The red emerges at a smaller angle psi than the violet
   //scattering is concentrated at theta_0 = 137.97 but also occurs for angles greater than that. 
   // This is what causes the white appearance `insidea the primary bow (and `outsidea the secondary bow).
   //Note the gap between 129 and 138. This is a region of negligible scattering (from higher orders) and appears 
   //as a dark space between the primary and secondary rainbows
   
	   
} // main

/* 

builtin functions:

radians(degrees)
degrees(radians)
sin(angle in radians)
asin, acos, pow(x,y), exp, log, sqrt, abs, sign, floor, ceil, mod, fract, mix(x,y,a)
length(vector x), distance(p0,p1), dot, cross, normalize, lessThan(vec,vec) equal(vec,vec)
*/

/*http://www.opengl.org/wiki/Compute_eye_space_from_window_space [good]
http://www.edxgraphics.com/blog/reconstructing-position-from-linear-depth
http://www.gamedev.net/topic/498615-how-can-i-get-the-world-coordinates-in-a-glsl-pixel-shader-by-using-the-depth-value/
http://mynameismjp.wordpress.com/2010/09/05/position-from-depth-3/
http://www.txutxi.com/?p=182 (GOOD)
vec4 ndcPos;
ndcPos.xy = ((2.0 * gl_FragCoord.xy) - (2.0 * viewport.xy)) / (viewport.zw) - 1;
ndcPos.z = (2.0 * gl_FragCoord.z - gl_DepthRange.near - gl_DepthRange.far) /
    (gl_DepthRange.far - gl_DepthRange.near);
ndcPos.w = 1.0;
 
vec4 clipPos = ndcPos / gl_FragCoord.w;
vec4 eyePos = invPersMatrix * clipPos;

**gl_FragCoord.w stores the reciprocal of the depth in view space (positive in front of the camera).

gl_FragCoord contains the window-space position of the current sample that this fragment represents. 
The Z component is the value that will be written to the depth buffer if the user does not override this (see below). 
The W component is special; it is 1/Wclip. That is, it is 1 divided by the W component of gl_Position output
 from the vertex or geometry shader.

 // // Screen size
uniform float screenWidth;
uniform float screenHeight;

...

// Converting (x,y,z) to range [0,1]
float x = gl_FragCoord.x/screenWidth;
float y = gl_FragCoord.y/screenHeight;
float z = gl_FragCoord.z; // Already in range [0,1]

// Converting from range [0,1] to NDC [-1,1]
float ndcx = x * 2.0 - 1.0;
float ndcy = y * 2.0 - 1.0;
float ndcz = z * 2.0 - 1.0;
vec3 ndc = vec3(ndcx, ndcy, ndcz);


 The fragment shader defines some uniform built-in values for the sake of convenience:


struct gl_DepthRangeParameters
{
    float near;
    float far;
    float diff;
};
uniform gl_DepthRangeParameters gl_DepthRange; 


This struct provides access to the glDepthRange near and far values.
 The diff value is the far value minus the near value.

	    
   // Converting (x,y,z) to range [0,1]
     float x = gl_FragCoord.x/screenWidth;
     float y = gl_FragCoord.y/screenHeight;
      float z = gl_FragCoord.z; // Already in range [0,1]

	  // texelFetch(uColorTexture, ivec2(gl_FragCoord.xy), 0);
	  // texture(uColorTexture, gl_FragCoord.xy/ textureSize(uColorTexture,0) )
	   


     // Converting from range [0,1] to NDC [-1,1]
     float ndcx = x * 2.0 - 1.0;
     float ndcy = y * 2.0 - 1.0;
     float ndcz = z * 2.0 - 1.0;
     vec3 ndc = vec3(ndcx, ndcy, ndcz);
	 

	 

HOW TO get Z values from Z buffer


Use glReadPixels with format = GL_DEPTH_COMPONENT, for example:
float depth;
glReadPixels(0, 0, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);


Will get the depth of pixel (0, 0).
 
gl_FragCoord.z is the depth value of the currently processed fragment,
 not the current value in the depth buffer. To get a depth texture, 
 create a texture with a depth format, render to an FBO with your depth texture 
 attached to GL_DEPTH_ATTACHMENT, and then bind the depth texture for sampling in your fragment shader. 
 1.Create and bind an FBO. Look up calls like glGenFramebuffers and glBindFramebuffer if you have not used FBOs before.
2.Create a texture or renderbuffer to be used as your color buffer, and attach it to the GL_COLOR_ATTACHMENT0 attachment point of your FBO with glFramebufferTexture2D or glFramebufferRenderbuffer. If you only care about the depth from this rendering pass, you can skip this and render without a color buffer.
3.Create a depth texture, and attach it to the GL_DEPTH_ATTACHMENT attachment point of the FBO.
4.Do your rendering that creates the depth you want to use.
5.Use glBindFramebuffer to switch back to the default framebuffer.
6.Bind your depth texture to a sampler used by your fragment shader.
7.Your fragment shader can now sample from the depth texture.

On question 2: gl_FragCoord.z is the depth value of the fragment that your shader is operating on, not the current value of the depth buffer at the fragment position

// Input texture
uniform sampler2D inputTexture;
// Frustum definition values
uniform float n; // znear
uniform float f; // zfar
uniform float l; // left
uniform float r; // right
uniform float b; // bottom
uniform float t; // top
// Quad/texture size in pixels
uniform float textureWidth;
uniform float textureHeight;

// z in non linear range [0,1]
float zpixel = inputTexture.r;

// conversion into NDC [-1,1]
float zndc = zpixel * 2.0 - 1.0;

// conversion into eye space
float zeye = 2*f*n / (zndc*(f-n)-(f+n));
// or zeye = uProjectionMatrix[3][2]/ ( zndc - uProjectionMatrix[2]2] ) 

// Converting from pixel coordinates to NDC
float xndc = gl_FragCoord.x/textureWidth * 2.0 - 1.0;
float yndc = gl_FragCoord.y/textureHeight * 2.0 - 1.0;

float xeye = -zeye*(xndc*(r-l)+(r+l))/(2.0*n);
float yeye = -zeye*(yndc*(t-b)+(t+b))/(2.0*n);

vec3 eyecoords = vec3(xeye, yeye, zeye);





	 */

