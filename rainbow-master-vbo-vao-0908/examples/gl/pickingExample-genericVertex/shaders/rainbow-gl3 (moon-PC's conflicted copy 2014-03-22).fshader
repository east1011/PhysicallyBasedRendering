
#version 130
// OpenGL defines a pixel's center, whose location is 0.5 pixels up and to the right of the lower left 
//corner of the pixel (which defines the pixel). In other words, pixel 0,0 goes from 0,0 to 1,1 
//and its sample-center is 0.5, 0.5.
// The GL will sample the texture for each pixel filled in by our polygon. So what texture coordinates are used per pixel?

//The texture coordinates (and all fragment information) is interpolated per fragment (pixel) at the pixel center.


#define MAX_ITERATIONS 100

/* http://http.developer.nvidia.com/GPUGems2/gpugems2_chapter27.html
So the first step in writing a filter shader is to understand how to compute the integer sample coordinate 
corresponding to the current texture coordinate and how to convert sample coordinates back into conventional
 texture coordinates. The following equations comply with Direct3D's texture-mapping specifications:

sampleCoord = floor(textureCoord * texSize);

textureCoord = (sampleCoord + 0.5) / texSize;

One additional tip is that it is not necessary to snap texture coordinates to the center of the nearest texel as long as texture lookups do not perform any type of filtering. This means setting filtering to GL_NEAREST in OpenGL or D3DTEXF_POINT in Direct3D.
 This way, coordinate snapping is performed automatically during texture lookups
*/

vec3 texIncrement = vec3( 1 / nRadii, 1 / nSpectralSamples, 1/ nThetas);
 
#define TC_XMINUS_YMINUS(coord) (coord - texIncrement.xy)

#define TC_XPLUS_YPLUS(coord) (coord + texIncrement.xy)

#define TC_XCENTER_YMINUS(coord) (coord - texIncrement.zy)

#define TC_XCENTER_YPLUS(coord) (coord + texIncrement.zy)

#define TC_XMINUS_YCENTER(coord) (coord - texIncrement.xz)

#define TC_XPLUS_YCENTER(coord) (coord + texIncrement.xz)



#define TC_XMINUS_YPLUS(coord) (coord - texIncrement.xw)

#define TC_XPLUS_YMINUS(coord) (coord + texIncrement.xw)

uniform float uWindowWidth; 
uniform float uWindowHeight; 

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

// Samplers are special types in OpenGL; they represent a texture
// that has been bound to the OpenGL context
// Arrays in GPU memory are called textures or texture samplers
// On the CPU, we usually talk about array indices, on the GPU, we will need texture coordinates to access values
// stored in the textures. Texture coordinates need to address texel centers. 
//Coordinates have to be normalized to the range [0,1] by [0,1], independent of the dimension [0,M] by [0,N] of the texture. 

//In the fragment shader, we have a new uniform variable of type sampler2D, texture, which basically represents a pointer to the texture data,with the added capacity of sampling it at arbitrary texture coordinates (i.e: not necessarily at the center of a texel). Depending on the sampling configuration of the texture (linear, bilinear, etc.),
// a different interpolation algorithm will be used by the GPU
// texelFetch (...) disables filtering, including mip mapping.

uniform sampler2D uTexUnit0;
uniform sampler2D uTexNormal;

//uniform sampler3D uDistTex;
uniform sampler3D uPhaseTex;
uniform sampler2D uScatTex;

// lights in eye space
uniform vec3 uLight1Pos;
uniform vec3 uLight2Pos;

// Sun direction in eye space

uniform vec3  uSunRayDir;
// eyePos and aabb box for rainbow
//uniform vec3 uEyePos; 
// in eye space, so it should be (0,0,0): see below


uniform vec4 uLight1Color, uLight2Color;

uniform vec4 uMaterialColor;
uniform vec3 uLightRGBColor;

uniform float uRadius;

uniform float uDropDensity;

// AABB bounding box in eyespace
uniform vec3 uAABBmin; 
uniform vec3 uAABBmax;


in vec2 vTexCoord;
in vec3 vNormal;
in vec3 vPosition;
in   mat3 vNTMat;

out vec4 fragColor;

/* The ways arrays are used: 
float[5] foo() { }
as a constructor of an array
float[5](3.4, 4.2, 5.0, 5.2, 1.1)

void foo(float[5])
float[5] a;

Arrays can have initializers formed from array constructors:
float a[5] = float[5](3.4, 4.2, 5.0, 5.2, 1.1);
float a[5] = float[](3.4, 4.2, 5.0, 5.2, 1.1); // same thing
*/ 



vec3 uEyePos = vec3(0,0,0);
float uTolerance = 1.0;

vec3[2]  uAABB = vec3[2](uAABBmin, uAABBmax); // This is OK in version 120

// parameters for rainbow computation

const  float pi = acos( -1.0);  // -1.0 const float in default; -1.0f = float

const int nSpectralSampleSteps = 120;
const int nThetaSteps = 100; // 100 steps between theta 130 deg and theta 142;

const int nRadiusSteps  = 2;

const int nRadii = nRadiusSteps +1;
const int nSpectralSamples = nSpectralSampleSteps +1;

const int nThetas = nThetaSteps +1; // 100 steps between theta 130 deg and theta 142;

float irainbow[ nSpectralSamples ];
float particlePhase[nSpectralSamples ];


const float  lambdaStart = 400;   // 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
const float lambdaEnd = 700;
     
const  float radiusStart =  1.0e-3;  // 
//const  float radiusEnd = 2.0e-3;           // 2 mm = 2 * e-3 m
const  float radiusEnd = 1.5e-3;           // 2 mm = 2 * e-3 m

const	float thetaStart = 130.0;
const	float thetaEnd = 142.0;

const float lambdaStep = ( lambdaEnd - lambdaStart) / float( nSpectralSampleSteps );
const float thetaStep = ( thetaEnd - thetaStart) / float (nThetaSteps);
const float radiusStep = (radiusEnd - radiusStart) / float (nRadiusSteps);

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

const float textureSize = 512.0; //size of the texture
const float texelSize = 1.0 / textureSize; //size of one texel 
 
vec4 texture2DBilinear( sampler2D textureSampler, vec2 uv )
{
    // in vertex shaders you should use texture2DLod instead of texture2D
    vec4 tl = texture2D(textureSampler, uv);
    vec4 tr = texture2D(textureSampler, uv + vec2(texelSize, 0));
    vec4 bl = texture2D(textureSampler, uv + vec2(0, texelSize));
    vec4 br = texture2D(textureSampler, uv + vec2(texelSize , texelSize));
    vec2 f = fract( uv.xy * textureSize ); // get the decimal part
    vec4 tA = mix( tl, tr, f.x ); // will interpolate the red dot in the image
    vec4 tB = mix( bl, br, f.x ); // will interpolate the blue dot in the image
    return mix( tA, tB, f.y ); // will interpolate the green dot in the image
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


    
int intersection_distance_no_if ( 

   in Ray ray,   in vec3[2] aabb,
   out float tmin, out float tmax,  out vec3 normalAtFront )
   
 {
    
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

float computeScatCrossSection( float radius, float lambda) {

	// get the index to the texture for (radius, lambda)

	vec2 index = vec2( (radius  - radiusStart) / (radiusEnd -  radiusStart), 
	               (lambda - lambdaStart) / ( lambdaEnd - lambdaStart) );
	
	float scatCrossSection = texture( uScatTex, index).r;

    return scatCrossSection;
	
}


float computePhaseFunction( float radius, float lambda, float theta) {
  	
	// get the index to the texture for (radius, lambda)
	// The texture function for 1D textures expects the texture coordinate to be normalized

	vec3 index = vec3 ( (radius  - radiusStart) / (radiusEnd - radiusStart), 
	               (lambda - lambdaStart) / ( lambdaEnd -  lambdaStart),
				   (theta - thetaStart ) / ( thetaEnd -  thetaStart) ) ;
	
	
	float phase  = texture( uPhaseTex, index).r;
	
    return phase;

}


float calculate_irradiance_of_sun (float lambda) {
    
    float R_Earth = 1.496e+8;
    float R_Sun = 6.95e+5;
    float h = 6.6261e-34;      // Planck's constant
    float c = 2.9979e+8;       // speed of light in vacuo
    float k = 1.3806e-23;      // Boltzmann's constant
    float T = 5782;            // Sun's temperature in Kelvin
    float Isun;
    
    Isun = pow( (R_Sun/R_Earth),2.0 )  * 2*pi*h*pow(c,2.0)  / ( pow(lambda, 5.0) * ( exp( h*c/ (lambda*k*T) ) -1 ) ) * 1.0e-9  ;
	   //Moon:  1.0e-9 is multiplied to convert the irradiance per meter to the irradiance per nanometer
	   // Original had 1.0e-10. 
    
    return Isun;

}

vec3 calculate_rainbowColor (float radius, float uDropDensity,  float phi_obs,float phi_sun,
                             float psi, float tau_N) {

    float Isun;
    
    float lambda_m;

   	float lambda;
	int j;

	

    for (lambda = lambdaStart, j = 0; j < nSpectralSamples; lambda += lambdaStep, j++) {
		
		lambda_m = lambda * 1.0e-9;				// convert lambda to meters from nanometers

        Isun = calculate_irradiance_of_sun(lambda_m); // isun: irradiance per nanometer

		// ignore the attenuation of the sun light while it passes through the drop volume

		float Lsun = Isun /  pi; // convert the irradiance to radiance 

        particlePhase[j] = computePhaseFunction( radius,  lambda, psi); 
		

		//float particlePhase = 1.230634e-002;

		float avgPhasePerVolume  =  uDropDensity * particlePhase[j];
		        
        irainbow[j] = ( ( Lsun* avgPhasePerVolume  / cos(phi_obs) ) / ( 1/cos(phi_sun) + 1/cos(phi_obs) ) )
		                * ( 1 - exp(-tau_N * ( 1/cos(phi_sun) + 1/cos(phi_obs) ) ) );
         
        

	}


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


		int i1 =  ( int(lambda) - 380 ) / int( lambdaStep*2 );
		int i2 =  ( int(lambda + lambdaStep*2) - 380 ) / int( lambdaStep*2 );
		
		xBar = ( cie_colour_match[ i1 * 3  + 0] + cie_colour_match[ i2 * 3  + 0] ) / 2.0; 
		yBar = ( cie_colour_match[ i1 * 3 + 1 ] + cie_colour_match[ i2 * 3 + 1 ] ) /2.0;
		zBar = ( cie_colour_match[ i1 * 3 + 2 ] + cie_colour_match[ i2 * 3 + 2 ] ) /2.0;

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





void main() {
  
  vec3 normalAtFront;
  vec3 diffuse, specular;

  float phi_obs, phi_sun;
  float tmin, tmax;


/*	
	float x = gl_FragCoord.x/uWindowWidth;
    float y = gl_FragCoord.y/uWindowHeight;

    float z = gl_FragCoord.z; // Already in range [0,1]

vPosition =(x,y,z,w): vPosition.xyz / vPosition.w = (x_n, y_n, z_n) is the coordinate value 
 in normalized device coordinates (NDC). This value ranges from -1 to 1 in each dimension.

 viewport transformation: [x_w, y_w, z_w, 1] = [ ] * [x_n, y_n, z_n, 1]==> z_w = 1/2 * z_n + 1/2 ranges from 0 to 1.
 // (x_w, y_w) ranges from ( -0.5, -0.5) to ( W-0.5, H-0.5).

gl_FragCoord gives the window coordinates. Dividing (gl_FragCoord.x, gl_FragCoord.y) + 0.5  by width and height 
gives the position ranging from 0 to 1

*/

  
	
/*  
	
   //  Phase debug	
	float x = gl_FragCoord.x/uWindowWidth;
    float y = gl_FragCoord.y/uWindowHeight;

    float z = gl_FragCoord.z; // Already in range [0,1]
 
    float theta = thetaStart + (thetaEnd - thetaStart) * x;

	
	float tlambda;

	int tj;
	bool zeroPhase = false;

    for (tlambda  = lambdaStart, tj = 0; tj  < nSpectralSamples; tlambda += lambdaStep, tj++) {
        particlePhase[tj] = computePhaseFunction( uRadius,  tlambda, theta); //  
		
    }

	//normalize_positive_vector( particlePhase, nSpectralSamples );

	tlambda  = lambdaStart + (lambdaEnd - lambdaStart) * y;
	int lambdaIndex;

	//gl_FragColor = vec4(x,y,0,1);
	//return;

	lambdaIndex = int ( ( tlambda - lambdaStart) / lambdaStep  );

	
	if ( lambdaIndex < 60 ) gl_FragColor = vec4(1,0,0,1);
	else gl_FragColor = vec4(0,1,0,1);
	return;
	
	
	normalize_positive_vector( particlePhase, nSpectralSamples);
			
	float phase = particlePhase[ lambdaIndex ];

	
	if (phase < 1.0e-2 )  gl_FragColor = vec4(1,0,0,1); // phase func value > 1.0e-4, so no red color should
	                                                     // appear, but does.
	if ( phase >= 1.0e-2)   gl_FragColor = vec4( 0, 1,  0,1);

    return;
	
	
	
	gl_FragColor = vec4(0, phase, 0,1);
	return;

*/


  /*
  //  spectrum debug

    float Isun;
    
    float lambda_m;

 
   	float lambda;
	int j;

	float tau_N = 1.0;

	normalAtFront = vec3(0,0,1);

	
    vec3 trayDir = normalize( vPosition - uEyePos );
	Ray ray = makeRay( uEyePos, trayDir);


	phi_obs = acos( dot( normalAtFront, -trayDir) );       //  rayDir in eyeSpace
    phi_sun = acos( dot( normalAtFront, 
	                         uSunRayDir )  );


   // for debug	
	float x = gl_FragCoord.x/uWindowWidth;
    float y = gl_FragCoord.y/uWindowHeight;

    float z = gl_FragCoord.z; // Already in range [0,1]
 
    float theta = thetaStart + (thetaEnd - thetaStart) * x;

    if ( !  ( theta >= 130 && theta <= 142 )    )  {
	
	   fragColor = vec4(0,0,0,1);
	   return;
	}	 
	

    for (lambda = lambdaStart, j = 0; j < nSpectralSamples; lambda += lambdaStep, j++) {
		
		lambda_m = lambda * 1.0e-9;				// convert lambda to meters from nanometers

        Isun = calculate_irradiance_of_sun(lambda_m); // isun: irradiance per nanometer

		// ignore the attenuation of the sun light while it passes through the drop volume

		float Lsun = Isun /  pi; // convert the irradiance to radiance 

        particlePhase [j] = computePhaseFunction( uRadius,  lambda, theta); 
		

		//float particlePhase = 1.230634e-002;

		float avgPhasePerVolume  =  uDropDensity * particlePhase[j];
		        
        irainbow[j] = ( ( Lsun* avgPhasePerVolume  / cos(phi_obs) ) / ( 1/cos(phi_sun) + 1/cos(phi_obs) ) )
		                * ( 1 - exp(-tau_N * ( 1/cos(phi_sun) + 1/cos(phi_obs) ) ) );
         
        

	} // for

	
	
	normalize_positive_vector( irainbow, nSpectralSamples );

	
   lambda = lambdaStart + ( lambdaEnd - lambdaStart) * y;
	
	int lambdaIndex;

	lambdaIndex = int(  (lambda - lambdaStart) /  lambdaStep ) ; 
	
	float rainbow = irainbow[ lambdaIndex ];	
    fragColor =  vec4( 0, rainbow,  0,1);

    return;
	
	

	
	// display the color for a particular deflection angle

	// convert spectral color to rgb color
	
	float X = 0, Y = 0, Z = 0;
	float XYZ;

	float xBar, yBar, zBar;

	
	// lambda: 400 ~ 700 into 120 steps with lamdbaStep - 2.5nm
	// cie_colour_match: from 380 to 780 nm

	
	float epsilon = 1.0e-3;
	

	for (lambda = lambdaStart, j =0; j < nSpectralSamples; lambda += lambdaStep, j++) {
		
		
		// cie_color_match[ (int (lambda) - 380) /  int( lambdaStep*2) ][k] 
		// => cie_color_match[ (int(lambaa) - 380) / int(lambdaStep*2) * 3 + k ] 
		
		if (  definitelyLessThan ( fract ( lambda / (lambdaStep * 2) ), 
		                                   (2 * lambdaStep ) / 5.0, epsilon ) ) {
		 
		 // lambda is at a cie sample lambda which is sampled every lambdaStep * 2
		  
		   int i = int (  lambda - 380 ) / int ( lambdaStep*2 );

		   xBar =  cie_colour_match[ i * 3  + 0]; 
		   yBar =  cie_colour_match[ i * 3 + 1 ];
		   zBar =  cie_colour_match[ i * 3 + 2 ];
		
		}  
		else { // lambda is between two points and interpolate between them

		  int i1 =  int (  lambda - lambdaStep - 380 ) / int ( lambdaStep*2 );
		  int i2 =  int (  lambda + lambdaStep  - 380 ) / int ( lambdaStep*2 );
		
		  xBar = ( cie_colour_match[ i1 * 3  + 0] + cie_colour_match[ i2 * 3  + 0] ) / 2.0; 
		  yBar = ( cie_colour_match[ i1 * 3 + 1 ] + cie_colour_match[ i2 * 3 + 1 ] ) /2.0;
		  zBar = ( cie_colour_match[ i1 * 3 + 2 ] + cie_colour_match[ i2 * 3 + 2 ] ) /2.0;
        }
	
		X += irainbow[j] * xBar * lambdaStep; 
		Y += irainbow[j] * yBar * lambdaStep;
		Z += irainbow[j] * zBar * lambdaStep;
		

    } // for
	
	
	XYZ = (X+Y+Z);
	vec3 xyzColor = vec3( X/XYZ, Y/XYZ, Z/XYZ );

	
	vec3 rgbColor = xyz_to_rgb(g_cs, xyzColor);

 	

	if (constrain_rgb(rgbColor)) {
		norm_rgb(rgbColor);
		//output << rgbColor << "(Approximation)" << endl;
	} else {
		norm_rgb(rgbColor);
		//output << rgbColor << endl;
	}


 // constrain_rgb(rgbColor);

  fragColor  = vec4( rgbColor, 1);

  return;
  */


  

 
// rainbow 
  float scatCrossSection = computeScatCrossSection( uRadius,  lambdaStart);
  
    // current position along the ray

  vec3 rayDir = normalize( vPosition - uEyePos );
      
  float psi = acos ( dot ( -rayDir, -uSunRayDir ) ) * 180.0 / pi;

  if (  psi >=  138 &&  psi <= 140    ) { 
  
    
      Ray ray = makeRay( uEyePos, rayDir);
	  
	   
	 	 
	 int hitType = intersection_distance_no_if( ray, uAABB, tmin, tmax, normalAtFront );	

      if ( hitType == 0 ) { // the AABB  is not hit => do not paint the current fragment
	  
	       fragColor   = vec4(0,1,0,1);
		   return;

           //discard; 
      }

	  else {
	  
	    // hit the surface and compute the volume color from the "hit point" to the eye
		
	   
       // tmin is the distance to the hit surface or the back face of the AABB

	  // get the angle between the sun ray and the front face of the AABB
      // get the angle between the observer ray and the front face of the AABB
    
        phi_obs = acos( dot( normalAtFront, -rayDir) );       //  rayDir in eyeSpace
        phi_sun = acos( dot( normalAtFront, 
	                         uSunRayDir )  );

	  // compute the optical depth normal to the AABB

	    float rayPathNormal =  (tmax - tmin) * cos (phi_obs);
	    float tau_N = uDropDensity * scatCrossSection  * rayPathNormal;	
	  	
    // compute the surface color at vPosition and attenuated the color through
    // volume space between tmin and tmax.

	     tau_N = 1.0;

         vec3 volumeColor = calculate_rainbowColor ( uRadius, uDropDensity, phi_obs,  phi_sun, psi, tau_N);  
	    
		 // with psi different, volume color should be different. why not

		 fragColor  = vec4( volumeColor, 1);	
       
		 return;
	  

        
	   }  // else

	 } // if  

  else { 



   // no volume color for rainbow:  The red emerges at a smaller angle psi than the violet
   //scattering is concentrated at theta_0 = 137.97 but also occurs for angles greater than that. 
   // This is what causes the white appearance `insidea the primary bow (and `outsidea the secondary bow).
   //Note the gap between 129 and 138. This is a region of negligible scattering (from higher orders) and appears 
   //as a dark space between the primary and secondary rainbows

	fragColor = vec4(0,0,0,1); 
	return;
	 
	// causes the fragment shader to stop execution, without writing anything (including depth) 
	// to the output buffer.

  }
  

  

	   
} // main

/* 

builtin functions:

radians(degrees)
degrees(radians)
sin(angle in radians)
asin, acos, pow(x,y), exp, log, sqrt, abs, sign, floor, ceil, mod, fract, mix(x,y,a)
length(vector x), distance(p0,p1), dot, cross, normalize, lessThan(vec,vec) equal(vec,vec)
*/
/*
	    
   // Converting (x,y,z) to range [0,1]
     float x = gl_FragCoord.x/screenWidth;
     float y = gl_FragCoord.y/screenHeight;
      float z = gl_FragCoord.z; // Already in range [0,1]

     // Converting from range [0,1] to NDC [-1,1]
     float ndcx = x * 2.0 - 1.0;
     float ndcy = y * 2.0 - 1.0;
     float ndcz = z * 2.0 - 1.0;
     vec3 ndc = vec3(ndcx, ndcy, ndcz);
	 
	 */
