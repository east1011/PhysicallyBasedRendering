#version 130


// opengl 3.0, 2008 => GLSL 1.3
// Opengl 3.1, 2009 => GLSL 1.3 (ATI and NVIDiA)
// Opengl 3.2, 2009 => GLSL 1.5 ( //): Compatibility Profile: Openframework uses this version

// Opengl 3.3, 2010 => GLSL 3.3 (NVIDIA)

 // As of March, 2010: Opengl 4.0, GLSL 4.0
 // My computer: GL 4.0 support


#define lights_no 8
#define textures_no 8
#define M_PI       3.14159265358979323846f
#define INV_PI     0.31830988618379067154f
#define INV_TWOPI  0.15915494309189533577f
#define INV_FOURPI 0.07957747154594766788f
#define TEX_CONST 0
#define TEX_IMAGEMAP 1
#define LIGHT_POINT 0
#define LIGHT_DISTANT 1

in vec3 vNormal;
in vec3 vPosition;
in vec2 vTexCoord;
in mat3 vTBNMat;


out vec4 fragColor;


// A variable of sampler can only be defined in one of two ways. 
//It can be defined as a function parameter or as a uniform variable.
//uniform sampler2D texture1;
//void Function(in sampler2D myTexture);


uniform mat4 uEyeMatrix;
uniform mat4 uInvEyeMatrix;

uniform mat4 uModelMatrixAABB;
uniform mat4 uModelViewMatrix;


uniform sampler3D uDistTex;
uniform sampler3D uPhaseTex;
uniform sampler2D uScatTex;


// eyePos and aabb box for rainbow
//uniform vec3 uEyePos; 
// in eye space, so it should be (0,0,0): see below


uniform float uRadius;
uniform vec3  uSunRayDir;
uniform float uDropDensity;

// AABB bounding box
uniform vec3 uAABBmin; 
uniform vec3 uAABBmax;


out vec4 phaseSpectrum;
out vec4 rainbowColorSpectrum;


vec3 uEyePos = vec3(0,0,0);
float uTolerance = 1.0;

vec3[2]  uAABB = vec3[2](uAABBmin, uAABBmax); // This is OK in version 120

// parameters for rainbow computation

const  float pi = acos( -1.0);  // -1.0 const float in default; -1.0f = float

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
		

struct Ray {
  vec3 origin;
  vec3 direction;
  vec3 inv_direction;
  int sign[3];
};



struct Light {
       int lightType;  // LIGHT_POINT =0; LIGHT_DISTANT =1
	   vec3 intensity; // intensity for point light
	   vec3 L; // radiance for directional light
        vec3 lightPos;
        vec3 lightDir;
		//float cutoff;
		//float exponent;
		};

//struct UberMaterial {  
//       int KdTextureClass; # TEX_CONST =1, TEX_IMAGEMAP =2
//       vec3  KdConstColor; 
//	    sampler2D KdTextureUnit;  // sampler is not allowed in fields
//	   vec3  diffuse;         
//}

struct PlasticMaterial {  
       int KdTextureClass; // TEX_CONST =1, TEX_IMAGEMAP =2
       int KdTextureUnit;  // index to the constTextures or the textures
	   int KsTextureClass; // TEX_CONST =1, TEX_IMAGEMAP =2
       int KsTextureUnit;  // index to the constTextures or the textures
	   int roughnessTextureClass; // TEX_CONST =1, TEX_IMAGEMAP =2
       int roughnessTextureUnit;  // index to the constTextures or the textures
	   int bumpMapTextureClass; // TEX_CONST =1, TEX_IMAGEMAP =2
       int bumpMapTextureUnit;  // index to the constTextures or the textures
};

uniform int uNumOfLights;
uniform sampler2D uTextures[ textures_no]; 
uniform vec3  uConstTextures[ textures_no ];


uniform PlasticMaterial uPlasticMaterial;
uniform Light uLights [lights_no];



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



// mat3x2 (vec2, vec2, vec2): in column-major order
// mat3(float, float, float, // first column
//     float, float, float, // second column
//     float, float, float); // third column

// m[1] = the second column

//calculates tangent space matrix from normal, vector in plane and texture coordinates
//needs normalized normal and position!

mat3 computeTangentFrame(vec3 normal, vec3 position, vec2 texCoord)

{
// To choose T axis, the standard method is to orient the tangent vector in the same direction
// as the texture coordinates.
// deltaPos1 = deltaUV1.x * T + deltaUV1.y * B
// detaPos2 = deltaUV2.x * T  + deltaUV2.y * B
// Plug in known values of deltaPos1, 2 and deltaUV1, 2, and create six equations for T and B


// TBN: from tangent space to model space.
// invTBN = tranpose(TBN): from model space to tangent space

 // compute the derivative of the position coordinates
    vec3 dpx = dFdx(position);
    vec3 dpy = dFdy(position);

 // compute the derivative of the texture coordinates
    vec2 dtx = dFdx(texCoord);
    vec2 dty = dFdy(texCoord);
    
	// solve the linear system
	// mat3 M = mat3( dpx, dpy, cross(dpx, dpy) );
	// mat2 invM = inverse( M );
	// vec3 tangent = normalize( invM * vec3( dtx.t, dty.t, 0.0) );
	// vec3 binormal = normalize( invM * vec3( dtx.s, dty.s, 0.0) );
	 
    vec3 tangent = normalize(dpx * dty.t - dpy * dtx.t);
	vec3 binormal = normalize(-dpx * dty.s + dpy * dtx.s);
   
    return mat3(tangent, binormal, normal);
}




struct BSDF {
  vec3 nn; // shading normal
  vec3 ng;  // geometric normal
  vec3 sn;   // tangent = dp/du
  vec3 tn;   // binormal = cross(nn, sn)
 
};


struct  Fresnel { // FresnelDielectric
  float eta_i;
  float eta_t;
};
  

struct  SpecularTransmission {
   vec3 T; // transmission scale factor
   float etai;
   float etat;
   Fresnel fresnel; // fresnel(ei, et)
};

struct Lambertian { // BSDF_REFLECTION | BSDF_DIFFUSE
   vec3 R;  // kd (reflectance) is assigned to R
};
 
 
struct Blinn  { // Blinn is a  MicrofaceDistribution
  float exponent;

};

struct Microfacet {
 vec3 R;  // ks is assigned to R (reflectance)
  Fresnel fresnel; // FresnelDielectric
  Blinn  distribution; //  Blinn((1.f / rough) ) => exponent = 1.0/rough;

};

float AbsCosTheta( vec3 w) { return abs(w.z); }

struct SpecularReflection {
 
  vec3 R; // reflection scaling factor
  Fresnel fresnel;
};


float distribution_D( Blinn dist, vec3 wh ) {

   float costhetah = AbsCosTheta( wh );
   return (dist.exponent+2) * INV_TWOPI * pow(costhetah, dist.exponent);
}


	
vec3 evaluteKd( PlasticMaterial uPlasticMaterial, int KdTextureClass, vec2 vTexCoord) {

 if ( KdTextureClass == TEX_CONST ) {
       return uConstTextures[ uPlasticMaterial.KdTextureUnit ];
  }
  else {
      return vec3( texture( uTextures[ uPlasticMaterial.KdTextureUnit ], vTexCoord) );
  }
}


	
vec3 evaluateKs( PlasticMaterial uPlasticMaterial, int KsTextureClass, vec2 vTexCoord) {
  if ( KsTextureClass == TEX_CONST ) {
       return uConstTextures[ uPlasticMaterial.KsTextureUnit ];
  }
  else {
      return vec3( texture( uTextures[ uPlasticMaterial.KsTextureUnit ], vTexCoord) );
  }
}
	

float  evaluateRoughness (PlasticMaterial uPlasticMaterial, int roughnessTextureClass, vec2 vTexCoord ) {

  if ( roughnessTextureClass == TEX_CONST ) {
       return uConstTextures[ uPlasticMaterial.roughnessTextureUnit ][0]; // return the first component of the 
	                                                                   // constant texture value
  }
  else {
      return texture( uTextures[ uPlasticMaterial.roughnessTextureUnit ], vTexCoord);
  }
}




// BxDF Utility Functions


float CosTheta( vec3 w) { return w.z; }

float SinTheta2( vec3 w) {
    return max(0.0, 1.0 - CosTheta(w)*CosTheta(w));
}


float SinTheta( vec3 w) {
    return sqrt (SinTheta2(w));
}


float CosPhi( vec3 w) {
    float sintheta = SinTheta(w);
    if (sintheta == 0.0) return 1.0;
    return clamp(w.x / sintheta, -1.0, 1.0);
}


float SinPhi( vec3 w) {
    float sintheta = SinTheta(w);
    if (sintheta == 0.0) return 0.0;
    return clamp(w.y / sintheta, -1.0, 1.0);
}


bool SameHemisphere( vec3 w, vec3 wp) {
    return w.z * wp.z > 0.0;
}



vec3  FrDiel(float cosi, float cost,  vec3 etai, vec3 etat) {
    vec3  Rparl = ((etat * cosi) - (etai * cost)) /
                     ((etat * cosi) + (etai * cost));
    vec3  Rperp = ((etai * cosi) - (etat * cost)) /
                     ((etai * cosi) + (etat * cost));
    return (Rparl*Rparl + Rperp*Rperp) / 2.0;
}


vec3 evaluateFresnel( Fresnel fresnel, float cosi) {
 // Compute Fresnel reflectance for dielectric
    cosi = clamp(cosi, -1.0, 1.0);

    // Compute indices of refraction for dielectric
    bool entering = cosi > 0.0;

    float ei = fresnel.eta_i, et = fresnel.eta_t;

    if (!entering) {
       // swap(ei, et);
	   float temp = et;
	   et = ei;
	   ei = temp;  // ei is the index of refraction of the incident medium
	}

    // Compute _sint_ using Snell's law
    float sint = ei/et * sqrt(max(0.0, 1.0 - cosi*cosi));
    if (sint >= 1.) {
        // Handle total internal reflection
        return vec3(1.0);
    }
    else {
        float cost = sqrt(max(0.0, 1.0 - sint*sint));
        return FrDiel( abs(cosi), cost, vec3(ei), vec3(et) );
    }

}
	

vec3 SpecularReflection_Sample_f( SpecularReflection specularReflection, vec3 wo, out vec3 wi  ) {
  wi = vec3( -wo.x, -wo.y, wo.z );

  return evaluateFresnel ( specularReflection.fresnel, CosTheta(wo) ) * specularReflection.R / AbsCosTheta(wi);

}

vec3 SpecularTransmission_Sample_f( SpecularTransmission specularTransmission , vec3 wo, out vec3 wi ) {

 // Figure out which $\eta$ is incident and which is transmitted
    bool entering = CosTheta(wo) > 0.0;

    float ei = specularTransmission.etai; float et = specularTransmission.etat;
    if (!entering) {
        //swap(ei, et);

		float temp;
		temp = et;
		et = ei;
		ei = temp;

		}

    // Compute transmitted ray direction
    float sini2 = SinTheta2(wo);
    float eta = ei / et;
    float sint2 = eta * eta * sini2;

    // Handle total internal reflection for transmission
    if (sint2 >= 1.0) return 0.0;
    float cost = sqrt(max(0.0, 1.0 - sint2));
    if (entering) cost = -cost;
    float sintOverSini = eta;
    
	wi = vec3( sintOverSini * -wo.x, sintOverSini * -wo.y, cost);
	
	    
    vec3 F = evaluateFresnel( specularTransmission.fresnel, CosTheta(wo) );
    return   (ei*ei)/(et*et) * (vec3(1.0) - F ) * specularTransmission.T / AbsCosTheta(wi);

}
 // T=  kt = transmission scale factor

 vec3 Lambertian_f  ( Lambertian lambertian, vec3 wo, vec3 wi ) {
  return lambertian.R * INV_PI;
}


float Microfacet_G( vec3 wo, vec3 wi, vec3 wh ) {

        float NdotWh = AbsCosTheta( wh);
        float NdotWo = AbsCosTheta( wo);
        float NdotWi = AbsCosTheta( wi);
        float WOdotWh = abs( dot(wo, wh));
        return min(1.0, min((2.0 * NdotWh * NdotWo / WOdotWh),
                            (2.0 * NdotWh * NdotWi / WOdotWh)));
    

}

vec3  Microfacet_f ( Microfacet  microFacet, vec3 wo, vec3 wi ) {
    float cosThetaO = AbsCosTheta( wo );
    float cosThetaI = AbsCosTheta( wi );
    if (cosThetaI == 0.0 || cosThetaO == 0.0) return vec3(0.0);
    vec3  wh = normalize( wi + wo ); // half-angle vector
    if (wh.x == 0. && wh.y == 0. && wh.z == 0.) return vec3(0.0);
 
    float cosThetaH = dot(wi, wh);
    vec3  F = evaluateFresnel( microFacet.fresnel, cosThetaH);

    return microFacet.R * distribution_D(microFacet.distribution, wh) * Microfacet_G(wo, wi, wh) * F /
               (4.0 * cosThetaI * cosThetaO);


}





vec3 WorldToLocal(BSDF bsdf,  vec3 v) {
       return vec3 (dot(v, bsdf.sn), dot(v, bsdf.tn), dot(v, bsdf.nn));

}

vec3 LocalToWorld (BSDF bsdf,  vec3 v) {
       return vec3 ( bsdf.sn.x * v.x + bsdf.tn.x * v.y + bsdf.nn.x * v.z,
                      bsdf.sn.y * v.x + bsdf.tn.y * v.y + bsdf.nn.y * v.z,
                      bsdf.sn.z * v.x + bsdf.tn.z * v.y + bsdf.nn.z * v.z);

}

vec3  bsdf_f( BSDF bsdf, vec3 vPos, vec3 vNormal, vec2 vTexCoord, vec3 woW, vec3 wiW) {

    vec3 wi = WorldToLocal( bsdf, wiW );
	vec3 wo = WorldToLocal( bsdf, woW);


	
    vec3 f  = vec3(0.0, 0.0, 0.0);

	
	
    vec3 kd =  evaluteKd( uPlasticMaterial, uPlasticMaterial.KdTextureClass, vTexCoord);
	
    if (kd != vec3(0,0,0) ) {
        Lambertian lambertian = Lambertian( kd );
		f += Lambertian_f  ( lambertian,  wo, wi );
    }

 	
    vec3 ks =  evaluateKs( uPlasticMaterial, uPlasticMaterial.KsTextureClass, vTexCoord);
	
    if (ks != vec3(0,0,0) ) {
        Fresnel fresnel =  Fresnel(1.0, 1.0);
        float rough = evaluateRoughness (uPlasticMaterial, uPlasticMaterial.roughnessTextureClass, vTexCoord );
        Microfacet microFacet  =  Microfacet( ks, fresnel,  Blinn( 1.f / rough ) );
        f +=  Microfacet_f ( microFacet, wo, wi );
    }

	

	return f;

} // bsdf_f()

vec3  Sample_L( Light light, vec3 vPos, out vec3 wi ) {

if (light.lightType == LIGHT_POINT ) {

  // change the coordinates of light position to the view space
    
	vec3 lightPosInViewSpace = vec3( uInvEyeMatrix * vec4(light.lightPos, 1.0) );

    wi = normalize( lightPosInViewSpace - vPos);  // vPos is in view space
	float dist = distance( lightPosInViewSpace, vPos );

    //return light.intensity /  ( dist * dist  );
	return light.intensity; 
}
else if (light.lightType == LIGHT_DISTANT ) {
   wi = vec3( uInvEyeMatrix * vec4( light.lightDir, 0) );
   return light.L;
}
else {
  wi = vec3(0,0,0);

  return vec3(0,0,0);
  }



} // Sample_L()

vec3 estimateDirect( Light light, BSDF bsdf, vec3 vPos,  vec3 vNormal, vec2 vTexCoord,  vec3 wo ) {

    vec3  wi; // reflection dir; 
	          // wo = direction to the viewer =  -raydir = -normalize( vec3(vPos) );
	vec3  raydir = vPos; 
	//wi =  raydir - vNormal * 2 * dot(raydir, vNormal );

	//wi = normalize( wi ); 
	
	// if the sphere is also transparent compute refraction ray (transmission)
	//	if (sphere->transparency) {
	//		T ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
	//		T cosi = -nhit.dot(raydir);
	//		T k = 1 - eta * eta * (1 - cosi * cosi);
	//		Vec3<T> refrdir = raydir * eta + nhit * (eta *  cosi - sqrt(k));

	//		refrdir.normalize();
     //			}

        
    // Add light's contribution to reflected radiance
    // Get the light energy the current light

	vec3 Li = Sample_L(light,  vPos, wi ); // out wi


	vec3 f = bsdf_f( bsdf,  vPos, vNormal, vTexCoord, wo, wi);

    vec3 L =  f * Li *  abs( dot(wi, vNormal) );
         
    return L;
}

// Lo(p, 企o) = Le(p, 企o) + int [f (p, 企o, 企i) Li(p, 企i) |cos 亥i | d企i. 
// We evaluate  Li(p, 企i) for directions to light sources and for the directions
// of perfect reflection and refraction. wi = the reflection direction by  the ray
// coming from the camera
// Spaw a tree of rays starting from each pixel on the camera, and compute the contribution
// of each ray with its weight.

// We consider only specular surfaces like glass, mirrors, and water; we do not  account for other types of indirect lighting
// effects like light bouncing off a wall and illuminating a room.


vec3 estimateDirectDiffuse( Light light, vec3 vPosition, vec3 vNormal,  vec2 vTexCoord, vec3 toV ) {   	

   vec3 lightDir;

   vec3 Li = Sample_L(light,  vPosition, lightDir ); // out lightDir

   if ( lightDir == vec3(0,0,0) ) {
      return vec4(0,0,0,1);
	}


   // compute diffuse color
   float diffuseCoefficient = max( 0.0, dot(vNormal, lightDir) ) ;

    vec3 kd =  evaluteKd( uPlasticMaterial, uPlasticMaterial.KdTextureClass, vTexCoord);
	
   
    Lambertian lambertian = Lambertian( kd );
	vec3 f = kd * INV_PI;
   

   vec3 diffuseColor = Li  * diffuseCoefficient * f;
  
   // compute specular reflection color

  //vec3 toV = -normalize( vec3(vPosition) );
  vec3 h = normalize( toV + lightDir  );

   
  float specularCoefficient  = pow( max(0.0, dot(h, vNormal) ), 30.0 );

  vec3 whiteLight = vec3( 0.6, 0.6, 0.6);

  vec3  specularColor = whiteLight * specularCoefficient;
     
  vec3  ambientColor = vec3(0.1, 0.1, 0.1);


  vec3 combinedColor = ambientColor + diffuseColor + specularColor;

  
  
  return combinedColor;

  }



void main() {

// prepare BSDF
// create the tangent frame at the current interpolated position vPosition

//mat3 TBN = computeTangentFrame(vNormal, vPosition, vTexCoord);

BSDF bsdf = BSDF( vTBNMat[2], vTBNMat[2], vTBNMat[0], vTBNMat[1] );

//struct BSDF {
//  vec3 nn; // shading normal
//  vec3 ng;  // geometric normal
//  vec3 sn;   // tangent = dp/du
//  vec3 tn;   // binormal = cross(nn, sn)
// 
//};

vec3 L = vec3(0,0,0);

//////////////
vec3 kd =  evaluteKd( uPlasticMaterial, uPlasticMaterial.KdTextureClass, vTexCoord);
fragColor= vec4(kd,1);
return ;

////////////


for (int i = 0; i < uNumOfLights; i++ ) {

  Light light = uLights[i];
  //vec3 lightPos = light.lightPos;
  vec3  wo  = normalize( - vPosition ); 
 
 // vec3 toLight = normalize( lightPos - vPosition);
  
  vec3 vNormal = normalize(vNormal);

  L += estimateDirect( light, bsdf,  vPosition, vNormal,  vTexCoord, wo );  	
 //L += estimateDirectDiffuse( light,  vPosition, vNormal,  vTexCoord, wo );  	
}
						   
fragColor = vec4( L, 1.0);

}

