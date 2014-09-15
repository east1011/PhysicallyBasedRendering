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

Ray makeRay( vec3 origin, vec3 direction) {

 vec3 inv_direction = vec3(1.0) / direction;

 return
    Ray( 
      origin, direction, 
	  inv_direction,
	  { if (inv_direction.x < 0)  1 else 0,
	    if (inv_direction.y < 0)  1 else 0,
		if (inv_direction.z < 0)  1 else 0  }
	 );
	  
}

vec3 posInAABB (vec3 posInEyeSpace, vec2[2] aabb) {

// ModelMatrixAABB * localPos = uEyeMatrix * posInEyeSpace;
 // localPos = inv(uModelMatrixAABB) * uEyeMattrix * posInEyeSpace 


	
	float AABBwidth = aabb[1].x - aabb[0].x;
	float AABBheight = aabb[1].y - aabb[0].y;
	float AABBdepth = aabb[1].z - aabb[0].z;

	vec3 localPoss = vec3( inverse( uModelMatrixAABB ) * uEyeMatrix * vec4(posInEyeSpace,0);

	vec3 scaledPos = vec3( localPos[0]  / AABBwidth,
	                        localPos[1] / AABBheight,
		    	            localPos[2] / AABBdepth );

	return scaledPos;
	
}


    
bool intersection_distance_no_if ( 

   in Ray ray,   in vec3[2] aabb,
   out float tmin, out tmax,  out vec3 normalAtFront )  {
 // ModelMatAABB * localPos = uEyeMatrix * posInEyeSpace;
 // localPos = inv(uModelMatrixAABB) * uEyeMatrix * posInEyeSpace  
   
   tmin = ( aabb[ ray.sign[0] ].x - ray.origin.x ) * ray.inv_direction.x;
   tmax = ( aabb[1- ray.sign[0] ].x - ray.origin.x ) * ray.inv_direction.x;
   tymin = ( aabb[ ray.sign[1]  ].y - ray.origin.y ) * ray.inv_direction.y;
   tymax = ( aabb[ 1 - ray.sign[1]  ].y - ray.origin.y ) * ray.inv_direction.y;
   tzmin = ( aabb[ ray.sign[2] ].z - ray.origin.z ) * ray.inv_direction.z;
   tzmax = ( aabb[ 1- ray.sign[2] ].z - ray.origin.z ) * ray.inv_direction.z;
   tmin = max( max( tmin, tymin), tzmin );
   tmax = min( min( tmax, tymax), tzmax );


  // compute the normal at the front face relative to the eye, where
  // the hit position is at tmin along the ray direction.
  // Assume that the front face is hit for now  	  
    
	normalAtFront = vec3( 0,0,1); 

	float t = tmin;

	// traverse the distance ﬁeld until we hit the surface or leave the AABB
	
	vec3 posInEyeSpace = ray.origin + ray.direction * tmin;

	for( int i = 0; i < MAX_ITERATIONS; ++i ) {

        //vec3 vb = query( level, pos );						// locate the cell’s voxel block
		//vec3 cp = toCellSpace( pos, level );				// pos into cell space
		//float d = reconstructDist( vb, cp ) - isovalue;

		//float d = getDist(pos) - isovalue;

		
		// Convert eyespace pos into the local position in AABB coordinate system;
		// The AABB axes are parallel to the world system (while they pass through the center of
		//  AABB) ;
		//  The conversion of pos into texture index requires it, 
		//  because distance fields are stored at positions which are relative to AABB 
		
		vec3 indexPos = posInAABB( posInEyeSpace,  aabb);
		

		float d = texture3D (uDistTex, indexPos).r;

        ///////////////////////////////////////////////

		t += d;

        // step along the ray
		posInEyeSpace = ray.origin + ray.direction * t;
        // update the current position
		
		if( d < uTolerance ) {
            // did we hit the isosurface?
			tmax = t;			 
			return true; // hit the surface
		 }

		
		if( t >= tmax  ) {
		   
		   return false; //  not hit the surface
		 } 
			
	} // endfor (MAX_ITERATION)
	
	
	return false; //  not hit the surface
	
 }

    
 


//Cvec3  xyz_to_rgb(struct colourSystem *cs,   const vec3 xyz) {

vec3  xyz_to_rgb(struct colourSystem cs,   const vec3 xyz) {
    double xr, yr, zr, xg, yg, zg, xb, yb, zb;
    double xw, yw, zw;
    double rx, ry, rz, gx, gy, gz, bx, by, bz;
    double rw, gw, bw;

	vec3 rgbColor;

    xr = cs.xRed;    yr = cs.yRed;    zr = 1 - (xr + yr);
    xg = cs.xGreen;  yg = cs.yGreen;  zg = 1 - (xg + yg);
    xb = cs.xBlue;   yb = cs.yBlue;   zb = 1 - (xb + yb);

    xw = cs.xWhite;  yw = cs.yWhite;  zw = 1 - (xw + yw);

    /* xyz -> rgb matrix, before scaling to white. */
    
    rx = (yg * zb) - (yb * zg);  ry = (xb * zg) - (xg * zb);  rz = (xg * yb) - (xb * yg);
    gx = (yb * zr) - (yr * zb);  gy = (xr * zb) - (xb * zr);  gz = (xb * yr) - (xr * yb);
    bx = (yr * zg) - (yg * zr);  by = (xg * zr) - (xr * zg);  bz = (xr * yg) - (xg * yr);

    /* White scaling factors.
       Dividing by yw scales the white luminance to unity, as conventional. */
       
    rw = ((rx * xw) + (ry * yw) + (rz * zw)) / yw;
    gw = ((gx * xw) + (gy * yw) + (gz * zw)) / yw;
    bw = ((bx * xw) + (by * yw) + (bz * zw)) / yw;

    /* xyz -> rgb matrix, correctly scaled to white. */
    
    rx = rx / rw;  ry = ry / rw;  rz = rz / rw;
    gx = gx / gw;  gy = gy / gw;  gz = gz / gw;
    bx = bx / bw;  by = by / bw;  bz = bz / bw;

    /* rgb of the desired point */
    
	rgbColor = vec3((rx * xyz[0]) + (ry * xyz[1]) + (rz * xyz[2]), 
				(gx * xyz[0]) + (gy * xyz[1]) + (gz * xyz[2]), 
				(bx * xyz[0]) + (by * xyz[1]) + (bz * xyz[2]));

	return rgbColor;
}



int constrain_rgb(inout vec3 rgb)
{
    double w;

    /* Amount of white needed is w = - min(0, *r, *g, *b) */
    
	w = if (0 < rgb[0])  0 else  rgb[0];
	w = if (w < rgb[1])  w  else  rgb[1];
	w = if (w < rgb[2])  w  else rgb[2];
    w = -w;

    /* Add just enough white to make r, g, b all positive. */
    
    if (w > 0) {
		rgb += w;
        return 1;                     // Colour modified to fit RGB gamut
    }

    return 0;                         /* Colour within RGB gamut */
}



double computeScatCrossSection( double radius, double lambda) {
	
    interval radius_interval, spectrum_interval;
	double scatCrossSection;

	
	// check if key falls within the min and max values of the table. Otherwise, return false;
	/*
    if ( ! ( radius  >=  radii[ 0] && radius <= radii[ nRadii-1] )  ) {
		return (-1.0);
		
	}
	
	if ( ! ( lambda  >=  spectralSamples[ 0] && lambda <= spectralSamples[ nSpectralSamples-1] )  ) {
		return (-1.0);
		
	}
	*/
	
	// get the index to the texture for (radius, lambda)

	vec2 index = vec2( (radius  - radiusStart) / (radiusEnd - radiusStart), 
	               (lambda - lambdaStart) / ( lambdaEnd - lambdaStart) );
	
	float scatCrossSection = texture3D( uScatTex, index);
	
    return scatCrossSection;
}


double computePhaseFunction( double radius, double lambda, double psi) {
	
  
    double  P = 0.0;
	


	// check if key falls within the min and max values of the table. Otherwise, return false;
	/*
    if ( ! ( radius  >=  radii[ 0] && radius <= radii[ nRadii-1] )  ){
		return ( -1.0);
		
	}  
	
	if ( ! ( lambda  >=  spectralSamples[ 0] && lambda <= spectralSamples[ nSpectralSamples-1] )  ) {
		return( -1.0);
		
	}
	if ( ! ( psi  >=  thetas[ 0] && psi <= thetas[ nThetas-1] )  ) {
		return(-1.0);
		
	}
	*/

	
	// get the index to the texture for (radius, lambda)

	vec3 index = vec3 ( (radius  - radiusStart) / (radiusEnd - radiusStart), 
	               (lambda - lambdaStart) / ( lambdaEnd - lambdaStart),
				   (theta - thetaStart ) / ( thetaEnd - thetaStart) ) ;
	
	float P = texture3D( uPhaseTex, index);

    return P;
}



vec3 computeSurfaceColorForCurrentFrag(vec3 vPosInEyeSpace, bool bHit, 
                                       float tmin, float tmax, float phi_obs, float phi_sun ) {

  if ( !bHit ) { // the AABB box is not hit => just return the black color
    return  vec3(0,0,0);
	
  }

  
// attenuate the color throughout the volume
      
  double sunTau; // optical depth
				 
  double rayPathNormal =  (tmax - tmin) * cos (phi_obs);
	
  double sunRayPath = rayPathNormal / cos (phi_sun);	 
			
  sunTau = uDropDensity * scatCrossSection * sunRayPath; // attenuate from the sphere to the hit
		                                              // point of the volume
    	  
 // attenuate the light color through volume: attenuation is the same for each wavelength
 
  vec3 attenuatedLightColor = uLightRGBColor * exp ( - sunTau ); 

  vec3 vViewDir = normalize(-vPosInEyeSpace);

  vec3 normal = normalize(vNormal);
  float nDotL = dot( normal, uSunRayDir); // use the directional light direction

  vec3 reflectionDir = normalize( 2.0 * normal * nDotL - uSunRayDir);
  float rDotV = max(0.0, dot(reflectionDir, vViewDir));
  float specularCoeff = pow(rDotV, 32.0);
  float diffuseCoeff = max(nDotL, 0.0);

   
  vec3  diffuseColor  = uTextureFlag * attenuatedLightColor * texture2D( uTexUnit0, vTexCoord0) * diffuseCoeff + 
               (1.0 - uTextureFlag) * attenuatedLightColor * uMaterialColor  * diffuseCoeff;

  vec3 whiteLight = vec3( 0.6, 0.6, 0.6);

  vec3  specularColor = whiteLight * specularEoeff;
     
  vec3  ambientColor = vec3(0.1, 0.1, 0.1);



  vec3 surfaceColor = ambientColor + diffuseColor + specularColor;
  
  
  // Attenuate the combined Color through the water drop volume from tmin to tmax
   double scatCrossSection = computeScatCrossSection( uRadius,  lambdaStart);

   	
  // attenuate the backgroundSurfaceColor  throughout the volume
                  
   double rayTau = uDropDensity * scatCrossSection  * ( tmax -tmin );
		   

   surfaceColor = surfaceColor  * exp ( - rayTau );   

   return surfaceColor;
  }



vec3 computeRainbowColorForCurrentFrag( in vec3 rayDir,
                            in vec3 uSunRay, bool bHit, in float tmin, in float tmax, in float phi_obs, in float phi_sun) {

  
    if ( !bHit ) { // the AABB  is not hit => just return the black color
      return  vec3(0,0,0);
	  
    }

  
	double scatCrossSection = computeScatCrossSection( uRadius,  lambdaStart);
	
	vec3  volumeColor;
	
   
	// attenuate the color throughout the volume
    		 
	double rayPathNormal =  (tmax - tmin) * cos (phi_obs);
	double tau_N = uDropDensity * scatCrossSection  * rayPathNormal;		
    			  
	  // compute rainbow volume color
	  // check if the psi angle is less than 130, in which case the rainbow is not computed

     double psi = acos ( dot ( -rayDir, -uSunRayDir ) ) * 180 / pi;

	 if ( psi <  130 || psi > 141.8 ) { // no volume color for rainbow:  The red emerges at a smaller angle psi than the violet
				      //scattering is concentrated at theta_0 = 137.97 but also occurs for angles greater than that. 
				      // This is what causes the white appearance `insidea the primary bow (and `outsidea the secondary bow).
				  //Note the gap between 129 and 138. This is a region of negligible scattering (from higher orders) and appears 
				  //as a dark space between the primary and secondary rainbows

				volumeColor = vec3(0,0,0);
			  }
	 else {
     	        volumeColor = calculate_rainbowColor ( uRadius, uDropDensity, phi_sun, phi_obs, psi, tau_N);  
				
		} 

     return volumeColor;

}

double calculate_irradiance_of_sun (double lambda) {
    
    double R_Earth = 1.496e+8;
    double R_Sun = 6.95e+5;
    double h = 6.6261e-34;      // Planck's constant
    double c = 2.9979e+8;       // speed of light in vacuo
    double k = 1.3806e-23;      // Boltzmann's constant
    double T = 5782;            // Sun's temperature in Kelvin
    double Isun;
    
    Isun = pow((R_Sun/R_Earth),2) * 2*PI*h*pow(c,2)  / ( pow(lambda, 5) * ( exp( h*c/(lambda*k*T) ) -1 ) ) * 1.0e-9  ;
	   //Moon:  1.0e-9 is multiplied to convert the irradiance per meter to the irradiance per nanometer
	   // Original had 1.0e-10. 
    
    return Isun;
}

vec3 calculate_rainbowColor (double radius, double uDropDensity, double phi_sun, double phi_obs,
                             double psi, double tau_N)
{
    double Isun;
    
    double lambda_m;

    double irainbow[ nSpectralSamples ];
     
	vec3 rgbColor = vec3(0.0, 0.0, 0.0);
	vec3 xyzColor = vec3(0.0, 0.0, 0.0);

	double X = 0, Y = 0, Z = 0, XYZ;

   	double lambda;
	int j;

    for (lambda = lambdaStart, j = 0; lambda < lambdaEnd; lambda += lambdaStep*2, j++) {
		
		lambda_m = lambda * 1.0e-9;				// convert lambda to meters from nanometers

        Isun = calculate_irradiance_of_sun(lambda_m); // isun: irradiance per nanometer

		// ignore the attenuation of the sun light while it passes through the drop volume

		double Lsun = Isun /  pi; // convert the irradiance to radiance 

        double particlePhase = computePhaseFunction( radius,  lambda, psi);

		double avgPhasePerVolume  =  uDropDensity * particlePhase;
		        
        irainbow[j] = ( ( Lsun* avgPhasePerVolume  / cos(phi_obs) ) / ( 1/cos(phi_sun) + 1/cos(phi_obs) ) )
		                * ( 1 - exp(-tau_N * ( 1/cos(phi_sun) + 1/cos(phi_obs) ) ) );

	}

	// convert spectral color to rgb color

	double xBar, yBar, zBar;

	/* cie_colour_match[(lambda - 380) / 5][0] = xBar
      cie_colour_match[(lambda - 380) / 5][1] = yBar
      cie_colour_match[(lambda - 380) / 5][2] = zBar
	  ==> cie_coulor_match [] from 380 nanometer to 780 nanometer
    */
	
	// 400 ~ 700 into 120 steps with lamdbaStep - 2.5nm
	for (lambda = lambdaStart, j =0; lambda < lambdaEnd; lambda += lambdaStep*2, j++) {
		
		//intensity = bb_spectrum(lambda);  // You already have the function that computes the irradiance of the sun

		xBar = cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep*2 ) + 0]; // it is assumed that lambdaStep is 5 nanometer
		yBar = cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep*2 ) +1 ];
		zBar = cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep*2 ) +2];

		X += irainbow[j] * xBar * lambdaStep*2; 
		Y += irainbow[j] * yBar * lambdaStep*2;
		Z += irainbow[j] * zBar * lambdaStep*2;
		
		
    }
	
		
	XYZ = (X+Y+Z);
	xyzColor = vec3( X/XYZ, Y/XYZ, Z/XYZ );
	rgbColor = xyz_to_rgb(g_cs, xyzColor);

	
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
  
  vec3 rayPos = vPosInEyeSpace;

    // current position along the ray
  vec3 rayDir = normalize( rayPos - uEyePos );
    // ray direction
	
  Ray ray = makeRay( uEyePos, rayDir);

  vec3 normalAtFront;
  
  vec3 diffuse, specular;

  float phi_obs, phi_sun;

  float tmin, tmax;

  vec3[2] aabb;

  // Get the camera space positions of the AABB box corners.

  aabb[0] = uInvEyeMatrix * uModelMatrixAABB * uAABB[0];
  aabb[1] = uInvEyeMatrix * uModelMatrixAABB * uAABB[1];

 

  bool bHit = intersection_distance_no_if( ray, aabb, tmin, tmax, normalAtFront );	

  // get the angle between the sun ray and the front face of the AABB
  // get the angle between the observer ray and the front face of the AABB
    
  phi_obs = acos( dot( normalAtFront, -rayDir) );       // in radians
  phi_sun = acos( dot( normalAtFront,  uSunRayDir) );          // in radians

  // compute the surface color at vPosInEyeSpace and attenuated the color through
  // volume space between tmin and tmax.

  vec3 surfaceColor = computeSurfaceColorForCurrentFrag( vPosInEyeSpace, bHit, tmin, tmax, phi_obs, phi_sun );
  
  vec3 volumeColor = computeRainbowColorForCurrentFrag(dirDir, uSunRayDir, bHit, tmin,  tmax, phi_obs, phi_sun );

  fragColor = vec4(surfaceColor, 1) + vec4(volumeColor,1);
}

/* 

builtin functions:

radians(degrees)
degrees(radians)
sin(angle in radians)
asin, acos, pow(x,y), exp, log, sqrt, abs, sign, floor, ceil, mod, fract, mix(x,y,a)
length(vector x), distance(p0,p1), dot, cross, normalize, lessThan(vec,vec) equal(vec,vec)
*/

