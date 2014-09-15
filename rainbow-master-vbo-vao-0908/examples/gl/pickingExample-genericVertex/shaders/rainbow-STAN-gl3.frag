
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

//uniform sampler2D uColorTex;
//uniform sampler2D uDepthTex;

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
in   mat3 vNTMat;

out vec4 fragColor;
out vec4 spect1;

// parameters for rainbow computation

float pi = acos( -1.0);  // -1.0 const float in default; -1.0f = float

const int nSpectralSampleSteps = 120;
const int nThetaSteps = 100; // 100 steps between theta 130 deg and theta 142;

const int nRadiusSteps  = 2;

const int nSpectralSamples = nSpectralSampleSteps +1;
const int nSamples = nSpectralSamples;


const int nRadii = nRadiusSteps +1;
const int nThetas = nThetaSteps +1; // 100 steps between theta 130 deg and theta 142;

//uniform float uPhaseFunction[nRadii * nSpectralSamples * nThetas];

//float irainbow[ nSpectralSamples ];
float particlePhase[nSpectralSamples ];

// uniform block for primary spectrums

const float CIE_Y_integral = 106.856895;
const float  lambdaStart = 400;   // 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
const float lambdaEnd = 700;
     
const  float radiusStart =  1.0e-3;  // 
const  float radiusEnd = 2.0e-3;           // 2 mm = 2 * e-3 m

const	float thetaStart = 130.0;
const	float thetaEnd = 142.0;

const float lambdaStep = ( lambdaEnd - lambdaStart) / float( nSpectralSamples -1);
const float thetaStep = ( thetaEnd - thetaStart) / float (nThetas -1);
const float radiusStep = (radiusEnd - radiusStart) / float (nRadii -1);



struct Ray {
  vec3 origin;
  vec3 direction;
  vec3 inv_direction;
  int sign[3];
};


float[nSpectralSamples] X;
float[nSpectralSamples]  Y;
float [nSpectralSamples]  Z;

float[nSpectralSamples] rgbRefl2SpectWhite;
float[nSpectralSamples] rgbRefl2SpectCyan;
float[nSpectralSamples] rgbRefl2SpectMagenta;
float[nSpectralSamples] rgbRefl2SpectYellow;
float[nSpectralSamples] rgbRefl2SpectRed;
float[nSpectralSamples] rgbRefl2SpectGreen;
float[nSpectralSamples] rgbRefl2SpectBlue;


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



float[nSpectralSamples] spectralColor(vec3 rgb) {
// primarySpectrums
// As background -- GLSL looks a lot like C, but compiles a bit different. 
//Things are very unrolled, and conditionals may be executed in parallel and switched at the end,
// that sort of thing. Depends on the hardware...

//You can use **loop indices or constants** [== contant-index-expression] to index into arrays. 
//The assignment in your loop is ok, but the access by tileID isn't.

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

    return r; 

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

          return r; 
		  
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

            return r; 
			
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

			 // return rgbRefl2SpectWhite;
		    for (int j =0; j < nSpectralSamples; j++) {
               r[j] *= .94;
            }

            return r; 
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

          return r; 
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

           return r; 

  } // if case (3.2)
 


} // spectralColor

vec3 XYZToRGB(vec3 xyz);

vec3 toXYZ(float[nSamples] spect);
 
vec3 toRGB(float[nSamples] spect)  {
        vec3  xyz;
        xyz = toXYZ(spect);
        return XYZToRGB(xyz);
    }

vec3  toXYZ(float[nSamples] spect)  {
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

vec3 XYZToRGB(vec3 xyz) {
    vec3 rgb;
    rgb[0] =  3.240479f*xyz[0] - 1.537150f*xyz[1] - 0.498535f*xyz[2];
    rgb[1] = -0.969256f*xyz[0] + 1.875991f*xyz[1] + 0.041556f*xyz[2];
    rgb[2] =  0.055648f*xyz[0] - 0.204043f*xyz[1] + 1.057311f*xyz[2];

	return rgb;
}//XYZToRGB


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


float avgPhaseFunction(float radius, float lambda, float psi) {

    float particlePhase  = computePhaseFunction( radius,  lambda, psi); 

	//float particlePhase = 1.230634e-002;

    float avgPhasePerVolume  =  uDropDensity * particlePhase;
    return avgPhasePerVolume;

}

float calculate_radiance_of_sun(float lambda) {
    
    float R_Earth = 1.496e+8;
    float R_Sun = 6.95e+5;
    float h = 6.6261e-34;      // Planck's constant
    float c = 2.9979e+8;       // speed of light in vacuo
    float k = 1.3806e-23;      // Boltzmann's constant
    float T = 5782;            // Sun's temperature in Kelvin
    float Isun;

			
		float lambda_m = lambda * 1.0e-9;				// convert lambda to meters from nanometers

        Isun = pow( (R_Sun/R_Earth),2.0 )  * 2*pi*h*pow(c,2.0)  / ( pow(lambda_m, 5.0) * ( exp( h*c/ (lambda_m*k*T) ) -1 ) ) * 1.0e-9  ;
	      //Moon:  1.0e-9 is multiplied to convert the irradiance per meter to the irradiance per nanometer
	      // Original had 1.0e-10.
		float Lsun = Isun / pi; //convert the irradiance to radiance 

		return Lsun;
		
} //calculate_radiance_of_sun()


float[nSpectralSamples] calculate_radiance_of_sun() {
    
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
		
		Lsun = calculate_radiance_of_sun( lambda );
		LSpectrumSun[ j ] = Lsun;		

		lambda += lambdaStep;
		
	}  // for 
    
    return LSpectrumSun;

} //calculate_radiance_of_sun()



float[nSpectralSamples] spectralColor(vec3  rgb );
vec3 toRGB(float[nSpectralSamples] spect);

float[nSpectralSamples]  LSpectrumIn( float[nSpectralSamples] LSpectrumSun, 
                                      vec3 currPos, vec3 dirToLight, float sigma_e) {

  // attenuate the sun light Lsun0 from the entering point to AABB to Pos

  float tmin, tmax;
  float[nSpectralSamples] LSpectIn;

  vec3 lightDirInBoxFrame = vec3( uInvModelMatrixAABB * uEyeMatrix * vec4( dirToLight, 0) );
  vec3 lightOriginInBoxFrame = vec3( uInvModelMatrixAABB * uEyeMatrix * vec4( currPos,1) );

  tmin = 0.0; tmax = 0.0;

  bool isHit =   rayOBBIntersect ( lightOriginInBoxFrame, lightDirInBoxFrame, 
                                  uAABBmin, uAABBmax, tmin, tmax);

 

  if ( isHit ) {
       // spect1 = vec4(tmin, tmax, 0.0, 1.0);
 // return LSpectrumSun;
	 // the light origin is at currPos [with t_i = 0, within the water view volume] goes to tmax 
	 // at which it exits AABB; tmin < 0 is the opposite direction. 

	 //   L(t_i, w_fromLight(t_i) ) = LSpectrumIn is given by
      //   tau(tmax_sun, t_i) Lsun( tmax_sun, dirFromLight(tmax_sun) ) 
      //   = exp( -sigma_e ( tmax_sun - t_i) )  Lsun( tmax_sun, dirFromLight(tmax_sun) ) 

	  for (int j=0; j < nSpectralSamples; j++ ) {
	       LSpectIn[j] = exp( - sigma_e * tmax ) * LSpectrumSun[j];
	  }
	
	 // return LSpectIn; for debugging

	 return LSpectrumSun;

  }
  else {
   // spect1 = vec4(-1,-1, 0.0, 0.0);
      return LSpectrumSun;
  }
  
} // LSpectrumIn()


float[nSpectralSamples] singleScatteringAndAttenuation(float radius, vec3 surfaceColor, float zEye, 
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

  bool isHit =  rayOBBIntersect ( uEyeOriginInBoxFrame, rayDirInBoxFrame, 
                                  uAABBmin, uAABBmax, tmin, tmax);
  
  
   //Loop iterations that only some fragments execute are non-uniform control flow.
   // In that case, if the accessed texture uses **mipmapping or anisotropic** filtering of any kind, 
   // then any texture function that is not "Lod" or "Grad" will retrieve undefined results

  if ( !isHit ) { // the AABB  is not hit? 	 getBackgroundColor:
  	 
	 // for debugging
  
      //spect1 = vec4(0,0,0,-1);
         	 
	  LSpectrumOut0 = spectralColor( surfaceColor);
	  return LSpectrumOut0; // this surface color is NOT attenuated through the water volume
	                       // and reaches the eye without attenuation
	              
  } // if (!isHit)

  else { //  hit the front at tmin and the back face of AABB box at tmax;
  

	 // if the background scene lies within the AABB box, the attenuation of the background color
	 //  begins at the intersection point between the ray and the background scene. The 
	 //  attenuation of the volume color begins also at the intersection point. 
	 //  Otherwise, the attenuation of background and volume color begins at the back face of the box.

	  LSpectrumOut0 = spectralColor( surfaceColor );
	  	  
  // compute the intensity of the light with wavelenth lambda scattered and attenuated 
  // in the direction of -rayDir
  
	  float tmaxAtten = tmax; // the t value at which attenuation begins

	  //float distFromEye = -zEye;

	 // if ( distFromEye >= tmin && distFromEye <= tmax ) {
	 //    tmaxAtten = distFromEye;
	 // }
	 // else {
	 //  tmaxAtten = tmax; 
	 // }
	  	    
     int N = 30;
     float deltaT  =  (tmax - tmin) / N;
   

  // get the surface reflection color which has been attenuated through the volume
  // assume that the surface color attenuated only through water volume and not through atmosphere

     // LSpectrumOut = exp( - sigma_e * ( tmaxAtten - tmin) )  * LSpectrumOut0;

     for (int j=0; j < nSpectralSamples; j++ ) {
	    LSpectrumOut[j] = exp( - sigma_e * ( tmaxAtten - tmin) )  * LSpectrumOut0[j];
	 }
 

  // get the light which has been scattered in the eye direction by each region in the volume 
  // starting from tmaxAtten.

     float[nSpectralSamples] LSpectIn;

     if ( isPointLight) { // point light
       vec3 lightPos = lightPosOrDir;

       for (float t = tmaxAtten; t >= tmin; t -= deltaT  ) { 
      
	    vec3 currPos =  vec3(0,0,0) +  rayDir * t;
		vec3 dirFromPointLight =  currPos - lightPos;

	    dirFromPointLight = normalize( dirFromPointLight );

		LSpectIn = LSpectrumIn( LSpectrumSun, currPos, -dirFromPointLight, sigma_e);
		
		// for debugging
		//return LSpectIn;


		// The light with given lambda is scattered in the direction of -rayDir with varying amount
		// depending on the scattering angle between the light direction and -rayDir, which varies
		// with each t. Consider only the scattering angles  within the range [ thetaStart, thetaEnd]. 
		//  The other scattering angles constributs less to the scattered intensity of lambda
		//   They are ignored for computational reasons. Some scattering spectrum is computed for
		// every ray direction (-rayDir)?? 
	  
	    float psi = acos ( dot ( -rayDir, dirFromPointLight) ) * 180.0 / pi;
		
		if (  psi >=  thetaStart  &&  psi <= thetaEnd    ) { // only the light scattered within these angles
		                                                     // contribute to the rainbow color.
															 // Otherwise, the background scene color is used 
			float lambda = lambdaStart;
			int   lambdaIndex = 0;	
		    for ( int j = 0; j < nSpectralSamples;  j++) {
		
		       LSpectrumOut[j] += exp( - sigma_e * ( t - tmin) )
		                      * sigma_s * avgPhaseFunction(radius, lambda, psi)
	                        *  LSpectIn[j];
               lambda += lambdaStep;
			   
			   
	        } // for
		   	   
        }   //if 
       
       } // for 

       return LSpectrumOut; 

     } // if (point light)
	
     else { // directional light

      vec3 dirFromSun  = lightPosOrDir;
	  float psi = acos ( dot ( -rayDir, dirFromSun ) ) * 180.0 / pi;

	  // for debugging
	 // spect1 = vec4(dirFromSun, 1.0);

	  //spect1 = vec4(dot ( -rayDir, dirFromSun ), acos(dot ( -rayDir, dirFromSun ) ), psi, 1.0 );
	  //return LSpectrumOut;

	 // consider only the cases where the light is scattered in the direction of the eye
	 // with scattering angle within [ thetaStart, thetaEnd]. These scattering angles contricute
	 // to the rainbow and its neighborhood. The cases with other scattering
	 // angles contribute to the further region outside of the rainbow and its neighbor.
	 // These regions are not considered and the background color are assumed to dominate. 


	  // SUM_{i=0, i=N-1} exp( - sigma_e ( t_i - tmin) ) sigma_s  f_p( dot(-rayDir, dirFromLight(t_i) ) ) L(t_i, dirFromLight(t_i) )
      //   L(t_i, dirFromLight(t_i) ) = LSpectrumIn is given by
      //   tau(tmin_sun, t_i) Lsun( tmin_sun, dirFromLight(tmin_sun) ) 
      //   = exp( -sigma_e ( t_i - tmin_sun) )  Lsun( tmin_sun, dirFromLight(tmin_sun) ) 
	  //  where tmin_sun is the point at which the sun light enters the volume to reach t_i.
	  
	  if (  !( psi >=  thetaStart  &&  psi <= thetaEnd )    ) {
	      return LSpectrumOut; // just return the background color
	  }  

	  else {	// psi is a rainbow scattering angle
		// LSpectrumOut = LSpectrumOut + SUM_{t = tmin, tmaxAtten} exp( - sigma_e * ( t - tmin) )  * sigma_s 
		//                                          * avgPhaseFunction( dot ( -rayDir, dirFromSun(t) ) ) 
		//                                          * LSpectIn(t, dirFromLight(t) )

		// for debugging

		//return LSpectrumOut;
		float[nSpectralSamples] LSpectIn;

        for (float t = tmaxAtten; t >= tmin; t-= deltaT  ) { 
      
	          vec3 currPos = vec3(0,0,0) +  rayDir * t;
		 		     
			  LSpectIn = LSpectrumIn(LSpectrumSun, currPos, -dirFromSun, sigma_e);
		  
		// for debugging
		//return LSpectIn;

		      float lambda = lambdaStart;
		    		   		   
		      for ( int j =0; j < nSpectralSamples; j++) {			  
				 
                 LSpectrumOut[j] += exp( - sigma_e * ( t - tmin) ) * sigma_s  
		                           * avgPhaseFunction(radius, lambda, psi) * LSpectIn[j];
                 lambda += lambdaStep;

				  
		     } // for each lambda	
			 				   
        } // for t		
			
		return 	  LSpectrumOut; 

      } // else ( psi is rainbow angle)

	    
    } // else (directional light) 	  


  } // else (the ray hits the volume )

} // singleScatteringAndAttenuation


float[nSpectralSamples] singleScatteringAndAttenuation_stan(float radius, vec3 surfaceColor, float zEye, 
                              vec3  dirFromSun,  vec3 rayDir, float sigma_s, float sigma_e )  {

  // check if the current ray through the fullscreen quad intersects the water volume AABB.
  // otherwise, the background color is used as the color at the pixel color of the ray.

  float tmin, tmax;
  float[nSpectralSamples] LSpectrumOut, LSpectrumOut0;
 

  vec3 rayDirInBoxFrame = vec3( uInvModelMatrixAABB * uEyeMatrix * vec4(rayDir, 0) );

  bool isHit =  rayOBBIntersect ( uEyeOriginInBoxFrame, rayDirInBoxFrame, 
                                  uAABBmin, uAABBmax, tmin, tmax);
  
  if ( !isHit ) { // the AABB  is not hit? 	 getBackgroundColor:
  	 
      // for debugging
	   
	  surfaceColor = vec3(0,0,0);
	     	 
	  LSpectrumOut0 = spectralColor( surfaceColor);
	  return LSpectrumOut0; // this surface color is NOT attenuated through the water volume
	                       // and reaches the eye without attenuation
	              
  } // if (!isHit)

  else { //  hit the front at tmin and the back face of AABB box at tmax;
  

	 // if the background scene lies within the AABB box, the attenuation of the background color
	 //  begins at the intersection point between the ray and the background scene. The 
	 //  attenuation of the volume color begins also at the intersection point. 
	 //  Otherwise, the attenuation of background and volume color begins at the back face of the box.

	  // for debugging
	  surfaceColor = vec3(0,0,0);

	  LSpectrumOut0 = spectralColor( surfaceColor );
	  
	  // attenuate the surface color	
     for (int j=0; j < nSpectralSamples; j++ ) {
	    LSpectrumOut[j] = exp( - sigma_e * ( tmax - tmin) )  * LSpectrumOut0[j];
	 }
   
   
  // compute and attenuate the intensity of the light scattered by the water volume
  // in the direction of -rayDir

     float scatCrossSection = computeScatCrossSection( radius,  lambdaStart);

     float psi = acos ( dot ( -rayDir, dirFromSun ) ) * 180.0 / pi;
   
     float phi_obs = acos( dot( uFrontNormalToAABB, -rayDir) );       //  rayDir in eyeSpace
     float phi_sun = acos( dot( uFrontNormalToAABB,  uSunRayDir )  );

	  // compute the optical depth normal to the AABB

	  float rayPathNormal =  (tmax - tmin) * cos (phi_obs);
	  float tau_N = uDropDensity * scatCrossSection  * rayPathNormal;	
	  	
    // compute the surface color at vPosition and attenuated the color through
    // volume space between tmin and tmax.

	  tau_N = 1.0;



	  float lambda = lambdaStart;

      for (int j = 0; j < nSpectralSamples; lambda += lambdaStep, j++) {
		
		float Lsun = calculate_radiance_of_sun(lambda); // Lsun: radiance per nanometer

		// ignore the attenuation of the sun light while it passes through the drop volume
				
        float particlePhase= computePhaseFunction( radius,  lambda, psi); 
		
		//float particlePhase = 1.230634e-002;

		float avgPhasePerVolume  =  uDropDensity * particlePhase;
		        
        LSpectrumOut[j] = ( ( Lsun* avgPhasePerVolume  / cos(phi_obs) ) / ( 1/cos(phi_sun) + 1/cos(phi_obs) ) )
		                * ( 1 - exp(-tau_N * ( 1/cos(phi_sun) + 1/cos(phi_obs) ) ) );
         
        
	  } // for j

      return LSpectrumOut;

  } // else (the ray hits the volume )

} // singleScatteringAndAttenuation_stan


vec3 calculate_rainbowColor (float radius, vec3 surfaceColor, float zEye, bool isPointLight,  vec3 lightPosOrDir, vec3 rayDir, 
                             float sigma_s, float sigma_e) {
	// compute the scattered light, which is attenuated along the direction 
	// -rayDir 
	
	float[nSpectralSamples]	irainbow;
		
	float[nSpectralSamples] LSpectrumSun = calculate_radiance_of_sun(); 



	// LSpectrumSun: radiances for wavelengths [nanometer]

	irainbow  = singleScatteringAndAttenuation(radius, surfaceColor, zEye, 
	                          isPointLight, lightPosOrDir, LSpectrumSun,
	                          rayDir, sigma_s, sigma_e );
    
	return toRGB(irainbow);
	  
}//calculate_rainbowColor

vec3 calculate_rainbowColor_stan (float radius, vec3 surfaceColor, float zEye,  vec3 dirFromSun, vec3 rayDir, 
                             float sigma_s, float sigma_e) {
	// compute the scattered light, which is attenuated along the direction 
	// -rayDir 
	
	float[nSpectralSamples]	irainbow;
		
	float[nSpectralSamples] LSpectrumSun = calculate_radiance_of_sun(); 



	// LSpectrumSun: radiances for wavelengths [nanometer]

	irainbow  = singleScatteringAndAttenuation_stan(radius, surfaceColor, zEye, 
	                          dirFromSun,  rayDir, sigma_s, sigma_e );
    
	return toRGB(irainbow);
	  
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



 // for debugging: dump the content of the FBO color texture 
  // get the pixel coordinates of the current fragment

	  //ivec2 pixelPos = ivec2( gl_FragCoord.xy);
	      
	 // fragColor =vec4(gl_FragCoord.x, gl_FragCoord.y, 0, 1);
	
	  // get the surface color and depth of the current fragment
	 
	  // vec3  surfaceColor = vec3( texelFetch(uColorTex, pixelPos, 0) ); // color

	   // for debugging
	   //fragColor = vec4(surfaceColor, 1.0);
	   //return;

	   // for debugging
	   //fragColor = vec4( toRGB( spectralColor(surfaceColor) ), 1.0 );
	   //return;


      //vec3  surfaceColor = vec3( texture(uColorTex, vTexCoord) ); // color
	  //float zWin  = texelFetch(uDepthTex, pixelPos, 0).r;  // zWin in non linear range [0,1]


      // conversion range ([0,1] into NDC range  [-1,1]
      //float zNDC = zWin * 2.0 - 1.0; // ( zWin = 1/2 zNDC + 1/2 by Viewport transformation: zWin in [0,1] )

	  //float zEye = -uProjectionMatrix[2][3]/ ( zNDC + uProjectionMatrix[2][2] ); 

      //if ( uProjectionMatrix[2][3] > 0.0 && uProjectionMatrix[2][2] > 0.0 ) { // good
	  //    fragColor = vec4(0.0, 1.0, 1.0, 1.0);
	  //}
      //else // bad
	  // fragColor = vec4(0.0, 1.0, 1.0, 1.0);

	
	   // for debugging
	  //spect1 = vec4( zWin, zNDC, zEye, 1.0);
	
	 // return;

	  // In our system, the third row of the projection matrix is different from a typical convention. 
	  // See p. 107 of Steven Gortler's book. This is a source of confusion. Be careful when you consult internet
	  // material which assumes the typical convention.

	  	   
	   // In this shader, depth test is enabled. so the incoming fragment is discarded if
	   // the depth comparison fails. The depth 
		//NDC (after projection), the depth from the eye is [-1,1] where -1 is the farthest, and 1 is
		// the film plane. But in the fragment shader, the viewport transformation
		// changes [-1,1] to [0,1], where 0 is the farthest, and 1 is the film plane.
	  // This shader renders the near plane of the view frustum, so that the rendered image covers the 
	  // full screen square. 
	  //
	  // ProjMat = [ 1/ a * tan(th/2) 0           0                0 
	  //             0                1/tan(th/2) 0                0
	  //             0                0          f+n/f-n        -2fn/f-n
	  //             0                0           -1               0 ]
	  //  zClip = [ ProjMat * zEye ]. z. 
	  // http://www.txutxi.com/?p=182
	  //http://www.opengl.org/wiki/Compute_eye_space_from_window_space

 
  vec3 surfaceColor = vec3(0.0, 0.0, 0.0);
  float zEye = 0.0;

  vec3 volumeColor;
  
   
 float scatCrossSection = computeScatCrossSection( uRadius,  lambdaStart);  // scatCrossSection is independent of
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
	     volumeColor = calculate_rainbowColor (uRadius, surfaceColor, zEye, isPointLight, lightPos,  rayDir, sigma_s, sigma_e); 
			  
	  //lightPos = uLight2Pos;
	  //volumeColor += calculate_rainbowColor (uRadius, surfaceColor, zEye, isPointLight, lightPos, rayDir,sigma_s, sigma_e); 
	
    }
	else {

	 // uSunRayDir is  already conveted into eye space in the main program.
	      vec3 dirFromSun  = -uSunRayDir;
	     volumeColor = calculate_rainbowColor_stan (uRadius, surfaceColor, zEye,  dirFromSun, rayDir, sigma_s, sigma_e);  
		 // volumeColor = calculate_rainbowColor (uRadius, surfaceColor, zEye,  isPointLight, dirFromSun, rayDir, sigma_s, sigma_e);  

	 	  
		  
	}

	
	fragColor  = vec4( volumeColor, 1);	
       
	// no volume color for rainbow:  The red emerges at a smaller angle psi than the violet
   //scattering is concentrated at theta_0 = 137.97 but also occurs for angles greater than that. 
   // This is what causes the white appearance `insidea the primary bow (and `outsidea the secondary bow).
   //Note the gap between 129 and 138. This is a region of negligible scattering (from higher orders) and appears 
   //as a dark space between the primary and secondary rainbows
   
	   
} // main




