
#version 130

uniform int uWindowWidth; 
uniform int uWindowHeight; 

// ***Opengl 3.2, 2009 => GLSL 1.5 ( //): Compatibility Profile: Openframework uses this version

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


in vec2 vTexCoord;
in vec3 vPosition;

out vec4 fragColor;
out vec4 phaseSpectrum;
out vec4 rainbowColorSpectrum;

// parameters for rainbow computation

float pi = acos( -1.0);  // -1.0 const float in default; -1.0f = float

const int nSpectralSampleSteps = 120;
const int nThetaSteps = 100; // 100 steps between theta 130 deg and theta 142;

const int nRadiusSteps  = 2;

const int nSpectralSamples = nSpectralSampleSteps +1;
const int nSamples = nSpectralSamples;


const int nRadii = nRadiusSteps +1;
const int nThetas = nThetaSteps +1; // 100 steps between theta 130 deg and theta 142;


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

float XSpectrum(int j ) {
  return  texelFetch(uXYZSpectrums, 0 * nSpectralSamples + j, 0).r; 
}

float YSpectrum(int j ) {
  return texelFetch(uXYZSpectrums, 1 * nSpectralSamples + j, 0).r; 
}

float ZSpectrum(int j) {
  return texelFetch(uXYZSpectrums, 2 * nSpectralSamples + j, 0).r; 
}

float  whiteSpectrum(int j ) {
  return texelFetch(uPrimarySpectrums, 0 * nSpectralSamples + j, 0).r; 
}

float cyanSpectrum(int j ) {
  return texelFetch(uPrimarySpectrums, 1 * nSpectralSamples + j, 0).r; 
}
float magentaSpectrum(int j ) {
 return texelFetch(uPrimarySpectrums, 2 * nSpectralSamples + j, 0).r; 
}
float  yellowSpectrum(int j ) {
  return  texelFetch(uPrimarySpectrums, 3 * nSpectralSamples + j, 0).r; 
}

float  redSpectrum(int j ) {
  return texelFetch(uPrimarySpectrums, 4 * nSpectralSamples + j, 0).r; 
}
float greenSpectrum(int j ) {
 return texelFetch(uPrimarySpectrums, 5 * nSpectralSamples + j, 0).r; 
}

float blueSpectrum(int j ) {
  return  texelFetch(uPrimarySpectrums, 6 * nSpectralSamples + j, 0).r; 
}

float[nSamples] spectralColorOptimal (vec3 rgb) {

//You can use **loop indices or constants** [== contant-index-expression] to index into arrays. 

float[nSamples]  r;
for (int j=0; j < nSamples; j++ ) {
     r[j] = 0.0;
}

for (int j =0; j < nSpectralSamples; j++ ) {

 float white = whiteSpectrum(j);
 float cyan = cyanSpectrum(j);
 float magenta= magentaSpectrum(j);
 float yellow= yellowSpectrum(j);
 float red = redSpectrum(j);
 float green = greenSpectrum(j);
 float blue  = blueSpectrum(j);


  // Convert reflectance spectrum to RGB
  
  if (rgb[0] <= rgb[1] && rgb[0] <= rgb[2]  && rgb[1] <= rgb[2] ) {
                          // rgb[0] <= rgb[1] << rgb[2], case (1.1)	  

   // Compute reflectance _SampledSpectrum_ with _rgb[0]_ as minimum

     r[j] += rgb[0] * white;
	 r[j] += (rgb[1] - rgb[0]) * cyan;		
	 r[j] += (rgb[2] - rgb[1]) * blue;
	
     r[j] *= .94;
    
  } // if case (1.1)
  
  
  if ( rgb[0] <= rgb[1] && rgb[0] <= rgb[2]  && rgb[2] <= rgb[1]   ) { 
                      // rgb[0] << rgb[2] << rgb[1]: case (1.2)

   // Compute reflectance _SampledSpectrum_ with _rgb[0]_ as minimum

    r[j] += rgb[0] * white;
	r[j] += (rgb[2] - rgb[0]) * cyan;
	r[j] += (rgb[1] - rgb[2]) * green;
			
    r[j] *= .94;
   		  
  } //  case (1.2)

  
  if (rgb[1] <= rgb[0] && rgb[1] <= rgb[2]  &&  rgb[0] <= rgb[2] ) { // case (2.1)
                      // rgb[1] <= rgb[0] <= rgb[2]: case (2.1)
	   
    // Compute reflectance _SampledSpectrum_ with _rgb[1]_ as minimum
	 r[j] += rgb[1] * white;
	 r[j] += (rgb[0] - rgb[1]) * magenta;
	 r[j] += (rgb[2] - rgb[0]) * blue;
		
	 r[j] *= .94;
    
				
  }// if case (2.1)

  
  if (  rgb[1] <= rgb[0] && rgb[1] <= rgb[2]  &&  rgb[2] <= rgb[0]  )  { 
                    // rgb[1] <= rfb[2] <= rgb[0]: case (2.2)
	      
    // Compute reflectance _SampledSpectrum_ with _rgb[1]_ as minimum
	 r[j] += rgb[1] * white;
	 r[j] += (rgb[2] - rgb[1]) * magenta;
	 r[j] += (rgb[0] - rgb[2]) * red;

    r[j] *= .94;

 
  } // if case (2.2)



  if ( rgb[2] <= rgb[0] &&  rgb[2] <= rgb[1] && rgb[0] <= rgb[1] )  { 
                     // rgb[2] <= rgb[0] <= rgb[1]: case (3.1) 

    // Compute reflectance _SampledSpectrum_ with _rgb[2]_ as minimum
	  r[j] += rgb[2] * white;
	  r[j] += (rgb[0] - rgb[2]) * yellow;
	  r[j] += (rgb[1] - rgb[0]) * green; 

     r[j] *= .94;
 
  } // if case (3.1) 


  if (  rgb[2] <= rgb[0] &&  rgb[2] <= rgb[1] && rgb[1] <= rgb[0]  )   { 
                    // rgb[2] <= rgb[1] <= rgb[0]: case (3.2)
			
    // Compute reflectance _SampledSpectrum_ with _rgb[2]_ as minimum
	  r[j] += rgb[2] * white;
	  r[j] += (rgb[1] - rgb[2]) * yellow;
	  r[j] += (rgb[0] - rgb[1]) * red;

      r[j] *= .94;


  } // if case (3.2)
 
 } // for j 

 return r;

} // spectralColorOptimal


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



void main() {
  // access textures in uniform contrl flow

 // access textures in uniform contrl flow

X = XSpectrum();
Y = YSpectrum();
Z = ZSpectrum();


 // for debugging: dump the content of the FBO color texture 
  // get the pixel coordinates of the current fragment

	  ivec2 pixelPos = ivec2( gl_FragCoord.xy);
	      
	 // fragColor =vec4(gl_FragCoord.x, gl_FragCoord.y, 0, 1);
	
	  // get the surface color and depth of the current fragment
	 
	   vec3  surfaceColor = vec3( texelFetch(uColorTex, pixelPos, 0) ); // color

	   float[nSpectralSamples] r;
	   r = spectralColorOptimal( surfaceColor );

	   fragColor = vec4( toRGB( r ), 1.0 );
	   return;


} // main




