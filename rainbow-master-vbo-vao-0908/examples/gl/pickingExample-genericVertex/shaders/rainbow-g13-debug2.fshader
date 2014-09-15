#version 130
uniform sampler1D uPrimarySpectrums;
uniform sampler1D uXYZSpectrums;

uniform vec4 uMaterialColor;
uniform mat4 uProjectionMatrix;

uniform sampler2D uColorTex;
uniform sampler2D uDepthTex;

in vec3 vPosition;
in vec2 vTexCoord;

out vec4 fragColor;

const int nSpectralSamples = 30;

const float  lambdaStart = 400;   // 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
const float lambdaEnd = 700;

const float CIE_Y_integral = 106.856895;

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


float[nSpectralSamples] spectralColorOptimal (vec3 rgb) {
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

vec3 toXYZ(float[nSpectralSamples] spect);
 
vec3 toRGB(float[nSpectralSamples] spect)  {
        vec3  xyz;
        xyz = toXYZ(spect);
        return XYZToRGB(xyz);
    }

vec3  toXYZ(float[nSpectralSamples] spect)  {
        vec3 xyz;
        xyz[0] = xyz[1] = xyz[2] = 0.f;
	

        for (int i = 0; i < nSpectralSamples; ++i) {
            xyz[0] += X[i] * spect[i];
            xyz[1] += Y[i] * spect[i];
            xyz[2] += Z[i] * spect[i];
        }

        float scale = float( lambdaEnd - lambdaStart) /
                          float(CIE_Y_integral * nSpectralSamples);


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

	  ivec2 pixelPos = ivec2( gl_FragCoord.xy);
	      
	 // fragColor =vec4(gl_FragCoord.x, gl_FragCoord.y, 0, 1);
	// fragColor = vec4(1,0,0,1);
	  //return;

	  // get the surface color and depth of the current fragment
	 
	   vec3  surfaceColor = vec3( texelFetch(uColorTex, pixelPos, 0) ); // color

      //vec3  surfaceColor = vec3( texture(uColorTex, vTexCoord) ); // color
	  //float zpixel  = texelFetch(uDepthTex, vTexCoord, 0 // z in non linear range [0,1]

	  float zpixel  = texture(uDepthTex, vTexCoord).r;  // z in non linear range [0,1]
     
      // conversion into NDC [-1,1]gl_FragCoord.z
      float zndc = zpixel * 2.0 - 1.0;

	  float zeye = uProjectionMatrix[3][2]/ ( zndc - uProjectionMatrix[2][2] ); 
	
	  // fragColor = vec4(surfaceColor, 1.0);
	  // return;

	   //float[nSpectralSamples] r = spectralColorOptimal( surfaceColor);
	   float[nSpectralSamples] r = spectralColor( surfaceColor);

	   vec3 rgb =  toRGB( r);
	   fragColor = vec4(rgb,1.0);


} // main
