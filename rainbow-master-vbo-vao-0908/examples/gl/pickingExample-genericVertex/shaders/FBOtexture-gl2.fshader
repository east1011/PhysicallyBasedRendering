#version 120

uniform int uWindowWidth; 
uniform int uWindowHeight; 

uniform sampler1D uPrimarySpectrums;
uniform sampler1D uXYZSpectrums;

uniform vec4 uMaterialColor;
uniform mat4 uProjectionMatrix;

uniform sampler2D uColorTex;
uniform sampler2D uDepthTex;

varying  vec3 vPosition;
varying  vec2 vTexCoord;

//out vec4 fragColor;

const int nSpectralSamples = 10;
const int nSamples =10;

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
     float pos = ( float( 0 * nSpectralSamples + j)  + 0.5 )/ float ( nSpectralSamples * 3);
     r[j] = texture1D(uXYZSpectrums,  pos).r; 
  }
  return r;
}

float[nSpectralSamples] YSpectrum() {
  float[nSpectralSamples]  r;
  for (int j=0; j < nSpectralSamples; j++ ) {
     float pos =( float( 1 * nSpectralSamples + j ) + 0.5 )/ float ( nSpectralSamples * 3);
     r[j] = texture1D( uXYZSpectrums, pos).r; 
  }
  return r;
}

float[nSpectralSamples] ZSpectrum() {
  float[nSpectralSamples]  r;
  for (int j=0; j < nSpectralSamples; j++ ) {
     float pos =( float( 2 * nSpectralSamples + j ) + 0.5 )/ float ( nSpectralSamples * 3);
     r[j] = texture1D( uXYZSpectrums, pos).r; 
  }
  return r;
}

float[nSpectralSamples] whiteSpectrum() {
  float[nSpectralSamples]  r;
  for (int j=0; j < nSpectralSamples; j++ ) {
     float pos =( float( 0 * nSpectralSamples + j ) + 0.5 )/ float ( nSpectralSamples * 7);
     r[j] = texture1D( uPrimarySpectrums, pos).r; 
  }
  return r;
}
float[nSpectralSamples] cyanSpectrum() {
  float[nSpectralSamples]  r;
  for (int j=0; j < nSpectralSamples; j++ ) {
     float pos =( float( 1 * nSpectralSamples + j ) + 0.5 )/ float ( nSpectralSamples * 7);
     r[j] = texture1D(uPrimarySpectrums, pos).r; 
  }
  return r;
}
float[nSpectralSamples] magentaSpectrum() {
  float[nSpectralSamples]  r;
  for (int j=0; j < nSpectralSamples; j++ ) {
     float pos =( float( 2 * nSpectralSamples + j ) + 0.5 )/ float ( nSpectralSamples * 7);
     r[j] = texture1D(uPrimarySpectrums, pos).r; 
  }
  return r;
}
float[nSpectralSamples] yellowSpectrum() {
  float[nSpectralSamples]  r;
  for (int j=0; j < nSpectralSamples; j++ ) {
     float pos =( float( 3 * nSpectralSamples + j ) + 0.5 )/ float ( nSpectralSamples * 7);
     r[j] = texture1D(uPrimarySpectrums, pos).r; 
  }
  return r;
}

float[nSpectralSamples] redSpectrum() {
  float[nSpectralSamples]  r;
  for (int j=0; j < nSpectralSamples; j++ ) {
     float pos =( float( 4 * nSpectralSamples + j ) + 0.5 )/ float ( nSpectralSamples * 7);
     r[j] = texture1D(uPrimarySpectrums, pos).r; 
  }
  return r;
}
float[nSpectralSamples] greenSpectrum() {
  float[nSpectralSamples]  r;
  for (int j=0; j < nSpectralSamples; j++ ) {
     float pos =( float( 5 * nSpectralSamples + j ) + 0.5 )/ float ( nSpectralSamples * 7);
     r[j] = texture1D(uPrimarySpectrums, pos).r; 
  }
  return r;
}

float[nSpectralSamples] blueSpectrum() {
  float[nSpectralSamples]  r;
  for (int j=0; j < nSpectralSamples; j++ ) {
     float pos =( float( 6 * nSpectralSamples + j ) + 0.5 )/ float ( nSpectralSamples * 7);
     r[j] = texture1D(uPrimarySpectrums, pos).r; 
  }
  return r;
}



float[nSamples] spectralColor(vec3 rgb) {
// primarySpectrums
// As background -- GLSL looks a lot like C, but compiles a bit different. 
//Things are very unrolled, and conditionals may be executed in parallel and switched at the end,
// that sort of thing. Depends on the hardware...

//You can use **loop indices or constants** [== contant-index-expression] to index into arrays. 
//The assignment in your loop is ok, but the access by tileID isn't.

  float[nSamples]  r;
  for (int j=0; j < nSamples; j++ ) {
     r[j] = 0.0;
  }

 
  // Convert reflectance spectrum to RGB
  
  if (rgb[0] <= rgb[1] && rgb[0] <= rgb[2]  && rgb[1] <= rgb[2] ) {
                          // rgb[0] <= rgb[1] << rgb[2], case (1.1)	  

   // Compute reflectance _SampledSpectrum_ with _rgb[0]_ as minimum

     for (int j =0; j < nSamples; j++) {
          r[j] += rgb[0] * rgbRefl2SpectWhite[j];

	 } // for

	 
    for (int j =0; j < nSamples; j++) {
               r[j] += (rgb[1] - rgb[0]) * rgbRefl2SpectCyan[j];		
			  
               r[j] += (rgb[2] - rgb[1]) * rgbRefl2SpectBlue[j];
	
    } // for
		 
		 
	for (int j =0; j < nSamples; j++) {
               r[j] *= .94;
    }

    return r; 

  } // if case (1.1)
  
  // return r;

  if ( rgb[0] <= rgb[1] && rgb[0] <= rgb[2]  && rgb[2] <= rgb[1]   ) { 
                      // rgb[0] << rgb[2] << rgb[1]: case (1.2)

   // Compute reflectance _SampledSpectrum_ with _rgb[0]_ as minimum

     for (int j =0; j < nSamples; j++) {
          r[j] += rgb[0] * rgbRefl2SpectWhite[j];

	 } // for

	      for (int j =0; j < nSamples; j++) {
                r[j] += (rgb[2] - rgb[0]) * rgbRefl2SpectCyan[j];
		
                r[j] += (rgb[1] - rgb[2]) * rgbRefl2SpectGreen[j];
			
          } // for

        
		  for (int j =0; j < nSamples; j++) {
               r[j] *= .94;
          }

          return r; 
		  
  } //  case (1.2)




  if (rgb[1] <= rgb[0] && rgb[1] <= rgb[2]  &&  rgb[0] <= rgb[2] ) { // case (2.1)
                      // rgb[1] <= rgb[0] <= rgb[2]: case (2.1)
	   
    // Compute reflectance _SampledSpectrum_ with _rgb[1]_ as minimum
	   for (int j =0; j < nSamples; j++) {
           r[j] += rgb[1] * rgbRefl2SpectWhite[j];

	   } // for

   
	        for (int j =0; j < nSamples; j++) {
                r[j] += (rgb[0] - rgb[1]) *rgbRefl2SpectMagenta[j];
			
                r[j] += (rgb[2] - rgb[0]) *rgbRefl2SpectBlue[j];
		
			} // for
		
		    for (int j =0; j < nSamples; j++) {
               r[j] *= .94;
            }

            return r; 
			
  }// if case (2.1)




  if (  rgb[1] <= rgb[0] && rgb[1] <= rgb[2]  &&  rgb[2] <= rgb[0]  )  { 
                    // rgb[1] <= rfb[2] <= rgb[0]: case (2.2)
	      
    // Compute reflectance _SampledSpectrum_ with _rgb[1]_ as minimum
	   for (int j =0; j < nSamples; j++) {
           r[j] += rgb[1] * rgbRefl2SpectWhite[j];

	   } // for

            for (int j =0; j < nSamples; j++) { 
			  r[j] += (rgb[2] - rgb[1]) * rgbRefl2SpectMagenta[j];

              r[j] += (rgb[0] - rgb[2]) * rgbRefl2SpectRed[j];
			
			} // for

			 // return rgbRefl2SpectWhite;
		    for (int j =0; j < nSamples; j++) {
               r[j] *= .94;
            }

            return r; 
  } // if case (2.2)



  if ( rgb[2] <= rgb[0] &&  rgb[2] <= rgb[1] && rgb[0] <= rgb[1] )  { 
                     // rgb[2] <= rgb[0] <= rgb[1]: case (3.1) 

    // Compute reflectance _SampledSpectrum_ with _rgb[2]_ as minimum
	   for (int j =0; j < nSamples; j++) {
            r[j] += rgb[2] * rgbRefl2SpectWhite[j];
		} // for
       
	      for (int j =0; j < nSamples; j++) {
            r[j] += (rgb[0] - rgb[2]) * rgbRefl2SpectYellow[j];
	       
            r[j] += (rgb[1] - rgb[0]) * rgbRefl2SpectGreen[j]; 

	      } // for
		 
		  for (int j =0; j < nSamples; j++) {
               r[j] *= .94;
          }

          return r; 
  } // if case (3.1) 


  



  if (  rgb[2] <= rgb[0] &&  rgb[2] <= rgb[1] && rgb[1] <= rgb[0]  )   { 
                    // rgb[2] <= rgb[1] <= rgb[0]: case (3.2)
			
    // Compute reflectance _SampledSpectrum_ with _rgb[2]_ as minimum
	   for (int j =0; j < nSamples; j++) {
            r[j] += rgb[2] * rgbRefl2SpectWhite[j];
		} // for
		   
	       for (int j =0; j < nSamples; j++) {

               r[j] += (rgb[1] - rgb[2]) * rgbRefl2SpectYellow[j];
			   r[j] += (rgb[0] - rgb[1]) * rgbRefl2SpectRed[j];
			   
		   } // for

		    // return rgbRefl2SpectWhite;
		   for (int j =0; j < nSamples; j++) {
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

	  vec2 pixelPos = vec2( gl_FragCoord.x / uWindowWidth, gl_FragCoord.y / uWindowHeight ) ;
	 	      

	  // get the surface color and depth of the current fragment
	 //gl_FragCoord contents:http://www.txutxi.com/?p=182

     //The first two values (x,y) contain the pixel¡¯s center coordinates where the fragment is being rendered. 
	 //For instance, with a frame buffer resolution of 800¡¿600, a fragment being rendered in the bottom-left corner
	 // would fall over the pixel position (0.5,0.5); 
	 //a fragment rendered into the most top-right corner would have coordinates (799.5, 599.5).
	 // For applications that use multi-sampling these values may fall elsewhere on the pixel area.

	  // vec3  surfaceColor = vec3( texelFetch(uColorTex, pixelPos, 0) ); // color

	  vec3  surfaceColor = vec3( texture2D(uColorTex, pixelPos) ); // color

      //vec3  surfaceColor = vec3( texture(uColorTex, vTexCoord) ); // color
	 
	  float zpixel  = texture2D(uDepthTex, pixelPos).r;  // z in non linear range [0,1]
     
      // conversion from range [0,10 to NDC [-1,1]
      float zndc = zpixel * 2.0 - 1.0;

	  float zeye = uProjectionMatrix[3][2]/ ( zndc - uProjectionMatrix[2][2] ); 
	 

	  //fragColor = vec4(surfaceColor, 1);
	 // return;

	 
	   float[nSamples] r = spectralColor( surfaceColor);

	   vec3 rgb =  toRGB( r);
	   gl_FragColor = vec4(rgb,1.0);

//fragColor = uMaterialColor;
//fragColor = vec4(0, 0, 200/255.0, 1.0);


}
