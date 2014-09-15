
#include <vector>

#include "ppm.h"
//#include "glsupport.h"
#include "shadertexture.h"
#include "error.h"

#include "tga.h"

//#include "phaseFunction.h"

//extern float g_phaseFunction[];

using namespace std;

//extern QOpenGLFunctions_3_0 *m_funcs;
extern bool g_Gl2Compatible;


bool LoadTGAFile(char *filename, TGAFILE *tgaFile, int & colorMode)
{
    FILE *filePtr;
    unsigned char ucharBad;
    short int sintBad;
    long imageSize;
    //int colorMode;
    unsigned char colorSwap;

    // Open the TGA file.
    filePtr = fopen(filename, "rb");
    if (filePtr == NULL)
    {
        return false;
    }

    // Read the two first bytes we don't need.
	// size_t fread( void *ptr, size_t size, size_t nitems, FILE *stream);

    fread(&ucharBad, sizeof(unsigned char), 1, filePtr);
    fread(&ucharBad, sizeof(unsigned char), 1, filePtr);

    // Which type of image gets stored in imageTypeCode.
    fread(&tgaFile->imageTypeCode, sizeof(unsigned char), 1, filePtr);

    // For our purposes, the type code should be 2 (uncompressed RGB image)
    // or 3 (uncompressed black-and-white images).
    if (tgaFile->imageTypeCode != 2 && tgaFile->imageTypeCode != 3)
    {
        fclose(filePtr);
        return false;
    }

    // Read 13 bytes of data we don't need.
    fread(&sintBad, sizeof(short int), 1, filePtr); // 1 = count
    fread(&sintBad, sizeof(short int), 1, filePtr);
    fread(&ucharBad, sizeof(unsigned char), 1, filePtr);
    fread(&sintBad, sizeof(short int), 1, filePtr);
    fread(&sintBad, sizeof(short int), 1, filePtr);

    // Read the image's width and height.
    fread(&tgaFile->imageWidth, sizeof(short int), 1, filePtr);
    fread(&tgaFile->imageHeight, sizeof(short int), 1, filePtr);

    // Read the bit depth.
    fread(&tgaFile->bitCount, sizeof(unsigned char), 1, filePtr);

    // Read one byte of data we don't need.
    fread(&ucharBad, sizeof(unsigned char), 1, filePtr);

    // Color mode -> 3 = BGR, 4 = BGRA.
    colorMode = tgaFile->bitCount / 8;
    imageSize = tgaFile->imageWidth * tgaFile->imageHeight * colorMode;

    // Allocate memory for the image data.
    tgaFile->imageData = (unsigned char*)malloc(sizeof(unsigned char)*imageSize);

    // Read the image data.
    fread(tgaFile->imageData, sizeof(unsigned char), imageSize, filePtr);

    // Change from BGR to RGB so OpenGL can read the image data.
    for (int imageIdx = 0; imageIdx < imageSize; imageIdx += colorMode)
    {
        colorSwap = tgaFile->imageData[imageIdx];
        tgaFile->imageData[imageIdx] = tgaFile->imageData[imageIdx + 2];
        tgaFile->imageData[imageIdx + 2] = colorSwap;
    }

    fclose(filePtr);
    return true;
}

ShaderImageTexture1D_RF_RF::ShaderImageTexture1D_RF_RF(const char* ppmFileName): tex() {
  int width;

  vector<float> floatData;

  read1DFloatTex(ppmFileName, width, floatData);

  glBindTexture(GL_TEXTURE_1D, tex); // tex ( member of the class) : tex Handle, set the current texture

  /* no mip map
  if (g_Gl2Compatible)
     glTexParameteri(GL_TEXTURE_1D, GL_GENERATE_MIPMAP, GL_TRUE);
  else  glGenerateMipmap(GL_TEXTURE_1D);
  */

  glTexImage1D(GL_TEXTURE_1D, 0, GL_R32F, width, 
               0, GL_RED, GL_FLOAT, &floatData[0]);
  
 
  //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  //glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  checkGlErrors();
}

ShaderImageTexture1D_RF_RF::ShaderImageTexture1D_RF_RF(const int width, const float *floatData) : tex() {
 

  glBindTexture(GL_TEXTURE_1D, tex); // tex ( member of the class) : tex Handle, set the current texture

  /* no mip map
  if (g_Gl2Compatible)
     glTexParameteri(GL_TEXTURE_1D, GL_GENERATE_MIPMAP, GL_TRUE);
  else  glGenerateMipmap(GL_TEXTURE_1D);
  */

  glTexImage1D(GL_TEXTURE_1D, 0, GL_R32F, width, 
               0, GL_RED, GL_FLOAT, floatData);
  
 
  //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  //glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  checkGlErrors();
}



 ShaderImageTexture2D_RGB_RGB::ShaderImageTexture2D_RGB_RGB(const char* ppmFileName, bool srgb): tex() {
  int width, height;

  vector<PackedPixel> pixData; // PackedPixel = struct { GLubyte r,g,b};
  TGAFILE tga;


  if( strstr(ppmFileName,".ppm") !=NULL){
    ppmRead(ppmFileName, width, height, pixData);

    glBindTexture(GL_TEXTURE_2D, tex); // this->tex : set the current texture this->tex, This was created by the new operator in this case.
	                                   // when this is created its field tex is created by its default constructor.
	
	checkGlErrors();
	// In GL 2: GL_GENERATE_MIPMAP​ is part of the texture object state and it is a flag (GL_TRUE​ or GL_FALSE​). If it is set to GL_TRUE​,  then whenever texture level 0 is updated, 	
	// the mipmaps will all be regenerated.
	// The image being input to a texture is in sRGB format:   This is done using the call glTexImage2D(GL TEXTURE 2D, 0, GL SRGB, twidth, weight, 0, GL RGB, GL UNSIGNED BYTE, pixdata).  
	// Whenever this texture is accessed in a fragment shader, the data is rst converted to linear [R; G;B] coordinates,  before given to the shader. 
   

    glTexImage2D(GL_TEXTURE_2D, 0, (!srgb) || g_Gl2Compatible ? GL_RGB8 : GL_SRGB, width, height,
               0, GL_RGB, GL_UNSIGNED_BYTE, &pixData[0]);

	checkGlErrors();

	if (g_Gl2Compatible) 
        glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
    else glGenerateMipmap(GL_TEXTURE_2D);
	
	
  // The MIN  filter is which quality to show when the texture is near the view, 
  // and the MAG filter is which  quality to show when the texture is far from the view.

    //The qualities are (in order from worst to best)
     //GL_NEAREST
     //GL_LINEAR
     //GL_LINEAR_MIPMAP_NEAREST
     //GL_LINEAR_MIPMAP_LINEAR

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
   // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    checkGlErrors();
  }

  else if(strstr(ppmFileName,".tga")!=NULL){

		int colorMode; 

		LoadTGAFile((char*)ppmFileName , &tga, colorMode);

		glBindTexture(GL_TEXTURE_2D, tex); // tex: tex Handle, set the current texture

		

		if ( colorMode == 4 ) {
		  
		

		  glTexImage2D(GL_TEXTURE_2D, 0, (!srgb) || g_Gl2Compatible ? GL_RGBA : GL_SRGB , tga.imageWidth, tga.imageHeight,
			   0, GL_RGBA, GL_UNSIGNED_BYTE, tga.imageData);
		  
	       if (g_Gl2Compatible)
	    		glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
		  else  glGenerateMipmap( GL_TEXTURE_2D);
		}

		else if ( colorMode == 3) { // colorMode == 3

		 
		  glTexImage2D(GL_TEXTURE_2D, 0, (!srgb) || g_Gl2Compatible ? GL_RGB : GL_SRGB , tga.imageWidth, tga.imageHeight,
			   0, GL_RGB, GL_UNSIGNED_BYTE, tga.imageData);

		 
		   if (g_Gl2Compatible)
	    		glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
		  else  glGenerateMipmap( GL_TEXTURE_2D);	  

		}
		else {
			Error( "16 bit tga texture file is not supported: %s", ppmFileName);
			cout <<  "16 bit tga texture file is not supported: " << ppmFileName << endl;
			messageFile <<  "16 bit tga texture file is not supported: " << ppmFileName << endl;
			
			
		}

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

		//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

		checkGlErrors();


	} // ".tga" file

  else { 
	  Error( "unsupported texture file: %s\n", ppmFileName );
	  cout <<  "unsupported texture file: " << ppmFileName << endl;
	  messageFile <<  "unsupported texture file: " << ppmFileName << endl;


  }


}



ShaderImageTexture2D_RGBF_RGBF::ShaderImageTexture2D_RGBF_RGBF(const char* ppmFileName): tex() { 
// needs to be reimplemented

  int width, height;

  vector<PackedPixel> pixData;

  ppmRead(ppmFileName, width, height, pixData);

  glBindTexture(GL_TEXTURE_2D, tex); // tex ( member of the class) : tex Handle, set the current texture

  
  // which determines the format that OpenGL will use to store the internal
  // texture data. They also take a format and type parameter indicating the
  // format and type of the data SUPPLIED BY  by the application

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, width, height,
               0, GL_RGB, GL_FLOAT, &pixData[0]);

  if (g_Gl2Compatible)
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
  else  glGenerateMipmap(GL_TEXTURE_2D);

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  checkGlErrors();

 

}



ShaderImageTexture2D_RF_RF::ShaderImageTexture2D_RF_RF(const char* texFileName): tex() {

   int height, width;
  vector<float> distData;

 
  read2DFloatTex(texFileName, height, width,  distData);
  
   
  glBindTexture(GL_TEXTURE_2D, tex); // tex ( member of the class ShaderImageTexture3D) : tex Handle, set the current texture

  /* no mip map
  if (g_Gl2Compatible)
    glTexParameteri(GL_TEXTURE_3D, GL_GENERATE_MIPMAP, GL_TRUE);
  else glGenerateMipmap(GL_TEXTURE_3D);
   */

  // We've got a texture name and we've created a texture with the glBindTexture() command, 
  // now we can specify some properties or parameters for the texture

  glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, width, height, 
               0, GL_RED, GL_FLOAT, &distData[0]);



  //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
 

  // GL_LINEAR Returns the weighted average of  the  four  texture
   //                        elements  that  are  closest  to  the center of the
   //                        pixel being textured.   These  can  include  border
   //                        texture   elements,
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  

  checkGlErrors();
}

// The tex member field is initialized by calling constructor GlTexture() with no arguments

ShaderImageTexture3D_RF_RF::ShaderImageTexture3D_RF_RF(const char* texFileName): tex() {
   int depth, height, width;
  vector<float> distData;

 
  read3DFloatTex(texFileName, depth, height, width,  distData);
  
  // let g_phaseFunction points to distData so that it could be transferred to shader as a uniform variable
  // texLength == nRadii; texLengthY == nSpectralSamples; texLengthZ = nThetas
 // for ( int i = 0; i < nRadii * nSpectralSamples * nThetas; i++ ) {
 //      g_phaseFunction[i] = distData[i];
 // }


  
  glBindTexture(GL_TEXTURE_3D, tex); // tex ( member of the class ShaderImageTexture3D) : tex Handle, set the current texture

  //if (g_Gl2Compatible)
  //  glTexParameteri(GL_TEXTURE_3D, GL_GENERATE_MIPMAP, GL_TRUE);
  // else  glGenerateMipmap(GL_TEXTURE_3D);

  // We've got a texture name and we've created a texture with the glBindTexture() command, 
  // now we can specify some properties or parameters for the texture

  glTexImage3D(GL_TEXTURE_3D, 0, GL_R32F, width, height, depth,
               0, GL_RED, GL_FLOAT, &distData[0]);


  //glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  // GL_LINEAR Returns the weighted average of  the  four  texture
   //                        elements  that  are  closest  to  the center of the
   //                        pixel being textured.   These  can  include  border
   //                        texture   elements,
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

  checkGlErrors();
}

// constructor Color Texture which were bound to FBO:
// To default-initialize an object of type T means:
//— if T is a (possibly cv-qualified) class type (Clause 9), the default constructor for T is called (and the initialization is ill-formed if T has no accessible default constructor);
//— if T is an array type, each element is default-initialized;
//— otherwise, no initialization is performed.

// Each class should initialize members that belong to IT, 
// after any parent class constructor call): 

// Here, the constructor is invoked by giving the name of the object to be constructed 
// rather than the name of the class (as in the case of using initialization lists to call the parent class's constructor).
// constant members and reference-valued members also should be always initialized.

Texture_FROM_FBO::Texture_FROM_FBO(GLuint texId): tex(texId)  { 
// needs to be reimplemented

   
  glBindTexture(GL_TEXTURE_2D, tex); // tex ( member of the class) : tex Handle, set the current texture

  /* no mip map
  if (g_Gl2Compatible)
    glTexParameteri(GL_TEXTURE_2D, GL_GENERATE_MIPMAP, GL_TRUE);
   else glGenerateMipmap(GL_TEXTURE_2D);
   */

  // which determines the format that OpenGL will use to store the internal
  // texture data. They also take a format and type parameter indicating the
  // format and type of the data SUPPLIED BY  by the application

  /*
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, width, height,
               0, GL_RGB, GL_FLOAT, &pixData[0]);

 // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  */

  checkGlErrors();

 

}

