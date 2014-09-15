#ifndef TEXTURE_H
#define TEXTURE_H

#include "glsupport.h"
#include "ppm.h"

typedef struct
{
    unsigned char imageTypeCode;
    short int imageWidth;
    short int imageHeight;
    unsigned char bitCount;
    unsigned char *imageData;
} TGAFILE;

//extern bool g_Gl2compatible;

bool LoadTGAFile(char *filename, TGAFILE *tgaFile, int & colorMode);

class ShaderTexture {
public:
  // Must return one of GL_SAMPLER_1D, GL_SAMPLER_2D, GL_SAMPLER_3D, GL_SAMPLER_CUBE,
  // GL_SAMPLER_1D_SHADOW, or GL_SAMPLER_2D_SHADOW, as its intended usage by GLSL shader
  virtual GLenum getSamplerType() const = 0;

  // Binds the texture. (The caller is responsible for setting the active texture unit)
  virtual void bind() const = 0;

  virtual ~ShaderTexture() {}
};

//----------------------------------------
// One concrete implementation of Texture
//----------------------------------------



class ShaderImageTexture1D_RF_RF : public ShaderTexture {
 
	GlTexture tex; // GlTexture constructor is called and generate the texture id

public:
  // Loades a PPM image files with three channels, and create
  // a 2D texture off it. if `srgb' is true, the image is assumed
  // to be in SRGB color space
  ShaderImageTexture1D_RF_RF(const int width, const float *primarySpectrum); // implemented in texture.cpp
  ShaderImageTexture1D_RF_RF( const char  *ppmFileName); // implemented in texture.cpp
  virtual GLenum getSamplerType() const {
    return GL_SAMPLER_1D;
  }

  virtual void bind() const {
    glBindTexture(GL_TEXTURE_1D, tex);
  }
};


class ShaderImageTexture2D_RGB_RGB : public ShaderTexture {
 
	GlTexture tex; // this is simply a declaration. The constructor of tex will be called when an object of type
	               // ShaderImageTexture2D_RGB_RGB is created.

public:
  // Loades a PPM image files with three channels, and create
  // a 2D texture off it. if `srgb' is true, the image is assumed
  // to be in SRGB color space
  ShaderImageTexture2D_RGB_RGB(const char* ppmFileName, bool srgb); // implemented in texture.cpp
 
 
  virtual GLenum getSamplerType() const {
    return GL_SAMPLER_2D;
  }

  virtual void bind() const {
    glBindTexture(GL_TEXTURE_2D, tex);
  }
};



class ShaderImageTexture2D_RGBF_RGBF : public ShaderTexture {
 
	GlTexture tex;

public:
  // Loades a PPM image files with three channels, and create
  // a 2D texture off it. if `srgb' is true, the image is assumed
  // to be in SRGB color space
  ShaderImageTexture2D_RGBF_RGBF(const char* ppmFileName); // implemented in  shadertexture.cpp
  
  virtual GLenum getSamplerType() const { 
    return GL_SAMPLER_2D;
  }

  virtual void bind() const {
    glBindTexture(GL_TEXTURE_2D, tex);
  }
};





class ShaderImageTexture2D_RF_RF : public ShaderTexture {
 
	GlTexture tex;

public:
  // Loades a PPM image files with three channels, and create
  // a 2D texture off it. if `srgb' is true, the image is assumed
  // to be in SRGB color space
  ShaderImageTexture2D_RF_RF(const char* ppmFileName); // implemented in  shadertexture.cpp
  
  virtual GLenum getSamplerType() const { 
    return GL_SAMPLER_2D;
  }

  virtual void bind() const {
    glBindTexture(GL_TEXTURE_2D, tex);
  }
};



class ShaderImageTexture3D_RF_RF : public ShaderTexture {
 
	GlTexture tex;

public:
  // Loades a PPM image files with three channels, and create
  // a 2D texture off it. if `srgb' is true, the image is assumed
  // to be in SRGB color space
  ShaderImageTexture3D_RF_RF(const char* ppmFileName); // implemented in  shadertexture.cpp
  
  virtual GLenum getSamplerType() const { 
    return GL_SAMPLER_3D;
  }

  virtual void bind() const {
    glBindTexture(GL_TEXTURE_3D, tex);
  }
};


class Texture_FROM_FBO : public ShaderTexture {
 
	GlTexture tex; // this is simply a declaration. The constructor of tex will be called when tex is initialized
	               // by a given id or when a texture is created from an image file by  ShaderImageTexture2D_RGB_RGB.

public:
  // Loades a PPM image files with three channels, and create
  // a 2D texture off it. if `srgb' is true, the image is assumed
  // to be in SRGB color space
  Texture_FROM_FBO(GLuint texId); // implemented in texture.cpp
 
 
  virtual GLenum getSamplerType() const {
    return GL_SAMPLER_2D;
  }

  virtual void bind() const {
    glBindTexture(GL_TEXTURE_2D, tex);
  }
};



#endif