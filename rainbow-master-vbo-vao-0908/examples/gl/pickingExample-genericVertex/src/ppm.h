#ifndef PPM_H
#define PPM_H

//#include "ofMain.h"

//#include <QtGui/QOpenGLFunctions_3_0>

#include "GL/glew.h"

void writePpmScreenshot(const int width, const int height, const char *filename);
void printFBO(GLuint g_fbo, int g_windowWidth, int g_windowHeight,  const char * filename );

// A 3-byte structure storing R,G,B value of a pixel
struct PackedPixel {
  GLubyte r,g,b;
};

// The image file is read into `pixels' and its dimension stored into `width'
// and `height'. Throws an exception on error.
void ppmRead(const char *filename, int& width, int& height, std::vector<PackedPixel>& pixels);
void read1DFloatTex(const char *ppmFileName, int& width,   std::vector<float>& floatData );

void read2DFloatTex(const char * ppmFileName, int& width, int& height, 
					std::vector<float>& floatData);

void read3DFloatTex(const char* texFileName, int& texLengthX, int& texLengthY, int& texLengthZ,  
					std::vector<float>& distData );



#endif
