#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* malloc, free, rand */

//#include <GL/glew.h>

#include "picker_0.h"
#include "ofMain.h"

//extern QOpenGLFunctions_3_0 *m_funcs;
extern bool g_Gl2Compatible;

using namespace std;
using namespace std::tr1;

Picker::Picker()
	:idCounter_(0),  srgbFrameBuffer_(!g_Gl2Compatible) {}
		

shared_ptr<Object> Picker::getRbtNodeAtXY(int x, int y) {
  // TODO

  //unsigned char pixel[3];
  //GLubyte pixel[3];

   int width = 1;
   int height = 1;

   //GLubyte *pixel = reinterpret_cast<GLubyte *> ( malloc( 3 * width * height) );

   //PackedPixel  pixel; // PackedPixel= a structure of r, g, b bytes
 // GLubyte pixel[3];
   unsigned char pixel[4] = {0};


  // read the pixel at (x,y), get the color of the pixel
  // get the object id corresponding to the color
  // return the Rbt node corressponding to the ide

   //vector<float> rcolors(4);
  // glReadPixels(x, y, 1, 1, GL_RGBA, GL_FLOAT, &rcolors[0]);

 // for debugging

  glReadBuffer(GL_BACK); // default
   
  glReadPixels(x,y, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixel); //  GLvoid *  pixel 
  // This read pixels from the back buffer

  /* More generally, do the following:
  GLint format, type;
  glGetIntegerv(GL_IMPLEMENTATION_COLOR_READ_FORMAT, &format);
  glGetIntegerv(GL_IMPLEMENTATION_COLOR_READ_TYPE, &type);
  glReadPixels(x, y, width, height, format, type, color);
  */
   
  checkGlErrors();

  std::cout << "picked position: x,y=" << x <<"," << y << endl;
  std::cout << "Picked [Integer] Color (sRGBA or RGB space depending on GL version)=" << (int) pixel[0] << "," << (int) pixel[1] << "," << (int) pixel[2] <<"," << (int) pixel[3] << endl;
   
  unsigned int id = colorToId( pixel ); // color to integer conversion: pixel color is in sRGB space

  std::cout << "Picked object ID =" << id << endl;

  //Cvec3 color = idToColor(id); // the color components will be printed

  
  return  find(id);
  
}

//------------------
// Helper functions
//------------------
//

int Picker::getidCounter() {
	  return  idCounter_;
}


//------------------
// Helper functions
//------------------
//
void Picker::addToMap(int id, shared_ptr<Object> node) {

	idToRbtNode_[id] = node;
}

shared_ptr<Object> Picker::find(int id) {
  
	IdToRbtNodeMap::iterator it = idToRbtNode_.find(id);
  
  if (it != idToRbtNode_.end())
    
	  return it->second;
  else
    return nullptr; // set to null
}




// encode 255^3  IDs  corresponding to each RGB color;

static const int NBITS = 8, N = 1 << NBITS, MASK = N-1;
//     1 << 8 == N=100000000; MASK = 126 - 1 = 255 = 011111111
//     id & MASK => the last 8 bits of id


Cvec4 Picker::idToColor(unsigned int id) { //Id to (R, G, B,A) color


 // assert( ("id must ge greater than 0 and less than 255^3 in Picker::idToColor", id > 0 && id < N * N * N) ); 
	
   assert( id > 0 && id < N * N * N ); 

  // Integer Id has 32 bits (8 * 4 bits): The frist byte becomes R, the second  byte G, the third  byte B,
  //                                      the last byte A.                                        
  
 // Cvec4 linearRGBAColor = Cvec4( ( id <<  3* NBITS ) & MASK, (id <<  2*NBITS) & MASK, 
//	                                (id <<  NBITS) & MASK, id & MASK ) / 255;



   Cvec4 linearRGBAColor = Cvec4( ( id >>   3* NBITS ) & MASK, (id >>  2* NBITS) & MASK, 
	   (id >> NBITS) & MASK, id & MASK ); 

   std::cout << "Id => int linearRGBA : " << id  <<"=>" << linearRGBAColor[0] << "," << linearRGBAColor[1] << ","  
	   << linearRGBAColor[2] << "," << linearRGBAColor[3] << endl;

  Cvec4 floatLinearRGBAColor = linearRGBAColor / 255.0;

  //std::cout << "Id => float linearRGBA : " << id  <<"=>" << floatLinearRGBAColor[0] << "," << floatLinearRGBAColor[1] << ","  
  //	  << floatLinearRGBAColor[2] << "," << floatLinearRGBAColor[3] << endl;


  //return floatLinearRGBAColor;
   
  
  // srgbFrameBuffer is a member of Picker class
  
  if (!srgbFrameBuffer_) { // g_Gl2Compatible == !srgbFrameBuffer; an older version Gl2 of OpenGL is used

   std::cout << "Id => int linearRGBA : " << id  <<"=>" << linearRGBAColor << endl;
    return linearRGBAColor; 
  }


  else { // !g_Gl2Compatible => GL 3 is used:
    // if GL3 is used, the framebuffer will be in SRGB format, and the color we supply is automatically converted to SRGB format in the 
	 // process of rendering.
	// So, to find the color in the framebuffer which corresponds to a given linear color, we need to convert this color into 
    // sRGB format.


	Cvec4 sRGBColor;
 
	for (int i = 0; i < 3; ++i) {
      sRGBColor[i] = floatLinearRGBAColor[i] <= 0.04045 ? floatLinearRGBAColor[i]/12.92 : pow( ( floatLinearRGBAColor[i] + 0.055 )/1.055, 2.4);
    }

	sRGBColor[3] = floatLinearRGBAColor[3]; // alpha not changed

	std::cout << "  Id => float linear RGBA : " << floatLinearRGBAColor  << endl;
	
	std::cout << " Id => sRGBA Color : " << id  <<"=>" << sRGBColor << endl;
    return sRGBColor;
  } 

}

   
 // (1)  Gamma-correction (for CRT monitor):
 // The strength of the electron beam determined the brightness of that part of the image. 
 // However, the strength of the beam did not vary linearly with the brightness of the image.
 // The easiest way to deal with that in the earliest days of TV was to simply modify the incoming image at the source.
// TV broadcasts sent image data that was non-linear in the opposite direction of the CRT's normal non-linearity.
// Thus, the final output was displayed linearly, as it was originally captured by the camera.

// The term for this process, de-linearizing an image to compensate for a non-linear display, is called gamma correction.

// You may be wondering why this matters. After all, odds are, you do not have a CRT-based monitor;
// you probably have some form of LCD, plasma, LED, or similar technology.
//   So what does the vagaries of CRT monitors matter to you?

// Because gamma correction is everywhere. It's in DVDs, video-tapes, and Blu-Ray discs. 
// Every digital camera does it. This is how it has been for a long time. 
// Because of that, you could not sell an LCD monitor that tried to do linear color reproduction; 
// nobody would buy it because all media for it (including your OS) was designed 
//   and written expecting CRT-style non-linear displays.

/*
  (2) When a texture uses one of the sRGB formats, and texture access functions they fetch a texel, 
      OpenGL automatically linearizes the color from the sRGB colorspace. This is exactly what we want. 
      The best part is that the linearisation cost is negligible. So there is no need to play with the data 
	  or otherwise manually linearize it.  OpenGL does it for us. Note that the shader does not change. 
	  It still uses a regular sampler2D, accesses it with a 2D texture coordinate and the texture function, etc.
	  The shader does not have to know or care  whether the image data is in the sRGB colorspace or a linear one.
	  It simply calls the texture function and expects it to return lRGB color values.
*/


// (3) convert srgb color space to a linear color space:
	  // The sRGB colorspace is an approximation of a gamma of 2.2. It is not exactly 2.2, 
	  // but it is close enough that you can display an sRGB image to the screen without gamma correction.

	  // OpenGL defines sRGB textures (OpenGL 2.1) and sRGB framebuffers (OpenGL 3.0) to take advantage of the sRGB color space.
	  // sRGB is not a linear color space. So, if we have 2 sRGB pixels A and B then A + B is not a valid operation.
	  // Each time a sRGB value is used for an operation, the value need to be converted to a linear color space value. 
	  // To convert a 8 bits sRGB pixel to a linear RGB pixel, the graphics card might need 16 bits per channel to be lossless. 

    // On the OpenGL pipeline we will find these conversions for three tasks: texel sampling, multisample resolution and blending. 
	// OpenGL gives control of these conversions only at blending stage with the calls glEnable(GL_FRAMEBUFFER_SRGB) and glDisable(GL_FRAMEBUFFER_SRGB).
	// When enabled, the values are linearized prior to their use for blending. A lack of conversions will generally darken the framebuffer. 
	// Unfortunatly, OpenGL ES 1 and 2 don't provide any support of sRGB. 

	
/*
(4) Most of the inputs to the rendering pipeline are in sRGB space, as well as (almost always) 
the final output. But we want most of the stuff in the middle to be in linear space, for the color computatation. 
So we end up doing a bunch of conversions between color spaces. Such conversions are done using 
the following functions, implemented in hardware (usually in approximate form), in shader code or in tools code:

f_{lin => sRGB} (x) = 12.92 x if x <= 0.00031308; 1.055x^{1/2.4} - 0.055 otherwise
f_{sRGB => lin }(x) = x/12.92 if x <= 0.04045; ([ x + 0.055 ]/1.055  )^ 2.4 otherwise
input, output in (0,1)

Note that the alpha value itself represents a physical coverage amount, thus it is always 
in linear space and never has to be converted.

*/

/*
void printBitString( string&  s ) {
for (int i = 0; i < s.length(); i++)
    for (char c = 1; c; c <<= 1) // little bits first
        std::cout << (s[i] & c ? "1" : "0");
}
*/

unsigned int Picker::colorToId(const unsigned char p[] ) {
  
	const int  NBITS = 8;
  
  //cout << "Picked [Integer] color (r,g,b,a)= in colorToId " << (unsigned int) p[0] << "," << 
  //	                (unsigned int) p[1] << "," << (unsigned int) p[2] << "," << (unsigned int) p[3] << endl;

  
  unsigned int id1 = (unsigned int) p[0]; 
  unsigned  int id2 = (unsigned int) p[1]; 
  unsigned int id3 = (unsigned int) p[2]; 
  unsigned int id4 = (unsigned int) p[3];
 
  cout << "id1 =" << id1 << endl;
  cout << "id2 =" << id2 << endl;
  cout << "id3 =" << id3 << endl;
  cout << "id4=" << id4 << endl;

  unsigned int  id;

  id = ( id1 << (3*NBITS) ) | ( id2 << (2* NBITS) ) | ( id3 << NBITS ) | id4;
    
 
  //cout << "id for the given color=" << (unsigned int) id << endl;

  return id;
}

 