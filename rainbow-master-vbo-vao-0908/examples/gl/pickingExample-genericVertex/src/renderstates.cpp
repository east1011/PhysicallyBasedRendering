#include <stdexcept>

#include "glsupport.h"
#include "renderstates.h"

using namespace std;

static const unsigned int kBlendBit = 1,
                          kCullFaceBit = 2;

RenderStates::RenderStates() // default constructor, 
  : glFront(GL_FILL)
  , glBack(GL_FILL)
  , glBlendSrcFactor(GL_ONE)
  , glBlendDstFactor(GL_ZERO)
  , glCullFaceMode(GL_BACK)
  , flags(kCullFaceBit) {}

RenderStates& RenderStates::polygonMode(GLenum face, GLenum mode) {
  switch (mode) {
  case GL_POINT:
  case GL_LINE:
  case GL_FILL:
    switch (face) {
    case GL_FRONT_AND_BACK:
      glFront = glBack = mode;
      return *this;
	  /* deprecated in opengl 3.2 onwards
    case GL_FRONT:
      glFront = mode;
      return *this;
    case GL_BACK:
      glBack = mode;
      return *this;
	  */
    default:
      ;
    }
  default:
    ;
  }

  throw invalid_argument("RenderStates::glPolygonMode: invalid argument");
}

RenderStates& RenderStates::blendFunc(GLenum sfactor, GLenum dfactor) {
  // future TODO: should check if sfactor and dfactor are valid enums
  glBlendSrcFactor = sfactor;
  glBlendDstFactor = dfactor;
  return *this;
}

RenderStates& RenderStates::cullFace(GLenum mode) {
  glCullFaceMode = mode;
  return *this;
}

RenderStates& RenderStates::enable(GLenum target) {
  switch (target) {
  case GL_BLEND:
    flags |= kBlendBit;
    return *this;
  case GL_CULL_FACE:
    flags |= kCullFaceBit;
    return *this;
  default:
    ;
  }
  throw invalid_argument("RenderStates::glEnable: unsupported target");
}

RenderStates& RenderStates::disable(GLenum target) {
  switch (target) {
  case GL_BLEND:
    flags &= ~kBlendBit;
    return *this;
  case GL_CULL_FACE:
    flags &= ~kCullFaceBit;
    return *this;
  default:
    ;
  }
  throw invalid_argument("RenderStates::glEnable: unsupported target");
}

void RenderStates::apply() const {
  static bool firstRun = false; //for some reasons, keep this way. Otherwise, 
                                // the main scene and the OFWindow image do not
                                // blend as expected.

  static RenderStates currentRs; // default constructor is called
  if (firstRun) {
    currentRs.captureFromGl();
    firstRun = false;
  }

  if (this->glFront != currentRs.glFront) {
    //::glPolygonMode(GL_FRONT, glFront);
	  glPolygonMode( GL_FRONT_AND_BACK, this->glFront);
    // A shader can look at the (transformed) normal vector to determine 
	 // whether it is rendering the front or back face. –  Ben Voigt
	  checkGlErrors(); 
	  
	  // If you use GL_FRONT => ENUM error;  GL_FRONT is deprecated in opengl 3.2 onwards
	
    currentRs.glFront = this->glFront;
  }
  if (glBack != currentRs.glBack) {
    ::glPolygonMode(GL_FRONT_AND_BACK, this->glFront);
	  checkGlErrors();
	
    currentRs.glBack = this->glBack;
  }
  if (glBlendSrcFactor != currentRs.glBlendSrcFactor || glBlendDstFactor != currentRs.glBlendDstFactor) {
    ::glBlendFunc(glBlendSrcFactor, glBlendDstFactor);
    currentRs.glBlendSrcFactor = glBlendSrcFactor;
    currentRs.glBlendDstFactor = glBlendDstFactor;
  }

  if (glCullFaceMode != currentRs.glCullFaceMode) {
    ::glCullFace(glCullFaceMode);
    currentRs.glCullFaceMode = glCullFaceMode;
  }

  if ((flags & kBlendBit) != (currentRs.flags & kBlendBit)) {
    if (flags & kBlendBit) {
      ::glEnable(GL_BLEND);
	  checkGlErrors();
	 
	}
    else {
      ::glDisable(GL_BLEND);
	  checkGlErrors();
	}
    currentRs.flags = (currentRs.flags & (~kBlendBit)) | (flags & kBlendBit);
  }

  if ((flags & kCullFaceBit) != (currentRs.flags & kCullFaceBit)) {
    if (flags & kCullFaceBit)
     { ::glEnable(GL_CULL_FACE);
	  checkGlErrors();
	 }
    else
     { ::glDisable(GL_CULL_FACE);
	  checkGlErrors();
	 }
    currentRs.flags = (currentRs.flags & (~kCullFaceBit)) | (flags & kCullFaceBit);
  }
}

void RenderStates::captureFromGl() {
  GLint values[2];

  ::glGetIntegerv(GL_POLYGON_MODE, values);
    checkGlErrors();
  glFront = values[0];
  glBack = values[1];

  ::glGetIntegerv(GL_BLEND_SRC_RGB, values); // Ignore glBlendFuncSeparate fo rnow
    checkGlErrors();
  glBlendSrcFactor = values[0];

  ::glGetIntegerv(GL_BLEND_DST_RGB, values);
    checkGlErrors();
  glBlendDstFactor = values[0];

  ::glGetIntegerv(GL_CULL_FACE_MODE, values);
    checkGlErrors();
  glCullFaceMode = values[0];

  flags = 0;
  if (::glIsEnabled(GL_BLEND))
    flags |= kBlendBit;

    checkGlErrors();

  if (::glIsEnabled(GL_CULL_FACE))
    flags |= kCullFaceBit;

  checkGlErrors();
}


