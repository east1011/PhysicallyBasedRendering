#ifndef MATERIAL_H
#define MATERIAL_H

#include <cstddef>
#include <vector>
#include <map>
#include <memory>
#include <stdexcept>
#if __GNUG__
#   include <tr1/memory>
#endif



//#include "glsupport.h"
#include "uniforms.h"
#include "renderstates.h"

class Geometry; // I do this instead of including "geometry.h", because it will include "material.h" causing a circular include
                // This forward, incomplete declaration is OK, because Geometry is only used as a type for reference or pointer

struct GlProgramDesc;

class MaterialShader {
public:

  const string shaderName; // vertex shader name + fragment shader name

  MaterialShader(const string shaderName, const std::string& vsFilename, const std::string& fsFilename);

  void draw( Geometry& geometry, const Uniforms& extraUniforms);

  Uniforms& getUniforms() { return uniforms_; }
  const Uniforms& getUniforms() const { return uniforms_; }

  RenderStates& getRenderStates() { return renderStates_; }
  const RenderStates& getRenderStates() const { return renderStates_; }

//protected: // by Moon Jung, 2014/8/12 debugging

  std::tr1::shared_ptr<GlProgramDesc> programDesc_;

  Uniforms uniforms_;

  RenderStates renderStates_;
};


#endif
