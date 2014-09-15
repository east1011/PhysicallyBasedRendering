#ifndef OBJECT_H
#define OBJECT_H

#include <vector>
#include <cassert>
#include <map>
#include <cmath>
#include <string>
#include <stdexcept>
#include <memory>

#if __GNUG__
#   include <tr1/memory>
#endif

//#include "cvec.h"
#include "material.h" 

using namespace std;

extern shared_ptr<MaterialShader> g_overridingMaterial;

// by Moon Jung
class Object {
public:
	shared_ptr<ofMatrix4x4> objectRbt;
	shared_ptr<ofMesh> mesh;	
	shared_ptr<MaterialShader> material;

	ofVec4f     objectColor;
	ofVec4f     pickColor;
	
	Object( shared_ptr<ofMatrix4x4> objectRbt, shared_ptr<ofMesh> mesh, shared_ptr<MaterialShader> material ) : 
	  objectRbt ( objectRbt ), mesh (mesh), material  (material) 
	{ }
	  
};

#endif
