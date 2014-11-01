#ifndef GEOMETRYMAKER_H
#define GEOMETRYMAKER_H

#include <cmath>

#include "matrix4.h"
#include "geometry.h"  // from pbrt's core
//--------------------------------------------------------------------------------
// Helpers for creating some special geometries such as plane, cubes, and spheres
//--------------------------------------------------------------------------------

typedef Matrix4 SgRbtNode;

// A generic vertex structure containing position, normal, and texture information
// Used by make* functions to pass vertex information to the caller

//  You can store different attributes of vertex into different vectors
//  as follows:
//	vector<ofVec3f> vertices;
//	vector<ofFloatColor> colors;
//	vector<ofVec3f> normals;
//	vector<ofVec2f> texCoords;
//	vector<ofIndexType> indices;

// OR: you can use an general attribute structure and an array of
// these structures as follows:

struct GenericVertex {
  Cvec3f pos;
  Cvec3f normal;
  Cvec2f tex;
  Cvec3f tangent, binormal;

  GenericVertex(
    float x, float y, float z,
    float nx, float ny, float nz,
    float tu, float tv,
    float tx, float ty, float tz,
    float bx, float by, float bz)
    : pos(x,y,z), normal(nx,ny,nz), tex(tu, tv), tangent(tx, ty, tz), binormal(bx, by, bz)
  {}
};

// inline functions and template functions can be defined in headers files
// without the danger of multiple function definition

inline void getPlaneVbIbLen(int& vbLen, int& ibLen) {
  vbLen = 4;
  ibLen = 6;
}

template<typename VtxOutIter, typename IdxOutIter>
void makePlane(float size, VtxOutIter vtxIter, IdxOutIter idxIter) {
  float h = size / 2.0;
  *vtxIter = GenericVertex(    -h, 0, -h, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex(-h, 0,  h, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex( h, 0,  h, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex( h, 0, -h, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, -1);
  *idxIter = 0;
  *(++idxIter) = 1;
  *(++idxIter) = 2;
  *(++idxIter) = 0;
  *(++idxIter) = 2;
  *(++idxIter) = 3;
}

template<typename VtxOutIter, typename IdxOutIter>
void makeNearPlaneQuad(float width, float height, float nearPlane, VtxOutIter vtxIter, IdxOutIter idxIter) {
  float w = width / 2.0;
  float h = height / 2.0;
  *vtxIter = GenericVertex(    -w,  -h, nearPlane, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex(-w,  h, nearPlane, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex( w,  h, nearPlane, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex( w, -h, nearPlane, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, -1);
  *idxIter = 0;
  *(++idxIter) = 1;
  *(++idxIter) = 2;
  *(++idxIter) = 0;
  *(++idxIter) = 2;
  *(++idxIter) = 3;
}

template<typename VtxOutIter, typename IdxOutIter>
void makeNDC(float size, VtxOutIter vtxIter, IdxOutIter idxIter) {
  float h = size / 2.0;
  *vtxIter = GenericVertex(    -h,  -h, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex(-h,  h, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex( h,  h, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 0, -1);
  *(++vtxIter) = GenericVertex( h, -h, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, -1);
  *idxIter = 0;
  *(++idxIter) = 1;
  *(++idxIter) = 2;
  *(++idxIter) = 0;
  *(++idxIter) = 2;
  *(++idxIter) = 3;
}


inline void getCubeVbIbLen(int& vbLen, int& ibLen) {
  vbLen = 24; // 3 * 8
  ibLen = 36; // 6 * 2 triangles * 3 index = 36
}




template<typename VtxOutIter, typename IdxOutIter>

void makeCube(float width, float height, float depth, VtxOutIter vtxIter, IdxOutIter idxIter) {
  float w = width / 2.0;
  float h = height / 2.0;
  float d = depth / 2.0;

#define DEFV(x, y, z, nx, ny, nz, tu, tv) { \
    *vtxIter = GenericVertex(x w, y h, z d, \
                             nx, ny, nz, tu, tv, \
                             tan[0], tan[1], tan[2], \
                             bin[0], bin[1], bin[2]); \
    ++vtxIter; \
}

  Cvec3f tan(0, 1, 0), bin(0, 0, 1); // tangent vector, binormial vector

  DEFV(+, -, -, 1, 0, 0, 0, 0); // facing +X
  DEFV(+, +, -, 1, 0, 0, 1, 0);
  DEFV(+, +, +, 1, 0, 0, 1, 1);
  DEFV(+, -, +, 1, 0, 0, 0, 1);

  tan = Cvec3f(0, 0, 1);
  bin = Cvec3f(0, 1, 0);
  DEFV(-, -, -, -1, 0, 0, 0, 0); // facing -X
  DEFV(-, -, +, -1, 0, 0, 1, 0);
  DEFV(-, +, +, -1, 0, 0, 1, 1);
  DEFV(-, +, -, -1, 0, 0, 0, 1);

  tan = Cvec3f(0, 0, 1);
  bin = Cvec3f(1, 0, 0);
  DEFV(-, +, -, 0, 1, 0, 0, 0); // facing +Y
  DEFV(-, +, +, 0, 1, 0, 1, 0);
  DEFV(+, +, +, 0, 1, 0, 1, 1);
  DEFV(+, +, -, 0, 1, 0, 0, 1);

  tan = Cvec3f(1, 0, 0);
  bin = Cvec3f(0, 0, 1);
  DEFV(-, -, -, 0, -1, 0, 0, 0); // facing -Y
  DEFV(+, -, -, 0, -1, 0, 1, 0);
  DEFV(+, -, +, 0, -1, 0, 1, 1);
  DEFV(-, -, +, 0, -1, 0, 0, 1);

  tan = Cvec3f(1, 0, 0);
  bin = Cvec3f(0, 1, 0);
  DEFV(-, -, +, 0, 0, 1, 0, 0); // facing +Z
  DEFV(+, -, +, 0, 0, 1, 1, 0);
  DEFV(+, +, +, 0, 0, 1, 1, 1);
  DEFV(-, +, +, 0, 0, 1, 0, 1);

  tan = Cvec3f(0, 1, 0);
  bin = Cvec3f(1, 0, 0);
  DEFV(-, -, -, 0, 0, -1, 0, 0); // facing -Z
  DEFV(-, +, -, 0, 0, -1, 1, 0);
  DEFV(+, +, -, 0, 0, -1, 1, 1);
  DEFV(+, -, -, 0, 0, -1, 0, 1);

#undef DEFV

  /* index list for triangles
  for (int v = 0; v < 24; v +=4) {
    *idxIter = v;
    *++idxIter = v + 1;
    *++idxIter = v + 2;
    *++idxIter = v;
    *++idxIter = v + 2;
    *++idxIter = v + 3;
    ++idxIter;
  }

  */
 // index list for rectangles 
  for (int v = 0; v < 24; v +=4) {
    *idxIter = v;
    *++idxIter = v + 1;
    *++idxIter = v + 2;
    *++idxIter = v;
    *++idxIter = v + 2;
    *++idxIter = v + 3;
    ++idxIter;
  }

  

}

// Function template instantiation

//A function template by itself is not a type, or a function, or any other entity. 
// No code is generated from a source file that contains only template definitions. 
//In order for any code to appear, a template must be instantiated: the template arguments must be determined 
//so that the compiler can generate an actual function (or class, from a class template). 

//#include <iostream>
 
//template<typename T>
//void f(T s)
//{
//    std::cout << s << '\n';
//}
 
//int main()
//{
//    f<double>(1); // instantiates and calls f<double>(double)
//    f<>('a'); // instantiates and calls f<char>(char)
//    f(7); // instantiates and calls f<int>(int)
//    void (*ptr)(std::string) = f; // instantiates f<string>(string)
//}

template<typename VtxOutIter, typename IdxOutIter>

void makeAABB(float width, float height, float depth, VtxOutIter vtxIter, IdxOutIter idxIter) {
  float w = width / 2.0;
  float h = height / 2.0;
  float d = depth / 2.0;

#define DEFV(x, y, z, nx, ny, nz, tu, tv) { \
    *vtxIter = GenericVertex(x w, y h, z d, \
                             nx, ny, nz, tu, tv, \
                             tan[0], tan[1], tan[2], \
                             bin[0], bin[1], bin[2]); \
    ++vtxIter; \
}

  Cvec3f tan(0, 1, 0), bin(0, 0, 1); // tangent vector, binormial vector

  DEFV(+, -, -, 1, 0, 0, 0, 0); // facing +X
  DEFV(+, +, -, 1, 0, 0, 1, 0);
  DEFV(+, +, +, 1, 0, 0, 1, 1);
  DEFV(+, -, +, 1, 0, 0, 0, 1);

  tan = Cvec3f(0, 0, 1);
  bin = Cvec3f(0, 1, 0);
  DEFV(-, -, -, -1, 0, 0, 0, 0); // facing -X
  DEFV(-, -, +, -1, 0, 0, 1, 0);
  DEFV(-, +, +, -1, 0, 0, 1, 1);
  DEFV(-, +, -, -1, 0, 0, 0, 1);

  tan = Cvec3f(0, 0, 1);
  bin = Cvec3f(1, 0, 0);
  DEFV(-, +, -, 0, 1, 0, 0, 0); // facing +Y
  DEFV(-, +, +, 0, 1, 0, 1, 0);
  DEFV(+, +, +, 0, 1, 0, 1, 1);
  DEFV(+, +, -, 0, 1, 0, 0, 1);

  tan = Cvec3f(1, 0, 0);
  bin = Cvec3f(0, 0, 1);
  DEFV(-, -, -, 0, -1, 0, 0, 0); // facing -Y
  DEFV(+, -, -, 0, -1, 0, 1, 0);
  DEFV(+, -, +, 0, -1, 0, 1, 1);
  DEFV(-, -, +, 0, -1, 0, 0, 1);

  tan = Cvec3f(1, 0, 0);
  bin = Cvec3f(0, 1, 0);
  DEFV(-, -, +, 0, 0, 1, 0, 0); // facing +Z
  DEFV(+, -, +, 0, 0, 1, 1, 0);
  DEFV(+, +, +, 0, 0, 1, 1, 1);
  DEFV(-, +, +, 0, 0, 1, 0, 1);

  tan = Cvec3f(0, 1, 0);
  bin = Cvec3f(1, 0, 0);
  DEFV(-, -, -, 0, 0, -1, 0, 0); // facing -Z
  DEFV(-, +, -, 0, 0, -1, 1, 0);
  DEFV(+, +, -, 0, 0, -1, 1, 1);
  DEFV(+, -, -, 0, 0, -1, 0, 1);

#undef DEFV

  /* index list for triangles
  for (int v = 0; v < 24; v +=4) {
    *idxIter = v;
    *++idxIter = v + 1;
    *++idxIter = v + 2;
    *++idxIter = v;
    *++idxIter = v + 2;
    *++idxIter = v + 3;
    ++idxIter;
  }

  */
 // index list for quads 
  for (int v = 0; v < 24; v +=4) {
    *idxIter = v;
    *++idxIter = v + 1;
    *++idxIter = v + 2;
    *++idxIter = v + 3;
    ++idxIter;
  }

  

}

// overload makeCube: in this case, the inputs are minCorner and maxCorner points of the box
template<typename VtxOutIter, typename IdxOutIter>

void  DEFV(char x, char y, char z, float nx, float ny, float nz, float tu, float tv, 
		   Cvec3f tan, Cvec3f bin,
		   Cvec3 minCorner, Cvec3 maxCorner,  VtxOutIter vtxIter, IdxOutIter idxIter) {
 float w, h, d;

	if ( x == '+' ) w = maxCorner[0];  
	if ( x == '-' ) w = minCorner[0];  
	if ( y == '+' )  h = maxCorner[1]; 
	if ( y == '-' )  h = minCorner[1]; 
	if ( z == '+')   d = maxCorner[2]; 
	if ( z == '-' )  d = minCorner[2]; 

    *vtxIter = GenericVertex( w,  h,  d, 
                             nx, ny, nz, tu, tv,   
                             tan[0], tan[1], tan[2],   
                             bin[0], bin[1], bin[2]);  
    ++vtxIter; 
 }


template<typename VtxOutIter, typename IdxOutIter>

void makeCube(Cvec3 minCorner, Cvec3 maxCorner, VtxOutIter vtxIter, IdxOutIter idxIter) {

float w; 
float h; 
float d; 


/*
#define DEFV(x, y, z, nx, ny, nz, tu, tv) { \
    *vtxIter = GenericVertex(x w, y h, z d, \
                             nx, ny, nz, tu, tv, \
                             tan[0], tan[1], tan[2], \
                             bin[0], bin[1], bin[2]); \
    ++vtxIter; \
}
*/

  Cvec3f tan(0, 1, 0), bin(0, 0, 1); // tangent vector, binormial vector

  DEFV('+', '-', '-', 1, 0, 0, 0, 0, tan, bin, minCorner, maxCorner, vtxIter, idxIter); // facing '+'X
  DEFV('+', '+', '-', 1, 0, 0, 1, 0, tan, bin, minCorner, maxCorner, vtxIter, idxIter);
  DEFV('+', '+', '+', 1, 0, 0, 1, 1, tan, bin, minCorner, maxCorner, vtxIter, idxIter);
  DEFV('+', '-', '+', 1, 0, 0, 0, 1, tan, bin, minCorner, maxCorner, vtxIter, idxIter);

  tan = Cvec3f(0, 0, 1);
  bin = Cvec3f(0, 1, 0);
  DEFV('-', '-', '-', -1, 0, 0, 0, 0, tan, bin, minCorner, maxCorner, vtxIter, idxIter); // facing '-'X
  DEFV('-', '-', '+', -1, 0, 0, 1, 0, tan, bin, minCorner, maxCorner, vtxIter, idxIter);
  DEFV('-', '+', '+', -1, 0, 0, 1, 1,tan, bin, minCorner, maxCorner, vtxIter, idxIter );
  DEFV('-', '+', '-', -1, 0, 0, 0, 1, tan, bin, minCorner, maxCorner, vtxIter, idxIter);

  tan = Cvec3f(0, 0, 1);
  bin = Cvec3f(1, 0, 0);
  DEFV('-', '+', '-', 0, 1, 0, 0, 0, tan, bin, minCorner, maxCorner, vtxIter, idxIter); // facing '+'Y
  DEFV('-', '+', '+', 0, 1, 0, 1, 0, tan, bin, minCorner, maxCorner, vtxIter, idxIter);
  DEFV('+', '+', '+', 0, 1, 0, 1, 1, tan, bin, minCorner, maxCorner, vtxIter, idxIter);
  DEFV('+', '+', '-', 0, 1, 0, 0, 1, tan, bin, minCorner, maxCorner, vtxIter, idxIter);

  tan = Cvec3f(1, 0, 0);
  bin = Cvec3f(0, 0, 1);
  DEFV('-', '-', '-', 0, -1, 0, 0, 0, tan, bin, minCorner, maxCorner, vtxIter, idxIter); // facing -Y
  DEFV('+', '-', '-', 0, -1, 0, 1, 0, tan, bin, minCorner, maxCorner, vtxIter, idxIter);
  DEFV('+', '-', '+', 0, -1, 0, 1, 1, tan, bin, minCorner, maxCorner, vtxIter, idxIter);
  DEFV('-', '-', '+', 0, -1, 0, 0, 1, tan, bin, minCorner, maxCorner, vtxIter, idxIter);

  tan = Cvec3f(1, 0, 0);
  bin = Cvec3f(0, 1, 0);
  DEFV('-', '-', '+', 0, 0, 1, 0, 0, tan, bin, minCorner, maxCorner, vtxIter, idxIter); // facing '+'Z
  DEFV('+', '-', '+', 0, 0, 1, 1, 0, tan, bin, minCorner, maxCorner, vtxIter, idxIter);
  DEFV('+', '+', '+', 0, 0, 1, 1, 1, tan, bin, minCorner, maxCorner, vtxIter, idxIter);
  DEFV('-', '+', '+', 0, 0, 1, 0, 1, tan, bin, minCorner, maxCorner, vtxIter, idxIter);

  tan = Cvec3f(0, 1, 0);
  bin = Cvec3f(1, 0, 0);
  DEFV('-', '-', '-', 0, 0, -1, 0, 0, tan, bin, minCorner, maxCorner, vtxIter, idxIter); // facing -Z
  DEFV('-', '+', '-', 0, 0, -1, 1, 0, tan, bin, minCorner, maxCorner, vtxIter, idxIter);
  DEFV('+', '+', '-', 0, 0, -1, 1, 1, tan, bin, minCorner, maxCorner, vtxIter, idxIter);
  DEFV('+', '-', '-', 0, 0, -1, 0, 1, tan, bin, minCorner, maxCorner, vtxIter, idxIter);


  for (int v = 0; v < 24; v +=4) {
    *idxIter = v;
    *++idxIter = v + 1;
    *++idxIter = v + 2;
    *++idxIter = v;
    *++idxIter = v + 2;
    *++idxIter = v + 3;
    ++idxIter;
  }

}



inline void getSphereVbIbLen(int slices, int stacks, int& vbLen, int& ibLen) {
  assert(slices > 1);
  assert(stacks >= 2);
  vbLen = (slices + 1) * (stacks + 1);
  ibLen = slices * stacks * 6;
}

template<typename VtxOutIter, typename IdxOutIter>
void makeSphere(float radius, int slices, int stacks, VtxOutIter vtxIter, IdxOutIter idxIter) {
  using namespace std;
  assert(slices > 1);
  assert(stacks >= 2);

  const double radPerSlice = 2 * CS175_PI / slices;
  const double radPerStack = CS175_PI / stacks;

  vector<double> longSin(slices+1), longCos(slices+1);
  vector<double> latSin(stacks+1), latCos(stacks+1);
  for (int i = 0; i < slices + 1; ++i) {
    longSin[i] = sin(radPerSlice * i);
    longCos[i] = cos(radPerSlice * i);
  }
  for (int i = 0; i < stacks + 1; ++i) {
    latSin[i] = sin(radPerStack * i);
    latCos[i] = cos(radPerStack * i);
  }

  for (int i = 0; i < slices + 1; ++i) {
    for (int j = 0; j < stacks + 1; ++j) {
      float x = longCos[i] * latSin[j];
      float y = longSin[i] * latSin[j];
      float z = latCos[j];

      Cvec3f n(x, y, z);
      Cvec3f t(-longSin[i], longCos[i], 0);
      Cvec3f b = cross(n, t);

      *vtxIter = GenericVertex(
        x * radius, y * radius, z * radius,
        x, y, z,
        1.0/slices*i, 1.0/stacks*j,
        t[0], t[1], t[2],
        b[0], b[1], b[2]);
      ++vtxIter;

      if (i < slices && j < stacks ) {
        *idxIter = (stacks+1) * i + j;
        *++idxIter = (stacks+1) * i + j + 1;
        *++idxIter = (stacks+1) * (i + 1) + j + 1;

        *++idxIter = (stacks+1) * i + j;
        *++idxIter = (stacks+1) * (i + 1) + j + 1;
        *++idxIter = (stacks+1) * (i + 1) + j;
        ++idxIter;
      }
    }
  }
}

inline void computeTangentBasis(
	// inputs
	int nverts,
	pbrt::Point *  vertices,
	float * uvs,
	Normal *  normals,
	// outputs
	std::vector<Vector> & tangents,
	std::vector<Vector> & bitangents) {

	for (unsigned int i=0; i< nverts; i+=3 ){

		// Shortcuts for vertices
		pbrt::Point  & v0 = vertices[i+0];
		pbrt::Point & v1 = vertices[i+1];
		pbrt::Point & v2 = vertices[i+2];

		// Shortcuts for UVs
		Cvec2f uv0 (uvs[ 2*i ], uvs[2*i+1] ); 
		Cvec2f uv1 ( uvs[2* i+ 2 ], uvs[2*i +3] );
		Cvec2f uv2  ( uvs[ 2* i+ 4 ], uvs[2*i+ 5] ); 


		// Edges of the triangle : postion delta
		Vector deltaPos1 = v1-v0;
		Vector deltaPos2 = v2-v0;

		// UV delta
		Cvec2f deltaUV1 = uv1-uv0;
		Cvec2f deltaUV2 = uv2-uv0;

		float r = 1.0f / (deltaUV1[0] * deltaUV2[1] - deltaUV1[1] * deltaUV2[0]);
		Vector tangent = (deltaPos1 * deltaUV2[1]   - deltaPos2 * deltaUV1[1])*r;
		Vector  bitangent = (deltaPos2 * deltaUV1[0]   - deltaPos1 * deltaUV2[0])*r;

		// Set the same tangent for all three vertices of the triangle.
		// They will be merged later, in vboindexer.cpp
		tangents.push_back(tangent);
		tangents.push_back(tangent);
		tangents.push_back(tangent);

		// Same thing for binormals
		bitangents.push_back(bitangent);
		bitangents.push_back(bitangent);
		bitangents.push_back(bitangent);

	}

	// See "Going Further"
	for (unsigned int i=0; i< nverts; i+=1 )
	{
		Normal  & n =  normals[i];  // n is the "reference" to normals[i]
		Vector  & t = tangents[i];
		Vector  & b = bitangents[i];
		
		// Gram-Schmidt orthogonalize: to make the tangent orthogonal to normal
		t = Normalize(t - Vector(n) * Dot( n, t));  
		// t - Vector(n) * Dot(n,t) can be zero, which means that n is almost parallel to t, so that the normalization fails.
		
		if (   ( t - Vector(n) * Dot( n, t)  ).Length()   < 1.0e-7 ) continue;

		// Calculate handedness
		if ( Dot(Cross(n, t), b) < 0.0f){
			t = t * -1.0f;
		}

	}


} // ComputeTangentBasis


#endif
