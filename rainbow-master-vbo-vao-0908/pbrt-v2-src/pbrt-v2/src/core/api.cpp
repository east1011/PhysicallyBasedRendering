
/*
pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

This file is part of pbrt.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

//#include "file":
//1.In the same directory as the file that contains the #include statement. 
//2.In the directories of any previously opened include files in the reverse order in which they were opened. 
// The search starts from the directory of the include file that was opened last and continues 
// through the directory of the include file that was opened first. 



// core/api.cpp*
#include "stdafx.h"
#include "api.h"
#include "parallel.h"
#include "paramset.h" // includes Reference class 
#include "spectrum.h"
#include "scene.h"
#include "renderer.h"
#include "film.h"
#include "volume.h"
#include "probes.h"
#include "material.h"

// API Additional Headers
#include "accelerators/bvh.h"
#include "accelerators/grid.h"
#include "accelerators/kdtreeaccel.h"
#include "cameras/environment.h"
#include "cameras/orthographic.h"
#include "cameras/perspective.h"
#include "film/image.h"
#include "filters/box.h"
#include "filters/gaussian.h"
#include "filters/mitchell.h"
#include "filters/sinc.h"
#include "filters/triangle.h"
#include "integrators/ambientocclusion.h"
#include "integrators/diffuseprt.h"
#include "integrators/dipolesubsurface.h"
#include "integrators/directlighting.h"
#include "integrators/emission.h"
#include "integrators/glossyprt.h"
#include "integrators/igi.h"
#include "integrators/irradiancecache.h"
#include "integrators/path.h"
#include "integrators/photonmap.h"
#include "integrators/single.h"
#include "integrators/useprobes.h"
#include "integrators/whitted.h"
#include "lights/diffuse.h"
#include "lights/distant.h"
#include "lights/goniometric.h"
#include "lights/infinite.h"
#include "lights/point.h"
#include "lights/projection.h"
#include "lights/spot.h"
#include "materials/glass.h"
#include "materials/kdsubsurface.h"
#include "materials/matte.h"
#include "materials/measured.h"
#include "materials/metal.h"
#include "materials/mirror.h"
#include "materials/mixmat.h"
#include "materials/plastic.h"
#include "materials/substrate.h"
#include "materials/subsurface.h"
#include "materials/shinymetal.h"
#include "materials/translucent.h"
#include "materials/uber.h"
#include "renderers/aggregatetest.h"
#include "renderers/createprobes.h"
#include "renderers/metropolis.h"
#include "renderers/samplerrenderer.h"
#include "renderers/surfacepoints.h"
#include "samplers/adaptive.h"
#include "samplers/bestcandidate.h"
#include "samplers/halton.h"
#include "samplers/lowdiscrepancy.h"
#include "samplers/random.h"
#include "samplers/stratified.h"
#include "shapes/cone.h"
#include "shapes/cylinder.h"
#include "shapes/disk.h"
#include "shapes/heightfield.h"
#include "shapes/hyperboloid.h"
#include "shapes/loopsubdiv.h"
#include "shapes/nurbs.h"
#include "shapes/paraboloid.h"
#include "shapes/sphere.h"
#include "shapes/trianglemesh.h"
#include "textures/bilerp.h"
#include "textures/checkerboard.h"
#include "textures/constant.h"
#include "textures/dots.h"
#include "textures/fbm.h"
#include "textures/imagemap.h"
#include "textures/marble.h"
#include "textures/mix.h"
#include "textures/scale.h"
#include "textures/uv.h"
#include "textures/windy.h"
#include "textures/wrinkled.h"
#include "volumes/exponential.h"
#include "volumes/homogeneous.h"
#include "volumes/volumegrid.h"
#include <map>
#if (_MSC_VER >= 1400)
#include <stdio.h>
#define snprintf _snprintf
#endif
using std::map;

// API Global Variables
Options PbrtOptions;

////////////////////////////////////////////////////


char g_stringFileName[100];

/////////////////////////////////////
//class OpenGLWindow;
//#include <QtGui/QWindow>
//#include <QtGui/QOpenGLFunctions_3_0>
//#include <QtGui/QResizeEvent>
//#include <QtGui/QGuiApplication>

#include "glsupport.h"
//#include "glu.h"
#include "materialshader.h"
#include "shadergeometry.h" // includes Object class
#include "geometrymaker.h"
#include "picker_0.h"
#include "ppm.h"

#include "uniforms.h"





int render_flag;


extern double g_frustMinFov;  // A minimal of 60 degree field of view
extern double g_frustFovY; // FOV in y direction (updated by updateFrustFovY)

extern  double g_frustNear;    // near plane
extern  double g_frustFar;    // far plane


extern int g_windowWidth;
extern int g_windowHeight;


extern GLdouble g_clearColor[4];
// --------- Scene
extern Cvec4 g_light1Pos, g_light2Pos;  // define two lights positions in world space
extern Cvec4 g_light1Color;
extern Cvec4 g_light2Color;
extern Matrix4 g_eyeRbt ; // ClASS METHOD CALL => AN INSTANCE OF MATRIX4


// shared_ptr< Picker > g_picker_ptr;
// g_picker_ptr.reset( new Picker() );

extern Picker g_picker; // used to define pick colors
extern shared_ptr< Object>  g_currentPickedObject;
// references are not objects, i.e. types regions of storage; so they are not actually "variables";
// you cannot point to references; so there are no references of references, array of references, pointers to
// references. They are just aliases of other variables. So, you cannot use vector < Object & > objectList

extern vector < shared_ptr<Object> >  g_objectPtrList;

shared_ptr<MaterialShader> g_uberMat; 
shared_ptr<MaterialShader> g_matteMat, g_plasticMat; 

shared_ptr<MaterialShader> g_currentMaterial; 


///
///////////////////////////////////

// API Local Classes
#define MAX_TRANSFORMS 2
#define START_TRANSFORM_BITS (1 << 0)
#define END_TRANSFORM_BITS   (1 << 1)
#define ALL_TRANSFORMS_BITS  ((1 << MAX_TRANSFORMS) - 1)
struct TransformSet {
	// TransformSet Public Methods
	Transform &operator[](int i) {
		Assert(i >= 0 && i < MAX_TRANSFORMS);
		return t[i];
	}
	const Transform &operator[](int i) const { Assert(i >= 0 && i < MAX_TRANSFORMS); return t[i]; }
	friend TransformSet Inverse(const TransformSet &ts) {
		TransformSet t2;
		for (int i = 0; i < MAX_TRANSFORMS; ++i)
			t2.t[i] = Inverse(ts.t[i]);
		return t2;
	}
	bool IsAnimated() const {
		for (int i = 0; i < MAX_TRANSFORMS-1; ++i)
			if (t[i] != t[i+1]) return true;
		return false;
	}
private:
	Transform t[MAX_TRANSFORMS];
};


struct RenderOptions {
	// RenderOptions Public Methods
	RenderOptions();
	Scene *MakeScene();
	Camera *MakeCamera() const;
	Renderer *MakeRenderer() const;

	void    GPURender(); // by Moon Jung, 2014/4/4

	// RenderOptions Public Data
	float transformStartTime, transformEndTime;
	string FilterName;
	ParamSet FilterParams;
	string FilmName;
	ParamSet FilmParams;

	string SamplerName;
	ParamSet SamplerParams;
	string AcceleratorName;
	ParamSet AcceleratorParams;
	string RendererName;
	string SurfIntegratorName, VolIntegratorName;
	ParamSet RendererParams;
	ParamSet SurfIntegratorParams, VolIntegratorParams;
	string CameraName;
	ParamSet CameraParams;
	TransformSet CameraToWorld;
	vector<Light *> lights;

	vector<Reference<Primitive> > primitives; // equiv to g_objectPtrList 

	mutable vector<VolumeRegion *> volumeRegions;
	map<string, vector<Reference<Primitive> > > instances;
	vector<Reference<Primitive> > *currentInstance;
};


RenderOptions::RenderOptions() {
	// RenderOptions Constructor Implementation
	transformStartTime = 0.f;
	transformEndTime = 1.f;
	FilterName = "box";
	FilmName = "image";
	SamplerName = "lowdiscrepancy";
	AcceleratorName = "bvh";
	RendererName = "sampler";
	SurfIntegratorName = "directlighting";
	VolIntegratorName = "emission";
	CameraName = "perspective";
	currentInstance = NULL;
}


struct GraphicsState { // the current attributes to be used when creating objects, e.g. the scene's lights,
	                   // geometry, participating media, etc. 
	                   // They are mangaged hierarchically by using stacks (vectors) by means of 
	                   // pbrtAttributeBegin() and pbtrAttributeEnd()
	                   // Transforms are treated independently of other attributes. 

	// Graphics State Methods
	GraphicsState();

	// Reference<T> x: is the same as T* x but handles the management of the pointer automatically
	// x.GetPtr() returns the pointer
	// x->member returns x's pointer->member

	Reference<Material> CreateMaterial(const ParamSet &params);

	// Graphics State
	map<string, Reference<Texture<float> > > floatTextures;
	map<string, Reference<Texture<Spectrum> > > spectrumTextures;
	ParamSet materialParams;
	string material;
	map<string, Reference<Material> > namedMaterials; // map is a table with keys
	string currentNamedMaterial;
	ParamSet areaLightParams;
	string areaLight;
	bool reverseOrientation;
};


GraphicsState::GraphicsState() {
	// GraphicsState  default Constructor Implementation
	material = "matte";
	reverseOrientation = false;
}


class TransformCache {
public:
	// TransformCache Public Methods
	void Lookup(const Transform &t, // this is curTransform[0]
		              Transform **tCached, // this is objectToWorld
		              Transform **tCachedInverse) {
			map<Transform, std::pair<Transform *, Transform *> >::iterator iter;

			iter = cache.find(t);

			if (iter == cache.end()) {
				Transform *tr = arena.Alloc<Transform>();

				*tr = t; // curTransform[0]=objectToWorld

				Transform *tinv = arena.Alloc<Transform>();
				*tinv = Transform(Inverse(t));

				cache[t] = std::make_pair(tr, tinv);
				iter = cache.find(t);
				PBRT_ALLOCATED_CACHED_TRANSFORM();
			}
			else
				PBRT_FOUND_CACHED_TRANSFORM();

			if (tCached) *tCached = iter->second.first;
			if (tCachedInverse) *tCachedInverse = iter->second.second;
	}
	void Clear() {
		arena.FreeAll();
		cache.erase(cache.begin(), cache.end());
	}
private:
	// TransformCache Private Data
	map<Transform, std::pair<Transform *, Transform *> > cache;
	MemoryArena arena;
};



// API Static Data
#define STATE_UNINITIALIZED  0
#define STATE_OPTIONS_BLOCK  1
#define STATE_WORLD_BLOCK    2
static int currentApiState = STATE_UNINITIALIZED;
static TransformSet curTransform;
static int activeTransformBits = ALL_TRANSFORMS_BITS;
static map<string, TransformSet> namedCoordinateSystems;
static RenderOptions *renderOptions = NULL;
static GraphicsState graphicsState;
static vector<GraphicsState> pushedGraphicsStates;
static vector<TransformSet> pushedTransforms;
static vector<uint32_t> pushedActiveTransformBits;
static TransformCache transformCache;

// API Macros
#define VERIFY_INITIALIZED(func) \
	if (currentApiState == STATE_UNINITIALIZED) { \
	Error("pbrtInit() must be before calling \"%s()\". " \
	"Ignoring.", func); \
	return; \
	} else /* swallow trailing semicolon */
#define VERIFY_OPTIONS(func) \
	VERIFY_INITIALIZED(func); \
	if (currentApiState == STATE_WORLD_BLOCK) { \
	Error("Options cannot be set inside world block; " \
	"\"%s\" not allowed.  Ignoring.", func); \
	return; \
	} else /* swallow trailing semicolon */
#define VERIFY_WORLD(func) \
	VERIFY_INITIALIZED(func); \
	if (currentApiState == STATE_OPTIONS_BLOCK) { \
	Error("Scene description must be inside world block; " \
	"\"%s\" not allowed. Ignoring.", func); \
	return; \
	} else /* swallow trailing semicolon */
#define FOR_ACTIVE_TRANSFORMS(expr) \
	for (int i = 0; i < MAX_TRANSFORMS; ++i) \
	if (activeTransformBits & (1 << i)) { expr }
#define WARN_IF_ANIMATED_TRANSFORM(func) \
	do { if (curTransform.IsAnimated()) \
	Warning("Animated transformations set; ignoring for \"%s\"" \
	"and using the start transform only", func); \
	} while (false)

// Object Creation Function Definitions
Reference<Shape> MakeShape(const string &name,
						   const Transform *object2world, const Transform *world2object,
						   bool reverseOrientation, const ParamSet &paramSet) {
							   Shape *s = NULL;

							   if (name == "sphere")
								   s = CreateSphereShape(object2world, world2object,
								   reverseOrientation, paramSet);
							   // Create remaining _Shape_ types
							   else if (name == "cylinder")
								   s = CreateCylinderShape(object2world, world2object, reverseOrientation,
								   paramSet);
							   else if (name == "disk")
								   s = CreateDiskShape(object2world, world2object, reverseOrientation,
								   paramSet);
							   else if (name == "cone")
								   s = CreateConeShape(object2world, world2object, reverseOrientation,
								   paramSet);
							   else if (name == "paraboloid")
								   s = CreateParaboloidShape(object2world, world2object, reverseOrientation,
								   paramSet);
							   else if (name == "hyperboloid")
								   s = CreateHyperboloidShape(object2world, world2object, reverseOrientation,
								   paramSet);
							   else if (name == "trianglemesh")
								   s = CreateTriangleMeshShape(object2world, world2object, reverseOrientation,
								   paramSet, &graphicsState.floatTextures);
							   else if (name == "heightfield")
								   s = CreateHeightfieldShape(object2world, world2object, reverseOrientation,
								   paramSet);
							   else if (name == "loopsubdiv")
								   s = CreateLoopSubdivShape(object2world, world2object, reverseOrientation,
								   paramSet);
							   else if (name == "nurbs")
								   s = CreateNURBSShape(object2world, world2object, reverseOrientation,
								   paramSet);
							   else
								   Warning("Shape \"%s\" unknown.", name.c_str());
							   paramSet.ReportUnused();
							   return s;
}


Reference<Material> MakeMaterial(const string &name,
								 const Transform &mtl2world,
								 const TextureParams &mp) {
									 Material *material = NULL;
									 if (name == "matte")
										 material = CreateMatteMaterial(name, mtl2world, mp);
									 else if (name == "plastic")
										 material = CreatePlasticMaterial(name, mtl2world, mp);
									 else if (name == "translucent")
										 material = CreateTranslucentMaterial(name, mtl2world, mp);
									 else if (name == "glass")
										 material = CreateGlassMaterial(name, mtl2world, mp);
									 else if (name == "mirror")
										 material = CreateMirrorMaterial(name, mtl2world, mp);
									 else if (name == "mix") {
										 string m1 = mp.FindString("namedmaterial1", "");
										 string m2 = mp.FindString("namedmaterial2", "");
										 Reference<Material> mat1 = graphicsState.namedMaterials[m1];
										 Reference<Material> mat2 = graphicsState.namedMaterials[m2];
										 if (!mat1) {
											 Error("Named material \"%s\" undefined.  Using \"matte\"",
												 m1.c_str());
											 mat1 = MakeMaterial("matte", curTransform[0], mp);
										 }
										 if (!mat2) {
											 Error("Named material \"%s\" undefined.  Using \"matte\"",
												 m2.c_str());
											 mat2 = MakeMaterial("matte", curTransform[0], mp);
										 }

										 material = CreateMixMaterial(name, mtl2world, mp, mat1, mat2);
									 }
									 else if (name == "metal")
										 material = CreateMetalMaterial(name, mtl2world, mp);
									 else if (name == "substrate")
										 material = CreateSubstrateMaterial(name, mtl2world, mp);
									 else if (name == "uber")
										 material = CreateUberMaterial(name, mtl2world, mp);
									 else if (name == "subsurface")
										 material = CreateSubsurfaceMaterial(name, mtl2world, mp);
									 else if (name == "kdsubsurface")
										 material = CreateKdSubsurfaceMaterial(name, mtl2world, mp);
									 else if (name == "measured")
										 material = CreateMeasuredMaterial(name, mtl2world, mp);
									 else if (name == "shinymetal")
										 material = CreateShinyMetalMaterial(name, mtl2world, mp);
									 else
										 Warning("Material \"%s\" unknown.", name.c_str());
									 mp.ReportUnused();
									 if (!material) Error("Unable to create material \"%s\"", name.c_str());
									 return material;
}


Reference<Texture<float> > MakeFloatTexture(const string &name,  const string &texClass, 
											const Transform &tex2world, const TextureParams &tp) {
												Texture<float> *tex = NULL;
												if ( texClass  == "constant")
													tex = CreateConstantFloatTexture(name, texClass, tex2world, tp);
												else if ( texClass == "scale")
													tex = CreateScaleFloatTexture(name, texClass, tex2world, tp);
												else if ( texClass  == "mix")
													tex = CreateMixFloatTexture(name, texClass, tex2world, tp);
												else if ( texClass  == "bilerp")
													tex = CreateBilerpFloatTexture(name, texClass, tex2world, tp);
												else if ( texClass  == "imagemap")
													tex = CreateImageFloatTexture(name, texClass, tex2world, tp);
												else if ( texClass  == "uv")
													tex = CreateUVFloatTexture(name, texClass, tex2world, tp);
												else if ( texClass  == "checkerboard")
													tex = CreateCheckerboardFloatTexture(name, texClass, tex2world, tp);
												else if ( texClass  == "dots")
													tex = CreateDotsFloatTexture(name, texClass, tex2world, tp);
												else if ( texClass  == "fbm")
													tex = CreateFBmFloatTexture(name, texClass, tex2world, tp);
												else if ( texClass == "wrinkled")
													tex = CreateWrinkledFloatTexture(name, texClass,tex2world, tp);
												else if ( texClass  == "marble")
													tex = CreateMarbleFloatTexture(name, texClass,tex2world, tp);
												else if ( texClass == "windy")
													tex = CreateWindyFloatTexture(name, texClass,tex2world, tp);
												else
													Warning("Float texture \"%s\" unknown.", texClass.c_str());
												tp.ReportUnused();
												return tex;
}

Reference<Texture<Spectrum> > MakeSpectrumTexture(const string &name, const string & texClass, 
												  const Transform &tex2world, const TextureParams &tp) {
													  Texture<Spectrum> *tex = NULL;
													  if ( texClass == "constant")
														  tex = CreateConstantSpectrumTexture(name, texClass,tex2world, tp);
													  else if (texClass == "scale")
														  tex = CreateScaleSpectrumTexture(name, texClass,tex2world, tp);
													  else if (texClass == "mix")
														  tex = CreateMixSpectrumTexture(name, texClass,tex2world, tp);
													  else if (texClass == "bilerp")
														  tex = CreateBilerpSpectrumTexture(name, texClass,tex2world, tp);
													  else if (texClass == "imagemap")
														  tex = CreateImageSpectrumTexture(name, texClass,tex2world, tp);
													  else if (texClass == "uv")
														  tex = CreateUVSpectrumTexture(name, texClass,tex2world, tp);
													  else if (texClass == "checkerboard")
														  tex = CreateCheckerboardSpectrumTexture(name, texClass,tex2world, tp);
													  else if (texClass == "dots")
														  tex = CreateDotsSpectrumTexture(name, texClass,tex2world, tp);
													  else if (texClass == "fbm")
														  tex = CreateFBmSpectrumTexture(name, texClass,tex2world, tp);
													  else if (texClass == "wrinkled")
														  tex = CreateWrinkledSpectrumTexture(name, texClass,tex2world, tp);
													  else if (texClass == "marble")
														  tex = CreateMarbleSpectrumTexture(name, texClass,tex2world, tp);
													  else if (texClass == "windy")
														  tex = CreateWindySpectrumTexture(name, texClass,tex2world, tp);
													  else
														  Warning("Spectrum texture \"%s\" unknown.", texClass.c_str());
													  tp.ReportUnused();
													  return tex;
}


Light *MakeLight(const string &name,
				 const Transform &light2world, const ParamSet &paramSet) {
					 Light *light = NULL;
					 if (name == "point")
						 light = CreatePointLight(name, light2world, paramSet); // light2world == currTransform[0]
					 else if (name == "spot")
						 light = CreateSpotLight(name, light2world, paramSet);
					 else if (name == "goniometric")
						 light = CreateGoniometricLight(name, light2world, paramSet);
					 else if (name == "projection")
						 light = CreateProjectionLight(name, light2world, paramSet);
					 else if (name == "distant")
						 light = CreateDistantLight(name, light2world, paramSet);
					 else if (name == "infinite" || name == "exinfinite")
						 light = CreateInfiniteLight(name, light2world, paramSet);
					 else
						 Warning("Light \"%s\" unknown.", name.c_str());
					 paramSet.ReportUnused();
					 return light;
}


AreaLight *MakeAreaLight(const string &name,
						 const Transform &light2world, const ParamSet &paramSet,
						 const Reference<Shape> &shape) {
							 AreaLight *area = NULL;
							 if (name == "area" || name == "diffuse")
								 area = CreateDiffuseAreaLight(name, light2world, paramSet, shape);
							 else
								 Warning("Area light \"%s\" unknown.", name.c_str());
							 paramSet.ReportUnused();
							 return area;
}


VolumeRegion *MakeVolumeRegion(const string &name,
							   const Transform &volume2world, const ParamSet &paramSet) {
								   VolumeRegion *vr = NULL;
								   if (name == "homogeneous")
									   vr = CreateHomogeneousVolumeDensityRegion(volume2world, paramSet);
								   else if (name == "volumegrid")
									   vr = CreateGridVolumeRegion(volume2world, paramSet);
								   else if (name == "exponential")
									   vr = CreateExponentialVolumeRegion(volume2world, paramSet);
								   else
									   Warning("Volume region \"%s\" unknown.", name.c_str());
								   paramSet.ReportUnused();
								   return vr;
}


SurfaceIntegrator *MakeSurfaceIntegrator(const string &name,
										 const ParamSet &paramSet) {
											 SurfaceIntegrator *si = NULL;
											 if (name == "whitted")
												 si = CreateWhittedSurfaceIntegrator(paramSet);
											 else if (name == "directlighting")
												 si = CreateDirectLightingIntegrator(paramSet);
											 else if (name == "path")
												 si = CreatePathSurfaceIntegrator(paramSet);
											 else if (name == "photonmap" || name == "exphotonmap")
												 si = CreatePhotonMapSurfaceIntegrator(paramSet);
											 else if (name == "irradiancecache")
												 si = CreateIrradianceCacheIntegrator(paramSet);
											 else if (name == "igi")
												 si = CreateIGISurfaceIntegrator(paramSet);
											 else if (name == "dipolesubsurface")
												 si = CreateDipoleSubsurfaceIntegrator(paramSet);
											 else if (name == "ambientocclusion")
												 si = CreateAmbientOcclusionIntegrator(paramSet);
											 else if (name == "useprobes")
												 si = CreateRadianceProbesSurfaceIntegrator(paramSet);
											 else if (name == "diffuseprt")
												 si = CreateDiffusePRTIntegratorSurfaceIntegrator(paramSet);
											 else if (name == "glossyprt")
												 si = CreateGlossyPRTIntegratorSurfaceIntegrator(paramSet);
											 else
												 Warning("Surface integrator \"%s\" unknown.", name.c_str());

											 paramSet.ReportUnused();
											 return si;
}


VolumeIntegrator *MakeVolumeIntegrator(const string &name,
									   const ParamSet &paramSet) {
										   VolumeIntegrator *vi = NULL;
										   if (name == "single")
											   vi = CreateSingleScatteringIntegrator(paramSet);
										   else if (name == "emission")
											   vi = CreateEmissionVolumeIntegrator(paramSet);
										   else
											   Warning("Volume integrator \"%s\" unknown.", name.c_str());
										   paramSet.ReportUnused();
										   return vi;
}


Primitive *MakeAccelerator(const string &name,
						   const vector<Reference<Primitive> > &prims,
						   const ParamSet &paramSet) {
							   Primitive *accel = NULL;
							   if (name == "bvh")
								   accel = CreateBVHAccelerator(prims, paramSet);
							   else if (name == "grid")
								   accel = CreateGridAccelerator(prims, paramSet);
							   else if (name == "kdtree")
								   accel = CreateKdTreeAccelerator(prims, paramSet);
							   else
								   Warning("Accelerator \"%s\" unknown.", name.c_str());
							   paramSet.ReportUnused();
							   return accel;
}


Camera *MakeCamera(const string &name,
				   const ParamSet &paramSet,
				   const TransformSet &cam2worldSet, float transformStart,
				   float transformEnd, Film *film) {
					   Camera *camera = NULL;
					   Assert(MAX_TRANSFORMS == 2);
					   Transform *cam2world[2];
					   transformCache.Lookup(cam2worldSet[0], &cam2world[0], NULL);
					   transformCache.Lookup(cam2worldSet[1], &cam2world[1], NULL);
					   AnimatedTransform animatedCam2World(cam2world[0], transformStart,
						   cam2world[1], transformEnd);
					   if (name == "perspective")
						   camera = CreatePerspectiveCamera(paramSet, animatedCam2World, film);
					   else if (name == "orthographic")
						   camera = CreateOrthographicCamera(paramSet, animatedCam2World, film);
					   else if (name == "environment")
						   camera = CreateEnvironmentCamera(paramSet, animatedCam2World, film);
					   else
						   Warning("Camera \"%s\" unknown.", name.c_str());
					   paramSet.ReportUnused();
					   return camera;
}


Sampler *MakeSampler(const string &name,
					 const ParamSet &paramSet, const Film *film, const Camera *camera) {
						 Sampler *sampler = NULL;
						 if (name == "adaptive")
							 sampler = CreateAdaptiveSampler(paramSet, film, camera);
						 else if (name == "bestcandidate")
							 sampler = CreateBestCandidateSampler(paramSet, film, camera);
						 else if (name == "halton")
							 sampler = CreateHaltonSampler(paramSet, film, camera);
						 else if (name == "lowdiscrepancy")
							 sampler = CreateLowDiscrepancySampler(paramSet, film, camera);
						 else if (name == "random")
							 sampler = CreateRandomSampler(paramSet, film, camera);
						 else if (name == "stratified")
							 sampler = CreateStratifiedSampler(paramSet, film, camera);
						 else
							 Warning("Sampler \"%s\" unknown.", name.c_str());
						 paramSet.ReportUnused();
						 return sampler;
}


Filter *MakeFilter(const string &name,
				   const ParamSet &paramSet) {
					   Filter *filter = NULL;
					   if (name == "box")
						   filter = CreateBoxFilter(paramSet);
					   else if (name == "gaussian")
						   filter = CreateGaussianFilter(paramSet);
					   else if (name == "mitchell")
						   filter = CreateMitchellFilter(paramSet);
					   else if (name == "sinc")
						   filter = CreateSincFilter(paramSet);
					   else if (name == "triangle")
						   filter = CreateTriangleFilter(paramSet);
					   else
						   Warning("Filter \"%s\" unknown.", name.c_str());
					   paramSet.ReportUnused();
					   return filter;
}


Film *MakeFilm(const string &name,
			   const ParamSet &paramSet, Filter *filter) {
				   Film *film = NULL;
				   if (name == "image")
					   film = CreateImageFilm(paramSet, filter);
				   else
					   Warning("Film \"%s\" unknown.", name.c_str());
				   paramSet.ReportUnused();
				   return film;
}



// API Function Definitions
void pbrtInit(const Options &opt) {
	PbrtOptions = opt;
	// API Initialization
	if (currentApiState != STATE_UNINITIALIZED)
		Error("pbrtInit() has already been called.");
	currentApiState = STATE_OPTIONS_BLOCK;
	renderOptions = new RenderOptions;
	graphicsState = GraphicsState();
	SampledSpectrum::Init();
}


void pbrtCleanup() {
	ProbesCleanup();
	// API Cleanup
	if (currentApiState == STATE_UNINITIALIZED)
		Error("pbrtCleanup() called without pbrtInit().");
	else if (currentApiState == STATE_WORLD_BLOCK)
		Error("pbrtCleanup() called while inside world block.");
	currentApiState = STATE_UNINITIALIZED;
	delete renderOptions;
	renderOptions = NULL;
}


void pbrtIdentity() {
	VERIFY_INITIALIZED("Identity");
	FOR_ACTIVE_TRANSFORMS(curTransform[i] = Transform();)
}


void pbrtTranslate(float dx, float dy, float dz) {
	VERIFY_INITIALIZED("Translate");
	FOR_ACTIVE_TRANSFORMS(curTransform[i] =
		curTransform[i] * Translate(Vector(dx, dy, dz));)
}


void pbrtTransform(float tr[16]) {
	VERIFY_INITIALIZED("Transform");
	FOR_ACTIVE_TRANSFORMS(curTransform[i] = Transform(Matrix4x4(
		tr[0], tr[4], tr[8], tr[12],
		tr[1], tr[5], tr[9], tr[13],
		tr[2], tr[6], tr[10], tr[14],
		tr[3], tr[7], tr[11], tr[15]));)
}


void pbrtConcatTransform(float tr[16]) {
	VERIFY_INITIALIZED("ConcatTransform");
	FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * Transform(
		Matrix4x4(tr[0], tr[4], tr[8], tr[12],
		tr[1], tr[5], tr[9], tr[13],
		tr[2], tr[6], tr[10], tr[14],
		tr[3], tr[7], tr[11], tr[15]));)
}


void pbrtRotate(float angle, float dx, float dy, float dz) {
	VERIFY_INITIALIZED("Rotate");
	FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * Rotate(angle, Vector(dx, dy, dz));)
}


void pbrtScale(float sx, float sy, float sz) {
	VERIFY_INITIALIZED("Scale");
	FOR_ACTIVE_TRANSFORMS(curTransform[i] = curTransform[i] * Scale(sx, sy, sz);)
}


void pbrtLookAt(float ex, float ey, float ez, float lx, float ly,
				float lz, float ux, float uy, float uz) {
					VERIFY_INITIALIZED("LookAt");
					FOR_ACTIVE_TRANSFORMS({ Warning("This version of pbrt fixes a bug in the LookAt transformation.\n"
						"If your rendered images unexpectedly change, add a \"Scale -1 1 1\"\n"
						"to the start of your scene file."); break; })
						FOR_ACTIVE_TRANSFORMS(curTransform[i] =
						curTransform[i] * LookAt( pbrt::Point(ex, ey, ez), pbrt::Point(lx, ly, lz), Vector(ux, uy, uz));)
						// LookAT() is the inverse to the cameraToWorld Transformation; So, curTransform in the context of 
						// camera positioning is the worldToCamera transformation.

                 /*
						Cvec3 p(ex,ey,ez);
					Cvec3 q(lx,ly,lz);
					Cvec3 u(ux,uy,uz);
					Cvec3 z=q-p;
					z=normalize(z);
					Cvec3 y=u;
					y=normalize(y);
					Cvec3 x;
					x=cross(y,z);
					*/

					//g_eyeRbt= Matrix4::makeRbt(x,y,z,p);
					//g_eyeRbt = Matrix4::makeTranslation( Cvec3(0, 0.25, 10) ); // ClASS METHOD CALL => AN INSTANCE OF MATRIX4
					
				//	Transform temp_Transfrom=LookAt( Point(ex, ey, ez), Point(lx, ly, lz), Vector(ux, uy, uz));;
				//	Matrix4x4 temp_Matrix=temp_Transfrom.GetMatrix(); GetMetrix() returns the Structure that contains m[4][4]
				//	Matrix4 temp(temp_Matrix.m[0][0] , temp_Matrix.m[0][1],temp_Matrix.m[0][2],temp_Matrix.m[0][3],temp_Matrix.m[1][0],temp_Matrix.m[1][1],temp_Matrix.m[1][2],temp_Matrix.m[1][3],temp_Matrix.m[2][0],temp_Matrix.m[2][1],temp_Matrix.m[2][2],temp_Matrix.m[2][3], temp_Matrix.m[3][0] ,temp_Matrix.m[3][1], temp_Matrix.m[3][2],temp_Matrix.m[3][3]);
				//	g_eyeRbt=temp;

}


void pbrtCoordinateSystem(const string &name) {
	VERIFY_INITIALIZED("CoordinateSystem");
	namedCoordinateSystems[name] = curTransform;
}


void pbrtCoordSysTransform(const string &name) {
	VERIFY_INITIALIZED("CoordSysTransform");
	if (namedCoordinateSystems.find(name) !=
		namedCoordinateSystems.end())
		curTransform = namedCoordinateSystems[name];
	else
		Warning("Couldn't find named coordinate system \"%s\"",
		name.c_str());
}


void pbrtActiveTransformAll() {
	activeTransformBits = ALL_TRANSFORMS_BITS;
}


void pbrtActiveTransformEndTime() {
	activeTransformBits = END_TRANSFORM_BITS;
}


void pbrtActiveTransformStartTime() {
	activeTransformBits = START_TRANSFORM_BITS;
}


void pbrtTransformTimes(float start, float end) {
	VERIFY_OPTIONS("TransformTimes");
	renderOptions->transformStartTime = start;
	renderOptions->transformEndTime = end;
}


void pbrtPixelFilter(const string &name, const ParamSet &params) {
	VERIFY_OPTIONS("PixelFilter");
	renderOptions->FilterName = name;
	renderOptions->FilterParams = params;
}


void pbrtFilm(const string &type, const ParamSet &params) {
	VERIFY_OPTIONS("Film");
	renderOptions->FilmParams = params;
	renderOptions->FilmName = type;
	g_windowWidth=params.FindOneInt("yresolution",512);
	g_windowHeight=params.FindOneInt("xresolution",512);
	//window->resize(g_windowWidth, g_windowHeight);
}


void pbrtSampler(const string &name, const ParamSet &params) {
	VERIFY_OPTIONS("Sampler");
	renderOptions->SamplerName = name;
	renderOptions->SamplerParams = params;
}


void pbrtAccelerator(const string &name, const ParamSet &params) {
	VERIFY_OPTIONS("Accelerator");
	renderOptions->AcceleratorName = name;
	renderOptions->AcceleratorParams = params;
}


void pbrtSurfaceIntegrator(const string &name, const ParamSet &params) {
	VERIFY_OPTIONS("SurfaceIntegrator");
	renderOptions->SurfIntegratorName = name;
	renderOptions->SurfIntegratorParams = params;
}


void pbrtVolumeIntegrator(const string &name, const ParamSet &params) {
	VERIFY_OPTIONS("VolumeIntegrator");
	renderOptions->VolIntegratorName = name;
	renderOptions->VolIntegratorParams = params;
}


void pbrtRenderer(const string &name, const ParamSet &params) {
	VERIFY_OPTIONS("Renderer");
	renderOptions->RendererName = name;
	renderOptions->RendererParams = params;
}


void pbrtCamera(const string &name, const ParamSet &params) {
	VERIFY_OPTIONS("Camera");
	renderOptions->CameraName = name;	
	renderOptions->CameraParams = params;
	renderOptions->CameraToWorld = Inverse(curTransform); // curTransform is the worldToCamera transformation
	namedCoordinateSystems["camera"] = renderOptions->CameraToWorld;

	g_frustMinFov=params.FindOneFloat("fov",60);
	// get frameaspectratio as well
	//window->updateFrustFovY();



}


void pbrtWorldBegin() {
	VERIFY_OPTIONS("WorldBegin");
	currentApiState = STATE_WORLD_BLOCK;
	for (int i = 0; i < MAX_TRANSFORMS; ++i)
		curTransform[i] = Transform();    // reset all the current transforms; camera transforms and world transforms are independent of each other

	activeTransformBits = ALL_TRANSFORMS_BITS;
	namedCoordinateSystems["world"] = curTransform;
	render_flag=0;
}


void pbrtAttributeBegin() {
	VERIFY_WORLD("AttributeBegin");
	pushedGraphicsStates.push_back(graphicsState);
	pushedTransforms.push_back(curTransform);
	pushedActiveTransformBits.push_back(activeTransformBits);
}


void pbrtAttributeEnd() {
	VERIFY_WORLD("AttributeEnd");
	if (!pushedGraphicsStates.size()) {
		Error("Unmatched pbrtAttributeEnd() encountered. "
			"Ignoring it.");
		return;
	}
	graphicsState = pushedGraphicsStates.back();
	pushedGraphicsStates.pop_back();
	curTransform = pushedTransforms.back();
	pushedTransforms.pop_back();
	activeTransformBits = pushedActiveTransformBits.back();
	pushedActiveTransformBits.pop_back();
}


void pbrtTransformBegin() {
	VERIFY_WORLD("TransformBegin");
	pushedTransforms.push_back(curTransform);
	pushedActiveTransformBits.push_back(activeTransformBits);
}


void pbrtTransformEnd() {
	VERIFY_WORLD("TransformEnd");
	if (!pushedTransforms.size()) {
		Error("Unmatched pbrtTransformEnd() encountered. "
			"Ignoring it.");
		return;
	}
	curTransform = pushedTransforms.back();
	pushedTransforms.pop_back();
	activeTransformBits = pushedActiveTransformBits.back();
	pushedActiveTransformBits.pop_back();
}



void pbrtTexture(const string &name, const string &type,
				 const string &texname, const ParamSet &params) { 
					 // Texture ¡°name¡±¡°type¡±¡°class¡± [ parameter-list ]
					 // texname ==  class: "imagemap", "constant", "windy",...
					
					 VERIFY_WORLD("Texture");

					 TextureParams tp(params, params, graphicsState.floatTextures,
						 graphicsState.spectrumTextures);
					 
					 if (type == "float")  { // float texture 
						 // Create _float_ texture and store in _floatTextures_
						 if (graphicsState.floatTextures.find(name) !=
							 graphicsState.floatTextures.end())
							 Info("Texture \"%s\" being redefined", name.c_str());
						 WARN_IF_ANIMATED_TRANSFORM("Texture");

						 Reference<Texture<float> > ft = MakeFloatTexture(name, texname,
							 curTransform[0], tp);

						 if (ft) graphicsState.floatTextures[name] = ft;
					 }

					 else if (type == "color" || type == "spectrum")  { // color texture: "color" and "spectrum" are synonyms here
						 // Create _color_ texture and store in _spectrumTextures_
						 if (graphicsState.spectrumTextures.find(name) != graphicsState.spectrumTextures.end())
							 Info("Texture \"%s\" being redefined", name.c_str());
						 WARN_IF_ANIMATED_TRANSFORM("Texture");
						 
					 		 

						 Reference<Texture<Spectrum> > st = MakeSpectrumTexture(name, texname,
							 curTransform[0], tp);

						 // st: texture object; it has RGBSpectrum *texels, width, height members
						 // for retrieving the texture data; class RGBSpectrum has toRGB( float *rgb) method.

						 if (st) graphicsState.spectrumTextures[name] = st;
					 }

					 else
						 Error("Texture type \"%s\" unknown.", type.c_str());
}



void pbrtMaterial(const string &name, const ParamSet &params) {
	VERIFY_WORLD("Material");

	graphicsState.material = name; // graphicsState.material contains the current material name
	graphicsState.materialParams = params; // graphicsState.materialParams contains the parameters of
	                                       // the current material

	graphicsState.currentNamedMaterial = "";
}


void pbrtMakeNamedMaterial(const string &name,
						   const ParamSet &params) {
							   VERIFY_WORLD("MakeNamedMaterial");
							   // error checking, warning if replace, what to use for transform?
							   TextureParams mp(params, graphicsState.materialParams,
								   graphicsState.floatTextures,
								   graphicsState.spectrumTextures);
							   string matName = mp.FindString("type");
							   WARN_IF_ANIMATED_TRANSFORM("MakeNamedMaterial");
							   if (matName == "") Error("No parameter string \"type\" found in MakeNamedMaterial");
							   else {
								   Reference<Material> mtl = MakeMaterial(matName, curTransform[0], mp);
								   if (mtl) graphicsState.namedMaterials[name] = mtl;
							   }
}



void pbrtNamedMaterial(const string &name) {
	VERIFY_WORLD("NamedMaterial");
	graphicsState.currentNamedMaterial = name;
}


void pbrtLightSource(const string &name, const ParamSet &params) {
	VERIFY_WORLD("LightSource");
	WARN_IF_ANIMATED_TRANSFORM("LightSource");
	Light *lt = MakeLight(name, curTransform[0], params);
	if (lt == NULL)
		Error("pbrtLightSource: light type \"%s\" unknown.", name.c_str());
	else
		renderOptions->lights.push_back(lt);
	/////////////////////////////////////////////////////////////////////////////

	/*
	Vector x;
	Point a(0,0,0);
	Spectrum b;
	float rgbColor[3] = { 1, 1, 1};

	a=params.FindOnePoint("from",a);
	b=params.FindOneSpectrum("L", Spectrum::FromRGB( rgbColor ));

	b.ToRGB( rgbColor );

	g_light1Color = Cvec4( rgbColor[0], rgbColor[1], rgbColor[2], 1.0);

	g_light1Pos[0]=a[0];
	g_light1Pos[1]=a[1];
	g_light1Pos[2]=a[2];

	*/
	/////////////////////////////////////////////////////////////////////////////////

}


void pbrtAreaLightSource(const string &name,
						 const ParamSet &params) {
							 VERIFY_WORLD("AreaLightSource");
							 graphicsState.areaLight = name;
							 graphicsState.areaLightParams = params;
}


void pbrtShape(const string &name, const ParamSet &params) { // params = "geometrical parameters"
	VERIFY_WORLD("Shape");
	Reference<Primitive> prim;

	AreaLight *area = NULL;
	if (!curTransform.IsAnimated()) {

		// Create primitive for static shape
		Transform *obj2world, *world2obj;
		transformCache.Lookup(curTransform[0], &obj2world, &world2obj);

		Reference<Shape> shape = MakeShape(name, obj2world, world2obj,
			graphicsState.reverseOrientation, params);

		//  TriangleMesh has members:
		//  protected:
        //    TriangleMesh Protected Data
        //    int ntris, nverts;
        //    int *vertexIndex;
        //     Point *p;
        //     Normal *n;
        //     Vector *s;
        //     float *uvs;
        //      Reference<Texture<float> > alphaTexture;


		if (!shape) return;
		
		Reference<Material> mtl = graphicsState.CreateMaterial(params);
		// graphicsState.CreateMaterial(params) creates a particular texture but is assigned
		// to the most general type Material.

		 // graphicsState.material, graphicsState.materialParams are the current material 
		 // name and parameters; mtl is the current material
		 // graphicsState contains material parameters, so it CreateMaterial () can use it
		 // to create material.

		// mtl refers to several textures which define the parameters of the material

		params.ReportUnused();

		// Possibly create area light for shape
		if (graphicsState.areaLight != "") {
			area = MakeAreaLight(graphicsState.areaLight, curTransform[0],
				graphicsState.areaLightParams, shape);
		}
		prim = new GeometricPrimitive(shape, mtl, area); // just init the members
		// once the primitive is created, the current material name and parameters are lost.
		// To keep the material name, it should be stored in mtl itself. Expand Material so that
		// it has a member name. 

		// 
	} else {
		// Create primitive for animated shape

		// Create initial _Shape_ for animated shape
		if (graphicsState.areaLight != "")
			Warning("Ignoring currently set area light when creating "
			"animated shape");
		Transform *identity;
		transformCache.Lookup(Transform(), &identity, NULL);
		Reference<Shape> shape = MakeShape(name, identity, identity,
			graphicsState.reverseOrientation, params);
		if (!shape) return;
		Reference<Material> mtl = graphicsState.CreateMaterial(params);
		params.ReportUnused();

		// Get _animatedWorldToObject_ transform for shape
		Assert(MAX_TRANSFORMS == 2);
		Transform *world2obj[2];
		transformCache.Lookup(curTransform[0], NULL, &world2obj[0]);
		transformCache.Lookup(curTransform[1], NULL, &world2obj[1]);
		AnimatedTransform
			animatedWorldToObject(world2obj[0], renderOptions->transformStartTime,
			world2obj[1], renderOptions->transformEndTime);
		Reference<Primitive> baseprim = new GeometricPrimitive(shape, mtl, NULL);
		if (!baseprim->CanIntersect()) {
			// Refine animated shape and create BVH if more than one shape created
			vector<Reference<Primitive> > refinedPrimitives;
			baseprim->FullyRefine(refinedPrimitives);
			if (refinedPrimitives.size() == 0) return;
			if (refinedPrimitives.size() > 1)
				baseprim = new BVHAccel(refinedPrimitives);
			else
				baseprim = refinedPrimitives[0];
		}
		prim = new TransformedPrimitive(baseprim, animatedWorldToObject);
	}
	// Add primitive to scene or current instance
	if (renderOptions->currentInstance) {
		if (area)
			Warning("Area lights not supported with object instancing");

		renderOptions->currentInstance->push_back(prim);
	}

	else {

		renderOptions->primitives.push_back(prim);

		if (area != NULL) {
			renderOptions->lights.push_back(area);
		}
	}

}


Reference<Material> GraphicsState::CreateMaterial(const ParamSet &params) {
	// GraphicsState has public members:
	
	// Graphics State
	// map<string, Reference<Texture<float> > > floatTextures;
	//  map<string, Reference<Texture<Spectrum> > > spectrumTextures;
	//  ParamSet materialParams;
	//  string material;
	//  map<string, Reference<Material> > namedMaterials; // map is a table with keys
	//  string currentNamedMaterial;
	//  ParamSet areaLightParams;
	//  string areaLight;
	//  bool reverseOrientation;

	TextureParams mp(params, materialParams, floatTextures, spectrumTextures); // just initialize members

	Reference<Material> mtl; // the default constructor of Reference is ivoked.
	
	if (currentNamedMaterial != "" &&
		namedMaterials.find(currentNamedMaterial) != namedMaterials.end())
		mtl = namedMaterials[graphicsState.currentNamedMaterial];

	if (!mtl)
		mtl = MakeMaterial(material, curTransform[0], mp);
	if (!mtl)
		mtl = MakeMaterial("matte", curTransform[0], mp);
	if (!mtl)
		Severe("Unable to create \"matte\" material?!");
	return mtl;
}


void pbrtReverseOrientation() {
	VERIFY_WORLD("ReverseOrientation");
	graphicsState.reverseOrientation =
		!graphicsState.reverseOrientation;
}


void pbrtVolume(const string &name, const ParamSet &params) {
	VERIFY_WORLD("Volume");
	WARN_IF_ANIMATED_TRANSFORM("Volume");
	VolumeRegion *vr = MakeVolumeRegion(name, curTransform[0], params);
	if (vr) renderOptions->volumeRegions.push_back(vr);
}


void pbrtObjectBegin(const string &name) {
	VERIFY_WORLD("ObjectBegin");
	pbrtAttributeBegin();
	if (renderOptions->currentInstance)
		Error("ObjectBegin called inside of instance definition");
	renderOptions->instances[name] = vector<Reference<Primitive> >();
	renderOptions->currentInstance = &renderOptions->instances[name];
}


void pbrtObjectEnd() {
	VERIFY_WORLD("ObjectEnd");
	if (!renderOptions->currentInstance)
		Error("ObjectEnd called outside of instance definition");
	renderOptions->currentInstance = NULL;
	pbrtAttributeEnd();
}


void pbrtObjectInstance(const string &name) {
	VERIFY_WORLD("ObjectInstance");
	// Object instance error checking
	if (renderOptions->currentInstance) {
		Error("ObjectInstance can't be called inside instance definition");
		return;
	}
	if (renderOptions->instances.find(name) == renderOptions->instances.end()) {
		Error("Unable to find instance named \"%s\"", name.c_str());
		return;
	}
	vector<Reference<Primitive> > &in = renderOptions->instances[name];
	if (in.size() == 0) return;
	if (in.size() > 1 || !in[0]->CanIntersect()) {
		// Refine instance _Primitive_s and create aggregate
		Reference<Primitive> accel =
			MakeAccelerator(renderOptions->AcceleratorName,
			in, renderOptions->AcceleratorParams);
		if (!accel) accel = MakeAccelerator("bvh", in, ParamSet());
		if (!accel) Severe("Unable to create \"bvh\" accelerator");
		in.erase(in.begin(), in.end());
		in.push_back(accel);
	}
	Assert(MAX_TRANSFORMS == 2);
	Transform *world2instance[2];
	transformCache.Lookup(curTransform[0], NULL, &world2instance[0]);
	transformCache.Lookup(curTransform[1], NULL, &world2instance[1]);
	AnimatedTransform animatedWorldToInstance(world2instance[0],
		renderOptions->transformStartTime,
		world2instance[1], renderOptions->transformEndTime);
	Reference<Primitive> prim =
		new TransformedPrimitive(in[0], animatedWorldToInstance);
	renderOptions->primitives.push_back(prim);
}


void pbrtWorldEnd() {
	
	VERIFY_WORLD("WorldEnd");
		// Ensure there are no pushed graphics states
	while (pushedGraphicsStates.size()) {
			Warning("Missing end to pbrtAttributeBegin()");
			pushedGraphicsStates.pop_back();
			pushedTransforms.pop_back();
		}
	while (pushedTransforms.size()) {
			Warning("Missing end to pbrtTransformBegin()");
			pushedTransforms.pop_back();
		}

		
		
	// renderOptions has all the information about primitives, volumeRegions, Cameras, lights,
	// and other information needed to create the scene, the camera, and the renderer.
    // The default constructor of renderOptions:
	    //transformStartTime = 0.f;
	    //transformEndTime = 1.f;
	   //FilterName = "box";
	   //FilmName = "image";
	   //SamplerName = "lowdiscrepancy";
	   //AcceleratorName = "bvh";
	   //RendererName = "sampler";
	   //SurfIntegratorName = "directlighting";
	  //VolIntegratorName = "emission";
	  //CameraName = "perspective";
	   //currentInstance = NULL;

	// For GPU rendering, the information about primitives, volumeRegions, Cameras, lights is
	// used.

	
	renderOptions->GPURender();


		// Clean up after rendering: This is copied from below
	
	graphicsState = GraphicsState();
	transformCache.Clear();
	currentApiState = STATE_OPTIONS_BLOCK;
	ProbesPrint(stdout);
	
	for (int i = 0; i < MAX_TRANSFORMS; ++i)
		curTransform[i] = Transform();

	activeTransformBits = ALL_TRANSFORMS_BITS;
	namedCoordinateSystems.erase(namedCoordinateSystems.begin(),
	                         	namedCoordinateSystems.end());
	ImageTexture<float, float>::ClearCache();
	ImageTexture<RGBSpectrum, Spectrum>::ClearCache();
	
	return; // For GPU rendering, the following code is ignored.

		
	// Create scene
	
	Scene *scene = renderOptions->MakeScene();
	//  MakeScene() sets the following membes of Scene:
	// lights = lts;
    // aggregate = accel;
    // volumeRegion = vr;
    // Scene Constructor Implementation
    //bound = aggregate->WorldBound();
    //if (volumeRegion) bound = Union(bound, volumeRegion->WorldBound());

	
	
	// create renderer
	Renderer *renderer = renderOptions->MakeRenderer();
	 // renderer has members: Render(), Li(), Transmittance(),
	 //  bool visualizeObjectIds;     Sampler *sampler;     Camera *camera;
     // SurfaceIntegrator *surfaceIntegrator;
    //  VolumeIntegrator *volumeIntegrator;methods.
	 // For GPU rendering, we will not use all information contained in renderer, only part of it.
	 // we will use the default sampler renderer, but this sampler renderer will be replaced
	 // by our own rendering using GPU shader.


	if (scene && renderer) renderer->Render(scene);
	TasksCleanup();
	delete renderer;
	delete scene;

		// Clean up after rendering
	graphicsState = GraphicsState();
	transformCache.Clear();
	currentApiState = STATE_OPTIONS_BLOCK;
	ProbesPrint(stdout);
	
	for (int i = 0; i < MAX_TRANSFORMS; ++i)
		curTransform[i] = Transform();

	activeTransformBits = ALL_TRANSFORMS_BITS;
	namedCoordinateSystems.erase(namedCoordinateSystems.begin(),
	                         	namedCoordinateSystems.end());
	ImageTexture<float, float>::ClearCache();
	ImageTexture<RGBSpectrum, Spectrum>::ClearCache();
	

}




Scene *RenderOptions::MakeScene() {
	// Initialize _volumeRegion_ from volume region(s)
	VolumeRegion *volumeRegion;
	if (volumeRegions.size() == 0)
		volumeRegion = NULL;
	else if (volumeRegions.size() == 1)
		volumeRegion = volumeRegions[0];
	else
		volumeRegion = new AggregateVolume(volumeRegions);

	Primitive *accelerator = MakeAccelerator(AcceleratorName,
		primitives, AcceleratorParams); // primitives: member of RenderOptions

	if (!accelerator)
		accelerator = MakeAccelerator("bvh", primitives, ParamSet());
	if (!accelerator)
		Severe("Unable to create \"bvh\" accelerator.");
	Scene *scene = new Scene(accelerator, lights, volumeRegion);
	// Erase primitives, lights, and volume regions from _RenderOptions_

	primitives.erase(primitives.begin(), primitives.end());
	lights.erase(lights.begin(), lights.end());
	volumeRegions.erase(volumeRegions.begin(), volumeRegions.end());
	return scene;
}


Renderer *RenderOptions::MakeRenderer() const {
	Renderer *renderer = NULL;
	Camera *camera = MakeCamera();
	if (RendererName == "metropolis") {
		renderer = CreateMetropolisRenderer(RendererParams, camera);
		RendererParams.ReportUnused();
		// Warn if no light sources are defined
		if (lights.size() == 0)
			Warning("No light sources defined in scene; "
			"possibly rendering a black image.");
	}
	// Create remaining _Renderer_ types
	else if (RendererName == "createprobes") {
		// Create surface and volume integrators
		SurfaceIntegrator *surfaceIntegrator = MakeSurfaceIntegrator(SurfIntegratorName,
			SurfIntegratorParams);
		if (!surfaceIntegrator) Severe("Unable to create surface integrator.");
		VolumeIntegrator *volumeIntegrator = MakeVolumeIntegrator(VolIntegratorName,
			VolIntegratorParams);
		if (!volumeIntegrator) Severe("Unable to create volume integrator.");
		renderer = CreateRadianceProbesRenderer(camera, surfaceIntegrator, volumeIntegrator, RendererParams);
		RendererParams.ReportUnused();
		// Warn if no light sources are defined
		if (lights.size() == 0)
			Warning("No light sources defined in scene; "
			"possibly rendering a black image.");
	}
	else if (RendererName == "aggregatetest") {
		renderer = CreateAggregateTestRenderer(RendererParams, primitives);
		RendererParams.ReportUnused();
	}
	else if (RendererName == "surfacepbrt::Points") {
		pbrt::Point pCamera = camera->CameraToWorld(camera->shutterOpen, pbrt::Point(0, 0, 0));
		renderer = CreateSurfacePointsRenderer(RendererParams, pCamera, camera->shutterOpen);
		RendererParams.ReportUnused();
	}
	else {
		if (RendererName != "sampler")
			Warning("Renderer type \"%s\" unknown.  Using \"sampler\".",
			RendererName.c_str());
		bool visIds = RendererParams.FindOneBool("visualizeobjectids", false);
		RendererParams.ReportUnused();
		Sampler *sampler = MakeSampler(SamplerName, SamplerParams, camera->film, camera);
		if (!sampler) Severe("Unable to create sampler.");
		// Create surface and volume integrators
		SurfaceIntegrator *surfaceIntegrator = MakeSurfaceIntegrator(SurfIntegratorName,
			SurfIntegratorParams);
		if (!surfaceIntegrator) Severe("Unable to create surface integrator.");
		VolumeIntegrator *volumeIntegrator = MakeVolumeIntegrator(VolIntegratorName,
			VolIntegratorParams);
		if (!volumeIntegrator) Severe("Unable to create volume integrator.");
		renderer = new SamplerRenderer(sampler, camera, surfaceIntegrator,
			volumeIntegrator, visIds);
		// Warn if no light sources are defined
		if (lights.size() == 0)
			Warning("No light sources defined in scene; "
			"possibly rendering a black image.");
	}
	return renderer;
}


Camera *RenderOptions::MakeCamera() const {
	Filter *filter = MakeFilter(FilterName, FilterParams);
	Film *film = MakeFilm(FilmName, FilmParams, filter);
	if (!film) Severe("Unable to create film.");
	Camera *camera = ::MakeCamera(CameraName, CameraParams,
		CameraToWorld, renderOptions->transformStartTime,
		renderOptions->transformEndTime, film);
	if (!camera) Severe("Unable to create camera.");
	return camera;
}


void RenderOptions::GPURender(  ) {
	// RenderOptions members:
	//string CameraName;
	//ParamSet CameraParams;
	//TransformSet CameraToWorld;
	//vector<Light *> lights;

	//vector<Reference<Primitive> > primitives; // equiv to g_objectPtrList 

	//mutable vector<VolumeRegion *> volumeRegions;


	// prepare to draw each primitive, by creating vbo and shader program for each primitive
	// AABB needs to be set as well: predefined for now.  

	
	//scene also has light info: Scene *scene = new Scene(accelerator, lights, volumeRegion);
	

	// TriangleMesh has members:  protected members are redeclared as pubic for GPU renderer
	// protected:
    // int ntris, nverts;
    // int *vertexIndex;
    // Point *p;
    // Normal *n;
    // Vector *s;
    // float *uvs;  // (u,v) parameters of the triangle mesh
    // Reference<Texture<float> > alphaTexture;
	
	// GeometricPrimitive has members:  => private  are redeclared as public for GPU renderer

	// private:
    //  Reference<Shape> shape;
    //  Reference<Material> material;
    //  AreaLight *areaLight;
	// 
    // Primitive Public Data
    // const uint32_t primitiveId;
    //  protected:
    //  static uint32_t nextprimitiveId;

	
	enum TexClass {
       TEX_CONST = 0,
	   TEX_IMAGEMAP = 1

	};

	
	enum LightType {
       LIGHT_POINT = 0,
	   LIGHT_DISTANT = 1

	};


	 
    // Loop over primitives and draw them

    for (uint32_t i = 0; i < primitives.size(); ++i) {


		Reference<Primitive> &prim = primitives[i];
				
		//GeometricPrimitive * triangleMesh = (GeometricPrimitive *)const_cast<Primitive * >( prim.GetPtr() );

		Reference<GeometricPrimitive> geoPrim = (GeometricPrimitive *) ( prim.GetPtr() );

		Reference<TriangleMesh> triMesh = (TriangleMesh *) ( geoPrim->shape.GetPtr() );
		
	    //  const Transform *ObjectToWorld: the current transformation for the current shape

		int objId = geoPrim->primitiveId;

	
	
		// struct VertexPNX : public VertexPN {
        //   VertexPNX(float x, float y, float z,
        //     float nx, float ny, float nz,
         //    float u, float v)
         //  : VertexPN(x, y, z, nx, ny, nz), x(u, v) {}

		

	   int nverts = triMesh->nverts; 
	   int nindices = triMesh->ntris * 3;

       vector<VertexPNX>  vtxTriangleMesh ( nverts);

       vector<unsigned short> idxTriangleMesh ( nindices );

	   // Fill vtxTriangleMesh with vertices of triMesh; it is assumed that trianglemesh has only vertices, normals, tex coords,
	    // not the tangent vector s (which is NULL, in this case)

	   // compute the tangent and binormal vectors on the shape
	   
	   std::vector<Vector> tangents, bitangents;

	   // localP is a pointer to POINT; triMesh->localP simply pass the value of the pointer: This is a call by value

	   computeTangentBasis( nverts, triMesh->localP, triMesh->uvs,  triMesh->n, tangents, bitangents);
	   

	   for ( int i=0; i < nverts; i++ ) {
		   vtxTriangleMesh[i] = VertexPNX( triMesh->localP[i][0],  triMesh->localP[i][1], triMesh->localP[i][2],
			                               triMesh->n[i][0], triMesh->n[i][1],triMesh->n[i][1],
										   triMesh->uvs[ 2* i], triMesh->uvs[ 2* i +1] );
		 //    vtxTriangleMesh[i] = VertexPNTBX( triMesh->localP[i][0],  triMesh->localP[i][1], triMesh->localP[i][2],
		// 	                               triMesh->n[i][0], triMesh->n[i][1],triMesh->n[i][2],		 
		//								    tangents[i][0], tangents[i][1], tangents[i][2],
		//									bitangents[i][0], bitangents[i][1], bitangents[i][2],
		//								    triMesh->uvs[ 2* i], triMesh->uvs[ 2* i +1] );

	   }


	   
	   for ( int i=0; i< nindices; i++ ) {
		   idxTriangleMesh[i] = triMesh->vertexIndex[i];
		   
	   }

	   // create a geometry 
	    char meshName[100];
		sprintf(meshName, "triMesh%d", objId );	// objId = geoPrim->primitiveId;
	   shared_ptr<Geometry> triangleMesh ( new SimpleIndexedGeometryPNX(meshName, &vtxTriangleMesh[0], 
		                                            &idxTriangleMesh[0], nverts, nindices, GL_TRIANGLES ) );
		
	//shared_ptr<Geometry> triangleMesh  ( new SimpleIndexedGeometryPNTBX( &vtxTriangleMesh[0], 
	//	                                            &idxTriangleMesh[0], nverts, nindices ) );
		

		// create the object reference frame of triangleMesh

		
		 
		Matrix4x4 mat =  triMesh->ObjectToWorld->GetMatrix(); 

		shared_ptr<Matrix4> TriangleMeshRbt ( new SgRbtNode( &mat.m[0][0], 16 ) ); // 
				

		// the kind of material name:  "matte", "plastic", "translucent", "glass", "mirror", "mix", 
		// "metal", "substrate", "uber", "subsurface", "kdsubsurface", "measured",  "shinymetal"
		//  


		
	
		
		Reference<Material> material = (Material *) ( geoPrim->material.GetPtr() );
		
		//NOTE: typedef RGBSpectrum Spectrum is active in pbrt.h rather than 
        // typedef SampledSpectrum Spectrum;

		if ( material->name == "matte" ) {
			  // MatteMaterial Private Data
              // Reference<Texture<Spectrum> > Kd;
              //  Reference<Texture<float> > sigma, bumpMap;
						
			// Each time material is introduced, a new shader program should be created
			MaterialShader matteMatShader ("matteMat", "shaders/normal-gl3.vshader", "shaders/matte-rainbow-gl3.fshader");
			//MaterialShader sets its member  programDesc_; Another member is uniforms_, which will be set by put() methods 

	       g_matteMat.reset(new MaterialShader( matteMatShader) ); // g_uberMat is a shared_ptr to the newly created object.
		                                                                 // the object created by new() is maintained until detroyed explictly

	        vector< shared_ptr<ShaderTexture> > imgTextures;
	        vector< Cvec3f > constTextures;

			int textureUnit = 0;
			int constTextureUnit = 0;
			

			g_currentMaterial = g_matteMat;

		    Reference<MatteMaterial> matteMat  = (MatteMaterial *) ( geoPrim->material.GetPtr() );
			string diffuseTexClass = matteMat->Kd->texClass;
			string sigmaTexClass  = matteMat->sigma->texClass;

			


			// texture types: constSpectrumTexture constFloatTexture or ImageSpectrumTexure ImageFloatTexture
			// diffuse texture is a spectrum texture, whereas sigma texture and bumpMap texture are float texture

			// check diffuse texture
			if ( diffuseTexClass == "constant" ) {
				
			
			    Reference<ConstantTexture<Spectrum>> constTexture = ( ConstantTexture<Spectrum> *) matteMat->Kd.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				RGBSpectrum value =constTexture->value;

				float rgbColor[3]; 
				value.ToRGB(rgbColor);

				// uMatteMaterial is a struct defined in the shader
				g_currentMaterial->getUniforms().put("uMattMaterial.KdTextureClass", TEX_CONST );
				constTextures.push_back( Cvec3f( rgbColor, 3) );
				g_currentMaterial->getUniforms().put("uMatteMaterial.KdTextureUnit", constTextureUnit++ );

				//g_currentMaterial->getUniforms().put("uMattMaterial.KdConstColor",  Cvec3f( rgbColor, 3) );
				
			    
				
			}

			else if ( diffuseTexClass =="imagemap" ) {

				

				Reference<ImageTexture<RGBSpectrum, Spectrum> > imageTexture = ( ImageTexture<RGBSpectrum, Spectrum> *) matteMat->Kd.GetPtr(); // Kd: Reference to Texture
				
				string filename = imageTexture->filename;
				
				//Reference<Texture<Spectrum> > Kd; 
			    // Texture class has m_texels ( RGBSpectrum *), an array of colors.

				// Add the texture to the shader either through the texture  filename or the texture data itself.
				// We use the latter here.

				// imageTexture->mipmap->pyradmid[0] contains the original image (base level) of texture

			    //RGBSpectrum * texels = imageTexture->mipmap->pyramid[0]->data;  //

				//int width = imageTexture->mipmap->width;
				//int height = imageTexture->mipmap->height;

		        //vector<float> pixData ( width * height * 3 );

				//for (int i=0; i < width * height; i++ ) {
                 //  texels[i].ToRGB( &pixData[3*i]);

				//}
				
				// the sensor in the camera is a linear responsor.
				// When you see the color recorded in the camera on the CRT, which is a non-linear device, it looks much darker than the original
				// To solve this problem, you gamma-encode (gamma correction) the original color, so that the reproduced color by TV looks behave in a linear manner.
				// => This means that the internal RGB color space (RGB 0..255 values) is a gamma encoded color space, not a linear color space.
				// Any digital camera that produces JPEG's will produce RGB pixel values that are gamma encoded (not linear). 
				// So viewing on any PC display will result in the original linear color information. 
				//  if we really had linear response between our RGB colors and luminance, we would need more than 8 bits per channel.

				// this non-linear encoding is not good 
				// if we're going to perform 3-d calculations to create computer-generated images of light sources.

				// In order to correctly create lighting effects, we need to:
                // 1.Convert from sRGB to linear color.
                // 2.Do the lighting accumulation in linear color space.
                //  3.Convert back to sRGB because that's the format the framebuffer needs.


				/* 
				The OpenGL extensions GL_EXT_texure_sRGB and GL_ARB_framebuffer_sRGB basically do steps 1 and 3 for you; 
				when you set a texture's internal type to sRGB, the GPU converts from sRGB to linear space during texel fetch. 
				When framebuffer_sRGB is enabled, the GPU converts from linear back to sRGB before writing your fragment out
				to the framebuffer. Thus your shader runs in linear space (which is fine because it has floating point precision) 
				while your imgTextures and framebuffer are sRGB like they've always been.*

                 The advantage of using these extensions on DirectX 10 hardware is that the conversion happens before texture filtering 
				 and after framebuffer blending - two operations you couldn't "fix" manually in your shader. 
				 So you get linear blending too, which makes the blend of colors look correct.

				*/

				g_currentMaterial->getUniforms().put("uMattMaterial.KdTextureClass", TEX_IMAGEMAP );
				
				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RGB_RGB>( new ShaderImageTexture2D_RGB_RGB( filename.c_str(), true  ) ) );

				//g_currentMaterial->getUniforms().put("uKdTextureUnit", shared_ptr<ShaderImageTexture2Df>(
				//	                                                              new ShaderImageTexture2Df( &pixData[0], width, height ) ) );
				g_currentMaterial->getUniforms().put("uMatteMaterial.KdTextureUnit", textureUnit++ );

				// put() =>  valueMap[name].reset(new TexturesValue(&value, 1) );

				//currentMaterialShader->getUniforms().put("uTexUnit0", shared_ptr<ShaderImageTexture2D>(
				//	                                          new ShaderImageTexture2D("image/Fieldstone.ppm", true)));



			}

			else { 

				Error( "unsupported texture class: %d \n", diffuseTexClass);

			}

			// check sigma texture
			if ( sigmaTexClass == "constant" ) {
				
			
			    Reference<ConstantTexture<float> > constTexture = ( ConstantTexture<float> *) matteMat->sigma.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				float  value =constTexture->value;

				g_currentMaterial->getUniforms().put("uMattMaterial.sigmaTextureClass", TEX_CONST );
				constTextures.push_back( Cvec3f( value, 0,0)  );
				g_currentMaterial->getUniforms().put("uMatteMaterial.sigmaTextureUnit", constTextureUnit++ );
				//g_currentMaterial->getUniforms().put("uMattMaterial.sigmaConstValue", value );
				
			    
				
			}

			else if ( sigmaTexClass =="imagemap" ) {  // imageFloat texture
				

				Reference<ImageTexture<float, float> > imageTexture = ( ImageTexture<float, float> *) matteMat->sigma.GetPtr(); // Kd: Reference to Texture
				string filename = imageTexture->filename;

						       
				g_currentMaterial->getUniforms().put("uMattMaterial.sigmaTextureClass", TEX_IMAGEMAP );
				//textures.push_back( shared_ptr<ShaderImageTexture2D_RGB_RGB>( new ShaderImageTexture2D_RGB_RGB( texels, width, height,false ) ) );

				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RF_RF>( new ShaderImageTexture2D_RF_RF( filename.c_str() ) ) );

			
				g_currentMaterial->getUniforms().put("uMatteMaterial.sigmaTextureUnit", textureUnit++ );
				
				// put() =>  valueMap[name].reset(new TexturesValue(&value, 1) );

				//currentMaterialShader->getUniforms().put("uTexUnit0", shared_ptr<ShaderImageTexture2D>(
				//	                                          new ShaderImageTexture2D("image/Fieldstone.ppm", true)));



			}

			else { 

				Error( "unsupported texture class: %d \n", sigmaTexClass);

			}

			if ( matteMat->bumpMap == NULL ) { // there is no bumpMapTexture => set bumpMapTextureClass uniform to -1. 
				
				g_currentMaterial->getUniforms().put("uMatteMaterial.bumpMapTextureClass", -1 );
			
				g_currentMaterial->getUniforms().put("uMatteMaterial.bumpMapTextureUnit", -1 );
			}

			else { 

			// check bumpMap texture

			
			   string bumpMapTexClass  = matteMat->bumpMap->texClass;  // material may NOT have bumpMap parameter; In that case, 
			                                                        // bumpMap is set to NULL. 
			
			  if ( bumpMapTexClass == "constant" ) {
				
			
			    Reference<ConstantTexture<float>> constTexture = ( ConstantTexture<float> *) matteMat->bumpMap.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				float value =constTexture->value;

				
				g_currentMaterial->getUniforms().put("uMattMaterial.bumpMapTextureClass", TEX_CONST );
				constTextures.push_back( Cvec3f(value, 0,0) );
				g_currentMaterial->getUniforms().put("uMattMaterial.bumpMapTextureUnit", constTextureUnit++  );
				
			    
				
			  }

			  else if ( bumpMapTexClass =="imagemap" ) {
				

				Reference<ImageTexture<float, float> > imageTexture = ( ImageTexture<float, float> *) matteMat->bumpMap.GetPtr(); // Kd: Reference to Texture
				string filename = imageTexture->filename;

				g_currentMaterial->getUniforms().put("uMattMaterial.bumpMapTextureClass", TEX_IMAGEMAP );
				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RF_RF>( new ShaderImageTexture2D_RF_RF( 
					                    filename.c_str()  ) ) );
				g_currentMaterial->getUniforms().put("uMatteMaterial.bumpMapTextureUnit", textureUnit++ );

				// put() =>  valueMap[name].reset(new imgTexturesValue(&value, 1) );

				//currentMaterialShader->getUniforms().put("uTexUnit0", shared_ptr<ShaderImageTexture2D>(
				//	                                          new ShaderImageTexture2D("image/Fieldstone.ppm", true)));



			   }

			else { 

				Error( "unsupported texture class: %d \n", bumpMapTexClass);

			}

			} // matteMat->bumpMap is not NULL

			// send the array of textures to the shader
			
			if ( imgTextures.size() != 0 ) {
			   g_currentMaterial->getUniforms().put( "uTextures",  &imgTextures[0], imgTextures.size() );
			}
			if ( constTextures.size() != 0 ) { 
			   g_currentMaterial->getUniforms().put( "uConstTexture",  &constTextures[0], constTextures.size() );
			}


		} // materail->name == "matte"


		else if ( material->name == "uber" ) {
			
			// Each time material is introduced, a new shader program should be created
			MaterialShader uberMatShader ("uberMat", "shaders/normal-gl3.vshader", "shaders/uber-rainbow-gl3.fshader");
			//MaterialShader sets its member  programDesc_; Another member is uniforms_, which will be set by put() methods 

	       g_uberMat.reset(new MaterialShader( uberMatShader) ); // g_uberMat is a shared_ptr to the newly created object.
		                                                                 // the object created by new() is maintained until detroyed explictly
		   // MaterialShader(uberMat) calls the  copy constructor


        // UberMaterial public  Data:
        // Reference<Texture<Spectrum> > Kd, Ks, Kr, Kt, opacity; Ks= glossy reflection coeff, Kr = specular reflection coeff
        // Reference<Texture<float> > roughness, eta, bumpMap; eta = local index of refraction

			//MaterialShader currentShader ("uber", "shaders/basic-gl3.vshader", "shaders/uber-rainbow-gl3.fshader"); 
			//currentMaterialShader.reset( new MaterialShader( currentShader) );
			vector< shared_ptr<ShaderTexture> > imgTextures;
	        vector< Cvec3f > constTextures;  // The constant textures of roughness, eta, bumpMap  have  one float value, but it should be
			                                 // stored in Cvec3 variable; it will be represented as (roughness, 0,0)

			int constTextureUnit =0;
			int textureUnit = 0;

			g_currentMaterial = g_uberMat; // assignment operator of shared_ptr

			Reference<UberMaterial> uberMat  = (UberMaterial *) ( geoPrim->material.GetPtr() );

			// uberMat->Kd, ->Ks, ->Kr, ->Kt is either an image texture or a constant texture
			// image texture has m_texels member which points to an array of RGBSpectrum or
			//  a value member which is of RGBSpectrum. ( In this implementation, Spectrum = RGBSpectrum, see pbrt.h

			string diffuseTexClass = uberMat->Kd->texClass;  // texClass may be "imagemap", "constant", etc.
				
			string glossyTexClass = uberMat ->Ks->texClass;  // "color Ks" [a,b,c]: which is a constant texture, whose value is stored in value member
			                                                 //  ConstantTexture. In this case, the value of Ks is returned by FindOneSpectrum("Ks", default)

			string reflectTexClass = uberMat->Kr->texClass; // "color Kr" [a,b,c]
			string transTexClass = uberMat->Kt->texClass;  // "color opacity" [a,b,c]

			// currentMaterialShader->getUniforms().put( "Kd", Texture ); => but a texture unit should be passed as the 
			// uniform variable, after  the texture data is loaded and bound to the texture unit
			// 
			string roughnessTexClass = uberMat->roughness->texClass;  // "color opacity" [a,b,c]
			string opacityTexClass = uberMat->opacity->texClass;  // "color opacity" [a,b,c]
			string etaTexClass = uberMat->eta->texClass;  // "color opacity" [a,b,c]
			//string bumpMapTexClass = uberMat->bumpMap->texClass;  // "color opacity" [a,b,c]
			
			// check diffuse texture

			if ( diffuseTexClass == "constant" ) {
				
				// add the material parameters Kd, sigma, bumpMap as uniform variables to the shader

			    Reference<ConstantTexture<Spectrum>> constTexture = ( ConstantTexture<Spectrum> *) uberMat->Kd.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				RGBSpectrum value =constTexture->value;

				float rgbColor[3]; 
				value.ToRGB(rgbColor);

				g_currentMaterial->getUniforms().put("uUberMaterial.KdTextureClass", TEX_CONST );
				constTextures.push_back( Cvec3f( rgbColor, 3) );
				g_currentMaterial->getUniforms().put("uUberMaterial.KdTextureUnit", constTextureUnit++  );
				
			    
			}

			else if ( diffuseTexClass =="imagemap" ) {
				

				Reference<ImageTexture<RGBSpectrum, Spectrum> > imageTexture = ( ImageTexture<RGBSpectrum, Spectrum> *) uberMat->Kd.GetPtr(); // Kd: Reference to Texture
				
				string filename = imageTexture->filename;
								

				g_currentMaterial->getUniforms().put("uUberMaterial.KdTextureClass", TEX_IMAGEMAP );

			    //textures.push_back( shared_ptr<ShaderImageTexture2D_RGBF_RGB>( new ShaderImageTexture2D_RGBF_RGB( &pixData[0], width, height, false ) ) );
			   
				 // for debugging

				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RGB_RGB>(
				  	new ShaderImageTexture2D_RGB_RGB( filename.c_str(), true ) ) );

				//imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RGB_RGB>(
				 // 	new ShaderImageTexture2D_RGB_RGB("sourceimages/sourceimages/Fieldstone.ppm", true) ) );

				// textureUnit is the index to the array of imgTextures which are used to define the same material

				g_currentMaterial->getUniforms().put("uUberMaterial.KdTextureUnit", textureUnit++ );
			}

			else { 

				Error( "unsupported texture class: %d \n", diffuseTexClass);

			}

			// glossy texture
			
			if ( glossyTexClass == "constant" ) {
				
				// add the material parameters Kd, sigma, bumpMap as uniform variables to the shader

			    Reference<ConstantTexture<Spectrum>> constTexture = ( ConstantTexture<Spectrum> *) uberMat->Ks.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				RGBSpectrum value =constTexture->value;

				float rgbColor[3]; 
				value.ToRGB(rgbColor);

				g_currentMaterial->getUniforms().put("uUberMaterial.KsTextureClass", TEX_CONST );
				constTextures.push_back(   Cvec3f( rgbColor, 3) );
				g_currentMaterial->getUniforms().put("uUberMaterial.KsTextureUnit",  constTextureUnit++ );
				
			   
			}

			else if ( glossyTexClass =="imagemap" ) {
				

				Reference<ImageTexture<RGBSpectrum, Spectrum> > imageTexture = ( ImageTexture<RGBSpectrum, Spectrum> *) uberMat->Ks.GetPtr(); // Kd: Reference to Texture
				
				string filename = imageTexture->filename;

				
				g_currentMaterial->getUniforms().put("uUberMaterial.KsTextureClass", TEX_IMAGEMAP );
				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RGB_RGB>(
					            new ShaderImageTexture2D_RGB_RGB( filename.c_str(), true ) ) );
				g_currentMaterial->getUniforms().put("uUberMaterial.KsTextureUnit", textureUnit++ );


			}

			else { 

				Error( "unsupported texture class: %d \n", glossyTexClass);

			}

	 
			// reflect texture  
	  		
			if ( reflectTexClass == "constant" ) {
				
				// add the material parameters Kd, sigma, bumpMap as uniform variables to the shader

			    Reference<ConstantTexture<Spectrum>> constTexture = ( ConstantTexture<Spectrum> *) uberMat->Kr.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				RGBSpectrum value =constTexture->value;

				float rgbColor[3]; 
				value.ToRGB(rgbColor);

				g_currentMaterial->getUniforms().put("uUberMaterial.KrTextureClass", TEX_CONST );
				constTextures.push_back(  Cvec3f( rgbColor, 3) );
				g_currentMaterial->getUniforms().put("uUberMaterial.KrTextureUnit",  constTextureUnit++ );
				
			
			}

			else if ( reflectTexClass =="imagemap" ) {
				

				Reference<ImageTexture<RGBSpectrum, Spectrum> > imageTexture = ( ImageTexture<RGBSpectrum, Spectrum> *) uberMat->Kr.GetPtr(); // Kd: Reference to Texture
				
				string filename = imageTexture->filename;

								
				g_currentMaterial->getUniforms().put("uUberMaterial.KrTextureClass", TEX_IMAGEMAP );
				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RGB_RGB>(
					                     new ShaderImageTexture2D_RGB_RGB( filename.c_str(), true  ) ) );
				g_currentMaterial->getUniforms().put("uUberMaterial.KrTextureUnit", textureUnit++ );

			}

			else { 

				Error( "unsupported texture class: %d \n", reflectTexClass);

			}

	// trans texture  
	  		
			if ( transTexClass == "constant" ) {
				
				// add the material parameters Kd, sigma, bumpMap as uniform variables to the shader

			    Reference<ConstantTexture<Spectrum>> constTexture = ( ConstantTexture<Spectrum> *) uberMat->Kt.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				RGBSpectrum value =constTexture->value;

				float rgbColor[3]; 
				value.ToRGB(rgbColor);

				g_currentMaterial->getUniforms().put("uUberMaterial.KtTextureClass", TEX_CONST );
				constTextures.push_back( Cvec3f( rgbColor, 3) );
				g_currentMaterial->getUniforms().put("uUberMaterial.KtTextureUnit",  constTextureUnit++  );
				
			    
			}

			else if ( transTexClass =="imagemap" ) {
				

				Reference<ImageTexture<RGBSpectrum, Spectrum> > imageTexture = ( ImageTexture<RGBSpectrum, Spectrum> *) uberMat->Kt.GetPtr(); // Kd: Reference to Texture
				string filename = imageTexture->filename;
				
				
				g_currentMaterial->getUniforms().put("uUberMaterial.KtTextureClass", TEX_IMAGEMAP );
				g_currentMaterial->getUniforms().put("uKtTextureUnit", shared_ptr<ShaderImageTexture2D_RGB_RGB>(
					                            new ShaderImageTexture2D_RGB_RGB( filename.c_str(), true  ) ) );
				g_currentMaterial->getUniforms().put("uUberMaterial.KtTextureUnit", textureUnit++ );

			}

			else { 

				Error( "unsupported texture class: %d \n", transTexClass);

			}
	

	// check roughness  texture
			if ( roughnessTexClass == "constant" ) {
				
			
			    Reference<ConstantTexture<float> > constTexture = ( ConstantTexture<float> *) uberMat->roughness.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				float  value =constTexture->value;

				g_currentMaterial->getUniforms().put("uUberMaterial.roughnessTextureClass", TEX_CONST );
				constTextures.push_back( Cvec3f(value, 0,0)  );
				g_currentMaterial->getUniforms().put("uUberMaterial.roughnessTextureUnit", constTextureUnit++  );
				
			    
				
			}

			else if ( roughnessTexClass =="imagemap" ) {  // imageFloat texture
				

				Reference<ImageTexture<float, float> > imageTexture = ( ImageTexture<float, float> *) uberMat->roughness.GetPtr(); // Kd: Reference to Texture
				
				string filename = imageTexture->filename;
				
				g_currentMaterial->getUniforms().put("uUberMaterial.roughnessTextureClass", TEX_IMAGEMAP );

				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RGB_RGB>( 
					     new ShaderImageTexture2D_RGB_RGB( filename.c_str(), true ) ) );

				g_currentMaterial->getUniforms().put("uUberMaterial.roughnessTextureUnit", textureUnit++ );
				// put() =>  valueMap[name].reset(new TexturesValue(&value, 1) );

				//currentMaterialShader->getUniforms().put("uTexUnit0", shared_ptr<ShaderImageTexture2D>(
				//	                                          new ShaderImageTexture2D("image/Fieldstone.ppm", true)));



			}

			else { 

				Error( "unsupported texture class: %d \n", roughnessTexClass);

			}

		// check opaticy  texture
			if ( opacityTexClass == "constant" ) {
				
			
			    Reference<ConstantTexture<float> > constTexture = ( ConstantTexture<float> *) uberMat->opacity.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				RGBSpectrum value =constTexture->value;

				float rgbColor[3]; 
				value.ToRGB(rgbColor);
				

				g_currentMaterial->getUniforms().put("uUberMaterial.opacityTextureClass", TEX_CONST );
				constTextures.push_back( Cvec3f(rgbColor,3) );

				g_currentMaterial->getUniforms().put("uUberMaterial.opacityTextureUnit", constTextureUnit++  );
				
			    
				
			}

			else if ( opacityTexClass =="imagemap" ) {  // imageFloat texture
				


				Reference<ImageTexture<RGBSpectrum, Spectrum> > imageTexture = ( ImageTexture<RGBSpectrum, Spectrum> *) uberMat->opacity.GetPtr(); // Kd: Reference to Texture
				string filename = imageTexture->filename;

				
				g_currentMaterial->getUniforms().put("uUberMaterial.opacityTextureClass", TEX_IMAGEMAP );
				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RGB_RGB>( 
					         new ShaderImageTexture2D_RGB_RGB( filename.c_str(), true  ) ) );

				g_currentMaterial->getUniforms().put("uUberMaterial.opacityTextureUnit", textureUnit++ );

				// put() =>  valueMap[name].reset(new TexturesValue(&value, 1) );

				//currentMaterialShader->getUniforms().put("uTexUnit0", shared_ptr<ShaderImageTexture2D>(
				//	                                          new ShaderImageTexture2D("image/Fieldstone.ppm", true)));



			}

			else { 

				Error( "unsupported texture class: %d \n", roughnessTexClass);

			}

			// check eta   texture
			if ( etaTexClass == "constant" ) {
				
			
			    Reference<ConstantTexture<float> > constTexture = ( ConstantTexture<float> *) uberMat->eta.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				float  value =constTexture->value;

				g_currentMaterial->getUniforms().put("uUberMaterial.etaTextureClass", TEX_CONST );
				constTextures.push_back( Cvec3f(value, 0,0) );
				g_currentMaterial->getUniforms().put("uUberMaterial.etaTextureUnit", constTextureUnit++  );
				
			    
				
			}

			else if ( etaTexClass =="imagemap" ) {  // imageFloat texture
				

				Reference<ImageTexture<float, float> > imageTexture = ( ImageTexture<float, float> *) uberMat->eta.GetPtr(); // Kd: Reference to Texture
				string filename = imageTexture->filename;
				
				g_currentMaterial->getUniforms().put("uUberMaterial.etaTextureClass", TEX_IMAGEMAP );
				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RF_RF>( 
					    new ShaderImageTexture2D_RF_RF( filename.c_str()  ) ) );
				g_currentMaterial->getUniforms().put("uUberMaterial.etaTextureUnit", textureUnit++ );
				// put() =>  valueMap[name].reset(new TexturesValue(&value, 1) );

				//currentMaterialShader->getUniforms().put("uTexUnit0", shared_ptr<ShaderImageTexture2D>(
				//	                                          new ShaderImageTexture2D("image/Fieldstone.ppm", true)));



			}

			else { 

				Error( "unsupported texture class: %d \n", roughnessTexClass);

			}

			// check bumpMap texture
			if ( uberMat->bumpMap == NULL ) { // there is no bumpMapTexture => set bumpMapTextureClass uniform to -1. 
				
				g_currentMaterial->getUniforms().put("uUberMaterial.bumpMapTextureClass", -1 );
			
				g_currentMaterial->getUniforms().put("uUberMaterial.bumpMapTextureUnit", -1 );
			}

			else {

			  string bumpMapTexClass  = uberMat->bumpMap->texClass;  // material may NOT have bumpMap parameter; In that case, 
			                                                        // bumpMap is set to NULL. 
			
			  if ( bumpMapTexClass == "constant" ) {
				
			
			    Reference<ConstantTexture<float>> constTexture = ( ConstantTexture<float> *) uberMat->bumpMap.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				float value =constTexture->value;

				
				g_currentMaterial->getUniforms().put("uUberMaterial.bumpMapTextureClass", TEX_CONST );
				constTextures.push_back( Cvec3f(value, 0,0) );
				g_currentMaterial->getUniforms().put("uUberMaterial.bumpMapTextureUnit", textureUnit++  );
			    
				
			  }

			  else if ( bumpMapTexClass =="imagemap" ) {
				

				Reference<ImageTexture<RGBSpectrum, Spectrum> > imageTexture = ( ImageTexture<RGBSpectrum, Spectrum> *) uberMat->bumpMap.GetPtr(); // Kd: Reference to Texture
				string filename = imageTexture->filename;

						       
				g_currentMaterial->getUniforms().put("uUberMaterial.bumpMapTextureClass", TEX_IMAGEMAP );
				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RF_RF>( new ShaderImageTexture2D_RF_RF( filename.c_str() ) ) );
				g_currentMaterial->getUniforms().put("uUberMaterial.bumpMapTextureUnit", textureUnit++ );
				// put() =>  valueMap[name].reset(new TexturesValue(&value, 1) );

				//currentMaterialShader->getUniforms().put("uTexUnit0", shared_ptr<ShaderImageTexture2D>(
				//	                                          new ShaderImageTexture2D("image/Fieldstone.ppm", true)));



  			  }

			  else { 

				Error( "unsupported texture class: %d \n", bumpMapTexClass);

			  }
			} // uberMat->bumpMap is not NULL

			
			// send the array of textures to the shader
			
			if ( imgTextures.size() != 0 ) {
			    g_currentMaterial->getUniforms().put( "uTextures",  &imgTextures[0], imgTextures.size() );
				  // put() simply stores the name and its value at valueMap, so that it will be used when drawing

				// for debugging, moon jung, 2014/4/22
			  // g_currentMaterial->getUniforms().put( "uTexture",  textures[0] );


			//	g_currentMaterial->getUniforms().put("uTexture", shared_ptr<ShaderImageTexture2D_RGB_RGB>(
			//	  	new ShaderImageTexture2D_RGB_RGB("sourceimages/sourceimages/Fieldstone.ppm", true) ) );
			}

			if (constTextures.size() != 0 ) {
			   g_currentMaterial->getUniforms().put( "uConstTextures",  &constTextures[0], constTextures.size() );
			}



		} // material->name =="uber"

		else if ( material->name == "plastic" ) {
			
			// Each time material is introduced, a new shader program should be created
			MaterialShader plasticMatShader ("plasticMat", "shaders/normal-gl3.vshader", "shaders/plastic-rainbow-gl3.fshader");
			//MaterialShader sets its member  programDesc_; Another member is uniforms_, which will be set by put() methods 

	       g_plasticMat.reset(new MaterialShader( plasticMatShader) ); // g_uberMat is a shared_ptr to the newly created object.
		                                                                 // the object created by new() is maintained until detroyed explictly
		   // MaterialShader(uberMat) calls the copy constructor


        // PlasticMaterial public  Data:
        // Reference<Texture<Spectrum> > Kd, Ks, 
        // Reference<Texture<float> > roughness,  bumpMap; 

			//MaterialShader currentShader ("uber", "shaders/basic-gl3.vshader", "shaders/uber-rainbow-gl3.fshader"); 
			//currentMaterialShader.reset( new MaterialShader( currentShader) );

			vector< shared_ptr<ShaderTexture> > imgTextures;

	        vector< Cvec3f > constTextures;  // The constant textures of roughness, eta, bumpMap  have  one float value, but it should be
			                                 // stored in Cvec3 variable; it will be represented as (roughness, 0,0)

			int constTextureUnit =0;
			int textureUnit = 0;

			g_currentMaterial = g_plasticMat;

			Reference<PlasticMaterial> plasticMat  = (PlasticMaterial *) ( geoPrim->material.GetPtr() );

			// uberMat->Kd, ->Ks, ->Kr, ->Kt is either an image texture or a constant texture
			// image texture has m_texels member which points to an array of RGBSpectrum or
			//  a value member which is of RGBSpectrum. ( In this implementation, Spectrum = RGBSpectrum, see pbrt.h

			string diffuseTexClass = plasticMat->Kd->texClass;  // texClass may be "imagemap", "constant", etc.
				
			string glossyTexClass = plasticMat ->Ks->texClass;  // "color Ks" [a,b,c]: which is a constant texture, whose value is stored in value member
			                                                 //  ConstantTexture. In this case, the value of Ks is returned by FindOneSpectrum("Ks", default)

			
			string roughnessTexClass = plasticMat->roughness->texClass;  // "color opacity" [a,b,c]
			
			//string bumpMapTexClass = plasticMat->bumpMap->texClass;  // "color opacity" [a,b,c]
			
			// check diffuse texture

			if ( diffuseTexClass == "constant" ) {
				
				// add the material parameters Kd, sigma, bumpMap as uniform variables to the shader

			    Reference<ConstantTexture<Spectrum>> constTexture = ( ConstantTexture<Spectrum> *) plasticMat->Kd.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				RGBSpectrum value =constTexture->value;

				float rgbColor[3]; 
				value.ToRGB(rgbColor);

				g_currentMaterial->getUniforms().put("uPlasticMaterial.KdTextureClass", TEX_CONST );
				constTextures.push_back( Cvec3f( rgbColor, 3) );
				g_currentMaterial->getUniforms().put("uPlasticMaterial.KdTextureUnit", constTextureUnit++  );
				
			    
			}

			else if ( diffuseTexClass =="imagemap" ) {
				

				Reference<ImageTexture<RGBSpectrum, Spectrum> > imageTexture = ( ImageTexture<RGBSpectrum, Spectrum> *) plasticMat->Kd.GetPtr(); // Kd: Reference to Texture
				string filename = imageTexture->filename;
				
				
				g_currentMaterial->getUniforms().put("uPlasticMaterial.KdTextureClass", TEX_IMAGEMAP );

			
				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RGB_RGB>(
					                                                              new ShaderImageTexture2D_RGB_RGB( filename.c_str(), true ) ) );
				g_currentMaterial->getUniforms().put("uPlasticMaterial.KdTextureUnit", textureUnit++ );
			}

			else { 

				Error( "unsupported texture class: %d \n", diffuseTexClass);

			}

			// glossy texture
			
			if ( glossyTexClass == "constant" ) {
				
				// add the material parameters Kd, sigma, bumpMap as uniform variables to the shader

			    Reference<ConstantTexture<Spectrum>> constTexture = ( ConstantTexture<Spectrum> *) plasticMat->Ks.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				RGBSpectrum value =constTexture->value;

				float rgbColor[3]; 
				value.ToRGB(rgbColor);

				g_currentMaterial->getUniforms().put("uPlasticMaterial.KsTextureClass", TEX_CONST );
				constTextures.push_back(   Cvec3f( rgbColor, 3) );
				g_currentMaterial->getUniforms().put("uPlasticMaterial.KsTextureUnit",  constTextureUnit++ );
				
			   
			}

			else if ( glossyTexClass =="imagemap" ) {
				

				Reference<ImageTexture<RGBSpectrum, Spectrum> > imageTexture = ( ImageTexture<RGBSpectrum, Spectrum> *) plasticMat->Ks.GetPtr(); // Kd: Reference to Texture
				string filename = imageTexture->filename;
				
				
				g_currentMaterial->getUniforms().put("uPlasticMaterial.KsTextureClass", TEX_IMAGEMAP );
				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RGB_RGB>(
					                                                              new ShaderImageTexture2D_RGB_RGB( filename.c_str(), true ) ) );
				g_currentMaterial->getUniforms().put("uPlasticMaterial.KsTextureUnit", textureUnit++ );


			}

			else { 

				Error( "unsupported texture class: %d \n", glossyTexClass);

			}

	 

	// check roughness  texture
			if ( roughnessTexClass == "constant" ) {
				
			
			    Reference<ConstantTexture<float> > constTexture = ( ConstantTexture<float> *) plasticMat->roughness.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				float  value =constTexture->value;

				g_currentMaterial->getUniforms().put("uPlasticMaterial.roughnessTextureClass", TEX_CONST );
				constTextures.push_back( Cvec3f(value, 0,0)  );
				g_currentMaterial->getUniforms().put("uPlasticMaterial.roughnessTextureUnit", constTextureUnit++  );
				
			    
				
			}

			else if ( roughnessTexClass =="imagemap" ) {  // imageFloat texture
				

				Reference<ImageTexture<float, float> > imageTexture = ( ImageTexture<float, float> *) plasticMat->roughness.GetPtr(); // Kd: Reference to Texture
				string filename = imageTexture->filename;
				
		        
				g_currentMaterial->getUniforms().put("uPlasticMaterial.roughnessTextureClass", TEX_IMAGEMAP );

				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RF_RF>( new ShaderImageTexture2D_RF_RF( filename.c_str() ) ) );

				g_currentMaterial->getUniforms().put("uPlasticMaterial.roughnessTextureUnit", textureUnit++ );
				// put() =>  valueMap[name].reset(new TexturesValue(&value, 1) );

				//currentMaterialShader->getUniforms().put("uTexUnit0", shared_ptr<ShaderImageTexture2D>(
				//	                                          new ShaderImageTexture2D("image/Fieldstone.ppm", true)));



			}

			else { 

				Error( "unsupported texture class: %d \n", roughnessTexClass);

			}

			

			// check bumpMap texture
			if ( plasticMat->bumpMap == NULL ) { // there is no bumpMapTexture => set bumpMapTextureClass uniform to -1. 
				
				g_currentMaterial->getUniforms().put("uUberMaterial.bumpMapTextureClass", -1 );
			
				g_currentMaterial->getUniforms().put("uUberMaterial.bumpMapTextureUnit", -1 );
			}

			else {

			  string bumpMapTexClass  = plasticMat->bumpMap->texClass;  // material may NOT have bumpMap parameter; In that case, 
			                                                        // bumpMap is set to NULL. 
			
			  if ( bumpMapTexClass == "constant" ) {
				
			
			    Reference<ConstantTexture<float>> constTexture = ( ConstantTexture<float> *) plasticMat->bumpMap.GetPtr(); // Kd: Reference to Texture

			    // get the constant color
				float value =constTexture->value;

				
				g_currentMaterial->getUniforms().put("uPlasticMaterial.bumpMapTextureClass", TEX_CONST );
				constTextures.push_back( Cvec3f(value, 0,0) );
				g_currentMaterial->getUniforms().put("uPlasticMaterial.bumpMapTextureUnit", textureUnit++  );
			    
				
			  }

			  else if ( bumpMapTexClass =="imagemap" ) {
				

				Reference<ImageTexture<RGBSpectrum, Spectrum> > imageTexture = ( ImageTexture<RGBSpectrum, Spectrum> *) plasticMat->bumpMap.GetPtr(); // Kd: Reference to Texture
				string filename = imageTexture->filename;
				
				g_currentMaterial->getUniforms().put("uPlasticMaterial.bumpMapTextureClass", TEX_IMAGEMAP );
				imgTextures.push_back( shared_ptr<ShaderImageTexture2D_RF_RF>( new ShaderImageTexture2D_RF_RF( filename.c_str()  ) ) );
				g_currentMaterial->getUniforms().put("uPlasticMaterial.bumpMapTextureUnit", textureUnit++ );
				// put() =>  valueMap[name].reset(new TexturesValue(&value, 1) );

				//currentMaterialShader->getUniforms().put("uTexUnit0", shared_ptr<ShaderImageTexture2D>(
				//	                                          new ShaderImageTexture2D("image/Fieldstone.ppm", true)));



  			  }

			  else { 

				Error( "unsupported texture class: %d \n", bumpMapTexClass);

			  }
			} // uberMat->bumpMap is not NULL

			
			// send the array of textures to the shader
			
			if ( imgTextures.size() != 0 ) {
			    g_currentMaterial->getUniforms().put( "uTextures",  &imgTextures[0], imgTextures.size() );
			}

			if (constTextures.size() != 0 ) {
			   g_currentMaterial->getUniforms().put( "uConstTextures",  &constTextures[0], constTextures.size() );
			}



		} // material->name =="plastic"

		else {
			
			std::cout << " unsupported material: " <<  material->name << endl;
			Error (" unsupported material %s\n", material->name  );

		}


	// process the lights: only consider positional light "point" and directional light "distant" now

	// lightDir = Normalize(LightToWorld(dir));     L = radiance;
	// lightPos = LightToWorld(Point(0,0,0)); 
	//  Intensity = intensity;

	
// normal drawing mode
// Declare an empty uniforms

//	When RiWorldBegin is invoked, all rendering options are frozen and cannot be changed until the picture is finished. 
//		The world-to-camera transformation is set to the current transformation and the current transformation is 
//		reinitialized to the identity. Inside an RiWorldBegin-RiWorldEnd block, 
//		the current transformation is interpreted to be the object-to-world transformation.
 
  // use the skyRbt as the eyeRbt
  // renderOptions->CameraToWorld = Inverse(curTransform); // used to transform the camera coordinates to the world coordinates
	// curTransform is interpreted as the transformation from the  world  space to the camera space

  g_currentMaterial->getUniforms().put("uNumOfLights", (int) lights.size() );   // lights is a member of RenderOptions class

  /* The camera is set so that the rainbow can be seen. The camera position of the scene file is overriden

  g_eyeRbt =  Matrix4( &CameraToWorld[0].GetMatrix().m[0][0], 16);
  
  cout << "g_eyeRbt  = CameraToWorld =\n " << g_eyeRbt << endl;
  messageFile <<  "g_eyeRbt = CameraToWorld = \n" << g_eyeRbt << endl;

  Matrix4 invEyeRbt = inv( g_eyeRbt );
  cout << "invEyeRbt =\n " << invEyeRbt << endl;
  messageFile <<  "invEyeRbt  =\n " << invEyeRbt << endl;

  */

  for (int i=0; i < lights.size(); i++ ) {

	 Light * light = lights[i];
	 if (light->lightType == "point" ) {
		     PointLight *pLight = (PointLight *) light;

			 char name[100];
			 std::string uniformName;

             sprintf(name, "uLights[%d].lightType", i);
             uniformName = name;
			 g_currentMaterial->getUniforms().put(uniformName, LIGHT_POINT); 
				
             sprintf(name, "uLights[%d].lightTPos", i);
             uniformName = name;

			 // this is the global position of the light position
             g_currentMaterial->getUniforms().put("uniformName", 
				                                Cvec3f( pLight->lightPos[0], pLight->lightPos[1],  pLight->lightPos[2] ) );
			 float rgbColor[3];
			 pLight->Intensity.ToRGB( rgbColor );

			 sprintf(name, "uLights[%d].intensity", i);
             uniformName = name;
			 g_currentMaterial->getUniforms().put(uniformName,  Cvec3f( rgbColor, 3 ) );
	 } // point light
	 else if (light->lightType =="distant" ) {
		     DistantLight *distLight = (DistantLight *) light;

		     char name[100];
			 std::string uniformName;

             sprintf(name, "uLights[%d].lightType", i);
             uniformName = name;
		     
			 g_currentMaterial->getUniforms().put(uniformName, LIGHT_DISTANT); 

			 sprintf(name, "uLights[%d].lightDir", i);
             uniformName = name;
				                               
             g_currentMaterial->getUniforms().put(uniformName,  
				                                Cvec3f( distLight->lightDir[0], distLight->lightDir[1],  distLight->lightDir[2] ) );
			 float rgbColor[3];
			 distLight->L.ToRGB( rgbColor );

			 sprintf(name, "uLights[%d].L", i);
             uniformName = name;
			 g_currentMaterial->getUniforms().put(uniformName,  Cvec3f( rgbColor,3 ) );
	 } // distant light

	 else { Error( "not supported light type: %s \n", light->lightType );
	 }


    } // for each light 
		
  
		// add each primitive to the global object list

        char objectName[100];
		sprintf(objectName, "triMesh%d", objId );	// objId = geoPrim->primitiveId;

		cout << "primitive: " << objectName << endl;

		shared_ptr<Object> object ( new Object( objectName, TriangleMeshRbt, triangleMesh, g_currentMaterial) );


			
		Cvec4 pickColor;

		pickColor = g_picker.idToColor(objId+1); // objId = 0 is AABB box, so all other objects start with 1
	
		object->pickColor = pickColor; // this st	ored pickColor will be sent to the shader Uniform variable when drawing 

		//store the objects in the object list


		g_objectPtrList.push_back( object );

		//add the object and the id to the list for picking
		g_picker.addToMap( objId+1, object );
		  

} // for each primtive

} // RenderOptions::GPURender(  )




