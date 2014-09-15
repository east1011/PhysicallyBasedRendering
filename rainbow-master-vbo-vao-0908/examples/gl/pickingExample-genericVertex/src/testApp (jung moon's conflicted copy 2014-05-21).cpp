#include "testApp.h"
#include "ofAppGLFWWindow.h"

ofPtr<ofAppBaseWindow> 		window;  // changed by Moon Jung, 2014/4/23

//#include "ofConstants.h"
//#include "ofGLUtils.h"

//#include "pickerwindow.h" 

//#define USE_PROGRAMMABLE_GL

// from PickingMat

#include "shadergeometry.h" // includes Object class

#include "materialshader.h"
#include "glsupport.h"


#include "picker_0.h"

#include "shadertexture.h"

#include "ppm.h"
#include <stdio.h>
#include "cvec.h"

#include "floatComparison.h"


// for rainbows
#include "spa.h"  //include the SPA header file
//#include "optical_depth.h"
#include "irradiance_sun.h"
//#include "phaseFunction.h"
#include "spectrum_to_rgb.h"

#include "phaseFunction.h"

#include <ctype.h>  // from pbrt for readFloatFile
#include <stdlib.h>  // from pbrt for readFloatFile


using namespace std; // string type is defined in namespace 

#define Assert(expr) \
    ((expr) ? (void)0 : \
        fprintf(stderr, "Assertion \"%s\" failed in %s, line %d", \
               #expr, __FILE__, __LINE__))




#include "stdafx.h"
#include "api.h"
#include "probes.h"
#include "parser.h"
#include "parallel.h"

bool ParseFile(const string &filename);
extern int g_argc;
extern char **g_argv;



// extern GLFWwindow* ofAppGLFWWindow::windowP; => glfwSawpBuffers(windowP)
// static ofPtr<ofAppBaseWindow> 		window; => window.display()


// --------- Materials

shared_ptr<MaterialShader> g_redDiffuseMat,
                            g_blueDiffuseMat,
							g_greenDiffuseMat,
                            g_bumpFloorMat,
                            g_arcballMat,
                            g_pickingMat,
                            g_lightMat,
							g_rainbowMat;

shared_ptr<MaterialShader> g_overridingMaterial;

// PickingMat

using namespace std;      // for string, vector, iostream, shared_ptr and other standard C++ stuff


bool reDrawWindowEvent = true;

bool g_Gl2Compatible = false; // use shader programs of type "*-gl3 *"
int  g_objId = 1; 

// set the parameters for the camera, which define the internal working of the camera.
// All projection happens relative to the origin 0. near and far determine the distance of the clipping planes
//	from the origin. left and right define the opening angle of the frustum together with the near plane. 
//	If you place a "near" and far plane in opposite directions of the origin, 
//	your projection space becomes shaped like an hourglass.
//the contents of the view volume are "projected" onto the near clipping plane.
// The projection is not meant to place the view. The projection is sort of the "lens" of your virtual camera.
//You can use it to tilt and shift your lens and set the "focal length".
//But it's not meant to "place your camera". This should happen through the modelview matrix.

// From the far to the "near" plane, objects will become larger as they approach 0, 
// for Z=0 they're in a singularity and blow up infinitely, and then getting closer to near
// they will become smaller but will be inverted, i.e. rotated by 180° around the Z axis, and depth values being turned around, 
// i.e. depth sorting by depth testing will reject fragments closer to near than to 0.
double  g_frustMinFov = 100.0;  // A minimal of 100 degree field of view for rainbow viewing

double  g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

double  g_frustNear = -0.1;    // near plane
double  g_frustFar = -1000.0;    // far plane

int  g_windowWidth; 
int  g_windowHeight; 


bool g_mouseClickDown = false;    // is the mouse button pressed
bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;

int g_prev_mouseClickX, g_prev_mouseClickY; // coordinates for mouse click event

int  g_pressed_button;

int g_activeShader = 0;

bool g_camera_mode = false;
bool g_backgroundBackedUp = false;

bool g_picking_mode = false; // mode during the process of picking
bool g_picked_mode = false;

bool g_rotation_mode  = false;



GLdouble g_clearColor[4];
GLuint g_fbo; 

bool g_debugMode = true;

// for rainbow drawing
// The coordinates of the AABB box relative to the  world coordinate system

//float g_phaseFunction[ nRadii * nSpectralSamples * nThetas ];

// THis includes the hemisphere for sky
//static const Cvec3 g_initAABBmax(90.3325, 75.8859, 90.3325);
//static  const Cvec3 g_initAABBmin(-90.3325, -14.4466, -90.3325);

// THis includes the hemisphere for sky

static const float rainbowVolumeHeight = 50; 
static const float rainbowVolumeWidth = 100;
static const float rainbowVolumeDepth = 10;

// The coordinates of the AABB, which will be used to compute the sizes of the view volume

//static const Cvec3 g_initAABBmax( 53.0225, 3.86467+ rainbowVolumeHeight , 53.02289 );
//static  const Cvec3 g_initAABBmin(-53.0226, -0.931774, -53.02225);

static const Cvec3 g_initAABBmax( rainbowVolumeWidth / 2.0,  rainbowVolumeHeight , 0.0);
static  const Cvec3 g_initAABBmin( -rainbowVolumeWidth / 2.0, 0.0, -rainbowVolumeDepth);


// The local coordinates of  the AABB box relative to the AABB coordinate system.
// Only these coordinates are important. 


  // relative to the coordinate system at g_initAABBcenter

static      	Cvec3 g_AABBmax, g_AABBmin; 
static const   Cvec3 g_AABBsize = g_initAABBmax - g_initAABBmin;
static           Cvec3 g_initEyeLocation;
static           Cvec3 g_AABBcenterInViewSpace;




shared_ptr<SgRbtNode> g_AABBRbt;  

// --------- Scene
// --------- Scene

Cvec4 g_light1Pos(2.0, 3.0, 14.0, 1.0), g_light2Pos(-2, -3.0, -5.0, 1.0);  // define two lights positions in world space


Cvec4 g_light1Color(1.0, 1.0, 1.0, 1.0);
Cvec4 g_light2Color(1.0, 1.0, 1.0, 1.0);

Matrix4 g_eyeRbt, g_invEyeRbt, g_cameraRotMat, g_invCameraRotMat;

Cvec3  g_sunRayDir ;

float  g_radius = 1.0e-3; // 1 mm

float  g_dropDensity =  30000;

Cvec3  g_lightRGBColor;

float  g_textureFlag = 1.0;  // shader is  supposed to use textures for rendering

//  void makeTranslationMatrix( const ofVec3f& );
//	void makeTranslationMatrix( float, float, float );
// shared_ptr< Picker > g_picker_ptr;
// g_picker_ptr.reset( new Picker() );

Picker g_picker; // default constructor

shared_ptr< Object>  g_currentPickedObject;

// references are not objects, i.e. types regions of storage; so they are not actually "variables";
// you cannot point to references; so there are no references of references, array of references, pointers to
// references. They are just aliases of other variables. So, you cannot use vector < Object & > objectList

vector < shared_ptr<Object> >  g_objectList;

std::ofstream messageFile("./messageFile.txt", std::ios::out);
std::ofstream textureFile("./textureDebug.txt", std::ios::out);

//--------------------------------------------------------------
testApp:: testApp() 

{
}

void testApp::setup(){

	
	// set up the sun light intensities and the sun direction

	setupSun();

	
	// set up the eye position and direction relative to the g_localAABB which is created at
	// initObjects(). Then create the eye Matrix

	// For now, It is assumed that localAABB is centered at the world origin
	// height / horizontalDistFromEyeToAABB  = tan (g_frustFovY / 2.0) 

	setupCamera();

	bool g_Gl2Compatible = false; // use shader programs of type "*-gl3 *"

    //ofSetBackgroundAuto( false); // It will make a single buffering by drawing to the front buffer
	ofSetBackgroundAuto( true);


	g_windowWidth = (float) ofGetWidth(); 
	g_windowHeight = (float) ofGetHeight(); 

	
	//g_windowWidth =  nThetas;  
	//g_windowHeight = nSpectralSamples;  

		
	initGLState();
    initMaterials();
	

	//initPBRTObjects();
	initObjects();
    initRainbowAABB();



}

void testApp::exit() {


}


void testApp::setupSun() {

spa_data spa;  //declare the SPA structure
    int result;
    Cvec3 sunRay;
    
    // Default values to calculate sun's position
    spa.year          = 2014;
    spa.month         = 1;
    spa.day           = 31;
    spa.hour          = 16;
    spa.minute        = 30;
    spa.second        = 30;
    spa.timezone      = 9.0;             // KST (GMT + 9)
    //spa.delta_ut1     =0;
    spa.delta_t       = 67;              // Doesn't it depend on location?  go to the related webpage to find the correct value.
    spa.longitude     = 126.9410634;     // Sogang University's longitude in decimal
    spa.latitude      = 37.5517132;      // Sogang University's latitude in decimal
    spa.elevation     = 200;             // m
    spa.pressure      = 820;             // mbar [milibar]
    spa.temperature   = 23;              // celcius
    //spa.slope         = 30;              // surface slope angle, used to compute the sun's incidence angle to the surface
    //spa.azm_rotation  = -10;             // surface azimuth angle, used to compute the sun's incidence angle to the surface
    spa.atmos_refract = 0.5667;          // ?
    spa.function      = SPA_ZA;          // find the zenith and azimuth angle
	//spa.function      = SPA_ALL;

	
	g_lightRGBColor = getSunLightRGBColor();

	messageFile << "lightRGBColor =" << g_lightRGBColor << endl;

    get_time_and_location(spa.year, spa.month, spa.day, spa.hour, spa.minute, spa.timezone, spa.longitude, spa.latitude);
    
    // calculate sun's position
    result = spa_calculate(&spa);   // input is valid when result eq 0
    
	if ( result !=0 ) {
		messageFile   << "error in the calculation of sun position" << endl;
		return; 
	}


    // calculate sunRay

    calculate_sunRay(sunRay, spa.zenith, spa.azimuth); // vector sunRay points to the sun
	
	messageFile  << "sunRay zenith angle = " << spa.zenith <<" sunRay azimuth" << spa.azimuth << endl;
	
	messageFile  << "sunRay in geocentric coord (world coord system) = (" << sunRay  << endl;

	g_sunRayDir = Cvec3(-sunRay[1], sunRay[2], -sunRay[0]); // rename the axes to 3D graphics convention : z => y, x => -z, y => -x => y=z, z = -x, x=-y
	
	// E.g: In the original coord system: (-10, 20, 5) [ azimth= south-west, polar = positive] =? (-20, 5, 10)
	
    
    messageFile  << "sunRay in Graphics coord = (" << g_sunRayDir << endl;

	
}

void testApp::setupCamera() {

	
	
// set the parameters for the camera, which define the internal working of the camera.
    g_frustMinFov = 100.0;  // A minimal of 100 degree field of view for rainbow viewing

    g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

    g_frustNear = -0.1;    // near plane
    g_frustFar = -200.0;    // far plane


	// set the camera location and direction so that it can see the rainbow well

	Cvec3 yAxis (0,1,0);
	Cvec3 zAxis (0,0,1);
	Cvec3 negZAxis (0,0,-1);

	// These axes define the geocentric coordinate system (North = negZAxis, East = xAxis, and Gravity Directions = yAxis)

	// set the camera ground direction so that it is equal to the ground projection of the sun ray direction

	Cvec3 groundCamDir = (-g_sunRayDir) - yAxis * dot(-g_sunRayDir, yAxis);
	Cvec3 upCamDir = -g_sunRayDir - groundCamDir;

	Assert ( upCamDir == yAxis * dot(-g_sunRayDir, yAxis) );


	
	groundCamDir.normalize();

	//Matrix4 rotMat;


	if ( norm2( negZAxis - groundCamDir ) < 1.0e-8 ) { // no need to rotate; The rotation matrix will be identity
		g_cameraRotMat =  Matrix4();
	}
	
	else {
		Cvec3 rotAxis = cross( negZAxis, groundCamDir); // rotAxis = zAxis x groundCamDir = toEarthDir

	    double  sinTheta = norm( rotAxis ); // length(zAxis x camDir) = norm(zAxis) * norm(groundCamDir) * sin (theta) = sin (theta)
	    double  theta = asin( sinTheta );
		
        g_cameraRotMat = Matrix4::makeAxisRotation( theta, normalize(rotAxis) ); // so that the camera dir = sun dir
	

     	messageFile << "groundCamDir:" ;
	    messageFile << groundCamDir << endl;

		messageFile << "rotation axis:";
	    messageFile << rotAxis << endl;

	    messageFile << "rot theta:";
	    messageFile << theta << endl;

	    messageFile << "rotated Axis:";
	    messageFile << Cvec3(  g_cameraRotMat * Cvec4( negZAxis, 0.0 ) ) << endl;

		messageFile << "g_rotMat :\n" <<  g_cameraRotMat << endl;
		


	    Assert ( norm2 ( groundCamDir - Cvec3(  g_cameraRotMat * Cvec4( negZAxis, 0.0 ) ) ) <= 1.0e-8 );
	   
	}

    
	double eyeHeight = 1.7;
	
	g_invCameraRotMat = inv( g_cameraRotMat);
	
    float distanceToEye =  10;  //  m
	
	// Cvec3 localCenter = Cvec3( g_invCameraRotMat * Cvec4( g_initAABBcenter,0) ); // (-13,25, -2) 

	float volumeDepth = g_initAABBmax[2] - g_initAABBmin[2];


	
	Cvec3 eyeLocation = Cvec3(0, eyeHeight, distanceToEye);

	// the position of the eye relative the world coordinate system.	
	//g_eyeRbt = g_cameraRotMat * Matrix4::makeTranslation( eyeLocation  );
	

	// for debugging

	g_eyeRbt =  Matrix4::makeTranslation( eyeLocation  );

	// The reference frame of the AABB box

	g_AABBcenterInViewSpace = Cvec3( 0, g_AABBsize[1]/2.0, -g_AABBsize[2]/ 2.0 );

	//shared_ptr<SgRbtNode> g_AABBRbt  ( new SgRbtNode ( g_cameraRotMat * 
	//	                                                Matrix4::makeTranslation( 	g_AABBcenterInViewSpace  ) ) ); 
	g_AABBRbt.reset(  new SgRbtNode ( g_cameraRotMat * 
		                                                Matrix4::makeTranslation( 	g_AABBcenterInViewSpace  ) ) ); 
	
	g_AABBmin = Cvec3(  *g_AABBRbt * Cvec4( -g_AABBsize[0]/2.0, -g_AABBsize[1]/2.0, -g_AABBsize[2]/2.0, 1.0 )  );
	g_AABBmax = Cvec3(  *g_AABBRbt * Cvec4( g_AABBsize[0]/2.0, g_AABBsize[1]/2.0, g_AABBsize[2]/2.0, 1.0 )  );
	
	g_invCameraRotMat = inv( g_cameraRotMat);
	messageFile << "g_invRotMat: \n" << g_invCameraRotMat << endl;

	Matrix4 g_invEyeRbt = inv( g_eyeRbt);
}



Cvec3  testApp::getSunLightRGBColor() {
// convert the sun light to its RGB representation
	
	double  sunIntensity, X = 0, Y = 0, Z = 0, XYZ;
	double lambda_m; // meter

	double xBar, yBar, zBar;

	//const double lambdaStart = 400;   // 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
    //const double lambdaEnd = 700;
     

	/* cie_colour_match[(lambda - 380) / 5][0] = xBar
      cie_colour_match[(lambda - 380) / 5][1] = yBar
      cie_colour_match[(lambda - 380) / 5][2] = zBar
	  ==> cie_coulor_match [] from 380 nanometer to 780 nanometer
    */
	double lambda;
	int j;
	
	for (lambda = lambdaStart, j=0; j < nSpectralSamples; lambda += lambdaStep, j++ ) {
	 // integration of cie-color-match curve over  nSpectralSampleSteps intervals,
    	      	
     //  cie is sampled every lambdaStep * 2
		  
	   int i1 = int ( (lambda - 380 ) / ( lambdaStep*2 ) );
	   int i2 = int ( ( lambda + lambdaStep - 380) / (lambdaStep*2 ) );

	   messageFile << "j =" << j  << " i1=" << i1 << "i2=" << i2 << endl;

	   if ( i1 == i2) { 
		   xBar =  cie_colour_match[ i1][ 0]; 
		   yBar =  cie_colour_match[ i1][ 1 ];
		   zBar =  cie_colour_match[ i1 ][ 2 ];
		
		}  
	   else { // j = i+1 => lambda is between two sample points 
		   
		  xBar = ( cie_colour_match[ i1][ 0] + cie_colour_match[ i2 ][ 0] ) / 2.0; 
		  yBar = ( cie_colour_match[ i1 ][ 1 ] + cie_colour_match[ i2][ 1 ] ) /2.0;
		  zBar = ( cie_colour_match[ i1][ 2 ] + cie_colour_match[ i2][ 2 ] ) /2.0;
        }

		//intensity = bb_spectrum(lambda);  // You already have the function that computes the irradiance of the sun
		/*
		if (  definitelyLessThan (  lambda / (lambdaStep * 2) - floor ( lambda / (lambdaStep * 2) ), 
		                                   (2 * lambdaStep ) / 10.0, epsilon ) ) {
		 
		 // lambda is at a cie sample lambda which is sampled every lambdaStep * 2
		  
		   int i = int (  lambda - 380 ) / int ( lambdaStep*2 );

		   xBar =  cie_colour_match[ i][ 0]; 
		   yBar =  cie_colour_match[ i][ 1 ];
		   zBar =  cie_colour_match[ i ][ 2 ];
		
		}  
		else { // lambda is between two points and interpolate between them

		  int i1 =  int (  lambda - lambdaStep - 380 ) / int ( lambdaStep*2 );
		  int i2 =  int (  lambda + lambdaStep  - 380 ) / int ( lambdaStep*2 );
		
		  xBar = ( cie_colour_match[ i1][ 0] + cie_colour_match[ i2 ][ 0] ) / 2.0; 
		  yBar = ( cie_colour_match[ i1 ][ 1 ] + cie_colour_match[ i2][ 1 ] ) /2.0;
		  zBar = ( cie_colour_match[ i1][ 2 ] + cie_colour_match[ i2][ 2 ] ) /2.0;
        }
		*/

		/*
		xBar =  ( cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep*2)] [0]
				         + cie_colour_match[ ( int(lambda + lambdaStep*2) - 380 )  / int( lambdaStep*2)] [0] ) /2.0; 
				 // lambdaStep is 2.5 nanometer
		yBar = ( cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep*2)] [1] 
		          + cie_colour_match[ ( int(lambda + lambdaStep*2) - 380 )  / int( lambdaStep*2)] [1] ) /2.0; 
		zBar = ( cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep*2)] [2] 
		           + cie_colour_match[ ( int(lambda + lambdaStep*2) - 380 )  / int( lambdaStep*2)] [2] ) /2.0; 

		//messageFile << "In dropvolume.cpp: int(lambda)- 380 =" << int(lambda) - 380 << endl;
		*/

		lambda_m = lambda * 1.0e-9;

        sunIntensity = calculate_irradiance_of_sun(lambda_m); // isun: the return value of irradiance is per nanometer, but the function uses
		                                                       // lambda_m in unit meter

		X += sunIntensity * xBar * lambdaStep; 
		Y += sunIntensity * yBar * lambdaStep;
		Z += sunIntensity * zBar * lambdaStep;
		
		
    }
	

	XYZ = (X+Y+Z);
	Cvec3 light_xyzColor = Cvec3( X/XYZ, Y/XYZ, Z/XYZ );
	
	return xyz_to_rgb(cs, light_xyzColor);

}




Matrix4 testApp::makeProjectionMatrix() {

	return Matrix4::makeProjection(
		g_frustFovY, g_windowWidth / static_cast <double> (g_windowHeight),
		g_frustNear, g_frustFar);
}


// update g_frustFovY from g_frustMinFov, g_windowWidth, and g_windowHeight
void testApp::updateFrustFovY() {

  if (g_windowWidth >= g_windowHeight)

    g_frustFovY = g_frustMinFov;

  else {
    const double RAD_PER_DEG = 0.5 * 3.14159 /180;

    g_frustFovY = atan2(sin(g_frustMinFov * RAD_PER_DEG) * g_windowHeight / g_windowWidth, cos(g_frustMinFov * RAD_PER_DEG)) / RAD_PER_DEG;
  }
}

//--------------------------------------------------------------
void testApp::update(){

}

static bool backgroundBackedUp = false;

void testApp::draw()  { // callback function for draw event


	//if ( !reDrawWindowEvent ) return; // if there is no need to redraw the window, skip draw() method.


	//initGLState();

	if (backgroundBackedUp) { // If background color has been backed up, this should be used
		// for ordinary drawing. 
		//Now set back the clear color in order to draw the real scene
		glClearColor( g_clearColor[0], g_clearColor[1],g_clearColor[2], g_clearColor[3] );   
	}

	//g_windowWidth = (float) ofGetWidth(); 
	//g_windowHeight = (float) ofGetHeight(); 
	//g_windowWidth = (float) nThetas;
	//g_windowHeight = (float) nSpectralSamples; 

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

	bool picking = false;

	drawStuff( picking);

	try {
	  checkGlErrors(" after drawStuff() in draw()" );
	}
	catch (const runtime_error & error ) {
      std::cout << error.what() << endl;
	  messageFile << error.what() << endl;

	}

	//reDrawWindowEvent = false;

	//glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)

} // draw()

void testApp::drawForPicking() {
 
  // We need to set the clear color to black, for pick rendering.
  // so let's save the original clear color

  //initGLState();

  if ( !backgroundBackedUp ) { // is the background color backed up? it should be backed up only once, in the 
	                          // beginning. Otherwise, the background color for picking will be backed up and 
	                          // used for normal scene rendering

  glGetDoublev(GL_COLOR_CLEAR_VALUE, g_clearColor); // glGetIntegerv() // glGetDoublev() integer or double values
  backgroundBackedUp = true;
  }


  // clear the background to transparent black color for picking

  
  glClearColor(0, 0, 0, 0); // this black background color is only meant to be used
                               // for the background of the back buffer for rendering in picking.


  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // CLEAR THE COLOR AND DEPTH BUFFER

  bool picking = true;

  drawStuff( picking);


  checkGlErrors(" drawStuff() in drawForPicking()");

  glFinish();
  checkGlErrors("glFinish() in drawForPicking()" );
}



void testApp::drawStuff( bool picking ) {

// normal drawing mode
// Declare an empty uniforms

  Uniforms extraUniforms;
  // build & send proj. matrix to vshader

  
  const Matrix4 projMatrix = makeProjectionMatrix();


  extraUniforms.put("uProjectionMatrix", projMatrix);

  // Uniforms.put() is not done for uScatTex, uDistTex  ???


  // use the skyRbt as the eyeRbt

 // const Matrix4 eyeRbt = g_eyeRbt;

  g_invEyeRbt = inv(g_eyeRbt);

 

  const Cvec3 eyeLight1 = Cvec3(g_invEyeRbt * g_light1Pos ); // g_light1 position in eye coordinates
  const Cvec3 eyeLight2 = Cvec3(g_invEyeRbt * g_light2Pos ); // g_light2 position in eye coordinates
  
  // send the eye space coordinates of lights to uniforms
  extraUniforms.put("uLight1Pos", eyeLight1);
  extraUniforms.put("uLight2Pos", eyeLight2);

  
  //safe_glUniform3f(curSS.h_uLight1Pos,  eyeLight1[0], eyeLight1[1], eyeLight1[2] );
  //safe_glUniform3f(curSS.h_uLight2Pos,  eyeLight2[0], eyeLight2[1], eyeLight2[2] );


  extraUniforms.put("uLight1Color", g_light1Color );
  extraUniforms.put("uLight2Color", g_light2Color );

  
  

 // extraUniforms.put("uTextureFlag", g_textureFlag);
  
  /*
  safe_glUniform4f(curSS.h_uLight1Color,   g_light1Color[0], g_light1Color[1], g_light1Color[2], g_light1Color[3] );
  safe_glUniform4f(curSS.h_uLight2Color,   g_light2Color[0], g_light2Color[1], g_light2Color[2], g_light2Color[3] );
  */
 //extraUniforms.put("uTextureFlag", g_textureFlag);
  
//  extraUniforms.put("uInvEyeMatrix", invEyeRbt);

  //Matrix4 MVMAABB = invEyeRbt * g_modelMatrixAABB; 

  // draw the objects in the objectList 
 
                                           // eye space coordinates to global coordinates
  // upload the sun direction relative to the eye space

 // extraUniforms.put("uSunRayDir", Cvec3( g_invEyeRbt * Cvec4(g_sunRayDir,0) ) );
  extraUniforms.put("uSunRayDir", g_sunRayDir );
  extraUniforms.put("uLightRGBColor", g_lightRGBColor);

  if ( !picking ) { // not picking mode but an ordinary drawing mode
	   
 
	// ordinary draw setting
	 glBindFramebuffer(GL_FRAMEBUFFER, 0);
		
	 glViewport(0,0,g_windowWidth, g_windowHeight); 
				
        // g_windowWidth and g_windowHeight is updated when the window is resized

     extraUniforms.put("uWindowWidth", g_windowWidth );
	 extraUniforms.put("uWindowHeight",g_windowHeight);
	  


   // g_overridingMaterial = g_pickingMat;

   for ( int i  = 0; i < g_objectList.size(); i++ ) { // (5) This loop goes over the list of objects g_objectList and draw each object in the list.
	                                                //     Where and how is this object list created?

//	 for ( int i  = 2; i < 3; i++ ) { // for debugging

     Matrix4 MVM = g_invEyeRbt * ( *( g_objectList[i] -> objectRbt ) ) ; // g_currentPickedObject->objectRbt may have been
	                                                              // changed by picking
    //Matrix4 NMVM = normalMatrix(MVM);
  
	 extraUniforms.put("uModelViewMatrix", MVM).put("uNormalMatrix", normalMatrix(MVM) );

	 if ( g_objectList[i]->objectName == "AABBrainbow" ) {

       extraUniforms.put("uEyeMatrix", g_eyeRbt); // This will be used in shader to convert the
       extraUniforms.put("uInvEyeMatrix", g_invEyeRbt); // This will be used in shader to convert the
       extraUniforms.put("uRadius", g_radius );
       extraUniforms.put("uDropDensity", g_dropDensity );

	    // AABBmin and AABBmax with respect to camera space: The inversion will be done in fragment shader


	   //Cvec3 AABBmin = Cvec3(  g_invEyeRbt * Cvec4( g_AABBmin, 1.0) );
	   // Cvec3 AABBmax =  Cvec3( g_invEyeRbt * Cvec4( g_AABBmax, 1.0) );

	   extraUniforms.put("uAABBmin", g_AABBmin );   
	   extraUniforms.put("uAABBmax", g_AABBmax );   
	 }

   
	  
	 if ( g_objectList[i]->objectName == "ground" ) {
		 
		 extraUniforms.put("uTextureFlag", g_textureFlag);
	 }
	  
	// Cvec4 pickColor = g_objectList[i]->pickColor;
	 //Cvec4 materialColor = g_objectList[i]->objectColor;

	 //std::cout << "Linear [float] Color for the Object to be Drawn For Picking=" << pickColor[0] << ","<< pickColor[1] << "," << 
	//	 pickColor[2]  << "," << pickColor[3] << endl;

	 //g_objectList[i]->material->getUniforms().put("uMaterialColor", pickColor);
	 // for debugging

	// g_overridingMaterial-> getUniforms().put("uMaterialColor", pickColor);

     // UniformextraUniforms.put( "uMaterialColor", pickColor ); // set material color. In the case of
	 // g_overridingMaterial, there are no uniform "uMaterialColor" assigned when
	 // when it is created.
	  


	 // draw
	 try {
      g_objectList[i]->draw( extraUniforms ); 
	  // This draw leads to  BufferedObjectGeometry::draw() which binds vertex buffers 
	  // and set vertex     attribute pointers, and then calls drawElements() */
	 
      checkGlErrors(" g_objectList[i]->draw( ) in drawStuff()" );
	  
	 }                                   


	 catch ( const runtime_error & error ) {
		 //std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	 }

 
    } // for each object


	//g_overridingMaterial.reset();

  } // not picking

  else { // picking mode

     g_overridingMaterial = g_pickingMat;

     for ( int i  = 0; i <g_objectList.size(); i++ ) { // (5) This loop goes over the list of objects g_objectList and draw each object in the list.
	                                                //     Where and how is this object list created?

     Matrix4 MVM = g_invEyeRbt * ( *(g_objectList[i]->objectRbt) ) ; // g_currentPickedObject->objectRbt may have been
	                                                              // changed by picking, object includes AABB
    //Matrix4 NMVM = normalMatrix(MVM);
  
	 extraUniforms.put("uModelViewMatrix", MVM).put("uNormalMatrix", normalMatrix(MVM) );

	 
    Cvec4 pickColor = g_objectList[i]->pickColor;

	//Cvec4 materialColor = g_objectList[i]->objectColor;
	
	//std::cout << "Linear [float] Color for the Object to be Drawn For Picking=" << pickColor[0] << ","<< pickColor[1] << "," << 
	//            pickColor[2]  << "," << pickColor[3] << endl;

    extraUniforms.put( "uMaterialColor", pickColor ); // set material color. In the case of
	                                                 // g_overridingMaterial, there are no uniform "uMaterialColor" assigned when
	                                                  // when it is created.


    checkGlErrors(" pick color in drawStuff( picking) " );



    // the following draw() method will use g_overridingMaterial->draw() if g_overridingMaterial is on

	 g_objectList[i]->draw( extraUniforms );  // the material Color is also a uniform value, but it is
	                                         // part of Material of this object. It will be taken care of 
	                                         // within the draw method.

     checkGlErrors(" g_objectList[i]->draw( ) in drawStuff(picking)" );
    }  // for

	// unset the overriding material
    g_overridingMaterial.reset();
	
  } // picking mode

} // drawStuff() for ordinary rendering or  picking

void testApp::readPixels() {

	unsigned char pixel[4] = {0};  // If you are using character types as numbers, use unsigned cha or signed char:
	// or unsigned int: you need to use it because you are dealing with bits in memory directly; or 
	// when doing manipulations such as bit masking or shifting on data,
	// for example when writing low-level code to read binary file formats such as audio files; 

  // read the pixel at (x,y), get the color of the pixel
  // get the object id corresponding to the color
  // return the Rbt node corressponding to the ide

   //vector<float> rcolors(4);
  // glReadPixels(x, y, 1, 1, GL_RGBA, GL_FLOAT, &rcolors[0]);

 // for debugging

  glReadBuffer(GL_BACK); // default: BACK is where we draw in double buffering

 // glReadBuffer(GL_FRONT);

// write the read pixel data to phaseFile

	for (int j = 0; j < g_windowHeight; j++ )
	 for ( int i = 0; i < g_windowWidth ; i++) 
	 {
		   glReadPixels(i,j, 1,1, GL_RGBA, GL_UNSIGNED_BYTE, pixel ); // GLvoid *pixels
		   // ostream& operator<< (unsigned int val);

		   textureFile << "pixel: " << "(" <<  j   << "," <<  i  << ") = "  
			                            << (int) pixel[0] << ", " <<  (int)  pixel[1] << ", " <<  (int) pixel[2]  << endl; 
		   cout << "pixel: " << "(" <<   j  << "," <<  i  << ") = " 
			                            << (int)  pixel[0] << ", " <<  (int) pixel[1] << ", " <<   (int)  pixel[2]  << endl; 

	}
	
    textureFile.close();

	
}

 void testApp::renderToFbo() {

 int   d_windowWidth =  g_windowWidth;
 int   d_windowHeight = g_windowHeight;
 


 // render to the FBO (the default framebuffer will be also rendered later)

 // create a FBO so that we can render to the texture associated with it as a debugg texture

	 Uniforms extraUniforms;     
  
	 
  const Matrix4 projMatrix = makeProjectionMatrix();


  extraUniforms.put("uProjMatrix", projMatrix);

  // Uniforms.put() is not done for uScatTex, uDistTex  ???


  // use the skyRbt as the eyeRbt

 // const Matrix4 eyeRbt = g_eyeRbt;
  g_invEyeRbt = inv(g_eyeRbt);

 

 // const Cvec3 eyeLight1 = Cvec3(g_invEyeRbt * g_light1Pos ); // g_light1 position in eye coordinates
 // const Cvec3 eyeLight2 = Cvec3(g_invEyeRbt * g_light2Pos ); // g_light2 position in eye coordinates
    

  //Matrix4 MVMAABB = invEyeRbt * g_modelMatrixAABB; 

                                           // eye space coordinates to global coordinates
  // upload the sun direction relative to the eye space

  extraUniforms.put("uSunRayDir", g_sunRayDir );
  extraUniforms.put("uLightRGBColor", g_lightRGBColor);

     extraUniforms.put("uWindowWidth", d_windowWidth );
     extraUniforms.put("uWindowHeight", d_windowHeight );

     glGenFramebuffers(1, &g_fbo); 
     glBindFramebuffer( GL_FRAMEBUFFER, g_fbo);  // a container for textures and an optional depth buffer
	//This is executed when drawing depending on the current mode: drawing or debug

	// The texture we're going to render to 
	
     GLuint debugTex; 
     glGenTextures(1, &debugTex);  
	
	// "Bind" the created texture : 
	// all future texture functions will modify this texture 
     glBindTexture(GL_TEXTURE_2D, debugTex);   
	// Give an empty image to OpenGL ( the last "0" ) 
	
	
	// Poor filtering. Needed ! 
	
	    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST); 
	    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); 
	    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE); 
	    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	   // glTexImage2D(GL_TEXTURE_2D, 0,GL_R32F, d_windowWidth, d_windowHeight, 0,GL_RGB,
		//         GL_FLOAT, NULL);   
		 glTexImage2D(GL_TEXTURE_2D, 0,GL_RGBA, d_windowWidth, d_windowHeight, 0,GL_RGBA,
		         GL_UNSIGNED_BYTE, NULL);   
	
	// If a fragment shader is active and it writes a value to the output variable gl_FragColor, 
	// then that value will 
	// be written into each of the buffers specified by bufs.

	// attach texture to framebuffer:  you could also use glFramebufferTexture2D with GL_TEXTURE_2D

	    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,  debugTex, 0); // 0 = mipmap level  
	
	    GLenum e = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
	    if (e != GL_FRAMEBUFFER_COMPLETE)     
	     	printf("There is a problem with the FBO\n");

	// // attach multiple colors
	//for (int i = 0; i < n; ++i) {
	//glFramebufferTexture(GL_DRAW_FRAMEBUFFER,  
	//GL_COLOR_ATTACHMENT0 + i , colorTex[i], 0); } 
	// check if everything is OK 
	
	
	

	   // glBindFramebuffer(GL_FRAMEBUFFER, g_fbo); // GL_FRAMEBUFFER is quiv to GL_DRAW_FRAMEBUFFER

		// define the index array for the outputs 
		// the buffer selection in glDrawBuffers is part of the framebuffer object state. 
		// Therefore, if this setting is constant, this function can be called only once when creating the framebuffer.
		// glDrawBuffers defines an array of buffers into which outputs from the fragment shader data will be written
		//  If a fragment shader writes a value to one or more user defined output variables, then the value of each variable
		// will be written into the buffer specified at a location within bufs
		// bufs \{ GL_FRONT_LEFT, ... GL_COLOR_ATTACHMENTn}


		GLuint attachments[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1,GL_COLOR_ATTACHMENT2 }; 
		glDrawBuffers(3,  attachments);

	    glViewport(0,0,d_windowWidth, d_windowHeight); // set the viewport for the new framebuffer

	// Render on the whole framebuffer, complete from the lower left corner to the upper right 
	// although only the default framebuffer will be visible on your screen, 
	//	you can read any framebuffer that is currently bound with a call to glReadPixels
	
	
   // g_overridingMaterial = g_pickingMat;

    for ( int i  = 0; i < g_objectList.size(); i++ ) { // (5) This loop goes over the list of objects g_objectList and draw each object in the list.
	                                                //     Where and how is this object list created?

    //for ( int i  = 0; i < 1; i++ ) { // for debugging

     Matrix4 MVM = g_invEyeRbt * ( *( g_objectList[i] -> objectRbt ) ) ; // g_currentPickedObject->objectRbt may have been
	                                                              // changed by picking
    //Matrix4 NMVM = normalMatrix(MVM);
  
	 extraUniforms.put("uModelViewMatrix", MVM).put("uNormalMatrix", normalMatrix(MVM) );

	 if ( g_objectList[i]->objectName == "AABBrainbow" ) {

       extraUniforms.put("uEyeMatrix", g_eyeRbt); // This will be used in shader to convert the
       extraUniforms.put("uInvEyeMatrix", g_invEyeRbt); // This will be used in shader to convert the
       extraUniforms.put("uRadius", g_radius );
       extraUniforms.put("uDropDensity", g_dropDensity );

     
       // AABBmin and AABBmax with respect to camera space: The inversion will be done in fragment shader


	   //Cvec3 AABBmin = Cvec3(  g_invEyeRbt * Cvec4( g_AABBmin, 1.0) );
	   //Cvec3 AABBmax =  Cvec3( g_invEyeRbt * Cvec4( g_AABBmax, 1.0) );

	   extraUniforms.put("uAABBmin", g_AABBmin );   
	   extraUniforms.put("uAABBmax", g_AABBmax );   
       //extraUniforms.put("uInvEyeMatrix", g_invEyeRbt);

	   //messageFile << "uAABBmin \n" << Cvec3( g_invEyeRbt * Cvec4( newAABBmin,1) ) << endl;
	   //messageFile << "uAABBmax \n" << Cvec3( g_invEyeRbt * Cvec4( newAABBmax,1) ) << endl;

	 }
	  
	 //if ( g_objectList[i]->objectName == "ground" ) {
		 
	//	 extraUniforms.put("uTextureFlag", g_textureFlag);
	 //}
	  
	
	 // draw
	 try {
      g_objectList[i]->draw( extraUniforms ); 
	  // This draw leads to  BufferedObjectGeometry::draw() which binds vertex buffers 
	  // and set vertex     attribute pointers, and then calls drawElements() */
	 
      checkGlErrors(" g_objectList[i]->draw( ) in drawStuff()" );
	  
	 }                                   


	 catch ( const runtime_error & error ) {
		 //std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		cout << error.what() << endl;

		throw; // A throw expression that has no operand re-throws the exception currently being handled

	 }

	
 
    } // for each object

	// After drawing, read buffer and print its content for debugging purpose

	  
   //vector<float> data( d_windowWidth * d_windowHeight) ; // GL_R32F format: 32 bit (float) for each pixel
	 unsigned char pixel[4];
	 
  // glReadPixels(x, y, 1, 1, GL_RGBA, GL_FLOAT, &rcolors[0]);

 // for debugging

 // glReadBuffer(GL_BACK); // the current buffer you're drawing to

    glReadBuffer(GL_COLOR_ATTACHMENT0); // fragColor = GL_COLOR_ATTACHMENT0

	checkGlErrors(" after glReadBuffer in debugMode" );
 // glReadPixels(x,y, width, height, GL_RGBA, GL_UNSIGNED_BYTE, &pixel[0]);
 
   // glReadPixels(0, 0, d_windowWidth, d_windowHeight, GL_RED, GL_FLOAT, &data[0]);

	checkGlErrors(" after glReadPixels in debugMode" );

	// write the read pixel data to phaseFile
	for (int j = 0; j < d_windowHeight; j++ )
	 for ( int i = 0; i < d_windowWidth ; i++) 
	 {
		   glReadPixels(i,j, 1,1, GL_RGBA, GL_UNSIGNED_BYTE, pixel ); // GLvoid *pixels
		   checkGlErrors(" after glReadPixels in debugMode" );

		   textureFile << "pixel: " << "(" <<  j   << "," <<  i  << ") = "  
			                            <<   pixel[0] << " " << pixel[1] << " " <<  pixel[2]   << endl; 
		   cout << "pixel: " << "(" <<   j  << "," <<  i  << ") = " 
			                            << pixel[0] << " " <<  pixel[1] <<" " << pixel[2] << endl; 

	}
	
    textureFile.close();


	glBindFramebuffer( GL_FRAMEBUFFER, 0); // back to the default buffer
	
} // renderToFbo



//If you're not using VAOs, then you would usually call glVertexAttribPointer
//(and the corresponding glEnableVertexAttribArray) right before rendering to setup the state properly. 
//If using VAOs though, you actually call it (and the enable function) inside the VAO creation code
//(which is usually part of some initialization or object creation), since its settings are stored inside the VAO and 
//all you need to do when rendering is bind the VAO and call a draw function.

/*
So in modern OpenGL using VAOs (which is recommended), it's usually similar to this workflow:

//initialization
  glGenVertexArrays
  glBindVertexArray

glGenBuffers
glBindBuffer
glBufferData

glVertexAttribPointer
glEnableVertexAttribArray

  glBindVertexArray(0)

  glDeleteBuffers //you can already delete it after the VAO is unbound, since the
//VAO still references it, keeping it alive (see comments below).

...

//rendering
  glBindVertexArray
glDrawWhatever

--------------------------------------------------------------------------------

*/

/* http://www.lighthouse3d.com/tutorials/glsl-tutorial/the-normal-matrix/ 

			this is needed when the modelview matrix contains a non-uniform scale.
			the correct matrix to transform the normal is the transpose of the inverse of the M matrix. 
			OpenGL computes this for us in the gl_NormalMatrix.
			
*/

//--------------------------------------------------------------
/*
void testApp::initObjects() {


	Options options;
	vector<string> filenames;
	// Process command-line arguments
	for (int i = 1; i < g_argc; ++i) {
		if (!strcmp(g_argv[i], "--ncores")) options.nCores = atoi(g_argv[++i]);
		else if (!strcmp(g_argv[i], "--outfile")) options.imageFile =g_argv[++i];
		else if (!strcmp(g_argv[i], "--quick")) options.quickRender = true;
		else if (!strcmp(g_argv[i], "--quiet")) options.quiet = true;
		else if (!strcmp(g_argv[i], "--verbose")) options.verbose = true;
		else if (!strcmp(g_argv[i], "--help") || !strcmp(g_argv[i], "-h")) {
			printf("usage: pbrt [--ncores n] [--outfile filename] [--quick] [--quiet] "
				"[--verbose] [--help] <filename.pbrt> ...\n");
			return ;
		}
		else filenames.push_back(g_argv[i]);
	}

	// Print welcome banner
	if (!options.quiet) {
		printf("pbrt version %s of %s at %s [Detected %d core(s)]\n",
			PBRT_VERSION, __DATE__, __TIME__, NumSystemCores());
		printf("Copyright (c)1998-2012 Matt Pharr and Greg Humphreys.\n");
		printf("The source code to pbrt (but *not* the book contents) is covered by the BSD License.\n");
		printf("See the file LICENSE.txt for the conditions of the license.\n");
		fflush(stdout);
	}
	pbrtInit(options);
	// Process scene description
	PBRT_STARTED_PARSING();
	if (filenames.size() == 0) {
		// Parse scene from standard input
		ParseFile("-");
	} else {
		// Parse scene from input files
		for (int i = 0; i < filenames.size(); i++)
			if (!ParseFile(filenames[i]))
				Error("Couldn't open scene file \"%s\"", filenames[i].c_str());
	}
	pbrtCleanup();




}

*/ 



void testApp::initGLState() {
	//glBindFramebuffer(GL_READ_FRAMEBUFFER, 0); 

	messageFile  << "OpenGL Version : " << glGetString(GL_VERSION) << endl; 

	glClearColor(0./255., 0./255., 0./255., 1.0); // sky color

	//glClearColor(128./255., 0./255., 0./255., 1.0); // sky color

	glClearDepth(0.);
	//glClearDepth(1.0); // the default value is 0.0


	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glCullFace(GL_BACK);
	//glEnable(GL_CULL_FACE);
	glDisable(GL_CULL_FACE);

	glEnable(GL_DEPTH_TEST);
	//glDisable(GL_DEPTH_TEST);

	glDepthFunc(GL_GREATER);
	//glDepthFunc(GL_LEQUAL);


	glColorMask( GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

	glDisable(GL_BLEND);
	glDisable(GL_ALPHA_TEST);

	// The default framebuffer has buffers like GL_FRONT​, GL_BACK​, GL_AUXi​, GL_ACCUM​, and so forth. 
	// FBOs do not have these. Instead, FBOs have a different set of images.
	// Each FBO image represents an attachment point, a location in the FBO where an image can be attached.
	//   FBOs have the following attachment points:

    //      GL_COLOR_ATTACHMENTi, GL_DEPTH_ATTACHMENT, GL_STENCIL_ATTACHMENT,GL_DEPTH_STENCIL_ATTACHMENT​:
		
	//Each framebuffer object has a set of attachment points that logical buffers can be attached to. 
	
	// There are three different ways to switch between framebuffers. 

	//The first one is to use several different FBOs, one for each combination of logical buffers 
	//	that you plan on using in you application. To change render targets you simply call glBindFramebufferEXT() 
	//	with the FBO containing the setup you wish to render to.
	// Another way is to use a single FBO and alter the attachments
	// The third way is to use a single FBO, but instead of altering attachments you call glDrawBuffer() or glDrawBuffers() 
	// to change which color attachment(s) the rendering goes to.

	// "renderbuffers = sub-buffers that constitute a framebuffer

	glBindFramebuffer(GL_FRAMEBUFFER, 0); // The default frame buffer 0 will be used both for reading and writing pixels.

	
	GLuint attachments[1] = { GL_BACK_LEFT }; 
	glDrawBuffers(1,  attachments);

	//glDrawBuffer(GL_BACK); // GL_BACK is the draw buffer for the zero framebuffer; for double-buffering context
	
	//// we  use glDrawBuffers(GLsizei n, const GLenum * buf) for drawing to multiple color attachments (targets)

	glReadBuffer(GL_FRONT); // default, the same with glDrawBuffer()

	if (!g_Gl2Compatible)
		glEnable(GL_FRAMEBUFFER_SRGB); // For opengl3.x, enable the sRGB framebuffer update and blending capability

	try {
		checkGlErrors("In initGLStates");
	}
	catch (const runtime_error & error ) {
		std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		//throw;
	}

/* HOW TO USE Framebuffer objects in addition to the default frame buffer

// Renderbuffer is simply a data storage object containing a single image of a renderable internal format. 
 // It is used to store OpenGL logical buffers
//  that do not have corresponding texture format, such as stencil or depth buffer.

 // create a FBO so that we can render to the texture associated with it
// The framebuffer, which regroups 0, 1, or more textures, and 0 or 1 depth buffer. 
// 	void glGenFramebuffers(GLsizei n, GLuint *ids);
	GLuint fbo; 

	glGenFramebuffers(1, &fbo); 
	glBindFramebuffer( GL_DRAW_FRAMEBUFFER, fbo);

	// The texture we're going to render to 
	
	GLuint colorTex; glGenTextures(1, &colorTex);  
	
	// "Bind" the created texture : 
	// all future texture functions will modify this texture 
	glBindTexture(GL_TEXTURE_2D, colorTex);   
	// Give an empty image to OpenGL ( the last "0" ) 
	
	glTexImage2D(GL_TEXTURE_2D, 0,GL_RGB, g_windowWidth, g_windowHeight, 0,GL_RGB, GL_UNSIGNED_BYTE, NULL);   
	// Poor filtering. Needed ! 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE); 
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	
	GLuint  depthTex; glGenTextures(1, &depthTex);
   	// "Bind" the reated texture : all future texture functions will modify this texture 
	glBindTexture(GL_TEXTURE_2D, depthTex);   
	// Give an empty image to OpenGL ( the last "0" ) 
	
	glTexImage2D(GL_TEXTURE_2D, 0,GL_DEPTH_COMPONENT, g_windowWidth, g_windowHeight, 0,GL_DEPTH_COMPONENT,
		         GL_FLOAT, NULL);   
	// Poor filtering. Needed ! 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE); 
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);


	
	// If a fragment shader is active and it writes a value to the output variable gl_FragColor, 
	// then that value will 
	// be written into each of the buffers specified by bufs.

	// attach texture to framebuffer
	glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0,  colorTex, 0); // 0 = mipmap level  
	glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_COMPONENT,  depthTex, 0); // 0 = mipmap level  
	
	// // attach multiple colors
	//for (int i = 0; i < n; ++i) {
	//glFramebufferTexture(GL_DRAW_FRAMEBUFFER,  
	//GL_COLOR_ATTACHMENT0 + i , colorTex[i], 0); } 
	// check if everything is OK 
	
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
	GLenum e = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
	if (e != GL_FRAMEBUFFER_COMPLETE)     
		printf("There is a problem with the FBO\n");
	
	glViewport(0,0,g_windowWidth, g_windowHeight); 
	// Render on the whole framebuffer, complete from the lower left corner to the upper right 
	// although only the default framebuffer will be visible on your screen, 
	//	you can read any framebuffer that is currently bound with a call to glReadPixels
	
	//Let’s assume we have a framebuffer object with n color attachments. 
	//To actually be able to render to all, or some, of these targets simultaneously
	// we must enable them in the OpenGL application with function glDrawBuffers:

//  set the color attachments for drawing  
	GLuint attachments[2] = { GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
	glDrawBuffers(2,  attachments);
	==>  if a framebuffer object is bound that is not the default framebuffer, then only GL_NONE
    and GL_COLOR_ATTACHMENTi are accepted, otherwise a
     GL_INVALID_ENUM error is generated.
    ==> For glDrawBuffers(), Only GL_NONE, GL_FRONT_LEFT, GL_FRONT_RIGHT, GL_BACK_LEFT, and GL_BACK_RIGHT are accepted.

	
	//layout (location = 0) out vec4 normalOut; 
	// layout (location = 1) out vec4 texCoordOut; 

	// the location is related to the indexes used as an argument in glDrawBuffers.
	// For instance, considering the above example for glDrawBuffers, 
	// location 0 will output to color attachment 1, and location 1 will output to color attachment 2.

	//glBindFramebuffer(GL_FRAMEBUFFER, fbo) ;  or glBindFramebuffer(GL_FRAMEBUFFER, 0)
	// glViewport ( 0, 0, backingWidth, backingHeight ); 
*/


}


// // set a texture reference
//	void setUniformTexture(const string & name, ofBaseHasTexture& img, int textureLocation);
//	void setUniformTexture(const string & name, ofTexture& img, int textureLocation);
//	void setUniformTexture(const string & name, int textureTarget, GLint textureID, int textureLocation);
	
void testApp::initMaterials() {


	// All the file names relative to the working directory where the project file is located.
	// get active uniforms and attributes and bind indices to them.

	MaterialShader diffuseMat( "basic+diffuse", "shaders/basic-gl3.vshader", "shaders/diffuse-gl3.fshader"); 
	
	

	// copy diffuse prototype and set red color
	g_redDiffuseMat.reset(new MaterialShader( diffuseMat ));  // MaterialShader(diffuse): uses the default copy constructor

	g_redDiffuseMat->getUniforms().put("uMaterialColor", Cvec4f(1, 0, 0,1.0));
	//g_redDiffuseMat->getUniforms().put("uTexUnit0", shared_ptr<ImageTexture>( new ImageTexture("reachup.ppm", true) ) );

	// copy diffuse prototype and set blue color
	g_blueDiffuseMat.reset(new MaterialShader(diffuseMat));
	g_blueDiffuseMat->getUniforms().put("uMaterialColor", Cvec4f(0, 0, 1,1.0));

	
	// normal mapping material

	try {
	  MaterialShader normalMat ("normal+normal", "shaders/normal-gl3.vshader", "shaders/normal-gl3.fshader");
	  g_bumpFloorMat.reset(new MaterialShader (normalMat) );
	
	
	  g_bumpFloorMat->getUniforms().put("uTexColor", shared_ptr<ShaderImageTexture2D_RGB_RGB>(new ShaderImageTexture2D_RGB_RGB("sourceimages/sourceimages/Fieldstone.ppm", true)));
	 g_bumpFloorMat->getUniforms().put("uTexNormal", shared_ptr<ShaderImageTexture2D_RGB_RGB>(new ShaderImageTexture2D_RGB_RGB("sourceimages/sourceimages/FieldstoneNormal.ppm", false)));
	  g_bumpFloorMat->getUniforms().put("uMaterialColor", Cvec4f(0.5, 0.5, 0.5,1.0));

	}
	catch (const runtime_error & error ) {
		messageFile << error.what() << endl;

	}

	
	
	MaterialShader solidMat("basic+solid", "shaders/basic-gl3.vshader", "shaders/solid-gl3.fshader");

	// copy solid prototype, and set to wireframed rendering
	g_arcballMat.reset(new MaterialShader(solidMat));


	g_arcballMat->getUniforms().put("uMaterialColor", Cvec4f(0.27f, 0.82f, 0.35f, 1.0));
	g_arcballMat->getRenderStates().polygonMode(GL_FRONT_AND_BACK, GL_LINE);

	// copy solid prototype, and set to color white
	g_lightMat.reset(new MaterialShader(solidMat));

	g_lightMat->getUniforms().put("uMaterialColor", Cvec4f(1, 1, 1,1));

	// pick shader

	g_pickingMat.reset(new MaterialShader("basic+pick", "shaders/basic-gl3.vshader", "shaders/pick-gl3.fshader") );

	

};



void testApp::initPBRTObjects() {


	Options options;
	vector<string> filenames;
	// Process command-line arguments
	 for (int i = 1; i < g_argc; ++i) {
		if (!strcmp(g_argv[i], "--ncores")) options.nCores = atoi(g_argv[++i]);
		else if (!strcmp(g_argv[i], "--outfile")) options.imageFile =g_argv[++i];
		else if (!strcmp(g_argv[i], "--quick")) options.quickRender = true;
		else if (!strcmp(g_argv[i], "--quiet")) options.quiet = true;
		else if (!strcmp(g_argv[i], "--verbose")) options.verbose = true;
		else if (!strcmp(g_argv[i], "--help") || !strcmp(g_argv[i], "-h")) {
			printf("usage: pbrt [--ncores n] [--outfile filename] [--quick] [--quiet] "
				"[--verbose] [--help] <filename.pbrt> ...\n");
			return ;
		}
		else filenames.push_back(g_argv[i]);
	}
	


	// Print welcome banner
	if (!options.quiet) {
		printf("pbrt version %s of %s at %s [Detected %d core(s)]\n",
			PBRT_VERSION, __DATE__, __TIME__, NumSystemCores());
		printf("Copyright (c)1998-2012 Matt Pharr and Greg Humphreys.\n");
		printf("The source code to pbrt (but *not* the book contents) is covered by the BSD License.\n");
		printf("See the file LICENSE.txt for the conditions of the license.\n");
		fflush(stdout);
	}
	pbrtInit(options);
	// Process scene description
	PBRT_STARTED_PARSING();
	if (filenames.size() == 0) {
		// Parse scene from standard input
		ParseFile("-");
	} else {
		// Parse scene from input files
		for (int i = 0; i < filenames.size(); i++)
			if (!ParseFile(filenames[i]))
				Error("Couldn't open scene file \"%s\"", filenames[i].c_str());
	}
	//pbrtCleanup();




}  // initPBRTObjects


void testApp::initObjects() {
  int ibLen, vbLen;


// 1: create a ground
	static const float groundY = -2.0;      // y coordinate of the ground
    static const float groundSize = 50.0;   // half the ground length

  // temporary storage for a plane geometry
  getPlaneVbIbLen( vbLen, ibLen); // get the sizes of the vertex buffer and the 
                                 // index buffer for a plane. The vertex buffer
                                 // size is the number of vertices.
  vector<VertexPNTBX> vtxGround (vbLen);
  vector<unsigned short> idxGround (ibLen);

  makePlane( groundSize * 2, vtxGround.begin(), idxGround.begin() );

  // vtxGround is an array of generic vertices (position, normal, tex coord, tangent, binormal)

  
  // A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
 /* VertexPN vtxGround[4] = {
		  
    VertexPN( -groundSize, groundY, groundSize, 0, 1, 0),
    VertexPN( groundSize, groundY,  groundSize, 0, 1, 0),
    VertexPN(  groundSize, groundY, -groundSize, 0, 1, 0),
    VertexPN(  -groundSize, groundY, -groundSize, 0, 1, 0),
  };
  
  unsigned short idxGround[] = {0, 1, 2, 0, 2, 3};
  */

  

  //shared_ptr<Geometry> groundGeometry  ( new Geometry( &vtxGround[0], &idxGround[0], 4, 6,  &sqTex[0], numOfTexCoords  )  );
  shared_ptr<Geometry> groundGeometry  ( new SimpleIndexedGeometryPNTBX( &vtxGround[0], &idxGround[0], vbLen, ibLen  )  );

 
//  SimpleIndexedGeometryPNTBX => init: vbo(new FormattedVbo( Vertex::FORMAT ) ), ibo(new FormattedIbo(size2IboFmt(sizeof(Index)))) 
//   where Vertex is either VertexPN, ...:

//  const VertexFormat VertexPN::FORMAT = VertexFormat(sizeof(VertexPN))
//	  .put("aPosition", 3, GL_FLOAT, GL_FALSE, offsetof(VertexPN, p))
//	  .put("aNormal", 3, GL_FLOAT, GL_FALSE, offsetof(VertexPN, n));

 

  // 2: create a cube
 // int ibLen, vbLen;

  getCubeVbIbLen(vbLen, ibLen);

  // Temporary storage for cube geometry

  // T = tangent, X = coordinates

  vector<VertexPNTBX> vtxCube(vbLen); // vtxCube:  a vector of VertexPNTBX type with length vbLen

  vector<unsigned short> idxCube(ibLen); // idxCube: 

  // makecube<VertexPNTBX, unsigned short>( vtxCube.begin(), idxCube.begin() );

  // vtxCube and idxCube are vertices whose stroages are already
  // allocated, so that you can access each element of a vertex by means
  // of iterator.

  makeCube(5,5,5, vtxCube.begin(), idxCube.begin() ); // fill vtxCube and idxCube, starting
                                                       // from their starting pointers => a list of vertex coords,
                                                       // tex coords, normal coords, tangent coords, binormal coords

  // create vbo and ibo for vtxCube and idxCube by using SimpleIndexedGeometryPNTBX, where
  // SimpleIndexedGeometry<VertexPNTBX, unsigned short> SimpleIndexedGeometryPNTBX;

 shared_ptr<Geometry> cubeGeometry (  new SimpleIndexedGeometryPNTBX( &vtxCube[0], &idxCube[0], vbLen, ibLen ) ); // reset(T * p)

  //shared_ptr<Geometry> cubeGeometry (  new Geometry( &vtxCube[0], &idxCube[0], vbLen, ibLen, &sqTex[0], numOfTexCoords  ) ); // reset(T * p)

 
  // 4: create a sphere
  
  int slices = 100;
  int stacks = 100;

  getSphereVbIbLen(slices, stacks, vbLen, ibLen);

  // Temporary storage for cube geometry

  vector<VertexPNTBX>  vtxSphere (vbLen);

  vector<unsigned short> idxSphere (ibLen);

  float radius = 3;

  makeSphere(radius, slices, stacks,  vtxSphere.begin(), idxSphere.begin() );

  //geometry.h: typedef SimpleIndexedGeometry<VertexPNTBX, unsigned short> SimpleIndexedGeometryPNTBX;

  shared_ptr<Geometry> sphereGeometry ( new SimpleIndexedGeometryPNTBX( &vtxSphere[0], &idxSphere[0], vbLen, ibLen ) ); // reset(T * p)
 
  // or  shared_ptr<Geometry> sphereGeometry =  new Geometry( &vtx[0], &idx[0], vbLen, ibLen ); // reset(T * p)
  // NOTE: here vtx and idx arrays are copied to the vertex buffer object internally, so this data need not be global.

  // create the object frames for the created geometry

  // The following does not work, because the  matrix created by SgRbtNode::makeTranslation() is
  //	destroyed outside of this function initObjects(). To avoid it, create the matrix by the "new" method.
 //	The objects created by new methods remain unless removed explicitly. This is handled by shared_ptr<>. 
																				   
//  shared_ptr<SgRbtNode> groundRbt ( &SgRbtNode::makeTranslation( Cvec3(0,0,0) ) ); 
 // shared_ptr<SgRbtNode>  cubeRbt ( &SgRbtNode::makeTranslation( Cvec3(5.0,0,0)  ) );
  // shared_ptr< SgRbtNode>  sphereRbt ( &SgRbtNode::makeTranslation( Cvec3(0,0,0) ) );
  
  
  shared_ptr<SgRbtNode> groundRbt ( new SgRbtNode( SgRbtNode::makeTranslation( Cvec3(0,0,0) ) ) ); // this uses the non-default copy constructor ?
  shared_ptr<SgRbtNode> cubeRbt ( new SgRbtNode( SgRbtNode::makeTranslation( Cvec3(5.0,0,0)  ) )  );
  shared_ptr< SgRbtNode> sphereRbt ( new SgRbtNode ( SgRbtNode::makeTranslation( Cvec3(-0.5,0,0) ) ) );
  
  
  shared_ptr<Object> ground ( new Object( "ground", groundRbt, groundGeometry, g_bumpFloorMat) );
  shared_ptr<Object> cube ( new Object( "cube", cubeRbt, cubeGeometry, g_redDiffuseMat) );
  shared_ptr<Object> sphere ( new Object( "sphere", sphereRbt, sphereGeometry, g_blueDiffuseMat) );

  // 1: create the object id color for picking purpose

  Cvec4 pickColor;
  unsigned int id; // 32 bit unsigned


  
  id = g_objId ++ ;
  id = id << 8;  // RGBA: fill A field with 0's, shifting id 100 to the B field
  id = id | 255 ; // id for  ground; make A field to be all 1's



  pickColor = g_picker.idToColor(id);

  ground->pickColor = pickColor; // this stored pickColor will be sent to the shader Uniform variable when drawing 


  // store the objects in the object list

  g_objectList.push_back( ground );

  // add the object and the id to the list for picking
  g_picker.addToMap( id, ground );



  // 2: the same for the cube
 
  id= g_objId ++;
  id = id << 8;  
  id = id | 255; // id for  cube 
  pickColor  = g_picker.idToColor(id );

  cube->pickColor = pickColor;
 

  g_objectList.push_back( cube );


  g_picker.addToMap( id, cube );

   
  
  
  // 4: the same for sphere
  id = g_objId ++ ;
  id = id << 8;
  id = id | 255; // id for  sphere 

  pickColor = g_picker.idToColor(id);

  sphere->pickColor  = pickColor;

  g_objectList.push_back( sphere );

  g_picker.addToMap( id, sphere );

  
 

};

 void testApp::initRainbowAABB() {
  int ibLen, vbLen;


  getCubeVbIbLen(vbLen, ibLen);

  // Temporary storage for cube geometry

  // T = tangent, X = coordinates

  vector<VertexPNTBX> vtxAABB(vbLen); // vtxCube:  a vector of VertexPNTBX type with length vbLen
                                      // what happens if I use VertexPNX ?

  vector<unsigned short> idxAABB(ibLen); // idxCube: 

 //Cvec3 AABBcenter = ( g_AABBmax + g_AABBmin ) / 2.0;

 // create a cube whose min coordinates are -SideLengths/2 and whose max coordinates are SideLengths/2,
 // relative the object coordinate system whose center passes through the center of AABB;
 // The object coordinate system is moved by AABBcenter from the world origin:

 // (t1, t1, t3, 1) = g_modelMatrixAABB * (0,0,0,1) where t = AABBcenter

 //g_modelMatrixAABB = Matrix4::makeTranslation( AABBcenter );

 //Matrix4 invModelMatrixAABB = inv( g_modelMatrixAABB);

//g_AABBmin = Cvec3( invModelMatrixAABB * Cvec4(g_AABBmin,1) );
//g_localAABBmax = Cvec3( invModelMatrixAABB * Cvec4(g_AABBmax,1) );

 //Cvec3 SideLengths = g_AABBmax - g_AABBmin;

 //messageFile << " SideLengths: "  << SideLengths  << endl;
 //messageFile << "g_AABBcenter="  << AABBcenter << endl;
 //messageFile << "localAABBmin =" << g_localAABBmin << endl;
 //messageFile << "localAABBmax =" << g_localAABBmax << endl;

//(std::ofstream &) ( messageFile << " SideLengths: " ) << SideLengths  << endl;
//(std::ofstream &) ( messageFile << "AABBcenter=" ) << AABBcenter << endl;
// why do I need this casting? (1) When messageFile << string is executed, the type of messageFile
//  ofstream is matched with its superclass "ostream"  which defines the << operator with respect to
// all kinds of input parameters. (2) When the second parameter Vector is encountered, the compiler
// tries to check if the superclass "ostream" has the << operator with Vector input. But it does not.
// When you cast (messageFile << "SideLengths:"), you don't have this issue, because std::ofstream has
// << operator with Vector parameter. The better way is to use std::ostream when you overload <<.


//messageFile <<  "g_modelMatrixAABB:"  << g_modelMatrixAABB << endl;
//messageFile << "invModelMatrixAABB:" << invModelMatrixAABB << endl;

 Cvec3 AABBsize = g_initAABBmax - g_initAABBmin;


 makeCube( AABBsize[0], AABBsize[1], AABBsize[2], vtxAABB.begin(), idxAABB.begin() ); // fill vtxAABB  starting
                                                       // from the  starting pointer with a GenericVertex which consists of  vertex coords,
                                                       // tex coords, normal coords, tangent coords, binormal coords

  // create vbo and ibo for vtxCube and idxCube by using SimpleIndexedGeometryPNTBX, where
  // SimpleIndexedGeometry<VertexPNTBX, unsigned short> SimpleIndexedGeometryPNTBX;

  // Upload the geometry data and bind it to the Vbo 
 shared_ptr<Geometry> AABBGeometry (  new SimpleIndexedGeometryPNTBX( &vtxAABB[0], &idxAABB[0], vbLen, ibLen ) ); 
 


  // rainbow shader

  MaterialShader rainbowMat ("basic+rainbow", "shaders/basic-gl3.vshader", "shaders/rainbow-gl3.fshader");


  g_rainbowMat.reset(new MaterialShader( rainbowMat) );

  g_rainbowMat->getUniforms().put("uScatTex", shared_ptr<ShaderImageTexture2D_RF_RF>(new ShaderImageTexture2D_RF_RF("rainbow/scattering2_120_fort.txt",  false)));

  g_rainbowMat->getUniforms().put("uPhaseTex", shared_ptr<ShaderImageTexture3D_RF_RF>(new ShaderImageTexture3D_RF_RF("rainbow/phaseFunction2_120_100_fort.txt",  false)));

  shared_ptr<Object> AABB ( new Object( "AABBrainbow",  g_AABBRbt, AABBGeometry, g_rainbowMat) );
  
  //shared_ptr<Object> AABB ( new Object( "AABBdiffuse",  g_AABBRbt, AABBGeometry, g_redDiffuseMat) );

  //shared_ptr<Object> cube ( new Object( cubeRbt, cubeGeometry, g_bumpFloorMat) );
  //shared_ptr<Object> sphere ( new Object( sphereRbt, sphereGeometry, g_bumpFloorMat) );

  // 1: create the object id color for picking purpose

  Cvec4 pickColor;
  unsigned int id; // 32 bit unsigned


  // 3: the same for the AABB
  id= g_objId++;
  id = id << 8;  
  id = id | 255; // id for  cube 
  pickColor  = g_picker.idToColor(id );

  AABB->pickColor = pickColor;


  g_objectList.push_back( AABB );


  g_picker.addToMap( id, AABB );

  
 }
 

//--------------------------------------------------------------
void testApp::windowResized(int w, int h){  // call back function for event processing
 
	// instance->windowW = w; => the new window size is already set before calling windowResized()

	// instance->windowH = h;
  g_windowHeight = h;
  g_windowWidth = w;
 
  
  glViewport(0, 0, g_windowWidth, g_windowHeight);



  cout   << "Size of window is now " << g_windowWidth << "x" << g_windowHeight << endl;
  cout << "results of ofGetHeight and ofGetWidth= " << ofGetWidth() << "x" << ofGetHeight() << endl;

  //cerr  << "Size of window is now " << g_windowWidth << "x" << g_windowHeight << endl;
  //system ("PAUSE");

  updateFrustFovY();
 
  // At every frame, display() callback is called to redraw the
	// scene. So, windowResized event handler needs not do anything to redraw the scene.

   reDrawWindowEvent = true;
	
}


//--------------------------------------------------------------
void testApp::mouseDragged(int x, int y, int button){ // call back function for event processing

 // This callback is called every time the pressed mouse is moved a little.
// x, y: the current mouse position after this little  motion is finished.


  /*cout  << "dragged mouse button " << button << "X = " << x << "Y = " << y << endl;
  // Here the moved button is 0, meaning that which button is pressed while moving the mouse is ignored.
 
  cout <<  "Control Key Pressed =" << ofGetKeyPressed( OF_KEY_CONTROL ) << endl;
  cout <<  "Shift  Key Pressed =" << ofGetKeyPressed( OF_KEY_SHIFT ) << endl;
  cout <<  "ALT Key Pressed =" << ofGetKeyPressed( OF_KEY_ALT ) << endl;
  */


  const double dx = x - g_prev_mouseClickX;
  const double dy = g_windowHeight - y - 1 - g_prev_mouseClickY;
  
  Matrix4 m;
  
  g_prev_mouseClickX = x;
  g_prev_mouseClickY = g_windowHeight - y - 1;

  if ( ofGetKeyPressed( OF_KEY_CONTROL ) && !ofGetKeyPressed( OF_KEY_ALT ) )  {
	  // translate the camera 
      switch ( button ) {
		  case OF_MOUSE_BUTTON_LEFT:
			 // cout  << " left button moves" << endl;
			  m = Matrix4::makeTranslation( Cvec3( dx, dy, 0.0) * 0.01 ); 
			  //m = Matrix4::makeTranslation( Cvec3( dx, 0.0, 0.0) * 0.01 ); 
			  break;
		  case OF_MOUSE_BUTTON_RIGHT:
			 // cout  << " right button moves" << endl;
			  m = Matrix4::makeTranslation( Cvec3( 0.0, 0.0, -dy) * 0.01 ); 
			  break;

		  default: 
			  m = Matrix4::makeTranslation( Cvec3( 0.0, 0.0, 0.0)  ); 
	  } // switch

	   g_eyeRbt *= m;

	   // window = ofPtr<ofAppBaseWindow>(new ofAppGLFWWindow());
	   //ofPtr<ofAppGLFWWindow>	appwindow =  window;
	   //appwindow  ->display();

	   reDrawWindowEvent = true;

  } // if (! g_rotation)

  else if ( ofGetKeyPressed( OF_KEY_CONTROL) && ofGetKeyPressed( OF_KEY_ALT) ) { //  rotate the camera

	  switch ( button ) {
		  case OF_MOUSE_BUTTON_LEFT:

			  //inline static ofMatrix4x4 newRotationMatrix( float angle, const ofVec3f& axis);
			  m = Matrix4::makeXRotation(-dy ) * Matrix4::makeZRotation( dx ); // // push the sphere on the top
			  // m = Matrix4::makeXRotation( dx); 
			  break;
		  case OF_MOUSE_BUTTON_RIGHT:
			  m = Matrix4::makeYRotation( dy ); // push the sphere on the side
			  break;

		  default: 
			  m = Matrix4::makeTranslation( Cvec3( 0.0, 0.0, 0.0)  ); 
	  } // switch

	  g_eyeRbt *= m;

	  //ofPtr<ofAppBaseWindow> 		window;  // changed by Moon Jung, 2014/4/23
	  
	   // window = ofPtr<ofAppBaseWindow>(new ofAppGLFWWindow());
	   //ofPtr<ofAppGLFWWindow>	appwindow =  window;
	   // 	  appwindow  ->display();

	  reDrawWindowEvent = true;
  } // else if
  
  

  else if ( g_picked_mode )  { // An object has been picked by pressing the left or right mouse

   if (   ofGetKeyPressed( OF_KEY_SHIFT ) && !ofGetKeyPressed( OF_KEY_ALT) )   { // translate the selected object
	 switch ( button ) {
	  case OF_MOUSE_BUTTON_LEFT:
		   //cout  << " left button moves" << endl;
		   m = Matrix4::makeTranslation( Cvec3( dx, dy, 0.0) * 0.01 ); 
		   //m = Matrix4::makeTranslation( Cvec3( dx, 0.0, 0.0) * 0.01 ); 
		   break;
	  case OF_MOUSE_BUTTON_RIGHT:
		  //cout  << " right button moves" << endl;
		   m = Matrix4::makeTranslation( Cvec3( 0.0, 0.0, -dy) * 0.01 ); 
		   break;

	  default: 
		   m = Matrix4::makeTranslation( Cvec3( 0.0, 0.0, 0.0)  ); 
	  } // switch

	  assert( ("g_currentPickedRbtNode should not be Null", g_currentPickedObject != nullptr) );

	  *(g_currentPickedObject->objectRbt) = m * *(g_currentPickedObject-> objectRbt) ;

    } // if ( translate )

    else if ( ofGetKeyPressed( OF_KEY_SHIFT) && ofGetKeyPressed( OF_KEY_ALT) ) { //  rotate the selected object
	   
		  switch ( button ) {
	       case OF_MOUSE_BUTTON_LEFT:

			//inline static ofMatrix4x4 newRotationMatrix( float angle, const ofVec3f& axis);
		     m = Matrix4::makeXRotation(-dy ) * Matrix4::makeZRotation( dx ); // push the sphere on the top
		  // m = Matrix4::makeXRotation( dx); 
		     break;
	       case OF_MOUSE_BUTTON_RIGHT:
		     m = Matrix4::makeYRotation( dy ); // push the sphere on the side.
		     break;

	       default: 
		     m = Matrix4::makeTranslation( Cvec3( 0.0, 0.0, 0.0)  ); 
	      } // switch

		  assert( ("g_currentPickedRbtNode should not be Null", g_currentPickedObject != nullptr) );

		  *(g_currentPickedObject->objectRbt) = m  *  *(g_currentPickedObject-> objectRbt) ;

       } // else if

  } // (g_picked_mode)
  
 
} // mouseDragged()





//--------------------------------------------------------------
void testApp::mouseReleased(int x, int y, int button){


} // mouseReleased()




//--------------------------------------------------------------
void testApp::mouseMoved(int x, int y ){

	// This callback is called every time the pressed mouse is moved a little.
	// x, y: the current mouse position after this little  motion is finished.


	
}


//--------------------------------------------------------------
void testApp::mousePressed(int x, int y, int button){
	
  //cout << "mouse x,y = " << x << "," << y << endl;

  g_prev_mouseClickX = x;
  g_prev_mouseClickY = g_windowHeight - y - 1;  // conversion from GLUT window-coordinate-system to OpenGL window-coordinate-system
  //g_prev_mouseClickY = y;

  /*
  cout  << "pressed mouse =" <<  button << endl;
  cout  << "shift  modifier key pressed? = " << ofGetKeyPressed(OF_KEY_SHIFT) << endl;
  */

  g_pressed_button = button; // the pressed button has meaning for moving the picked object or the camera later

  
  if ( ofGetKeyPressed( OF_KEY_SHIFT ) ) { // shift +  any mouse press leads to pickin an object
    cout  << "draw for picking " << endl;
	// for temp debugginhg
	
	drawForPicking();      //    In the picking mode, you render the simplified version of the scene to the back buffer,
		                       //     while maintaining the scene of the screen intact, because it is in the front buffer.
		                      
	   //  Read pixel operation needed for picking is performed with respect to the back buffer, not to the front buffer

	// for temp debugging

	g_currentPickedObject = g_picker.getRbtNodeAtXY( g_prev_mouseClickX, g_prev_mouseClickY );

	
	 // g_currentPickedRbtNode points to one of g_objectRbt[], which has been picked.

     // g_currentPickedRbtNode will be NULL, if no object has been actually picked.

   if ( g_currentPickedObject == nullptr ) {

		  cout << "no object has been picked" << endl;
		  g_picked_mode = false;

	 }
   else {
          cout << "an object has been picked" << endl;

		  g_picked_mode = true;
     }

  } // g_picking_mode

  // after drawForPicking(), normal draw is called, which is in fact called every frame.

} // mousePressEvent()




//--------------------------------------------------------------
void testApp::keyPressed  (int key){ 
//	OF_KEY_LEFT_CONTROL,  OF_KEY_LEFT_SHIFT,
// OF_KEY_LEFT_ALT,
//	OF_MOUSE_BUTTON_LEFT
 /* 
  cout << "pressed key=" << key<< endl;
  cout << "OF_KEY_SHIFT= " << OF_KEY_SHIFT << endl;
  cout << "OF_KEY_CONTROL= " << OF_KEY_CONTROL << endl;
  cout << "OF_KEY_ALT= " << OF_KEY_ALT << endl;
*/

  switch (key) {

  case OF_KEY_SHIFT:  // debug mode, render to texture and print the content
	  renderToFbo();
	  break;

  case 'r':  // debug mode, read pixels
	  readPixels();
	  break;

  case OF_KEY_ESC: 
    exitApp();                                  // ESC 
  case OF_KEY_CONTROL: 
    g_camera_mode = true;
	break;

  case 'h': 
    cout << " ============== H E L P ==============\n\n"
    << "h\t\thelp menu\n"
    << "s\t\tsave screenshot\n"
    << "f\t\tToggle flat shading on/off.\n"
    << "ESC\t\t exit\n"
    << "Cntrl\t\t  camera moving mode\n"
	<< "Alt\t\t rotate the selected object or the camera\n"
    << "Shift\t\t  picking mode\n" << endl;
    break;

  case 's': 
    glFlush();
    writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
	break;


  default: break;

 } // switch
	
} // keyPressed()

//--------------------------------------------------------------
void testApp::keyReleased(int key){ 


  //cout   << "released key=" << key << endl;
  switch (key) {

   case OF_KEY_CONTROL: 
	  g_camera_mode = false;
	  break;
  
   case OF_KEY_SHIFT: 
	   g_picked_mode  = false;
	   break;

   default: 
	   break;
  }

}



//--------------------------------------------------------------
void testApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void testApp::dragEvent(ofDragInfo dragInfo){ 

}

