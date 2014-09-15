#include "ofApp.h"


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
#include "glut.h"

// for rainbows
#include "spa.h"  //include the SPA header file
//#include "optical_depth.h"
#include "irradiance_sun.h"
//#include "phaseFunction.h"
#include "spectrum_to_rgb.h"

#include "phaseFunction.h"

#include <ctype.h>  // from pbrt for readFloatFile
#include <stdlib.h>  // from pbrt for readFloatFile


// Where is shared_ptr? 
// 1.If your C++ implementation supports C++11 (or at least the C++11 shared_ptr),
//	then std::shared_ptr will be defined in <memory>.


// 2.If your C++ implementation supports the C++ TR1 library extensions,
// then std::tr1::shared_ptr will likely be in <memory> (Microsoft Visual C++)
//
// or <tr1/memory> (g++'s libstdc++). 
// Boost also provides a TR1 implementation that you can use.
// 3.Otherwise, you can obtain the Boost libraries and use boost::shared_ptr, 
// which can be found in <boost/shared_ptr.hpp>.


#include <memory> // for shared_ptr
using namespace std; // string type is defined in namespace 
using namespace std::tr1; // for shared_ptr

/*
std::shared_ptr is a smart pointer that retains shared ownership of an object through a pointer. Several shared_ptr objects may own the same object. The object is destroyed and its memory deallocated when either of the following happens: 
 the last remaining shared_ptr owning the object is destroyed. 
 the last remaining shared_ptr owning the object is assigned another pointer via operator= or reset(). 

The object is destroyed using delete-expression or a custom deleter that is supplied to shared_ptr during construction. 

A shared_ptr can share ownership of an object while storing a pointer to another object. This feature can be used to point to member objects while owning the object they belong to. 

A shared_ptr may also own no objects, in which case it is called empty. 

*/

#define Assert(expr) \
    ((expr) ? (void)0 : \
        fprintf(stderr, "Assertion \"%s\" failed in %s, line %d", \
               #expr, __FILE__, __LINE__))



/* for debugging: exclude pbrt for the time being, by Moon Jung, 2014/8/7
#include "stdafx.h"
#include "api.h"
#include "probes.h"
#include "parser.h"
#include "parallel.h"

bool ParseFile(const string &filename);
extern int g_argc;
extern char **g_argv;

*/

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
							g_rainbowMat,
							g_drawTextureMat;

shared_ptr<MaterialShader> g_overridingMaterial;

shared_ptr<Object> g_drawTextureObj; // used  simply to draw a texture

//Just avoid any using namespace directives in the headers. 
//	(And use with care in C++ files, if at all. Inside functions is preferred.)
// as long as they  * appear after all #includes.

using namespace std;      // for string, vector, iostream, shared_ptr and other standard C++ stuff
using namespace std::tr1; // for shared_ptr<T>

bool g_redrawWindowEvent = true; // used in ofAppGLFWWindow.cpp

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
double  g_frustMinFov,  g_frustFovY;
double  g_frustNear, g_frustFar;

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

GLuint g_fbo, g_colorTex, g_depthTex; 

bool g_debugMode = true;

// for rainbow drawing
// The coordinates of the AABB box relative to the  world coordinate system

//float g_phaseFunction[ nRadii * nSpectralSamples * nThetas ];

// THis includes the hemisphere for sky
//static const Cvec3 g_initAABBmax(90.3325, 75.8859, 90.3325);
//static  const Cvec3 g_initAABBmin(-90.3325, -14.4466, -90.3325);

// THis includes the hemisphere for sky

// for National Modern Musium
//static const float rainbowVolumeHeight = 21.3; 
//static const float rainbowVolumeWidth = 30;
//static const float rainbowVolumeDepth = 5;

// for Asian Game 
static const float rainbowVolumeHeight = 30; 
static const float rainbowVolumeWidth = 40;
static const float rainbowVolumeDepth = 10;

// The coordinates of the AABB, which will be used to compute the sizes of the view volume

static      	Cvec3 g_AABBmax, g_AABBmin; 
static const   Cvec3 g_AABBsize (rainbowVolumeWidth, rainbowVolumeHeight, rainbowVolumeDepth  );
static           Cvec3 g_initEyeLocation;
static           Cvec3 g_AABBcenter;

shared_ptr<SgRbtNode> g_AABBRbt, g_quadRbt;
SgRbtNode g_SceneRbt;  

// --------- Scene
// --------- Scene

Cvec4 g_light1Pos, g_light2Pos;  // define two lights positions in world space, which is set right in front of the water
                                 // volume

Cvec4 g_light1Color, g_light2Color; 

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

vector < shared_ptr<Object> >  g_objectPtrList;
shared_ptr< Object> g_rainbowObj; // default constructor is called to create an empty shared_ptr object

shared_ptr< Geometry > g_geometry;

std::ofstream messageFile("./messageFile.txt", std::ios::out);
std::ofstream textureFile("./textureDebug.txt", std::ios::out);

//--------------------------------------------------------------
ofApp:: ofApp() // default constructor

{
}






void ofApp::setup(){
	// drop down menu settup

	
	ofBackground(150);

	
    gui = new ofxUICanvas(300,300); // this->gui variable exists

    gui->addLabel("DROPDOWN MENU", OFX_UI_FONT_MEDIUM);
    gui->addSpacer();
    gui->addLabel("'1' TO ADD TO LIST", OFX_UI_FONT_SMALL);
    gui->addLabel("'2' TO DELETE FROM LIST", OFX_UI_FONT_SMALL);
    gui->addLabel("'3' TO DELETE ALL IN LIST", OFX_UI_FONT_SMALL);
    gui->addSpacer();
    vector<string> names;
    names.push_back("render background to system framebuffer");   
	names.push_back("dump the system framebuffer to file");   
	names.push_back("render background to FBO");
	names.push_back("dump FBO to file");  
	names.push_back("render background to system framebuffer using FBO texture");
	names.push_back("render background from FBO file");
	names.push_back("render background from  system framebuffer file");
	names.push_back("render rainbow using FBO background");

    gui->setWidgetFontSize(OFX_UI_FONT_SMALL);
    gui->addToggle("SHOW ACTIVE", false);
    ddl = gui->addDropDownList("DROPDOWN", names); // this->ddl variable exists

   // ddl->setAllowMultiple(true);
//    ddl->setAutoClose(true);
    gui->autoSizeToFitWidgets(); 
//    gui->setDrawWidgetPadding(true);

    ofAddListener(gui->newGUIEvent, this, &ofApp::guiEvent);

	



	// set up the sun light intensities and the sun direction
//	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ALPHA );
	setupSun();

	
	// set up the eye position and direction relative to the g_localAABB which is created at
	// initObjects(). Then create the eye Matrix

	// For now, It is assumed that localAABB is centered at the world origin
	// height / horizontalDistFromEyeToAABB  = tan (g_frustFovY / 2.0) 


	setupCamera();

	bool g_Gl2Compatible = false; // use shader programs of type "*-gl3 *"

	//Sets the background clearing function to be auto (default) or not. 
	//	If non-auto, then background clearing will not occur per frame (at the start of draw)
	//but rather, whenever ofBackground is called.
	
    //ofSetBackgroundAuto( false); // It will make a single buffering by drawing to the front buffer
	
	//ofBackground(128,200,255); // sky color

	ofSetBackgroundAuto( true); // which is default

	//------------------------------------------------------------
    // void ofAppGLFWWindow::disableSetupScreen(){
	//   bEnableSetupScreen = false;
    // };

	//_description: _

    // Every update/draw cycle, the function ofSetupScreen is called. 
	// That function sets the perspective, coordinate system, and some other openGL parameters.
	// If you need to use your own parameters, the call to that function can be disabled with ofDisableSetupScreen.

	//ofDisableSetupScreen(); SetupScreen() is needed to draw the widgets. But the camera parameters set by
	// this function will be ignored by draw() function defined in this class.

	g_windowWidth = (float) ofGetWidth(); 
	g_windowHeight = (float) ofGetHeight(); 

	
	//g_windowWidth =  nThetas;  
	//g_windowHeight = nSpectralSamples;  

		
	initGLState();

    initMaterials();
	

	//initPBRTObjects();

	initObjects(); // initObjects() is for testing. use it or initPBRTObjects()
	//initAABB();
	//initDrawTextureMat();

   // initRainbow();



}

void ofApp::exit() {


}


void ofApp::setupSun() {

spa_data spa;  //declare the SPA structure
    int result;
    Cvec3 sunRay;
    
    // Default values to calculate sun's position
    spa.year          = 2014;
    spa.month         = 11;
    spa.day           = 1;
    spa.hour          = 16;
    spa.minute        = 30;
    spa.second        = 00;
    spa.timezone      = 9.0;             // KST (GMT + 9): Korean or Japan time, 9 hours ahead
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
	cout  << "lightRGBColor =" << g_lightRGBColor << endl;

    get_time_and_location(spa.year, spa.month, spa.day, spa.hour, spa.minute, spa.timezone, spa.longitude, spa.latitude);
    
    // calculate sun's position
    result = spa_calculate(&spa);   // input is valid when result eq 0
    
	if ( result !=0 ) {
		messageFile   << "error in the calculation of sun position" << endl;
		cout   << "error in the calculation of sun position" << endl;
		return; 
	}


    // calculate sunRay

    calculate_sunRay(sunRay, spa.zenith, spa.azimuth); // vector sunRay points to the sun
	
	messageFile  << "sunRay zenith angle = " << spa.zenith <<" sunRay azimuth" << spa.azimuth << endl;
	
	messageFile  << "sunRay in geocentric coord (world coord system) = (" << sunRay  << endl;

	cout  << "sunRay zenith angle = " << spa.zenith <<" sunRay azimuth" << spa.azimuth << endl;
	
	cout  << "sunRay in geocentric coord (world coord system) = (" << sunRay  << endl;
	g_sunRayDir = Cvec3(-sunRay[1], sunRay[2], -sunRay[0]); // rename the axes to 3D graphics convention : z => y, x => -z, y => -x => y=z, z = -x, x=-y
	
	// E.g: In the original coord system: (-10, 20, 5) [ azimth= south-west, polar = positive] =? (-20, 5, 10)
	
    
    messageFile  << "sunRay in Graphics coord = (" << g_sunRayDir << endl;
	cout   << "sunRay in Graphics coord = (" << g_sunRayDir << endl;

	
}

void ofApp::setupCamera() {

	
	
// set the parameters for the camera, which define the internal working of the camera.
    g_frustMinFov = 60.0;  // A minimal of 100 degree field of view for rainbow viewing

 //   g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)
//   It means there is very high precision at the near plane, but very little precision at the far plane. 
//	If the range [-n, -f] is getting larger, it causes a depth precision problem (z-fighting); 
//  a small change of ze around the far plane does not affect on zn value. 
//  The distance between n and f should be short as possible to minimize the depth buffer precision problem.

    g_frustNear = -10.1;    // near plane
    g_frustFar = -100.0;    // far plane


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

	// set up the reference frame of the scene: In this case, the scene is located on the east side of the 
	// geocentric reference frame (the z axis = the north, y = the gravity direction, x = the east side
	
	g_SceneRbt = Matrix4::makeAxisRotation( - PI / 2, Cvec3(0,1,0) ); //for NationalModernArt
	
	//g_SceneRbt = Matrix4::makeAxisRotation( 0.0, Cvec3(0,1,0) ); // for Incheon Asian Game

    // set the camera rotation so that it points to the scene to the east of the geocentric frame
	// Ignore the above calculation

	g_cameraRotMat = g_SceneRbt;  

	double eyeHeight = 1.7;
	
	g_invCameraRotMat = inv( g_cameraRotMat);
	
   
	float distanceToEye =  42 ;  // the observer on the ground in Asian Game

	
	Cvec3 eyeLocation = Cvec3(0, eyeHeight, distanceToEye);

	// the position of the eye relative to the world coordinate system.


	g_eyeRbt =   g_cameraRotMat * Matrix4::makeTranslation( eyeLocation  );
	

	
	g_invCameraRotMat = inv( g_cameraRotMat);
	messageFile << "g_invRotMat: \n" << g_invCameraRotMat << endl;

	Matrix4 g_invEyeRbt = inv( g_eyeRbt);


	// setup light sources for Asian Game
	float  distanceToLight =  42 + 24 + 5 + 300;
	float  heightToLight = 40 + 15;
	float  separationLight = 5;

	g_light1Pos = g_SceneRbt * Cvec4( separationLight/2.0, heightToLight, distanceToLight, 1.0);
	g_light2Pos = g_SceneRbt * Cvec4(-separationLight/2.0, heightToLight, distanceToLight, 1.0);  // define two lights positions in world space

    g_light1Color =  Cvec4(1.0, 1.0, 1.0, 1.0);
    g_light2Color =  Cvec4(1.0, 1.0, 1.0, 1.0);// set up light sources


	
}



Cvec3  ofApp::getSunLightRGBColor() {
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
	
	for (lambda = lambdaStart, j=0; j < nLambdas; lambda += lambdaStep, j++ ) {
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




Matrix4 ofApp::makeProjectionMatrix() {

	return Matrix4::makeProjection(
		g_frustFovY, g_windowWidth / static_cast <double> (g_windowHeight),
		g_frustNear, g_frustFar);
}


// update g_frustFovY from g_frustMinFov, g_windowWidth, and g_windowHeight
void ofApp::updateFrustFovY() {

  if (g_windowWidth >= g_windowHeight)

    g_frustFovY = g_frustMinFov;

  else {
    const double RAD_PER_DEG = 0.5 * 3.14159 /180;

    g_frustFovY = atan2(sin(g_frustMinFov * RAD_PER_DEG) * g_windowHeight / g_windowWidth, cos(g_frustMinFov * RAD_PER_DEG)) / RAD_PER_DEG;
  }
}

//--------------------------------------------------------------
void ofApp::update(){

}

static bool backgroundBackedUp = false;

void ofApp::draw() {}


void ofApp::myOwnDraw(bool *isFrameBufferRedrawn)  { // callback function for draw event: 
	// This callback function is registered at the system initialization, 
	// so called earlier than the other draw callbacks

	// for debugging
	//return;


	//if ( !g_redrawWindowEvent ) {
	//	*isFrameBufferRedrawn = false;
	//	return;
	//}

	
	cout << "I am here in myOwnDraw() to draw" << endl;
	//renderToFBO();
	renderToSysBuffer();



	//drawFBOTexture();

	try {
	  checkGlErrors();
	}
	catch (const runtime_error & error ) {
      std::cout << error.what() << endl;
	  messageFile << error.what() << endl;

	}

	//Different GL implementations buffer commands in several different locations,
	//including network buffers and the graphics accelerator itself. glFlush empties 
	//all of these buffers, causing all issued commands to be executed as quickly
	//as they are accepted by the actual rendering engine. Though this execution may not be
	//completed in any particular time period, it does complete in finite time


	//If you dot glFlush(), the current frame written to the buffer may not be immediately drawn to the screen.

	glFlush();

	// for debugging

	*isFrameBufferRedrawn = true;

	g_redrawWindowEvent = false;

	//glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)
	// This action will be done in ofAppGLFWWindow.cpp after ofNotifyDraw() which calls
	// ofApp::draw()

} // draw()

void ofApp::drawForPicking() {
 
   
  drawPseudoColors();
  
  checkGlErrors();

  glFinish();
  checkGlErrors();
}

void ofApp::drawPseudoColors() {

  glEnable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);

  glClearColor(0, 0, 0, 0); // this background color is only meant to be used
                               // for the background of the back buffer for rendering in picking.
   
  glBindFramebuffer(GL_FRAMEBUFFER, 0); // 0 = the system-provided buffer
  

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // CLEAR THE COLOR AND DEPTH BUFFER => does it apply to FBOs too?


// Declare an empty uniforms

  Uniforms extraUniforms;
  // build & send proj. matrix to vshader

  
  const Matrix4 projMatrix = makeProjectionMatrix();


  extraUniforms.put("uProjectionMatrix", projMatrix);
   
  g_invEyeRbt = inv(g_eyeRbt);

  cout << "current eye frame:" << "\n" << g_eyeRbt << endl;

  

  extraUniforms.put("uWindowWidth", g_windowWidth );
  extraUniforms.put("uWindowHeight",g_windowHeight);
	
  
  g_overridingMaterial = g_pickingMat;

  for ( int i  = 0; i <g_objectPtrList.size() - 1 ; i++ ) { // (5) This loop goes over the list of objects g_objectPtrList and draw each object in the list.
	                                                //     The last object is a rainbow, which is not to be picked up. So skip over it
		                                            //     during picking.

     Matrix4 MVM = g_invEyeRbt * (  *(g_objectPtrList[i]->objectRbt) ) ; // g_currentPickedObject->objectRbt may have been
	                                                              // changed by picking, object includes AABB
    //Matrix4 NMVM = normalMatrix(MVM);
  
	 extraUniforms.put("uModelViewMatrix", MVM).put("uNormalMatrix", normalMatrix(MVM) );

	 
    Cvec4 pickColor = g_objectPtrList[i]->pickColor; // sRGB color is written to the framebuffer; linear RGB when read from the buffer


	//Cvec4 materialColor = g_objectPtrList[i]->objectColor;
	
	//std::cout << "Linear [float] Color for the Object to be Drawn For Picking=" << pickColor[0] << ","<< pickColor[1] << "," << 
	//            pickColor[2]  << "," << pickColor[3] << endl;

    extraUniforms.put( "uMaterialColor", pickColor ); // set material color. In the case of
	                                                 // g_overridingMaterial, there are no uniform "uMaterialColor" assigned when
	                                                  // when it is created.


    checkGlErrors();



    // the following draw() method will use g_overridingMaterial->draw() if g_overridingMaterial is on

	 g_objectPtrList[i]->draw( extraUniforms );  // the material Color is also a uniform value, but it is
	                                         // part of Material of this object. It will be taken care of 
	                                         // within the draw method.

     checkGlErrors();
  }  // for

	// unset the overriding material
  g_overridingMaterial.reset();
	
  
} // drawPseudoColor() for  picking

void ofApp::setupFBO() {

	// framebuffer object with value zero is reserved to represent the default framebuffer 
    //  provided by the windowing system.
    // FBO uses a nonzero framebuffer object fbo:
    // GLuint fbo;
	 glGenFramebuffers(1, &g_fbo);
	 
	 try {
	  glBindFramebuffer(GL_FRAMEBUFFER, g_fbo);  
	  checkGlErrors();
	 }
	 
	  catch ( const runtime_error & error ) {
		 //std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	 }

    /* glClear on FBO should be called after it is completely specified. Move to later. See below

    try {

      glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	  checkGlErrors();
	  // GL_INVALID_FRAMEBUFFER_OPERATION is  caused by calling glClear when there's no current draw framebuffer
	  // (as will be the case for a windowless OpenGL context before you bind an FBO, for example).


	 }
	 
	  catch ( const runtime_error & error ) {
		 //std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	 }
	 */
	  
	 try {
	  glViewport(0,0,g_windowWidth, g_windowHeight);   
	  checkGlErrors();
	 }
	 
	  catch ( const runtime_error & error ) {
		 //std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	 }

				
   
	// The texture we're going to render colors to 

	 //GLuint colorTex;
	 glGenTextures(1, &g_colorTex);

	// "Bind" the generated texture as the current texture: 
	// all future texture functions will modify this current texture 
	 glBindTexture(GL_TEXTURE_2D, g_colorTex);
	
	 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	 glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	 glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	 glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
 
	 glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, g_windowWidth, g_windowHeight, 0, GL_RGBA,
		GL_UNSIGNED_BYTE, NULL);

	// attach the texture to the framebuffer
	 glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, g_colorTex, 0); // 0 = mipmap level  

	 glBindTexture(GL_TEXTURE_2D,0);
	 // depth texture
	 	
	//GLuint depthTex;
	glGenTextures(1, &g_depthTex);
	glBindTexture(GL_TEXTURE_2D, g_depthTex);
		
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, g_windowWidth, g_windowHeight, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
	
	// attach the depth texture to the framebuffer
	glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, g_depthTex, 0);
	glBindTexture(GL_TEXTURE_2D,0); 

	//GLuint attachments[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
	GLuint attachments[1] = { GL_COLOR_ATTACHMENT0 };
	glDrawBuffers(1, attachments);		// drawing color and depth all together? - need to check
	//DrawBuffer is per-framebuffer state, and will default to COLOR_ATTACHMENT0 for non-zero FBOs.
	 
	// the buffer selection in glDrawBuffers is part of the framebuffer object state. 
	// Therefore, if this setting is constant, this function can be called only once when creating the framebuffer.
	// glDrawBuffers defines an array of buffers into which outputs from the fragment shader data will be written
	//  If a fragment shader writes a value to one or more user defined output variables, then the value of each variable
	// will be written into the buffer specified at a location within bufs
	// bufs \{ GL_FRONT_LEFT, ... GL_COLOR_ATTACHMENTn}

	

	GLenum e = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (e != GL_FRAMEBUFFER_COMPLETE)
		cout << "There is a problem with the FBO\n" << endl;

	 //Now that FBO is completely specified, clear the buffer with  the original clear color in order to draw the real scene
      
	//You should use glClear (...), modern GPUs use color buffer, depth buffer and stencil
	 //buffer compression. Clearing the buffer is extremely cheap because these buffers are 
	 //often hierarchically tiled, and the process of clearing the buffer amounts to flipping 
	 //one or two bits in each tile. It even helps with many early fragment tests and general
	 //frame buffer throughput if you regularly clear the buffers. On Tile-Based Deferred 
	 //Rendering GPUs (e.g. PowerVR SGX - all iOS devices) it is equally important for similar reasons. 

	 
	 try {

      glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	  checkGlErrors();
	  // GL_INVALID_FRAMEBUFFER_OPERATION is  caused by calling glClear when there's no current draw framebuffer
	  // (as will be the case for a windowless OpenGL context before you bind an FBO, for example).


	 }
	 
	  catch ( const runtime_error & error ) {
		 //std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	 }
	 
} // setupFBO

void ofApp::renderToFBO() {
	setupFBO();
	drawBackgroundStuff();
}

void ofApp::renderToSysBuffer() {

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	
	//glClearDepth(0.);

	//glClearDepth(1.0); // the default value

	//glDepthFunc(GL_GREATER);
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	float * bgPtr = ofBgColorPtr();

	// The following causes a blank screen, omitted by Moon Jung, 2014/8/4
	// it calls 	renderer->clear(r,g,b,a); which calls glClearcolor
	
	ofClear(bgPtr[0]*255,bgPtr[1]*255,bgPtr[2]*255, bgPtr[3]*255);
	
	//glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	drawBackgroundStuff();

}


void ofApp::drawBackgroundStuff()  {

    glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
  //setupFBO();

// Declare an empty uniforms

  Uniforms extraUniforms;
  // build & send proj. matrix to vshader
  
  const Matrix4 projMatrix = makeProjectionMatrix();

  extraUniforms.put("uProjectionMatrix", projMatrix);
   
  g_invEyeRbt = inv(g_eyeRbt);

  cout << "current eye frame:" << "\n" << g_eyeRbt << endl;


  const Cvec3 eyeLight1 = Cvec3(g_invEyeRbt * g_light1Pos ); // g_light1 position in eye coordinates
  const Cvec3 eyeLight2 = Cvec3(g_invEyeRbt * g_light2Pos ); // g_light2 position in eye coordinates
  
  // send the eye space coordinates of lights to uniforms
  extraUniforms.put("uLight1Pos", eyeLight1);
  extraUniforms.put("uLight2Pos", eyeLight2);

  extraUniforms.put("uLight1Color", g_light1Color );
  extraUniforms.put("uLight2Color", g_light2Color );

 
  // upload the sun direction relative to the eye space

  extraUniforms.put("uSunRayDir", Cvec3( g_invEyeRbt * Cvec4(g_sunRayDir,0) ) );
 // extraUniforms.put("uSunRayDir", g_sunRayDir );
  extraUniforms.put("uLightRGBColor", g_lightRGBColor);

  extraUniforms.put("uWindowWidth", g_windowWidth );
  extraUniforms.put("uWindowHeight",g_windowHeight);
	

  	  
	// (1) Draw the background objects to  FBO =  a container for textures and an optional depth buffer
	// FBO: there are no visible color buffer bitplanes, only a single "off-screen" color image attachment,
	//	so there is no sense of front and back buffers or SWAPPING.   
 
	// draw all the objects except for the last, which is a rainbow
    
   for ( int i  = 0; i < g_objectPtrList.size() - 1; i++ ) { // (5) This loop goes over the list of objects g_objectPtrList and draw each object in the list.
	                                                //     Where and how is this object list created?

     Matrix4 MVM = g_invEyeRbt * (  *( g_objectPtrList[i] -> objectRbt ) ) ; // g_currentPickedObject->objectRbt may have been
	                                                              // changed by picking
  
	 extraUniforms.put("uModelViewMatrix", MVM).put("uNormalMatrix", normalMatrix(MVM) );
   
	 // draw
	 try {
      g_objectPtrList[i]->draw( extraUniforms ); 
	  // This draw leads to  BufferedObjectGeometry::draw() which binds vertex buffers 
	  // and set vertex     attribute pointers, and then calls drawElements() */
	 
	  // unbind the current vao for the current mesh

	  //glBindVertexArray(0);

	  
	 }                                   
	 catch ( const runtime_error & error ) {
		 //std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	 }

	 
    } // for each object except for rainbow


} // drawBackgroundStuff() for ordinary rendering


void ofApp::renderRainbowUsingFBOTexture() {

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

   Uniforms extraUniforms;

      
    glBindFramebuffer(GL_FRAMEBUFFER, 0); // 0 = the system-provided buffer => unbind the current FBO
	//AFAIK, you can have one draw framebuffer bound at a time.
	// And that's the one that glClear should operate on. But this can be done only once?. 
	// No glClear should be done for every frame, if you do not intend to draw on only part of the screen

	glDrawBuffer(GL_BACK); // GL_BACK is the draw buffer for the zero framebuffer; for double-buffering context

	//glClearColor( g_clearColor[0], g_clearColor[1],g_clearColor[2], g_clearColor[3] );  
	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glViewport(0,0,g_windowWidth, g_windowHeight); 
		

    int rainbowIndex = g_objectPtrList.size() - 1;

	//extraUniforms.put("uColorTex", shared_ptr<ShaderImageTexture2D_RGB_RGB>(new ShaderImageTexture2D_RGB_RGB("sourceimages/sourceimages/Fieldstone.ppm", true)));

	// bind the "colorTex" and "depthTex" which are attached to FBO to the rainbow material shader
	extraUniforms.put("uRadius", g_radius );
    extraUniforms.put("uDropDensity", g_dropDensity );   
   
	// bind the FBO textures

	extraUniforms.put("uColorTex", 
		  shared_ptr<Texture_FROM_FBO>(new Texture_FROM_FBO( g_colorTex) ));
	extraUniforms.put("uDepthTex", 
		   shared_ptr<Texture_FROM_FBO>(new Texture_FROM_FBO(g_depthTex)));



    extraUniforms.put("uEyeMatrix", g_eyeRbt); // This will be used in shader to convert the
    extraUniforms.put("uInvEyeMatrix", g_invEyeRbt); // This will be used in shader to convert the
   	   
	
	Matrix4 uInvModelMatrixAABB = inv( *g_AABBRbt );

	cout << "current AABB frame:" << "\n" << *g_AABBRbt << endl;

	extraUniforms.put("uInvModelMatrixAABB", uInvModelMatrixAABB );

	Cvec3 uEyeOriginInBoxFrame = Cvec3( uInvModelMatrixAABB * 
		        Cvec4(g_eyeRbt[3], g_eyeRbt[7], g_eyeRbt[11], g_eyeRbt[15]) );
	  
	extraUniforms.put("uEyeOriginInBoxFrame", uEyeOriginInBoxFrame);

	g_AABBmin = Cvec3(  -g_AABBsize[0]/2.0, -g_AABBsize[1] / 2.0, -g_AABBsize[2]/2.0 );
	g_AABBmax = Cvec3(  g_AABBsize[0]/2.0, g_AABBsize[1] / 2.0, g_AABBsize[2]/2.0 );
	   

	extraUniforms.put("uAABBmin", g_AABBmin  );   
	extraUniforms.put("uAABBmax", g_AABBmax );   
	

	// draw rainbow
	 try {

      Matrix4 MVM = g_invEyeRbt * (  *( g_objectPtrList[ rainbowIndex] -> objectRbt ) ) ; // g_currentPickedObject->objectRbt may have been
	                                                              // changed by picking
  
	  extraUniforms.put("uModelViewMatrix", MVM).put("uNormalMatrix", normalMatrix(MVM) );
   
      g_objectPtrList[rainbowIndex]->draw( extraUniforms ); 

	  // This draw leads to  BufferedObjectGeometry::draw() which binds vertex buffers 
	  // and set vertex     attribute pointers, and then calls drawElements() 
	 
      checkGlErrors();
	  
	 }                                   


	 catch ( const runtime_error & error ) {
		 //std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	 }

} //  renderRainbowUsingFBOTexture() 



void ofApp::drawFBOTexture() {
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	// read the content of FBO to verify FBO is written well
    Uniforms extraUniforms;
      
    glBindFramebuffer(GL_FRAMEBUFFER, 0); // 0 = the system-provided buffer => unbind the current FBO
	//AFAIK, you can have one draw framebuffer bound at a time.
	// And that's the one that glClear should operate on. But this can be done only once?. 
	// No glClear should be done for every frame, if you do not intend to draw on only part of the screen

	glDrawBuffer(GL_BACK); // GL_BACK is the draw buffer for the zero framebuffer; for double-buffering context

	//glClearColor( g_clearColor[0], g_clearColor[1],g_clearColor[2], g_clearColor[3] );  
	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glViewport(0,0,g_windowWidth, g_windowHeight); 
		

	extraUniforms.put("uColorTex", 
		  shared_ptr<Texture_FROM_FBO>(new Texture_FROM_FBO( g_colorTex) ));
	extraUniforms.put("uDepthTex", 
		   shared_ptr<Texture_FROM_FBO>(new Texture_FROM_FBO(g_depthTex)));

	try {
      g_drawTextureObj->draw( extraUniforms ); 

	  // This draw leads to  BufferedObjectGeometry::draw() which binds vertex buffers 
	  // and set vertex     attribute pointers, and then calls drawElements() 
	 
      checkGlErrors();
	  
	 }                                   


	 catch ( const runtime_error & error ) {
		 //std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	 }

} //  drawFBOTexture() 


void ofApp::drawTextureFromFile( GLuint g_fbo) {

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	Uniforms  extraUniforms;

	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	extraUniforms.put("uTexColor", shared_ptr<ShaderImageTexture2D_RGB_RGB>(new ShaderImageTexture2D_RGB_RGB("colorFBO.ppm", false)));

    g_drawTextureObj->draw( extraUniforms );


}
void ofApp::drawTextureFromFile() {

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	Uniforms  extraUniforms;

	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	extraUniforms.put("uTexColor", shared_ptr<ShaderImageTexture2D_RGB_RGB>(new ShaderImageTexture2D_RGB_RGB("systembuffer.ppm", true)));
    g_drawTextureObj->draw( extraUniforms );


}



void ofApp::readPixels() {

	//unsigned char pixel[4] = {0};  // If you are using character types as numbers, use unsigned cha or signed char:
	// or unsigned int: you need to use it because you are dealing with bits in memory directly; or 
	// when doing manipulations such as bit masking or shifting on data,
	// for example when writing low-level code to read binary file formats such as audio files; 

 unsigned char *pixel = new unsigned char [ g_windowWidth * g_windowHeight * 4 ];

  // read the pixel at (x,y), get the color of the pixel
  // get the object id corresponding to the color
  // return the Rbt node corressponding to the ide

   //vector<float> rcolors(4);
  // glReadPixels(x, y, 1, 1, GL_RGBA, GL_FLOAT, &rcolors[0]);

 // for debugging
  glBindFramebuffer(GL_FRAMEBUFFER, g_fbo);  

  glReadBuffer(GL_COLOR_ATTACHMENT0); 
  // glReadBuffer(GL_BACK); /// default: BACK is where we draw in double buffering/
 // glReadBuffer(GL_FRONT);

// write the read pixel data to phaseFile
  /*
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
	
	*/ 
  glReadPixels(0,0, g_windowWidth, g_windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, pixel ); // GLvoid *pixels

  for (int j = 0; j < g_windowHeight; j++ )
	  for ( int i = 0; i < g_windowWidth ; i++) 
	 {
		  
		   // ostream& operator<< (unsigned int val);

		   cout << "pixel: " << "(" <<  j   << "," <<  i  << ") = "  
			                            << (int) pixel[j + i* g_windowWidth + 0] << ", " <<  (int)  pixel[j + i* g_windowWidth + 1] << ", " << 
										(int) pixel[j + i* g_windowWidth + 2]  << endl; 
		 

	}
	

    //textureFile.close();

	
} // readPixels()


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
void ofApp::initObjects() {


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





// // set a texture reference
//	void setUniformTexture(const string & name, ofBaseHasTexture& img, int textureLocation);
//	void setUniformTexture(const string & name, ofTexture& img, int textureLocation);
//	void setUniformTexture(const string & name, int textureTarget, GLint textureID, int textureLocation);
	
void ofApp::initMaterials() {


	// All the file names relative to the working directory where the project file is located.
	// get active uniforms and attributes and bind indices to them.

	MaterialShader diffuseMat( "basic+diffuse", "shaders/basic-gl3.vshader", "shaders/diffuse-gl3.fshader"); 
	
	// copy diffuse prototype and set red color
	g_redDiffuseMat.reset(new MaterialShader( diffuseMat ));  // MaterialShader(diffuse): uses the default copy constructor

	g_redDiffuseMat->getUniforms().put("uMaterialColor", Cvec4f(1, 0, 0,1.0));
	//g_redDiffuseMat->getUniforms().put("uTexUnit0", shared_ptr<ImageTexture>( new ImageTexture("reachup.ppm", true) ) );

	// copy diffuse prototype and set blue color
	g_blueDiffuseMat.reset(new MaterialShader(diffuseMat));
	g_blueDiffuseMat->getUniforms().put("uMaterialColor", Cvec4f(0, 1, 0, 1));

	
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
 
	//g_arcballMat->getUniforms().put("uMaterialColor", Cvec4f(0.27f, 0.82f, 0.35f, 1.0));
	g_arcballMat->getUniforms().put("uMaterialColor", Cvec4f(0.27f, 0.27f, 0.35f, 0.1));
	g_arcballMat->getRenderStates().polygonMode(GL_FRONT_AND_BACK, GL_LINE);
	
	//g_arcballMat->getRenderStates().polygonMode(GL_FRONT_AND_BACK, GL_POINT);
	// copy solid prototype, and set to color white
	g_lightMat.reset(new MaterialShader(solidMat));
	g_lightMat->getUniforms().put("uMaterialColor", Cvec4f(1, 1, 1,1));
	// pick shader
	g_pickingMat.reset(new MaterialShader("basic+pick", "shaders/basic-gl3.vshader", "shaders/pick-gl3.fshader") );

};


/*
void ofApp::initPBRTObjects() {


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
			if (!ParseFile(filenames[i])) { 
				cout << "Couldn't open scene file:" << filenames[i].c_str() << endl;
				Error("Couldn't open scene file \"%s\"", filenames[i].c_str());
			}
	}
	pbrtCleanup();




}  // initPBRTObjects

*/

void ofApp::initObjects() {
  int ibLen, vbLen;


// 1: create a ground
	//static const float groundY = -2.0;      // y coordinate of the ground
    static const float groundSize = 100.0;   // half the ground length

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
		  
    VertexPN( -groundSize, groundY, groundSize, 0, 1, 0), where groundY =0
    VertexPN( groundSize, groundY,  groundSize, 0, 1, 0),
    VertexPN(  groundSize, groundY, -groundSize, 0, 1, 0),
    VertexPN(  -groundSize, groundY, -groundSize, 0, 1, 0),
  };
  
  unsigned short idxGround[] = {0, 1, 2, 0, 2, 3};
  */

  

  //shared_ptr<Geometry> groundGeometry  ( new Geometry( &vtxGround[0], &idxGround[0], 4, 6,  &sqTex[0], numOfTexCoords  )  );
  shared_ptr<Geometry> groundGeometry  ( new SimpleIndexedGeometryPNTBX("ground", &vtxGround[0], &idxGround[0], vbLen, ibLen, GL_TRIANGLES  )  );

 
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

  vector<VertexPNTBX> vtxCube1(vbLen); // vtxCube:  a vector of VertexPNTBX type with length vbLen

  vector<unsigned short> idxCube1(ibLen); // idxCube: 

  // makecube<VertexPNTBX, unsigned short>( vtxCube.begin(), idxCube.begin() );

  // vtxCube and idxCube are vertices whose stroages are already
  // allocated, so that you can access each element of a vertex by means
  // of iterator.

  makeCube(30,6.4, 25, vtxCube1.begin(), idxCube1.begin() ); // fill vtxCube and idxCube, starting
                                                       // from their starting pointers => a list of vertex coords,
                                                       // tex coords, normal coords, tangent coords, binormal coords

  // create vbo and ibo for vtxCube and idxCube by using SimpleIndexedGeometryPNTBX, where
  // SimpleIndexedGeometry<VertexPNTBX, unsigned short> SimpleIndexedGeometryPNTBX;

 shared_ptr<Geometry> cubeGeometry1 (  new SimpleIndexedGeometryPNTBX("cube1",  &vtxCube1[0], &idxCube1[0], vbLen, ibLen, GL_TRIANGLES ) ); // reset(T * p)

  //shared_ptr<Geometry> cubeGeometry (  new Geometry( &vtxCube[0], &idxCube[0], vbLen, ibLen, &sqTex[0], numOfTexCoords  ) ); // reset(T * p)

 
 // 3: create another cube
 // int ibLen, vbLen;

  getCubeVbIbLen(vbLen, ibLen);

  // Temporary storage for cube geometry

  // T = tangent, X = coordinates

  vector<VertexPNTBX> vtxCube2(vbLen); // vtxCube:  a vector of VertexPNTBX type with length vbLen

  vector<unsigned short> idxCube2(ibLen); // idxCube: 

  // makecube<VertexPNTBX, unsigned short>( vtxCube.begin(), idxCube.begin() );

  // vtxCube and idxCube are vertices whose stroages are already
  // allocated, so that you can access each element of a vertex by means
  // of iterator.

  makeCube(20, 12.7, 25, vtxCube2.begin(), idxCube2.begin() ); // fill vtxCube and idxCube, starting
                                                       // from their starting pointers => a list of vertex coords,
                                                       // tex coords, normal coords, tangent coords, binormal coords

  // create vbo and ibo for vtxCube and idxCube by using SimpleIndexedGeometryPNTBX, where
  // SimpleIndexedGeometry<VertexPNTBX, unsigned short> SimpleIndexedGeometryPNTBX;

 shared_ptr<Geometry> cubeGeometry2 (  new SimpleIndexedGeometryPNTBX("cube2",  &vtxCube2[0], &idxCube2[0], vbLen, ibLen, GL_TRIANGLES ) ); // reset(T * p)

  //shared_ptr<Geometry> cubeGeometry (  new Geometry( &vtxCube[0], &idxCube[0], vbLen, ibLen, &sqTex[0], numOfTexCoords  ) ); // reset(T * p)


  // 4: create a sphere
  
 /*
  int slices = 100;
  int stacks = 100;

  getSphereVbIbLen(slices, stacks, vbLen, ibLen);

  // Temporary storage for cube geometry

  vector<VertexPNTBX>  vtxSphere (vbLen);

  vector<unsigned short> idxSphere (ibLen);

  float radius = 3;

  makeSphere(radius, slices, stacks,  vtxSphere.begin(), idxSphere.begin() );

  //geometry.h: typedef SimpleIndexedGeometry<VertexPNTBX, unsigned short> SimpleIndexedGeometryPNTBX;

  shared_ptr<Geometry> sphereGeometry ( new SimpleIndexedGeometryPNTBX( &vtxSphere[0], &idxSphere[0], vbLen, ibLen, GL_TRIANGLES ) ); // reset(T * p)
 
  // or  shared_ptr<Geometry> sphereGeometry =  new Geometry( &vtx[0], &idx[0], vbLen, ibLen ); // reset(T * p)
  // NOTE: here vtx and idx arrays are copied to the vertex buffer object internally, so this data need not be global.

  */


  // create the object frames for the created geometry

  // The following does not work, because the  matrix created by SgRbtNode::makeTranslation() is
  //	destroyed outside of this function initObjects(). To avoid it, create the matrix by the "new" method.
 //	The objects created by new methods remain unless removed explicitly. This is handled by shared_ptr<>. 
																				   
//  shared_ptr<SgRbtNode> groundRbt ( &SgRbtNode::makeTranslation( Cvec3(0,0,0) ) ); 
 // shared_ptr<SgRbtNode>  cubeRbt ( &SgRbtNode::makeTranslation( Cvec3(5.0,0,0)  ) );
  // shared_ptr< SgRbtNode>  sphereRbt ( &SgRbtNode::makeTranslation( Cvec3(0,0,0) ) );
  
   
	
  shared_ptr<SgRbtNode> groundRbt ( new SgRbtNode( g_SceneRbt * SgRbtNode::makeTranslation( Cvec3(0,0,0) ) ) ); // this uses the non-default copy constructor ?
  shared_ptr<SgRbtNode> cubeRbt1 ( new SgRbtNode( g_SceneRbt * SgRbtNode::makeTranslation( Cvec3(0, 3.2, -12.5)  ) )  );
 // shared_ptr< SgRbtNode> sphereRbt ( new SgRbtNode ( SgRbtNode::makeTranslation( Cvec3(-0.5,0,0) ) ) );
   shared_ptr<SgRbtNode> cubeRbt2 ( new SgRbtNode( g_SceneRbt * SgRbtNode::makeTranslation( Cvec3(25, 6.35, -12.5)  ) )  );
  

  
	   
  shared_ptr<Object> ground ( new Object( "ground",  groundRbt, groundGeometry, g_bumpFloorMat) );
  shared_ptr<Object> cube1 ( new Object( "cube1", cubeRbt1, cubeGeometry1, g_blueDiffuseMat) );
  shared_ptr<Object> cube2 ( new Object( "cube2", cubeRbt2, cubeGeometry2, g_redDiffuseMat) );
  //shared_ptr<Object> sphere ( new Object( "sphere", sphereRbt, sphereGeometry, g_blueDiffuseMat) );

  // 1: create the object id color for picking purpose 

  Cvec4 pickColor;
  unsigned int id; // 32 bit unsigned


  
  id = g_objId ++ ; 

  id = id << 8;  // RGBA: fill A field with 0's, shifting the value of id left 8 bits
  id = id | 255 ; // id for  ground; make A field to be all 1's



  pickColor = g_picker.idToColor(id); // id => linearTrans(idColor)

  ground->pickColor = pickColor; // the   pickColor = linearTrans(idColor) will be sent to the shader Uniform variable when drawing 
                                 // That is, by fragColor = linearTrans(id). WHen written to the framebuffer, sRGBTrans(linearTrans(idColor)) =
                                 //  idColor will written. When read by glReadPixel, this value is 
                                 // read. 
  
  cout << "ground object id=" << id <<"," << "object pseudo color=" << pickColor << endl;


  // store the objects in the object list

  g_objectPtrList.push_back( ground );

  // add the object and the id to the list for picking
  g_picker.addToMap( id, ground );


  
  // 2: the same for the cube
 
  id= g_objId ++;
  id = id << 8;  
  id = id | 255; // id for  cube 
  pickColor  = g_picker.idToColor(id );

  cube1->pickColor = pickColor;
 
  cout << "cube1 object id=" << id <<"," << "object pseudo color=" << pickColor << endl;

  g_objectPtrList.push_back( cube1 );


  g_picker.addToMap( id, cube1 );

   
  
  
  // 4: the same for sphere
  id = g_objId ++ ;
  id = id << 8;
  id = id | 255; // id for  sphere 

  pickColor = g_picker.idToColor(id);

  //sphere->pickColor  = pickColor;
  cube2->pickColor  = pickColor;

  cout << "cube2 object id=" << id <<"," << "object pseudo color=" << pickColor << endl;

  g_objectPtrList.push_back( cube2 );

  g_picker.addToMap( id, cube2 );

  
 

};

 void ofApp::initAABB() {
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


 // vtxAABB will have vertex coordinates relative the center of the box

// makeAABB( g_AABBsize[0], g_AABBsize[1], g_AABBsize[2], vtxAABB.begin(), idxAABB.begin() ); // fill vtxAABB  starting
                                                       // from the  starting pointer with a GenericVertex which consists of  vertex coords,
                                                       // tex coords, normal coords, tangent coords, binormal coords

  makeCube( g_AABBsize[0], g_AABBsize[1], g_AABBsize[2], vtxAABB.begin(), idxAABB.begin() ); 

  // create vbo and ibo for vtxCube and idxCube by using SimpleIndexedGeometryPNTBX, where
  // SimpleIndexedGeometry<VertexPNTBX, unsigned short> SimpleIndexedGeometryPNTBX;

  // Upload the geometry data and bind it to the Vbo 
 //shared_ptr<Geometry> AABBGeometry (  new SimpleIndexedGeometryPNTBX("AABB", &vtxAABB[0], &idxAABB[0], vbLen, ibLen, GL_QUADS ) ); 
 // GL_QAUDS has been removed from the core profile => ENUM error is caused
  shared_ptr<Geometry> AABBGeometry (  new SimpleIndexedGeometryPNTBX("AABB", &vtxAABB[0], &idxAABB[0], vbLen, ibLen, GL_TRIANGLES ) );

 
// for National Modern Musium
// 	g_AABBRbt.reset(  new SgRbtNode (  g_SceneRbt * 
//		                       Matrix4::makeTranslation( 	Cvec3( 0, g_AABBsize[1]/2 + 6.4, -g_AABBsize[2]/2.0 ) ) ) );
	 
	// for Asian Game [the world frame is assumed right at the front of the AABB box

 	g_AABBRbt.reset(  new SgRbtNode (  g_SceneRbt * 
		                       Matrix4::makeTranslation( 	Cvec3( 0, g_AABBsize[1]/2, -g_AABBsize[2]/2.0) ) ) );
	  
  shared_ptr<Object> AABB ( new Object( "AABB",  g_AABBRbt, AABBGeometry, g_arcballMat) );

   
  Cvec4 pickColor;
  unsigned int id; // 32 bit unsigned

  
  id= 0 ; // set the  AABB box id to zero. // AABB box can be picked and moved.
  id = id << 8;  
  id = id | 255; // 
  pickColor  = g_picker.idToColor(id ); // pickColor is an sRGB non-linear color

  cout << "AABB object id=" << id <<"," << "object pseudo color=" << pickColor << endl;

  AABB->pickColor = pickColor;


  g_objectPtrList.push_back( AABB );


  g_picker.addToMap( id, AABB );

 
 }
 
  void ofApp::initRainbow() {
  int ibLen, vbLen;

  
  static const float quadSize = 0.5;   // half the quad length

  // temporary storage for a plane geometry
  getPlaneVbIbLen( vbLen, ibLen); // get the sizes of the vertex buffer and the 
                                 // index buffer for a plane. The vertex buffer
                                 // size is the number of vertices.
  vector<VertexPNTBX> vtxQuad (vbLen);
  vector<unsigned short> idxQuad (ibLen);

  makeQuad( quadSize * 2, vtxQuad.begin(), idxQuad.begin() );

  
  //shared_ptr<Geometry> groundGeometry  ( new Geometry( &vtxGround[0], &idxGround[0], 4, 6,  &sqTex[0], numOfTexCoords  )  );
  shared_ptr<Geometry> quadGeometry  ( new SimpleIndexedGeometryPNTBX("rainbow", &vtxQuad[0], &idxQuad[0], vbLen, ibLen, GL_TRIANGLES  )  );

    // rainbow shader

  MaterialShader rainbowMat ("basic+rainbow", "shaders/FullScreenQuad-gl3.vshader", "shaders/rainbow-gl3.fshader");
  
 // MaterialShader initializes  its member programDesc_ which has a member program, 
  // by programDesc_ ( GlProgramLibrary::getSingleton().getProgramDesc(vsFilename, fsFilename) ), which in turn
  // sets the output variable for shaders, e.g.:
  //if (!g_Gl2Compatible) {
		// index 0, 1, 2 refers to the index of the array attachments in:
		// GLuint attachments[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1,GL_COLOR_ATTACHMENT2 }; 
		// glDrawBuffers(3,  attachments);

   //   glBindFragDataLocation(program, 0, "fragColor");
  //	  glBindFragDataLocation(program, 1, "phaseSpectrum");
  //	  glBindFragDataLocation(program, 2, "rainbowColorSpectrum");
  //	}



  g_rainbowMat.reset(new MaterialShader( rainbowMat) );

  g_rainbowMat->getUniforms().put("uScatTex", shared_ptr<ShaderImageTexture2D_RF_RF>(new ShaderImageTexture2D_RF_RF("rainbow/scattering2_120_fort.txt")));


  g_rainbowMat->getUniforms().put("uPhaseTex", shared_ptr<ShaderImageTexture3D_RF_RF>(new ShaderImageTexture3D_RF_RF("rainbow/phaseFunction2_120_100_fort.txt")));

 
  g_quadRbt.reset(  new SgRbtNode (  Matrix4::Matrix4()  ) );

  
//  g_rainbowObj.reset(   new Object( "rainbow",  g_quadRbt, quadGeometry, g_rainbowMat) );
//   g_rainbowObj.reset(   new Object( "rainbow",  NULL, NULL, NULL) );
   shared_ptr< Object > rainbowObj;
   rainbowObj.reset(   new Object( "rainbow",  g_quadRbt, quadGeometry, g_rainbowMat) );

 //  shared_ptr< Object > rainbow (  new Object( "rainbow",  g_quadRbt, quadGeometry, g_rainbowMat) );
  //  Construct an instance of shared_ptr which takes ownership of
  // pointer  new Object( "rainbow",  g_quadRbt, quadGeometry, g_rainbowMat) 
    
  g_objectPtrList.push_back( rainbowObj );
 //g_objectPtrList.push_back( NULL );



 } // initRainbow()
 
  void ofApp::initDrawTextureMat() {
  int ibLen, vbLen;

  getPlaneVbIbLen(vbLen, ibLen);
  
  static const float quadSize = 0.5;   // half the quad length

  // temporary storage for a plane geometry
  getPlaneVbIbLen( vbLen, ibLen); // get the sizes of the vertex buffer and the 
                                 // index buffer for a plane. The vertex buffer
                                 // size is the number of vertices.
  vector<VertexPNTBX> vtxQuad (vbLen);
  vector<unsigned short> idxQuad (ibLen);

  makeQuad( quadSize * 2, vtxQuad.begin(), idxQuad.begin() );

  
  //shared_ptr<Geometry> groundGeometry  ( new Geometry( &vtxGround[0], &idxGround[0], 4, 6,  &sqTex[0], numOfTexCoords  )  );
  shared_ptr<Geometry> quadGeometry  ( new SimpleIndexedGeometryPNTBX("drawTexture", &vtxQuad[0], &idxQuad[0], vbLen, ibLen, GL_TRIANGLES  )  );

    // rainbow shader

  MaterialShader drawTextureMat ("FBOToSysBuffer", "./shaders/FullScreenQuad-gl3.vshader", "./shaders/FBOtexture-gl3.fshader");
 

  g_drawTextureMat.reset(new MaterialShader( drawTextureMat) );

  
  g_quadRbt.reset(  new SgRbtNode (  Matrix4::Matrix4()  ) );

  g_drawTextureObj.reset ( new Object( "drawTexture",  g_quadRbt, quadGeometry, g_drawTextureMat) );

 } // initDrawTextureMat()
 

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){  // call back function for event processing
 
	// instance->windowW = w; => the new window size is already set before calling windowResized()

	// instance->windowH = h;
  g_windowHeight = h;
  g_windowWidth = w;
   
  glViewport(0, 0, g_windowWidth, g_windowHeight);
  
  cout   << "WINDOW RESIZED: Size of window is now " << g_windowWidth << "x" << g_windowHeight << endl;
  cout << "results of ofGetHeight and ofGetWidth= " << ofGetWidth() << "x" << ofGetHeight() << endl;

  //cerr  << "Size of window is now " << g_windowWidth << "x" << g_windowHeight << endl;
  //system ("PAUSE");

  updateFrustFovY();
 
  // At every frame, display() callback is called to redraw the
	// scene. So, windowResized event handler needs not do anything to redraw the scene.

   g_redrawWindowEvent = true;
	
}


//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){ // call back function for event processing

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

	   g_redrawWindowEvent = true;

  } // if (! g_rotation)

  else if ( ofGetKeyPressed( OF_KEY_CONTROL) && ofGetKeyPressed( OF_KEY_ALT) ) { //  rotate the camera

	  switch ( button ) {
		  case OF_MOUSE_BUTTON_LEFT:

			  //inline static ofMatrix4x4 newRotationMatrix( float angle, const ofVec3f& axis);
			  if ( abs(dy) > abs(dx)  ) {
				  cout << " rotate about the x axis" << endl;

			      m = Matrix4::makeXRotation(dy * 0.01); // unit= degree; pitch motion: up and down
			  }
			  else {
				  cout << " rotate about the y axis" << endl;
				  m = Matrix4::makeYRotation( dx * 0.01 ); // yaw motion: left and right
			  }
			   
			  break;
		  case OF_MOUSE_BUTTON_RIGHT:
			  cout << " rotate about the z axis" << endl;
			  m = Matrix4::makeZRotation( -dx * 0.01 ); //  rolling about the camera view direction
			  break;

		  default: 
			  m = Matrix4::makeTranslation( Cvec3( 0.0, 0.0, 0.0)  ); 
	  } // switch

	  g_eyeRbt *= m;

	  //ofPtr<ofAppBaseWindow> 		window;  // changed by Moon Jung, 2014/4/23
	  
	   // window = ofPtr<ofAppBaseWindow>(new ofAppGLFWWindow());
	   //ofPtr<ofAppGLFWWindow>	appwindow =  window;
	   // 	  appwindow  ->display();

	  g_redrawWindowEvent = true;
  } // else if
  
  

  else if (   ofGetKeyPressed( OF_KEY_SHIFT ) && !ofGetKeyPressed( OF_KEY_ALT) )   { // translate the selected object


     if ( g_currentPickedObject == nullptr ) { // no object is pickefd up, so no need to move
		 return;
	 }
	  

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

	  *(g_currentPickedObject->objectRbt) = *(g_currentPickedObject-> objectRbt) * m;

	  if ( g_currentPickedObject->objectName =="AABB" ) {
		  
		  g_AABBRbt = (g_currentPickedObject->objectRbt);
		  	 
	  }

	   g_redrawWindowEvent = true;

    } // if ( translate )

    else if ( ofGetKeyPressed( OF_KEY_SHIFT) && ofGetKeyPressed( OF_KEY_ALT) ) { //  rotate the selected object
	   
          if ( g_currentPickedObject == nullptr ) { // no object is picked up, so no need to move
		     return;
	      }

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

		  *(g_currentPickedObject->objectRbt) =  *(g_currentPickedObject-> objectRbt) * m;

		  if ( g_currentPickedObject->objectName =="AABB" ) {
		  
		     g_AABBRbt = (g_currentPickedObject->objectRbt);
		  	 
	      }
		  
		  g_redrawWindowEvent = true;
       } // else if

   
} // mouseDragged()





//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

	// no object is selected now, once the mouse is released.

	g_currentPickedObject = nullptr;

} // mouseReleased()




//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

	// This callback is called every time the pressed mouse is moved a little.
	// x, y: the current mouse position after this little  motion is finished.


	
}


//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){
	
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
	
	drawForPicking();    
	//    In the picking mode, you render the simplified version of the scene to the back buffer,
   //     while maintaining the scene of the screen intact, because it is in the front buffer.
		                      
   //  Read pixel operation needed for picking is performed with respect to the back buffer,
	//  not to the front buffer
		
	g_currentPickedObject = g_picker.getRbtNodeAtXY( g_prev_mouseClickX, g_prev_mouseClickY );

	 // g_currentPickedRbtNode points to one of g_objectRbt[], which has been picked.

     // g_currentPickedRbtNode will be NULL, if no object has been actually picked.

   if ( g_currentPickedObject == nullptr ) {

		  cout << "no object has been picked" << endl;
		  //g_picked_mode = false;

	 }
   else {
          cout << "an object has been picked" << endl;

		 // g_picked_mode = true;
     }

  } // g_picking_mode

  // after drawForPicking(), normal draw is called, which is in fact called every frame.

} // mousePressEvent()




//--------------------------------------------------------------
void ofApp::keyPressed  (int key){ 
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

  case 'f':  // debug mode, render to texture and print the content
	  renderToFBO();
	  break;

  case 'r':  // debug mode, read pixels
	  readPixels();
	  break;

  case OF_KEY_ESC: 
    exit();                                  // ESC 
  //case OF_KEY_CONTROL: 
  //  g_camera_mode = true;
  //	break;

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
void ofApp::keyReleased(int key){ 


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
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}


void ofApp::guiEvent(ofxUIEventArgs &e)
{
    string name = e.widget->getName();
    if(name == "SHOW ACTIVE")
    {
        ofxUIToggle *toggle = (ofxUIToggle *) e.widget;
        ddl->setShowCurrentSelected(toggle->getValue()); 
    }
    else if(name == "DROPDOWN")
    {
        ofxUIDropDownList *ddlist = (ofxUIDropDownList *) e.widget;
        vector<ofxUIWidget * > &selected = ddlist->getSelected(); 

		for(int i = 0; i < selected.size(); i++) {  
          string action_name = selected[i]->getName(); 
		  //cout << action_name << endl;
		  //break;
	
	     if ( action_name == "render background to system framebuffer") {
		    renderToSysBuffer(); // to system framebuffer
		    
	     }
	     else if ( action_name == "dump the system framebuffer to file" ) { // system buffer dump
		     glFlush();
             writePpmScreenshot(g_windowWidth, g_windowHeight, "systembuffer.ppm");
	     }
	     else if ( action_name == "render background to FBO") { // render to FBO

		     renderToFBO();
	     }
	     else if ( action_name == "dump FBO to file" ) { // FBO dump
             printFBO(g_fbo, g_windowWidth, g_windowHeight, "colorFBO.ppm" );
	
         }
	     else if (action_name == "render background to system framebuffer using FBO texture" ) { // render to system buffer using FBO texture
		     drawFBOTexture();
		 }
	     else if (action_name == "render background from FBO file") { // render to system buffer using dumped FBO texture
		     drawTextureFromFile(g_fbo); // this will give a non gamma corrected image
		                          // because the FBO image is not in a sRGB format.

		 }
	     else if (action_name  == "render background from system framebuffer file") { // render to system buffer using dumped FBO texture
		    drawTextureFromFile();
		 } 
	     else if (action_name  == "render rainbow using FBO background") { // render to system buffer using dumped FBO texture
	      	renderRainbowUsingFBOTexture();
		 } 
    
      } // for


	} // else if ("DROPDOWN")

	
} // ofApp::guiEvent

void ofApp::initGLState() {



	glBindFramebuffer(GL_FRAMEBUFFER, 0); 
	// Once a FBO is bound, all OpenGL operations affect onto the current bound framebuffer object.

	/* A Framebuffer is a collection of buffers that can be used as the destination for rendering.
		OpenGL has two kinds of framebuffers: the Default Framebuffer, which is provided by 
		the OpenGL Context; and user-created framebuffers called Framebuffer Objects (FBOs). 
		The buffers for default framebuffers are part of the context and usually represent a 
		window or display device. The buffers for FBOs reference images from either Textures 
		or Renderbuffers; they are never directly visible. */

	messageFile  << "OpenGL Version : " << glGetString(GL_VERSION) << endl; 

	//glClearColor(128./255,200./255,1,0); // sky color
	
	//glClearColor(0./255., 0./255., 0./255., 1.0); // black  color
	
	// backup the original background color, which is used for the background color after picking
	// which uses another background color

	glGetDoublev(GL_COLOR_CLEAR_VALUE, g_clearColor); // glGetIntegerv() // glGetDoublev() integer or double values

	checkGlErrors();
	
	//glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // how pixel data is written to or read from

	//glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glCullFace(GL_BACK);
	//glEnable(GL_CULL_FACE);
	glDisable(GL_CULL_FACE);

	glEnable(GL_DEPTH_TEST);
	//glDisable(GL_DEPTH_TEST);
	glClearDepth(0.);

	//glClearDepth(1.0); // the default value

	glDepthFunc(GL_GREATER);
	//glDepthFunc(GL_LEQUAL);


	glColorMask( GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

	glDisable(GL_BLEND);
	checkGlErrors();
	
	//Also note that using the constants GL_ONE (source) and GL_ZERO (destination) 
	//gives the same results as when blending is disabled

	/*
	glDisable(GL_ALPHA_TEST); 
	// fixed-pipe line alpha test deprecated in opengl 3.2 onward,
	// which we use in this program. 

	try {
		checkGlErrors();
	}
	catch (const runtime_error & error ) {
		std::cout << error.what() << endl;
		messageFile << error.what() << endl;
		//throw;
	}
	*/

	// the alpha test to accept or reject a fragment based on its alpha value
	//If enabled, the test compares the incoming alpha value with a reference value. 
	// The fragment is accepted or rejected depending on the result of the comparison. Both the reference value and the comparison function are set with glAlphaFunc(). 
	//By default, the reference value is zero, the comparison function is GL_ALWAYS[ always accept the fragment ], and the alpha test is disabled

	// The default framebuffer has buffers like GL_FRONT​, GL_BACK​, GL_AUXi​, GL_ACCUM​, and so forth. 
	// FBOs do not have these. Instead, FBOs have a different set of images.
	// Each FBO image represents an attachment point, a location in the FBO where an image can be attached.
	//   FBOs have the following attachment points:

    //      GL_COLOR_ATTACHMENTi, GL_DEPTH_ATTACHMENT, GL_STENCIL_ATTACHMENT,GL_DEPTH_STENCIL_ATTACHMENT​:
		
	//Each framebuffer object has a set of attachment points that logical buffers can be attached to. 
	
	// There are three different ways to switch between framebuffers. 

	//The first one is to use several different FBOs, one for each combination of logical buffers 
	//	that you plan on using in you application. To change render targets you simply call glBindFramebuffer() 
	//	with the FBO containing the setup you wish to render to.
	// Another way is to use a single FBO and alter the attachments
	// The third way is to use a single FBO, but instead of altering attachments you call glDrawBuffer() or glDrawBuffers() 
	// to change which color attachment(s) the rendering goes to.

	// for debugging: glDrawBuffer(GL_BACK); // GL_BACK is the draw buffer for the zero framebuffer; for double-buffering context
		
	
	//if (!g_Gl2Compatible)  glEnable(GL_FRAMEBUFFER_SRGB); // For opengl3.x,  request an sRGB frame buffer using the call glEnable(GL FRAMEBUFFER SRGB).
                                       // Then we can pass linear [R, G, B] values out from the fragment shader. 
	                                   //  and they will be gamma corrected into the sRGB format before begin sent to the framebuffer.
		// Any writes to images that are not in the sRGB format should not be affected. 
		//	So if you're writing to a floating-point image, nothing should happen.
		//	Thus, you should be able to just turn it on and leave it that way; 
		// OpenGL will know when you're rendering to an sRGB framebuffer.							   
	

	try {
		checkGlErrors();
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
	==>  if a framebuffer object is bound that is not the default framebuffer, then only
	 GL_NONE, GL_FRONT_LEFT, GL_FRONT_RIGHT, GL_BACK_LEFT, and GL_BACK_RIGHT, GL_COLOR_ATTACHMENTi
	 are accepted. 
	 GL_NONE: you will render the depth map, but color map:
	 
	 The depth in a GL_DEPTH_COMPONENT texture is written automatically right?
	 Yep, assuming you have depth writes enabled, and the texture bound as the depth attachment of an FBO.


   ## specifying outout variables:
	//layout (location = 0) out vec4 normalOut; 
	// layout (location = 1) out vec4 texCoordOut; 

	// the location is related to the indexes used as an argument in glDrawBuffers.
	// For instance, considering the above example for glDrawBuffers, 
	// normalOut = X will write X to location 0, which is to color attachment 1, and
	// texCoordOut = X will write X to  location 1, which is  color attachment 2.

	OR: use the following instead of layout ...:
	       glBindFragDataLocation(program, 0, "normalOut");
  //	  glBindFragDataLocation(program, 1, "texCoordOut");
  //	More than one varying out variable is bound to the same color
  //    number => the program fails to link. If a shader statically
  // assigns a location to a varying out variable in the shader text, 
  // that location is used and any location assigned with glBindFragDataLocation 
  // will be ignored. 
	
*/

/* The X, Y and Z values in glVertex() are given in object coordinates. Object coordinates are different coordinate from what glClearDepth() and glDepthRange() use. glClearDepth() uses the range 0..1 for depth values, as used in depth buffer.•Object Coordinates are transformed by the ModelView matrix to produce Eye Coordinates
•Eye Coordinates are transformed by the Projection matrix to produce Clip Coordinates.
•Clip Coordinate X, Y, and Z are divided by Clip Coordinate W to produce Normalized Device Coordinates.
•Normalized Device Coordinates are scaled and translated by the viewport parameters to produce Window Coordinates.
Effectively your projection near is mapped to 0.0 and far is mapped to 1.0, if you use default depth range. Thus glClearDepth(0.0) means "clear depth buffer with near plane" and glClearDepth(1.0) means "clear depth buffer with far plane.

 Depth testing happens with the window Z coordinate, which is in 0 to 1 range. This is the value that gets written to the depth buffer if depth and other tests pass. Further away fragments have larger value, in usual setups.

 If you have specified non-default values to glDepthRange(), then depth values will be mapped to your choice of depth range, however this range still needs to be in 0..1. glDepthRange(1.0, 1.0) makes all rendered fragments to go to the far plane. Notice that near and far given to glDepthRange() are different from near and far used in projections.

 Typically you clear with "far plane", 1.0, and render with depth test less / lequal. If you want to render the away facing surfaces instead, you can clear with "near plane", 0.0 and render with depth test greater / gequal. You also need to pay attention to have your backface culling setup properly. 
 
 The above descriptiop holds when you use the [left-handed] normalized device coordinates, where 
 nearer means smaller z_n.
 In OUR implementation, z_e = f maps to z_n = -1, and z_e = n maps to z_n = 1, p. 106, Gortler's book
 That is, we need the typical right-handed coordinate system for NDC. This is mapped to window coordinates
 : z_n = -1 = far => z_w = 0; z_n = 1 = near => z_w = 1. So, we must explicitly 
 tells Opengl to use this system by saying glDepthFunc(GL_GREATER), and glCLearDepth (0.0), where
0.0  is  the Widows coordinate depth of the far plane. 
 */
 

} // initGLStates

/* 
With program scope (3), separately compiled modules require access to 
function names, structures, and objects, often placed in global namespace.
Global namespace includes all names from global declarations;
that is, declarations appearing outside of function definitions. 
What is the proper way to deal with global namespaces? This problem can be significant, 
especially with applications and libraries. Let's investigate how namespaces can help.

namespace Blue {      // original namespace definition
  int j;
  void print(int);
 }

 namespace Blue {      // namespace extension
  char ch;
  char buffer[20];
 }
 . . .

 namespace BigBlue {     // equivalent to the above
  int j;
  void print(int);
  char ch;
  char buffer[20];
 }
 
When you create a namespace, avoid placing include files inside the namespace definition.
 // myfile.C
 namespace math {
 #include "geometry.h"      // not a good idea
  Point origin = { 3, 5 };
  . . .
 }

This approach creates link errors with other modules that do not include geometry.h in the same namespace. Instead, define a namespace in the include file itself, and use namespace extensions for modules that need to add to it.
 // geometry.h
 namespace math {      // define namespace
  struct Point { double x, y; };
  double slope(Point, Point);
 }

 // myfile.C
 #include "geometry.h"
 namespace math {       // namespace extension
  Point origin = { 3, 5 };
  . . .
 }

  namespace Blue {          // define namespace Blue
  int j;              // j is a member of namespace Blue
  void print(int);         // print() is a member of namespace Blue
 }

 using namespace Blue;        // global using directive

 void sub1() {
  j = 0;              // legal - Blue::j
  print(j);            // legal - Blue::print()
  . . .
 }

 void sub2() {
  j = 1;             // legal - Blue::j
  print(j);            // legal - Blue::print()
  . . .
 }

The global using directive provides access to j and print()
from namespace Blue without a qualifier and without the scope operator.
Note that j in sub1() and sub2() are not local variables.


Here are several examples of using declarations.
 namespace Black {           // define namespace Black
  int j;
  void print(int);
  char ch;
 }

 namespace White {           // define namespace White
  int j;
  void print(int);
  double vision;
 }
 using White::vision;         // global using declaration

 void sub1() {
  using Black::print;        // local using declaration
  ch = 'a';             // illegal - ch not defined
  vision = 7.65;           // legal - White::vision
  print(5);             // legal - Black::print()
  . . .
 }

 void sub2() {
  print(5);             // illegal - print() not defined
  vision = 3.45;          // legal - White::vision
  . . .
 }

 namespace A_Very_Long_Library_Name {
  struct Point { double x, y; };
  Point origin = { 10, 10 };
 }

 namespace ALN = A_Very_Long_Library_Name;   // alias

 void sub() {
  using ALN::Point;              // using declaration
  Point first = ALN::origin;         // namespace member
  . . .
 }

 Listing 3.14 geometry.h definitions
 #ifndef GEOMETRYH
 #define GEOMETRYH
 // geometry.h - geometry definitions

 namespace Geometry {       // namespace definition
  struct Point {
    double x, y;
  };
  double slope(Point, Point);
 }
 #endif

Namespace Geometry includes a Point type and a slope() function. Additional definitions appear in geometry.C.

Listing 3.15 geometry.C implementations
 // geometry.C - geometry implementations
 #include "geometry.h"

 namespace Geometry {       // namespace extension
  Point origin = { 0, 0 };
 }

 double Geometry::slope(Point a1, Point a2) {
  double dy = a2.y - a1.y;
  double dx = a2.x - a1.x;
  if (dx == 0)
    throw "slope(): undefined slope";
  return dy / dx;
 }



*/

