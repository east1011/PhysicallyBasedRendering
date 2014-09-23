#include "ofApp.h"
#include "spectrum.h"


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



// for debugging: exclude pbrt for the time being, by Moon Jung, 2014/8/7
#include "stdafx.h"
#include "api.h"
#include "probes.h"
#include "parser.h"
#include "parallel.h"




bool ParseFile(const string &filename);
extern int g_argc;
extern char **g_argv;


// Just avoid any using namespace directives in the headers. 
// (And use with care in C++ files, if at all. Inside functions is preferred.)
// as long as they  * appear after all #includes.

using namespace std;      // for string, vector, iostream, shared_ptr and other standard C++ stuff
using namespace std::tr1; // for shared_ptr<T>



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
						   g_rainbowAndSceneMat,
						   g_rainbowOnlyMat,
						   g_AABBRainbowOnlyMat,
						   g_fullScreenMat,
						   g_NDCMat;

shared_ptr<MaterialShader> g_overridingMaterial;

shared_ptr<Object> g_NDCMatObj; // used simply to draw a texture
shared_ptr<Object> g_fullScreenMatObj, g_rainbowAndSceneMatObj, g_rainbowOnlyMatObj; // used simply to draw a texture
shared_ptr<Object> g_AABBRainbowOnlyMatObj, g_oldAABBRainbowOnlyMatObj;;

Matrix4 g_projectionMatrix; 

int  g_objId = 1; 


double g_frustMinFov, g_frustFovY;
double g_frustNear, g_frustFar;

int g_windowWidth; 
int g_windowHeight; 


bool g_mouseClickDown = false;    // is the mouse button pressed
bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;

int g_prev_mouseClickX, g_prev_mouseClickY; // coordinates for mouse click event

int g_pressed_button;

int g_activeShader = 0;

bool g_camera_mode = false;
bool g_backgroundBackedUp = false;

bool g_picking_mode = false; // mode during the process of picking
bool g_picked_mode = false;

bool g_rotation_mode = false;



GLdouble g_clearColor[4];

GLuint g_fboScene, g_colorTex, g_dataTex, g_depthTex; 
GLuint g_fboRainbow, g_colorTex1, g_dataTex1, g_depthTex1, g_dataTex2;  ;

bool g_debugMode = true;

// for rainbow drawing
// The coordinates of the AABB box relative to the  world coordinate system

//float g_phaseFunction[ nRadii * nSpectralSamples * nThetas ];

// This includes the hemisphere for sky
//static const Cvec3 g_initAABBmax(90.3325, 75.8859, 90.3325);
//static const Cvec3 g_initAABBmin(-90.3325, -14.4466, -90.3325);

// This includes the hemisphere for sky

// for National Modern Musium
//static const float rainbowVolumeHeight = 21.3; 
//static const float rainbowVolumeWidth = 30;
//static const float rainbowVolumeDepth = 5;

// for Asian Game 
//static const float rainbowVolumeHeight = 30; 
//static const float rainbowVolumeWidth = 40;
//static const float rainbowVolumeDepth = 10;

// for testing by friends 
static const float rainbowVolumeHeight = 50; 
//static const float rainbowVolumeWidth = 60;
static const float rainbowVolumeWidth = 120;
static const float rainbowVolumeDepth = 10;



// The coordinates of the AABB, which will be used to compute the sizes of the view volume

static Cvec3 g_AABBmax, g_AABBmin; 
static Cvec3 g_AABBsize (rainbowVolumeWidth, rainbowVolumeHeight, rainbowVolumeDepth);
static Cvec3 g_initEyeLocation;
static Cvec3 g_AABBcenter;

shared_ptr<SgRbtNode> g_AABBRbt, g_quadRbt;
SgRbtNode g_SceneRbt;  

// OLD AABB Specification
static const Cvec3 g_initAABBmax(48.8672218, 60.4031639, 35.5040970);
static const Cvec3 g_initAABBmin(-60.4031792, -9.44913101, -60.4031944);
static const Cvec3 g_initAABBcenter = (g_initAABBmin + g_initAABBmax) / 2.0;

// The local coordinates of  the AABB box relative to the AABB coordinate system.
// Only these coordinates are important. 

static const Cvec3 g_localAABBmin = g_initAABBmin - g_initAABBcenter;
// relative to the coordinate system at g_initAABBcenter

static const Cvec3 g_localAABBmax = g_initAABBmax - g_initAABBcenter;


// --------- Scene
// --------- Scene

Cvec4 g_light1Pos, g_light2Pos;  // define two lights positions in world space, which is set right in front of the water
                                 // volume

Cvec4 g_light1Color, g_light2Color; 

Matrix4 g_eyeRbt, g_invEyeRbt, g_cameraRotMat, g_invCameraRotMat;

Cvec3 g_sunRayDir ;

float g_radius = 1.0e-3; // 1 mm, 1.5e-3, 2.03-3

float g_dropDensity = 30000;

Cvec3 g_lightRGBColor;


Picker g_picker; // default constructor

shared_ptr< Object> g_currentPickedObject;

vector < shared_ptr<Object> > g_objectPtrList;

shared_ptr<Geometry> g_geometry;

std::ofstream messageFile("./messageFile.txt", std::ios::out);
std::ofstream textureFile("./textureDebug.txt", std::ios::out);



bool g_redrawWindowEvent = true; // used in ofAppGLFWWindow.cpp
bool g_drawnToBackbuffer = false;

int g_renderMode = 1; // The background scene mode

GLint g_savedFramebuffer;

bool g_Gl2Compatible = false; // use shader programs of type "*-gl3 *"

//--------------------------------------------------------------
ofApp:: ofApp() // default constructor

{
}



void ofApp::setup(){
	// drop down menu settup
	
	//ofBackground(255, 0,0); // for debugging

	ofBackground(128,200,255); // background (buffer clear color) = sky color
     // => glClearColor(128/255, 200/255, 255/255);
	ofSetBackgroundAuto(true); // which is default

	//ofEnableDepthTest();
	//ofSetDepthTest(true);
	//ofEnableBlendMode(OF_BLENDMODE_ALPHA);
	//ofEnableAlphaBlending();
	
	/* This GUI causes the rainbow to have a strange color, perhapse some kind of
	   odd blending seems to happen

    gui = new ofxUICanvas(300,300); // this->gui variable exists

    gui->addLabel("DROPDOWN MENU", OFX_UI_FONT_MEDIUM);
    gui->addSpacer();
    gui->addLabel("'1' TO ADD TO LIST", OFX_UI_FONT_SMALL);
    gui->addLabel("'2' TO DELETE FROM LIST", OFX_UI_FONT_SMALL);
    gui->addLabel("'3' TO DELETE ALL IN LIST", OFX_UI_FONT_SMALL);
    gui->addSpacer();
    vector<string> names;
    names.push_back("render Scene");   
	names.push_back("dump Scene Sysbuffer to Image File");   
	names.push_back("renderToFBO and Dump");  
	names.push_back("dump SceneFBO to Image File");  
	names.push_back("render Rainbow");
	names.push_back("render Rainbow Only");
	names.push_back("render AABBRainbow Only");
	names.push_back("render OldAABBRainbow Only");
	names.push_back("draw Scene from FBO File");
	names.push_back("draw Scene from SysFramebuffer");


    gui->setWidgetFontSize(OFX_UI_FONT_SMALL);
    gui->addToggle("SHOW ACTIVE", false);
    ddl = gui->addDropDownList("DROPDOWN", names); // this->ddl variable exists

    //ddl->setAllowMultiple(true);
    //ddl->setAutoClose(true);
    gui->autoSizeToFitWidgets(); 
    //gui->setDrawWidgetPadding(true);

    ofAddListener(gui->newGUIEvent, this, &ofApp::guiEvent);

	*/
	



	// set up the sun light intensities and the sun direction
	//glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ALPHA );
	setupSun();

	
	// set up the eye position and direction relative to the g_localAABB which is created at
	// initObjects(). Then create the eye Matrix

	// For now, It is assumed that localAABB is centered at the world origin
	// height / horizontalDistFromEyeToAABB  = tan (g_frustFovY / 2.0) 

	setupCamera();
	//setupOldCamera();

	//bool g_Gl2Compatible = false; // use shader programs of type "*-gl3 *"

	// Sets the background clearing function to be auto (default) or not. 
	// If non-auto, then background clearing will not occur per frame (at the start of draw)
	// but rather, whenever ofBackground is called.
	
    //ofSetBackgroundAuto( false); // It will make a single buffering by drawing to the front buffer


	//------------------------------------------------------------
    //void ofAppGLFWWindow::disableSetupScreen(){
	//	bEnableSetupScreen = false;
    //};

	//_description: _

    // Every update/draw cycle, the function ofSetupScreen is called. 
	// That function sets the perspective, coordinate system, and some other openGL parameters.
	// If you need to use your own parameters, the call to that function can be disabled with ofDisableSetupScreen.

	// ofDisableSetupScreen(); SetupScreen() is needed to draw the widgets. But the camera parameters set by
	// this function will be ignored by draw() function defined in this class.

	g_windowWidth = ofGetWidth(); 
	g_windowHeight = ofGetHeight(); 

	updateFrustFovY(); // It uses g_windowWidth and g_windowHeight
	
	//g_windowWidth = nThetas;  
	//g_windowHeight = nSpectralSamples;  

		
	initGLState();

    initMaterials();
	

	initPBRTObjects();

	//initObjects(); // initObjects() is for testing. use it or initPBRTObjects()
	initAABB();

	// shader objects
	initFullScreen(); 
	initNDC();

    initRainbow();
	initAABBRainbow();
	
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
	//spa.month         = 1;
    //spa.day           = 31;
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

    //get_time_and_location(spa.year, spa.month, spa.day, spa.hour, spa.minute, spa.timezone, spa.longitude, spa.latitude);
    
    // calculate sun's position
    result = spa_calculate(&spa);   // input is valid when result eq 0
    
	if ( result !=0 ) {
		messageFile   << "error in the calculation of sun position" << endl;
		cout   << "error in the calculation of sun position" << endl;
		return; 
	}


    // calculate sunRay

    calculate_sunRay(sunRay, spa.zenith, spa.azimuth); // vector sunRay points to the sun
	
	messageFile << "sunRay zenith angle = " << spa.zenith <<" sunRay azimuth" << spa.azimuth << endl;
	
	messageFile << "sunRay in geocentric coord (world coord system) = " << sunRay  << endl;

	cout << "sunRay zenith angle = " << spa.zenith <<" sunRay azimuth" << spa.azimuth << endl;
	
	cout << "sunRay in geocentric coord (world coord system) = " << sunRay  << endl;
	g_sunRayDir = Cvec3(-sunRay[1], sunRay[2], -sunRay[0]); 
	// rename the axes to 3D graphics convention : z (up) => y, x (north) => -z, y(west) => -x
	//                         |
	//                         | x 
	//                   <------
	//                          y
	
	// E.g: In the original coord system: (-10, 20, 5) [ azimth= south-west, polar = positive] =? (-20, 5, 10)
	
    
    messageFile  << "sunRay in Graphics coord = " << g_sunRayDir << endl;
	cout   << "sunRay in Graphics coord = " << g_sunRayDir << endl;

	
}

void ofApp::setupCamera() {

	// set the parameters for the camera, which define the internal working of the camera.
    g_frustMinFov = 90.0;  // old: A minimal of 100 degree field of view for rainbow viewing
	//g_frustMinFov = 60.0;  //new

    g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

	//It means there is very high precision at the near plane, but very little precision at the far plane. 
	//If the range [-n, -f] is getting larger, it causes a depth precision problem (z-fighting); 
	//a small change of ze around the far plane does not affect on zn value. 
	//The distance between n and f should be short as possible to minimize the depth buffer precision problem.

    g_frustNear = -0.01;		// near plane
    g_frustFar = -200.0;		// far plane


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
	
	g_SceneRbt = Matrix4::makeAxisRotation( - PI / 3.0, Cvec3(0,1,0) ); //for NationalModernArt
	
	//g_SceneRbt = Matrix4::makeAxisRotation( 0.0, Cvec3(0,1,0) ); // for Incheon Asian Game

    // set the camera rotation so that it points to the scene to the east of the geocentric frame
	// Ignore the above calculation

	g_cameraRotMat = g_SceneRbt;  
	 

	double eyeHeight = 1.7;
   
	float distanceToEye = 42 ;  // the observer on the ground in Asian Game
	
	Cvec3 eyeLocation = Cvec3(0, eyeHeight, distanceToEye);

	// the position of the eye relative to the world coordinate system.

	g_eyeRbt = g_cameraRotMat * Matrix4::makeTranslation( eyeLocation );
	
	cout << "camera setup: " << endl;
	cout << g_eyeRbt << endl;

	g_invCameraRotMat = inv( g_cameraRotMat);
	messageFile << "g_invRotMat: \n" << g_invCameraRotMat << endl;

	Matrix4 g_invEyeRbt = inv(g_eyeRbt);


	// setup light sources for Asian Game
	//float distanceToLight = 42 + 24 + 5;
	float distanceToLight = 300;
	//float heightToLight = 40 + 15;
	float heightToLight = 30;
	float separationLight = 5;

	g_light1Pos = g_SceneRbt * Cvec4( separationLight/2.0, heightToLight, distanceToLight, 1.0);
	g_light2Pos = g_SceneRbt * Cvec4(-separationLight/2.0, heightToLight, distanceToLight, 1.0);  // define two lights positions in world space

    g_light1Color = Cvec4(1.0, 1.0, 1.0, 1.0);
    g_light2Color = Cvec4(1.0, 1.0, 1.0, 1.0);// set up light sources


	
} // setupCamera()

void ofApp::setupOldCamera() {

	
	
	// set the parameters for the camera, which define the internal working of the camera.
    g_frustMinFov = 100.0;  // A minimal of 100 degree field of view for rainbow viewing

    g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

    g_frustNear = -0.01;		// near plane
    g_frustFar = -200.0;    // far plane


	// set the camera location and direction so that it can see the rainbow well

	Cvec3 yAxis (0,1,0);
	Cvec3 zAxis (0,0,1);

	// set the camera ground direction so that it is equal to the ground projection of the sun ray direction

	Cvec3 groundCamDir = (-g_sunRayDir) - yAxis * dot(-g_sunRayDir, yAxis);
	Cvec3 upCamDir = -g_sunRayDir - groundCamDir;

	Assert ( upCamDir == yAxis * dot(-g_sunRayDir, yAxis) );


	
	groundCamDir.normalize();

	//Matrix4 rotMat;

	Cvec3 negZAxis (0,0,-1);

	if ( norm2( negZAxis - groundCamDir ) < 1.0e-8 ) { // no need to rotate; The rotation matrix will be identity
		g_cameraRotMat =  Matrix4();
	}
	
	else {
		Cvec3 rotAxis = cross( negZAxis, groundCamDir); // rotAxis = zAxis x groundCamDir

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
	    messageFile << Cvec3( g_cameraRotMat * Cvec4( negZAxis, 0.0 ) ) << endl;

		messageFile << "g_rotMat :\n" << g_cameraRotMat << endl;
		


	    Assert ( norm2 ( groundCamDir - Cvec3( g_cameraRotMat * Cvec4( negZAxis, 0.0 ) ) ) <= 1.0e-8 );
	   
	}

    static const Cvec3 g_initAABBmax(48.8672218, 60.4031639, 35.5040970);
    static const Cvec3 g_initAABBmin(-60.4031792, -9.44913101, -60.4031944);
    static const Cvec3 g_initAABBcenter = (g_initAABBmin + g_initAABBmax) / 2.0;


	double eyeHeight = 1.7;
	
	g_invCameraRotMat = inv(g_cameraRotMat);

    Cvec3 localCenter = Cvec3( g_invCameraRotMat * Cvec4(g_initAABBcenter,0) );

	float volumeDepth = g_initAABBmax[2] - g_initAABBmin[2];

	//Cvec3 eyeTranslation = localCenter + Cvec3 ( 0, eyeHeight, volumeDepth );
	Cvec3 eyeTranslation = Cvec3 ( 0, eyeHeight, g_initAABBmax[2] * 2 );
	// the position of the eye relative the world coordinate system.	
	g_eyeRbt = g_cameraRotMat * Matrix4::makeTranslation( eyeTranslation );
	g_invCameraRotMat = inv(g_cameraRotMat);
	messageFile << "g_invRotMat: \n" << g_invCameraRotMat << endl;

	Matrix4 g_invEyeRbt = inv(g_eyeRbt);
} // setupOldCamera()


Cvec3  ofApp::getSunLightRGBColor() {
	// convert the sun light to its RGB representation

	double sunIntensity, X = 0, Y = 0, Z = 0, XYZ;
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

		//  messageFile << "j =" << j  << " i1=" << i1 << "i2=" << i2 << endl;

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

		sunIntensity = calculate_irradiance_of_sun(lambda_m);	// Isun: the return value of irradiance is per nanometer, but the function uses
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
	cout << "g_frustFovY=" << g_frustFovY << endl;
	double aspectRatio =  static_cast <double>(g_windowWidth) / static_cast <double> (g_windowHeight);
	cout << "aspectRatio = " << aspectRatio  << endl;

	return Matrix4::makeProjection(g_frustFovY, g_windowWidth / static_cast <double> (g_windowHeight),
								   g_frustNear, g_frustFar);
}


// update g_frustFovY from g_frustMinFov, g_windowWidth, and g_windowHeight
void ofApp::updateFrustFovY() {

	if (g_windowWidth >= g_windowHeight)

		g_frustFovY = g_frustMinFov;

	else {
		const double RAD_PER_DEG = 0.5 * 3.14159 /180;

		g_frustFovY = atan2(sin(g_frustMinFov * RAD_PER_DEG) * g_windowHeight /static_cast <double> (g_windowWidth), 
							cos(g_frustMinFov * RAD_PER_DEG)) 
					  / RAD_PER_DEG;
	}
}

//--------------------------------------------------------------
void ofApp::update(){

}

static bool backgroundBackedUp = false;


void ofApp::draw() {
	myOwnDraw(); 
}

void ofApp::myOwnDraw() {

    // callback function for draw event: 
	// This callback function is registered at the system initialization, 
	// so called earlier than the other draw callbacks

	// for debugging
	//return;

	//return; // for debugging

	if ( !g_redrawWindowEvent ) return;
		
	cout << "I am here in myOwnDraw() to draw according to the renderMode:" << endl;
	//renderToFBO();

	if ( g_renderMode == 1 ) { // background scene only
		cout << "renderSceneToSysBuffer():"  << endl;
		try {
			renderSceneToSysBuffer();
		}
	
		catch (const runtime_error & error ) {
			std::cout << error.what() << endl;
			cout  <<"renderToSysBuffer:" <<  error.what() << endl;
		}
	}

    if ( g_renderMode == 2 ) { // rainbow only 
		cout << "renderRainbowOnlyToScreen():"  << endl;
		try {
			renderRainbowOnlyToScreen();

		}
	
		catch (const runtime_error & error ) {
			std::cout << "renderRainbowOnlyToScreen: " << error.what() << endl;
	     
		}

	}
	
	 if ( g_renderMode == 3 ) { // background scene + rainbow 
		cout << "renderRainbowAndScreenToScreen():"  << endl;
		try {
			renderRainbowAndSceneToScreen();
		
		}
	
		catch (const runtime_error & error ) {
			std::cout << error.what() << endl;
			cout  <<"renderRainbowAndScreenToScreen():" <<  error.what() << endl;
		}
		
	}


    if ( g_renderMode == 4 ) { // AABBrainbow only 
		cout << "renderAABBRainbowOnlyToScreen():"  << endl;
		try {
			renderAABBRainbowOnlyToScreen();
	
		}
	
		catch (const runtime_error & error ) {
			std::cout << error.what() << endl;
			cout  <<"renderAABBRainbowOnlyToScreen():" <<  error.what() << endl;
		}

	} 

		
	// g_renderMode == 0 => do nothing

	g_redrawWindowEvent = false; //for debugging

	g_drawnToBackbuffer = true;


	//glutSwapBuffers();   show the back buffer (where we rendered stuff)
	// This action will be done in ofAppGLFWWindow.cpp after ofNotifyDraw() which calls
	// ofApp::draw()

} // myOwnDraw()



void ofApp::renderToFBOAndDump(int renderMode ) {
	
	
	if ( renderMode == 1 ) { // scene  
		cout << "renderToFBOAndDump: scene  only:"  << endl;
		try {
			renderSceneToFBO();

			glBindFramebuffer(GL_FRAMEBUFFER, g_fboScene);
  
			checkGlErrors();
			readSceneFBOPixels("sceneFBO.txt");

			glBindFramebuffer(GL_FRAMEBUFFER, g_savedFramebuffer); // default system framebuffer
		}
	
		catch (const runtime_error & error ) {
			std::cout << "renderToFBOAndDump: scene only " << error.what() << endl;
	     
		}

	}
	
	   
    if ( renderMode == 2 ) { // rainbow only 
		cout << "renderToFBOAndDump: rainbow only"  << endl;
		try {
			renderRainbowOnlyToFBO();
			glBindFramebuffer(GL_FRAMEBUFFER, g_fboRainbow);
  
			checkGlErrors();
			readRainbowFBOPixels("rainbowOnlyFBO.txt");
			glBindFramebuffer(GL_FRAMEBUFFER, g_savedFramebuffer); // default system framebuffer
		}
	
		catch (const runtime_error & error ) {
			std::cout << "renderToFBOAndDump: rainbow only " << error.what() << endl;
	     
		}

	}
	
	 if ( renderMode == 3 ) { // background scene + rainbow 
		cout << "renderToFBOAndDump: RainbowAndScreen"  << endl;
		try {
			renderRainbowAndSceneToFBO();

			glBindFramebuffer(GL_FRAMEBUFFER, g_fboRainbow);
  
			checkGlErrors();
			readRainbowFBOPixels("rainbowAndSceneFBO.txt");

			glBindFramebuffer(GL_FRAMEBUFFER, g_savedFramebuffer); // default system framebuffer
			cout << "FBO dumped to rainbowAndSceneFBO.txt" << endl;
		
		}
	
		catch (const runtime_error & error ) {
			std::cout << error.what() << endl;
			cout  <<"renderToFBOAndDump:RainbowAndScreen  "<<  error.what() << endl;
		}

	}


} // renderToFBOAndDump()


void ofApp::renderSceneToFBO() {

	glBindFramebuffer(GL_FRAMEBUFFER, g_fboScene);

	checkGlErrors();
	glEnable(GL_DEPTH_TEST);
	checkGlErrors();
	glDisable(GL_BLEND);
	checkGlErrors();

	// generate and bind a framebuffer
	attachTexturesToSceneFBO();
	
	 
	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ); // for the FBO g_fboScene
	checkGlErrors();
	// GL_INVALID_FRAMEBUFFER_OPERATION is  caused by calling glClear when there's no current draw framebuffer
	// (as will be the case for a windowless OpenGL context before you bind an FBO, for example).
	 

	drawBackgroundStuff();  // draw to FBO g_fboScene

	glBindFramebuffer(GL_FRAMEBUFFER, g_savedFramebuffer);	// unbind the framebuffer and bind it to the default system buffer
	//glBindFramebuffer(GL_FRAMEBUFFER, 0);					// unbind the framebuffer and bind it to the default system buffer
	checkGlErrors();

}//renderToSceneFBO()




void ofApp::renderRainbowOnlyToScreen() {
	

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	
	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ); // for the default framebuffer
	checkGlErrors();
	 
	Uniforms extraUniforms;

	glDrawBuffer(GL_BACK); // GL_BACK is the draw buffer for the zero framebuffer; for double-buffering context
	checkGlErrors();
	
   
	const Matrix4 projMatrix = makeProjectionMatrix();

	// The basic transformations
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


	extraUniforms.put("uWindowWidth", g_windowWidth );
	extraUniforms.put("uWindowHeight",g_windowHeight);
  
	extraUniforms.put("uRadius", g_radius );
	extraUniforms.put("uDropDensity", g_dropDensity );  

	extraUniforms.put("uEyeMatrix", g_eyeRbt); // This will be used in shader to convert the
	extraUniforms.put("uInvEyeMatrix", g_invEyeRbt); // This will be used in shader to convert the
   	   
	cout << "current AABB frame (center):" << "\n" << *g_AABBRbt << endl;

	Matrix4 uInvModelMatrixAABB = inv( *g_AABBRbt );

	
	extraUniforms.put("uModelMatrixAABB",  *g_AABBRbt );
	extraUniforms.put("uInvModelMatrixAABB", uInvModelMatrixAABB );

	Cvec3 uEyeOriginInBoxFrame = Cvec3( uInvModelMatrixAABB * Cvec4(g_eyeRbt[3], g_eyeRbt[7], g_eyeRbt[11], g_eyeRbt[15]) );
	  
	extraUniforms.put("uEyeOriginInBoxFrame", uEyeOriginInBoxFrame);

	g_AABBmin = Cvec3( -g_AABBsize[0]/2.0, -g_AABBsize[1] / 2.0, -g_AABBsize[2]/2.0 );
	g_AABBmax = Cvec3( g_AABBsize[0]/2.0, g_AABBsize[1] / 2.0, g_AABBsize[2]/2.0 );
	   
	//g_AABBmin = Cvec3( *g_AABBRbt * Cvec4(  -g_AABBsize[0]/2.0, -g_AABBsize[1] / 2.0, -g_AABBsize[2]/2.0, 1.0 ) );
	//g_AABBmax = Cvec3( *g_AABBRbt  * Cvec4(  g_AABBsize[0]/2.0, g_AABBsize[1] / 2.0, g_AABBsize[2]/2.0, 1.0 )  );
	   
	extraUniforms.put("uAABBmin", g_AABBmin  );   
	extraUniforms.put("uAABBmax", g_AABBmax );   
	
	Cvec4 localFrontNormal = Cvec4( 0.0, 0.0, 1.0, 0.0);
	Cvec3 frontNormalToAABB = Cvec3( g_invEyeRbt * (*g_AABBRbt) * localFrontNormal );


	cout << "frontNormalToAABB:=" << localFrontNormal <<" " << (*g_AABBRbt) * localFrontNormal  << " " << frontNormalToAABB << endl;

	// The orientation of g_AABBRbt and that of g_EyeRbt are the same. 

	// draw rainbow to the default buffer
	try {

		//Matrix4 MVM = g_invEyeRbt * (  *( g_rainbowOnlyMatObj -> objectRbt ) ) ;		// g_currentPickedObject->objectRbt may have been
		//																			// changed by picking
		Matrix4 MVM = Matrix4();

		extraUniforms.put("uModelViewMatrix", MVM).put("uNormalMatrix", normalMatrix(MVM) );
		// debugging

		g_rainbowOnlyMatObj->draw( extraUniforms ); 

		// This draw leads to  BufferedObjectGeometry::draw() which binds vertex buffers 
		// and set vertex     attribute pointers, and then calls drawElements() 
	  
	}                                   


	catch ( const runtime_error & error ) {
		//std::cout << error.what() << endl;
		messageFile << "uModelViewMatrix in rainbow" <<error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	}

	
} //  renderRainbowOnlyToScreen() 



void ofApp::renderAABBRainbowOnlyToScreen() {
	

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	
	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ); // for the default framebuffer
	checkGlErrors();
	 
	Uniforms extraUniforms;

	glDrawBuffer(GL_BACK); // GL_BACK is the draw buffer for the zero framebuffer; for double-buffering context
	checkGlErrors();
	
   
	const Matrix4 projMatrix = makeProjectionMatrix();

	// The basic transformations
	extraUniforms.put("uProjectionMatrix", projMatrix);
   
	g_invEyeRbt = inv(g_eyeRbt);

	cout << "At renderAABBRainbow: current eye frame:" << "\n" << g_eyeRbt << endl;

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


	extraUniforms.put("uWindowWidth", g_windowWidth );
	extraUniforms.put("uWindowHeight",g_windowHeight);
  
	extraUniforms.put("uRadius", g_radius );
	extraUniforms.put("uDropDensity", g_dropDensity );  

	extraUniforms.put("uEyeMatrix", g_eyeRbt); // This will be used in shader to convert the
	extraUniforms.put("uInvEyeMatrix", g_invEyeRbt); // This will be used in shader to convert the
   	   
	cout << "At renderAABBRainbow: current AABB frame (center):" << "\n" << *g_AABBRbt << endl;

	Matrix4 uInvModelMatrixAABB = inv( *g_AABBRbt );

	

	extraUniforms.put("uInvModelMatrixAABB", uInvModelMatrixAABB );

	Cvec3 uEyeOriginInBoxFrame = Cvec3( uInvModelMatrixAABB * Cvec4(g_eyeRbt[3], g_eyeRbt[7], g_eyeRbt[11], g_eyeRbt[15]) );
	  
	extraUniforms.put("uEyeOriginInBoxFrame", uEyeOriginInBoxFrame);

	g_AABBmin = Cvec3( -g_AABBsize[0]/2.0, -g_AABBsize[1] / 2.0, -g_AABBsize[2]/2.0 );
	g_AABBmax = Cvec3( g_AABBsize[0]/2.0, g_AABBsize[1] / 2.0, g_AABBsize[2]/2.0 );
	   
	//g_AABBmin = Cvec3( *g_AABBRbt * Cvec4(  -g_AABBsize[0]/2.0, -g_AABBsize[1] / 2.0, -g_AABBsize[2]/2.0, 1.0 ) );
	//g_AABBmax = Cvec3( *g_AABBRbt  * Cvec4(  g_AABBsize[0]/2.0, g_AABBsize[1] / 2.0, g_AABBsize[2]/2.0, 1.0 )  );
	   
	extraUniforms.put("uAABBmin", g_AABBmin );   
	extraUniforms.put("uAABBmax", g_AABBmax );   
	
	Cvec4 localFrontNormal = Cvec4( 0.0, 0.0, 1.0, 0.0);
	Cvec3 frontNormalToAABB = Cvec3( g_invEyeRbt * (*g_AABBRbt) * localFrontNormal );


	cout << "frontNormalToAABB:=" << localFrontNormal <<" " << (*g_AABBRbt) * localFrontNormal  << " " << frontNormalToAABB << endl;

	// The orientation of g_AABBRbt and that of g_EyeRbt are the same. 

	// draw rainbow to the default buffer
	try {

		Matrix4 MVM = g_invEyeRbt * (  *( g_AABBRainbowOnlyMatObj -> objectRbt ) ) ;		// g_currentPickedObject->objectRbt may have been
																						// changed by picking
  
		extraUniforms.put("uModelViewMatrix", MVM).put("uNormalMatrix", normalMatrix(MVM) );
   
		g_AABBRainbowOnlyMatObj->draw( extraUniforms ); 

		// This draw leads to  BufferedObjectGeometry::draw() which binds vertex buffers 
		// and set vertex     attribute pointers, and then calls drawElements() 
	 
     
	  
	}                                   


	catch ( const runtime_error & error ) {
		//std::cout << error.what() << endl;
		messageFile << "uModelViewMatrix in rainbow" <<error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	}

	
} //  renderAABBRainbowOnlyToScreen() 




void ofApp::renderRainbowAndSceneToScreen( ) {

	renderSceneToFBO();


	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
		
	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	checkGlErrors();
	 
	glDrawBuffer(GL_BACK); // GL_BACK is the draw buffer for the zero framebuffer; for double-buffering context
	checkGlErrors();
	
   

	Uniforms extraUniforms;

  
	const Matrix4 projMatrix = makeProjectionMatrix();

	// The basic transformations
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

	extraUniforms.put("uSunRayDir", Cvec3( g_invEyeRbt * Cvec4(g_sunRayDir,0) ));
	// extraUniforms.put("uSunRayDir", g_sunRayDir );


	extraUniforms.put("uEyeMatrix", g_eyeRbt); // This will be used in shader to convert the
	extraUniforms.put("uInvEyeMatrix", g_invEyeRbt); // This will be used in shader to convert the
   	   
	cout << "current AABB frame (center):" << "\n" << *g_AABBRbt << endl;

	Matrix4 uInvModelMatrixAABB = inv( *g_AABBRbt );

	
	extraUniforms.put("uModelMatrixAABB",  *g_AABBRbt );
	extraUniforms.put("uInvModelMatrixAABB", uInvModelMatrixAABB );

	Cvec3 uEyeOriginInBoxFrame = Cvec3( uInvModelMatrixAABB * Cvec4(g_eyeRbt[3], g_eyeRbt[7], g_eyeRbt[11], g_eyeRbt[15]) );
	  
	extraUniforms.put("uEyeOriginInBoxFrame", uEyeOriginInBoxFrame);

	g_AABBmin = Cvec3( -g_AABBsize[0]/2.0, -g_AABBsize[1] / 2.0, -g_AABBsize[2]/2.0);
	g_AABBmax = Cvec3( g_AABBsize[0]/2.0, g_AABBsize[1] / 2.0, g_AABBsize[2]/2.0 );
	   

	extraUniforms.put("uAABBmin", g_AABBmin  );   
	extraUniforms.put("uAABBmax", g_AABBmax );   
	
	Cvec4 localFrontNormal = Cvec4( 0.0, 0.0, 1.0, 0.0);
	Cvec3 frontNormalToAABB = Cvec3( g_invEyeRbt * (*g_AABBRbt) * localFrontNormal );

	cout << "frontNormalToAABB:=" << localFrontNormal <<" " << (*g_AABBRbt) * localFrontNormal  << " " << frontNormalToAABB << endl;

	// bind the FBO textures

	extraUniforms.put("uColorTex", 
			shared_ptr<Texture_FROM_FBO>(new Texture_FROM_FBO( g_colorTex) ));	// g_colorTex is the bound texture
																				// which was attached to g_fboScene
	extraUniforms.put("uDepthTex", 
			shared_ptr<Texture_FROM_FBO>(new Texture_FROM_FBO(g_depthTex)));


	extraUniforms.put("uRadius", g_radius );
	extraUniforms.put("uDropDensity", g_dropDensity );   

	

	// draw rainbow to the default buffer
	try {
      
		//Matrix4 MVM = g_invEyeRbt * (  *( g_rainbowAndSceneMatObj -> objectRbt ) ) ; 
	   
		// for debugging: try g_fullScreenMatObj
		//Matrix4 MVM = g_invEyeRbt * (  *( g_fullScreenMatObj-> objectRbt ) ) ;	// g_currentPickedObject->objectRbt may have been
		//																		// changed by picking
		Matrix4 MVM = Matrix4();

		extraUniforms.put("uModelViewMatrix", MVM).put("uNormalMatrix", normalMatrix(MVM) );
   
		g_rainbowAndSceneMatObj->draw( extraUniforms ); 
		//g_fullScreenMatObj->draw( extraUniforms ); 

		// This draw leads to  BufferedObjectGeometry::draw() which binds vertex buffers 
		// and set vertex attribute pointers, and then calls drawElements() 
	 
     
	  
	}                                   


	catch ( const runtime_error & error ) {
		//std::cout << error.what() << endl;
		messageFile << "uModelViewMatrix in rainbow" <<error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	}

		
} //  renderRainbowAndSceenToScreen() 

void ofApp::renderRainbowOnlyToFBO() {
	
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	
	// render to  FBO g_fboRainbow 

	glBindFramebuffer(GL_FRAMEBUFFER, g_fboRainbow);
	// attach new textures to g_fboRainbow
	attachTexturesToRainbowFBO(); // attach textures g_colorTex1, g_dataTex1, and g_depthTex1, to g_fboRainbow


	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ); 
	checkGlErrors();
	 
    Uniforms extraUniforms;
     
    const Matrix4 projMatrix = makeProjectionMatrix();

	 // The basic transformations
    extraUniforms.put("uProjectionMatrix", projMatrix);   

    cout << "current eye frame:" << "\n" << g_eyeRbt << endl;

	 g_invEyeRbt = inv(g_eyeRbt);
   
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

	extraUniforms.put("uRadius", g_radius );
    extraUniforms.put("uDropDensity", g_dropDensity );   
   
	
    extraUniforms.put("uEyeMatrix", g_eyeRbt); // This will be used in shader to convert the
    extraUniforms.put("uInvEyeMatrix", g_invEyeRbt); // This will be used in shader to convert the
   	   
	cout << "current AABB frame (center):" << "\n" << *g_AABBRbt << endl;

	Matrix4 uInvModelMatrixAABB = inv( *g_AABBRbt );
	
	extraUniforms.put("uModelMatrixAABB",  *g_AABBRbt );
	extraUniforms.put("uInvModelMatrixAABB", uInvModelMatrixAABB );

	Cvec3 uEyeOriginInBoxFrame = Cvec3( uInvModelMatrixAABB * Cvec4(g_eyeRbt[3], g_eyeRbt[7], g_eyeRbt[11], g_eyeRbt[15]) );
	  
	extraUniforms.put("uEyeOriginInBoxFrame", uEyeOriginInBoxFrame);

	//g_AABBmin = Cvec3( *g_AABBRbt * Cvec4(  -g_AABBsize[0]/2.0, -g_AABBsize[1] / 2.0, -g_AABBsize[2]/2.0, 1.0 ) );
	//g_AABBmax = Cvec3( *g_AABBRbt  * Cvec4(  g_AABBsize[0]/2.0, g_AABBsize[1] / 2.0, g_AABBsize[2]/2.0, 1.0 )  );
	   
	g_AABBmin = Cvec3( -g_AABBsize[0]/2.0, -g_AABBsize[1] / 2.0, -g_AABBsize[2]/2.0 );
	g_AABBmax = Cvec3( g_AABBsize[0]/2.0, g_AABBsize[1] / 2.0, g_AABBsize[2]/2.0 );
	   

	extraUniforms.put("uAABBmin", g_AABBmin);   
	extraUniforms.put("uAABBmax", g_AABBmax);   

	Cvec4 localFrontNormal = Cvec4( 0.0, 0.0, 1.0, 0.0);
	Cvec3 frontNormalToAABB = Cvec3( g_invEyeRbt * (*g_AABBRbt) * localFrontNormal );


	cout << "frontNormalToAABB:=" << localFrontNormal <<" " << (*g_AABBRbt) * localFrontNormal  << " " << frontNormalToAABB << endl;

	
	// draw rainbow to the default buffer
	try {

		//Matrix4 MVM = g_invEyeRbt * (  *( g_rainbowOnlyMatObj -> objectRbt ) ) ;		// g_currentPickedObject->objectRbt may have been
		//																			//changed by picking
		Matrix4 MVM = Matrix4();
		extraUniforms.put("uModelViewMatrix", MVM).put("uNormalMatrix", normalMatrix(MVM) );
   
		g_rainbowOnlyMatObj->draw( extraUniforms ); 

		// This draw leads to  BufferedObjectGeometry::draw() which binds vertex buffers 
		// and set vertex     attribute pointers, and then calls drawElements() 
	 
     
	  
	}                                   


	catch ( const runtime_error & error ) {
		//std::cout << error.what() << endl;
		messageFile << "uModelViewMatrix in rainbow" <<error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	}

	 glBindFramebuffer(GL_FRAMEBUFFER, g_savedFramebuffer);
	
} //  renderRainbowOnlyToFBO() 




void ofApp::renderRainbowAndSceneToFBO() {
	
	// render the background scene to FBO

	renderSceneToFBO(); // render to FBO g_fboScene which contains textures g_colorTex, g_dataTex, and g_depthTex


	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	
	// render to another FBO g_fboRainbow using the textures g_colorTex and g_depthTex

	glBindFramebuffer(GL_FRAMEBUFFER, g_fboRainbow);
	// attach new textures to g_fboRainbow
	attachTexturesToRainbowFBO(); // attach textures g_colorTex1, g_dataTex1, and g_depthTex1, to g_fboRainbow


	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ); 
	checkGlErrors();
	 
    Uniforms extraUniforms;
     
    g_projectionMatrix = makeProjectionMatrix();

	// The basic transformations
    extraUniforms.put("uProjectionMatrix", g_projectionMatrix);   

    cout << "current eye frame:" << "\n" << g_eyeRbt << endl;

	g_invEyeRbt = inv(g_eyeRbt);
   
    const Cvec3 eyeLight1 = Cvec3(g_invEyeRbt * g_light1Pos ); // g_light1 position in eye coordinates
    const Cvec3 eyeLight2 = Cvec3(g_invEyeRbt * g_light2Pos ); // g_light2 position in eye coordinates
  
	// send the eye space coordinates of lights to uniforms
    extraUniforms.put("uLight1Pos", eyeLight1);
    extraUniforms.put("uLight2Pos", eyeLight2);

    extraUniforms.put("uLight1Color", g_light1Color );
    extraUniforms.put("uLight2Color", g_light2Color );

	// bind the FBO textures

	extraUniforms.put("uColorTex", 
		  shared_ptr<Texture_FROM_FBO>(new Texture_FROM_FBO( g_colorTex) )); // g_colorTex is the bound texture
	                                                                             // which was attached to g_fboScene
	extraUniforms.put("uDepthTex", 
		   shared_ptr<Texture_FROM_FBO>(new Texture_FROM_FBO(g_depthTex)));

	// upload the sun direction relative to the eye space

    extraUniforms.put("uSunRayDir", Cvec3( g_invEyeRbt * Cvec4(g_sunRayDir,0) ) );
	//extraUniforms.put("uSunRayDir", g_sunRayDir );

	extraUniforms.put("uRadius", g_radius );
    extraUniforms.put("uDropDensity", g_dropDensity );   
   
	
    extraUniforms.put("uEyeMatrix", g_eyeRbt); // This will be used in shader to convert the
    extraUniforms.put("uInvEyeMatrix", g_invEyeRbt); // This will be used in shader to convert the
   	   
	cout << "current AABB frame (center):" << "\n" << *g_AABBRbt << endl;

	Matrix4 uInvModelMatrixAABB = inv( *g_AABBRbt );
	
	extraUniforms.put("uModelMatrixAABB",  *g_AABBRbt );
	extraUniforms.put("uInvModelMatrixAABB", uInvModelMatrixAABB );

	Cvec3 uEyeOriginInBoxFrame = Cvec3( uInvModelMatrixAABB * Cvec4(g_eyeRbt[3], g_eyeRbt[7], g_eyeRbt[11], g_eyeRbt[15]) );
	  
	extraUniforms.put("uEyeOriginInBoxFrame", uEyeOriginInBoxFrame);

	g_AABBmin = Cvec3( -g_AABBsize[0]/2.0, -g_AABBsize[1] / 2.0, -g_AABBsize[2]/2.0 );
	g_AABBmax = Cvec3( g_AABBsize[0]/2.0, g_AABBsize[1] / 2.0, g_AABBsize[2]/2.0 );
	   
	//g_AABBmin = Cvec3( *g_AABBRbt * Cvec4(  -g_AABBsize[0]/2.0, -g_AABBsize[1] / 2.0, -g_AABBsize[2]/2.0, 1.0 ) );
	//g_AABBmax = Cvec3( *g_AABBRbt  * Cvec4(  g_AABBsize[0]/2.0, g_AABBsize[1] / 2.0, g_AABBsize[2]/2.0, 1.0 )  );
	   
	extraUniforms.put("uAABBmin", g_AABBmin);   
	extraUniforms.put("uAABBmax", g_AABBmax);   

	Cvec4 localFrontNormal = Cvec4(0.0, 0.0, 1.0, 0.0);
	Cvec3 frontNormalToAABB = Cvec3( g_invEyeRbt * (*g_AABBRbt) * localFrontNormal );


	cout << "frontNormalToAABB:=" << localFrontNormal <<" " << (*g_AABBRbt) * localFrontNormal  << " " << frontNormalToAABB << endl;

	// draw rainbow to the default buffer
	try {

		//Matrix4 MVM = g_invEyeRbt * (  *( g_rainbowAndSceneMatObj -> objectRbt ) ) ;		// g_currentPickedObject->objectRbt may have been
		//																				// changed by picking

		Matrix4 MVM = Matrix4();
  
		extraUniforms.put("uModelViewMatrix", MVM).put("uNormalMatrix", normalMatrix(MVM) );
   
		g_rainbowAndSceneMatObj->draw( extraUniforms ); 

		// This draw leads to  BufferedObjectGeometry::draw() which binds vertex buffers 
		// and set vertex attribute pointers, and then calls drawElements() 
	 
     
	  
	}                                   


	catch ( const runtime_error & error ) {
		//std::cout << error.what() << endl;
		messageFile << "uModelViewMatrix in rainbow" <<error.what() << endl;
		cout << error.what() << endl;

		//throw; // A throw expression that has no operand re-throws the exception currently being handled

	}

	glBindFramebuffer(GL_FRAMEBUFFER, g_savedFramebuffer);

} //  renderRainbowAndSceneToFBO() 



void ofApp::renderSceneToSysBuffer() {
	//GLint window_fbo =0; 
	//glGetIntegerv(GL_FRAMEBUFFER_BINDING, &window_fbo); 
	//glBindFramebuffer(GL_FRAMEBUFFER, window_fbo);
	//checkGlErrors();
	// GL_INVALID_OPERATION is the error you get when multiple combinations of parameters that depend on different state are in conflict.
	// If it were just a missing enum, you should get GL_INVALID_ENUM.

	// GL_INVALID_OPERATION is generated if framebuffer is not zero or the name of a framebuffer previously returned 
	// from a call to glGenFramebuffers. 

	/*
	I feel it's necessary to point out here that the call to glBindFramebuffer(GL_FRAMEBUFFER, 0); 
		does not return rendering to the main framebuffer
		although it would appear to work for machines that run Windows, Unix(Mac) or Linux.
		Desktops and laptops have no concept of a main default system buffer. 
		This idea started with handheld devices. 
		When you make an openGL bind call with zero as the parameter
		then what you are doing is setting this function to NULL.
		It's how you disable this function. It's the same with glBindTexture(GL_TEXTURE_2D, 0);
	It is possible that on some handheld devices that the driver automatically activates
		the main system framebuffer when you set the framebuffer to NULL without activating another.
		This would be a choice made by the manufacturer and is not something that you should count on,
		this is not part of the openGL ES spec.
		For desktops and laptops, this is absolutely necessary 
		since disabling the framebuffer is required to return to normal openGL rendering. 

        On an iOS device, you should make the following call, glBindFramebuffer(GL_FRAMEBUFFER, viewFramebuffer);,
		providing that you named your system framebuffer 'viewFramebuffer'.
		Look for through your initialization code for the following call, 
		glGenFramebuffers(1, &viewFramebuffer); 
	    Whatever you have written at the end there is what you bind to when returning to your main system buffer.


        If you are using GLKit then you can use the following call, [((GLKView *) self.view) bindDrawable];
		The 'self.view' may be slightly different depending on your particular startup code.
        Also, for iOS, you could use, glBindFramebuffer(GL_FRAMEBUFFER, 2); 
		but this is likely not going to be consistent across future devices released by Apple. 
		They may change the default value of '2' to be '3' or something else in the future 
		so you'd want to use the actual name instead of an integer value.
	*/

	//glClearDepth(0.);

	//glClearDepth(1.0); // the default value

	//glDepthFunc(GL_GREATER);

	checkGlErrors();

	glEnable(GL_DEPTH_TEST);
	checkGlErrors();
	glDisable(GL_BLEND);
	checkGlErrors();

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	checkGlErrors();

	drawBackgroundStuff();

}//renderToSysBuffer()


void ofApp::drawBackgroundStuff()  {

	//glEnable(GL_DEPTH_TEST);
	//checkGlErrors();
	
	//glDisable(GL_BLEND);
	//checkGlErrors();
	

	Uniforms extraUniforms;
	// build & send proj. matrix to vshader
  
	const Matrix4 projMatrix = makeProjectionMatrix();
	g_projectionMatrix =  projMatrix; 

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

	extraUniforms.put("uWindowWidth", g_windowWidth);
	extraUniforms.put("uWindowHeight",g_windowHeight);
	

  	  
	// (1) Draw the background objects to  FBO =  a container for textures and an optional depth buffer
	// FBO: there are no visible color buffer bitplanes, only a single "off-screen" color image attachment,
	// so there is no sense of front and back buffers or SWAPPING.   
 
	// draw all the objects except for the last, which is a rainbow
    
	for ( int i  = 0; i < g_objectPtrList.size(); i++ ) {	// (5) This loop goes over the list of objects g_objectPtrList and draw each object in the list.
															//     Where and how is this object list created?

		Matrix4 MVM = g_invEyeRbt * (  *( g_objectPtrList[i] -> objectRbt ) ) ; // g_currentPickedObject->objectRbt may have been
																				// changed by picking
		try {
			extraUniforms.put("uModelViewMatrix", MVM).put("uNormalMatrix", normalMatrix(MVM) );
		}
		catch ( const runtime_error & error ) {
			//std::cout << error.what() << endl;
			messageFile << "uModelViewMatrix:" << error.what() << endl;
			cout << error.what() << endl;			
		}

		// compute the depth of the vertcies in eye space for debugging

     
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
			messageFile << "object->draw():" <<  error.what() << endl;
			cout << error.what() << endl;

			//throw; // A throw expression that has no operand re-throws the exception currently being handled

		}

	 
	} // for each object except for rainbow


} // drawBackgroundStuff() for ordinary rendering



void ofApp::drawForPicking() {
 
   
	drawPseudoColors();
  
  
	glFinish();
	checkGlErrors();
}

void ofApp::drawPseudoColors() {

	glEnable(GL_DEPTH_TEST);
	checkGlErrors();
	glDisable(GL_BLEND);
	checkGlErrors();
	glClearColor(0, 0, 0, 0);	// this background color is only meant to be used
								// for the background of the back buffer for rendering in picking.
	checkGlErrors();
	//GLint window_fbo =0; 
	// glGetIntegerv(GL_FRAMEBUFFER_BINDING, &window_fbo);  
	//glBindFramebuffer(GL_FRAMEBUFFER, window_fbo); // 0 = the system-provided buffer
  

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // CLEAR THE COLOR AND DEPTH BUFFER => does it apply to FBOs too?
	checkGlErrors();

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

	for ( int i  = 0; i <g_objectPtrList.size() ; i++ ) {// (5) This loop goes over the list of objects g_objectPtrList and draw each object in the list.
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

		extraUniforms.put( "uMaterialColor", pickColor );	// set material color. In the case of
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

void ofApp::generateFBO() {
	glGenFramebuffers(1, &g_fboScene);
	checkGlErrors();
	// glBindFramebuffer(GL_FRAMEBUFFER, g_fboScene);  
	
	glGenFramebuffers(1, &g_fboRainbow);
	checkGlErrors();
}


void ofApp::attachTexturesToSceneFBO() {

	// framebuffer object with value zero is reserved to represent the default framebuffer 
    // provided by the windowing system.
    // FBO uses a nonzero framebuffer object fbo:
    //GLuint fbo;
	 
	//glBindFramebuffer(GL_FRAMEBUFFER, g_fboScene);  
	//checkGlErrors();
	

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
	  
	 
	glViewport(0,0,g_windowWidth, g_windowHeight);   
	checkGlErrors();
	 
	 	  
	// The texture we're going to render colors to 

	//GLuint colorTex;
	glGenTextures(1, &g_colorTex);
	checkGlErrors();
	// "Bind" the generated texture as the current texture: 
	// all future texture functions will modify this current texture 

	glBindTexture(GL_TEXTURE_2D, g_colorTex);

	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	checkGlErrors();
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	checkGlErrors();
 
	//Components of texels are named after color components. 
	//Because this parameter does not end in “_INTEGER”, OpenGL knows that the data we are uploading is 
	//either a floating-point value or a normalized integer value (which converts to a float when accessed by the shader).
	//The parameter GL_UNSIGNED_BYTE says that each component that we are uploading is stored in an 8-bit unsigned byte. 
	//glTexImage2D(GL_TEXTURE_2D, 0,  (!srgb) || g_Gl2Compatible ? GL_RGB : GL_SRGB, tga.imageWidth, tga.imageHeight,
	//0, GL_RGBA, GL_UNSIGNED_BYTE, tga.imageData);

	bool srgb = true;

	// do not use SRGB format because we are doing light calculation

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, g_windowWidth, g_windowHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);

	//  glTexImage2D(GL_TEXTURE_2D, 0, (!srgb) || g_Gl2Compatible ? GL_RGBA : GL_SRGB, 
	//	            g_windowWidth, g_windowHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	checkGlErrors();

	// attach the texture to the framebuffer
	glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, g_colorTex, 0); // 0 = mipmap level  
	checkGlErrors();

	glBindTexture(GL_TEXTURE_2D,0);  // why unbind the texture object?
	checkGlErrors();
	// depth texture
	 	

	//GLuint dataTex;
	glGenTextures(1, &g_dataTex);
	checkGlErrors();
	// "Bind" the generated texture as the current texture: 
	// all future texture functions will modify this current texture 

	glBindTexture(GL_TEXTURE_2D, g_dataTex);

	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	checkGlErrors();
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	checkGlErrors();
 
	//Components of texels are named after color components. 
	//Because this parameter does not end in “_INTEGER”, OpenGL knows that the data we are uploading is 
	//either a floating-point value or a normalized integer value (which converts to a float when accessed by the shader).
	//The parameter GL_UNSIGNED_BYTE says that each component that we are uploading is stored in an 8-bit unsigned byte. 

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, g_windowWidth, g_windowHeight, 0, GL_RGBA,
	GL_FLOAT, NULL);
	checkGlErrors();

	// attach the texture to the framebuffer
	glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, g_dataTex, 0); // 0 = mipmap level  
	checkGlErrors();

	glBindTexture(GL_TEXTURE_2D,0);  // why unbind the texture object?
	checkGlErrors();
	// depth texture

	//GLuint depthTex;
	glGenTextures(1, &g_depthTex);
	checkGlErrors();
	glBindTexture(GL_TEXTURE_2D, g_depthTex);
	checkGlErrors();
		
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, g_windowWidth, g_windowHeight, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	checkGlErrors();

	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
	
	// attach the depth texture to the framebuffer
	glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, g_depthTex, 0);
	checkGlErrors();
	glBindTexture(GL_TEXTURE_2D,0); 
	checkGlErrors();

	//GLuint attachments[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
	GLuint attachments[2] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
	glDrawBuffers(2, attachments);		// drawing color and depth all together? - need to check
	checkGlErrors();
	//DrawBuffer is per-framebuffer state, and will default to COLOR_ATTACHMENT0 for non-zero FBOs.
	 
	// the buffer selection in glDrawBuffers is part of the framebuffer object state. 
	// Therefore, if this setting is constant, this function can be called only once when creating the framebuffer.
	// glDrawBuffers defines an array of buffers into which outputs from the fragment shader data will be written
	// If a fragment shader writes a value to one or more user defined output variables, then the value of each variable
	// will be written into the buffer specified at a location within bufs
	// bufs \{ GL_FRONT_LEFT, ... GL_COLOR_ATTACHMENTn}

	

	GLenum e = glCheckFramebufferStatus(GL_FRAMEBUFFER);
		checkGlErrors();
	if (e != GL_FRAMEBUFFER_COMPLETE)
		cout << "There is a problem with the FBO\n" << endl;

	//Now that FBO is completely specified, clear the buffer with  the original clear color in order to draw the real scene
      
	//You should use glCZKlear (...), modern GPUs use color buffer, depth buffer and stencil
	//buffer compression. Clearing the buffer is extremely cheap because these buffers are 
	//often hierarchically tiled, and the process of clearing the buffer amounts to flipping 
	//one or two bits in each tile. It even helps with many early fragment tests and general
	//frame buffer throughput if you regularly clear the buffers. On Tile-Based Deferred 
	//Rendering GPUs (e.g. PowerVR SGX - all iOS devices) it is equally important for similar reasons. 


	 
} // attachTexturesToSceneFBO


void ofApp::attachTexturesToRainbowFBO() {

	// framebuffer object with value zero is reserved to represent the default framebuffer 
	// provided by the windowing system.
	// FBO uses a nonzero framebuffer object fbo:
	//GLuint fbo;
	 
	//glBindFramebuffer(GL_FRAMEBUFFER, g_fboScene);  
	//checkGlErrors();
	

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
	GLint maxColorAttachments;

	glGetIntegerv( GL_MAX_COLOR_ATTACHMENTS, &maxColorAttachments);

	cout << "maxColorAttachments =" << maxColorAttachments << endl;
	  	 
	glViewport(0,0,g_windowWidth, g_windowHeight);   
	checkGlErrors();
	 
	 	  
	// The texture we're going to render colors to 

	//GLuint colorTex;
	glGenTextures(1, &g_colorTex1);
	checkGlErrors();
	// "Bind" the generated texture as the current texture: 
	// all future texture functions will modify this current texture 

	glBindTexture(GL_TEXTURE_2D, g_colorTex1);

	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	checkGlErrors();
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	checkGlErrors();
 
	bool srgb = true;

	// do not use SRGB format because we are doing light calculation

	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, g_windowWidth, g_windowHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);

	//glTexImage2D(GL_TEXTURE_2D, 0, (!srgb) || g_Gl2Compatible ? GL_RGBA : GL_SRGB, g_windowWidth, g_windowHeight, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	checkGlErrors();

	// attach the texture to the framebuffer
	glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, g_colorTex1, 0); // 0 = mipmap level  
	checkGlErrors();

	glBindTexture(GL_TEXTURE_2D,0);  // why unbind the texture object?
	checkGlErrors();

	  
	// depth texture
	 	
	
	//GLuint depthTex;
	glGenTextures(1, &g_depthTex1);
	checkGlErrors();
	glBindTexture(GL_TEXTURE_2D, g_depthTex1);
	checkGlErrors();
		
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, g_windowWidth, g_windowHeight, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	checkGlErrors();

	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_FUNC, GL_LEQUAL);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
	
	// attach the depth texture to the framebuffer
	glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, g_depthTex1, 0);
	checkGlErrors();
	glBindTexture(GL_TEXTURE_2D,0); 
	checkGlErrors();



	//GLuint dataTex;
	glGenTextures(1, &g_dataTex1);
	checkGlErrors();
	// "Bind" the generated texture as the current texture: 
	// all future texture functions will modify this current texture 

	glBindTexture(GL_TEXTURE_2D, g_dataTex1);

	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	checkGlErrors();
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	checkGlErrors();
 
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, g_windowWidth, g_windowHeight, 0, GL_RGBA,
	GL_FLOAT, NULL);
	checkGlErrors();

	// attach the texture to the framebuffer
	glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1,g_dataTex1, 0); // 0 = mipmap level  
	checkGlErrors();

	glBindTexture(GL_TEXTURE_2D,0);  // why unbind the texture object?
	checkGlErrors();


	  
	  
	//GLuint dataTex;
	glGenTextures(1, &g_dataTex2);
	checkGlErrors();
	// "Bind" the generated texture as the current texture: 
	// all future texture functions will modify this current texture 

	glBindTexture(GL_TEXTURE_2D, g_dataTex2);

	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	checkGlErrors();
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	checkGlErrors();
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	checkGlErrors();
 
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, g_windowWidth, g_windowHeight, 0, GL_RGBA,
	GL_FLOAT, NULL);
	checkGlErrors();

	// attach the texture to the framebuffer
	glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2,g_dataTex2, 0); // 0 = mipmap level  
	checkGlErrors();

	glBindTexture(GL_TEXTURE_2D,0);  // why unbind the texture object?
	checkGlErrors();


	GLuint attachments[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
	//GLuint attachments[2] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
	//,  Max=8: GL_COLOR_ATTACHMENT1,GL_COLOR_ATTACHMENT1,    GL_COLOR_ATTACHMENT1,  
	//	GL_COLOR_ATTACHMENT5,  GL_COLOR_ATTACHMENT6,  GL_COLOR_ATTACHMENT7};
	  
	glDrawBuffers(3, attachments);		// drawing color and depth all together? - need to check
	checkGlErrors();
	//DrawBuffer is per-framebuffer state, and will default to COLOR_ATTACHMENT0 for non-zero FBOs.
	 
	// the buffer selection in glDrawBuffers is part of the framebuffer object state. 
	// Therefore, if this setting is constant, this function can be called only once when creating the framebuffer.
	// glDrawBuffers defines an array of buffers into which outputs from the fragment shader data will be written
	// If a fragment shader writes a value to one or more user defined output variables, then the value of each variable
	// will be written into the buffer specified at a location within bufs
	// bufs \{ GL_FRONT_LEFT, ... GL_COLOR_ATTACHMENTn}

	

	GLenum e = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	checkGlErrors();
	if (e != GL_FRAMEBUFFER_COMPLETE)
		cout << "There is a problem with the FBO\n" << endl;

	//Now that FBO is completely specified, clear the buffer with  the original clear color in order to draw the real scene
      
	//You should use glClear (...), modern GPUs use color buffer, depth buffer and stencil
	//buffer compression. Clearing the buffer is extremely cheap because these buffers are 
	//often hierarchically tiled, and the process of clearing the buffer amounts to flipping 
	//one or two bits in each tile. It even helps with many early fragment tests and general
	//frame buffer throughput if you regularly clear the buffers. On Tile-Based Deferred 
	//Rendering GPUs (e.g. PowerVR SGX - all iOS devices) it is equally important for similar reasons. 


	 
} // attachTexturesToRainbowFBO


void ofApp::drawTextureFromFBOFile() {

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	Uniforms  extraUniforms;
	//GLint window_fbo =0; 
    //glGetIntegerv(GL_FRAMEBUFFER_BINDING, &window_fbo); 
	//glBindFramebuffer(GL_FRAMEBUFFER, window_fbo);

	extraUniforms.put("uTexColor", shared_ptr<ShaderImageTexture2D_RGB_RGB>(new ShaderImageTexture2D_RGB_RGB("colorFBO.ppm", false)));

    g_NDCMatObj->draw( extraUniforms );


}
void ofApp::drawTextureFromSysBufferFile() {

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	checkGlErrors();

	Uniforms  extraUniforms;
	//GLint window_fbo =0; 
    //glGetIntegerv(GL_FRAMEBUFFER_BINDING, &window_fbo); 
	//glBindFramebuffer(GL_FRAMEBUFFER, window_fbo);

	extraUniforms.put("uTexColor", shared_ptr<ShaderImageTexture2D_RGB_RGB>(new ShaderImageTexture2D_RGB_RGB("systembuffer.ppm", true)));
    g_NDCMatObj->draw( extraUniforms );


}



void ofApp::readRainbowFBOPixels(char * fileName) {

	std::ofstream sceneFBOFile (fileName);
	
	unsigned char *pixels = new unsigned char [ g_windowWidth * g_windowHeight * 4 ];
	
	float *depths = new float [g_windowHeight * g_windowWidth];
	float *floatPixels = new float [g_windowHeight * g_windowWidth * 4];
	float *floatPixels2 = new float [g_windowHeight * g_windowWidth * 4];
	//float *depth = new float [1];
	//float *floatPixel = new float [4];
 
 
 
	//glReadPixels(0,0, g_windowWidth, g_windowHeight, GL_DEPTH_COMPONENT, GL_FLOAT, depths ); 
	// The near plane is drawn, so that the depth is all 1. It is verified.
	
	// read the color buffer 

	glReadBuffer(GL_COLOR_ATTACHMENT0); // of the current framebuffer
	
	glReadPixels(0,0, g_windowWidth, g_windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, pixels); // fragColor

	glReadBuffer(GL_COLOR_ATTACHMENT1); // spect1
	glReadPixels(0,0, g_windowWidth, g_windowHeight, GL_RGBA, GL_FLOAT, floatPixels); //

	glReadBuffer(GL_COLOR_ATTACHMENT2); // spect2
	glReadPixels(0,0, g_windowWidth, g_windowHeight, GL_RGBA, GL_FLOAT, floatPixels2); //
	//glReadPixels simply returns bytes in the order R, G, B, R, G, B, ... 
	//(based on your setting of GL_RGB) from the bottom left of the screen going up to the top right. From the OpenGL documentation:
	//From there you can either reference a pixel's component location 
	//with data[(py * width + px) * 3 + component] where px and py are the pixel locations you want to look up, 
	//and component being the R, G, or B components of the pixel.

	/*
	GLint colorReadType;

	glGetIntegerv(GL_IMPLEMENTATION_COLOR_READ_TYPE, &colorReadType);
	GLint colorReadFormat;

	glGetIntegerv(GL_IMPLEMENTATION_COLOR_READ_TYPE, &colorReadFormat);

	cout << "colorReadFormat=" << colorReadFormat << " colorReadType=" << colorReadType << endl;
	*/
	sceneFBOFile << "The RGBA value of each pixel (i,j) may be background colors or value of some variables in the case of true fragments." << endl;
   
	    
	for (int j = 0; j < g_windowHeight; j++ )
	// print jth row
		for ( int i = 0; i < g_windowWidth ; i++) { // print ith column
	
			//ostream& operator<< (unsigned int val);
			//glReadPixels(i,j, 1,1, GL_DEPTH_COMPONENT, GL_FLOAT, depth ); 
			//float zWin0 =	 depths[i + j* g_windowWidth ]; 
			//float zWin1 =	 depth[0]; 
			sceneFBOFile << endl;

			sceneFBOFile  << "At " << "(" <<  i   << "," <<  j  << "):" ;  
	    
			//float zNDC = zWin0 * 2.0 - 1.0; // ( zWin = 1/2 zNDC + 1/2 by Viewport transformation: zWin in [0,1] )
	   	       
			//float zEye = -g_projectionMatrix(2,3) / ( zNDC + g_projectionMatrix(2,2) ); 
		
	  
			//sceneFBOFile  << "-zEye= " << -zEye <<  endl;  // -200 = FAR PLANE
	
			float red = floatPixels[ (j* g_windowWidth + i) *4 + 0 ]; 
			float green = floatPixels[ (j* g_windowWidth + i) *4 + 1 ]; 
			float blue = floatPixels[ (j* g_windowWidth + i) *4 + 2 ]; 
			float alpha= floatPixels[ (j* g_windowWidth + i) *4 + 3 ]; 
	  
			float zNDC = red * 2.0 - 1.0; // ( zWin = 1/2 zNDC + 1/2 by Viewport transformation: zWin in [0,1] )
	   	       
			float zEye = -g_projectionMatrix(2,3) / ( zNDC + g_projectionMatrix(2,2) ); 
			
			sceneFBOFile << "(height) = " << red<< endl;


			// the names of variables should be changed in each case
			//sceneFBOFile << "-Eye= " << zEye <<  endl;  // -200 = FAR PLANE, -5 = near plane
			//sceneFBOFile << "(zWin, zNDC, zEye, A) = " << Cvec4f( red, green, blue, alpha) << endl;
	
			red = floatPixels2[ (j* g_windowWidth + i) *4 + 0 ]; 
			green = floatPixels2[ (j* g_windowWidth + i) *4 + 1 ]; 
			blue = floatPixels2[ (j* g_windowWidth + i) *4 + 2 ]; 
			alpha= floatPixels2[ (j* g_windowWidth + i) *4 + 3 ]; 
	  
			//sceneFBOFile << "(zEye, tmin, tmax, A) = " << Cvec4f( red, green, blue, alpha) << endl;


		}

 
	sceneFBOFile.close();
	cout <<"background Scene output dumped" << endl;

	// inbind the currently bound framebuffer

	//glBindFramebuffer(GL_FRAMEBUFFER, g_savedFramebuffer);
	//glBindFramebuffer(GL_FRAMEBUFFER, 0);
} // readRainbowFBOPixels()


void ofApp::readSceneFBOPixels(char * fileName) {

	std::ofstream sceneFBOFile (fileName);
	
	unsigned char *pixels = new unsigned char [ g_windowWidth * g_windowHeight * 4 ];
 
	float *depths = new float [  g_windowHeight * g_windowWidth];
	float *floatPixels = new float [ g_windowHeight *g_windowWidth *  4];
	float *floatPixels2 = new float [ g_windowHeight *g_windowWidth *  4];
	//float *depth = new float [1];
	//float *floatPixel = new float [ 4];
 
 
 

	// The near plane is drawn, so that the depth is all 1. It is verified.
 
	// read the color buffer 

	glReadBuffer(GL_COLOR_ATTACHMENT0); // of the current framebuffer
  
	glReadPixels(0,0, g_windowWidth, g_windowHeight, GL_RGBA, GL_UNSIGNED_BYTE, pixels ); // fragColor

	glReadPixels(0,0, g_windowWidth, g_windowHeight, GL_DEPTH_COMPONENT, GL_FLOAT, depths ); 

	//glReadBuffer(GL_COLOR_ATTACHMENT1); // spect1
	//glReadPixels(0,0, g_windowWidth, g_windowHeight, GL_RGBA, GL_FLOAT, floatPixels ); //

	//glReadBuffer(GL_COLOR_ATTACHMENT2); // spect2
	//glReadPixels(0,0, g_windowWidth, g_windowHeight, GL_RGBA, GL_FLOAT, floatPixels2 ); //
	// glReadPixels simply returns bytes in the order R, G, B, R, G, B, ... 
	// (based on your setting of GL_RGB) from the bottom left of the screen going up to the top right. From the OpenGL documentation:
	// From there you can either reference a pixel's component location 
	// with data[(py * width + px) * 3 + component] where px and py are the pixel locations you want to look up, 
	// and component being the R, G, or B components of the pixel.

	/*
	GLint colorReadType;

	glGetIntegerv(GL_IMPLEMENTATION_COLOR_READ_TYPE, &colorReadType);
	GLint colorReadFormat;

	glGetIntegerv(GL_IMPLEMENTATION_COLOR_READ_TYPE, &colorReadFormat);

	cout << "colorReadFormat=" << colorReadFormat << "  colorReadType =" << colorReadType << endl;
	*/
	sceneFBOFile << " The RGBA value of each pixel (i,j) may be background colors or value of some variables in the case of true fragments." << endl;
   
	    
	for (int j = 0; j < g_windowHeight; j++ )
			   // print jth row
		for ( int i = 0; i < g_windowWidth ; i++) { // print ith column
	
			// ostream& operator<< (unsigned int val);
			//glReadPixels(i,j, 1,1, GL_DEPTH_COMPONENT, GL_FLOAT, depth ); 
			float zWin =	 depths[i + j* g_windowWidth ]; 
			//float zWin1 =	 depth[0]; 
			sceneFBOFile << endl;

			sceneFBOFile  << "At " << "(" <<  i   << "," <<  j  << "):" ;  
	    
			//float zNDC = zWin0 * 2.0 - 1.0; // ( zWin = 1/2 zNDC + 1/2 by Viewport transformation: zWin in [0,1] )
	   	       
			//float zEye = -g_projectionMatrix(2,3) / ( zNDC + g_projectionMatrix(2,2) ); 
		
	  
			//sceneFBOFile  << "-zEye= " << -zEye <<  endl;  // -200 = FAR PLANE
			/*
			float red =	 pixels[ (j* g_windowWidth + i) *4 + 0 ]; 
			float green = pixels[ (j* g_windowWidth + i) *4 + 1 ]; 
			float blue = pixels[ (j* g_windowWidth + i) *4 + 2 ]; 
			float alpha= pixels[ (j* g_windowWidth + i) *4 + 3 ]; 
			*/

			float zNDC = zWin * 2.0 - 1.0; // ( zWin = 1/2 zNDC + 1/2 by Viewport transformation: zWin in [0,1] )
	   	       
			float zEye = -g_projectionMatrix(2,3) / ( zNDC + g_projectionMatrix(2,2) ); 
		
	  
			sceneFBOFile << "( g_projectionMatrix(2,2), g_projectionMatrix(2,3), zEye)= " << Cvec3( g_projectionMatrix(2,2), g_projectionMatrix(2,3), zEye) << endl;  // -200 = FAR PLANE, -5 = near plane

			/*
			sceneFBOFile << "(zWin, zNDC, zEye, A) = " << Cvec4f(red, green, blue, alpha) << endl;
	
			red = floatPixels2[ (j* g_windowWidth + i) *4 + 0 ]; 
			green = floatPixels2[ (j* g_windowWidth + i) *4 + 1 ]; 
			blue = floatPixels2[ (j* g_windowWidth + i) *4 + 2 ]; 
			alpha= floatPixels2[ (j* g_windowWidth + i) *4 + 3 ]; 
	  
			sceneFBOFile << "(zEye, tmin, tmax, A) = " << Cvec4f( red, green, blue, alpha ) << endl;
			*/

		}

 
    sceneFBOFile.close();
	cout <<"background Scene output dumped" << endl;

	// inbind the currently bound framebuffer

	//glBindFramebuffer(GL_FRAMEBUFFER, g_savedFramebuffer);
	//glBindFramebuffer(GL_FRAMEBUFFER, 0);
} // readSceneFBOPixels()



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





//// set a texture reference
//void setUniformTexture(const string & name, ofBaseHasTexture& img, int textureLocation);
//void setUniformTexture(const string & name, ofTexture& img, int textureLocation);
//void setUniformTexture(const string & name, int textureTarget, GLint textureID, int textureLocation);
	
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
		messageFile << "normal shader:" << error.what() << endl;

	}

	
	
	MaterialShader solidMat("basic+solid", "shaders/basic-gl3.vshader", "shaders/solid-gl3.fshader");

	// copy solid prototype, and set to wireframed rendering
	g_arcballMat.reset(new MaterialShader(solidMat));
 
	//g_arcballMat->getUniforms().put("uMaterialColor", Cvec4f(0.27f, 0.82f, 0.35f, 1.0));
	g_arcballMat->getUniforms().put("uMaterialColor", Cvec4f(1.0f, 0.0f, 0.0f, 0.1));
	
	g_arcballMat->renderStates_ = g_arcballMat->renderStates_.polygonMode(GL_FRONT_AND_BACK, GL_LINE);

	//g_arcballMat->getRenderStates().polygonMode(GL_FRONT_AND_BACK, GL_LINE);
	
	//g_arcballMat->getRenderStates().polygonMode(GL_FRONT_AND_BACK, GL_POINT);
	// copy solid prototype, and set to color white
	g_lightMat.reset(new MaterialShader(solidMat));
	g_lightMat->getUniforms().put("uMaterialColor", Cvec4f(1, 1, 1,1));
	// pick shader
	g_pickingMat.reset(new MaterialShader("basic+pick", "shaders/basic-gl3.vshader", "shaders/pick-gl3.fshader") );

};



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


void ofApp::initObjects() {
	int ibLen, vbLen;


	// 1: create a ground
	//static const float groundY = -2.0;      // y coordinate of the ground
    static const float groundSize = 200.0;   // half the ground length

	// temporary storage for a plane geometry
	getPlaneVbIbLen(vbLen, ibLen);	// get the sizes of the vertex buffer and the 
									// index buffer for a plane. The vertex buffer
									// size is the number of vertices.
	vector<VertexPNTBX> vtxGround (vbLen);
	vector<unsigned short> idxGround (ibLen);

	makePlane( groundSize * 2, vtxGround.begin(), idxGround.begin() );

	// vtxGround is an array of generic vertices (position, normal, tex coord, tangent, binormal)

  
	// A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
	/* VertexPN vtxGround[4] = {VertexPN( -groundSize, groundY, groundSize, 0, 1, 0), where groundY =0
								VertexPN( groundSize, groundY,  groundSize, 0, 1, 0),
								VertexPN(  groundSize, groundY, -groundSize, 0, 1, 0),
								VertexPN(  -groundSize, groundY, -groundSize, 0, 1, 0),
								};
  
	unsigned short idxGround[] = {0, 1, 2, 0, 2, 3};
	*/

  

	//shared_ptr<Geometry> groundGeometry ( new Geometry( &vtxGround[0], &idxGround[0], 4, 6,  &sqTex[0], numOfTexCoords ) );
	shared_ptr<Geometry> groundGeometry ( new SimpleIndexedGeometryPNTBX("ground", &vtxGround[0], &idxGround[0], vbLen, ibLen, GL_TRIANGLES ) );

 
	// SimpleIndexedGeometryPNTBX => init: vbo(new FormattedVbo( Vertex::FORMAT ) ), ibo(new FormattedIbo(size2IboFmt(sizeof(Index)))) 
	// where Vertex is either VertexPN, ...:

	//const VertexFormat VertexPN::FORMAT = VertexFormat(sizeof(VertexPN))
	//										.put("aPosition", 3, GL_FLOAT, GL_FALSE, offsetof(VertexPN, p))
	//										.put("aNormal", 3, GL_FLOAT, GL_FALSE, offsetof(VertexPN, n));

 

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

	//makeCube(30,6.4, 25, vtxCube1.begin(), idxCube1.begin() );	// fill vtxCube and idxCube, starting
	//															// from their starting pointers => a list of vertex coords,
	//															// tex coords, normal coords, tangent coords, binormal coords

	// for debugging
	makeCube(30, 26.4, 25, vtxCube1.begin(), idxCube1.begin() );	// fill vtxCube and idxCube, starting
																	// from their starting pointers => a list of vertex coords,
																	// tex coords, normal coords, tangent coords, binormal coords
	// create vbo and ibo for vtxCube and idxCube by using SimpleIndexedGeometryPNTBX, where
	// SimpleIndexedGeometry<VertexPNTBX, unsigned short> SimpleIndexedGeometryPNTBX;

	shared_ptr<Geometry> cubeGeometry1 ( new SimpleIndexedGeometryPNTBX("cube1", &vtxCube1[0], &idxCube1[0], vbLen, ibLen, GL_TRIANGLES ) ); // reset(T * p)

	//shared_ptr<Geometry> cubeGeometry ( new Geometry( &vtxCube[0], &idxCube[0], vbLen, ibLen, &sqTex[0], numOfTexCoords ) ); // reset(T * p)

 
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

	//makeCube(20, 12.7, 25, vtxCube2.begin(), idxCube2.begin() );	// fill vtxCube and idxCube, starting
	//																// from their starting pointers => a list of vertex coords,
	//																// tex coords, normal coords, tangent coords, binormal coords

	// for debugging
	makeCube(20, 62.7, 25, vtxCube2.begin(), idxCube2.begin() ); // fill vtxCube and idxCube, starting
	// create vbo and ibo for vtxCube and idxCube by using SimpleIndexedGeometryPNTBX, where
	// SimpleIndexedGeometry<VertexPNTBX, unsigned short> SimpleIndexedGeometryPNTBX;

	shared_ptr<Geometry> cubeGeometry2 ( new SimpleIndexedGeometryPNTBX("cube2", &vtxCube2[0], &idxCube2[0], vbLen, ibLen, GL_TRIANGLES ) ); // reset(T * p)

	//shared_ptr<Geometry> cubeGeometry ( new Geometry( &vtxCube[0], &idxCube[0], vbLen, ibLen, &sqTex[0], numOfTexCoords  ) ); // reset(T * p)


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

	shared_ptr<Geometry> sphereGeometry ( new SimpleIndexedGeometryPNTBX( "sphere", &vtxSphere[0], &idxSphere[0], vbLen, ibLen, GL_TRIANGLES ) ); // reset(T * p)
 
	// or  shared_ptr<Geometry> sphereGeometry =  new Geometry( &vtx[0], &idx[0], vbLen, ibLen ); // reset(T * p)
	// NOTE: here vtx and idx arrays are copied to the vertex buffer object internally, so this data need not be global.
	


	// create the object frames for the created geometry

	// The following does not work, because the  matrix created by SgRbtNode::makeTranslation() is
	// destroyed outside of this function initObjects(). To avoid it, create the matrix by the "new" method.
	// The objects created by new methods remain unless removed explicitly. This is handled by shared_ptr<>. 
																				   
	//shared_ptr<SgRbtNode> groundRbt ( &SgRbtNode::makeTranslation( Cvec3(0,0,0) ) ); 
	//shared_ptr<SgRbtNode> cubeRbt ( &SgRbtNode::makeTranslation( Cvec3(5.0,0,0) ) );
	//shared_ptr< SgRbtNode> sphereRbt ( &SgRbtNode::makeTranslation( Cvec3(0,0,0) ) );
  
   
	
	shared_ptr<SgRbtNode> groundRbt ( new SgRbtNode( g_SceneRbt * SgRbtNode::makeTranslation( Cvec3(0,0,0) ) ) ); // this uses the non-default copy constructor ?
	//shared_ptr<SgRbtNode> cubeRbt1 ( new SgRbtNode( g_SceneRbt * SgRbtNode::makeTranslation( Cvec3(0, 3.2, -12.5) ) )  );
	shared_ptr<SgRbtNode> cubeRbt1 ( new SgRbtNode( g_SceneRbt * SgRbtNode::makeTranslation( Cvec3(0, 26.4/2.0, -25.0/2.0) ) ) );
	shared_ptr< SgRbtNode> sphereRbt ( new SgRbtNode ( SgRbtNode::makeTranslation( Cvec3(-0.5,0,0) ) ) );
	//shared_ptr<SgRbtNode> cubeRbt2 ( new SgRbtNode( g_SceneRbt * SgRbtNode::makeTranslation( Cvec3(25, 6.35, -12.5) ) ) );
  
	shared_ptr<SgRbtNode> cubeRbt2 ( new SgRbtNode( g_SceneRbt * SgRbtNode::makeTranslation( Cvec3(25, 62.7/2.0, -25.0/2.0) ) ) );
  
	   
	shared_ptr<Object> ground ( new Object( "ground",  groundRbt, groundGeometry, g_bumpFloorMat) );
	shared_ptr<Object> cube1 ( new Object( "cube1", cubeRbt1, cubeGeometry1, g_blueDiffuseMat) );
	shared_ptr<Object> cube2 ( new Object( "cube2", cubeRbt2, cubeGeometry2, g_redDiffuseMat) );
	shared_ptr<Object> sphere ( new Object( "sphere", sphereRbt, sphereGeometry, g_blueDiffuseMat) );

	cout << endl << "cube1 RBT=" << cubeRbt1 << endl;


	// 1: create the object id color for picking purpose 

	Cvec4 pickColor;
	unsigned int id; // 32 bit unsigned


	id = g_objId ++; 
	//id = 0; 

	id = id << 8;  // RGBA: fill A field with 0's, shifting the value of id left 8 bits
	id = id | 255 ; // id for  ground; make A field to be all 1's



	pickColor = g_picker.idToColor(id); // id => linearTrans(idColor)

	ground->pickColor = pickColor;	// the pickColor = linearTrans(idColor) will be sent to the shader Uniform variable when drawing 
									// That is, by fragColor = linearTrans(id). WHen written to the framebuffer, sRGBTrans(linearTrans(idColor)) =
									// idColor will written. When read by glReadPixel, this value is 
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

   
  
  
	// 3: another cube
	id = g_objId ++;
	id = id << 8;
	id = id | 255; // id for sphere 

	pickColor = g_picker.idToColor(id);

	//sphere->pickColor = pickColor;
	cube2->pickColor = pickColor;

	cout << "cube2 object id=" << id <<"," << "object pseudo color=" << pickColor << endl;

	// for debugging

	g_objectPtrList.push_back( cube2 );

	g_picker.addToMap( id, cube2 );

	// 4: sphere 
	id = g_objId ++;
	id = id << 8;
	id = id | 255; // id for  sphere 

	pickColor = g_picker.idToColor(id);

	//sphere->pickColor  = pickColor;
	sphere->pickColor  = pickColor;

	cout << "sphere object id=" << id <<"," << "object pseudo color=" << pickColor << endl;

	// for debugging

	g_objectPtrList.push_back( sphere );

	g_picker.addToMap( id, sphere );

  
 

} //initObjects()


 void ofApp::initAABB() {
	int ibLen, vbLen;


	getCubeVbIbLen(vbLen, ibLen);

 
	vector<VertexPNTBX> vtxAABB(vbLen); // vtxCube:  a vector of VertexPNTBX type with length vbLen
										// what happens if I use VertexPNX ?

	vector<unsigned short> idxAABB(ibLen); // idxCube: 

	makeCube( g_AABBsize[0], g_AABBsize[1], g_AABBsize[2], vtxAABB.begin(), idxAABB.begin() ); 

 
	shared_ptr<Geometry> AABBGeometry ( new SimpleIndexedGeometryPNTBX("AABB", &vtxAABB[0], &idxAABB[0], vbLen, ibLen, GL_TRIANGLES ) );

 
	// for National Modern Musium
	//g_AABBRbt.reset( new SgRbtNode ( g_SceneRbt * Matrix4::makeTranslation( Cvec3( 0, g_AABBsize[1]/2 + 6.4, -g_AABBsize[2]/2.0 ) ) ) );
	g_AABBRbt.reset( new SgRbtNode ( g_SceneRbt * Matrix4::makeTranslation( Cvec3( 0, g_AABBsize[1]/2, -g_AABBsize[2]/2.0 ) ) ) );
	shared_ptr<Object> AABB ( new Object( "AABB", g_AABBRbt, AABBGeometry, g_arcballMat) );

  
	cout << "At initAABB: current AABB frame (center):" << "\n" << *g_AABBRbt << endl;

	Cvec4 pickColor;
	unsigned int id; // 32 bit unsigned

  
	id= 0; // set the  AABB box id to zero. // AABB box can be picked and moved.
	id = id << 8;  
	id = id | 255; // 
	pickColor  = g_picker.idToColor(id ); // pickColor is an sRGB non-linear color

	cout << "At initAABB: AABB object id=" << id <<"," << "object pseudo color=" << pickColor << endl;

	AABB->pickColor = pickColor;


	g_objectPtrList.push_back( AABB ); //for debugging


	g_picker.addToMap( id, AABB );

 
 } // initAABB

 // RAINBOW INIT
void ofApp::initRainbow() {
	int ibLen, vbLen;

	// create the near plane quad, which is the front of the view frustrum.
	//g_frustMinFov = 60.0;  
	//g_frustFovY = g_frustMinFov; 
	//g_frustNear = -1.1;     // near plane
	//g_frustFar = -100.0;    // far plane
	//y = tan(60/2 * PI / 180) * z => y = tan(60/2) * (-g_frustNear) is the height of the 
	// near quad of the near plane. 


	float quadZLocation = (float) g_frustNear + (-0.1); // a little bit further than the near plane; the near plane itself will not be seen

	float y = std::tan( g_frustFovY * 0.5 * PI / 180 ) * (-quadZLocation);


	static const float quadHeight = 2.0 * y;   // 

	// temporary storage for a plane geometry
	getPlaneVbIbLen( vbLen, ibLen); // get the sizes of the vertex buffer and the 
									// index buffer for a plane. The vertex buffer
									// size is the number of vertices.
	vector<VertexPNTBX> vtxQuad (vbLen);
	vector<unsigned short> idxQuad (ibLen);

	// create the quad at g_frustNear
	float aspectRatio = g_windowWidth / static_cast <double> (g_windowHeight);
	float quadWidth = quadHeight * aspectRatio;

	

	makeNearPlaneQuad( quadWidth, quadHeight, quadZLocation, vtxQuad.begin(), idxQuad.begin() );

	
	//shared_ptr<Geometry> groundGeometry ( new Geometry( &vtxGround[0], &idxGround[0], 4, 6,  &sqTex[0], numOfTexCoords ) );
	shared_ptr<Geometry> nearPlaneQuadGeometry ( new SimpleIndexedGeometryPNTBX("rainbow", &vtxQuad[0], &idxQuad[0], vbLen, ibLen, GL_TRIANGLES ) );

	// rainbow shader

	//MaterialShader rainbowAndSceneMat ("rainbow+rainbow", "shaders/rainbow-gl3.vshader", "shaders/rainbow-gl3.fshader");
   
	MaterialShader rainbowAndSceneMat ("rainbow+rainbow", "shaders/rainbow-gl3.vshader", "shaders/rainbow-RGB-STAN-gl3.frag");
	//MaterialShader rainbowOnlyMat ("rainbow+rainbow", "shaders/rainbow-gl3.vshader", "shaders/rainbow-old-org-gl3.fshader");
	//MaterialShader rainbowMat ("basic+rainbow", "shaders/FullScreenQuad-gl3.vshader", "shaders/FBOtexture-gl3.fshader");
	//MaterialShader rainbowOnlyMat ("rainbow+rainbow", "shaders/rainbow-gl3.vshader", "shaders/rainbow-8-22-gl3.fshader");
 
	//MaterialShader rainbowOnlyMat ("rainbow+rainbow", "shaders/rainbow-gl3.vshader", "shaders/rainbow-STAN-gl3.fshader");

	//MaterialShader rainbowOnlyMat ("rainbow+rainbow", "shaders/rainbow-gl3.vshader", "shaders/rainbow-old-gl3.fshader");


	// MaterialShader initializes  its member programDesc_ which has a member program, 
	// by programDesc_ ( GlProgramLibrary::getSingleton().getProgramDesc(vsFilename, fsFilename) ), which in turn
	// sets the output variable for shaders, e.g.:
	//if (!g_Gl2Compatible) {
	//	index 0, 1, 2 refers to the index of the array attachments in:
	//	GLuint attachments[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1,GL_COLOR_ATTACHMENT2 }; 
	//	glDrawBuffers(3,  attachments);

	//	glBindFragDataLocation(program, 0, "fragColor");
	//	glBindFragDataLocation(program, 1, "phaseSpectrum");
	//	glBindFragDataLocation(program, 2, "rainbowColorSpectrum");
	//}

	int primaryWidth = PrimarySpectrum::getPrimarySpectrumsWidth();
	int XYZWidth = PrimarySpectrum::getXYZSpectrumsWidth();
	float *primarySpectrums = new float[ primaryWidth ];
	float *XYZSpectrums = new float[ XYZWidth ];

	PrimarySpectrum::initPrimarySpectrums(XYZSpectrums, primarySpectrums);

	// rainbow mat
	g_rainbowAndSceneMat.reset(new MaterialShader( rainbowAndSceneMat) );
	g_rainbowAndSceneMat->getUniforms().put("uXYZSpectrums", 
			shared_ptr<ShaderImageTexture1D_RF_RF>(new ShaderImageTexture1D_RF_RF(XYZWidth, XYZSpectrums)));

	g_rainbowAndSceneMat->getUniforms().put("uPrimarySpectrums", 
			shared_ptr<ShaderImageTexture1D_RF_RF>(new ShaderImageTexture1D_RF_RF(primaryWidth, primarySpectrums)));

	g_rainbowAndSceneMat->getUniforms().put("uScatTex", shared_ptr<ShaderImageTexture2D_RF_RF>(new ShaderImageTexture2D_RF_RF("rainbow/scattering2_120_fort.txt")));


	g_rainbowAndSceneMat->getUniforms().put("uPhaseTex", shared_ptr<ShaderImageTexture3D_RF_RF>(new ShaderImageTexture3D_RF_RF("rainbow/phaseFunction2_120_100_fort.txt")));

	// quadRbt for the near plane quad is the same as the eye frame, because the full screen quad (near plane) is specified 
	// relative to the eye frame. 

	// When the camera is changed, g_rainbowAndSceneMatObj->objectRbt should be changed as well.
	// To avoid it, make MVM = g_invEyeRbt * g_eyeRbt = Identity.
	// Or Use g_NDCMatObj to use the full screen quad.

	std::shared_ptr<SgRbtNode> nearPlaneQuadRbt ( new SgRbtNode ( g_eyeRbt) );

	// rainbow material obj
  
	g_rainbowAndSceneMatObj.reset( new Object( "rainbowAndScene", nearPlaneQuadRbt, nearPlaneQuadGeometry, g_rainbowAndSceneMat ) );



	// rainbow only mat

	MaterialShader rainbowOnlyMat ("rainbow+rainbow", "shaders/rainbow-gl3.vshader", "shaders/rainbow-RGB-only-STAN-gl3.frag");
	g_rainbowOnlyMat.reset(new MaterialShader( rainbowOnlyMat) );
   
	g_rainbowOnlyMat->getUniforms().put("uXYZSpectrums", 
			shared_ptr<ShaderImageTexture1D_RF_RF>(new ShaderImageTexture1D_RF_RF(XYZWidth, XYZSpectrums)));

	g_rainbowOnlyMat->getUniforms().put("uPrimarySpectrums", 
			shared_ptr<ShaderImageTexture1D_RF_RF>(new ShaderImageTexture1D_RF_RF(primaryWidth, primarySpectrums)));

	g_rainbowOnlyMat->getUniforms().put("uScatTex", shared_ptr<ShaderImageTexture2D_RF_RF>(new ShaderImageTexture2D_RF_RF("rainbow/scattering2_120_fort.txt")));


	g_rainbowOnlyMat->getUniforms().put("uPhaseTex", shared_ptr<ShaderImageTexture3D_RF_RF>(new ShaderImageTexture3D_RF_RF("rainbow/phaseFunction2_120_100_fort.txt")));

  
	g_rainbowOnlyMatObj.reset( new Object( "rainbowOnly", nearPlaneQuadRbt, nearPlaneQuadGeometry, g_rainbowOnlyMat ) );


 } // initRainbow()
 
void ofApp::initAABBRainbow() {
	int ibLen, vbLen;

	getCubeVbIbLen(vbLen, ibLen);

	// Temporary storage for cube geometry

	// T = tangent, X = coordinates

	vector<VertexPNTBX> vtxAABB(vbLen); // vtxCube:  a vector of VertexPNTBX type with length vbLen
										// what happens if I use VertexPNX ?

	vector<unsigned short> idxAABB(ibLen); // idxCube: 

	makeCube( g_AABBsize[0], g_AABBsize[1], g_AABBsize[2], vtxAABB.begin(), idxAABB.begin() ); 

	// create vbo and ibo for vtxCube and idxCube by using SimpleIndexedGeometryPNTBX, where
	// SimpleIndexedGeometry<VertexPNTBX, unsigned short> SimpleIndexedGeometryPNTBX;

	// Upload the geometry data and bind it to the Vbo 
	//shared_ptr<Geometry> AABBGeometry ( new SimpleIndexedGeometryPNTBX("AABB", &vtxAABB[0], &idxAABB[0], vbLen, ibLen, GL_QUADS ) ); 
	// GL_QAUDS has been removed from the core profile => ENUM error is caused
	shared_ptr<Geometry> AABBGeometry ( new SimpleIndexedGeometryPNTBX("AABBRainbow", &vtxAABB[0], &idxAABB[0], vbLen, ibLen, GL_TRIANGLES ) );

 
	// for National Modern Musium: g_AABBRbt is set in initAABB()

	//g_AABBRbt.reset( new SgRbtNode ( g_SceneRbt * 
	//							Matrix4::makeTranslation( Cvec3( 0, g_AABBsize[1]/2 + 6.4, -g_AABBsize[2]/2.0 ) ) ) );
	// for debugging

  
	//g_AABBRbt.reset( new SgRbtNode ( g_eyeRbt * 
	//							Matrix4::makeTranslation( Cvec3( 0, g_AABBsize[1]/2 + 6.4, -g_AABBsize[2]/2.0 ) ) ) );

	cout << " At initAABBRainbow: AABBRbt =" << "\n" << *g_AABBRbt << endl;

	// for Asian Game [the world frame is assumed right at the front of the AABB box

	//g_AABBRbt.reset( new SgRbtNode ( g_SceneRbt * 
	//	                       Matrix4::makeTranslation( Cvec3( 0, g_AABBsize[1]/2, -g_AABBsize[2]/2.0) ) ) );
	
	// rainbow shader

 
	//MaterialShader rainbowOnlyMat ("rainbow+rainbow", "shaders/rainbow-gl3.vshader", "shaders/rainbow-only-gl3.fshader");
	MaterialShader rainbowOnlyMat ("rainbow+rainbow", "shaders/rainbow-gl3.vshader", "shaders/rainbow-only-STAN-gl3.frag");
	g_AABBRainbowOnlyMat.reset( new MaterialShader( rainbowOnlyMat) );

	g_AABBRainbowOnlyMat->getUniforms().put("uScatTex", shared_ptr<ShaderImageTexture2D_RF_RF>(new ShaderImageTexture2D_RF_RF("rainbow/scattering2_120_fort.txt")));

	g_AABBRainbowOnlyMat->getUniforms().put("uPhaseTex", shared_ptr<ShaderImageTexture3D_RF_RF>(new ShaderImageTexture3D_RF_RF("rainbow/phaseFunction2_120_100_fort.txt")));

	g_AABBRainbowOnlyMatObj.reset ( new Object( "AABBRainbow", g_AABBRbt, AABBGeometry, g_AABBRainbowOnlyMat) );
  
 
}//initAABBRainbow()
 

void ofApp::initFullScreen() {
  int ibLen, vbLen;

	// create the near plane quad, which is the front of the view frustrum.
	//g_frustMinFov = 60.0;  
	//g_frustFovY = g_frustMinFov; 
	//g_frustNear = -1.1;    // near plane
	//g_frustFar = -100.0;    // far plane
	//y = tan( 60/2 * PI / 180) * z => y = tan(60/2) * (-g_frustNear) is the height of the 
	// near quad of the near plane. 

	float quadZLocation = (float) g_frustNear + (-0.1); // a little bit further than the near plane; the near plane itself will not be seen

	static const float y = std::tan( g_frustFovY * 0.5 * PI / 180.0 ) * (-quadZLocation);


	static const float quadHeight = 2.0 * y;

 

	// temporary storage for a plane geometry
	getPlaneVbIbLen( vbLen, ibLen); // get the sizes of the vertex buffer and the 
									// index buffer for a plane. The vertex buffer
									// size is the number of vertices.
	vector<VertexPNTBX> vtxQuad (vbLen);
	vector<unsigned short> idxQuad (ibLen);

	// create the quad at g_frustNear

	float aspectRatio = g_windowWidth / static_cast <double> (g_windowHeight);
	float quadWidth = aspectRatio * quadHeight;

	//float quadZLocation = (float) g_frustNear + (-0.1); // a little bit further than the near plane; the near plane itself will not be seen
  
 
	makeNearPlaneQuad( quadWidth, quadHeight, quadZLocation, vtxQuad.begin(), idxQuad.begin() );

	cout << "quad Width=" << quadWidth << "quad Height=" << quadHeight << endl;

	//shared_ptr<Geometry> groundGeometry ( new Geometry( &vtxGround[0], &idxGround[0], 4, 6,  &sqTex[0], numOfTexCoords ) );
	shared_ptr<Geometry> nearPlaneQuadGeometry  ( new SimpleIndexedGeometryPNTBX( "fullscreen", &vtxQuad[0], &idxQuad[0], vbLen, ibLen, GL_TRIANGLES ) );
   
	MaterialShader fullScreenMat ("fullScreen", "shaders/FullScreenQuad-gl3.vshader", "shaders/FBOtexture-gl3.fshader");

	int primaryWidth = PrimarySpectrum::getPrimarySpectrumsWidth();
	int XYZWidth = PrimarySpectrum::getXYZSpectrumsWidth();
	float *primarySpectrums = new float[ primaryWidth ];
	float *XYZSpectrums = new float[ XYZWidth ];

	PrimarySpectrum::initPrimarySpectrums(XYZSpectrums, primarySpectrums);

	g_fullScreenMat.reset(new MaterialShader( fullScreenMat));
   

	g_fullScreenMat->getUniforms().put("uXYZSpectrums", 
			shared_ptr<ShaderImageTexture1D_RF_RF>(new ShaderImageTexture1D_RF_RF(XYZWidth, XYZSpectrums)));

	g_fullScreenMat->getUniforms().put("uPrimarySpectrums", 
			shared_ptr<ShaderImageTexture1D_RF_RF>(new ShaderImageTexture1D_RF_RF(primaryWidth, primarySpectrums)));

  
	std::shared_ptr<SgRbtNode> nearPlaneQuadRbt ( new SgRbtNode ( g_eyeRbt) );
    
	g_fullScreenMatObj.reset( new Object( "fullscreen", nearPlaneQuadRbt, 
											nearPlaneQuadGeometry, g_fullScreenMat) );
	// print the vertices of the qaud

	for (int i =0; i < vtxQuad.size(); i++ ) {
		cout << " local quad[i]= " << vtxQuad[i].p << endl; 

	}
 
	for (int i =0; i < vtxQuad.size(); i++ ) {
		cout << " global quad[i]= " << ( *nearPlaneQuadRbt) *
										Cvec4f(vtxQuad[i].p,  1.0) << endl; 

	}


	const Matrix4 projMatrix = makeProjectionMatrix();
   
	for (int i =0; i < vtxQuad.size(); i++ ) {
		cout << " projected quad[i]= " <<  projMatrix * Cvec4f(vtxQuad[i].p, 1) << endl; 

	}


 } // initFullScreenQuad()

void ofApp::initNDC() {
	int ibLen, vbLen;

	 
	getPlaneVbIbLen(vbLen, ibLen);
  
	static const float quadSize = 2.0; 

	// temporary storage for a plane geometry
	getPlaneVbIbLen( vbLen, ibLen); // get the sizes of the vertex buffer and the 
									// index buffer for a plane. The vertex buffer
									// size is the number of vertices.
	vector<VertexPNTBX> vtxQuad (vbLen);
	vector<unsigned short> idxQuad (ibLen);

	makeNDC( quadSize, vtxQuad.begin(), idxQuad.begin() );

  
	//shared_ptr<Geometry> groundGeometry ( new Geometry( &vtxGround[0], &idxGround[0], 4, 6, &sqTex[0], numOfTexCoords )  );
	shared_ptr<Geometry> NDCGeometry ( new SimpleIndexedGeometryPNTBX("NDC", &vtxQuad[0], &idxQuad[0], vbLen, ibLen, GL_TRIANGLES)  );

		
	MaterialShader NDCMat ("NDC", "./shaders/NDC-gl3.vshader", "./shaders/FBOtexture-gl3.fshader");
     

	g_NDCMat.reset(new MaterialShader( NDCMat) );

	int primaryWidth = PrimarySpectrum::getPrimarySpectrumsWidth();
    int XYZWidth = PrimarySpectrum::getXYZSpectrumsWidth();
    float *primarySpectrums = new float[ primaryWidth ];
    float *XYZSpectrums = new float[ XYZWidth ];

    PrimarySpectrum::initPrimarySpectrums(XYZSpectrums, primarySpectrums);

    g_NDCMat->getUniforms().put("uXYZSpectrums", 
			shared_ptr<ShaderImageTexture1D_RF_RF>(new ShaderImageTexture1D_RF_RF(XYZWidth, XYZSpectrums)));

    g_NDCMat->getUniforms().put("uPrimarySpectrums", 
			shared_ptr<ShaderImageTexture1D_RF_RF>(new ShaderImageTexture1D_RF_RF(primaryWidth, primarySpectrums)));
	
    // NDCRbt is the identity matrix
    shared_ptr<SgRbtNode> NDCRbt ( new SgRbtNode ( Matrix4() ) );

	g_NDCMatObj.reset ( new Object("NDC",  NDCRbt, NDCGeometry, g_NDCMat) );

} // initNDC()
 

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){  // call back function for event processing
 
	// instance->windowW = w; => the new window size is already set before calling windowResized()

	// instance->windowH = h;
	g_windowHeight = h;
	g_windowWidth = w;
   
	glViewport(0, 0, g_windowWidth, g_windowHeight);
  
	cout << "WINDOW RESIZED: Size of window is now " << g_windowWidth << "x" << g_windowHeight << endl;
	cout << "results of ofGetHeight and ofGetWidth= " << ofGetWidth() << "x" << ofGetHeight() << endl;

	//cerr << "Size of window is now " << g_windowWidth << "x" << g_windowHeight << endl;
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


	/*cout << "dragged mouse button " << button << "X = " << x << "Y = " << y << endl;
	// Here the moved button is 0, meaning that which button is pressed while moving the mouse is ignored.
 
	cout << "Control Key Pressed =" << ofGetKeyPressed( OF_KEY_CONTROL ) << endl;
	cout << "Shift  Key Pressed =" << ofGetKeyPressed( OF_KEY_SHIFT ) << endl;
	cout << "ALT Key Pressed =" << ofGetKeyPressed( OF_KEY_ALT ) << endl;
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
				// cout << " left button moves" << endl;
				m = Matrix4::makeTranslation( Cvec3( dx, dy, 0.0) * 0.01 ); 
				//m = Matrix4::makeTranslation( Cvec3( dx, 0.0, 0.0) * 0.01 ); 
				break;
			case OF_MOUSE_BUTTON_RIGHT:
				// cout << " right button moves" << endl;
				m = Matrix4::makeTranslation( Cvec3( 0.0, 0.0, -dy) * 0.01 ); 
				break;

			default: 
				m = Matrix4::makeTranslation( Cvec3( 0.0, 0.0, 0.0) ); 
		} // switch

		g_eyeRbt *= m;

		//window = ofPtr<ofAppBaseWindow>(new ofAppGLFWWindow());
		//ofPtr<ofAppGLFWWindow> appwindow = window;
		//appwindow ->display();

		g_redrawWindowEvent = true;

	} // if (! g_rotation)

	else if ( ofGetKeyPressed(OF_KEY_CONTROL) && ofGetKeyPressed(OF_KEY_ALT) ) { //  rotate the camera

		switch ( button ) {
			case OF_MOUSE_BUTTON_LEFT:

				//inline static ofMatrix4x4 newRotationMatrix( float angle, const ofVec3f& axis);
				if ( abs(dy) > abs(dx) ) {
					cout << " rotate about the x axis" << endl;

					m = Matrix4::makeXRotation( dy * 0.01 ); // unit= degree; pitch motion: up and down
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
				m = Matrix4::makeTranslation( Cvec3(0.0, 0.0, 0.0) ); 
		} // switch

		g_eyeRbt *= m;

		//ofPtr<ofAppBaseWindow> window;  // changed by Moon Jung, 2014/4/23
	  
		//window = ofPtr<ofAppBaseWindow>(new ofAppGLFWWindow());
		//ofPtr<ofAppGLFWWindow> appwindow = window;
		//appwindow ->display();

		g_redrawWindowEvent = true;
	} // else if
  
  

	else if ( ofGetKeyPressed(OF_KEY_SHIFT) && !ofGetKeyPressed(OF_KEY_ALT) )   { // translate the selected object


		if ( g_currentPickedObject == nullptr ) { // no object is pickefd up, so no need to move
			return;
		}
	  

		switch ( button ) {
		case OF_MOUSE_BUTTON_LEFT:
			//cout  << " left button moves" << endl;
			m = Matrix4::makeTranslation( Cvec3(dx, dy, 0.0) * 0.01 ); 
			//m = Matrix4::makeTranslation( Cvec3( dx, 0.0, 0.0) * 0.01 ); 
			break;
		case OF_MOUSE_BUTTON_RIGHT:
			//cout  << " right button moves" << endl;
			m = Matrix4::makeTranslation( Cvec3(0.0, 0.0, -dy) * 0.01 ); 
			break;

		default: 
			m = Matrix4::makeTranslation( Cvec3(0.0, 0.0, 0.0) ); 
		} // switch

		assert( ("g_currentPickedRbtNode should not be Null", g_currentPickedObject != nullptr) );

		*(g_currentPickedObject->objectRbt) = *(g_currentPickedObject-> objectRbt) * m;

		if ( g_currentPickedObject->objectName =="AABB" ) {
		  
			g_AABBRbt = (g_currentPickedObject->objectRbt);
		  	 
		}

		g_redrawWindowEvent = true;

	} // if ( translate )

    else if ( ofGetKeyPressed(OF_KEY_SHIFT) && ofGetKeyPressed(OF_KEY_ALT) ) { //  rotate the selected object
	   
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
			m = Matrix4::makeTranslation( Cvec3(0.0, 0.0, 0.0) ); 
		} // switch

		assert( ("g_currentPickedRbtNode should not be Null", g_currentPickedObject != nullptr) );

		*(g_currentPickedObject->objectRbt) = *(g_currentPickedObject-> objectRbt) * m;

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

  
	if ( ofGetKeyPressed(OF_KEY_SHIFT) ) { // shift +  any mouse press leads to pickin an object
		cout  << "draw for picking " << endl;
		// for temp debugging
	
		drawForPicking();    
		// In the picking mode, you render the simplified version of the scene to the back buffer,
		// while maintaining the scene of the screen intact, because it is in the front buffer.
		                      
		// Read pixel operation needed for picking is performed with respect to the back buffer,
		// not to the front buffer
		
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
void ofApp::keyPressed (int key){ 
	// OF_KEY_LEFT_CONTROL, OF_KEY_LEFT_SHIFT,
	// OF_KEY_LEFT_ALT,
	// OF_MOUSE_BUTTON_LEFT
	/* 
	cout << "pressed key=" << key<< endl;
	cout << "OF_KEY_SHIFT= " << OF_KEY_SHIFT << endl;
	cout << "OF_KEY_CONTROL= " << OF_KEY_CONTROL << endl;
	cout << "OF_KEY_ALT= " << OF_KEY_ALT << endl;
	*/

	switch (key) {

	case '1': g_renderMode =1;
		g_redrawWindowEvent = true;
		break;
	case '2': g_renderMode = 2;
		g_redrawWindowEvent = true;
		break;
	case '3': g_renderMode = 3;
		g_redrawWindowEvent = true;
		break;
	case '4': g_renderMode = 4;
		g_redrawWindowEvent = true;
		break;
  
	case 'q':  renderToFBOAndDump( 1 );
		break;
	case 'w':  renderToFBOAndDump( 2 );
		break;
	case 'e':  renderToFBOAndDump( 3 );
		break;
 

	case OF_KEY_ESC: 
	exit();			// ESC 
	// case OF_KEY_CONTROL: 
	//g_camera_mode = true;
	//break;

	case 'h': 
		cout << " ============== H E L P ==============\n\n"
		<< "h             help menu\n"
		<< "s             save screenshot\n"
		<< "f             Toggle flat shading on/off.\n"
		<< "ESC           exit\n"
		<< "Cntrl         camera moving mode\n"
		<< "Alt           rotate the selected object or the camera\n"
		<< "Shift         picking mode\n" << endl;
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

	//cout << "released key=" << key << endl;
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
    if(name == "SHOW ACTIVE") {
        ofxUIToggle *toggle = (ofxUIToggle *) e.widget;
        ddl->setShowCurrentSelected(toggle->getValue()); 
    }
    else if(name == "DROPDOWN") {
        ofxUIDropDownList *ddlist = (ofxUIDropDownList *) e.widget;
        vector<ofxUIWidget * > &selected = ddlist->getSelected(); 

		g_renderMode = 0; // neither background scene or rainbow; can be changed
		                  // below.

		for(int i = 0; i < selected.size(); i++) {  
			  string action_name = selected[i]->getName(); 
			  cout << action_name << endl;
		  


			if ( action_name == "render Scene") {
				checkGlErrors();
				renderSceneToSysBuffer(); // to system framebuffer

				checkGlErrors();
				g_drawnToBackbuffer = true;
				g_renderMode = 1; // background scene only
		 
		    
			}
			else if ( action_name == "dump Scene Sysbuffer to Image File" ) { // system buffer dump
				//glFlush();
				glReadBuffer(GL_FRONT);
				writePpmScreenshot(g_windowWidth, g_windowHeight, "systembuffer.ppm");

				g_drawnToBackbuffer = true;
			}
	   
			else if ( action_name == "renderToFBO and Dump" ) { // FBO dump
				checkGlErrors();
				renderToFBOAndDump( g_renderMode );
				checkGlErrors();
	
			}
			else if ( action_name == "dump SceneFBO to Image File" ) { // FBO dump
				checkGlErrors();
				printFBO(g_fboScene, g_windowWidth, g_windowHeight, "colorFBO.ppm" );
				checkGlErrors();
	
			}
			else if ( action_name == "render Rainbow") { 
				checkGlErrors();
				renderRainbowAndSceneToScreen();
				checkGlErrors();
				g_renderMode = 3; // background + rainbow
	     
			}
	    
         
			else if ( action_name == "render Rainbow Only") { // render to FBO
				checkGlErrors();
				renderRainbowOnlyToScreen();
				checkGlErrors();

				g_renderMode = 2; // rainbow only
	     
			}
	    
			else if ( action_name == "render AABBRainbow Only") { // render to FBO
				checkGlErrors();
				renderAABBRainbowOnlyToScreen();
				checkGlErrors();

				g_renderMode = 4; // AABBrainbow only
	     
			}
	      
          
			else if (action_name == "draw Scene from FBO File") { // render to system buffer using dumped FBO texture
				drawTextureFromFBOFile();	// this will give a non gamma corrected image
											// because the FBO image is not in a sRGB format.
				g_drawnToBackbuffer = true;

			}
			else if (action_name  == "draw Scene from SysFramebuffer") { // render to system buffer using dumped FBO texture
				drawTextureFromSysBufferFile();
				g_drawnToBackbuffer = true;
			} 
	    
		
		} // for


	} // else if ("DROPDOWN")

	
} // ofApp::guiEvent

void ofApp::initGLState() {


	//glBindFramebuffer(GL_FRAMEBUFFER, 0); the framebuffer is already bound to the system-provided one

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
	
	glGetIntegerv( GL_FRAMEBUFFER_BINDING, &g_savedFramebuffer );

	cout << "the current bound framebuffer = " << g_savedFramebuffer << endl;

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


	// whether or not the color data will be clampled:
    //http://stackoverflow.com/questions/11211698/framebuffer-object-with-float-texture-clamps-values
	glClampColor(GL_CLAMP_READ_COLOR, GL_FALSE);
    //glClampColor(GL_CLAMP_VERTEX_COLOR, GL_FALSE);
    //glClampColor(GL_CLAMP_FRAGMENT_COLOR, GL_FALSE);

	//Note that this is just for the reading of the color. What gets written by the fragment shader will always be unclamped.

	generateFBO(); // generate an FBO and set it to g_fboScene 

	if (!g_Gl2Compatible)  glEnable(GL_FRAMEBUFFER_SRGB); 

	// For opengl3.x,  request an sRGB frame buffer using the call glEnable(GL FRAMEBUFFER SRGB).
    // Then we can pass linear [R, G, B] values out from the fragment shader. 
	// and they will be gamma corrected into the sRGB format before begin sent to the framebuffer.
	// Any writes to images that are not in the sRGB format should not be affected. 
	// So if you're writing to a floating-point image, nothing should happen.
	// Thus, you should be able to just turn it on and leave it that way; 
	// OpenGL will know when you're rendering to an sRGB framebuffer.							   
	

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
	// If enabled, the test compares the incoming alpha value with a reference value. 
	// The fragment is accepted or rejected depending on the result of the comparison. Both the reference value and the comparison function are set with glAlphaFunc(). 
	// By default, the reference value is zero, the comparison function is GL_ALWAYS[ always accept the fragment ], and the alpha test is disabled

	// The default framebuffer has buffers like GL_FRONT​, GL_BACK​, GL_AUXi​, GL_ACCUM​, and so forth. 
	// FBOs do not have these. Instead, FBOs have a different set of images.
	// Each FBO image represents an attachment point, a location in the FBO where an image can be attached.
	// FBOs have the following attachment points:

    // GL_COLOR_ATTACHMENTi, GL_DEPTH_ATTACHMENT, GL_STENCIL_ATTACHMENT,GL_DEPTH_STENCIL_ATTACHMENT​:
		
	// Each framebuffer object has a set of attachment points that logical buffers can be attached to. 
	
	// There are three different ways to switch between framebuffers. 

	// The first one is to use several different FBOs, one for each combination of logical buffers 
	// that you plan on using in you application. To change render targets you simply call glBindFramebuffer() 
	// with the FBO containing the setup you wish to render to.
	// Another way is to use a single FBO and alter the attachments
	// The third way is to use a single FBO, but instead of altering attachments you call glDrawBuffer() or glDrawBuffers() 
	// to change which color attachment(s) the rendering goes to.

	// for debugging: glDrawBuffer(GL_BACK); // GL_BACK is the draw buffer for the zero framebuffer; for double-buffering context
		
	
	

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
	

	GLint param[1];
	
	glGetFramebufferAttachmentParameteriv(GL_FRAMEBUFFER, GL_FRONT_LEFT,
	GL_FRAMEBUFFER_ATTACHMENT_COLOR_ENCODING, param);
	cout << "color encoding of sys buffer =" << param[0] << endl;

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
