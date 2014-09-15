#include "ofMain.h"
#include "ofApp.h"
#include "ofGLProgrammableRenderer.h"


int g_argc;
char **g_argv;

//========================================================================
int main(int argc, char **argv ){

	g_argc=argc;
	g_argv=argv;

	// to use a default GPU shader
	ofSetCurrentRenderer(ofGLProgrammableRenderer::TYPE);
	// => ofSetCurrentRenderer(ofPtr<ofBaseRenderer>(new ofGLProgrammableRenderer),setDefaults);

	// otherwise the fixed pipeline will be used

    //#define USE_PROGRAMMABLE_RENDERER will has the same effect

	
	ofSetLogLevel(OF_LOG_VERBOSE);


	ofSetupOpenGL(1024,768, OF_WINDOW);			// setup the GL context: glfwMakeContextCurrent(windowP);
	//	window = ofPtr<ofAppBaseWindow>(new ofAppGLFWWindow()) =>
	//  ofSetupOpenGL(ofPtr<ofAppBaseWindow> windowPtr, int w, int h, int screenMode):
	// window = windowPtr;

	
	// this kicks off the running of my app
	// can be OF_WINDOW or OF_FULLSCREEN
	// pass in width and height too:
	ofRunApp( new ofApp()); // create an object of ofApp and run it

}
