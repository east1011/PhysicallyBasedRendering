#pragma once

#include "ofMain.h"  // #include "GL/glew.h"
#include "ofxUI.h"

#include "shadergeometry.h"


#include "matrix4.h"


class ofApp : public ofBaseApp{
	
	public:

		ofApp();

		void setup();
		void update();
		void draw();  // This declaration requires draw() to be implemented by ofApp.cpp,
		              // even though the parent class implements this virtual class
		              // by providing the body.
		void exit();


		void keyPressed  (int key);
		void keyReleased (int key);

		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
		
		ofxUICanvas *gui; 
        ofxUIDropDownList *ddl;
        void guiEvent(ofxUIEventArgs &e);

		void updateFrustFovY();
		void initMaterials();

		
		void drawForPicking();
		
		void renderToFBOAndDump(int renderMode );

		void renderRainbowAndSceneToScreen(); 
		void renderRainbowOnlyToScreen(); 
		void renderAABBRainbowAndSceneToScreen();
		void renderAABBRainbowOnlyToScreen();
		void renderSceneToScreen();
		void renderSceneToFBO();
		void renderRainbowAndSceneToFBO();
		void renderRainbowOnlyToFBO();
		void renderRainbowOnly();
		void drawTextureFromFBOFile();
		void drawTextureFromSysBufferFile();
		void dumpRainbowFBOToFile();
		void dumpRainbowOnlyFBOToFile();
		void dumpSceneFBOToFile();
		void attachTexturesToRainbowFBO();
		void attachTexturesToSceneFBO();
		void readSceneFBOPixels(char *fileName);
		void readRainbowFBOPixels(char *fileName);

		void pick();

		// void display(); nop. we need to use display() function in ofAppGLFWWindow.cpp
		void drawBackgroundStuff();
		void renderSceneToSysBuffer();
		void initNDC();
		void initFullScreen();
		void drawPseudoColors();
		
		void initGLState();
		
		void initPBRTObjects();
		void initObjects();
		void initRainbow();
		void initAABB();
		void initAABBRainbow();
		void setupSun();
		void setupCamera();
		void setupOldCamera();
		void generateFBO();
		
		void setupFloatTexturesForVolume();
		void  getSunLightRGBColor(Cvec3& g_lightRGBColor, Cvec3& g_lightXYZColor);		
		Matrix4 makeProjectionMatrix();

		virtual void myOwnDraw(); 

		
			
};
