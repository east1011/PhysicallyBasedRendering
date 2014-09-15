#pragma once

#include "ofMain.h"  // #include "GL/glew.h"

#include "materialshader.h"  // also #include "GL/glew.h"
#include "geometrymaker.h"


class testApp : public ofBaseApp{
	
	public:

		testApp();

		void setup();
		void update();
		void draw();
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
		
		
		void updateFrustFovY();
		void initMaterials();

		void drawForPicking();
		void renderToFbo();
		void readPixels();
		void pick();
		// void display(); nop. we need to use display() function in ofAppGLFWWindow.cpp
		void drawObjects();
		void drawPseudoColors();
		
		void initGLState();
		
		void initPBRTObjects();
		void initObjects();
		void initRainbow();
		void initAABB();
		void setupSun();
		void setupCamera();
		void setupFBO();
		void setupFloatTexturesForVolume();
		Cvec3  getSunLightRGBColor();		
		Matrix4 makeProjectionMatrix();
			
};

