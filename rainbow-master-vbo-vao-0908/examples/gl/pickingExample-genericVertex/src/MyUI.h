#pragma once

#include "ofMain.h"
#include "ofxUI.h"
#include "ofApp.h"



class MyUI : public ofBaseApp {

public:
	void setup();
    void update();
    void draw();
    void exit();

    void keyPressed(int key);
    void keyReleased(int key);
    void mouseMoved(int x, int y );
    void mouseDragged(int x, int y, int button);
    void mousePressed(int x, int y, int button);
    void mouseReleased(int x, int y, int button);
    void windowResized(int w, int h);
    void dragEvent(ofDragInfo dragInfo);
    void gotMessage(ofMessage msg);
    
    void drawGrid(float x, float y);
    
	void setGUI1();
	void setGUI2();
	void setGUI3();
	void setGUI4();
	void setGUI5();
	
	ofxUISuperCanvas *gui1;
	ofxUISuperCanvas *gui2;
	ofxUISuperCanvas *gui3;
    ofxUISuperCanvas *gui4;
    ofxUISuperCanvas *gui5;
    

	bool hideGUI;
	
	float red, green, blue;
	bool bdrawGrid;
	bool bdrawPadding;
	
	void guiEvent(ofxUIEventArgs &e);
    ofxUISuperCanvas *gui_event;


    ofxUIMovingGraph *mg;
    ofxUIDropDownList *ddl;
    ofxUIToggleMatrix *tm;
    
    float *buffer;
    ofImage *img;



	ofxUITextInput *textInput; 
	ofxUITextInput *dx; 
	ofxUITextInput *dy; 
	ofxUITextInput *dz; 
	ofxUITextInput *distanceToLight;
	ofxUITextInput *heightToLight;
	ofxUITextInput *separationLight;
	ofxUITextInput *PositionX;
	ofxUITextInput *PositionY;
	ofxUITextInput *PositionZ;
	ofxUITextInput *LookX;
	ofxUITextInput *LookY;
	ofxUITextInput *LookZ;
	ofxUITextInput *UpX;
	ofxUITextInput *UpY;
	ofxUITextInput *UpZ; 

	ofPoint position;



	float  g_distanceToLight;
	float  g_heightToLight;
	float  g_separationLight;

	float g_Lookat[4][4];
	void setg_lightPos();


};



