#include "MyUI.h"
#include "api.h"
//--------------------------------------------------------------
void MyUI::setup(){
	
  //  ofSetCircleResolution(120);
  //  red = 233; blue = 233; green = 233;
  //  hideGUI = false;
  //  bdrawGrid = false;
//	bdrawPadding = false;

//    ddl = NULL;
//    textInput = NULL;

    
	setGUI1();

}

//--------------------------------------------------------------
void MyUI::update(){
   // mg->addPoint(buffer[0]);
   // for(int i = 0; i < 256; i++) { buffer[i] = ofNoise(i/100.0, ofGetElapsedTimef()); }
}

//--------------------------------------------------------------
void MyUI::draw(){
    ofBackground(red, green, blue, 255);
	
	ofPushStyle();
	ofEnableBlendMode(OF_BLENDMODE_ALPHA);
    

    
	ofPopStyle();
    
    ofSetRectMode(OF_RECTMODE_CENTER);
}

void MyUI::guiEvent(ofxUIEventArgs &e)
{
	ofApp::setRedrawWindowEvent(true);
	string name = e.getName();
	int kind = e.getKind();
	cout << "got event from: " << name << endl;
   

	if(name=="Light position"){

		if(gui_event!=NULL)
				gui_event->disable();
		gui_event=new ofxUISuperCanvas("Input field");
		
		gui_event->addLabel(name);
		gui_event->addLabel("distanceToLight");
		distanceToLight= gui_event->addTextInput("distanceToLight","0");
		distanceToLight->setAutoClear(false);
		gui_event->addLabel("heightToLight");
		heightToLight= gui_event->addTextInput("heightToLight","0");
		heightToLight->setAutoClear(false);
		gui_event->addLabel("separationLight");
		separationLight= gui_event->addTextInput("separationLight","0");
		separationLight->setAutoClear(false);
	
	}
	//Pbrtlookat
	else if(name=="Pbrtlookat"){
		if(gui_event!=NULL)
				gui_event->disable();
		gui_event=new ofxUISuperCanvas("Input field");

		gui_event->addLabel(name);
		gui_event->addLabel("Position x,y,z");
		PositionX= gui_event->addTextInput("PositionX","0");
		PositionX->setAutoClear(false);
		PositionY= gui_event->addTextInput("PositionY","0");
		PositionY->setAutoClear(false);
		PositionZ= gui_event->addTextInput("PositionZ","0");
		PositionZ->setAutoClear(false);

		gui_event->addLabel("Look x,y,z");
		LookX= gui_event->addTextInput("LookX","0");
		LookX->setAutoClear(false);
		LookY= gui_event->addTextInput("LookY","0");
		LookY->setAutoClear(false);
		LookZ= gui_event->addTextInput("LookZ","0");
		LookZ->setAutoClear(false);
		

		gui_event->addLabel("Up x,y,z");
		UpX= gui_event->addTextInput("UpX","0");
		UpX->setAutoClear(false);
		UpY= gui_event->addTextInput("UpY","0");
		UpY->setAutoClear(false);
		UpZ= gui_event->addTextInput("UpZ","0");
		UpZ->setAutoClear(false);



	}
	else if(name=="Remove input field"){
		if(gui_event!=NULL)
			gui_event->disable();
	}
	else if(name=="PositionX"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float PositionX;
			
			istringstream(ti->getTextString()) >> PositionX; 
			g_Lookat[0][0]=PositionX;
			cout<<PositionX;
			cout << "ON ENTER: ";
			
        }
		//istringstream(s) >> f; 
		
	}else if(name=="PositionY"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float PositionY;
			
			istringstream(ti->getTextString()) >> PositionY; 
			cout<<PositionY;
			g_Lookat[0][1]=PositionY;
            cout << "ON ENTER: ";
        }
		//istringstream(s) >> f; 
		
	}
	else if(name=="PositionZ"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float PositionZ;
			
			istringstream(ti->getTextString()) >> PositionZ; 
			cout<<PositionZ;
			g_Lookat[0][2]=PositionZ;
            cout << "ON ENTER: ";
        }
		//istringstream(s) >> f; 
		
	}
	else if(name=="LookX"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float LookX;
			
			istringstream(ti->getTextString()) >> LookX; 
			cout<<LookX;
			g_Lookat[1][0]=LookX;
            cout << "ON ENTER: ";
        }
		//istringstream(s) >> f; 
		
	}
	else if(name=="LookY"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float LookY;
			
			istringstream(ti->getTextString()) >> LookY; 
			cout<<LookY;
			g_Lookat[1][1]=LookY;
            cout << "ON ENTER: ";
        }
		//istringstream(s) >> f; 
		
	}
	else if(name=="LookZ"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float LookZ;
			
			istringstream(ti->getTextString()) >> LookZ; 
			cout<<LookZ;
			g_Lookat[1][2]=LookZ;
            cout << "ON ENTER: ";
        }
		//istringstream(s) >> f; 
		
	}
	else if(name=="UpX"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float UpX;
			
			istringstream(ti->getTextString()) >> UpX; 
			cout<<UpX;
			g_Lookat[2][0]=UpX;
            cout << "ON ENTER: ";
        }
		//istringstream(s) >> f; 
		
	}
	else if(name=="UpY"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float UpY;
			
			istringstream(ti->getTextString()) >> UpY; 
			g_Lookat[2][1]=UpY;
			cout<<UpY;
            cout << "ON ENTER: ";
        }
		//istringstream(s) >> f; 
		
	}
	else if(name=="UpZ"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float UpZ;
			
			istringstream(ti->getTextString()) >> UpZ; 
			g_Lookat[2][2]=UpZ;
			cout<<UpZ;
		//	guiLookAt(g_Lookat[0][0],g_Lookat[0][1],g_Lookat[0][2],g_Lookat[1][0],g_Lookat[1][1],g_Lookat[1][2],g_Lookat[2][0],g_Lookat[2][1],g_Lookat[2][2]);
            cout << "ON ENTER: ";
        }
		//istringstream(s) >> f; 
		
	}
	else if(name=="distanceToLight"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float distanceToLight;
			
			istringstream(ti->getTextString()) >> g_distanceToLight; 
			setg_lightPos();
            cout << "ON ENTER: ";
        }
		//istringstream(s) >> f; 
		
	}
	else if(name=="heightToLight"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {	
			istringstream(ti->getTextString()) >> g_heightToLight; 
			setg_lightPos();
            cout << "ON ENTER: ";
        }
		//istringstream(s) >> f; 
		
	}
	else if(name=="separationLight"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			istringstream(ti->getTextString()) >> g_separationLight; 
			setg_lightPos();
		
            cout << "ON ENTER: ";
        }
		//istringstream(s) >> f; 
		
	}


	else if(name=="Camera translation"){
		if(gui_event!=NULL)
				gui_event->disable();
		gui_event=new ofxUISuperCanvas("Input field");

		gui_event->addLabel(name);
		gui_event->addLabel("dx (Camera translation)");
		dx= gui_event->addTextInput("dx (Camera translation)","0");
		gui_event->addLabel("dy (Camera translation)");
		dy= gui_event->addTextInput("dy (Camera translation)","0");
		gui_event->addLabel("dz (Camera translation)");
		dz= gui_event->addTextInput("dz (Camera translation)","0");
		dz->setAutoClear(false);
		dx->setAutoClear(false);
		dy->setAutoClear(false);
	}
	else if(name=="Camera rotation"){
		if(gui_event!=NULL)
				gui_event->disable();
		gui_event=new ofxUISuperCanvas("Input field");

		gui_event->addLabel(name);
		gui_event->addLabel("dx (Camera rotation)");
		dx= gui_event->addTextInput("dx (Camera rotation)","0");
		
		gui_event->addLabel("dy (Camera rotation)");
		dy= gui_event->addTextInput("dy (Camera rotation)","0");
		
		gui_event->addLabel("dz (Camera rotation)");
		dz= gui_event->addTextInput("dz (Camera rotation)","0");
		dz->setAutoClear(false);
		dx->setAutoClear(false);
		dy->setAutoClear(false);
		
	}

	else if(name=="dx (Camera translation)"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float dx_temp;
			
			istringstream(ti->getTextString()) >> dx_temp; 
			Matrix4 m = Matrix4::makeTranslation( Cvec3( dx_temp, 0.0,0.0) * 0.01 ); 
			ofApp::changeG_Rbt(m);
			cout<<dx_temp;
            cout << "ON ENTER: ";


        }
		//istringstream(s) >> f; 
		
	}
	else if(name=="dy (Camera translation)"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float dy_temp;
			istringstream(ti->getTextString()) >> dy_temp; 
			Matrix4 m = Matrix4::makeTranslation( Cvec3( 0.0, dy_temp,0.0) * 0.01 ); 
			ofApp::changeG_Rbt(m);
            cout << "ON ENTER: ";
        }
		
	}
	else if(name=="dz (Camera translation)"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float dz_temp;
			istringstream(ti->getTextString()) >> dz_temp; 
			Matrix4 m = Matrix4::makeTranslation( Cvec3( 0.0, 0.0,dz_temp) * 0.01 ); 
			ofApp::changeG_Rbt(m);
            cout << "ON ENTER: ";
        }
		
	}

	else if(name=="dx (Camera rotation)"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float dx_temp;
			istringstream(ti->getTextString()) >> dx_temp; 
			Matrix4 m=Matrix4::makeYRotation( dx_temp* 0.01 );
			ofApp::changeG_Rbt(m);
            cout << "ON ENTER: ";
        }
		
	}
	else if(name=="dy (Camera rotation)"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float dy_temp;
			istringstream(ti->getTextString()) >> dy_temp; 
			Matrix4 m=Matrix4::makeXRotation( dy_temp* 0.01 );
			ofApp::changeG_Rbt(m);
            cout << "ON ENTER: ";
        }
		
	}
	else if(name=="dz (Camera rotation)"){
		ofxUITextInput *ti = (ofxUITextInput *) e.widget;
		if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
			float dz_temp;
			istringstream(ti->getTextString()) >> dz_temp; 
			Matrix4 m=Matrix4::makeZRotation( dz_temp* 0.01 );
			ofApp::changeG_Rbt(m);
            cout << "ON ENTER: ";
        }
		
	}

	if(gui_event!=NULL){
		gui_event->setPosition(212*3,0);
		gui_event->autoSizeToFitWidgets();
		ofAddListener(gui_event->newGUIEvent,this,&MyUI::guiEvent);
	}


}

//--------------------------------------------------------------
void MyUI::exit()
{

}

//--------------------------------------------------------------
void MyUI::keyPressed(int key){
 
	cout<<"asf";
}



void MyUI::setGUI1()
{
	gui_event=NULL;
	ddl = NULL;
    textInput = NULL;
	gui1 = new ofxUISuperCanvas("PANEL 1: BASICS");
  

	vector<string> items;
    items.push_back("Pbrtlookat");items.push_back("Light position");items.push_back("Camera translation"); items.push_back("Camera rotation"); 
	items.push_back("Remove input field"); 
	gui1->addDropDownList("DROP DOWN LIST", items);

	


    gui1->autoSizeToFitWidgets();
	ofAddListener(gui1->newGUIEvent,this,&MyUI::guiEvent);
}


//--------------------------------------------------------------
void MyUI::keyReleased(int key){

}

//--------------------------------------------------------------
void MyUI::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void MyUI::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void MyUI::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void MyUI::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void MyUI::windowResized(int w, int h){

}

//--------------------------------------------------------------
void MyUI::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void MyUI::dragEvent(ofDragInfo dragInfo){ 

}


void MyUI::setg_lightPos(){
//	g_light1Pos = g_SceneRbt * Cvec4( g_separationLight/2.0, g_heightToLight, g_distanceToLight, 1.0);
//	g_light2Pos = g_SceneRbt * Cvec4(-g_separationLight/2.0, g_heightToLight, g_distanceToLight, 1.0); 
}