#include "MyUI.h"

//--------------------------------------------------------------
void MyUI::setup(){
	
  //  ofSetCircleResolution(120);
  //  red = 233; blue = 233; green = 233;
  //  hideGUI = false;
  //  bdrawGrid = false;
//	bdrawPadding = false;

//    ddl = NULL;
//    textInput = NULL;
    //img = new ofImage();
   // img->loadImage("nerd_me.png");
    //buffer = new float[256];
    //for(int i = 0; i < 256; i++) { buffer[i] = ofNoise(i/100.0); }
    
	setGUI1();
	/*setGUI2();
    setGUI3();
    setGUI4();
    setGUI5();*/
    
    //gui1->loadSettings("gui1Settings.xml");
 /*   gui2->loadSettings("gui2Settings.xml");
    gui3->loadSettings("gui3Settings.xml");
    gui4->loadSettings("gui4Settings.xml");
    gui5->loadSettings("gui5Settings.xml");*/
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
    
    
	if(bdrawGrid)
	{
		ofSetColor(255, 255, 255, 25);
		drawGrid(8,8);
	}
    
	ofPopStyle();
    
    ofSetRectMode(OF_RECTMODE_CENTER);
}

void MyUI::guiEvent(ofxUIEventArgs &e)
{
	ofApp::setRedrawWindowEvent(true);
	string name = e.getName();
	int kind = e.getKind();
	cout << "got event from: " << name << endl;
    if(kind == OFX_UI_WIDGET_NUMBERDIALER)
    {
        ofxUINumberDialer *n = (ofxUINumberDialer *) e.widget;
        cout << n->getValue() << endl;
    }
	
    if(name == "SAMPLER")
    {
        ofxUIImageSampler *is = (ofxUIImageSampler *) e.widget;
        ofColor clr = is->getColor();
        red = clr.r;
        blue = clr.b;
        green = clr.g;
    }
	else if(name == "BUTTON")
	{
		ofxUIButton *button = (ofxUIButton *) e.getButton();
		bdrawGrid = button->getValue();
	}
	else if(name == "TOGGLE")
	{
		ofxUIToggle *toggle = (ofxUIToggle *) e.getToggle();
		bdrawGrid = toggle->getValue();
        if(textInput != NULL)
        {
            textInput->setFocus(bdrawGrid);
        }
	}
    else if(name == "RADIO VERTICAL")
    {
        ofxUIRadio *radio = (ofxUIRadio *) e.widget;
        cout << radio->getName() << " value: " << radio->getValue() << " active name: " << radio->getActiveName() << endl; 
    }
    else if(name == "TEXT INPUT")
    {
        ofxUITextInput *ti = (ofxUITextInput *) e.widget;
        if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_ENTER)
        {
            cout << "ON ENTER: ";
        }
        else if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_FOCUS)
        {
            cout << "ON FOCUS: ";
        }
        else if(ti->getInputTriggerType() == OFX_UI_TEXTINPUT_ON_UNFOCUS)
        {
            cout << "ON BLUR: ";
        }
        string output = ti->getTextString();
        cout << output << endl;
    }
}

//--------------------------------------------------------------
void MyUI::exit()
{
    gui1->saveSettings("gui1Settings.xml");
    gui2->saveSettings("gui2Settings.xml");
    gui3->saveSettings("gui3Settings.xml");
    gui4->saveSettings("gui4Settings.xml");
    gui5->saveSettings("gui5Settings.xml");
    
	delete gui1;
	delete gui2;
    delete gui3;
    delete gui4;
    delete gui5;
	delete[] buffer;
    delete img;
}

//--------------------------------------------------------------
void MyUI::keyPressed(int key){
    if(gui2->hasKeyboardFocus())
    {
        return;
    }
	switch (key)
	{
		case 't':
        {
            if(textInput != NULL)
            {
                textInput->setTextString(ofGetTimestampString());
            }
        }
			break;

		case 'T':
        {
            if(tm != NULL)
            {
                int cols = tm->getColumnCount();
                int rows = tm->getRowCount();
                for(int row = 0; row < rows; row++)
                {
                    for(int col = 0; col < cols; col++)
                    {
                        cout << tm->getState(row, col) << "\t";
                    }
                    cout << endl;
                }
            }
        }
			break;

		case 'd':
        {
            if(ddl != NULL)
            {
                vector<ofxUIWidget *> selected = ddl->getSelected();
                for(vector<ofxUIWidget *>::iterator it = selected.begin(); it != selected.end(); ++it)
                {
                    ofxUILabelToggle *lt = (ofxUILabelToggle *) (*it);
                    cout << lt->getName() << endl;
                }
            }
        }
			break;
            
        case 'D':
        {
            if(ddl != NULL)
            {
                vector<string> names = ddl->getSelectedNames();
                for(vector<string>::iterator it = names.begin(); it != names.end(); ++it)
                {
                    cout << (*it) << endl;
                }
            }
        }
			break;
            
		case 'r':
        {
            if(textInput != NULL)
            {
                textInput->setFocus(!textInput->isFocused());
            }
        }
			break;
            
		case 'f':
			ofToggleFullscreen();
			break;
            
        case 'F':
        {
            if(tm != NULL)
            {
                tm->setDrawOutlineHighLight(!tm->getDrawOutlineHighLight());
//                tm->setDrawPaddingOutline(!tm->getDrawPaddingOutline());
            }
        }
			break;
            
		case 'h':
            gui1->toggleVisible();
            gui2->toggleVisible();
            gui3->toggleVisible();
            gui4->toggleVisible();
            gui5->toggleVisible();
			break;
            
		case 'p':
			bdrawPadding = !bdrawPadding;
			gui1->setDrawWidgetPaddingOutline(bdrawPadding);
			gui2->setDrawWidgetPaddingOutline(bdrawPadding);
			gui3->setDrawWidgetPaddingOutline(bdrawPadding);
			gui4->setDrawWidgetPaddingOutline(bdrawPadding);
			gui5->setDrawWidgetPaddingOutline(bdrawPadding);
			break;
            
		case '[':
			gui1->setDrawWidgetPadding(false);
			gui2->setDrawWidgetPadding(false);
			gui3->setDrawWidgetPadding(false);
			gui4->setDrawWidgetPadding(false);
			gui5->setDrawWidgetPadding(false);
			break;
            
		case ']':
			gui1->setDrawWidgetPadding(true);
			gui2->setDrawWidgetPadding(true);
			gui3->setDrawWidgetPadding(true);
			gui4->setDrawWidgetPadding(true);
			gui5->setDrawWidgetPadding(true);
			break;
			
        case '1':
            gui1->toggleVisible();
            break;
            
        case '2':
            gui2->toggleVisible();
            break;
            
        case '3':
            gui3->toggleVisible();
            break;
            
        case '4':
            gui4->toggleVisible();
            break;

        case '5':
            gui5->toggleVisible();
            break;

		default:
			break;
	}
}

void MyUI::drawGrid(float x, float y)
{
    float w = ofGetWidth();
    float h = ofGetHeight();
    
    for(int i = 0; i < h; i+=y)
    {
        ofLine(0,i,w,i);
    }
    
    for(int j = 0; j < w; j+=x)
    {
        ofLine(j,0,j,h);
    }
}

void MyUI::setGUI1()
{
    vector<string> names;
	names.push_back("RAD1");
	names.push_back("RAD2");
	names.push_back("RAD3");
	
	gui1 = new ofxUISuperCanvas("PANEL 1: BASICS");
    gui1->addSpacer();
    gui1->addLabel("Press 'h' to Hide GUIs", OFX_UI_FONT_SMALL);
    
    gui1->addSpacer();
	gui1->addLabel("H SLIDERS");
	gui1->addSlider("RED", 0.0, 255.0, &red)->setTriggerType(OFX_UI_TRIGGER_ALL);
	gui1->addSlider("GREEN", 0.0, 255.0, &green)->setTriggerType(OFX_UI_TRIGGER_BEGIN|OFX_UI_TRIGGER_CHANGE|OFX_UI_TRIGGER_END);
	gui1->addSlider("BLUE", 0.0, 255.0, &blue)->setTriggerType(OFX_UI_TRIGGER_BEGIN|OFX_UI_TRIGGER_CHANGE);
    
    gui1->addSpacer();
    gui1->addLabel("V SLIDERS");
	gui1->addSlider("0", 0.0, 255.0, 150, 17, 160);
	gui1->setWidgetPosition(OFX_UI_WIDGET_POSITION_RIGHT);
	gui1->addSlider("1", 0.0, 255.0, 150, 17, 160);
	gui1->addSlider("2", 0.0, 255.0, 150, 17, 160);
	gui1->addSlider("3", 0.0, 255.0, 150, 17, 160);
	gui1->addSlider("4", 0.0, 255.0, 150, 17, 160);
	gui1->addSlider("5", 0.0, 255.0, 150, 17, 160);
	gui1->addSlider("6", 0.0, 255.0, 150, 17, 160);
	gui1->addSlider("7", 0.0, 255.0, 150, 17, 160);
	gui1->addSlider("8", 0.0, 255.0, 150, 17, 160);
	gui1->setWidgetPosition(OFX_UI_WIDGET_POSITION_DOWN);
    
    gui1->addSpacer();
	gui1->addRadio("RADIO HORIZONTAL", names, OFX_UI_ORIENTATION_HORIZONTAL);
	gui1->addRadio("RADIO VERTICAL", names, OFX_UI_ORIENTATION_VERTICAL);
    
    gui1->addSpacer();
    gui1->setWidgetFontSize(OFX_UI_FONT_SMALL);
	gui1->addButton("BUTTON", false);
	gui1->addToggle( "TOGGLE", false);
    
    gui1->addSpacer();
    gui1->addLabel("RANGE SLIDER");
	gui1->addRangeSlider("RSLIDER", 0.0, 255.0, 50.0, 100.0);
    
    string textString = "This widget is a text area widget. Use this when you need to display a paragraph of text. It takes care of formatting the text to fit the block.";
    gui1->addSpacer();
    
    gui1->addTextArea("textarea", textString, OFX_UI_FONT_SMALL);
    
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