#pragma once

#include "ofPoint.h"
#include "ofEvents.h"
#include "ofBaseTypes.h"

class ofBaseApp : public ofBaseSoundInput, public ofBaseSoundOutput{

	public:
        ofBaseApp() {
            mouseX = mouseY = 0;
        }

		virtual ~ofBaseApp(){
		}

		virtual void setup(){}
		virtual void update(){}
		virtual void draw(){}
		virtual void exit(){}

		virtual void windowResized(int w, int h){}

		virtual void keyPressed( int key ){}
		virtual void keyReleased( int key ){}

		virtual void mouseMoved( int x, int y ){}
		virtual void mouseDragged( int x, int y, int button ){}
		virtual void mousePressed( int x, int y, int button ){}
		virtual void mouseReleased(int x, int y, int button ){}
		
		virtual void dragEvent(ofDragInfo dragInfo) { }
		virtual void gotMessage(ofMessage msg){ }
	
		virtual void windowEntry ( int state ) { }
		
		int mouseX, mouseY;			// for processing heads

		virtual void setup(ofEventArgs & args){
			setup();
		}
		virtual void update(ofEventArgs & args){
			update();
		}
		virtual void draw(ofEventArgs & args){
			draw();
		}
		virtual void exit(ofEventArgs & args){
			exit();
		}

		virtual void windowResized(ofResizeEventArgs & resize){
			windowResized(resize.width,resize.height);
		}

		virtual void keyPressed( ofKeyEventArgs & key ){
			keyPressed(key.key);
		}
		virtual void keyReleased( ofKeyEventArgs & key ){
			keyReleased(key.key);
		}

		virtual void mouseMoved( ofMouseEventArgs & mouse ){
			mouseX=mouse.x;
			mouseY=mouse.y;
			mouseMoved(mouse.x,mouse.y);
		}
		virtual void mouseDragged( ofMouseEventArgs & mouse ){
			mouseX=mouse.x;
			mouseY=mouse.y;
			mouseDragged(mouse.x,mouse.y,mouse.button);
		}
		virtual void mousePressed( ofMouseEventArgs & mouse ){
			mouseX=mouse.x;
			mouseY=mouse.y;
			mousePressed(mouse.x,mouse.y,mouse.button);
		}
		virtual void mouseReleased(ofMouseEventArgs & mouse){
			mouseX=mouse.x;
			mouseY=mouse.y;
			mouseReleased(mouse.x,mouse.y,mouse.button);
		}
		virtual void windowEntry(ofEntryEventArgs & entry){
			windowEntry(entry.state);
		}
		virtual void dragged(ofDragInfo & drag){
			dragEvent(drag);
		}
		virtual void messageReceived(ofMessage & message){
			gotMessage(message);
		}

		virtual void myOwnDraw(bool *isWindowRedrawn) = 0; // by Moon Jung, 2014/8/9
		//For a virtual function you need to provide implementation in the base class. However derived class can override this implementation with its own implementation. Normally , 
		//	for pure virtual functions implementation is not provided. 
        //You can actually provide implementations of pure virtual functions in C++. 
		//The only difference is all pure virtual functions must be implemented
	    // by derived classes before the class can be instantiated.

		// Base class pointer can point to derived class object. 
		// In this case, using base class pointer if we call some function
		//	which is in both classes, then base class function is invoked.
		//	But if we want to invoke derived class function using base class pointer, 
		//	it can be achieved by defining the function as virtual in base class, 
		// this is how virtual functions support runtime polymorphism.

};



