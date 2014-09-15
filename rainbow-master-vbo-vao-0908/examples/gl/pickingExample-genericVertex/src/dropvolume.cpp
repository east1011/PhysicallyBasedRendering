
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <cassert>

//#include "phaseFunction.h"
#include  "dropvolume.h"

// define arrays to be used externally; they are externally declared in dropvolume.h
double Cabs[ nRadii ] [nSpectralSamples];
double FCscat[nRadii][nSpectralSamples];


double particlePhase [nRadii] [nSpectralSamples] [nThetas];
double  radii[ nRadii ];
double  spectralSamples [ nSpectralSamples ];
double  thetas[ nThetas ];


DropVolume::DropVolume( double near, double elevation, double width, double height, double depth,  
		           int dropDensity):
		   mnear (near), melevation (elevation), mwidth (width), mheight (height), mdepth (depth),
		   mdropDensity (dropDensity)  	
 {
		 
		mFrontNormal = Cvec3(0,0,1);
		mBackNormal = Cvec3(0,0,-1);
		mTopNormal = Cvec3( 0,1,0);
		mBottomNormal = Cvec3(0,-1,0);
		mLeftNormal = Cvec3(-1,0,0);
		mRightNormal = Cvec3(0,0,1);

		mFrontCenter = Cvec3(0, mheight/2, 0);
		mBackCenter = Cvec3(0, mheight/2.0, -mdepth);

		mTopCenter = Cvec3( 0, mheight, -mdepth / 2.0);
		mBottomCenter = Cvec3(0, 0, -mdepth/2.0);
		mLeftCenter = Cvec3( - mwidth /2.0, mheight/2.0, - mdepth/2.0);
		mRightCenter = Cvec3( mwidth /2.0, mheight/2.0,  -mdepth/2.0);

		front.leftBottom = Cvec3( -mwidth/2.0, 0, 0);
		front.leftTop = Cvec3( -mwidth/2.0, mheight, 0);
		front.rightBottom = Cvec3(mwidth / 2.0, 0,0);
		front.rightTop = Cvec3( mwidth/2.0, mheight, 0);

		back.leftBottom = Cvec3( -mwidth/2.0, 0, -mdepth);
		back.leftTop = Cvec3( -mwidth/2.0, mheight, -mdepth);
		back.rightBottom = Cvec3(mwidth / 2.0, 0, -mdepth);
		back.rightTop = Cvec3( mwidth/2.0, mheight, -mdepth);

		top.leftBottom = Cvec3( -mwidth/2.0, mheight, 0);
		top.leftTop = Cvec3( -mwidth/2.0, mheight, -mdepth);
		top.rightBottom = Cvec3(mwidth / 2.0, mheight,0);
		top.rightTop = Cvec3( mwidth/2.0, mheight, -mdepth);

		bottom.leftBottom = Cvec3( -mwidth/2.0, 0, 0);
		bottom.leftTop = Cvec3( -mwidth/2.0, 0, -mdepth);
		bottom.rightBottom = Cvec3(mwidth / 2.0, 0, 0);
		bottom.rightTop = Cvec3( mwidth/2.0, 0, -mdepth);


		left.leftBottom = Cvec3( -mwidth/2.0, 0, -mdepth);
		left.leftTop = Cvec3( -mwidth/2.0, mheight, -mdepth);
		left.rightBottom = Cvec3(-mwidth / 2.0, 0, 0);
		left.rightTop = Cvec3( -mwidth/2.0, mheight, 0);

		
		right.leftBottom = Cvec3( mwidth/2.0,  0, -mdepth);
		right.leftTop = Cvec3( mwidth/2.0,   mheight,-mdepth);
		right.rightBottom = Cvec3(mwidth / 2.0,  0,  0 );
		left.rightTop = Cvec3(mwidth/2.0,  mheight,  0);
		

		// check the coordinates
		//
		Assert(   norm( mFrontCenter - ( front.leftBottom + front.rightBottom + front.leftTop + front.rightTop ) / 4.0 )   < 1.0e-8   ) ;
		Assert( norm(mTopCenter - ( top.leftBottom + top.rightBottom + top.leftTop + top.rightTop ) / 4.0) < 1.0e-8 ) ;
		Assert( norm(mBackCenter - ( back.leftBottom + back.rightBottom + back.leftTop + back.rightTop ) / 4.0) < 1.0e-8 ) ;
		Assert( norm(mBottomCenter - ( bottom.leftBottom + bottom.rightBottom + bottom.leftTop + bottom.rightTop ) / 4.0) < 1.0e-8 ) ;
		Assert( norm(mLeftCenter - ( left.leftBottom + left.rightBottom + left.leftTop + left.rightTop ) / 4.0) < 1.0e-8 ) ;
		Assert( norm(mRightCenter - ( right.leftBottom + right.rightBottom + right.leftTop + right.rightTop ) / 4.0) < 1.0e-8 ) ;
 }

					   					   				   

double mix(const double &a, const double &b, const double &mix)
{
	return b * mix + a * ( double(1) - mix);
}


void readScatCrossFile(const char *filename, int nRadii, int nSpectralSamples ) {

	//extern FILE *debugFile;

	FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Unable to open file \"%s\"", filename);
        return;
    }

    int c; 

    bool inNumber = false;

	bool atHeaderLine = true;
	bool inDataLine = false;

    char curNumber[32];
    int curNumberPos = 0; 
    int lineNumber = 1;

	int intValues[2];
	double realValues[3];
	int intValueIndex =0;
	int realValueIndex =0;



	int n = 0; // data entry index

	int radius_i =0;
	int wavelength_i =0;
	
    while ( ( c = getc(f) ) != EOF )  { 

	
		//fprintf(debugFile, "char = %c  ", c); by Moon Jung
		
		// whether the  end of line or not, the followig is executed

		if ( inNumber ) { // during parsing of a number; otherwise, beginning of a number
			// char is part of a number

            if (isdigit(c) || c == '.' || c == 'e' || c == '-' || c == '+') {
                curNumber[curNumberPos++] = c; // curNumber is an array of character

				inDataLine = true;			  
				}
            			
            else { // non-number character ( '/n' and space  )  => the end of a number: 
				
                curNumber[curNumberPos++] = '\0';
				
				//fprintf(debugFile, "number =%s \n", curNumber); //by Moon Jung

				// the first line of numbers is supposed to be the header line of the file 
		        // the numbers reached while atHeaderLine is false ARE HEADER NUMBERS
				
				if ( atHeaderLine ) { // the numbers for the header line are read
						
						intValues[ intValueIndex ++ ] =  atoi( curNumber );
						//cout << "number =" << atoi(curNumber) << " ending at" << " curNumberPos =" << curNumberPos << endl;
				}
				else { // non header line numbersr
					    realValues[ realValueIndex ++ ] =  atof(curNumber);
						//cout << "number = " << atof(curNumber) << " ending at" << " curNumberPos =" << curNumberPos << endl;

				}

				
			

                Assert( curNumberPos < (int)sizeof(curNumber) );

                inNumber = false;

                curNumberPos = 0;
							
               } // a number is parsed

        } // if ( inNumber )

		// not yet part of a number, but the beginning of a line

        else { // The beginning of a number

			//cout << "beginning char=" << c << endl << endl;
			 
            if (isdigit(c) || c == '.' || c == '-' || c == '+') { // ascii char == integer 

                inNumber = true;
				curNumber[curNumberPos++] = c;
            }
           
			else if (c == '#') {
                while ((c = getc(f)) != '\n' && c != EOF)  // skip all the characters except for 
					                                       // the end of a line or the end of file
                    ;
                ++lineNumber;
            }


            else if (  c != 32 )  { // other characters are supposed to be a space

				
				cout << "not  space  char=" << c << endl << endl;
                fprintf(stderr, "Unexpected text found at line %d of float file \"%s\n",
                        lineNumber, filename);
				
            }

			// c == 32 
			continue;
			

        } // the beginning of the a number


	
        if (c == '\n') { //  the end of a line is read, and  the end of line specially treated. 
			             
		  ++lineNumber;  

		  if ( inDataLine ) { // a line of numbers is read; # comment is skipped

			  intValueIndex =0;
			  realValueIndex =0;

			   // check if the read header line numbers are equal to the provided sizes of the arrays
			   
			  if ( atHeaderLine ) { // we are at the header line
				  // vector push_back() eliminates duplication by default. 

				  Assert( intValues[0] == nRadii && intValues[1] == nSpectralSamples  );
				  atHeaderLine = false;
				  
			  }

			  else { // in pure data lines => store the read data to the associated  arrays
				     // phase(i,j,k) is stored in hash table, with key = (i,j,k)
				  
				  // n: data entry index

				   
				   if ( n == 0 ) { // the first data entry
					   
					   radius_i = 0;
					   wavelength_i = 0;
	
					   radii[ radius_i ] = realValues[0];
					   spectralSamples [ wavelength_i ] = realValues[1];

					  
				   } // n is zero


				   else { //  n is not zero
					   
				       if (  n %  nSpectralSamples == 0 ) { // n == k * nSpectralSamples, k>=1 
					       radius_i ++;                                              
					       wavelength_i = 0;

					       radii [ radius_i  ] = realValues[0];				  
				       }

				       if (   n %  nSpectralSamples  != 0 ) { // n != k * nSpectralSamples, k >= 1
					    					 						
						   wavelength_i ++;  // wavelength_i is reset to zero when n == k * nSpectralSamples

						   if ( radius_i == 0 )  {
							   spectralSamples [ wavelength_i ] = realValues[1];		
						   }

				       }
				  
				  	
					
				   } // n is not zero 			  	

				   //  FCscat
				   FCscat[ radius_i  ][ wavelength_i ]  = realValues[2];
				  				
				   n++; // data entry index
			  }
			  		  

		  } // if ( inDataLine )

	 
		} // if ( end of line )

	
    }  // while 

    fclose(f);
    
}


void readPhaseFile(const char *filename, int nRadii, int nSpectralSamples, int nThetas ) {

	//extern FILE *debugFile;

	FILE *f = fopen(filename, "r");
    if (!f) {
        fprintf(stderr, "Unable to open file \"%s\"", filename);
        return;
    }

    int c;
    bool inNumber = false;

	bool atHeaderLine = true;
	bool inDataLine = false;

    char curNumber[32];
    int curNumberPos = 0; 
    int lineNumber = 1;


	int intValues[3];
	double realValues[4];

	int intValueIndex =0;
	int realValueIndex =0;

	
	int n = 0; // data entry index

	int radius_i =0;
	int wavelength_i =0;
	int theta_i =0;


    while ( ( c = getc(f) ) != EOF ) { 

		//fprintf(debugFile, "char = %c  ", c); by Moon Jung

        if (inNumber) { // during parsing of a number
			// char is part of a number
            if (isdigit(c) || c == '.' || c == 'e' || c == '-' || c == '+' ) {
                curNumber[curNumberPos++] = c;
				inDataLine = true;			  
				}
            			
            else { // non-number character ( '/n' included )  => the end of a number: 
				
                curNumber[curNumberPos++] = '\0';
				
				//fprintf(debugFile, "number =%s \n", curNumber); //by Moon Jung

				// the first line of numbers is supposed to be the header line of the file 
		        // the numbers reached while atHeaderLine is false ARE HEADER NUMBERS
				
				if ( atHeaderLine ) { // the numbers for the header line are read
						
						intValues[ intValueIndex++ ] =  atoi( curNumber );
				}
				else { // non header line numbers
					    realValues[ realValueIndex++ ] =  atof(curNumber);
				}

				

                Assert( curNumberPos < (int)sizeof(curNumber) );
                inNumber = false;
                curNumberPos = 0;
							
               } // a number is parsed

        } // if ( inNumber )

		// not yet part of a number, but the beginning of a line

        else {
			// The beginning of a number
            if (isdigit(c) || c == '.' || c == '-' || c == '+') {

                inNumber = true;
                curNumber[curNumberPos++] = c;
            }
            else if (c == '#') {
                while ((c = getc(f)) != '\n' && c != EOF)  // skip all the characters except for 
					                                       // the end of a line or the end of file
                    ;
                ++lineNumber;
            }
            else if ( !isspace(c)  )  { // other characters are supposed to be a space
				fprintf(stderr, "char=%s", c);

                fprintf(stderr, "Unexpected text found at line %d of float file \"%s\"",
                        lineNumber, filename);
            }
        }

		// process the end of line

        if (c == '\n') {
		  ++lineNumber;  

		  if ( inDataLine ) {

			   // check if the read header line numbers are equal to the provided sizes of the arrays
			  intValueIndex =0;
			  realValueIndex =0; 

			  if ( atHeaderLine ) {
				  Assert( intValues[0] == nRadii && intValues[1] == nSpectralSamples && intValues[2] == nThetas );
				  atHeaderLine = false;
			  }

			  else { // in pure data lines => store the read data to the associated  arrays
				     // phase(i,j,k) is stored in hash table, with key = (i,j,k)
				  
				  // n: data entry index
				  			  

				   if ( n == 0 ) { // the first data entry
					
					   radius_i = 0;
					   wavelength_i = 0;
					   theta_i = 0;
					
					   radii[ radius_i ] = realValues[0];
					   spectralSamples[ wavelength_i ] = realValues[1];
					   thetas[ theta_i ] = realValues[2];
					   					   				   
					   
				   } // n is zero


				   else { //  n is not zero

				     if (   n %  ( nSpectralSamples * nThetas ) == 0  ) { 
						 radius_i ++;
					     wavelength_i = 0;
					     theta_i = 0;

                         radii [ radius_i  ] = realValues[0]; 						 		 
					  					  
				     }
				  
				  				
				    if (  n % ( nSpectralSamples * nThetas ) != 0 ) { // within the subarray of nSpectralSamples x nThetas 

					   if ( n % nThetas == 0 ) {
						  wavelength_i ++; 
						  theta_i = 0;

						  if ( radius_i == 0 ) {
							  spectralSamples[ wavelength_i ] = realValues[1];
						  }

					   }
					  
					   if ( n % nThetas != 0 ) {

						   theta_i ++;

						   if ( radius_i == 0 && wavelength_i ==0 ) {
							   thetas[ theta_i ] = realValues[2];
						   }

					   }
					} //  if (  n % ( nSpectralSamples * nThetas ) != 0 )

				
				   } // n is not zero 			  	
				    
				  
				   particlePhase[ radius_i ][ wavelength_i ][ theta_i ] = realValues[3];				
				   n++; // data entry index
			  }

		  
		      inDataLine = false; // once a line of data is read, set the flag inDataLine false because the next line may not be 
		                          // data line

		  } // if ( inDataLine )

		} // if ( end of line )

	

    }  // while (! end of file)

    fclose(f);
    
}

void setupVolume() {

	string fileName1 ="rainbow/scattering.txt";

	readScatCrossFile(fileName1.c_str(),  nRadii, nSpectralSamples );  

	string fileName2 ="rainbow/phaseFunction.txt";


	readPhaseFile(fileName2.c_str(), nRadii,  nSpectralSamples, nThetas );

	writeScatCrossFile();
	writePhaseFile();

	writeSpectrumFile();
	writeRadiusFile();
	writeThetaFile();


}

void writeRadiusFile()  { 
	
  //std::ofstream output ("output.txt");

    std::ofstream file;
    
    //task3_file << std::setw(15) << "radisu[mm]" << std::setw(20) <<  "wavelength[nm]" <<  std::setw(15) << "scat angle" << std::setw(15) << "phase" << std::endl;
    
  
	string fileName;

	    
    fileName = "rainbow/radii.txt";
	file.open( fileName.c_str() );
    file << std::setiosflags (std::ios::scientific );

	file  << std::setprecision(6) << "radius[mm]"  << std::endl;
	file << nRadii <<  std::endl;

    for (int i=0; i< nRadii; i++ ) {
		              
       		file <<  std::setprecision(6) << radii[i]  << std::endl;
			     
    }
   
	file.close();

}

void writeSpectrumFile()  { 
	
  
    std::ofstream file;
    
    //task3_file << std::setw(15) << "radisu[mm]" << std::setw(20) <<  "wavelength[nm]" <<  std::setw(15) << "scat angle" << std::setw(15) << "phase" << std::endl;
  
	string fileName;
    
    fileName = "rainbow/spectrum.txt";
	file.open( fileName.c_str() );
    file << std::setiosflags (std::ios::scientific );

	file  << std::setprecision(6) << "spectrum"  << std::endl;
	file << nSpectralSamples <<  std::endl;

    for (int i=0; i< nSpectralSamples; i++ ) {
		              
       		file <<  std::setprecision(6) << spectralSamples[i]  << std::endl;
			     
    }
   
	file.close();

}

void writeThetaFile()  { 
	
  
    std::ofstream file;
    
    //task3_file << std::setw(15) << "radisu[mm]" << std::setw(20) <<  "wavelength[nm]" <<  std::setw(15) << "scat angle" << std::setw(15) << "phase" << std::endl;
        
	string fileName;
		    
    fileName = "rainbow/theta.txt";
	file.open( fileName.c_str() );
    file << std::setiosflags (std::ios::scientific );

	file  << std::setprecision(6) << "theta"  << std::endl;
	file << nThetas <<  std::endl;

    for (int i=0; i< nThetas; i++ ) {
		              
       		file <<  std::setprecision(6) << thetas[i]  << std::endl;
			     
    }
   
	file.close();

}


void writeScatCrossFile()  { 
	
  
    std::ofstream file;
    
    //task3_file << std::setw(15) << "radisu[mm]" << std::setw(20) <<  "wavelength[nm]" <<  std::setw(15) << "scat angle" << std::setw(15) << "phase" << std::endl;
    
    double radius, lambda;

	string fileName;

	radius = radiusStart;
    
    fileName = "rainbow/scattering2.txt";
	file.open( fileName.c_str() );
    file << std::setiosflags (std::ios::scientific );

	file  << std::setprecision(6) << "radius[mm]" << std::setw(15) <<  "wavelength[nm]" <<  std::setw(15) << "scatCrossSection" << std::endl;
	file << nRadii << " " << nSpectralSamples <<  std::endl;

    for (int i=0; i< nRadii; i++ ) {

        //string fileName = "C:/Users/moon/pbrt-v2-src/pbrt-v2/pbrt-scenes/rainbow/scatradius";
		//fileName += to_string( i+1);

      
      //  file <<  std::setw(15) << "wavelength[m]" << std::setw(15) << "Cabs" << std::endl;
        
        lambda = lambdaStart;
		for (int j =0; j < nSpectralSamples; j++ )  {
            
			//FCscat[ i ][ j ] = Fscat_crossSection(lambda, radius );
                        
			file <<  std::setprecision(6) << radius << std::setw(15) << lambda  << std::setw(15) << FCscat[ i ][ j ] << std::endl;
			                       
            lambda += lambdaStep;
            
		}

		
        radius += radiusStep;
      
    }
   
	file.close();

}


void writePhaseFile() {
	 
	std::ofstream file;

  
	
	double radius, lambda;
   
	double theta;
	 
	
	radius = radiusStart;


	string fileName = "rainbow/phaseFunction2.txt";
    file.open( fileName.c_str() );

	file << std::setiosflags (std::ios::scientific );

	file  << "radius[mm]" << std::setw(15) <<  "wavelength[nm]" <<  std::setw(15) << "scatteringAngle" << std::setw(15) << "particlePhase" << std::endl;
	file << nRadii << " " << nSpectralSamples << " " << nThetas << std::endl;

    for ( int i=0; i  < nRadii; i++ ) {
        
        lambda = lambdaStart;
        for (int j=0; j < nSpectralSamples ; j++ ) {          
           
		   theta = thetaStart;
           
		   for (int k =0; k < nThetas; k++ ) {
                
		       // particlePhase[ i ][ j ][ k ] = particlePhaseFunction (lambda, radius, theta, FCscat[ i ][ j ] );
                                           
				file << std::setw(15) << std::setprecision( 6) << radius << std::setw(15) << lambda << std::setw(15) <<  
					theta * 180./pi << std::setw(15) <<  particlePhase[ i ][ j ] [ k ]  << std::endl;
                

			    theta += thetaStep;
           } // for k

		 

           lambda += lambdaStep;

            
        } // for j
        radius += radiusStep;

	} // for i

	file.close();
	}


bool intersectLineFace( Cvec3 rayorig, Cvec3 raydir, Cvec3 faceCenter, Cvec3 faceNormal, Face face,  double *t) {
 
	double normalDistToCenter = dot ( faceCenter, faceNormal);
	double normalDistToRayOrig = dot( rayorig, faceNormal );
	double distFromRayOrigToFace = normalDistToCenter - normalDistToRayOrig;

	double normalRatio = dot( raydir, faceNormal);

	// distFromRayOrigToFace = distToFaceAlongRay * normalRatio

	double distToFaceAlongRay = distFromRayOrigToFace / normalRatio;

	Cvec3 hitPoint = rayorig + raydir * distToFaceAlongRay;
	
	// check if hitPoint is within the face, at  corners of the face, or outside of the face
    // Project the corner points of the face and the hitPoint along the normal direction to the face
	// and test if the projected hitPoint is within the projected face.
	
	Cvec3 leftEdgeDir = normalize( face.leftTop - face.leftBottom );
	Cvec3 bottomEdgeDir = normalize( face.rightBottom - face.leftBottom );

	double hitPointAlongLeftEdge = dot( hitPoint, leftEdgeDir );
	double hitPointAlongBottomEdge = dot( hitPoint, bottomEdgeDir );

	double leftBottomAlongLeftEdge = dot(face.leftBottom, leftEdgeDir);
	
	double leftTopAlongLeftEdge = dot(face.leftTop, leftEdgeDir);

	double leftBottomAlongBottomEdge = dot( face.leftBottom, bottomEdgeDir );
	double rightBottomAlongBottomEdge = dot( face.rightBottom, bottomEdgeDir );


	if ( hitPointAlongLeftEdge >= leftBottomAlongLeftEdge && hitPointAlongLeftEdge <= leftTopAlongLeftEdge &&
		 hitPointAlongBottomEdge >= leftBottomAlongBottomEdge && hitPointAlongBottomEdge <= rightBottomAlongBottomEdge ) {
		  
			 return true;
	}
    else return false;
	
}

// This is the main trace function. It takes a ray as argument (defined by its origin
// and direction). We test if this ray intersects any of the geometry in the scene.
// If the ray intersects an object, we compute the intersection point, the normal
// at the intersection point, and shade this point using this information.
// Shading depends on the surface property (is it transparent, reflective, diffuse).
// The function returns a color for the ray. If the ray intersects an object that
// is the color of the object at the intersection point, otherwise it returns
// the background color.


DropVolume DropVolume::rotateVolume( Matrix4 rotMat, DropVolume volumeBox) {

	DropVolume g_volumeBox = volumeBox;

	g_volumeBox.front.leftBottom = Cvec3( rotMat * Cvec4( volumeBox.front.leftBottom, 1.0) ) ;
	g_volumeBox.front.rightBottom = Cvec3( rotMat * Cvec4( volumeBox.front.rightBottom, 1.0) ) ;

	g_volumeBox.front.leftTop = Cvec3( rotMat * Cvec4( volumeBox.front.leftTop, 1.0) ) ;
	g_volumeBox.front.rightTop = Cvec3( rotMat * Cvec4( volumeBox.front.rightTop, 1.0) ) ;


	g_volumeBox.back.leftBottom = Cvec3( rotMat * Cvec4( volumeBox.back.leftBottom, 1.0) ) ;
	g_volumeBox.back.rightBottom = Cvec3( rotMat * Cvec4( volumeBox.back.rightBottom, 1.0) ) ;

	g_volumeBox.back.leftTop = Cvec3( rotMat * Cvec4( volumeBox.back.leftTop, 1.0) ) ;
	g_volumeBox.back.rightTop = Cvec3( rotMat * Cvec4( volumeBox.back.rightTop, 1.0) ) ;


	g_volumeBox.left.leftBottom = Cvec3( rotMat * Cvec4( volumeBox.left.leftBottom, 1.0) ) ;
	g_volumeBox.left.rightBottom = Cvec3( rotMat * Cvec4( volumeBox.left.rightBottom, 1.0) ) ;

	g_volumeBox.left.leftTop = Cvec3( rotMat * Cvec4( volumeBox.left.leftTop, 1.0) ) ;
	g_volumeBox.left.rightTop = Cvec3( rotMat * Cvec4( volumeBox.left.rightTop, 1.0) ) ;

	
	g_volumeBox.right.leftBottom = Cvec3( rotMat * Cvec4( volumeBox.right.leftBottom, 1.0) ) ;
	g_volumeBox.right.rightBottom = Cvec3( rotMat * Cvec4( volumeBox.right.rightBottom, 1.0) ) ;

	g_volumeBox.right.leftTop = Cvec3( rotMat * Cvec4( volumeBox.right.leftTop, 1.0) ) ;
	g_volumeBox.right.rightTop = Cvec3( rotMat * Cvec4( volumeBox.right.rightTop, 1.0) ) ;

	
	g_volumeBox.top.leftBottom = Cvec3( rotMat * Cvec4( volumeBox.top.leftBottom, 1.0) ) ;
	g_volumeBox.top.rightBottom = Cvec3( rotMat * Cvec4( volumeBox.top.rightBottom, 1.0) ) ;

	g_volumeBox.top.leftTop = Cvec3( rotMat * Cvec4( volumeBox.top.leftTop, 1.0) ) ;
	g_volumeBox.top.rightTop = Cvec3( rotMat * Cvec4( volumeBox.top.rightTop, 1.0) ) ;

	
	g_volumeBox.bottom.leftBottom = Cvec3( rotMat * Cvec4( volumeBox.bottom.leftBottom, 1.0) ) ;
	g_volumeBox.bottom.rightBottom = Cvec3( rotMat * Cvec4( volumeBox.bottom.rightBottom, 1.0) ) ;

	g_volumeBox.bottom.leftTop = Cvec3( rotMat * Cvec4( volumeBox.bottom.leftTop, 1.0) ) ;
	g_volumeBox.bottom.rightTop = Cvec3( rotMat * Cvec4( volumeBox.bottom.rightTop, 1.0) ) ;
	
	return g_volumeBox;


}


Cvec3  calculate_background_surfaceColor( DropVolume  volumeBox, double t0, double t1  )  {
	double radius = 1.0e-3; // 1mm {-3} meter
	
    for (double lambda = lambdaStart, j = 0; lambda < lambdaEnd; lambda += lambdaStep, j++) {
		
		double lambda_m = lambda * 1.0e-9;				// convert lambda to meters from nanometers

        double isun = calculate_irradiance_of_sun(lambda_m); // isun: irradiance per nanometer

        double scatCrossSection = computeScatCrossSection( radius,  lambda);
        
		rayTau = volumeBox.mdropDensity * scatCrossSection n * ( t1 -t0 );
		
		* exp ( - rayTau ); 


}


 
SurfVolColor  DropVolume::rayIntersect(  const std::vector< Sphere * > &spheres, DropVolume volumeBox, 
								  Cvec3 rayorig, Cvec3 raydir, Cvec3 g_sunRay) {


	//if (raydir.length() != 1) std::cerr << "Error " << raydir << std::end
	
	// find the intersections with the AABB volume box: Rays hit the front or top face of the volume and then
	// one of the remaining faces.  
	// It is assumed that one intersection at a corner of the box or two intersection at two faces, or no intersection occurs.
	 
    int hitType; // 0: no hit, 1: hit at one point, and 2: hit at two points

	
	 
	double tnear = INFINITY;
    
	double t0, t1;


	const Sphere *sphere = NULL;
	

	if ( intersectLineFace( rayorig, raydir, volumeBox.mFrontCenter, volumeBox.mFrontNormal, volumeBox.front, &t0 )) { // first hit to front face
		 
		if ( intersectLineFace( rayorig, raydir, volumeBox.mBackCenter, volumeBox.mBackNormal, volumeBox.back, &t1 ) ) {
			if ( (t1 - t0) < 1.0e-8 )  hitType = 2;
			else hitType = 1; // no genuine second hit => hit only a corner of volume
		} 
		else if (intersectLineFace( rayorig, raydir, volumeBox.mLeftCenter, volumeBox.mLeftNormal, volumeBox.left, &t1 ) ) {
			if ( (t1 - t0) <= 1.0e-8 )  hitType = 2;
			else hitType = 1;
			
		}
		else if (intersectLineFace( rayorig, raydir, volumeBox.mRightCenter, volumeBox.mRightNormal, volumeBox.right, &t1 ) ) {
			if ( (t1 - t0) <= 1.0e-8 )  hitType = 2;
			else hitType = 1;
			
		}
		else if ( intersectLineFace( rayorig, raydir, volumeBox.mBottomCenter, volumeBox.mBottomNormal, volumeBox.bottom, &t1 ) ) {
			if ( (t1 - t0) <= 1.0e-8 )  hitType = 2;
			else hitType = 1;
		}
		else if ( intersectLineFace( rayorig, raydir, volumeBox.mTopCenter, volumeBox.mTopNormal, volumeBox.top, &t1 ) ) {
			if ( (t1 - t0) <= 1.0e-8 )  hitType = 2;
			else hitType = 1;
			
		}
		else { // no 2nd hit 
			hitType = 1;  // it means that the ray hits one of the corners of the face. 
			
		}
				
	} // if: first hit the front face

	else if ( intersectLineFace( rayorig, raydir, volumeBox.mTopCenter, volumeBox.mTopNormal, volumeBox.front, &t0 )) { // first hit to top face
		 
		if ( intersectLineFace( rayorig, raydir, volumeBox.mBackCenter, volumeBox.mBackNormal, volumeBox.back, &t1 ) ) {
			if ( (t1 - t0) < 1.0e-8 )  hitType = 2;
			else hitType = 1; // no genuine second hit => hit only a corner of volume
		} 
		else if (intersectLineFace( rayorig, raydir, volumeBox.mLeftCenter, volumeBox.mLeftNormal, volumeBox.left, &t1 ) ) {
			if ( (t1 - t0) < 1.0e-8 )  hitType = 2;
			else hitType = 1; // no second hit => hit only a corner of volume
			
		}
		else if (intersectLineFace( rayorig, raydir, volumeBox.mRightCenter, volumeBox.mRightNormal, volumeBox.right, &t1 ) ) {
			if ( (t1 - t0) < 1.0e-8 )  hitType = 2;
			else hitType = 1; // no second hit => hit only a corner of volume
			
		}
		else if ( intersectLineFace( rayorig, raydir, volumeBox.mBottomCenter, volumeBox.mBottomNormal, volumeBox.bottom, &t1 ) ) {
			if ( (t1 - t0) < 1.0e-8 )  hitType = 2;
			else hitType = 1; // no genuine second hit => hit only a corner of volume
			
		}
		
		else { // pass by the volume
			hitType = 1;  // it means that the ray hits one of the corners of the face. 
			
		}
				
	} // else if
	
	else { // no hit
		 hitType = 0;
	}

	// find intersection of this ray with the sphere in the scene

	for (unsigned i = 0; i < spheres.size(); ++i) {

		double l0 = INFINITY, l1 = INFINITY;

		if (spheres[i]->intersect(rayorig, raydir, &l0, &l1)) {
			if (l0 < 0) l0 = l1;
			if (l0 < tnear) {
				tnear = l0;
				sphere = spheres[i];
			}
		}
	}

	// check the intersection with the drop volume
	   //phi_obs = PI - acos(normal_waterbox.dot(raydir));       // in radians

   
	double phi_obs = acos( dot( volumeBox.mFrontNormal, -raydir) );       // in radians
    double phi_sun = acos( dot( volumeBox.mFrontNormal, g_sunRay) );          // in radians
    
	// COMPUTE THE SURFACE COLOR and VOLUMEL COLOR

	Cvec3 surfaceColor, volumeColor;

    // if there's no intersection return the background color which has attenuated by the drop volume if the ray passes through
	// volume

	if (!sphere) { // no surface intersection
		Cvec3 backgroundSurfaceColor = Cvec3( 0 );
		
		volumeColor = Cvec3(0); // to be used when there is no volume intersection

		if ( hitType == 2 ) { // light passes through the volume
			 double rayTau; // optical depth
	
			 double scat
			
			 // attenuate throughout the volume
             surfaceColor  = calculate_background_surfaceColor (volumeBox, t0, t1);
			

			 
			 // check if the psi angle is less than 130, in which case the rainbow is not computed

			 double psi = acos ( dot ( -raydir, -g_sunRay ) ) * 180 / PI;
			
			 if ( psi <  130 ) { 
				volumeColor = Cvec3(0) * exp( -rayTau) ;
			 }
			
             else { // within the rainbow field of view
				double tau_N = 	volumeBox.mdropDensity * volumeBox.mscatCrossSection * volumeBox.mdepth;

				volumeColor = calculate_rainbowColor ( volumeBox,  raydir,  g_sunRay);  
			}
		} // passes through volume
			 
		
     } // no surface hit

	else { // a surface has been hit

	   Cvec3 surfaceDiffuseColor, surfaceSpecularColor; // color of the surface of the object intersected by the ray
	   Cvec3 surfaceColor;
	   
	   Cvec3 volumeColor = Cvec3(0); // to be used when there is  no volumen intersection

	   Cvec3 phit = rayorig + raydir * tnear; // point of intersection
	
	   Cvec3 nhit = phit - sphere->center; // normal at the intersection point
	   nhit= normalize( nhit ); // normalize normal direction
	 
	
	   Cvec3 halfVec =  -raydir + g_sunRay;
	   halfVec = normalize( halfVec );
	   
       double  specularCoefficient  = pow( Max(0.0, dot(halfVec,  nhit) ), 30.0 );
	   Cvec3 surfaceDiffuseCoefficient = sphere->surfaceColor *   Max( double(0), dot( nhit, g_sunRay ) );
	  
	
	   // process the attenuation through volume

	   if ( hitType == 2 ) { // light passes through the volume
			 double sunTau; // optical depth
			 Cvec3 attenuatedLightColor;
		
			 Cvec3 rayEnterPoint = rayorig + raydir * t0; // the point at which the ray hit the volume
			 
			 double rayPathNormal =  norm( phit - rayEnterPoint ) * cos (phi_obs);

			 double sunRayPath = rayPathNormal / cos (phi_sun);
			 
			 sunTau = volumeBox.mdropDensity * volumeBox.mscatCrossSection * sunRayPath; // attenuate from the sphere to the hit
				                                                                         // point of the volume
             attenuatedLightColor = g_light_rgbColor * exp ( - sunTau );

			 surfaceSpecularColor = attenuatedLightColor  * specularCoefficient;
		       

	         for (int i=0; i < 3; i++ ) {
		          surfaceDiffuseColor[i] =  surfaceDiffuseCoefficient[i] * attenuatedLightColor[i];
	         }
	
			 double rayTau; // optical depth

					 
			 double rayPath =  norm( phit - rayEnterPoint );
		
			 rayTau = volumeBox.mdropDensity * volumeBox.mscatCrossSection * rayPath; // attenuate from the sphere to the hit
				                                                                         // point of the volume
			 double normalRayPath = rayPath * cos (phi_obs);

			 double tau_N = volumeBox.mdropDensity * volumeBox.mscatCrossSection * normalRayPath;

             surfaceColor = ( surfaceDiffuseColor + surfaceSpecularColor ) * exp ( - rayTau ) ;
			
			 // check if the psi angle is less than 130, in which case the rainbow is not computed

			 double psi = acos ( dot ( -raydir, -g_sunRay ) ) * 180 / PI;

			 if ( psi <  130 ) { 
				volumeColor = Cvec3(0) * exp( -rayTau) ;
			 }
			
             else { // within the rainbow field of view

				
				volumeColor = calculate_rainbow ( volumeBox,  raydir,  g_sunRay); 
			} // within the rainbow field of view


	   } // if (hitType ==2)

	   else { // surface hit but no volume  hit =>  get the surface color without attenuation.
		   		   						 
			 surfaceSpecularColor = g_light_rgbColor * specularCoefficient;
	         
	         for (int i=0; i < 3; i++ ) {
		          surfaceDiffuseColor[i] =  surfaceDiffuseCoefficient[i] * g_light_rgbColor[i];
	         }
				 
             surfaceColor = ( surfaceDiffuseColor + surfaceSpecularColor );
			 
	   } // not hitType ==2
       
	  
	} //  surface hit

	SurfVolColor surfaceVolumeColor;

	surfaceVolumeColor.surfaceColor = surfaceColor;
	surfaceVolumeColor.volumeColor = volumeColor;

	return surfaceVolumeColor;

}


void setSunLightRGBColor() {
// convert the sun light to its RGB representation
	
	double lambda, sunIntensity, X = 0, Y = 0, Z = 0, XYZ;
	double lambda_m; // meter

	double xBar, yBar, zBar;

	/* cie_colour_match[(lambda - 380) / 5][0] = xBar
      cie_colour_match[(lambda - 380) / 5][1] = yBar
      cie_colour_match[(lambda - 380) / 5][2] = zBar
	  ==> cie_coulor_match [] from 380 nanometer to 780 nanometer
    */
	
	for (lambda = lambdaStart; lambda < lambdaEnd; lambda += lambdaStep) {
		
		//intensity = bb_spectrum(lambda);  // You already have the function that computes the irradiance of the sun

		xBar = cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep)] [0]; // it is assumed that lambdaStep is 5 nanometer
		yBar = cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep)] [1];
		zBar = cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep)] [2];

		lambda_m = lambda * 1.0e-9;

        sunIntensity = calculate_irradiance_of_sun(lambda_m); // isun: irradiance per nanometer, but lambda_m is meter

		X += sunIntensity * xBar * lambdaStep; 
		Y += sunIntensity * yBar * lambdaStep;
		Z += sunIntensity * zBar * lambdaStep;
		
		
    }
	
	XYZ = (X+Y+Z);
	g_light_xyzColor = Cvec3( X/XYZ, Y/XYZ, Z/XYZ );
	g_light_rgbColor = xyz_to_rgb(cs, g_light_xyzColor);

	//output << "light_xyzColor: " << g_light_xyzColor << ", light_rgbColor: " << g_light_rgbColor << endl;
}

///////////////////////////////////////////////////////////
//  calcuate radiance of rainbow for a pixel
//////////////////////////////////////////////////////////

// it is assumed that the scattering angle psi between g_gunRay and newRayDir0
// goes from 130 to 180 in one degree step, with nThetas = 60;
// the phase function is precomputed only for those angles.
// This assumption is satisfied by having the field of view 100 deg
// If the resolution in the vertical direciton is 1000, the difference
// between the adjacent rays is 0.1 deg.
// The field of view in the horizontal direction should be also considered
// to the same as the vertical one, because the phase functions is defined
// for the angle from 130 to 180 deg. 

Cvec3  calculate_rainbow (DropVolume volumeBox, const Cvec3  &raydir, const Cvec3 &g_sunRay)
{
    double isun = 0.0;
    double P;
             // m
	//get_waterDensity_and_waterDepth(waterDensity, waterDepth);
    
    double lambda_m,  psi;


    double irainbow[ nSpectralSamples ];
    
    double phi_sun, phi_obs;
    
    
	Cvec3 rgbColor = Cvec3(0.0, 0.0, 0.0);
	Cvec3 xyzColor = Cvec3(0.0, 0.0, 0.0);

	double X = 0, Y = 0, Z = 0, XYZ;

    psi = 180 * acos( dot(g_sunRay, raydir) ) / PI;             //  -g_sunRay dot -rayDir should be computed, but it is equal to g_sunRay dot rayDir 

    //phi_obs = PI - acos(normal_waterbox.dot(raydir));       // in radians
	phi_obs = acos( dot( volumeBox.mFrontNormal, -raydir) );       // in radians
    phi_sun = acos( dot( volumeBox.mFrontNormal, g_sunRay) );          // in radians
    
   

	double lambda;
	double radius = 1.0e-3; // 1mm {-3} meter

	int j;
	//double lambdaStart = 400;   // 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
    //double lambdaEnd = 700;
	//double lambdaStep = 5.0; // 5 manometer

	
    for (lambda = lambdaStart, j = 0; lambda < lambdaEnd; lambda += lambdaStep, j++) {
		
		lambda_m = lambda * 1.0e-9;				// convert lambda to meters from nanometers

        isun = calculate_irradiance_of_sun(lambda_m); // isun: irradiance per nanometer

        P = computePhaseFunction( radius,  lambda_m, psi);
        
        irainbow[j] = ( ( isun*P / cos(phi_obs) ) / ( 1/cos(phi_sun) + 1/cos(phi_obs) ) ) * ( 1 - exp(-tau_N * ( 1/cos(phi_sun) + 1/cos(phi_obs) ) ) );
		
		//std::cout << "irainbow : " << irainbow << std::endl;

		//output << "psi: " << psi << ", phi_obs: " << phi_obs << ", phi_sun:" << phi_sun << endl;
		//output << "irainbow[" << int( lambda)  << "]=" << irainbow[j] << endl;

	}
	
	double xBar, yBar, zBar;

	/* cie_colour_match[(lambda - 380) / 5][0] = xBar
      cie_colour_match[(lambda - 380) / 5][1] = yBar
      cie_colour_match[(lambda - 380) / 5][2] = zBar
	  ==> cie_coulor_match [] from 380 nanometer to 780 nanometer
    */
	
	for (lambda = lambdaStart, j =0; lambda < lambdaEnd; lambda += lambdaStep, j++) {
		
		//intensity = bb_spectrum(lambda);  // You already have the function that computes the irradiance of the sun

		xBar = cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep)] [0]; // it is assumed that lambdaStep is 5 nanometer
		yBar = cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep)] [1];
		zBar = cie_colour_match[ ( int(lambda) - 380 ) / int( lambdaStep)] [2];

		X += irainbow[j] * xBar * lambdaStep; 
		Y += irainbow[j] * yBar * lambdaStep;
		Z += irainbow[j] * zBar * lambdaStep;
		
		
    }
	
		
	XYZ = (X+Y+Z);
	xyzColor = Cvec3( X/XYZ, Y/XYZ, Z/XYZ );
	rgbColor = xyz_to_rgb(cs, xyzColor);

	//output << "xyzColor: " << xyzColor << ", rgbColor: " << rgbColor << endl;

	if (constrain_rgb(rgbColor)) {
		norm_rgb(rgbColor);
		//output << rgbColor << "(Approximation)" << endl;
	} else {
		norm_rgb(rgbColor);
		//output << rgbColor << endl;
	}
	
    return rgbColor;
}

