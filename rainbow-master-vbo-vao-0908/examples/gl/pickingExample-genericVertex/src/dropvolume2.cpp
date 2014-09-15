
#include <iostream>
#include <fstream>
#include <complex> // for complex numbers
#include <cmath> // std::abs is here

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cassert>


#include  "dropvolume.h"

extern Cvec3  g_lightRGBColor; // defined in testApp.cpp

using namespace std;

extern std::ofstream messageFile;


void setupVolume() {

	string fileName1 ="rainbow/scattering2_120_fort.txt";

	readScatCrossFile(fileName1.c_str(),  nRadii, nSpectralSamples );  

	string fileName2 ="rainbow/phaseFunction2_120_100_fort.txt";



	readPhaseFile(fileName2.c_str(), nRadii,  nSpectralSamples, nThetas );

	writeScatCrossFile();
	writePhaseFile();

	writeSpectrumFile();
	writeRadiusFile();
	writeThetaFile();


}

Cvec3 rayIntersectVolume(   double uDropDensity, double radius,
								  Cvec3 rayorig, Cvec3 raydir, Cvec3 sunRay) {

	// find the intersections with the AABB volume box: The ray from the camera  hits the front or top face of the volume and then
	// one of the remaining faces.  
	// It is assumed that one intersection at a corner of the box or two intersection at two faces, or no intersection occurs.
	 
    int hitType; // 0: no hit, 1: hit at one point, and 2: hit at two points
			 
	double tnear = INFINITY;
    
	double t0, t1;
	
	double phi_sun, phi_obs, psi;

	const Sphere *sphere = NULL;
	
	// ray hits the front face or the top face of the drop volume

	
	
   // COMPUTE THE SURFACE COLOR and VOLUMEL COLOR

	double scatCrossSection = computeScatCrossSection( radius,  lambdaStart);
	
	Cvec3 surfaceColor, volumeColor;
	
    // if there's no intersection return the background color which has attenuated by the drop volume if the ray passes through
	// volume

	if (!sphere) { // no surface intersection

		messageFile << "no sphere hit" << endl;

		Cvec3 surfaceColor = Cvec3( 0 );
		
		
		if ( hitType == 2 ) { //  no surface hit and  only volume hit
								
			 // attenuate the backgroundSurfaceColor  throughout the volume
                  
		      double rayTau = uDropDensity * scatCrossSection  * ( t1 -t0 );

			  double rayPathNormal =  (t1 - t0)  * cos (phi_obs);
				
	          double tau_N = uDropDensity * scatCrossSection  * rayPathNormal;

			  messageFile << "tau_N ( better be 1.0) =" << tau_N << endl;

		      surfaceColor = surfaceColor * exp ( - rayTau ); 
			  
			  // compute rainbow volume color
			  // check if the psi angle is less than 130, in which case the rainbow is not computed

			  double psi = acos ( dot ( -raydir, -sunRay ) ) * 180 / PI;

			  if ( psi <  130 || psi > 141.8 ) { // no volume color for rainbow:  The red emerges at a smaller angle psi than the violet
				      //scattering is concentrated at theta_0 = 137.97 but also occurs for angles greater than that. 
				      // This is what causes the white appearance `insidea the primary bow (and `outsidea the secondary bow).
				  //Note the gap between 129 and 138. This is a region of negligible scattering (from higher orders) and appears 
				  //as a dark space between the primary and secondary rainbows

				volumeColor = Cvec3(0);
			  }
			  else {
                 messageFile <<" within rainbow angle: " << psi <<endl;
			     volumeColor = calculate_rainbowColor ( radius, uDropDensity, phi_sun, phi_obs, psi, tau_N);  
			  }

			
		} // passes through volume
		else  { // hitType ==0 or 1 => no surface hit and no volume hit

			
			volumeColor = Cvec3(0);

		}
	
     } // no surface hit

	else { // a surface has been hit

     

	   Cvec3 surfaceDiffuseColor, surfaceSpecularColor; // color of the surface of the object intersected by the ray
	   Cvec3 surfaceColor;
	   
	   
	   Cvec3 phit = rayorig + raydir * tnear; // point of intersection
	
	   Cvec3 nhit = phit - sphere->center; // normal at the intersection point
	   nhit= normalize( nhit ); // normalize normal direction
	 
	
	   Cvec3 halfVec =  normalize( -raydir + sunRay );
	  	   
       double  specularCoefficient  = pow( Max(0.0, dot(halfVec,  nhit) ), 30.0 );
	   Cvec3 surfaceDiffuseCoefficient = sphere->surfaceColor *   Max( double(0), dot( nhit, sunRay ) );
	  
	
	   // process the attenuation of the surface color through volume

	   double scatCrossSection = computeScatCrossSection( radius, lambdaStart); // scattering cross section of a drop is the same
                                                                                // regardless of wavelengths, because absorption is zero


	   if ( hitType == 2 ) { // volume hit and surface hit
			 double sunTau; // optical depth
			 Cvec3 attenuatedLightColor;
		
			 Cvec3 rayEnterPoint = rayorig + raydir * t0; // the point at which the ray hit the volume
			 
			 double rayPathNormal =  norm( phit - rayEnterPoint ) * cos (phi_obs);

			 double sunRayPath = rayPathNormal / cos (phi_sun);	 
			

			 sunTau = uDropDensity * scatCrossSection * sunRayPath; // attenuate from the sphere to the hit
				                                                                         // point of the volume

			 // attenuate the light color through volume: attenuation is the same for each wavelength
             attenuatedLightColor = g_lightRGBColor * exp ( - sunTau ); 

			 // compute the specular surface color
			 surfaceSpecularColor = attenuatedLightColor  * specularCoefficient;
		       
			 // compute the diffuse surface color
	         for (int i=0; i < 3; i++ ) {
		          surfaceDiffuseColor[i] =  surfaceDiffuseCoefficient[i] * attenuatedLightColor[i];
	         }
	
			// attenuate both surface and volume color from the sphere hit point to the entry  point of the volume	

			 double rayPathToEye  =  norm( phit - rayEnterPoint );
		
			 double rayToEyeTau = uDropDensity * scatCrossSection * rayPathToEye; 
			

             surfaceColor = ( surfaceDiffuseColor + surfaceSpecularColor ) * exp ( - rayToEyeTau ) ;
			
			 // compute the volume color

			 // check if the psi angle is less than 130, in which case the rainbow is not computed

			 double psi = acos ( dot ( -raydir, -sunRay ) ) * 180 / PI;
			 
			 // process only the primary rainbow: consider scattering angles from 140 to 137
			 
			 
			 if ( psi <=  130 || psi >= 142 ) { //  secondary bow for rainbow angles 50~53 deg=> scattering angle  130 ~ 127 deg
				volumeColor = Cvec3(0);
			 }
			
             else { // within the rainbow field of view, passing through the volume and hitting the surface
				messageFile <<" within rainbow angle: " << psi <<endl;

			    Cvec3 rayEnterPoint = rayorig + raydir * t0; // the point at which the ray hit the volume
			 
	            double rayPathNormal =  norm( phit - rayEnterPoint ) * cos (phi_obs);
				
	            double tau_N = uDropDensity * scatCrossSection  * rayPathNormal;
				
				messageFile << "tau_N ( better be 1.0) =" << tau_N << endl;


				volumeColor = calculate_rainbowColor (radius, uDropDensity, phi_sun, phi_obs, psi, tau_N ); 

			} // within the rainbow field of view


	   } // if (hitType ==2)

	   else { // surface hit but no volume  hit =>  get the surface color without attenuation.
		   		   						 
			 surfaceSpecularColor = g_lightRGBColor * specularCoefficient;
	         
	         for (int i=0; i < 3; i++ ) {
		          surfaceDiffuseColor[i] =  surfaceDiffuseCoefficient[i] * g_lightRGBColor[i];
	         }
				 
             surfaceColor = ( surfaceDiffuseColor + surfaceSpecularColor );

			 volumeColor = Cvec3(0);
			 
	   } // not hitType ==2
       
	 
	} //  surface hit

	return surfaceColor + volumeColor;


}

Cvec3 calculate_rainbowColor (double radius, double uDropDensity, double phi_sun, double phi_obs, double psi, double tau_N)
{
    double Isun;

             // m
	//get_waterDensity_and_waterDepth(waterDensity, waterDepth);
    
    double lambda_m;


    double irainbow[ nSpectralSamples ];
    
      
	Cvec3 rgbColor = Cvec3(0.0, 0.0, 0.0);
	Cvec3 xyzColor = Cvec3(0.0, 0.0, 0.0);

	double X = 0, Y = 0, Z = 0, XYZ;

   	double lambda;
	int j;

//	double radius = 1.0e-3; // 1mm {-3} meter

			
	//double lambdaStart = 400;   // 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
    //double lambdaEnd = 700;
	//double lambdaStep = 5.0; // 5 manometer
		
	messageFile  << "psi: " << psi << ", phi_obs: " << phi_obs << ", phi_sun:" << phi_sun << endl;

    for (lambda = lambdaStart, j = 0; j < nSpectralSamples; lambda += lambdaStep, j++) {
		
		lambda_m = lambda * 1.0e-9;				// convert lambda to meters from nanometers

        Isun = calculate_irradiance_of_sun(lambda_m); // isun: irradiance per nanometer
		double Lsun = Isun /  pi;

		messageFile << "Lsun at  " << lambda_m <<":" << Lsun << endl;

        double particlePhase = computePhaseFunction( radius,  lambda, psi);
		double avgPhasePerVolume  =  uDropDensity * particlePhase;

		messageFile << "particlePhase =" << particlePhase << endl;
		messageFile << "avgPhasePerVolume =" << avgPhasePerVolume << endl;
		
		        
        irainbow[j] = ( ( Lsun* avgPhasePerVolume  / cos(phi_obs) ) / ( 1/cos(phi_sun) + 1/cos(phi_obs) ) ) * ( 1 - exp(-tau_N * ( 1/cos(phi_sun) + 1/cos(phi_obs) ) ) );
		
		
		messageFile  << "irainbow[" << int( lambda)  << "]=" << irainbow[j] << endl;

	}

	// convert spectral color to rgb color

	double xBar, yBar, zBar;

	/* cie_colour_match[(lambda - 380) / 5][0] = xBar
      cie_colour_match[(lambda - 380) / 5][1] = yBar
      cie_colour_match[(lambda - 380) / 5][2] = zBar
	  ==> cie_coulor_match [] from 380 nanometer to 780 nanometer
    */
	
	// 400 ~ 700 into 120 steps with lamdbaStep - 2.5nm
	for (lambda = lambdaStart, j =0; j < nSpectralSamples; lambda += lambdaStep, j++) {
		
		//intensity = bb_spectrum(lambda);  // You already have the function that computes the irradiance of the sun
		  
	   int i1 = int ( (lambda - 380 ) /  (lambdaStep*2 ) );
	   int i2 = int ( ( lambda + lambdaStep - 380) / (lambdaStep*2 ) );

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

	   	X += irainbow[j] * xBar * lambdaStep; 
		Y += irainbow[j] * yBar * lambdaStep;
		Z += irainbow[j] * zBar * lambdaStep;
		
		
    }
	
		
	XYZ = (X+Y+Z);
	xyzColor = Cvec3( X/XYZ, Y/XYZ, Z/XYZ );
	rgbColor = xyz_to_rgb(cs, xyzColor); // cs: a static variable

	
  /*
	if (constrain_rgb(rgbColor)) {
		norm_rgb(rgbColor);
		//output << rgbColor << "(Approximation)" << endl;
	} else {
		norm_rgb(rgbColor);
		//output << rgbColor << endl;
	}
 */
	constrain_rgb(rgbColor);

    return rgbColor;
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

            if (isdigit(c) || c == '.' || c == 'e' ||  c == 'E' ||c == '-' || c == '+') {
				if ( c == 'E' ) c ='e';

                curNumber[curNumberPos++] = c; // curNumber is an array of character

				inDataLine = true;			  
				}
            			
            else { // non-number character ( '/n' and space  )  => the end of a number: 
				
                curNumber[curNumberPos++] = '\0';

				//messageFile << "parsed number=";
				//messageFile << curNumber << endl;
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
            if (isdigit(c) || c == '.' || c == 'e' ||  c == 'E' ||c == '-' || c == '+' ) {
				if ( c == 'E' ) c ='e';
                curNumber[curNumberPos++] = c;
				inDataLine = true;			  
				}
            			
            else { // non-number character ( '/n' included )  => the end of a number: 
				
                curNumber[curNumberPos++] = '\0';
				
				//messageFile << "parsed number=";
				//messageFile << curNumber << endl;

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
				    
				  
				   particlePhase[ radius_i ][ wavelength_i ][ theta_i ] = realValues[3] / (4 * pi);				
				   n++; // data entry index
			  }

		  
		      inDataLine = false; // once a line of data is read, set the flag inDataLine false because the next line may not be 
		                          // data line

		  } // if ( inDataLine )

		} // if ( end of line )

	

    }  // while (! end of file)

    fclose(f);
    
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
    
    fileName = "rainbow/scattering3.txt";
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


	string fileName = "rainbow/phaseFunction3.txt";
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
					theta << std::setw(15) <<  particlePhase[ i ][ j ] [ k ]  << std::endl;
                

			    theta += thetaStep;
           } // for k

		 

           lambda += lambdaStep;

            
        } // for j
        radius += radiusStep;

	} // for i

	file.close();
	}

