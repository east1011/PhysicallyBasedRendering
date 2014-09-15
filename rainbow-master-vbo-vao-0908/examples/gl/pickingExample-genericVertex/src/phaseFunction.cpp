
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <iomanip>
#include <stdlib.h>     /* exit, EXIT_FAILURE */

#include "phaseFunction.h"

/*
Violet: 400 ~ 420
INDIGO: 420 ~ 400
BLUE: 440 ~ 490e
GREEN: 490 ~ 570
YELLOW: 570 ~ 585
ORANGE: 585 ~ 620
RED: 620 ~ 670

*/
extern std::ofstream messageFile;


// define arrays to be used externally; they are externally declared in dropvolume.h
double Cabs[ nRadii][nLambdas];
double FCscat[nRadii][nLambdas];


double particlePhase [nRadii][nLambdas][nThetas];
double  radii[ nRadii ];
double  spectralSamples [ nLambdas ];
double  thetas[ nThetas ];


interval  binarySearch( double *table, interval indexInterval, double key) {


   // if interval pair contains consecutive indices, it means that key is found to fall within the values at the indices
   
   if ( indexInterval.start + 1  == indexInterval.end  ) return indexInterval;

   int midIndex  = (  indexInterval.start + indexInterval.end )  / 2; 
   
   interval pair;

   if ( key < table[ midIndex ] ) {
	          
			   pair.start = indexInterval.start;
			   pair.end = midIndex;

	           return binarySearch ( table, pair, key );
    }
    else {
                 pair.start = midIndex;
			     pair.end = indexInterval.end;
		         return binarySearch( table, pair, key );
	}
		 
}

double computeScatCrossSection( double radius, double lambda) {
	
    interval radius_interval, spectrum_interval;
	double scatCrossSection;

	// check if key falls within the min and max values of the table. Otherwise, return false;
	
    if ( ! ( radius  >=  radii[ 0] && radius <= radii[ nRadii-1] )  ) {
		messageFile << "computeScatCrossSection: invalid radius: " << radius << std::endl;
		exit( EXIT_FAILURE);
	}
	
	if ( ! ( lambda  >=  spectralSamples[ 0] && lambda <= spectralSamples[ nLambdas-1] )  ) {
		messageFile << "computeScatCrossSection: invalid wavelength" << lambda << std::endl;
		exit( EXIT_FAILURE);
	}
	
	
	/*
	radius_interval.start =0;
	radius_interval.end = nRadii  - 1;

	spectrum_interval.start = 0;
	spectrum_interval.end = nLambdas -1;
	
    radius_interval = binarySearch( radii, radius_interval,  radius);
	spectrum_interval = binarySearch( spectralSamples, spectrum_interval, lambda );
	
	// check if key falls within the min and max values of the table. Otherwise, return false;
	
  

	// trilinear interpolation with respect to radius, wavelength, and theta
    // bilinear interpolation: f (x,y) = wx * wy * f(x0, y0) + wx * (1 - wy) * f(x0, y1) +
	//                                  (1-wx) * wy * f(x1, x0) + (1-wx)(1-wy)f(x1,y1)
	// trilinear interpolation: f(x,y,z) = wx wy wz f(x0,y0,z0) + wxwy(1-wz)f(x0,y0,z1) +
	//                                     wx(1-wy)wz f(x0,y1,z0) + wx(1-wy)(1-wz) f(x0,y1,z1) +
	//                                     (1-wx)wywzf(x1,y0,z0) + (1-wx)wy(1-wz)f(x1,y0,z1) +
	//                                     (1-wx)(1-wy)wz f(x1,y1,z0) + (1-wx)(1-wy)(1-wz)f(x1,y1,z1) 
    
	int x0 = radius_interval.start;
	int x1 = radius_interval.end;
	int y0 = spectrum_interval.start;
	int y1 = spectrum_interval.end;
	*/

	int i = int ( ( radius - radiusStart ) / radiusStep );	  
	int j = int ( (lambda - lambdaStart ) /  lambdaStep );

	double wx = ( radius - radii[ i ] ) / ( radii[i+1] - radii[i] );
	double wy = ( lambda - spectralSamples[ j] ) / ( spectralSamples[ j+1] - spectralSamples[j] );
	
	scatCrossSection = wx * wy*  FCscat[i][j] + wx*  (1-wy) *  FCscat[i][j+1] + 
		(1-wx) * wy *  FCscat[i+1][j] +  (1- wx) * (1 -wy) *  FCscat[i+1][j+1];
	


    return scatCrossSection;
}


double computePhaseFunction( double radius, double lambda, double theta) {
	
  
    double  P = 0.0;
	


	// check if key falls within the min and max values of the table. Otherwise, return false;
	
    if ( ! ( radius  >=  radii[ 0] && radius <= radii[ nRadii-1] )  ){
		messageFile << "computePhase: invalid radius" << radius << std::endl;
		exit( EXIT_FAILURE);
	}  
	
	if ( ! ( lambda  >=  spectralSamples[ 0] && lambda <= spectralSamples[ nLambdas-1] )  ) {
		messageFile << "computePhase: invalid wavelength" << lambda << std::endl;
		exit( EXIT_FAILURE);
	}
	if ( ! ( theta  >=  thetas[ 0] && theta <= thetas[ nThetas-1] )  ) {
		messageFile << "compputePhase: invalid scattering angle" << theta << std::endl;
		exit( EXIT_FAILURE);
	}

	int i = int ( ( radius - radiusStart ) / radiusStep );	  
	int j = int ( (lambda - lambdaStart ) /  lambdaStep );
	int k = int ( ( theta - thetaStart ) /  thetaStep );

	 
	/*
	interval radius_interval, spectrum_interval, theta_interval;
	
	radius_interval.start =0;
	radius_interval.end = nRadii  - 1;

	spectrum_interval.start = 0;
	spectrum_interval.end = nLambdas -1;
	theta_interval.start =0;
	theta_interval.end = nThetas -1;

    radius_interval = binarySearch( radii, radius_interval,  radius);
	spectrum_interval = binarySearch( spectralSamples, spectrum_interval, lambda );
	theta_interval =    binarySearch( thetas, theta_interval, psi);

	// check if key falls within the min and max values of the table. Otherwise, return false;
	
  

	// trilinear interpolation with respect to radius, wavelength, and theta
    // bilinear interpolation: f (x,y) = wx * wy * f(x0, y0) + wx * (1 - wy) * f(x0, y1) +
	//                                  (1-wx) * wy * f(x1, x0) + (1-wx)(1-wy)f(x1,y1)
	// trilinear interpolation: f(x,y,z) = wx wy wz f(x0,y0,z0) + wxwy(1-wz)f(x0,y0,z1) +
	//                                     wx(1-wy)wz f(x0,y1,z0) + wx(1-wy)(1-wz) f(x0,y1,z1) +
	//                                     (1-wx)wywzf(x1,y0,z0) + (1-wx)wy(1-wz)f(x1,y0,z1) +
	//                                     (1-wx)(1-wy)wz f(x1,y1,z0) + (1-wx)(1-wy)(1-wz)f(x1,y1,z1) 
    
	int x0 = radius_interval.start;
	int x1 = radius_interval.end;
	int y0 = spectrum_interval.start;
	int y1 = spectrum_interval.end;
	int z0 = theta_interval.start;
	int z1 = theta_interval.end;

	*/


	double wx = ( radius - radii[ i ] ) / ( radii[i+1] - radii[i] );
	double wy = ( lambda - spectralSamples[ j] ) / ( spectralSamples[ j+1] - spectralSamples[j] );
	double wz = ( theta - thetas[ k ] ) / ( thetas[k+1] - thetas[k] );

	P = wx * wy* wz* particlePhase[i][j][k] + wx* wy* (1-wz) * particlePhase[i][j][k+1] + 
		wx* (1-wy) * wz * particlePhase[i][j+1][k] + wx * (1- wy) * (1 -wz) * particlePhase[i][j+1][k+1] +
		(1-wx) * wy * wz * particlePhase[i+1][j][k] + (1-wx) * wy * (1-wz) * particlePhase[i+1][j][k+1] + 
		(1-wx) * (1-wy) * wz * particlePhase[i+1][j+1][k] + (1-wx) * (1-wy) * (1 - wz) * particlePhase[i+1][j+1][k+1];
	


    return P;
}


void get_lambda_and_psi_for_P (double &l, double &p) {
    
    std::string input;
	// typedef basic_string<char, char_traits<char>, allocator<char> > 	string;


    // get numbers
    while (true) {
        std::cout << "Please enter lambda (0.00000040 ~ 0.00000069): ";
        getline(std::cin, input);
        
        // This code converts from string to number safely.
        std::stringstream myStream(input);
        if (myStream >> l)
            break;
        std::cout << "Invalid number, please try again" << std::endl;
    }
    std::cout << "You entered: " << l << std::endl;
    
    while (true) {
        std::cout << "Please enter psi (0 ~ 180): ";
        getline(std::cin, input);
        
        // This code converts from string to number safely.
        std::stringstream myStream(input);
        if (myStream >> p)
            break;
        std::cout << "Invalid number, please try again" << std::endl;
    }
    std::cout << "You entered: " << p << std::endl;
    
}

