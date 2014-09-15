//
//  irradiance_sun.h
//  Solar Position Algorithm
//
//  Created by jd on 2014. 2. 1..
//  Copyright (c) 2014 jd. All rights reserved.
//

#ifndef Solar_Position_Algorithm_irradiance_sun_h
#define Solar_Position_Algorithm_irradiance_sun_h

#include <cmath>

#include "irradiance_sun.h"
#include "cvec.h"

double calculate_irradiance_of_sun (double lambda) {
    
    double R_Earth = 1.496e+8;
    double R_Sun = 6.95e+5;
    double h = 6.6261e-34;      // Planck's constant
    double c = 2.9979e+8;       // speed of light in vacuo
    double k = 1.3806e-23;      // Boltzmann's constant
    double T = 5782;            // Sun's temperature in Kelvin
    double Isun;
    
    Isun = pow((R_Sun/R_Earth),2) * 2*PI*h*pow(c,2)  / ( pow(lambda, 5) * ( exp( h*c/(lambda*k*T) ) -1 ) ) * 1.0e-9  ;
	   //Moon:  1.0e-9 is multiplied to convert the irradiance per meter to the irradiance per nanometer
	   // Original had 1.0e-10. 
    
    return Isun;
}

#endif
