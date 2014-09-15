//
//  irradiance_sun.h
//  Solar Position Algorithm
//
//  Created by jd on 2014. 2. 1..
//  Copyright (c) 2014 jd. All rights reserved.
//

#ifndef irradiance_sun_h
#define irradiance_sun_h

static double pi = acos( -1.0e0);  // -1.0 const double in default; -1.0f = float

#define PI 3.1415926535897932384626433832795028841971

double calculate_irradiance_of_sun (double lambda);

#endif
