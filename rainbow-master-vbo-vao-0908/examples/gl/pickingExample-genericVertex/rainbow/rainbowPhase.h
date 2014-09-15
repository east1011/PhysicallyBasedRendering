// rainbow/rainbowPhase.h

#ifndef PBRT_RAINBOW_RAINBOWPHASE_H
#define PBRT_RAINBOW_RAINBOWPHASE_H

#include "spectrum.h" // this includes "spectrum.h" which is needed in rainbowPhase.h

Spectrum particle_phasefunction( double radius, double theta, const Spectrum FCscatSpec);
double lambda_phasefunction(double lambda, double radius, double theta, double FCscat);



#endif // PBRT_RAINBOW_RAINBOWPHASE_H