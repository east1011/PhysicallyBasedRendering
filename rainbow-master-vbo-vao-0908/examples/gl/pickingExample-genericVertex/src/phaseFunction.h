
#ifndef Solar_Position_Algorithm_phaseFunction_h
#define Solar_Position_Algorithm_phaseFunction_h




const int nLambdasteps = 120;
const int nRadiusSteps = 2;
const int nThetaSteps = 100; // 1 degree in each step between theta 130 deg and theta 142;

const int nLambdas = nLambdasteps +1;
const int nRadii = nRadiusSteps +1;
const int nThetas = nThetaSteps +1;

const double lambdaStart = 400;   // 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
const double lambdaEnd = 700;
     
const  double radiusStart =  0.05e-3;  // 0.05 mm = 0.05 * e-3 m;
const  double radiusEnd = 2e-3;           // 2 mm = 2 * e-3 m

const	double thetaStart = 130.0;
const	double thetaEnd = 142.0;

// the step size of each sample; each sample is considered the center point of the step
const double lambdaStep = ( lambdaEnd - lambdaStart) / double( nLambdas );
const double thetaStep = ( thetaEnd - thetaStart) / double (nThetas);
const double radiusStep = (radiusEnd - radiusStart) / double (nRadii);


extern  double Cabs[ nRadii][nLambdas];
extern double FCscat[nRadii][nLambdas];


extern double particlePhase [nRadii][nLambdas][nThetas];
extern double  radii[ nRadii ];
extern double  spectralSamples [ nLambdas ];
extern double  thetas[ nThetas ];

// externals whose contents are defined in dropvolume.cpp

/*
extern double Cabs[ nRadii ] [nLambdas];
extern double FCscat[nRadii][nLambdas];


extern double particlePhase [nRadii] [nLambdas] [nThetas];
extern double  radii[ nRadii ];
extern double  spectralSamples [ nLambdas ];
extern double  thetas[ nThetas ];

*/

typedef struct {
	int start;
	int end;
} interval;


					   		
double computeScatCrossSection( double radius, double lambda);

double computePhaseFunction( double radius, double lambda, double psi);

interval  binarySearch(double *table, interval pair, double key);
//void get_lambda_and_psi_for_P (double &l, double &p);

#endif