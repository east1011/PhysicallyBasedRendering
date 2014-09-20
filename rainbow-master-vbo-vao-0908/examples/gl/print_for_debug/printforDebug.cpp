//#include "cvec.h"
//#include "spa.h"
//#include "matrix4.h"
//
//#define IlluminantCxWhite   0.3101	    /* For NTSC television */
//#define IlluminantCyWhite	0.3162
//#define IlluminantD65xWhite 0.3127		/* For EBU and SMPTE */
//#define IlluminantD65yWhite 0.3291	
//#define IlluminantExWhite 	0.33333333	/* CIE equal-energy illuminant */
//#define IlluminantEyWhite	0.33333333 
//#define GAMMA_REC709		0	        /* Rec. 709 */
//
//#define PI 3.1415926535897932384626433832795028841971
//
//int g_windowWidth = 1024;
//int g_windowHeight = 768;
//float g_dropDensity = 30000;
//
//Cvec3 g_sunRayDir;
//Matrix4 g_cameraRotMat;
//Matrix4 g_eyeRbt;
//Matrix4 g_invEyeRbt;
//Matrix4 g_SceneRbt;
//Cvec4 g_light1Pos;
//float  g_radius;
//Cvec3 g_AABBsize;
//Matrix4 g_AABBRbt;
//
//static const float rainbowVolumeHeight = 40; 
//static const float rainbowVolumeWidth =100;
//static const float rainbowVolumeDepth = 10;
//
//
//struct colourSystem {
//	//uint name[15];     	    /* Colour system name */
//	float xRed, yRed,	    	/* Red x, y */
//		  xGreen, yGreen,  	    /* Green x, y */
//		  xBlue, yBlue,    	    /* Blue x, y */
//		  xWhite, yWhite,  	    /* White point x, y */
//		  gamma;   	    		/* Gamma correction for system */
//};
//
//const int nSpectralSamples = 61;
//const float lambdaStart = 400;		// 400 nm = 400 * e-9 m = 0.4 * e-6 m = 0.4 um
//const float lambdaEnd = 700;
//const float radiusStart =  1.0e-3;  
//const float radiusEnd = 2.0e-3;		// 2 mm = 2 * e-3 m
//const float thetaStart = 130.0;
//const float thetaEnd = 142.0;
//float pi = acos(-1.0);  // -1.0 const float in default; -1.0f = float
//
//const float lambdaStep = ( lambdaEnd - lambdaStart) / float( nSpectralSamples -1); //lambdaStep = (700-400)/(61-1)
//
//float ISpectrumSun[nSpectralSamples] = {
////wavelength (nm)	normalized irradiance (arbitrary units)
//// 700 - 400 = 300 / 5 = 60 + 1 = nSpectralSamples= 61 
//// The original wavelength starts from 380 to 700 with step 5 nm
////380
////	0.4651826, 	0.3830807, 	0.4409676, 	0.4415811, 	
////400
//	0.5531539, 	0.7498946, 0.746196, 0.7945927,
////420
//	0.8004899, 	0.7938542, 	0.7373845, 	0.7613856,   0.8335401,  0.9062264,   0.9455949, 0.9749124,
////460
//	0.9885056, 	0.9996937, 	0.9703307, 	0.9874594, 	1, 	0.9923742, 	0.9242011, 	0.9640443, 	0.9617534,
////505
//	0.9306258, 	0.9463357, 	0.9275853, 	0.8597925, 	0.9149771, 	0.9065616, 	0.9197563, 	0.9100043,
////545
//	0.8943484, 0.8951805, 0.8915439, 0.8640333,  	0.8611063,   0.8411431,  	0.8347283,  0.82272,
//// 585	
//   0.8433171, 	0.7949293,  	0.7905713,  0.7947434,  0.7932764,  0.7953064, 	0.769975, 	0.7614744,
////625
//	0.7494112,  0.7195148,  	0.7222514,  0.7319365,  0.7231518,  0.6976414,  0.6917001,  0.6692451,
////665
//	0.7032121,  0.6983164,  	0.6871358,  0.6830449, 	0.6650319, 	0.5508773, 	0.5994822, 	0.6185095
//};
//
//const int nLambdasForColorMatch = 81;
//const int nRGB = 3;
//float cie_colour_match[ nLambdasForColorMatch * nRGB ] = { 
//	    0.0014,0.0000,0.0065, 0.0022,0.0001,0.0105, 0.0042,0.0001,0.0201,
//        0.0076,0.0002,0.0362, 0.0143,0.0004,0.0679, 0.0232,0.0006,0.1102,
//        0.0435,0.0012,0.2074, 0.0776,0.0022,0.3713, 0.1344,0.0040,0.6456,
//        0.2148,0.0073,1.0391, 0.2839,0.0116,1.3856, 0.3285,0.0168,1.6230,
//        0.3483,0.0230,1.7471, 0.3481,0.0298,1.7826, 0.3362,0.0380,1.7721,
//        0.3187,0.0480,1.7441, 0.2908,0.0600,1.6692, 0.2511,0.0739,1.5281,
//        0.1954,0.0910,1.2876, 0.1421,0.1126,1.0419, 0.0956,0.1390,0.8130,
//        0.0580,0.1693,0.6162, 0.0320,0.2080,0.4652, 0.0147,0.2586,0.3533,
//        0.0049,0.3230,0.2720, 0.0024,0.4073,0.2123, 0.0093,0.5030,0.1582,
//        0.0291,0.6082,0.1117, 0.0633,0.7100,0.0782, 0.1096,0.7932,0.0573,
//        0.1655,0.8620,0.0422, 0.2257,0.9149,0.0298, 0.2904,0.9540,0.0203,
//        0.3597,0.9803,0.0134, 0.4334,0.9950,0.0087, 0.5121,1.0000,0.0057,
//        0.5945,0.9950,0.0039, 0.6784,0.9786,0.0027, 0.7621,0.9520,0.0021,
//        0.8425,0.9154,0.0018, 0.9163,0.8700,0.0017, 0.9786,0.8163,0.0014,
//        1.0263,0.7570,0.0011, 1.0567,0.6949,0.0010, 1.0622,0.6310,0.0008,
//        1.0456,0.5668,0.0006, 1.0026,0.5030,0.0003, 0.9384,0.4412,0.0002,
//        0.8544,0.3810,0.0002, 0.7514,0.3210,0.0001, 0.6424,0.2650,0.0000,
//        0.5419,0.2170,0.0000, 0.4479,0.1750,0.0000, 0.3608,0.1382,0.0000,
//        0.2835,0.1070,0.0000, 0.2187,0.0816,0.0000, 0.1649,0.0610,0.0000,
//        0.1212,0.0446,0.0000, 0.0874,0.0320,0.0000, 0.0636,0.0232,0.0000,
//        0.0468,0.0170,0.0000, 0.0329,0.0119,0.0000, 0.0227,0.0082,0.0000,
//        0.0158,0.0057,0.0000, 0.0114,0.0041,0.0000, 0.0081,0.0029,0.0000,
//        0.0058,0.0021,0.0000, 0.0041,0.0015,0.0000, 0.0029,0.0010,0.0000,
//        0.0020,0.0007,0.0000, 0.0014,0.0005,0.0000, 0.0010,0.0004,0.0000,
//        0.0007,0.0002,0.0000, 0.0005,0.0002,0.0000, 0.0003,0.0001,0.0000,
//        0.0002,0.0001,0.0000, 0.0002,0.0001,0.0000, 0.0001,0.0000,0.0000,
//		0.0001,0.0000,0.0000, 0.0001,0.0000,0.0000, 0.0000,0.0000,0.0000};
//
//colourSystem g_cs = {0.67,  0.33,  0.21,  0.71,  0.14, 0.08, IlluminantCxWhite, IlluminantCyWhite, GAMMA_REC709};
//
//float LSpectrumSun_Lee[nSpectralSamples];
//float LSpectrumSun_Plank[nSpectralSamples];
//
//void calculate_radiance_of_sun_Lee();
//void calculate_radiance_of_sun_Plank();
//float calculate_radiance_of_sun(float lambda);
//
//Cvec3  xyz_to_rgb( colourSystem cs, Cvec3 xyz);
//bool constrain_rgb(Cvec3 rgb);
//void norm_rgb( Cvec3 rgb);
//Cvec3 oldToRGB( float irainbow[] );
//float max(float LSpectrumSun_Plank[]);
//void normalize_radiance_of_sun_Plank (float LSpectrumSun_Plank[], float L_Plank_max);
//
//
//
//
//void calculate_radiance_of_sun_Lee() {
//
//    float R_Earth = 1.496e+8;
//    float R_Sun = 6.95e+5;
//    float h = 6.6261e-34;      // Planck's constant
//    float c = 2.9979e+8;       // speed of light in vacuo
//    float k = 1.3806e-23;      // Boltzmann's constant
//    float T = 5782;            // Sun's temperature in Kelvin
//    float Lsun;
//
//	float lambda = lambdaStart;
//	
//    for ( int j = 0; j < nSpectralSamples;  j++) {
//		
//		//LSpectrumSun_Lee[ j ] = ISpectrumSun[j] / pi;		
//		LSpectrumSun_Lee[ j ] = ISpectrumSun[j];	
//
//	}  // for 
//
//} //calculate_radiance_of_sun_Lee()
//
//
//void calculate_radiance_of_sun_Plank() {
//
//    float R_Earth = 1.496e+8;
//    float R_Sun = 6.95e+5;
//    float h = 6.6261e-34;      // Planck's constant
//    float c = 2.9979e+8;       // speed of light in vacuo
//    float k = 1.3806e-23;      // Boltzmann's constant
//    float T = 5782;            // Sun's temperature in Kelvin
//    float Lsun;
//
//	float lambda = lambdaStart;
//	
//    for ( int j = 0; j < nSpectralSamples;  j++) {
//		
//		Lsun = calculate_radiance_of_sun(lambda);
//
//		LSpectrumSun_Plank[ j ] = Lsun;
//		lambda += lambdaStep;
//		
//	}  // for 
//    
//} //calculate_radiance_of_sun_Plank()
//
//
//float calculate_radiance_of_sun(float lambda) {
//    // http://www.oceanopticsbook.info/view/light_and_radiometry/level_2/blackbody_radiation
//	// Radiation in thermodynamic equilibrium is isotropic and unpolarized ( eq. (4) )
//    float R_Earth = 1.496e+11;	//The radious of the sphere at the center of the sun passing through the earth
//    float R_Sun = 6.95e+8;		//the radious of the sun's surface
//    float h = 6.6261e-34;      // Planck's constant
//    float c = 2.9979e+8;       // speed of light in vacuo
//    float k = 1.3806e-23;      // Boltzmann's constant
//    float T = 5782;            // Sun's temperature in Kelvin
//    float Isun;
//
//			
//	float lambda_m = lambda * 1.0e-9;				// convert lambda to meters from nanometers
//	// blackbody radiation spectrum [W/m^2/nm]
//	//Isun = pow( (R_Sun/R_Earth), (float)2.0 ) * 2 * pi * h * pow(c,(float)2.0)/( pow(lambda_m, (float)5.0) * ( exp(  h*c/ ( lambda_m*k*T ) ) -1 ) ) * 1.0e-9  ;
//	//    //Moon:  1.0e-9 is multiplied to convert the irradiance per meter to the irradiance per nanometer
//	//    // Original had 1.0e-10.
//	//float Lsun = Isun / pi; //convert the irradiance to radiance 
//
//	float Lsun =  (R_Sun/R_Earth)*(R_Sun/R_Earth) * 2 * h * (c * c)/
//					( pow(lambda_m, (float)5.0) * ( exp( h*c/(lambda_m*k*T) ) -1 ) ) * 1.0e-9  ;
//
//	return Lsun;
//		
//} //calculate_radiance_of_sun()
//
//Cvec3  xyz_to_rgb( colourSystem cs, Cvec3 xyz) {
//
//    float xr, yr, zr, xg, yg, zg, xb, yb, zb;
//    float xw, yw, zw;
//    float rx, ry, rz, gx, gy, gz, bx, by, bz;
//    float rw, gw, bw;
//
//	xr = cs.xRed;    yr = cs.yRed;    zr = 1 - (xr + yr);
//    xg = cs.xGreen;  yg = cs.yGreen;  zg = 1 - (xg + yg);
//    xb = cs.xBlue;   yb = cs.yBlue;   zb = 1 - (xb + yb);
//
//    xw = cs.xWhite;  yw = cs.yWhite;  zw = 1 - (xw + yw);
//
//    // xyz -> rgb matrix, before scaling to white. 
//    
//    rx = (yg * zb) - (yb * zg);  ry = (xb * zg) - (xg * zb);  rz = (xg * yb) - (xb * yg);
//    gx = (yb * zr) - (yr * zb);  gy = (xr * zb) - (xb * zr);  gz = (xb * yr) - (xr * yb);
//    bx = (yr * zg) - (yg * zr);  by = (xg * zr) - (xr * zg);  bz = (xr * yg) - (xg * yr);
//
//    // White scaling factors.
//    //  Dividing by yw scales the white luminance to unity, as conventional. 
//       
//    rw = ((rx * xw) + (ry * yw) + (rz * zw)) / yw;
//    gw = ((gx * xw) + (gy * yw) + (gz * zw)) / yw;
//    bw = ((bx * xw) + (by * yw) + (bz * zw)) / yw;
//
//    // xyz -> rgb matrix, correctly scaled to white. 
//    
//    rx = rx / rw;  ry = ry / rw;  rz = rz / rw;
//    gx = gx / gw;  gy = gy / gw;  gz = gz / gw;
//    bx = bx / bw;  by = by / bw;  bz = bz / bw;
//
//    // rgb of the desired point 
//    
//	Cvec3 rgbColor = Cvec3(	(rx * xyz[0]) + (ry * xyz[1]) + (rz * xyz[2]), 
//							(gx * xyz[0]) + (gy * xyz[1]) + (gz * xyz[2]), 
//							(bx * xyz[0]) + (by * xyz[1]) + (bz * xyz[2]) );
//
//	return rgbColor;
//
//} // xyz_to_rgb
//
//bool constrain_rgb(Cvec3 rgb) {
//    float w;
//
//    // Amount of white needed is w = - min(0, *r, *g, *b) 
//	    
//	if (0.0 < rgb[0]) w = 0.0; else w = rgb[0];
//	if  (w < rgb[1])  w =w; else w =  rgb[1];
//	if  (w < rgb[2])  w =w; else w = rgb[2];
//    w = -w;
//
//    // Add just enough white to make r, g, b all positive. 
//    
//    if (w > 0.0) {
//		rgb += w;
//        return true;                     // Colour modified to fit RGB gamut 
//    }
//    return false;                         // Colour within RGB gamut 
//}
//
//void norm_rgb( Cvec3 rgb) {
//	#define Max(a, b)   (((a) > (b)) ? (a) : (b))
//	float greatest = Max(rgb[0], Max(rgb[1], rgb[2]));
//    
//    if (greatest > 0) {
//		rgb[0] /= greatest;
//		rgb[1] /= greatest;
//		rgb[2] /= greatest;
//    }
//	#undef Max
//}
//
//Cvec3 oldToRGB( float irainbow[] ) {
//
//	// convert spectral color to rgb color
//	
//	float X = 0;
//	float Y = 0;
//	float Z = 0;
//	float XYZ;
//
//	float xBar, yBar, zBar;
//	
//	for (int j= 0; j < nSpectralSamples; j++) {
//		
//		int k = j + 4;
//
//		xBar =  cie_colour_match[ k * 3 + 0 ]; 
//		yBar =  cie_colour_match[ k * 3 + 1 ];
//		zBar =  cie_colour_match[ k * 3 + 2 ];
//
//		X += irainbow[j] * xBar * lambdaStep; 
//		Y += irainbow[j] * yBar * lambdaStep;
//		Z += irainbow[j] * zBar * lambdaStep;
//		
//    }
//	
//	XYZ = (X+Y+Z);
//	Cvec3 xyzColor = Cvec3( X/XYZ, Y/XYZ, Z/XYZ );
//
//	Cvec3 rgbColor = xyz_to_rgb(g_cs, xyzColor);
//
//	if (constrain_rgb(rgbColor)) {  // Colour modified to fit RGB gamut 
//		norm_rgb(rgbColor);
//		//output << rgbColor << "(Approximation)" << endl;
//	} else {
//		norm_rgb(rgbColor);
//		//output << rgbColor << endl;
//	}
//	//constrain_rgb(rgbColor);
//
//	return rgbColor;
//
//} // oldToRGB()
//
//float max(float LSpectrumSun_Plank[]) {
//
//	float Plank_max = 0;
//
//	for (int i = 0; i < nSpectralSamples; i++) {
//		if (Plank_max < LSpectrumSun_Plank[i]) Plank_max = LSpectrumSun_Plank[i];
//	}
//	return Plank_max;
//}
//
//void normalize_radiance_of_sun_Plank (float LSpectrumSun_Plank[], float L_Plank_max) {
//
//	for (int i = 0; i < nSpectralSamples; i++) {
//		LSpectrumSun_Plank[i] = LSpectrumSun_Plank[i]/L_Plank_max;
//	}
//}
//
//int* test3()
//{
//    static int array[10];
//    return array;
//}
//
////float* LSpectrumIn ( float LSpectrumSun[], Cvec3 currPos, Cvec3 dirToLight, float sigma_e) {
////
////	// attenuate the sun light Lsun0 from the entering point to AABB to Pos
////
////	float tmin, tmax;
////	float LSpectIn[nSpectralSamples];
////
////	Matrix4 uInvModelMatrixAABB = inv( g_AABBRbt );
////
////	Cvec3 lightDirInBoxFrame = Cvec3( uInvModelMatrixAABB * g_eyeRbt * Cvec4( dirToLight, 0) );
////	Cvec3 lightOriginInBoxFrame = Cvec3( uInvModelMatrixAABB * g_eyeRbt * Cvec4( currPos,1) );
////
////	tmin = 0.0; tmax = 0.0;
////
////	bool isHit = rayOBBIntersect ( lightOriginInBoxFrame, lightDirInBoxFrame, 
////									uAABBmin, uAABBmax, tmin, tmax);
////
////	static float array[10];
////    return array;
////	//return LSpectrumSun;
////  
////} // LSpectrumIn()
//
//float computePhaseFunction( float radius, float lambda, float theta) {
//  	
//    float z0 =   (radius  - radiusStart) / (radiusEnd - radiusStart );
//	float y0 =  (lambda - lambdaStart) / (lambdaEnd - lambdaStart );
//	float x0 =  (theta - thetaStart ) / ( thetaEnd - thetaStart );
//	
//	float phase = texture(uPhaseTex, Cvec3(x0,y0,z0) ).r;
//	
//    return phase;
//
//}
//
//float avgPhaseFunction(float radius, float lambda, float psi) {
//
//    float particlePhase  = computePhaseFunction( radius,  lambda, psi); 
//
//	//float particlePhase = 1.230634e-002;
//
//    float avgPhasePerVolume  =  g_dropDensity * particlePhase;
//    return avgPhasePerVolume;
//
//}
//
//Cvec3 singleScatteringAndAttenuation(float radius, Cvec3 surfaceColor, float zEye, 
//                              bool isPointLight, Cvec3 lightPosOrDir, float LSpectrumSun[],
//	                          Cvec3 rayDir, float sigma_s, float sigma_e )  {
//
//	float tmin, tmax;
//	float LSpectrumOut[nSpectralSamples], LSpectrumOut0[nSpectralSamples];
//
//
//	
//	g_SceneRbt = Matrix4::makeAxisRotation( - PI / 3.0, Cvec3(0,1,0) );
//	g_AABBsize = (rainbowVolumeWidth, rainbowVolumeHeight, rainbowVolumeDepth);
//	g_AABBRbt = g_SceneRbt * Matrix4::makeTranslation( 	Cvec3( 0, g_AABBsize[1]/2 + 6.4, -g_AABBsize[2]/2.0 ) ) ) );
//	Matrix4 uInvModelMatrixAABB = inv( g_AABBRbt );
//	
//	Cvec3 rayDirInBoxFrame = Cvec3( uInvModelMatrixAABB * g_eyeRbt * Cvec4(rayDir, 0) );
//
//	//bool isHit = rayOBBIntersect ( uEyeOriginInBoxFrame, rayDirInBoxFrame, uAABBmin, uAABBmax, tmin, tmax );
//	bool isHit = false;
//	if ( !isHit ) { // the AABB is not hit? 	 getBackgroundColor:
//  		 
//		return surfaceColor;
//	  	  					
//	} // if (!isHit) // NOT HIT
//
//	else { //  hit the front at tmin and the back face of AABB box at tmax;
//	  
//		float tmaxAtten = tmax; // the t value at which attenuation begins
//
//		float distFromEye = -zEye;
//
//		if ( distFromEye < tmin ) { // the surface is before the AABB box => return the surfaceColor
//			
//			return surfaceColor;
//		
//		}
//
//		if ( distFromEye >= tmin && distFromEye <= tmax ) {
//			tmaxAtten = distFromEye;
//		}
//		else { // the surface is behind the AABB box
//			tmaxAtten = tmax; 
//		}
//	  	    
//		int N = 30;
//		float deltaT = (tmaxAtten - tmin) / float(N);
//   
//   
//		// get the light which has been scattered in the eye direction by each region in the volume 
//		// starting from tmaxAtten.
//
//		float LSpectIn[nSpectralSamples];
//
//		for (int j = 0; j < nSpectralSamples; j++ ) {
//			LSpectrumOut[j] = 0.0;
//		}
//
//
//
//		// LIGHT DIRECTIONS
//		if ( isPointLight ) { // point light
//      
//			Cvec3 lightPos = lightPosOrDir;
//
//			for (float t = tmaxAtten; t >= tmin; t -= deltaT) {  // for each position on the ray
//      
//				Cvec3 currPos = Cvec3(0,0,0) + rayDir * t;
//
//				Cvec3 dirFromPointLight = currPos - lightPos;
//
//				dirFromPointLight = normalize( dirFromPointLight );
//
//				LSpectIn = LSpectrumIn( LSpectrumSun, currPos, -dirFromPointLight, sigma_e);
//		
//				// The light with given lambda is scattered in the direction of -rayDir with varying amount
//				// depending on the scattering angle between the light direction and -rayDir, which varies
//				// with each t. Consider only the scattering angles  within the range [ thetaStart, thetaEnd]. 
//				// The other scattering angles constributs less to the scattered intensity of lambda
//				// They are ignored for computational reasons. Some scattering spectrum is computed for
//				// every ray direction (-rayDir)?? 
//	  
//				// compute the scattered light scattered  along the ray
//				float psi = acos ( dot ( -rayDir, dirFromPointLight) ) * 180.0 / pi;
//		
//				if (  psi >=  thetaStart  &&  psi <= thetaEnd  ) {	// only the light [at the current position on the ray]
//																	// whose scattering angle is  within these angles
//																	// contribute to the rainbow color along the current ray.
//																	// Otherwise, no color is contributed.  
//					float lambda = lambdaStart;
//					int   lambdaIndex = 0;	
//					for ( int j = 0; j < nSpectralSamples;  j++) {
//		
//						LSpectrumOut[j] += exp( - sigma_e * ( t - tmin) )
//										* sigma_s * avgPhaseFunction(radius, lambda, psi)
//										* LSpectIn[j];
//						lambda += lambdaStep;
//			   
//			   
//					} // for
//		   	   
//				}   //if 
//
//				else continue; // go to the next position
//
//			} // for all positions on the current ray
//
//			return oldToRGB( LSpectrumOut );
//			//return toRGB( LSpectrumOut );
//
//		} // if (point light)
//		
//		else { // directional light
//
//			Cvec3 dirFromSun  = lightPosOrDir;
//			float psi = acos ( dot ( -rayDir, dirFromSun ) ) * 180.0 / pi;
//
//			// for debugging
//			// spect1 = vec4(dirFromSun, 1.0);
//
//			//spect1 = vec4(dot ( -rayDir, dirFromSun ), acos(dot ( -rayDir, dirFromSun ) ), psi, 1.0 );
//			//return LSpectrumOut;
//
//			// consider only the cases where the light is scattered in the direction of the eye
//			// with scattering angle within [ thetaStart, thetaEnd]. These scattering angles contricute
//			// to the rainbow and its neighborhood. The cases with other scattering
//			// angles contribute to the further region outside of the rainbow and its neighbor.
//			// These regions are not considered and the background color are assumed to dominate. 
//
//
//			// SUM_{i=0, i=N-1} exp( - sigma_e ( t_i - tmin) ) sigma_s  f_p( dot(-rayDir, dirFromLight(t_i) ) ) L(t_i, dirFromLight(t_i) )
//			// L(t_i, dirFromLight(t_i) ) = LSpectrumIn is given by
//			// tau(tmin_sun, t_i) Lsun( tmin_sun, dirFromLight(tmin_sun) ) 
//			// = exp( -sigma_e ( t_i - tmin_sun) )  Lsun( tmin_sun, dirFromLight(tmin_sun) ) 
//			// where tmin_sun is the point at which the sun light enters the volume to reach t_i.
//	  
//			if ( !( psi >= thetaStart  &&  psi <= thetaEnd ) ) {
//				
//				return surfaceColor; // just return the background color
//		 
//			}  
//
//			else {	// psi is a rainbow scattering angle
//				// LSpectrumOut = LSpectrumOut + SUM_{t = tmin, tmaxAtten} exp( - sigma_e * ( t - tmin) )  * sigma_s 
//				//                                          * avgPhaseFunction( dot ( -rayDir, dirFromSun(t) ) ) 
//				//                                          * LSpectIn(t, dirFromLight(t) )
//
//		
//
//				//return LSpectrumOut;
//				float LSpectIn[nSpectralSamples];
//
//				float incrementalSpect[nSpectralSamples];
//				
//				for (float t = tmaxAtten; t >= tmin; t-= deltaT  ) { 
//
//					//LSpectrumOut = LSpectrumOut +  exp( - sigma_e * ( t - tmin) )  * sigma_s 
//					//									 * avgPhaseFunction( dot ( -rayDir, dirFromSun(t) ) * LSpectrumIn
//				  
//					Cvec3 currPos = Cvec3(0,0,0) +  rayDir * t;
//		 		     
//					LSpectIn = LSpectrumIn(LSpectrumSun, currPos, -dirFromSun, sigma_e);
//		
//
//					float lambda = lambdaStart;
//		    		   		   
//					for ( int j =0; j < nSpectralSamples; j++) {			  
//				 
//						incrementalSpect[j] = exp( - sigma_e * ( t - tmin) ) * sigma_s  
//													* avgPhaseFunction(radius, lambda, psi) * LSpectIn[j];
//
//				 
//						//incrementalSpect[j] = exp( - sigma_e * ( t - tmin) ) * sigma_s  
//						//                          * avgPhaseFunction(radius, lambda, psi) * LSpectrumSun[j];
//						//incrementalSpect[j] = exp( - sigma_e * ( t - tmin) ) * sigma_s  
//						//                          * avgPhaseFunction(radius, lambda, psi); 
//						//LSpectrumOut[j] += exp( - sigma_e * ( t - tmin) ) * sigma_s  
//						//							* avgPhaseFunction(radius, lambda, psi) * LSpectIn[j];
//			   
//						lambda += lambdaStep;				  
//					} // for each lambda	
//
//			
//					for ( int j =0; j < nSpectralSamples; j++) {			  
//						LSpectrumOut[j] += incrementalSpect[j];						  			  
//					} // for each lambda	
//
//
//		 				   
//				} // for t: all positions on the current ray	
//		
//				// for debugging		
//
//				return oldToRGB( LSpectrumOut ); // the light accumulated along the current ray
//				//return toRGB( LSpectrumOut );
//				
//
//
//			} // else ( psi is rainbow angle)
//
//	    
//		} // else (DIRECTIONal light) 	  
//
//
//	} // else (the ray HITs the volume )
//
//} // singleScatteringAndAttenuation
//
//
//Cvec3 calculate_rainbowColor (float radius, Cvec3 surfaceColor, float zEye,  bool isPointLight,  Cvec3 lightPosOrDir, Cvec3 rayDir, 
//							  float sigma_s, float sigma_e, float LSpectrumSun_Plank[]) {
//	// compute the scattered light, which is attenuated along the direction 
//	// -rayDir 
//	
//	Cvec3 rainbowRGB;
//		
//	//float[nSpectralSamples] LSpectrumSun = calculate_radiance_of_sun_Lee(); 
//	//float LSpectrumSun[nSpectralSamples] = calculate_radiance_of_sun_Plank(); 
//	
//	rainbowRGB = singleScatteringAndAttenuation(radius,  surfaceColor, zEye, isPointLight, 
//												lightPosOrDir, LSpectrumSun_Plank, rayDir, sigma_s, sigma_e );
//    
//	return rainbowRGB;
//	
//}//calculate_rainbowColor
//
//Matrix4 makeProjectionMatrix() {
//
//	double g_frustFovY = 100.0;
//	double g_frustNear = -5.0;    // near plane
//    double g_frustFar = -200.0;    // far plane
//
//	return Matrix4::makeProjection(
//		g_frustFovY, g_windowWidth / (double) (g_windowHeight),
//		g_frustNear, g_frustFar);
//
//}
//
//float computeScatCrossSection( float radius, float lambda ) {
//
//	// get the index to the texture for (radius, lambda)
//	Cvec2 index = Cvec2((lambda - lambdaStart) / (lambdaEnd - lambdaStart),
//						(radius - radiusStart) / (radiusEnd - radiusStart));
//	
//	float scatCrossSection = texture( uScatTex, index).r;
//		
//    return scatCrossSection;
//	
//}
//
//void main () {
//	
//	calculate_radiance_of_sun_Lee();
//	ofstream outFile0("Lee_out.txt");
//	for ( int i = 0; i < nSpectralSamples;  i++) {
//		
//		outFile0 << LSpectrumSun_Lee[i] << endl;
//		
//	}  // for
//	outFile0.close();
//	
//	float L_Plank_max = max(LSpectrumSun_Plank);
//	normalize_radiance_of_sun_Plank(LSpectrumSun_Plank, L_Plank_max);
//
//	ofstream outFile2("Plank_out_after.txt");
//	for ( int i = 0; i < nSpectralSamples;  i++) {
//		
//		outFile2 << LSpectrumSun_Plank[i] << endl;
//		
//	}  // for
//	outFile2.close();
//
//
//	spa_data spa;  //declare the SPA structure
//		int result;
//		Cvec3 sunRay;
//    
//		// Default values to calculate sun's position
//		spa.year          = 2014;
//		spa.month         = 11;
//		spa.day           = 1;
//		//spa.month         = 1;
//		//spa.day           = 31;
//		spa.hour          = 16;
//		spa.minute        = 30;
//		spa.second        = 00;
//		spa.timezone      = 9.0;             // KST (GMT + 9): Korean or Japan time, 9 hours ahead
//		//spa.delta_ut1     =0;
//		spa.delta_t       = 67;              // Doesn't it depend on location?  go to the related webpage to find the correct value.
//		spa.longitude     = 126.9410634;     // Sogang University's longitude in decimal
//		spa.latitude      = 37.5517132;      // Sogang University's latitude in decimal
//		spa.elevation     = 200;             // m
//		spa.pressure      = 820;             // mbar [milibar]
//		spa.temperature   = 23;              // celcius
//		//spa.slope         = 30;              // surface slope angle, used to compute the sun's incidence angle to the surface
//		//spa.azm_rotation  = -10;             // surface azimuth angle, used to compute the sun's incidence angle to the surface
//		spa.atmos_refract = 0.5667;          // ?
//		spa.function      = SPA_ZA;          // find the zenith and azimuth angle
//		//spa.function      = SPA_ALL;
//
//	calculate_sunRay(sunRay, spa.zenith, spa.azimuth); // vector sunRay points to the sun
//	g_sunRayDir = Cvec3(-sunRay[1], sunRay[2], -sunRay[0]); 
//	Cvec3 yAxis (0,1,0);
//
//	Cvec3 groundCamDir = (-g_sunRayDir) - yAxis * dot(-g_sunRayDir, yAxis);
//	Cvec3 negZAxis (0,0,-1);
//	Cvec3 rotAxis = cross( negZAxis, groundCamDir);
//	double  sinTheta = norm( rotAxis );
//	double  theta = asin( sinTheta );
//	g_cameraRotMat = Matrix4::makeAxisRotation( theta, normalize(rotAxis) ); // so that the camera dir = sun dir
//	float distanceToEye =  42 ;  // the observer on the ground in Asian Game
//	double eyeHeight = 1.7;
//	Cvec3 eyeLocation = Cvec3(0, eyeHeight, distanceToEye);
//	g_eyeRbt = g_cameraRotMat * Matrix4::makeTranslation( eyeLocation );
//	g_invEyeRbt = inv(g_eyeRbt);
//
//	float  distanceToLight =  42 + 24 + 5 ;
//	float  heightToLight = 40 + 15;
//	float  separationLight = 5;
//	g_SceneRbt = Matrix4::makeAxisRotation( - PI / 3.0, Cvec3(0,1,0) );  //?????????????????? SgRbtNode ??????????????????
//	g_light1Pos = g_SceneRbt * Cvec4( separationLight/2.0, heightToLight, distanceToLight, 1.0);
//
//	const Cvec3 eyeLight1 = Cvec3(g_invEyeRbt * g_light1Pos ); // g_light1 position in eye coordinates
//	Cvec3 lightPos = eyeLight1;
//
//	g_radius = 1.0e-3;
//	bool isPointLight = false;
//	
//	
//
//
//
//
//	Cvec3 gl_FragCoord;
//	
//	Cvec2 pixelPos = Cvec2( gl_FragCoord[0], gl_FragCoord[1] );
//
//	Cvec3 surfaceColor = Cvec3( texelFetch(uColorTex, pixelPos, 0) ); // color
//  
//	float zWin  = texelFetch(uDepthTex, pixelPos, 0).r;  // zWin in non linear range [0,1]
//	// conversion range ([0,1] into NDC range  [-1,1]
//      
//	float zNDC = zWin * 2.0 - 1.0; // ( zWin = 1/2 zNDC + 1/2 by Viewport transformation: zWin in [0,1] )
//
//	// This uProjectionMatrix should be the same as the projectionMatrix when the background scene was rendered
//	const Matrix4 projMatrix = makeProjectionMatrix();
//	float alpha = projMatrix(2,2);
//	float beta = projMatrix(3,2); // glsl matrix is accessed with column first order
//	float zEye = - beta / ( zNDC + alpha ); 
//
//	Cvec3 volumeColor;
//  
//   
//	float scatCrossSection = computeScatCrossSection( g_radius, lambdaStart);	// scatCrossSection is independent of
//																				// lambda, so use lambdaStart
//  																			    
//	Cvec3 rayDir = normalize( vPosition ); // the direction from pixel to fragment position vPosition
//  
//     
//	float sigma_s = g_dropDensity * scatCrossSection;	// scattering coefficient [1/m]
//	float sigma_e = sigma_s; // assume that the water drop not absorb light, the extinction coeff = scat coeff. 
//
//	// compute the surface color at vPosition and attenuate the color through
//	// volume space between tmin and tmax.
//
//	bool isPointLight = false;
//	//Cvec3 lightPosOrDir;
//
//	Cvec3 dirFromSun  = -Cvec3( g_invEyeRbt * Cvec4(g_sunRayDir,0) );
//	
//
//	g_light1Pos = g_SceneRbt * Cvec4( separationLight/2.0, heightToLight, distanceToLight, 1.0);
//
//
//	Cvec3 volumeColor;
//
//	volumeColor = calculate_rainbowColor (g_radius,  surfaceColor, zEye, isPointLight, dirFromSun, 
//										  rayDir, sigma_s, sigma_e, LSpectrumSun_Plank);
//}