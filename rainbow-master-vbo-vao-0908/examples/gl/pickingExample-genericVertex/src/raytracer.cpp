/*
 A very basic raytracer example.
 Copyright (C) 2012  www.scratchapixel.com
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 - changes 02/04/13: fixed flag in ofstream causing a bug under Windows,
 added default values for M_PI and INFINITY
 - changes 24/05/13: small change to way we compute the refraction direction
 vector (eta=ior if we are inside and 1/ior if we are outside the sphere)
 
 Compile with the following command: c++ -o raytracer -O3 -Wall raytracer.cpp
 
 */


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

#include  "raytracer.h"
#include "spa.h"

std::ofstream messageFile("./messageFile.txt", std::ios::out);
	 

void render(std::vector< Sphere * > &spheres, DropVolume volumeBox,  Cvec3 eyePos, Cvec3 sunRay )
{
	unsigned width = 1024, height = 720;

	Cvec3  *image = new Cvec3 [width * height], *pixel = image;

	
	double fov = 100; // 50 degree with respect to the sun light direction (the view direction)
	                  // ==> The light deflected from drops by 130 deg up to 180 deg will contribute to the image on thnae
	                  //     retina

	
	double aspectratio = width / double(height);	

	// set up the image plane at z = -1, where z is the z axis of the camera coordinate system rotMat

	double imageY  = tan(PI * 0.5 * fov / double(180));
    double imageX = imageY * aspectratio;

	double radius = 1.0e-3; // 1.0 mm


	// Trace rays
	double xStep = 2 * imageX / double(width);
	double yStep = 2 * imageY / double (height);
	
	double y, x;

	/*
	// initialize the image array to black color

	for (unsigned j = 0; j < height; ++j) {
		
		for (unsigned i = 0; i < width; ++i, ++pixel) {
			*pixel =  Cvec3(0);
        
			
		} // for x
        // go to the next step
        
	} // for y


	// compute the color along the central vertical line on the image plane

	x = 0.0;
	y = imageY;

	for (unsigned j = 0; j < height; ++j) {
		
								
	  Cvec3 raydir(x, y - yStep/2.0, -1); // the ray direction with respect to the camera coordinate system
	  raydir.normalize();
		

	  double psi = acos ( dot ( -raydir, -sunRay ) ) * 180 / PI;			
		
         
	  messageFile <<  "In raytracer.cpp: ray point (" << x  << ", " << y - yStep/2.0 << ")" << endl;
	  messageFile <<  "ray dir = ";
	  messageFile << raydir << endl;

	  messageFile << "dot between globalRay and sunRay=" << dot ( -raydir, -sunRay ) << endl;

	  messageFile << "scattering angle=";
	  messageFile << psi << endl;

	  if ( psi >  130 &&  psi < 141.981 ) { 

	 // calculate the surface color while considering the attenuation through drops
					  
	    SurfVolColor color  = volumeBox.rayIntersect(  spheres, volumeBox, radius, eyePos, raydir, sunRay); 
			
		   // hitType  = 0: hit the background (including the ground) ;  =1: hit the sphere; =2: hit the volume and sphere; 
		   // =3: hit the  volume and sphere; 

	    image[width * j + width/2]  = ( color.surfaceColor + color.volumeColor) ;

		messageFile << "image[" <<  width/2   << "]" << "[" << j << "]" ;
        messageFile << image[width *j + width/2] << endl;


            ////////////////////////////////////////////////////////////////////
            // TODO: add Volume rendered result to each pixels. Vec3<T> type.
            ///////////////////////////////////////////////////////////////////
            
			//output << "pixel (" << x << ", " << y << ")" << endl;

	  }	
	
        // go to the next step
      y = y - yStep;

	} // for y

	*/


	y = imageY;
	 
	for (unsigned j = 0; j < height; ++j) {
		
		
		x = -imageX;

		for (unsigned i = 0; i < width; ++i, ++pixel) {
						
			Cvec3 raydir(x + xStep/2.0, y - yStep/2.0, -1); // the ray direction with respect to the camera coordinate system
			raydir.normalize();

			// get the glocal ray dir
			//Cvec3 globalRayDir = Cvec3( rotMat * Cvec4( raydir, 0.0) );

			double psi = acos ( dot ( -raydir, -sunRay ) ) * 180 / PI;			
		
         
			messageFile <<  "In raytracer.cpp: pixel  (" << i << ", " << j << ")" << endl;

			messageFile <<  "In raytracer.cpp: ray point (" << x << ", " << y << ")" << endl;
			messageFile <<  "ray dir = ";
			messageFile << raydir << endl;

			messageFile << "dot between globalRay and sunRay=" << dot ( -raydir, -sunRay ) << endl;

			messageFile << "scattering angle=";
			messageFile << psi << endl;

		 	// calculate the surface color while considering the attenuation through drops
					  
			SurfVolColor color  = volumeBox.rayIntersect(  spheres, volumeBox, radius, eyePos, raydir, sunRay); 
			
		   // hitType  = 0: hit the background (including the ground) ;  =1: hit the sphere; =2: hit the volume and sphere; 
		   // =3: hit the  volume and sphere; 

			*pixel =  color.surfaceColor + color.volumeColor ;

            ////////////////////////////////////////////////////////////////////
            // TODO: add Volume rendered result to each pixels. Vec3<T> type.
            ///////////////////////////////////////////////////////////////////
            
			//output << "pixel (" << x << ", " << y << ")" << endl;

            // go to the next pixel
            x = x + xStep;

			
		} // for x
        // go to the next step
        y = y - yStep;

	} // for y
   

	// Save result to a PPM image (keep these flags if you compile under Windows)
	std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);

	ofs << "P6\n" << width << " " << height << "\n255\n";
	for (unsigned i = 0; i < width * height; ++i) {

		ofs << (unsigned char)(std::min( double(1), image[i][0]) * 255) <<
		(unsigned char)(std::min(  double(1), image[i][1]) * 255) <<
		(unsigned char)(std::min( double (1), image[i][2]) * 255);
	}
	ofs.close();
	delete [] image;
}



//////////////////////////////////////////////////////////////////////////
//	THIS CODE USES RIGHT HAND COORDINATE!
/////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    spa_data spa;  //declare the SPA structure
    int result;
    Cvec3 sunRay;
    
    // Default values to calculate sun's position
    spa.year          = 2014;
    spa.month         = 1;
    spa.day           = 31;
    spa.hour          = 16;
    spa.minute        = 30;
    spa.second        = 30;
    spa.timezone      = 9.0;             // KST (GMT + 9)
    //spa.delta_ut1     =0;
    spa.delta_t       = 67;              // Doesn't it depend on location?  go to the related webpage to find the correct value.
    spa.longitude     = 126.9410634;     // Sogang University's longitude in decimal
    spa.latitude      = 37.5517132;      // Sogang University's latitude in decimal
    spa.elevation     = 200;             // m
    spa.pressure      = 820;             // mbar [milibar]
    spa.temperature   = 23;              // celcius
    //spa.slope         = 30;              // surface slope angle, used to compute the sun's incidence angle to the surface
    //spa.azm_rotation  = -10;             // surface azimuth angle, used to compute the sun's incidence angle to the surface
    spa.atmos_refract = 0.5667;          // ?
    spa.function      = SPA_ZA;          // find the zenith and azimuth angle
	//spa.function      = SPA_ALL;

   
	
    /////////////////////////////////////////////
    // Add objects
    /////////////////////////////////////////////
    
	
	std::vector<Sphere *> spheres;

	// position, radius, surface color, reflectivity, transparency, emission color
	spheres.push_back(new Sphere( Cvec3(0, -100, -2000), 1000, Cvec3(0.2), 0, 0.0)); // surface color =(2,2,2), reflection =0, transmission - 0.0, emission color = default =(0,0,0)
	spheres.push_back(new Sphere( Cvec3 (0, 0, -200), 100, Cvec3 (1.00, 0.32, 0.36), 1, 0.5));
	spheres.push_back(new Sphere( Cvec3 (5, -10, -150), 100, Cvec3(0.90, 0.76, 0.46), 1, 0.0));
	spheres.push_back(new Sphere( Cvec3(5, 0, -250), 100, Cvec3 (0.65, 0.77, 0.97), 1, 0.0));
	spheres.push_back(new Sphere( Cvec3 (-5.5, -10, -150), 100, Cvec3 (0.90, 0.90, 0.90), 1, 0.0));


	setSunLightRGBColor();


    get_time_and_location(spa.year, spa.month, spa.day, spa.hour, spa.minute, spa.timezone, spa.longitude, spa.latitude);
    
    // calculate sun's position
    result = spa_calculate(&spa);   // input is valid when result eq 0
    
	if ( result !=0 ) {
		messageFile   << "error in the calculation of sun position" << endl;
		exit(-1);
	}


    // calculate sunRay

    calculate_sunRay(sunRay, spa.zenith, spa.azimuth); // vector sunRay points to the sun
	messageFile  << "sunRay zenith angle = " << spa.zenith <<" sunRay azimuth" << spa.azimuth << endl;
	
	messageFile  << "sunRay in geocentric coord = (" << sunRay[0] << ", " << sunRay[1] << ", " << sunRay[2] << ")\n" << endl;

	g_sunRay = Cvec3(-sunRay[1], sunRay[2], -sunRay[0]); // rename the axes to 3D graphics convention : z => y, x => -z, y => -x => y=z, z = -x, x=-y
	// E.g: In the original coord system: (-10, 20, 5) [ azimth= south-west, polar = positive] =? (-20, 5, 10)
	
    
    messageFile  << "sunRay in Graphics coord = (" << g_sunRay[0] << ", " << g_sunRay[1] << ", " << g_sunRay[2] << ")\n" << endl;

	int dropDensity = 30000;
    //int dropDensity = 1000;
	double width=60; double height=30; double depth=5;
	
	DropVolume volumeBox (width, height, depth,  dropDensity);

	
	 // Establish the camera reference systtem relative to which the view volume will be considered

	// Find the rotation matrix R such that groundCamDir = R * (0,0,-1)
	//  (1) Find the axis   that is perpendicular both groundCamDir and (0,0,-1)
	//  (2) Find the angle theta between these two vectors
	//  (3) Find the rotation matrix which rotates the vector (0,0,-1) by theta about the axis u.
	
	// Set up the camera coordinate system so that the sun is right behind the camera. The camera points in the direction 
	// which is the ground projection of the sun direction.
	

	//Cvec3 camDir (0,0,-1); //  camDir is initially set to point horizontally to the negative Z direction
	                        // it can be changed. 

	// get the ground projection of camDir: camDir = upCamDir + groundCamDir = (camDir dot yAxis) yAxis + groundCamDir
	// The camera reference system consists of groundCamDir, yAxis, and the axis orthogonal to these axes.

	Cvec3 yAxis (0,1,0);
	Cvec3 zAxis (0,0,1);

	Cvec3 groundCamDir = (-g_sunRay) - yAxis * dot(-g_sunRay, yAxis);
	Cvec3 upCamDir = -g_sunRay - groundCamDir;

	Assert ( upCamDir == yAxis * dot(-g_sunRay, yAxis) );


	double eyeHeight = 1.7;
	double eyeX = 0.0;
	double eyeZ = 15;

	Cvec3 eyePos (eyeX,  eyeHeight, eyeZ) ; // the eye of the observer relative to the camera coordinate system.

	groundCamDir.normalize();

	Matrix4 rotMat;

	Cvec3 negZAxis (0,0,-1);

	if ( norm2( negZAxis - groundCamDir ) < 1.0e-8 ) { // no need to rotate; The rotation matrix will be identity
		rotMat =  Matrix4();
	}
	
	else {
		Cvec3 rotAxis = cross( negZAxis, groundCamDir); // rotAxis = zAxis x groundCamDir

	    double  sinTheta = norm( rotAxis ); // length(zAxis x camDir) = norm(zAxis) * norm(groundCamDir) * sin (theta) = sin (theta)
	    double  theta = asin( sinTheta );
		
        rotMat = Matrix4::makeAxisRotation( theta, normalize(rotAxis) );
	

     	messageFile << "groundCamDir:" ;
	    messageFile << groundCamDir << endl;

		messageFile << "rotation axis:";
	    messageFile << rotAxis << endl;

	    messageFile << "rot theta:";
	    messageFile << theta << endl;

	    messageFile << "rotated Axis:";
	    messageFile << Cvec3(  rotMat * Cvec4( negZAxis, 0.0 ) ) << endl;

	   Assert ( norm2 ( groundCamDir - Cvec3(  rotMat * Cvec4( negZAxis, 0.0 ) ) ) <= 1.0e-8 );
	   
	}

	// find the coordinates (ie. centers) of the spheres relative to the camera coordinate system

	// rotMat represents the camera coordinate system: globalCoord = robtMat * eyeCoord
	// eyeCoord = inv(rotMat) * globalCoord
	Matrix4 invRotMat = inv( rotMat);

	for (unsigned i = 0; i < spheres.size(); ++i) {
		
		spheres[i]->center = Cvec3( invRotMat * Cvec4( spheres[i]->center, 0.0) ) ;
	}

	    ///////////////////////////////////////////
    // Rendering starts
    //////////////////////////////////////////
    

	Cvec3 localSunRay ( invRotMat * Cvec4(g_sunRay, 0.0) );

	
    render( spheres, volumeBox, eyePos, localSunRay);
    
	
	return 0;
}
