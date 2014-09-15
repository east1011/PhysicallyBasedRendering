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


#ifndef raytracer_h
#define raytracer_h

#include  "dropvolume.h"
//#include "phaseFunction.h"


#define Assert(expr) \
    ((expr) ? (int) 0 : \
         fprintf(stderr, "Assertion \"%s\" failed in %s, line %d", \
               #expr, __FILE__, __LINE__))


void render(std::vector< Sphere * > &spheres, DropVolume volumeBox,  Cvec3 eyePos, Cvec3 sunRay);



#endif