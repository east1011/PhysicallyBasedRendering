//
//  sphere.h
//  rainbow_raytracer
//
//  Created by jd on 2014. 2. 4..
//  Copyright (c) 2014. All rights reserved.
//

#ifndef rainbow_raytracer_sphere_h
#define rainbow_raytracer_sphere_h

#include "cvec.h"

class Sphere
{
public:
	Cvec3  center;                         /// position of the sphere
	double radius, radius2;                      /// sphere radius and radius^2
	Cvec3  surfaceColor, emissionColor;    /// surface color and emission (light)
	double transparency, reflection;             /// surface transparency and reflectivity

	Sphere(const Cvec3 &c, const double &r, const Cvec3 &sc = 0,
           const double &refl = 0, const double &transp = 0, const Cvec3  &ec = 0) :
      center(c), radius(r), radius2(r * r), surfaceColor(sc), emissionColor(ec),
      transparency(transp), reflection(refl)
	  {  }


	// compute a ray-sphere intersection using the geometric solution
	bool intersect(const Cvec3 &rayorig, const Cvec3 &raydir,  double *t0 = NULL, double *t1 = NULL) const
	{
		Cvec3 l = center - rayorig;
		double tca = dot(l, raydir);
		if (tca < 0) return false;
		double d2 = dot(l, l) - tca * tca;
		if (d2 > radius2) return false;
		double thc = sqrt(radius2 - d2);
		if (t0 != NULL && t1 != NULL) {
			*t0 = tca - thc;
			*t1 = tca + thc;
		}
        
		return true;
	}
};


#endif
