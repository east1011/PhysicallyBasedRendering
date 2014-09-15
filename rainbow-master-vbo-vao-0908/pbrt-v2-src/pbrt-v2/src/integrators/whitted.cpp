
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// integrators/whitted.cpp*
#include "stdafx.h"
#include "integrators/whitted.h"
#include "intersection.h"
#include "paramset.h"

// WhittedIntegrator [ surface integrator ] Method Definitions
Spectrum WhittedIntegrator::Li(const Scene *scene,
        const Renderer *renderer, const RayDifferential &ray,
        const Intersection &isect, const Sample *sample, RNG &rng,
        MemoryArena &arena) const {
    Spectrum L(0.);
    // Compute emitted and reflected light at ray intersection point

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena); // return BSDF of the material of the hit primitive


    // Initialize common variables for Whitted integrator
    const pbrt::Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    Vector wo = -ray.d;

    // Compute emitted light at the intersection point. The emiited light exists if 
	// the intersection occured at an area light source; otherwise, it is zero. 

    L += isect.Le(wo);

    // Add contribution of each light source
    for (uint32_t i = 0; i < scene->lights.size(); ++i) {

        Vector wi;
        float pdf;
        VisibilityTester visibility;

        Spectrum Li = scene->lights[i]->Sample_L(p, isect.rayEpsilon,
            LightSample(rng), ray.time, &wi, &pdf, &visibility);

		// if light is a delta light, then "reflection" direction wi is set to the light direction with pdf =1

        if (Li.IsBlack() || pdf == 0.f) continue;

        Spectrum f = bsdf->f(wo, wi); 
		//
		// compute the fraction of light to be reflected from wi to wo.
		// ONLY non-delta BXDF, e.g. Lambertian (diffuse reflection)  contributes to f here.  
		
		

        if (!f.IsBlack() && visibility.Unoccluded(scene)) 
			// bsdf indicates that some of the incident
			// light from direction wi is in fact scattered to the direction wo

            L += f * Li * AbsDot(wi, n) *
                 visibility.Transmittance(scene, renderer,
                                          sample, rng, arena) / pdf;
    }

	// consider the specular reflection and the specular transmission;
	// you can get this component only when the light sources are NOT delta distribution.
	 // there's no possible way that the peaks of the two delta distributions (the glass 
     // and the point light) match.  it is not possible 
     // to get a highlight of a point light source in a mirror.

    if (ray.depth + 1 < maxDepth) {
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample,
                             arena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample,
                              arena);
    }
    return L;
}


WhittedIntegrator *CreateWhittedSurfaceIntegrator(const ParamSet &params)
{
    int maxDepth = params.FindOneInt("maxdepth", 5);
    return new WhittedIntegrator(maxDepth);
}


