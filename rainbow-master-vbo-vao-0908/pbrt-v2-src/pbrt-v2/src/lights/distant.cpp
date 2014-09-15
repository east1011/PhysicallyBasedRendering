
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


// lights/distant.cpp*
#include "stdafx.h"
#include "lights/distant.h"
#include "paramset.h"
#include "montecarlo.h"

// DistantLight Method Definitions
DistantLight::DistantLight(const string lightType, const Transform &light2world,
        const Spectrum &radiance, const Vector &dir)
    : Light(lightType, light2world) {
    lightDir = Normalize(LightToWorld(dir));
    L = radiance;
}


Spectrum DistantLight::Sample_L(const pbrt::Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    *wi = lightDir;
    *pdf = 1.f;
    visibility->SetRay(p, pEpsilon, *wi, time);
    return L;
}


Spectrum DistantLight::Power(const Scene *scene) const {
    pbrt::Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return L * M_PI * worldRadius * worldRadius;
}


DistantLight *CreateDistantLight(const string lightType, const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    pbrt::Point from = paramSet.FindOnePoint("from", pbrt::Point(0,0,0));
    pbrt::Point to = paramSet.FindOnePoint("to", pbrt::Point(0,0,1));
    Vector dir = from-to;
    return new DistantLight(lightType, light2world, L * sc, dir);
}


float DistantLight::Pdf(const pbrt::Point &, const Vector &) const {
    return 0.;
}


Spectrum DistantLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    // Choose point on disk oriented toward infinite light direction
    pbrt::Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(lightDir, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(ls.uPos[0], ls.uPos[1], &d1, &d2);
    pbrt::Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);

    // Set ray origin and direction for infinite light ray
    *ray = Ray(Pdisk + worldRadius * lightDir, -lightDir, 0.f, INFINITY,
               time);
    *Ns = (Normal)ray->d;

    *pdf = 1.f / (M_PI * worldRadius * worldRadius);
    return L; // DistanceLight::L, which is returned by Sample_L(), not by Le(), which is set to zero
	          // in delta lights. 
}


