
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


// textures/marble.cpp*
#include "stdafx.h"
#include "textures/marble.h"

// MarbleTexture Method Definitions
Texture<float> *CreateMarbleFloatTexture(const string& name, const string & texClass,const Transform &tex2world,
        const TextureParams &tp) {
    return NULL;
}



MarbleTexture *CreateMarbleSpectrumTexture(const string& name, const string & texClass,const Transform &tex2world,
        const TextureParams &tp) {
    // Initialize 3D texture mapping _map_ from _tp_
    TextureMapping3D *map = new IdentityMapping3D(tex2world);
    return new MarbleTexture( name, texClass, tp.FindInt("octaves", 8),
        tp.FindFloat("roughness", .5f),
        tp.FindFloat("scale", 1.f),
        tp.FindFloat("variation", .2f),
        map);
}

