
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_TEXTURES_IMAGEMAP_H
#define PBRT_TEXTURES_IMAGEMAP_H

// textures/imagemap.h*
#include "pbrt.h"
#include "texture.h"
#include "mipmap.h"
#include "paramset.h"
#include <map>

// TexInfo Declarations
struct TexInfo {
    TexInfo(const string &f, bool dt, float ma, ImageWrap wm, float sc, float ga)
        : filename(f), doTrilinear(dt), maxAniso(ma), wrapMode(wm), scale(sc), gamma(ga) { }
    string filename;
    bool doTrilinear;
    float maxAniso;
    ImageWrap wrapMode;
    float scale, gamma;
    bool operator<(const TexInfo &t2) const {
        if (filename != t2.filename) return filename < t2.filename;
        if (doTrilinear != t2.doTrilinear) return doTrilinear < t2.doTrilinear;
        if (maxAniso != t2.maxAniso) return maxAniso < t2.maxAniso;
        if (scale != t2.scale) return scale < t2.scale;
        if (gamma != t2.gamma) return gamma < t2.gamma;
        return wrapMode < t2.wrapMode;
    }
};



// ImageTexture Declarations
template <typename Tmemory, typename Treturn>
    class ImageTexture : public Texture<Treturn> {
public:

	//int m_width, m_height;     // by Moon Jung, 2014/4/5/
	//RGBSpectrum *m_texels; // by Moon Jung, 2014/4/5; m_texels could be float??

    // ImageTexture Public Methods
    ImageTexture(const string& name, const string & texClass,TextureMapping2D *m, const string &filename, bool doTri,
                 float maxAniso, ImageWrap wm, float scale, float gamma); // =>  this-filename = filename

    Treturn Evaluate(const DifferentialGeometry &) const;
    ~ImageTexture();

    static void ClearCache() {
        typename std::map<TexInfo, MIPMap<Tmemory> *>::iterator iter;
        iter = textures.begin();
        while (iter != textures.end()) {
            delete iter->second;
            ++iter;
        }
        textures.erase(textures.begin(), textures.end());
    }

	//private:
public: // changed by Moon Jung, 2014/4/19

    // ImageTexture Private Methods
    static MIPMap<Tmemory> *GetTexture(const string &filename,
        bool doTrilinear, float maxAniso, ImageWrap wm, float scale, float gamma);
   
	static void convertIn(const RGBSpectrum &from, RGBSpectrum *to,
                          float scale, float gamma) {
        *to = Pow(scale * from, gamma);
    }

    static void convertIn(const RGBSpectrum &from, float *to,
                          float scale, float gamma) {
        *to = powf(scale * from.y(), gamma); // get the brightness of the RGB color.
		                                     // We need to use a different convertIn() in order to treat GL_RED
		                                     // texture.
    }

    static void convertOut(const RGBSpectrum &from, Spectrum *to) {
        float rgb[3];
        from.ToRGB(rgb);
        *to = Spectrum::FromRGB(rgb);
    }
    static void convertOut(float from, float *to) {
        *to = from;
    }

    // ImageTexture Private Data
	string filename; 
    MIPMap<Tmemory> *mipmap; // mipmap has BlockedArray<T> **pyramid member to hold the texture data
	                         // mipmap->pyradmid[0] contains the original image.

    TextureMapping2D *mapping;
    static std::map<TexInfo, MIPMap<Tmemory> *> textures; // texture cache: this is a class variable, a sort of global variable
	                                                      // declaration, which should be initialized in imagemap.cpp file.

	
};


ImageTexture<float, float> *CreateImageFloatTexture(const string& name, const string & texClass,
													const Transform &tex2world,
        const TextureParams &tp);
ImageTexture<RGBSpectrum, Spectrum> *CreateImageSpectrumTexture(const string& name, const string & texClass,
																const Transform &tex2world,
        const TextureParams &tp);

#endif // PBRT_TEXTURES_IMAGEMAP_H
