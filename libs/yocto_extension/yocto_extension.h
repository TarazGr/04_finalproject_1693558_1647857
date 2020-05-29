//
// # Yocto/Extension: Tiny Yocto/GL extension
//
//

//
// LICENSE:
//
// Copyright (c) 2020 -- 2020 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

#ifndef _YOCTO_EXTENSION_H_
#define _YOCTO_EXTENSION_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_image.h>
#include <yocto/yocto_math.h>

#include <atomic>
#include <future>
#include <memory>

// -----------------------------------------------------------------------------
// ALIASES
// -----------------------------------------------------------------------------
namespace yocto::extension {

// Namespace aliases
namespace ext = yocto::extension;
namespace img = yocto::image;

// Math defitions
using math::bbox3f;
using math::byte;
using math::frame3f;
using math::identity3x4f;
using math::ray3f;
using math::rng_state;
using math::vec2f;
using math::vec2i;
using math::vec3b;
using math::vec3f;
using math::vec3i;
using math::vec4f;
using math::vec4i;
using math::zero2f;
using math::zero3f;

}  // namespace yocto::pathtrace

// -----------------------------------------------------------------------------
// HIGH LEVEL API
// -----------------------------------------------------------------------------
namespace yocto::extension {
//⟨HairBSDF Constants⟩ ≡ 
static const int pMax = 3;
static const float SqrtPiOver8 = 0.626657069f;

//for hair computations
inline float safeSqrt(float x);
inline float safeASin(float x);
inline float Phi(int p, float gammaO, float gammaT);
inline float Logistic(float x, float s);
inline float LogisticCDF(float x, float s);
inline float TrimmedLogistic(float x, float s, float a, float b);
inline float I0(float x);
inline float LogI0(float x);

//General utility functions
inline float safeSqrt(float x) { if(x >= -0.0001) return std::sqrt(std::max(float(0), x)); 
    else return 0; };
inline float safeASin(float x) { if(x >= -1.0001 && x <= 1.0001) return std::asin(yocto::math::clamp(x, -1.0, 1.0)); 
    else return 0; }

//for hair
inline float Phi(int p, float gammaO, float gammaT) {
    return 2 * p * gammaT - 2 * gammaO + p * yocto::math::pif;
}

inline float Logistic(float x, float s) {
    x = std::abs(x);
    return std::exp(-x / s) / (s * yocto::math::pow2(1 + std::exp(-x / s)));
}

inline float LogisticCDF(float x, float s) {
    return 1 / (1 + std::exp(-x / s));
}

inline float TrimmedLogistic(float x, float s, float a, float b) {
    return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
}

inline float I0(float x) {
    float   val   = 0;
    float   x2i   = 1;
    int64_t ifact = 1;
    int     i4    = 1;
    // I0(x) \approx Sum_i x^(2i) / (4^i (i!)^2)
    for (int i = 0; i < 10; ++i) {
        if (i > 1) ifact *= i;
        val += x2i / (i4 * yocto::math::pow2(ifact));
        x2i *= x * x;
        i4 *= 4;
    }
    return val;
}

inline float LogI0(float x) {
    if (x > 12)
        return x + 0.5 * (-std::log(2 * yocto::math::pif) + std::log(1 / x) + 1 / (8 * x));
    else
        return std::log(I0(x));
}


}  // namespace yocto::pathtrace

#endif
