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
#include <yocto_pathtrace/yocto_pathtrace.h>

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

using math::clamp;
using math::max;
using math::pif;
using math::pow2;

}  // namespace yocto::extension

// -----------------------------------------------------------------------------
// HIGH LEVEL API
// -----------------------------------------------------------------------------

namespace yocto::extension {
//⟨HairBSDF Constants⟩ ≡
static const int   pMax        = 3;
static const float SqrtPiOver8 = 0.626657069f;

struct hair {
  float h       = 0;
  float gammaO  = 0;
  float eta     = 0;
  float beta_m  = 0;
  float beta_n  = 0;
  float alpha   = 0;
  vec3f sigma_a = zero3f;

  std::vector<float> v = std::vector<float>();

  float s = 0;

  vec3f sin2kAlpha = zero3f;
  vec3f cos2kAlpha = zero3f;
};

inline float SafeSqrt(float x);
inline float SafeASin(float x);
inline float Phi(int p, float gammaO, float gammaT);
inline float Logistic(float x, float s);
inline float LogisticCDF(float x, float s);
inline float TrimmedLogistic(float x, float s, float a, float b);
inline float I0(float x);
inline float LogI0(float x);

inline float SafeSqrt(float x) {
  if (x >= -0.0001)
    return sqrt(max(float(0), x));
  else
    return 0;
}
inline float SafeASin(float x) {
  if (x >= -1.0001 && x <= 1.0001)
    return asin(clamp(x, -1.0, 1.0));
  else
    return 0;
}

inline float Phi(int p, float gammaO, float gammaT) {
  return 2 * p * gammaT - 2 * gammaO + p * pif;
}

inline float Logistic(float x, float s) {
  x = abs(x);
  return exp(-x / s) / (s * pow2(1.0f + exp(-x / s)));
}

inline float LogisticCDF(float x, float s) {
  return 1.0f / (1.0f + exp(-x / s));
}

inline float TrimmedLogistic(float x, float s, float a, float b) {
  return Logistic(x, s) / (LogisticCDF(b, s) - LogisticCDF(a, s));
}

inline float I0(float x) {
  auto    val   = 0.0f;
  auto    x2i   = 1.0f;
  int64_t ifact = 1;
  auto    i4    = 1;
  for (auto i = 0; i < 10; i++) {
    if (i > 1) ifact *= i;
    val += x2i / (i4 * pow2(ifact));
    x2i *= x * x;
    i4 *= 4;
  }
  return val;
}

inline float LogI0(float x) {
  if (x > 12)
    return x + 0.5 * (-log(2 * pif) + log(1 / x) + 1 / (8 * x));
  else
    return log(I0(x));
}

static uint32_t Compact1By1(uint32_t x) {
  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
  x &= 0x55555555;
  // x = --fe --dc --ba --98 --76 --54 --32 --10
  x = (x ^ (x >> 1)) & 0x33333333;
  // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
  x = (x ^ (x >> 2)) & 0x0f0f0f0f;
  // x = ---- ---- fedc ba98 ---- ---- 7654 3210
  x = (x ^ (x >> 4)) & 0x00ff00ff;
  // x = ---- ---- ---- ---- fedc ba98 7654 3210
  x = (x ^ (x >> 8)) & 0x0000ffff;
  return x;
}

static vec2f DemuxFloat(float f) {
  uint64_t v       = f * (1ull << 32);
  uint32_t bits[2] = {Compact1By1(v), Compact1By1(v >> 1)};
  return {bits[0] / float(1 << 16), bits[1] / float(1 << 16)};
}
static float SampleTrimmedLogistic(float u, float s, float a, float b) {
  auto k = LogisticCDF(b, s) - LogisticCDF(a, s);
  auto x = -s * log(1 / (u * k + LogisticCDF(a, s)) - 1);
  return clamp(x, a, b);
}

}  // namespace yocto::extension

#endif
