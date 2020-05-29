//
// Implementation for Yocto/Extension.
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

#include "yocto_extension.h"

#include <yocto_pathtrace/yocto_pathtrace.h>

#include <atomic>
#include <deque>
#include <future>
#include <memory>
#include <mutex>
using namespace std::string_literals;

// -----------------------------------------------------------------------------
// MATH FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::extension {
// import math symbols for use
using math::abs;
using math::acos;
using math::atan2;
using math::clamp;
using math::cos;
using math::exp;
using math::flt_max;
using math::fmod;
using math::fresnel_conductor;
using math::fresnel_dielectric;
using math::identity3x3f;
using math::invalidb3f;
using math::log;
using math::make_rng;
using math::max;
using math::min;
using math::pif;
using math::pow;
using math::pow2;
using math::sample_discrete_cdf;
using math::sample_discrete_cdf_pdf;
using math::sample_uniform;
using math::sample_uniform_pdf;
using math::sin;
using math::sqrt;
using math::zero2f;
using math::zero2i;
using math::zero3f;
using math::zero3i;
using math::zero4f;
using math::zero4i;

}  // namespace yocto::extension

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EXTENSION
// -----------------------------------------------------------------------------
namespace yocto::extension {

float Np(float phi, int p, float s, float gammaO, float gammaT) {
  float dphi = phi - Phi(p, gammaO, gammaT);
  while (dphi > pif) dphi -= 2 * pif;
  while (dphi < -pif) dphi += 2 * pif;
  return TrimmedLogistic(dphi, s, -pif, pif);
}

static float Mp(float cosThetaI, float cosThetaO, float sinThetaI,
    float sinThetaO, float v) {
  float a  = cosThetaI * cosThetaO / v;
  float b  = sinThetaI * sinThetaO / v;
  float mp = (v <= .1)
                 ? (exp(LogI0(a) - b - 1 / v + 0.6931f + log(1 / (2 * v))))
                 : (exp(-b) * I0(a)) / (sinh(1 / v) * 2 * v);
  return mp;
}

static std::array<vec3f, pMax + 1> Ap(float cosThetaO, float eta, vec3f normal,
    vec3f outging, float h, const vec3f& T) {
  std::array<vec3f, pMax + 1> ap;
  //⟨Compute p = 0 attenuation at initial cylinder intersection⟩
  // float cosGammaO = safeSqrt(1 - h * h);
  // float cosTheta = cosThetaO * cosGammaO;
  float f = fresnel_dielectric(eta, normal, outging);
  // float f = FrDielectric(cosTheta, 1.f, eta);
  ap[0] = {f, f, f};  // Spectrum
  //⟨Compute p = 1 attenuation term⟩
  ap[1] = pow2(1 - f) * T;
  //⟨Compute attenuation terms up to p = pMax⟩
  for (int p = 2; p < pMax; ++p) ap[p] = ap[p - 1] * T * f;
  //⟨Compute attenuation term accounting for remaining orders of scattering⟩
  ap[pMax] = ap[pMax - 1] * f * T / (vec3f(1.f) - T * f);
  return ap;
}

//⟨HairBSDF Method Definitions⟩ ≡
hair eval_hair(const yocto::pathtrace::material* material) {
  auto hdata    = hair{};
  hdata.h       = material->h;
  hdata.gammaO  = safeASin(material->h);
  hdata.eta     = material->eta;
  hdata.sigma_a = material->color;
  hdata.beta_m  = material->beta_m;
  hdata.beta_n  = material->beta_n;
  hdata.alpha   = material->alpha;

  auto v          = std::vector<float>(pMax + 1);
  auto s          = 0;
  auto sin2kAlpha = zero3f;
  auto cos2kAlpha = zero3f;

  //⟨Compute longitudinal variance from βm⟩ //roughness
  v[0] = pow2(0.726f * material->beta_m + 0.812f * pow2(material->beta_m) +
              3.7f * pow(material->beta_m, 20));
  v[1] = .25 * v[0];
  v[2] = 4 * v[0];
  for (int p = 3; p <= pMax; ++p) v[p] = v[2];

  //⟨Compute azimuthal logistic scale factor from βn⟩
  s = SqrtPiOver8 *
      (0.265f * material->beta_n + 1.194f * pow2(material->beta_n) +
          5.372f * pow(material->beta_n, 22));

  //⟨Compute α terms for hair scales⟩
  sin2kAlpha[0] = sin(material->alpha);
  cos2kAlpha[0] = SafeSqrt(1 - pow2(sin2kAlpha[0]));
  for (int i = 1; i < 3; ++i) {
    sin2kAlpha[i] = 2 * cos2kAlpha[i - 1] * sin2kAlpha[i - 1];
    cos2kAlpha[i] = pow2(cos2kAlpha[i - 1]) - pow2(sin2kAlpha[i - 1]);
  }
  return hdata;
}

}  // namespace yocto::extension
