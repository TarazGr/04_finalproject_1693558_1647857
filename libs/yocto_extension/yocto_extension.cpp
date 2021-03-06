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
#include <numeric>
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
using math::pi;
using math::pif;
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

hair hair_bsdf(const yocto::pathtrace::material* material, vec2f uv) {
  auto hdata    = hair();
  hdata.h       = -1 + 2 * uv.y;
  hdata.gammaO  = SafeASin(hdata.h);
  hdata.eta     = 1.55f;
  hdata.beta_m  = 0.25f;
  hdata.beta_n  = 0.3f;
  hdata.alpha   = 2.0f;
  hdata.sigma_a = SigmaAFromReflectance(material->color, hdata.beta_n);

  // Compute longitudinal variance from βm
  hdata.v.push_back(Sqr(0.726f * hdata.beta_m + 0.812f * Sqr(hdata.beta_m) +
                        3.7f * Pow<20>(hdata.beta_m)));
  hdata.v.push_back(.25f * hdata.v[0]);
  hdata.v.push_back(4.0f * hdata.v[0]);
  for (auto p = 3; p <= pMax; p++) hdata.v.push_back(hdata.v[2]);
  // Compute azimuthal logistic scale factor from βn
  hdata.s = SqrtPiOver8 * (0.265f * hdata.beta_n + 1.194f * Sqr(hdata.beta_n) +
                              5.372f * Pow<22>(hdata.beta_n));
  // Compute α terms for hair scales
  hdata.sin2kAlpha.x = sin(2.0f);
  hdata.cos2kAlpha.x = SafeSqrt(1 - Sqr(hdata.sin2kAlpha.x));
  for (auto i = 1; i < 3; i++) {
    hdata.sin2kAlpha[i] = 2 * hdata.cos2kAlpha[i - 1] * hdata.sin2kAlpha[i - 1];
    hdata.cos2kAlpha[i] = Sqr(hdata.cos2kAlpha[i - 1]) -
                          Sqr(hdata.sin2kAlpha[i - 1]);
  }
  return hdata;
}

float Mp(float cosThetaI, float cosThetaO, float sinThetaI, float sinThetaO,
    float v) {
  auto a  = cosThetaI * cosThetaO / v;
  auto b  = sinThetaI * sinThetaO / v;
  auto mp = (v <= 0.1f)
                ? (exp(LogI0(a) - b - 1 / v + 0.6931f + log(1 / (2 * v))))
                : (exp(-b) * I0(a)) / (sinh(1 / v) * 2 * v);
  return mp;
}

std::vector<vec3f> Ap(float cosThetaO, float eta, vec3f normal, vec3f outging,
    float h, const vec3f& T) {
  std::vector<vec3f> ap = std::vector<vec3f>();
  // Compute p = 0 attenuation at initial cylinder intersection
  auto cosGammaO = SafeSqrt(1.0f - h * h);
  auto cosTheta  = cosThetaO * cosGammaO;
  auto f         = FrDielectric(cosTheta, 1.0f, eta);
  ap.push_back(vec3f{f});  // Spectrum
  // Compute p = 1 attenuation term
  ap.push_back(Sqr(1.0f - f) * T);
  // Compute attenuation terms up to p = pMax
  for (auto p = 2; p < pMax; p++) ap.push_back(ap[p - 1] * T * f);
  // Compute attenuation term accounting for remaining orders of scattering
  ap.push_back(ap[pMax - 1] * f * T / (vec3f{1.0f} - T * f));
  return ap;
}

float Np(float phi, int p, float s, float gammaO, float gammaT) {
  auto dphi = phi - Phi(p, gammaO, gammaT);
  while (dphi > pif) dphi -= 2 * pif;
  while (dphi < -pif) dphi += 2 * pif;
  return TrimmedLogistic(dphi, s, -pif, pif);
}

vec3f eval_hair(const hair& bsdf, const vec3f& normal, const vec3f& outgoing,
    const vec3f& incoming) {
  // Compute hair coordinate system terms related to wo
  auto sinThetaO = outgoing.x;
  auto cosThetaO = SafeSqrt(1.0f - Sqr(sinThetaO));
  auto phiO      = atan2(outgoing.z, outgoing.y);
  // Compute hair coordinate system terms related to wi
  auto sinThetaI = incoming.x;
  auto cosThetaI = SafeSqrt(1.0f - Sqr(sinThetaI));
  auto phiI      = atan2(incoming.z, incoming.y);
  // Compute cos θt for refracted ray
  auto sinThetaT = sinThetaO / bsdf.eta;
  auto cosThetaT = SafeSqrt(1.0f - Sqr(sinThetaT));
  // Compute γt for refracted ray
  auto etap      = sqrt(bsdf.eta * bsdf.eta - Sqr(sinThetaO)) / cosThetaO;
  auto sinGammaT = bsdf.h / etap;
  auto cosGammaT = SafeSqrt(1 - Sqr(sinGammaT));
  auto gammaT    = SafeASin(sinGammaT);
  // Compute the transmittance T of a single path through the cylinder
  auto T = exp(-bsdf.sigma_a * (2.0f * cosGammaT / cosThetaT));
  // Evaluate hair BSDF
  auto               phi = phiI - phiO;
  std::vector<vec3f> ap  = Ap(cosThetaO, bsdf.eta, normal, outgoing, bsdf.h, T);
  auto               fsum = zero3f;  // the evaluated spectrum
  for (auto p = 0; p < pMax; p++) {
    // Compute sin θo and cos θo terms accounting for scales
    float sinThetaOp, cosThetaOp;
    if (p == 0) {
      sinThetaOp = sinThetaO * bsdf.cos2kAlpha.y -
                   cosThetaO * bsdf.sin2kAlpha.y;
      cosThetaOp = cosThetaO * bsdf.cos2kAlpha.y +
                   sinThetaO * bsdf.sin2kAlpha.y;
    }
    // Handle remainder of p values for hair scale tilt
    else if (p == 1) {
      sinThetaOp = sinThetaO * bsdf.cos2kAlpha.x +
                   cosThetaO * bsdf.sin2kAlpha.x;
      cosThetaOp = cosThetaO * bsdf.cos2kAlpha.x -
                   sinThetaO * bsdf.sin2kAlpha.x;
    } else if (p == 2) {
      sinThetaOp = sinThetaO * bsdf.cos2kAlpha.z +
                   cosThetaO * bsdf.sin2kAlpha.z;
      cosThetaOp = cosThetaO * bsdf.cos2kAlpha.z -
                   sinThetaO * bsdf.sin2kAlpha.z;
    } else {
      sinThetaOp = sinThetaO;
      cosThetaOp = cosThetaO;
    }
    // Handle out-of-range cos θi from scale adjustment
    cosThetaOp = abs(cosThetaOp);
    fsum += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, bsdf.v[p]) *
            ap[p] * Np(phi, p, bsdf.s, bsdf.gammaO, gammaT);
  }
  // Compute contribution of remaining terms after pMax
  fsum += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, bsdf.v[pMax]) *
          ap[pMax] / (2.0f * pif);
  if (abs(incoming.z) > 0) fsum /= abs(incoming.z);
  return fsum;
}

std::vector<float> ComputeApPdf(const hair& bsdf, float cosThetaO,
    const vec3f& normal, const vec3f& outgoing) {
  // Compute array of Ap values for cosThetaO
  auto sinThetaO = SafeSqrt(1.0f - cosThetaO * cosThetaO);
  // Compute cos θt for refracted ray
  auto sinThetaT = sinThetaO / bsdf.eta;
  auto cosThetaT = SafeSqrt(1.0f - Sqr(sinThetaT));
  // Compute γt for refracted ray
  auto etap      = sqrt(bsdf.eta * bsdf.eta - Sqr(sinThetaO)) / cosThetaO;
  auto sinGammaT = bsdf.h / etap;
  auto cosGammaT = SafeSqrt(1.0f - Sqr(sinGammaT));
  auto gammaT    = SafeASin(sinGammaT);
  // Compute the transmittance T of a single path through the cylinder
  auto               T  = exp(-bsdf.sigma_a * (2.0f * cosGammaT / cosThetaT));
  std::vector<vec3f> ap = Ap(cosThetaO, bsdf.eta, normal, outgoing, bsdf.h, T);
  // Compute Ap PDF from individual Ap terms
  std::vector<float> apPdf = std::vector<float>(pMax + 1);
  auto               first = ap.begin();
  auto               last  = ap.end();
  auto               sumY  = 0.0f;
  auto               op    = [](float s, const vec3f& ap) { return s + ap.y; };
  for (; first != last; first++) {
    sumY = op(sumY, *first);
  }
  for (auto i = 0; i <= pMax; i++) {
    apPdf[i] = ap[i].y / sumY;
  }
  return apPdf;
}

std::pair<vec3f, float> sample_hair(const hair& bsdf, const vec3f& normal,
    const vec3f& outgoing, const vec2f& rng) {
  // Compute hair coordinate system terms related to wo
  auto sinThetaO = outgoing.x;
  auto cosThetaO = SafeSqrt(1.0f - Sqr(sinThetaO));
  auto phiO      = atan2(outgoing.z, outgoing.y);
  // Derive four random samples from u2
  std::vector<vec2f> u = {DemuxFloat(rng.x), DemuxFloat(rng.y)};
  // Determine which term p to sample for hair scattering
  std::vector<float> apPdf = ComputeApPdf(bsdf, cosThetaO, normal, outgoing);
  int                p;
  for (p = 0; p < pMax; p++) {
    if (u[0].x < apPdf[p]) break;
    u[0].x -= apPdf[p];
  }
  // Sample Mp to compute θi
  // Rotate sin θo and cos θo to account for hair scale tilt
  float sinThetaOp, cosThetaOp;
  if (p == 0) {
    sinThetaOp = sinThetaO * bsdf.cos2kAlpha.y - cosThetaO * bsdf.sin2kAlpha.y;
    cosThetaOp = cosThetaO * bsdf.cos2kAlpha.y + sinThetaO * bsdf.sin2kAlpha.y;
  } else if (p == 1) {
    sinThetaOp = sinThetaO * bsdf.cos2kAlpha.x + cosThetaO * bsdf.sin2kAlpha.x;
    cosThetaOp = cosThetaO * bsdf.cos2kAlpha.x - sinThetaO * bsdf.sin2kAlpha.x;
  } else if (p == 2) {
    sinThetaOp = sinThetaO * bsdf.cos2kAlpha.z + cosThetaO * bsdf.sin2kAlpha.z;
    cosThetaOp = cosThetaO * bsdf.cos2kAlpha.z - sinThetaO * bsdf.sin2kAlpha.z;
  } else {
    sinThetaOp = sinThetaO;
    cosThetaOp = cosThetaO;
  }
  u[1].x        = max(u[1].x, float(1e-5));
  auto cosTheta = 1 + bsdf.v[p] *
                          log(u[1].x + (1 - u[1].x) * exp(-2.0f / bsdf.v[p]));
  auto sinTheta  = SafeSqrt(1.0f - Sqr(cosTheta));
  auto cosPhi    = cos(2.0f * pif * u[1].y);
  auto sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
  auto cosThetaI = SafeSqrt(1 - Sqr(sinThetaI));
  // Sample Np to compute ∆φ
  auto  etap      = sqrt(bsdf.eta * bsdf.eta - Sqr(sinThetaO)) / cosThetaO;
  auto  sinGammaT = bsdf.h / etap;
  auto  cosGammaT = SafeSqrt(1.0f - Sqr(sinGammaT));
  auto  gammaT    = SafeASin(sinGammaT);
  float dphi;
  if (p < pMax)
    dphi = Phi(p, bsdf.gammaO, gammaT) +
           SampleTrimmedLogistic(u[0].y, bsdf.s, -pif, pif);
  else
    dphi = 2.0f * pif * u[0].y;
  // Compute wi from sampled hair scattering angles
  auto phiI     = phiO + dphi;
  auto incoming = vec3f{
      sinThetaI, cosThetaI * cos(phiI), cosThetaI * sin(phiI)};
  // Compute PDF for sampled hair scattering direction wi
  auto pdf = 0.0f;
  for (auto p = 0; p < pMax; p++) {
    // Compute sin θo and cos θo terms accounting for scales
    float sinThetaOp, cosThetaOp;
    if (p == 0) {
      sinThetaOp = sinThetaO * bsdf.cos2kAlpha.y -
                   cosThetaO * bsdf.sin2kAlpha.y;
      cosThetaOp = cosThetaO * bsdf.cos2kAlpha.y +
                   sinThetaO * bsdf.sin2kAlpha.y;
    }
    // Handle remainder of p values for hair scale tilt
    else if (p == 1) {
      sinThetaOp = sinThetaO * bsdf.cos2kAlpha.x +
                   cosThetaO * bsdf.sin2kAlpha.x;
      cosThetaOp = cosThetaO * bsdf.cos2kAlpha.x -
                   sinThetaO * bsdf.sin2kAlpha.x;
    } else if (p == 2) {
      sinThetaOp = sinThetaO * bsdf.cos2kAlpha.z +
                   cosThetaO * bsdf.sin2kAlpha.z;
      cosThetaOp = cosThetaO * bsdf.cos2kAlpha.z -
                   sinThetaO * bsdf.sin2kAlpha.z;
    } else {
      sinThetaOp = sinThetaO;
      cosThetaOp = cosThetaO;
    }
    // Handle out-of-range cos θi from scale adjustment
    cosThetaOp = abs(cosThetaOp);
    pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, bsdf.v[p]) *
           apPdf[p] * Np(dphi, p, bsdf.s, bsdf.gammaO, gammaT);
  }
  pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, bsdf.v[pMax]) *
         apPdf[pMax] * (1.0f / (2.0f * pif));
  return {incoming, pdf};
}
}  // namespace yocto::extension
