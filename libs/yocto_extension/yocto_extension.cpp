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
using math::pi;
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
float Mp(float cosThetaI, float cosThetaO, float sinThetaI, float sinThetaO,
    float v) {
  float a  = cosThetaI * cosThetaO / v;
  float b  = sinThetaI * sinThetaO / v;
  float mp = (v <= .1)
                 ? (exp(LogI0(a) - b - 1 / v + 0.6931f + log(1 / (2 * v))))
                 : (exp(-b) * I0(a)) / (sinh(1 / v) * 2 * v);
  return mp;
}

std::vector<vec3f> Ap(float cosThetaO, float eta, vec3f normal, vec3f outging,
    float h, const vec3f& T) {
  std::vector<vec3f> ap = std::vector<vec3f>(pMax + 1);
  // Compute p = 0 attenuation at initial cylinder intersection
  // float cosGammaO = safeSqrt(1 - h * h);
  // float cosTheta = cosThetaO * cosGammaO;
  float f = fresnel_dielectric(eta, normal, outging);
  // float f = FrDielectric(cosTheta, 1.f, eta);
  ap[0] = vec3f{f};  // Spectrum
  // Compute p = 1 attenuation term
  ap[1] = pow2(1 - f) * T;
  // Compute attenuation terms up to p = pMax
  for (int p = 2; p < pMax; p++) ap[p] = ap[p - 1] * T * f;
  // Compute attenuation term accounting for remaining orders of scattering
  ap[pMax] = ap[pMax - 1] * f * T / (vec3f{1} - T * f);
  return ap;
}

std::vector<float> ComputeApPdf(const hair& bsdf, float cosThetaO,
    const vec3f& normal, const vec3f& outgoing) {
  // Compute array of Ap values for cosThetaO
  float sinThetaO = SafeSqrt(1 - cosThetaO * cosThetaO);
  // Compute cos θt for refracted ray
  float sinThetaT = sinThetaO / bsdf.eta;
  float cosThetaT = SafeSqrt(1 - pow2(sinThetaT));
  // Compute γt for refracted ray
  float etap = std::sqrt(bsdf.eta * bsdf.eta - pow2(sinThetaO)) / cosThetaO;
  float sinGammaT = bsdf.h / etap;
  float cosGammaT = SafeSqrt(1 - pow2(sinGammaT));
  float gammaT    = SafeASin(sinGammaT);
  // Compute the transmittance T of a single path through the cylinder
  vec3f              T  = exp(-bsdf.sigma_a * (2 * cosGammaT / cosThetaT));
  std::vector<vec3f> ap = Ap(cosThetaO, bsdf.eta, normal, outgoing, bsdf.h, T);
  // Compute Ap PDF from individual Ap terms
  std::vector<float> aPdf = std::vector<float>(pMax + 1);
  return aPdf;
}
hair hair_bsdf(const yocto::pathtrace::material* material, vec2f uv) {
  auto hdata    = hair{};
  hdata.h       = -1 + 2 * uv.y;
  hdata.gammaO  = SafeASin(hdata.h);
  hdata.eta     = material->eta;
  hdata.sigma_a = material->color;
  hdata.beta_m  = material->beta_m;
  hdata.beta_n  = material->beta_n;
  hdata.alpha   = material->alpha;

  //⟨Compute longitudinal variance from βm⟩ //roughness
  hdata.v[0] = pow2(0.726f * material->beta_m +
                    0.812f * pow2(material->beta_m) +
                    3.7f * pow(material->beta_m, 20));
  hdata.v[1] = .25 * hdata.v[0];
  hdata.v[2] = 4 * hdata.v[0];
  for (int p = 3; p <= pMax; ++p) hdata.v[p] = hdata.v[2];

  //⟨Compute azimuthal logistic scale factor from βn⟩
  hdata.s = SqrtPiOver8 *
            (0.265f * material->beta_n + 1.194f * pow2(material->beta_n) +
                5.372f * pow(material->beta_n, 22));

  //⟨Compute α terms for hair scales⟩
  hdata.sin2kAlpha[0] = sin(material->alpha);
  hdata.cos2kAlpha[0] = SafeSqrt(1 - pow2(hdata.sin2kAlpha[0]));
  for (int i = 1; i < 3; ++i) {
    hdata.sin2kAlpha[i] = 2 * hdata.cos2kAlpha[i - 1] * hdata.sin2kAlpha[i - 1];
    hdata.cos2kAlpha[i] = pow2(hdata.cos2kAlpha[i - 1]) -
                          pow2(hdata.sin2kAlpha[i - 1]);
  }
  return hdata;
}

vec3f eval_hair(const hair& bsdf, const vec3f& normal, const vec3f& outgoing,
    const vec3f& incoming) {
  // Compute hair coordinate system terms related to wo
  float sinThetaO = outgoing.x;
  float cosThetaO = SafeSqrt(1 - pow2(sinThetaO));
  float phiO      = atan2(outgoing.z, outgoing.y);
  // Compute hair coordinate system terms related to wi
  float sinThetaI = incoming.x;
  float cosThetaI = SafeSqrt(1 - pow2(sinThetaO));
  float phiI      = atan2(incoming.z, incoming.y);
  // Compute cos θt for refracted ray
  float sinThetaT = sinThetaO / bsdf.eta;
  float cosThetaT = SafeSqrt(1 - pow2(sinThetaT));
  // Compute γt for refracted ray
  float etap      = sqrt(bsdf.eta * bsdf.eta - pow2(sinThetaO)) / cosThetaO;
  float sinGammaT = bsdf.h / etap;
  float cosGammaT = SafeSqrt(1 - pow2(sinGammaT));
  float gammaT    = SafeASin(sinGammaT);
  // Compute the transmittance T of a single path through the cylinder
  vec3f T = exp(-bsdf.sigma_a * (2 * cosGammaT / cosThetaT));
  // Evaluate hair BSDF
  float              phi = phiI - phiO;
  std::vector<vec3f> ap  = Ap(cosThetaO, bsdf.eta, normal, outgoing, bsdf.h, T);
  vec3f              fsum(0);  // calcola l'assoribimento //spectrum
  for (int p = 0; p < pMax; p++) {
    // Compute sin θi and cos θi terms accounting for scales
    float sinThetaIp, cosThetaIp;
    if (p == 0) {
      sinThetaIp = sinThetaI * bsdf.cos2kAlpha.y +
                   cosThetaI * bsdf.sin2kAlpha.y;
      cosThetaIp = cosThetaI * bsdf.cos2kAlpha.y -
                   sinThetaI * bsdf.sin2kAlpha.y;
    } else if (p == 1) {
      sinThetaIp = sinThetaI * bsdf.cos2kAlpha.x +
                   cosThetaI * bsdf.sin2kAlpha.x;
      cosThetaIp = cosThetaI * bsdf.cos2kAlpha.x -
                   sinThetaI * bsdf.sin2kAlpha.x;
    } else if (p == 2) {
      sinThetaIp = sinThetaI * bsdf.cos2kAlpha.z +
                   cosThetaI * bsdf.sin2kAlpha.z;
      cosThetaIp = cosThetaI * bsdf.cos2kAlpha.z -
                   sinThetaI * bsdf.sin2kAlpha.z;
    } else {
      sinThetaIp = sinThetaI;
      cosThetaIp = cosThetaI;
    }
    // Handle out-of-range cos θi from scale adjustment
    cosThetaIp = abs(cosThetaIp);
    fsum += Mp(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, bsdf.v[p]) *
            ap[p] * Np(phi, p, bsdf.s, bsdf.gammaO, gammaT);
  }
  // Compute contribution of remaining terms after pMax
  fsum += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, bsdf.v[pMax]) *
          ap[pMax] / (2 * pif);

  // WARNING: calcolare il AbscosTheta di wi. Per ora ho deciso che e'
  // l'absolute value di cosThetaI, ma sono in dubbio if (AbsCosTheta(wi) > 0)
  // fsum /= AbsCosTheta(wi);
  if (abs(cosThetaI) > 0) fsum /= abs(cosThetaI);
  return fsum;
}

vec3f sample_hair(const hair& bsdf, const vec3f& normal, const vec3f& outgoing,
    const vec2f& rng) {
  // Compute hair coordinate system terms related to wo
  float sinThetaO = outgoing.x;
  float cosThetaO = SafeSqrt(1 - pow2(sinThetaO));
  float phiO      = atan2(outgoing.z, outgoing.y);
  // Derive four random samples from u2
  std::vector<vec2f> u = {DemuxFloat(rng.x), DemuxFloat(rng.y)};
  // Determine which term p to sample for hair scattering
  std::vector<float> apPdf = ComputeApPdf(bsdf, cosThetaO, normal, outgoing);
  int                p;
  for (p = 0; p < pMax; ++p) {
    if (u[0][0] < apPdf[p]) break;
    u[0][0] -= apPdf[p];
  }
  // Rotate $\sin \thetao$ and $\cos \thetao$ to account for hair scale tilt
  float sinThetaOp, cosThetaOp;
  if (p == 0) {
    sinThetaOp = sinThetaO * bsdf.cos2kAlpha[1] -
                 cosThetaO * bsdf.sin2kAlpha[1];
    cosThetaOp = cosThetaO * bsdf.cos2kAlpha[1] +
                 sinThetaO * bsdf.sin2kAlpha[1];
  } else if (p == 1) {
    sinThetaOp = sinThetaO * bsdf.cos2kAlpha[0] +
                 cosThetaO * bsdf.sin2kAlpha[0];
    cosThetaOp = cosThetaO * bsdf.cos2kAlpha[0] -
                 sinThetaO * bsdf.sin2kAlpha[0];
  } else if (p == 2) {
    sinThetaOp = sinThetaO * bsdf.cos2kAlpha[2] +
                 cosThetaO * bsdf.sin2kAlpha[2];
    cosThetaOp = cosThetaO * bsdf.cos2kAlpha[2] -
                 sinThetaO * bsdf.sin2kAlpha[2];
  } else {
    sinThetaOp = sinThetaO;
    cosThetaOp = cosThetaO;
  }
  // Sample Mp to compute θi
  // taken from
  // https://github.com/mmp/pbrt-v3/blob/master/src/materials/hair.cpp
  u[1][0]        = max(u[1][0], float(1e-5));
  float cosTheta = 1 + bsdf.v[p] *
                           log(u[1][0] + (1 - u[1][0]) * exp(-2 / bsdf.v[p]));
  float sinTheta  = SafeSqrt(1 - pow2(cosTheta));
  float cosPhi    = cos(2 * pi * u[1][1]);
  float sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
  float cosThetaI = SafeSqrt(1 - pow2(sinThetaI));
  // Sample Np to compute ∆φ
  float etap      = sqrt(bsdf.eta * bsdf.eta - pow2(sinThetaO)) / cosThetaO;
  float sinGammaT = bsdf.h / etap;
  float cosGammaT = SafeSqrt(1 - pow2(sinGammaT));
  float gammaT    = SafeASin(sinGammaT);
  float dphi;
  if (p < pMax)
    dphi = Phi(p, bsdf.gammaO, gammaT) +
           SampleTrimmedLogistic(u[0][1], bsdf.s, -pif, pif);
  else
    dphi = 2 * pif * u[0][1];
  // Compute wi from sampled hair scattering angles
  float phiI     = phiO + dphi;
  auto  incoming = vec3f{
      sinThetaI, cosThetaI * cos(phiI), cosThetaI * sin(phiI)};
  // Compute PDF for sampled hair scattering direction wi
}
}  // namespace yocto::extension
