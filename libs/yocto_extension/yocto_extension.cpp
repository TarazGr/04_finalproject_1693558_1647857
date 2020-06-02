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

float Np(float phi, int p, float s, float gammaO, float gammaT) {
  auto dphi = phi - Phi(p, gammaO, gammaT);
  while (dphi > pif) dphi -= 2 * pif;
  while (dphi < -pif) dphi += 2 * pif;
  return TrimmedLogistic(dphi, s, -pif, pif);
}
float Mp(float cosThetaI, float cosThetaO, float sinThetaI, float sinThetaO,
    float v) {
  auto a  = cosThetaI * cosThetaO / v;
  auto b  = sinThetaI * sinThetaO / v;
  auto mp = (v <= .1) ? (exp(LogI0(a) - b - 1 / v + 0.6931f + log(1 / (2 * v))))
                      : (exp(-b) * I0(a)) / (sinh(1 / v) * 2 * v);
  if (!isinf(mp) && !isnan(mp))
    return mp;
  else
    return 1.0f;
}

std::vector<vec3f> Ap(float cosThetaO, float eta, vec3f normal, vec3f outging,
    float h, const vec3f& T) {
  std::vector<vec3f> ap = std::vector<vec3f>();
  // Compute p = 0 attenuation at initial cylinder intersection
  auto cosGammaO = SafeSqrt(1.0f - h * h);
  auto cosTheta  = cosThetaO * cosGammaO;
  // auto f = fresnel_dielectric(eta, normal, outging);
  auto f = FrDielectric(cosTheta, 1.0f, eta);
  ap.push_back(vec3f{f});  // Spectrum
  // Compute p = 1 attenuation term
  ap.push_back((1.0f - f) * (1.0f - f) * T);
  // Compute attenuation terms up to p = pMax
  for (auto p = 2; p < pMax; p++) ap.push_back(ap[p - 1] * T * f);
  // Compute attenuation term accounting for remaining orders of scattering
  ap.push_back(ap[pMax - 1] * f * T / (vec3f{1.0f} - T * f));
  return ap;
}

std::vector<float> ComputeApPdf(const hair& bsdf, float cosThetaO,
    const vec3f& normal, const vec3f& outgoing) {
  // Compute array of Ap values for cosThetaO
  auto sinThetaO = SafeSqrt(1 - cosThetaO * cosThetaO);
  // Compute cos θt for refracted ray
  auto sinThetaT = sinThetaO / bsdf.eta;
  auto cosThetaT = SafeSqrt(1 - sinThetaT * sinThetaT);
  // Compute γt for refracted ray
  auto etap = sqrt(bsdf.eta * bsdf.eta - sinThetaO * sinThetaO) / cosThetaO;
  auto sinGammaT = bsdf.h / etap;
  auto cosGammaT = SafeSqrt(1 - sinGammaT * sinGammaT);
  // Compute the transmittance T of a single path through the cylinder
  auto               T  = exp(-bsdf.sigma_a * (2 * cosGammaT / cosThetaT));
  std::vector<vec3f> ap = Ap(cosThetaO, bsdf.eta, normal, outgoing, bsdf.h, T);
  // Compute Ap PDF from individual Ap terms
  std::vector<float> apPdf = std::vector<float>(pMax + 1);
  auto               sumY  = 0.0f;
  auto               first = ap.begin();
  auto               last  = ap.end();
  while (first != last) {
    sumY += (bsdf.s + (*first).y);
    ++first;
  }
  for (auto i = 0; i <= pMax; i++) {
    apPdf[i] = ap[i].y / sumY;
  }
  return apPdf;
}

float sample_hair_pdf(const hair& bsdf, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming, const vec2f& rng) {
  auto sinThetaO = outgoing.x;
  auto cosThetaO = SafeSqrt(1 - sinThetaO * sinThetaO);
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
  // Sample Mp to compute θi
  // taken from
  // https://github.com/mmp/pbrt-v3/blob/master/src/materials/hair.cpp
  u[1].x        = max(u[1].x, float(1e-5));
  auto cosTheta = 1 +
                  bsdf.v[p] * log(u[1].x + (1 - u[1].x) * exp(-2 / bsdf.v[p]));
  auto sinTheta  = SafeSqrt(1 - cosTheta * cosTheta);
  auto cosPhi    = cos(2 * pif * u[1].y);
  auto sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
  auto cosThetaI = SafeSqrt(1 - sinThetaI * sinThetaI);
  // Sample Np to compute ∆φ
  auto  etap = sqrt(bsdf.eta * bsdf.eta - sinThetaO * sinThetaO) / cosThetaO;
  auto  sinGammaT = bsdf.h / etap;
  auto  cosGammaT = SafeSqrt(1 - sinGammaT * sinGammaT);
  auto  gammaT    = SafeASin(sinGammaT);
  float dphi;
  if (p < pMax)
    dphi = Phi(p, bsdf.gammaO, gammaT) +
           SampleTrimmedLogistic(u[0].y, bsdf.s, -pif, pif);
  else
    dphi = 2 * pif * u[0].y;
  // Compute PDF for sampled hair scattering direction wi
  auto pdf = 0.0f;
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
    pdf += Mp(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, bsdf.v[p]) *
           apPdf[p] * Np(dphi, p, bsdf.s, bsdf.gammaO, gammaT);
  }
  pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, bsdf.v[pMax]) *
         apPdf[pMax] * (1 / (2 * pif));
  return pdf;
  //// performs the same computation was we just implemented for hair_sample
  // auto sinThetaO = outgoing.x;
  // auto cosThetaO = SafeSqrt(1 - sinThetaO * sinThetaO);
  // auto phiO      = atan2(outgoing.z, outgoing.y);
  //// Compute hair coordinate system terms related to wi
  // auto sinThetaI = incoming.x;
  // auto cosThetaI = SafeSqrt(1 - sinThetaI * sinThetaI);
  // auto phiI      = atan2(incoming.z, incoming.y);
  //// Compute γt for refracted ray
  // auto etap = sqrt(bsdf.eta * bsdf.eta - sinThetaO * sinThetaO) / cosThetaO;
  // auto sinGammaT = bsdf.h / etap;
  // auto gammaT    = SafeASin(sinGammaT);

  // std::vector<float> apPdf = ComputeApPdf(bsdf, cosThetaO, normal, outgoing);
  // auto               phi   = phiI - phiO;
  // auto               pdf   = 0.0f;
  // for (auto p = 0; p < pMax; p++) {
  //  // Compute sin θi and cos θi terms accounting for scales
  //  float sinThetaOp, cosThetaOp;
  //  if (p == 0) {
  //    sinThetaOp = sinThetaO * bsdf.cos2kAlpha.y -
  //                 cosThetaO * bsdf.sin2kAlpha.y;
  //    cosThetaOp = cosThetaO * bsdf.cos2kAlpha.y +
  //                 sinThetaO * bsdf.sin2kAlpha.y;
  //  }
  //  // Handle remainder of p values for hair scale tilt
  //  else if (p == 1) {
  //    sinThetaOp = sinThetaO * bsdf.cos2kAlpha.x +
  //                 cosThetaO * bsdf.sin2kAlpha.x;
  //    cosThetaOp = cosThetaO * bsdf.cos2kAlpha.x -
  //                 sinThetaO * bsdf.sin2kAlpha.x;
  //  } else if (p == 2) {
  //    sinThetaOp = sinThetaO * bsdf.cos2kAlpha.z +
  //                 cosThetaO * bsdf.sin2kAlpha.z;
  //    cosThetaOp = cosThetaO * bsdf.cos2kAlpha.z -
  //                 sinThetaO * bsdf.sin2kAlpha.z;
  //  } else {
  //    sinThetaOp = sinThetaO;
  //    cosThetaOp = cosThetaO;
  //  }
  //  // Handle out-of-range cos θi from scale adjustment
  //  cosThetaOp = abs(cosThetaOp);
  //  pdf += Mp(cosThetaI, cosThetaOp, sinThetaI, sinThetaOp, bsdf.v[p]) *
  //         apPdf[p] * Np(phi, p, bsdf.s, bsdf.gammaO, gammaT);
  //}
  // pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, bsdf.v[pMax]) *
  //       apPdf[pMax] * (1 / (2 * pif));
  // return pdf;
}

hair hair_bsdf(const yocto::pathtrace::material* material, vec2f uv) {
  auto hdata    = hair();
  hdata.h       = -1 + 2 * uv.y;
  hdata.gammaO  = SafeASin(hdata.h);
  hdata.eta     = 1.55f;
  //hdata.sigma_a = material->color;
  hdata.sigma_a = SigmaAFromReflectance(material->color, 0.3f);
  hdata.beta_m  = 0.25f;
  hdata.beta_n  = 0.3f;
  hdata.alpha   = 2.0f;
  //⟨Compute longitudinal variance from βm⟩ //roughness
  /*hdata.v.push_back(
      (0.726f * 0.25f + 0.812f * 0.25f * 0.25f + 3.7f * Pow<20>(0.25f)) *
      (0.726f * 0.25f + 0.812f * 0.25f * 0.25f + 3.7f * Pow<20>(0.25f)));*/
  hdata.v.push_back(
    (0.726f * hdata.beta_m + 0.812f * (hdata.beta_m*hdata.beta_m) + 3.7f * Pow<20>(hdata.beta_m)) *
    (0.726f * hdata.beta_m + 0.812f * (hdata.beta_m*hdata.beta_m) + 3.7f * Pow<20>(hdata.beta_m))
  );
  hdata.v.push_back(.25f * hdata.v[0]);
  hdata.v.push_back(4.0f * hdata.v[0]);
  for (auto p = 3; p <= pMax; p++) hdata.v.push_back(hdata.v[2]);
  //⟨Compute azimuthal logistic scale factor from βn⟩
  /*hdata.s = SqrtPiOver8 *
            (0.265f * 0.3f + 1.194f * 0.3f * 0.3f + 5.372f * Pow<22>(0.3f));*/
  hdata.s = SqrtPiOver8 *
            (0.265f * hdata.beta_n + 1.194f *(hdata.beta_n*hdata.beta_n) + 5.372f * Pow<22>(hdata.beta_n));
  //⟨Compute α terms for hair scales⟩
  hdata.sin2kAlpha.x = sin(2.0f * pif / 180.0f);
  hdata.cos2kAlpha.x = SafeSqrt(1 - hdata.sin2kAlpha.x * hdata.sin2kAlpha.x);
  for (auto i = 1; i < 3; i++) {
    hdata.sin2kAlpha[i] = 2 * hdata.cos2kAlpha[i - 1] * hdata.sin2kAlpha[i - 1];
    hdata.cos2kAlpha[i] = hdata.cos2kAlpha[i - 1] * hdata.cos2kAlpha[i - 1] -
                          hdata.sin2kAlpha[i - 1] * hdata.sin2kAlpha[i - 1];
  }
  return hdata;
}

vec3f eval_hair(const hair& bsdf, const vec3f& normal, const vec3f& outgoing,
    const vec3f& incoming) {
  // Compute hair coordinate system terms related to wo
  auto sinThetaO = outgoing.x;
  auto cosThetaO = SafeSqrt(1 - sinThetaO * sinThetaO);
  auto phiO      = atan2(outgoing.z, outgoing.y);
  // Compute hair coordinate system terms related to wi
  auto sinThetaI = incoming.x;
  auto cosThetaI = SafeSqrt(1 - sinThetaI * sinThetaI);
  auto phiI      = atan2(incoming.z, incoming.y);
  // Compute cos θt for refracted ray
  auto sinThetaT = sinThetaO / bsdf.eta;
  auto cosThetaT = SafeSqrt(1 - sinThetaT * sinThetaT);
  // Compute γt for refracted ray
  auto etap = sqrt(bsdf.eta * bsdf.eta - sinThetaO * sinThetaO) / cosThetaO;
  auto sinGammaT = bsdf.h / etap;
  auto cosGammaT = SafeSqrt(1 - sinGammaT * sinGammaT);
  auto gammaT    = SafeASin(sinGammaT);
  // Compute the transmittance T of a single path through the cylinder
  auto T = exp(-bsdf.sigma_a * (2 * cosGammaT / cosThetaT));
  // Evaluate hair BSDF
  auto               phi = phiI - phiO;
  std::vector<vec3f> ap  = Ap(cosThetaO, bsdf.eta, normal, outgoing, bsdf.h, T);
  auto               fsum = zero3f;  // calcola l'assoribimento //spectrum
  for (auto p = 0; p < pMax; p++) {
    // Compute sin θi and cos θi terms accounting for scales
    float sinThetaOp, cosThetaOp;
    if (p == 0) {
      sinThetaOp = sinThetaO * bsdf.cos2kAlpha.y -
                   cosThetaO * bsdf.sin2kAlpha.y;
      cosThetaOp = cosThetaO * bsdf.cos2kAlpha.y +
                   sinThetaO * bsdf.sin2kAlpha.y;
    } else if (p == 1) {
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
          ap[pMax] / (2 * pif);
  if (abs(incoming.z) > 0) fsum /= abs(incoming.z);
  // if (fsum == zero3f) printf("fsum: {%f, %f, %f}\n", fsum.x, fsum.y, fsum.z);
  return fsum;
}

std::pair<vec3f, float> sample_hair(const hair& bsdf, const vec3f& normal,
    const vec3f& outgoing, const vec2f& rng) {
  // Compute hair coordinate system terms related to wo
  auto sinThetaO = outgoing.x;
  auto cosThetaO = SafeSqrt(1 - sinThetaO * sinThetaO);
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
  // Sample Mp to compute θi
  // taken from
  // https://github.com/mmp/pbrt-v3/blob/master/src/materials/hair.cpp
  u[1].x        = max(u[1].x, float(1e-5));
  auto cosTheta = 1 +
                  bsdf.v[p] * log(u[1].x + (1 - u[1].x) * exp(-2 / bsdf.v[p]));
  auto sinTheta  = SafeSqrt(1 - cosTheta * cosTheta);
  auto cosPhi    = cos(2 * pif * u[1].y);
  auto sinThetaI = -cosTheta * sinThetaOp + sinTheta * cosPhi * cosThetaOp;
  auto cosThetaI = SafeSqrt(1 - sinThetaI * sinThetaI);
  // Sample Np to compute ∆φ
  auto  etap = sqrt(bsdf.eta * bsdf.eta - sinThetaO * sinThetaO) / cosThetaO;
  auto  sinGammaT = bsdf.h / etap;
  auto  cosGammaT = SafeSqrt(1 - sinGammaT * sinGammaT);
  auto  gammaT    = SafeASin(sinGammaT);
  float dphi;
  if (p < pMax)
    dphi = Phi(p, bsdf.gammaO, gammaT) +
           SampleTrimmedLogistic(u[0].y, bsdf.s, -pif, pif);
  else
    dphi = 2 * pif * u[0].y;
  // Compute wi from sampled hair scattering angles
  auto phiI     = phiO + dphi;
  auto incoming = vec3f{
      sinThetaI, cosThetaI * cos(phiI), cosThetaI * sin(phiI)};
  // Compute PDF for sampled hair scattering direction wi
  auto pdf = 0.0f;
  for (auto p = 0; p < pMax; p++) {
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
    pdf += Mp(cosThetaIp, cosThetaO, sinThetaIp, sinThetaO, bsdf.v[p]) *
           apPdf[p] * Np(dphi, p, bsdf.s, bsdf.gammaO, gammaT);
  }
  pdf += Mp(cosThetaI, cosThetaO, sinThetaI, sinThetaO, bsdf.v[pMax]) *
         apPdf[pMax] * (1 / (2 * pif));
  return {incoming, pdf};
}
}  // namespace yocto::extension
