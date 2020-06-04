# Hair Shading
- Veronica Birindelli, 1647857
- Gabriel Radu Taranciuc - 1693558

In this project we integrated new functionalities into yocto_pathtrace, implementing a hair scattering model based on the implementation described by Matt Pharr in prbt (https://www.pbrt.org/hair.pdf).
Some details that were omitted in the paper were instead taken from their source code (https://github.com/mmp/pbrt-v3).

## Theory background

To perform all the correct calculations, the paper gives some geometric descriptions of the prblem at hand. In particular, it provides tools to measure incident and outgoing directions when a ray intersects a hair on a given point.
The hair is assumed to be a curve, which is the shape obtained when a circle is swept along a Bézier curve, generating a cylinder. 
They consider two different kind of scattering:

	- Scattering on the longitudinal plane
		It considers the scattering on the side view of the curve, or more formally the length of the curve. The surface normals at a certain point are given by the surface normals of the circular cross-section at that point. The plane on which they lay is called the normal plane.
	- Scattering on the azimuthal plane
		It consider scattering on the curve-width direction.

### Directions at the intersection point
When a ray hits a hair, the directions of the incoming and outgoing rays are measured by two angles, θ and φ.
	- θ: is called the longitudinal angle, and measures the offset of a ray w.r.s. to the normal plane
	- φ; is called the azimuthal angle, and is found by projecting a direction ω into the normal plain, and computing its angle with the y axis
Another useful parameter is the h parameter. It parametrizes the circle's diameter of the hair section. Given a h for a ray that has intersected the curve, the angle γ between the direction ω and  the surface normal can be calculated. Trigonometry tells us that sin γ = h.

### Scattering from hair
The paper assumes the cuticle can be modeled as a rough dielectric cylinder, with scales all angled at the same angle α. The hair interior is treated as a homogeneous medium that only absorbes light. Scattering is assumed to be modeled accurately by a BSDF, assuming the light enters and exits at the same place.
The parameter p measures the number of path segments that light arriving at the hair follows inside the hair before being scattered back out. p=0 means R, p=1 means TT, p=2 means TRT, p=3 means TRRT and so on. 

The hair BSDF (in the code denoted by fsum) will be written as a sum over the p-terms. It will also be divided into terms that depend only on angles θ or φ (the latter is given by φo-φi): Mp(θo,θi), the longitudinal scattering function; Ap(ωo), the attenuation function; and Np(φ), the azimuthal scattering function. 
In the implementation just the first few terms of the sum are calculated explicitly, all the higher-order terms will be represented by a single term. The constant Pmax will control exacty how many therms there will be before switching over. As the paper suggested, we used a Pmax value of 3.
