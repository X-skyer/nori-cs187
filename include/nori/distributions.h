#ifndef __DISTRIBUTIONS_H__
#define __DISTRIBUTIONS_H__

/***
	This header file contains many distribution functions for microfacet models.
	This results in cleaner code sharing among the BSDFs implements
*/

#include<nori/common.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

enum class Distributions
{
	BECKMANN,
	GGX,
	PHONG
};

class Distribution
{
public:
	Distribution() {}
	Distribution(const std::string& distribution, float alpha);

	Vector3f sample(const Vector2f& _sample, float& pdf) const;
	Vector3f sample(const Vector2f& _sample, float& pdf, float alpha) const;
	
	float D(const Vector3f& m) const;
	float D(const Vector3f& m, float alpha) const;

	float pdf(const Vector3f& m) const;
	float pdf(const Vector3f& m, float alpha) const;

	float G(const Vector3f& wi, const Vector3f& wo, const Vector3f& wm) const;
	float G(const Vector3f& wi, const Vector3f& wo, const Vector3f& wm, float alpha) const;

private:
	
	Distributions m_distribution;
	float smithBeckmannG1(const Vector3f& v, const Vector3f& m) const;
	float smithBeckmannG1(const Vector3f& v, const Vector3f& m, float alpha) const;

	float smithPhongG1(const Vector3f& v, const Vector3f& m) const;
	float smithPhongG1(const Vector3f& v, const Vector3f& m, float alpha) const;

	float smithGgxG1(const Vector3f& v, const Vector3f& m) const;
	float smithGgxG1(const Vector3f& v, const Vector3f& m, float alpha) const;

	float m_alpha;
};

NORI_NAMESPACE_END


#endif
