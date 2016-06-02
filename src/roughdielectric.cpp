/*
This file is part of Nori, a simple educational ray tracer

Copyright (c) 2015 by Wenzel Jakob

Nori is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License Version 3
as published by the Free Software Foundation.

Nori is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/distributions.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class RoughDielectric : public BSDF {
public:
	RoughDielectric(const PropertyList &propList) 
	{
		/* Interior IOR (default: BK7 borosilicate optical glass) */
		m_intIOR = propList.getFloat("intIOR", 1.5046f);

		/* Exterior IOR (default: air) */
		m_extIOR = propList.getFloat("extIOR", 1.000277f);

		m_alpha = propList.getFloat("alpha");

		// get distribution type
		const std::string dist_type = propList.getString("distribution");
		m_distribution = Distribution(dist_type, m_alpha);
		
		m_type = BsdfType::BSDF_ROUGHDIELECTRIC;
	}

	virtual Color3f eval(const BSDFQueryRecord &) const {
		/* Discrete BRDFs always evaluate to zero in Nori */
		return Color3f(0.0f);
	}

	virtual float pdf(const BSDFQueryRecord &) const {
		/* Discrete BRDFs always evaluate to zero in Nori */
		return 0.0f;
	}

	virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample, float optional_u) const 
	{
		// Sample a microfacet normal direction
		Vector3f m = m_distribution.sample(sample, bRec.pdf);
		float Fr = F(bRec.wi, m);

		// check if reflection or refraction
		// if TIR, Fr = 1, and hence this branch is always taken.
		if (optional_u < Fr)
		{
			bRec.wo = 2.0f * m.dot(bRec.wi) * m - bRec.wi;
			bRec.measure = EDiscrete;
		}
		else
		{
			// do refraction
			bool entering = Frame::cosTheta(bRec.wi) > 0.0f;
			float eta_i = m_extIOR;
			float eta_t = m_intIOR;

			if (!entering)
				std::swap(eta_i, eta_t);

			float eta = eta_i / eta_t;
			float c = bRec.wi.dot(m);

			float sign = bRec.wi.z() > 0.0f ? 1.0f : -1.0f;
			bRec.wo = (eta * c - sign * sqrtf(1.0f + eta * (c * c - 1.0f)))*m - eta * bRec.wi;
		}

		return m.dot(bRec.wi) * m_distribution.G(bRec.wi, bRec.wo, m) / (bRec.wi.z() * m.z());
	}

	virtual std::string toString() const {
		return tfm::format(
			"Rough Dielectric[\n"
			"  intIOR = %f,\n"
			"  extIOR = %f\n"
			"]",
			m_intIOR, m_extIOR);
	}
private:

	// Microfacet fresnel
	float F(const Vector3f& wi, const Vector3f& m) const
	{
		float eta_i = m_extIOR;
		float eta_t = m_intIOR;

		if (Frame::cosTheta(wi) < 0.0f)
			std::swap(eta_i, eta_t);

		float c = wi.dot(m);
		float term = (eta_t * eta_t / eta_i * eta_i) - 1.0f + c * c;
		if (term < 0.0f) return 1.0f;

		float g = sqrtf(term);
		float g_minus_c2 = (g - c) * (g - c);
		float g_plus_c2 = (g + c) * (g + c);
	}

	Distribution m_distribution;
	float m_intIOR, m_extIOR;
	float m_alpha;
};

NORI_REGISTER_CLASS(RoughDielectric, "roughdielectric");
NORI_NAMESPACE_END
