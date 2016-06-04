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
#include <nori/complexior.h>
#include <nori/distributions.h>

NORI_NAMESPACE_BEGIN

Color3f fresnel_conductor(float cosI, const Color3f& eta, const Color3f& k)
{
	Color3f tmp = (eta*eta + k*k) * cosI*cosI;
	Color3f Rparl2 = (tmp - (2.f * eta * cosI) + 1) / (tmp + (2.f * eta * cosI) + 1);
	Color3f tmp_f = eta*eta + k*k;
	Color3f Rperp2 =
		(tmp_f - (2.f * eta * cosI) + cosI*cosI) /
		(tmp_f + (2.f * eta * cosI) + cosI*cosI);
	return (Rparl2 + Rperp2) / 2.f;
}

class RoughConductor : public BSDF {
public:
	RoughConductor(const PropertyList &propList) {
		/* RMS surface roughness */
		m_alpha = propList.getFloat("alpha", 0.1f);
		m_type = BsdfType::BSDF_MICROFACET;
		const std::string& cior = propList.getString("cior");
		const std::string& dist = propList.getString("distribution", "beckmann");
		m_tex_filename = propList.getString("filename", "none");

		if (m_tex_filename != "none")
		{
			m_hasTexture = true;
			m_texture = Texture(m_tex_filename);
		}

		m_ior = cior_lookup(cior);
		m_distribution = Distribution(dist, m_alpha);
	}

	/// Evaluate the BRDF for the given pair of directions
	/// Always assume that BSDFQueryRecord has directions in the local frame
	virtual Color3f eval(const BSDFQueryRecord &bRec) const {

		if (bRec.measure != ESolidAngle || Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0)
			return Color3f(0.0f);

		Normal3f w_h = (bRec.wi + bRec.wo).normalized();
		float D = m_distribution.D(w_h);
		Color3f F = fresnel_conductor(w_h.dot(bRec.wi), m_ior.eta, m_ior.k);
		float G = m_distribution.G(bRec.wi, bRec.wo, w_h);

		return F * D * G / (4 * (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo)));
	}

	/// Evaluate the sampling density of \ref sample() wrt. solid angles
	virtual float pdf(const BSDFQueryRecord &bRec) const {

		if (Frame::cosTheta(bRec.wo) <= 0.0f || Frame::cosTheta(bRec.wi) <= 0.0f) return 0.0f;
		Normal3f w_h = (bRec.wi + bRec.wo).normalized();
		float jacobian = 0.25f / (w_h.dot(bRec.wo));
		return m_distribution.pdf(w_h) * jacobian;
	}

	/// Sample the BRDF
	virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample, float optional_u) const {

		if (Frame::cosTheta(bRec.wi) <= 0)
			return Color3f(0.0f);

		bRec.measure = ESolidAngle;		

		float pdf;
		Normal3f w_h = m_distribution.sample(_sample, pdf);
		if (pdf == 0.0f) return 0.0f;
		bRec.wo = 2.0f * w_h.dot(bRec.wi) * w_h - bRec.wi;
		float jacobian = 0.25f / fabsf(w_h.dot(bRec.wo));
		bRec.pdf = pdf * jacobian;
		return eval(bRec) / bRec.pdf;
	}

	virtual std::string toString() const {
		return tfm::format(
			"Microfacet[\n"
			"  alpha = %f,\n",
			m_alpha
		);
	}
private:
	float m_alpha;
	ComplexIor m_ior;
	Distribution m_distribution;
	std::string m_tex_filename;
};

NORI_REGISTER_CLASS(RoughConductor, "roughconductor");
NORI_NAMESPACE_END
