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
		const std::string dist_type = propList.getString("distribution", "beckmann");
		m_distribution = Distribution(dist_type, m_alpha);
		
		m_type = BsdfType::BSDF_ROUGHDIELECTRIC;
	}

	virtual Color3f eval(const BSDFQueryRecord& bRec) const {
		/*
		/* Discrete BRDFs always evaluate to zero in Nori */
		/*
		// Reflection
		float sign = bRec.wi.z() < 0.0f ? -1.0f : 1.0f;
		Vector3f w_hr = sign * (bRec.wo + bRec.wi).normalized();
		Color3f fr = fresnel(w_hr.dot(bRec.wi), m_extIOR, m_intIOR) * m_distribution.G(bRec.wi, bRec.wo, w_hr) * m_distribution.D(w_hr) * 0.25f / (bRec.wi.z() * bRec.wo.z());

		// Transmission
		float eta_i = m_extIOR, eta_t = m_intIOR;
		if (bRec.wi.z() < 0.0f)
			std::swap(eta_i, eta_t);

		Vector3f w_ht = -(eta_i * bRec.wi + eta_t * bRec.wo).normalized();

		float term = fabsf(bRec.wi.dot(w_ht)) * fabsf(bRec.wo.dot(w_ht)) / (fabsf(bRec.wi.z()) * fabsf(bRec.wo.z()));
		float fr_t = 1.0f - fresnel(w_ht.dot(bRec.wi), m_extIOR, m_intIOR);
		float denom = eta_i * bRec.wi.dot(w_ht) + eta_t * bRec.wo.dot(w_ht);
		Color3f ft = term * eta_t * eta_t * fr_t * m_distribution.G(bRec.wi, bRec.wo, w_ht) * m_distribution.D(w_ht) / (denom * denom);
		
		return fr;// +ft;*/
		
		//return 0.0f;
		Color3f fr, ft;
				
		float sign = bRec.wi.z() < 0.0f ? -1.0f : 1.0f;
		Normal3f w_h = sign * (bRec.wi + bRec.wo).normalized();
		Vector3f half_vector = bRec.wo + bRec.wi;
		if (half_vector.norm() < 1e-3f)
		{
			fr = 0.0f;
		}

		if (w_h.isZero()) fr = 0.0f;
		float D_r = m_distribution.D(w_h);
		Color3f F_r = fresnel(w_h.dot(bRec.wi), m_extIOR, m_intIOR);
		float G_r = m_distribution.G(bRec.wi, bRec.wo, w_h);
		
		fr = F_r * D_r * G_r / (4 * fabsf(Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo)));
		if (isnan(fr.x()) || isnan(fr.y()) || isnan(fr.z())) fr = Color3f(0.0f);		// might arise due 0/0

		// Transmission
		float eta_i = m_extIOR, eta_t = m_intIOR;
		if (bRec.wi.z() < 0.0f)
			std::swap(eta_i, eta_t);

		Vector3f w_ht = -(eta_i * bRec.wi + eta_t * bRec.wo).normalized();

		float term = fabsf(bRec.wi.dot(w_ht)) * fabsf(bRec.wo.dot(w_ht)) / (fabsf(bRec.wi.z()) * fabsf(bRec.wo.z()));
		float fr_t = 1.0f - fresnel(w_ht.dot(bRec.wi), m_extIOR, m_intIOR);
		float denom = eta_i * bRec.wi.dot(w_ht) + eta_t * bRec.wo.dot(w_ht);
		float G_t = m_distribution.G(bRec.wi, bRec.wo, w_ht);
		float D_t = m_distribution.D(w_ht);
		ft = term * eta_t * eta_t * fr_t * G_t * D_t / (denom * denom);
		if (isnan(ft.x()) || isnan(ft.y()) || isnan(ft.z())) ft = Color3f(0.0f);

		Color3f total = fr + ft;
		if (!total.isValid())
		{
			std::cout << "INSIDE EVAL : " << std::endl;
			std::cout << "--------------------------" << std::endl;
			std::cout << " fr : " << fr << std::endl;
			std::cout << " F : " << F_r << std::endl;
			std::cout << " D : " << D_r << std::endl;
			std::cout << " G : " << G_r << std::endl;
			std::cout << "--------------------------" << std::endl;			
			std::cout << " ft : " << ft;
			std::cout << " F : " << fr_t << std::endl;
			std::cout << " D : " << D_t << std::endl;
			std::cout << " G : " << G_t << std::endl;
		}
		return total;
	}

	virtual float pdf(const BSDFQueryRecord& bRec) const {
		
		/*
		// reflection pdf
		

		// transmission pdf
		float eta_i = m_extIOR, eta_t = m_intIOR;
		if (bRec.wi.z() < 0.0f)
			std::swap(eta_i, eta_t);

		Vector3f w_ht = -(eta_i * bRec.wi + eta_t * bRec.wo).normalized();
		float denom = eta_i * bRec.wi.dot(w_ht) + eta_t * bRec.wo.dot(w_ht);
		float jacobian_t = eta_t * eta_t * bRec.wo.dot(w_ht) / (denom * denom);
		float pdf_t = m_distribution.D(w_ht) * jacobian_t;

		return pdf_r;// +pdf_t;
		*/
		
		//return 0.0f;

		float sign = bRec.wi.z() < 0.0f ? -1.0f : 1.0f;
		Vector3f w_hr = sign * (bRec.wo + bRec.wi).normalized();
		float jacobian_r = 0.25f / fabsf(bRec.wo.dot(w_hr));
		float pdf_r = m_distribution.D(w_hr) * jacobian_r;
		
		float eta_i = m_extIOR, eta_t = m_intIOR;
		if (bRec.wi.z() < 0.0f)
			std::swap(eta_i, eta_t);

		Vector3f w_ht = -(eta_i * bRec.wi + eta_t * bRec.wo).normalized();
		float denom = eta_i * bRec.wi.dot(w_ht) + eta_t * bRec.wo.dot(w_ht);
		float jacobian_t = eta_t * eta_t * bRec.wo.dot(w_ht) / (denom * denom);
		float pdf_t = m_distribution.D(w_ht) * jacobian_t;

		return pdf_t + pdf_r;
	}

	virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample, float optional_u) const 
	{
		float pdf;
		Normal3f w_h = m_distribution.sample(_sample, pdf);
		if (pdf == 0.0f) return 0.0f;
		
		// compute Fresnel
		float Fr = fresnel(w_h.dot(bRec.wi), m_extIOR, m_intIOR);

		if (optional_u < Fr)
		{
			bRec.wo = 2.0f * w_h.dot(bRec.wi) * w_h - bRec.wi;
			bRec.measure = ESolidAngle;

			// check for sampled direction
			Vector3f half_vector = bRec.wo + bRec.wi;
			if (half_vector.norm() < 1e-3f)
			{
				return 0.0f;
			}

			float jacobian = 0.25f / fabsf(w_h.dot(bRec.wo));
			bRec.pdf = pdf * jacobian;
			Color3f e = eval(bRec);
			Color3f f = eval(bRec) / bRec.pdf;
			if (!f.isValid())
			{
				std::cout << "Catching here in REFL BRANCH" << std::endl;
				std::cout << "Wi : " << std::endl;
				std::cout << bRec.wi << std::endl;
				std::cout << "Wo : " << std::endl;
				std::cout << bRec.wo << std::endl;
				std::cout << "e " << e << std::endl;
				std::cout << "pdf : " << bRec.pdf << std::endl;
				std::cout << "jacobian" << std::endl;
				Vector3f w_hr = (bRec.wi + bRec.wo);
				std::cout << w_hr.x() << "," << w_hr.y() << w_hr.z() << std::endl;
				std::cout << "Norm : " << w_hr.norm();
			}
			return f;
		}
		else
		{
			bool entering = Frame::cosTheta(bRec.wi) > 0.0f;
			float eta_i = m_extIOR;
			float eta_t = m_intIOR;

			if (!entering)
				std::swap(eta_i, eta_t);

			float eta = eta_i / eta_t;
			float c = bRec.wi.dot(w_h);

			float sign = bRec.wi.z() > 0.0f ? 1.0f : -1.0f;
			float discriminant = 1.0f + eta * (c * c - 1.0f);
			if (discriminant < 0.0f)
			{
				// TIR BRANCH
				// FOR NOW ZERO
				return Color3f(0.0f);
			}

			bRec.wo = (eta * c - sign * sqrtf(1.0f + eta * (c * c - 1.0f)))*w_h - eta * bRec.wi;

			Vector3f w_ht = -(eta_i * bRec.wi + eta_t * bRec.wo).normalized();
			float denom = eta_i * bRec.wi.dot(w_ht) + eta_t * bRec.wo.dot(w_ht);
			float jacobian_t = eta_t * eta_t * fabsf(bRec.wo.dot(w_ht)) / (denom * denom);
			float pdf_t = m_distribution.D(w_ht) * jacobian_t;
			if (pdf_t == 0.0f) return 0.0f;

			bRec.pdf = pdf_t;
			Color3f e = eval(bRec);
			Color3f f = e / pdf_t;
			if (!f.isValid())
			{
				std::cout << "Catching here in REFR BRANCH" << std::endl;
				std::cout << "Wi : " << std::endl;
				std::cout << bRec.wi << std::endl;
				std::cout << "Wo : " << std::endl;
				std::cout << bRec.wo << std::endl;
				std::cout << "e " << e << std::endl;
				std::cout << "pdf : " << bRec.pdf << std::endl;
				std::cout << "jacobian" << jacobian_t << std::endl;
				Vector3f w_hr = (bRec.wi + bRec.wo);
				std::cout << w_hr.x() << "," << w_hr.y() << w_hr.z() << std::endl;
				std::cout << "Norm : " << w_hr.norm() << std::endl;
				std::cout << "F : " << f << std::endl;
			}
			return f;
		}
		
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

		float term_a = c * (g + c) - 1.0f;
		float term_b = c * (g - c) + 1.0f;
		float term_a2 = term_a * term_a;
		float term_b2 = term_b * term_b;

		return 0.5f * (g_minus_c2 / g_plus_c2) * (1.0f + term_a2 / term_b2);
	}

	Distribution m_distribution;
	float m_intIOR, m_extIOR;
	float m_alpha;
};

NORI_REGISTER_CLASS(RoughDielectric, "roughdielectric");
NORI_NAMESPACE_END
