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

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

		m_type = BsdfType::BSDF_DIELECTRIC;
    }

    virtual Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    virtual float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample, float optional_u) const {
		
		float Fr = fresnel(Frame::cosTheta(bRec.wi), m_extIOR, m_intIOR);

		// do refraction
		bool entering = Frame::cosTheta(bRec.wi) > 0.0f;
		float eta_i = m_extIOR;
		float eta_t = m_intIOR;

		if (!entering)
			std::swap(eta_i, eta_t);

		// check if reflection or refraction
		if (sample.x() < Fr)
		{
			bRec.wo = Vector3f(
				-bRec.wi.x(),
				-bRec.wi.y(),
				bRec.wi.z()
			);
			bRec.measure = EDiscrete;

			/* Relative index of refraction: no change */
			bRec.eta = eta_i;

			// The wi can be under the surface and hence negative. however, we need to take only the absolute value.
			return (Color3f(1.0f) / fabsf(Frame::cosTheta(bRec.wi)));
		}
		else
		{
			// Compute transmitted direction
			float sini2 = 1.0f - (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wi));
			float eta = eta_i / eta_t;
			float sint2 = eta * eta * sini2;

			if (sint2 >= 1.0f)
			{
				bRec.wi = Vector3f(0.0f);
				return 0.0f;
			}

			float cost = sqrtf(std::max(0.0f, 1.0f - sint2));
			if (entering) cost = -cost;
			float sintOverSini = eta;

			// Set all the values
			bRec.wo = Vector3f(sintOverSini * -bRec.wi.x(), sintOverSini * -bRec.wi.y(), cost);
			bRec.eta = eta_t;
			bRec.measure = EDiscrete;
			bRec.pdf = 1.0f;

			// The (1-fr) term disappears because the probability of sampling this refraction is also (1-fr)
			// Hence the numerator term and the denominator term cancel out.
			return (square(eta_t) / square(eta_i)) / (fabsf(Frame::cosTheta(bRec.wo)));
		}
    }

    virtual std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }
private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
