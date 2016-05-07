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
		
		// Compute the reflection coefficient.
		float cosThetaI = Frame::cosTheta(bRec.wi);
		float fR = fresnel(cosThetaI, m_extIOR, m_intIOR);

		const Vector3f n(0, 0, 1.0f);

		// Compute the transmitted direction.
		float eta1 = m_extIOR;
		float eta2 = m_intIOR;

		// Check for correct internal and outer IORs
		if (cosThetaI < 0.0f)
			std::swap(eta1, eta2);

		// Compute the refracted direction.
		bRec.measure = EDiscrete;
		bRec.eta = eta2;

		// Compute refracted direction.
		// There is a chance of total internal reflection.
		// TODO: handle it.
		if (fR == 1.0f)
			return 0.0f;
		bRec.wo = -(eta1 / eta2) * (bRec.wi - n * bRec.wi.z()) - n * sqrt(1.0f - (square(eta1 / eta2) * (1.0f - square(bRec.wi.z()))));
		return (square(eta1) / square(eta2)) * (1.0f - fR) / fabsf(Frame::cosTheta(bRec.wo));		
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
