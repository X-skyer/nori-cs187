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
#include <nori/complexior.h>
#include <string>

NORI_NAMESPACE_BEGIN

extern Color3f fresnel_conductor(float cosI, const Color3f& eta, const Color3f& k);

/// Ideal mirror BRDF
/// But gets the data from actual materials
class SmoothConductor : public BSDF {
public:
	SmoothConductor(const PropertyList& props) {
		m_type = BsdfType::BSDF_MIRROR;
		const std::string material = props.getString("cior", "Au");
		m_tex_filename = props.getString("filename", "none");

		if (m_tex_filename != "none")
		{
			m_hasTexture = true;
			m_texture = Texture(m_tex_filename);
		}

		m_ior = cior_lookup(material);
	}

	virtual Color3f eval(const BSDFQueryRecord &) const {
		/* Discrete BRDFs always evaluate to zero in Nori */
		return Color3f(0.0f);
	}

	virtual float pdf(const BSDFQueryRecord &) const {
		/* Discrete BRDFs always evaluate to zero in Nori */
		return 0.0f;
	}

	virtual Color3f sample(BSDFQueryRecord &bRec, const Point2f &, float optional_u) const {
		if (Frame::cosTheta(bRec.wi) <= 0)
			return Color3f(0.0f);

		// Reflection in local coordinates
		bRec.wo = Vector3f(
			-bRec.wi.x(),
			-bRec.wi.y(),
			bRec.wi.z()
		);
		bRec.measure = EDiscrete;

		/* Relative index of refraction: no change */
		bRec.eta = 1.0f;
		bRec.pdf = 1.0f;

		Color3f Fr = fresnel_conductor(Frame::cosTheta(bRec.wo), m_ior.eta, m_ior.k);

		Color3f tex(1.0f);
		if (m_hasTexture)
			tex = m_texture.getval(bRec.uv.x(), bRec.uv.y());

		return Fr * tex / Frame::cosTheta(bRec.wi);
	}

	virtual std::string toString() const {
		return tfm::format("SmoothConductor[ Material : %s\neta=[%f,%f,%f]\nk=[%f,%f,%f]\n]", m_ior.name, m_ior.eta.x(), m_ior.eta.y(), m_ior.eta.z(), m_ior.k.x(), m_ior.k.y(), m_ior.k.z());
	}
private:
	ComplexIor m_ior;
	std::string m_tex_filename;
};

NORI_REGISTER_CLASS(SmoothConductor, "smoothconductor");
NORI_NAMESPACE_END
