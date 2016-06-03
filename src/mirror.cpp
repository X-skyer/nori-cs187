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

/// Ideal mirror BRDF
class Mirror : public BSDF {
public:
    Mirror(const PropertyList &props) {
		m_type = BsdfType::BSDF_MIRROR;
		m_tex_filename = props.getString("filename", "none");
		if (m_tex_filename != "none")
		{
			m_hasTexture = true;
			m_texture = Texture(m_tex_filename);
		}
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

		Color3f f(1.0f);
		if (m_hasTexture)
			f = m_texture.getval(bRec.uv.x(), bRec.uv.y());

        return f / Frame::cosTheta(bRec.wi);
    }

    virtual std::string toString() const {
        return "Mirror[]";
    }
private:
	std::string m_tex_filename;
};

NORI_REGISTER_CLASS(Mirror, "mirror");
NORI_NAMESPACE_END
