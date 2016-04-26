/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Romain Prévost

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

#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

class AreaEmitter : public Emitter {
public:
    AreaEmitter(const PropertyList &props) {
		m_type = EmitterType::EMITTER_AREA;
		m_radiance = props.getColor("radiance");
    }

    virtual std::string toString() const {
        return tfm::format(
                "AreaLight[\n"
                "  radiance = %s,\n"
                "]",
                m_radiance.toString());
    }

	virtual Color3f eval(const EmitterQueryRecord & lRec) const {
        if(!m_mesh)
            throw NoriException("There is no shape attached to this Area light!");

		// Return radiance only from the outside.
		//if (lRec.n.dot(lRec.wi) < 0.0f) return m_radiance;
		//else return 0.0f;
		return m_radiance;
    }

    virtual Color3f sample(EmitterQueryRecord & lRec, const Point2f & sample) const {
        if(!m_mesh)
            throw NoriException("There is no shape attached to this Area light!");

		// Sample the underlying mesh for a position and normal.
		m_mesh->samplePosition(sample, lRec.p, lRec.n);

		// Construct the EmitterQueryRecord structure.
		lRec.wi = (lRec.p - lRec.ref).normalized();
		lRec.emitter = this;
		lRec.dist = (lRec.p - lRec.ref).norm();
		float pA = pdf(lRec);

		// convert pdf to solid angle measure
		// pW = pA * r^2 / cos(theta)
		float pW = pA * lRec.dist * lRec.dist / fabsf(lRec.n.dot(-lRec.wi));

		if (pW < 0.0f)
		{
			std::cout << "Caught here" << std::endl;
		}

		lRec.pdf = pW;
		
		// Return the appropriately weighted radiance term back
		if(lRec.pdf != 0.0f) return eval(lRec) / lRec.pdf;
		else return 0.0f;
    }

	// Returns probability with respect to Area
    virtual float pdf(const EmitterQueryRecord &lRec) const {
        if(!m_mesh)
            throw NoriException("There is no shape attached to this Area light!");

		return m_mesh->pdf();
    }


    virtual Color3f samplePhoton(Ray3f &ray, const Point2f &sample1, const Point2f &sample2) const {
        throw NoriException("To implement...");
    }

	// Get the parent mesh
	void setParent(NoriObject *parent)
	{
		auto type = parent->getClassType();
		if (type == EMesh)
			m_mesh = static_cast<Mesh*>(parent);
	}


protected:
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaEmitter, "area")
NORI_NAMESPACE_END