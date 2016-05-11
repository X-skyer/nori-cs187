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

	// We don't assume anything about the visibility of points specified in 'ref' and 'p' in the EmitterQueryRecord.
	virtual Color3f eval(const EmitterQueryRecord & lRec) const {
        if(!m_mesh)
            throw NoriException("There is no shape attached to this Area light!");

		// This function call can be done by bsdf sampling routines.
		// Hence the ray was already traced for us - i.e a visibility test was already performed.
		// Hence just check if the associated normal in emitter query record and incoming direction are not backfacing
		if (lRec.n.dot(lRec.wi) < 0.0f)
			return m_radiance;
		else return 0.0f;
    }

    virtual Color3f sample(EmitterQueryRecord & lRec, const Point2f & sample, float optional_u) const {
        if(!m_mesh)
            throw NoriException("There is no shape attached to this Area light!");

		// Sample the underlying mesh for a position and normal.
		m_mesh->samplePosition(sample, lRec.p, lRec.n, optional_u);

		// Construct the EmitterQueryRecord structure.
		lRec.wi = (lRec.p - lRec.ref).normalized();
		lRec.emitter = this;
		lRec.dist = (lRec.p - lRec.ref).norm();
		lRec.pdf = pdf(lRec);
		
		// Return the appropriately weighted radiance term back
		// NOTE: We are not checking visibility here. It's the integrator's responsibility to check for the shadow ray test.
		if(lRec.pdf != 0.0f || fabsf(lRec.pdf) != INFINITY) return eval(lRec) / lRec.pdf;
		else return 0.0f;
    }

	// Returns probability with respect to solid angle given by all the information inside the emitterqueryrecord.
	// Assumes all information about the intersection point is already provided inside.
	// WARNING: Use with care. Malformed EmitterQueryRecords can result in undefined behavior. Plus no visibility is considered.
    virtual float pdf(const EmitterQueryRecord &lRec) const {
        if(!m_mesh)
            throw NoriException("There is no shape attached to this Area light!");

		Vector3f inv_wi = -lRec.wi;
		float costheta_here = fabsf(lRec.n.dot(inv_wi));
		float pW = m_mesh->pdf() * lRec.dist * lRec.dist / costheta_here;
		if (isnan(pW) || fabsf(pW) == INFINITY)
			return 0.0f;
		return pW;
    }


    virtual Color3f samplePhoton(Ray3f &ray, const Point2f &sample1, const Point2f &sample2, float another_rand_u) const {
		if (!m_mesh)
			throw NoriException("There is no shape attached to this area light");

		// global pdf
		float _pdf = 0.0f;

		// Sample a point to emit from the mesh
		EmitterQueryRecord eRec;
		m_mesh->samplePosition(sample1, eRec.p, eRec.n, another_rand_u);
		
		// Sample a direction
		Vector3f wi = Warp::squareToCosineHemisphere(sample2);

		Frame my_frame(eRec.n);
		Vector3f xfm_wi = my_frame.toWorld(wi);

		_pdf = (1.0f / m_mesh->pdf()) * (Warp::squareToCosineHemispherePdf(wi));

		// Get the ray out.
		ray = Ray3f(eRec.p, xfm_wi, Epsilon, INFINITY);

		// Return power
		return M_PI * m_mesh->totalSurfaceArea() * m_radiance;
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