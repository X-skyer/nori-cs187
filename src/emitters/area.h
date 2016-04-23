#ifndef __AREA_H__
#define __AREA_H__

#pragma once
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

class AreaEmitter : public Emitter
{
public:

	AreaEmitter(const PropertyList& props);

	~AreaEmitter();

	Color3f sample(EmitterQueryRecord &lRec, const Point2f &sample) const;

	float pdf(const EmitterQueryRecord &lRec) const;

	Color3f eval(const EmitterQueryRecord &lRec) const;

	std::string toString() const;

private:
	Mesh*		m_mesh;								// underlying mesh the emitter points to
	Color3f		m_radiance;							// radiance emitted by each face of the mesh
	Transform	m_localToWorld, m_worldToLocal;		// transforms to be used to transform rays.
};

NORI_NAMESPACE_END

#endif
