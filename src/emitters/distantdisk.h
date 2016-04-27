#ifndef __DISTANT_DISK_H__
#define __DISTANT_DISK_H__

#pragma once
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

/**
	This class specifies a distant disk light which emits light in a cone of directions.
*/

class DistantDisk : public Emitter
{
public:

	DistantDisk(const PropertyList& props);

	~DistantDisk();

	Color3f sample(EmitterQueryRecord &lRec, const Point2f &sample, float optional_u) const;

	float pdf(const EmitterQueryRecord &lRec) const;

	Color3f eval(const EmitterQueryRecord &lRec) const;

	std::string toString() const;

private:
	Color3f m_radiance;
	float   m_thetaA;
	float   m_cosThetaMax;
	Transform m_localToWorld;
	Transform m_worldToLocal;
};


NORI_NAMESPACE_END
#endif

