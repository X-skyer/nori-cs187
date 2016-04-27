#ifndef __POINT_LIGHT_H__
#define __POINT_LIGHT_H__

#pragma once
#include <nori/emitter.h>

NORI_NAMESPACE_BEGIN

/**
	This class defines an isotropic point light source that is defined by its position and power.
*/

class PointLight : public Emitter
{
public:

	PointLight(const PropertyList& prop);

	// Methods that have to be implemented for all types of emitters.
	Color3f sample(EmitterQueryRecord &lRec, const Point2f &sample, float optional_u) const;

	float pdf(const EmitterQueryRecord &lRec) const;

	Color3f eval(const EmitterQueryRecord &lRec) const;

	std::string toString() const;

private:
	Color3f		m_power;
	Vector3f	m_position;
};

NORI_NAMESPACE_END

#endif