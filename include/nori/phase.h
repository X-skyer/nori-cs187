#ifndef __PHASE_H__
#define __PHASE_H__

#include <nori/common.h>
#include <nori/object.h>

NORI_NAMESPACE_BEGIN

// Can be used for sampling medium information
struct MediumQueryRecord
{
	Vector3f wi;
	Vector3f wo;
	const PhaseFunction* phase;
	Point3f p;
	float pdf;

	MediumQueryRecord(const Vector3f& wo) : wo(wo)
	{		
	}

	MediumQueryRecord(const Vector3f& _wo, const Vector3f& _wi, const Vector3f& _pos, const PhaseFunction* _p)
	{
		wo = _wo;
		wi = _wi;
		p = _pos;
		phase = _p;
	}
};


class PhaseFunction : public NoriObject
{
public:
	virtual Color3f sample(MediumQueryRecord& rec, const Vector2f& _sample) const = 0;
	virtual float    pdf(MediumQueryRecord& rRec) const = 0;
	virtual Color3f  eval(MediumQueryRecord& rRec) const = 0;
	virtual EClassType getClassType() const { return EPhaseFunction; }
};

NORI_NAMESPACE_END

#endif
