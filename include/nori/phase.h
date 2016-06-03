#ifndef __PHASE_H__
#define __PHASE_H__

#include <nori/common.h>
#include <nori/object.h>
#include <nori/medium.h>

NORI_NAMESPACE_BEGIN

struct PhaseFunctionSamplingRecord
{
	const MediumSamplingRecord& mRec;
	Vector3f p;
	Vector3f wi;
	Vector3f wo;
	PhaseFunctionSamplingRecord(const MediumSamplingRecord& rec, const Vector3f& _p, const Vector3f& _wi, const Vector3f& _wo) :
		mRec(rec), p(_p), wi(_wi), wo(_wo)
	{

	}
};

class PhaseFunction : public NoriObject
{
public:
	// evaluate phase function for a sample
	virtual float eval(const PhaseFunctionSamplingRecord& pRec) const = 0;

	// return value/pdf
	virtual float sample(PhaseFunctionSamplingRecord& pRec, const Point2f& sample) const = 0;

	// return only phase function value, but return pdf in pdf parameter
	virtual float sample(PhaseFunctionSamplingRecord& pRec, float& pdf, const Point2f& sample) const = 0;

	// get pdf for a sample
	virtual float pdf(const PhaseFunctionSamplingRecord& pRec) const = 0;

	virtual EClassType getClassType() const { return EPhaseFunction; }

	virtual std::string toString() const = 0;

	void setParent(NoriObject* obj)
	{
		if (obj->getClassType() == EMedium)
		{
			m_medium = static_cast<Medium*>(obj);
		}
	}

	Medium* m_medium;
};

NORI_NAMESPACE_END

#endif
