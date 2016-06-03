#ifndef __MEDIUM_H__
#define __MEDIUM_H__

#include <nori/common.h>
#include <nori/object.h>

NORI_NAMESPACE_BEGIN

// Medium sampling record is used for housekeeping information about media effects
struct MediumSamplingRecord
{
	float					t;					// how much the ray has travelled in the medium so far
	Vector3f				p;					// the point where event was recorded
	Vector3f				wi;					// direction of incoming particle
	Color3f					transmittance;		// transmittance from [min_t, t] or [min_t, max_t]	 if medium was successfully sampled
	Color3f					m_sigmaA;
	Color3f					m_sigmaS;
	float					pdf_success;		// pdf of successfully sampling the meidum
	float					pdf_failure;		// pdf of not successfully sampling the medium
	const Medium*			medium;				// underlying medium that was sampled
	const PhaseFunction*	phase_function;
};

class Medium : public NoriObject
{
public:

	~Medium();

	virtual bool sample_distance(const Ray3f& ray, MediumSamplingRecord& mRec, const Point2f& sample) const = 0;

	virtual void eval(const Ray3f& ray, MediumSamplingRecord& mRec) const = 0;

	virtual Color3f eval_transmittance(const Ray3f& ray) const = 0;

	Color3f get_sigmaT() const { return m_sigmaT; }
	Color3f get_sigmaA() const { return m_sigmaA; }
	Color3f get_sigmaS() const { return m_sigmaS; }
	
	virtual EClassType getClassType() const { return EMedium; }

	virtual std::string toString() const = 0;

	void addChild(NoriObject* obj);

	PhaseFunction*	m_phase_funtion;

protected:
	Color3f			m_sigmaA;
	Color3f			m_sigmaS;
	Color3f			m_sigmaT;

};

NORI_NAMESPACE_END

#endif