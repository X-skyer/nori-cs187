#include <nori/phase.h>
#include <nori/warp.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN


class Isotropic : public PhaseFunction
{
public:

	Isotropic(const PropertyList& props)
	{
	}

	float eval(const PhaseFunctionSamplingRecord& pRec) const
	{
		return INV_FOURPI;
	}

	float pdf(const PhaseFunctionSamplingRecord& pRec) const
	{
		return INV_FOURPI;
	}

	float sample(PhaseFunctionSamplingRecord& pRec, const Point2f& _sample) const
	{
		pRec.wo = Warp::squareToUniformSphere(_sample);
		return 1.0f;
	}

	float sample(PhaseFunctionSamplingRecord& pRec, float& pdf, const Point2f& _sample) const
	{
		pRec.wo = Warp::squareToUniformSphere(_sample);
		pdf = Warp::squareToUniformSpherePdf(pRec.wo);
		return 1.0f;
	}


	std::string toString() const
	{
		return tfm::format("Iostropic Phase Function[]\n");
	}

private:

};

// henyey-greenstein function

class PhaseHG : public PhaseFunction
{
public:

	PhaseHG(const PropertyList& props)
	{
		m_g = props.getFloat("g", 0.0f);
	}

	float eval(const PhaseFunctionSamplingRecord& pRec) const
	{
		float temp = 1.0f + m_g * m_g + 2.0f * m_g * pRec.wi.dot(pRec.wo);
		return INV_FOURPI * (1 - m_g * m_g) / (temp * std::sqrt(temp));
	}

	float pdf(const PhaseFunctionSamplingRecord& pRec) const
	{
		return eval(pRec);
	}

	float sample(PhaseFunctionSamplingRecord& pRec, const Point2f& _sample) const
	{
		float cos_theta;
		if (fabsf(m_g) < 1e-3f)
		{
			cos_theta = 1.0f - 2 * _sample.x();
		}
		else
		{
			float sqrTerm = (1 - m_g * m_g) / (1 - m_g + 2 * m_g * _sample.x());
			cos_theta = (1 + m_g * m_g - sqrTerm * sqrTerm) / (2 * m_g);
		}

		float sin_theta = sqrtf(1.0f - cos_theta * cos_theta);
		float phi = 2.0f * M_PI * _sample.y();
		float sin_phi = sinf(phi);
		float cos_phi = cosf(phi);
		
		pRec.wo = Frame(-pRec.wi).toWorld(Vector3f(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta));
		return 1.0f;
	}

	float sample(PhaseFunctionSamplingRecord& pRec, float& pdf, const Point2f& sample) const
	{
		PhaseHG::sample(pRec, sample);
		pdf = PhaseHG::eval(pRec);
		return 1.0f;
	}

	std::string toString() const
	{
		return tfm::format("Henyey-Greenstein Phase Function[ g : %f\n]", m_g);
	}

private:
	float m_g;
};

NORI_REGISTER_CLASS(Isotropic, "isotropic");
NORI_REGISTER_CLASS(PhaseHG, "hg");
NORI_NAMESPACE_END