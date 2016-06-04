#include <nori/medium.h>
#include <nori/color.h>
#include <nori/proplist.h>
#include <memory>

NORI_NAMESPACE_BEGIN

class HomogeneousMedium : public Medium
{
public:

	HomogeneousMedium(const PropertyList& props)
	{
		m_sigmaA = props.getColor("sigmaA", Color3f(0.0f, 0.0f, 0.0f));
		m_sigmaS = props.getColor("sigmaS", Color3f(0.0f, 0.0f, 0.0f));
		m_sigmaT = m_sigmaA + m_sigmaS;
	}

	bool sample_distance(const Ray3f& ray, MediumSamplingRecord& mRec, const Point2f& sample) const
	{
		// choose the component of sigmaT to sample
		int component = floor(sample.x() * 3.0f);
		if (component == 3) component = 2;

		float sampled_sigmaT = m_sigmaT[component];

		// sample the distance based on this component only.
		float sampled_distance = -log(1 - sample.y()) / sampled_sigmaT;
		float distance_to_surface = ray.maxt - ray.mint;
		bool status;
		if (sampled_distance > distance_to_surface)
		{
			// failure because sampled a point more than the nearest surface 
			status = false;
			sampled_distance = distance_to_surface;
		}
		else
		{
			// success
			mRec.p = ray(sampled_distance);
			mRec.m_sigmaA = m_sigmaA;
			mRec.m_sigmaS = m_sigmaS;
			mRec.t = ray.mint + sampled_distance;
			mRec.medium = this;
			status = true;
		}

		// compute probabilities
		mRec.pdf_failure = 0;
		mRec.pdf_success = 0;
		for (int i = 0; i < 3; i++)
		{
			float temp = std::exp(-m_sigmaT[i] * sampled_distance);
			mRec.pdf_failure += temp;
			mRec.pdf_success += temp * m_sigmaT[i];
		}

		mRec.pdf_failure /= 3.0f;
		mRec.pdf_success /= 3.0f;
		mRec.transmittance = eval_transmittance(sampled_distance);
		mRec.phase_function = m_phase_funtion;
		return status;
	}

	// We assume that ray extends only within the medium
	void eval(const Ray3f& ray, MediumSamplingRecord& mRec) const
	{
		float distance = ray.maxt - ray.mint;

		mRec.pdf_failure = 0;
		mRec.pdf_success = 0;
		for (int i = 0; i < 3; i++)
		{
			float temp = std::exp(-m_sigmaT[i] * distance);
			mRec.pdf_failure += temp;
			mRec.pdf_success += temp * m_sigmaT[i];
		}

		mRec.pdf_failure /= 3.0f;
		mRec.pdf_success /= 3.0f;

		mRec.m_sigmaA = m_sigmaA;
		mRec.m_sigmaS = m_sigmaS;
		mRec.medium = this;
		mRec.transmittance = eval_transmittance(distance);
		mRec.phase_function = m_phase_funtion;
	}

	// For now we assume medium is fully surrounding everything
	Color3f eval_transmittance(const Ray3f& ray) const
	{
		Color3f ret(0.0f);
		float t = ray.maxt - ray.mint;

		for (int i = 0; i < 3; i++)
		{
			ret[i] = std::exp(-m_sigmaT[i] * t);
		}

		return ret;
	}

	Color3f eval_transmittance(float t) const
	{
		Color3f ret(0.0f);
		for (int i = 0; i < 3; i++)
		{
			ret[i] = std::exp(-m_sigmaT[i] * t);
		}
		return ret;
	}

	std::string toString() const
	{
		return tfm::format("Homogeneous Medium : [\n sigmaA : (%f,%f,%f)\n sigmaS : (%f,%f,%f) \n sigmaT : (%f,%f,%f", m_sigmaA[0], m_sigmaA[1], m_sigmaA[2], m_sigmaS[0], m_sigmaS[1], m_sigmaS[2], m_sigmaT[0], m_sigmaT[1], m_sigmaT[2]);
	}

private:
};

NORI_REGISTER_CLASS(HomogeneousMedium, "homogeneous");
NORI_NAMESPACE_END