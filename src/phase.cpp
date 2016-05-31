#include <nori/phase.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class IsotropicPhaseFunction : public PhaseFunction
{
public:
	
	IsotropicPhaseFunction(const PropertyList& props)
	{
		// Empty constructor
	}

	IsotropicPhaseFunction()
	{
		// Empty Destructor
	}

	// Assuming that the calling function fills in wo and position of scattering event.
	Color3f sample(MediumQueryRecord& mRec, const Vector2f& _sample) const
	{
		mRec.wi = Warp::squareToUniformSphere(_sample);
		mRec.pdf = INV_FOURPI;
		mRec.phase = this;
		return eval(mRec);
	}

	float pdf(MediumQueryRecord& mRec) const
	{
		return Warp::squareToUniformSpherePdf(Vector3f(0.0f));
	}

	Color3f eval(MediumQueryRecord& mRec) const
	{
		return Color3f(INV_FOURPI);
	}

	std::string toString() const
	{
		return tfm::format("Isotropic Phase Function : []");
	}
};

// Henyey-Greenstein Phase function
class AnisotropicPhaseFunction : public PhaseFunction
{
public:
	AnisotropicPhaseFunction(const PropertyList& props)
	{
		float g = props.getFloat("g", 0.0f);		
	}

	Color3f sample(MediumQueryRecord& mRec, const Vector2f& _sample) const
	{
		mRec.wi = Warp::squareToUniformSphere(_sample);
		mRec.pdf = INV_FOURPI;
		mRec.phase = this;
		return eval(mRec);
	}

	float pdf(MediumQueryRecord& mRec) const
	{
		return INV_FOURPI;
	}
	
	Color3f eval(MediumQueryRecord& mRec) const
	{
		// compute cos(theta)
		float costheta = -mRec.wi.dot(mRec.wo);
		float val = INV_FOURPI * (1.0f - m_g * m_g) / std::powf((1.0f + m_g * m_g - (2.0f * m_g * costheta)), 1.5f);
		return Color3f(val);
	}

	std::string toString() const
	{
		return tfm::format("Anisotropic Phase Function : \n[\
							g = %f]", m_g);
	}
	
private:
	float m_g;	
};


NORI_REGISTER_CLASS(IsotropicPhaseFunction, "isotropic");
NORI_REGISTER_CLASS(AnisotropicPhaseFunction, "aniostropic");
NORI_NAMESPACE_END