#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

class Plastic : public BSDF
{
public:

	Plastic(const PropertyList& props)
	{

	}



	std::string toString() const
	{
		return tfm::format("Plastic BSDF :");
	}

private:
	Color3f m_kd;
	Color3f m_ks;
	float m_intIOR;
	float m_extIOR;
};

NORI_NAMESPACE_END