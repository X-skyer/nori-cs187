#ifndef __NULL_H__
#define __NULL_H__

#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN


// A null bsdf is used for enclosing media inside some shape
// Light does not interact with it any way.

class Null : public BSDF
{
public:
	Null(const PropertyList& props)
	{
		m_type = BsdfType::BSDF_NULL;
	}

	Color3f eval(const BSDFQueryRecord &bRec) const
	{
		return 0.0f;
	}

	float pdf(const BSDFQueryRecord &bRec) const
	{
		return 0.0f;
	}
	
	Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample, float optional_u = 0.0f) const
	{
		bRec.wo = -bRec.wi;
		bRec.eta = 1.0f;
	}

	std::string toString() const
	{
		return tfm::format("Null BSDF");
	}

private:
};

NORI_NAMESPACE_END

#endif