#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectMatsIntegrator : public Integrator
{
public:
	DirectMatsIntegrator(const PropertyList& prop)
	{
		// Empty constructor for now.
	}

	~DirectMatsIntegrator()
	{
		// Empty destructor for now.	
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
	{
		return 0.0f;
	}

	std::string toString() const
	{
		return "DirectMatsIntegrator[]";
	}

private:

};


NORI_REGISTER_CLASS(DirectMatsIntegrator, "direct_mats");
NORI_NAMESPACE_END