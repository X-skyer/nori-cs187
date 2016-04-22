#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectEmsIntegrator : public Integrator
{
public:
	DirectEmsIntegrator(const PropertyList& prop)
	{
		// Empty constructor for now.
	}
	
	~DirectEmsIntegrator()
	{
		// Empty destructor for now.	
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
	{

	}
	
	std::string toString() const
	{
		return "DirectEmsIntegrator[]";
	}

private:

};


NORI_REGISTER_CLASS(DirectEmsIntegrator, "direct_ems");
NORI_NAMESPACE_END