#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathIntegratorMts : public Integrator
{
public:

	PathIntegratorMts(const PropertyList& props)
	{
		m_rrStart = props.getInteger("rrStart", 3);
		m_maxDepth = props.getInteger("maxDepth", -1);
	}

	~PathIntegratorMts() {}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
	{
		Intersection isect;

		Color3f L(0.0f);
		int depth = 0;

		// Check if ray misses the scene
		if (!scene->rayIntersect(ray, isect))
			return scene->getBackground(ray, sampler->next2D());
	}

	std::string toString() const
	{
		return tfm::format("PathIntegratorMts[\nrrStart = %d\n]", m_rrStart);
	}

private:
	int m_rrStart;				// from which bounce should russian roulette start.
	int m_maxDepth;				// Fixed length cutoff

};

NORI_NAMESPACE_END