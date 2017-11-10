#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/medium.h>
#include <nori/phase.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class VolPath : public Integrator
{
public:

	VolPath(const PropertyList& props)
	{
		m_maxDepth = props.getInteger("maxDepth", -1);
		m_rrStart = props.getInteger("rrStart", 5);
	}

	Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
	{
		Color3f L(0.0f);
		Color3f throughput(1.0f);
		
		Intersection isect;
		Ray3f tray = ray;
		int depth = 0;

		while (true)
		{
			// check for intersections with the scene
			bool foundIntersection = scene->rayIntersect(tray, isect);
			if (foundIntersection)
				tray.maxt = isect.t;

			const Medium* medium = ray.medium;
			MediumSamplingRecord mRec;
			if (medium != NULL)
				medium->sample(tray, mRec, sampler->next2D());

			if (mRec.medium != NULL)
			{
				// perform medium sampling.
			}
			else
			{
				// Surface branch was taken
			}
		}

		return L;
	}

	std::string toString() const
	{
		return tfm::format("VolPath Integrator ["
			"Max Depth : %d\n"
			"RR Start : %d\n",
			m_maxDepth, m_rrStart
		);
	}

private:
	int m_maxDepth;
	int m_rrStart;
};

NORI_REGISTER_CLASS(VolPath, "volpath");
NORI_NAMESPACE_END