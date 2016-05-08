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
		m_rrStart = props.getInteger("rrStart", 5);
		m_maxDepth = props.getInteger("maxDepth", -1);
	}

	~PathIntegratorMts() {}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
	{
		Intersection isect;

		Color3f L(0.0f);
		Color3f throughput(1.0f);
		int depth = 0;
		Ray3f traced_ray = ray;

		while (depth < m_maxDepth || m_maxDepth == -1)
		{
			// Check if ray misses the scene
			if (!scene->rayIntersect(traced_ray, isect))
			{
				L += throughput * scene->getBackground(traced_ray, sampler->next2D());
				break;
			}

			// Check if direct emitter intersection
			if (isect.mesh->isEmitter())
			{
				EmitterQueryRecord eRec;
				eRec.ref = traced_ray.o;
				eRec.wi = traced_ray.d;
				eRec.n = isect.shFrame.n;
				L += throughput * isect.mesh->getEmitter()->eval(eRec);

				// Assume for now we dont bounce off the light sources.
				break;
			}

			const BSDF* bsdf = isect.mesh->getBSDF();			

			// Sample a reflection ray
			BSDFQueryRecord bRec(isect.toLocal(-traced_ray.d));
			Color3f f = bsdf->sample(bRec, sampler->next2D());
			Vector3f reflected_dir = isect.toWorld(bRec.wo);

			float cos_theta = fabsf(Frame::cosTheta(bRec.wo));
			throughput *= f * cos_theta;

			// Check if we've reached a zero throughput. No point in proceeding further.
			if (throughput.isZero())
				break;

			// Check for russian roulette
			if (depth > m_rrStart)
			{
				if (sampler->next1D() < 0.5f)
					break;
				else throughput *= 2.0f;
			}
			else if (depth > m_maxDepth && m_maxDepth != -1)
			{
				// forcibly terminate
				break;
			}

			// Propogate
			traced_ray = Ray3f(isect.p, reflected_dir, Epsilon, INFINITY);
			depth++;
		}
		
		return L;
	}

	std::string toString() const
	{
		return tfm::format("PathIntegratorMts[\nrrStart = %d\n]", m_rrStart);
	}

private:
	int m_rrStart;				// from which bounce should russian roulette start.
	int m_maxDepth;				// Fixed length cutoff

};

NORI_REGISTER_CLASS(PathIntegratorMts, "path_mats");
NORI_NAMESPACE_END