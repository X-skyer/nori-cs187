#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathIntegratorMats : public Integrator
{
public:
	PathIntegratorMats(const PropertyList& props)
	{
		m_rrStart = props.getInteger("rrStart", 10);
		m_maxDepth = props.getInteger("maxDepth", -1);
	}

	~PathIntegratorMats() {}

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
				L += throughput * scene->getBackground(traced_ray);
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
				//break;
			}

			const BSDF* bsdf = isect.mesh->getBSDF();

			// Sample a reflection ray
			BSDFQueryRecord bRec(isect.toLocal(-traced_ray.d));
			bRec.uv = isect.uv;

			Color3f f = bsdf->sample(bRec, sampler->next2D(), sampler->next1D());
			if (!f.isValid())
			{
				std::cout << "Catching here" << std::endl;
			}
			Vector3f reflected_dir = isect.toWorld(bRec.wo);

			float cos_theta = fabsf(Frame::cosTheta(bRec.wo));
			throughput *= f * cos_theta;

			// Check if we've fa
			if (throughput.isZero())
				break;

			// Check for russian roulette
			if (depth > m_rrStart)
			{
				float prob_cont = throughput.getLuminance();
				if (sampler->next1D() < prob_cont)
					throughput *= (1.0f / prob_cont);
				else break;
			}
			else if (depth > m_maxDepth && m_maxDepth != -1)
			{
				// forcibly terminate
				break;
			}

			// Propogate
			traced_ray = Ray3f(isect.p, reflected_dir, Epsilon, INFINITY);
			depth++;

			if (!L.isValid())
			{
				std::cout << "Catching" << std::endl;
			}
		}

		return L;
	}

	std::string toString() const
	{
		return tfm::format("PathIntegratorMats[\nrrStart = %d\n]", m_rrStart);
	}

private:
	int m_rrStart;				// from which bounce should russian roulette start.
	int m_maxDepth;				// Fixed length cutoff


};

NORI_REGISTER_CLASS(PathIntegratorMats, "path_mats")
NORI_NAMESPACE_END