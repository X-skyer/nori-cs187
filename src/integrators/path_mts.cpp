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

			// Else do NEE first
			const BSDF* bsdf = isect.mesh->getBSDF();

			/*
			EmitterQueryRecord eRec;
			
			// Sample one light among many.
			float pdf = 1.0f / scene->getLights().size();
			const Emitter* random_emitter = scene->getRandomEmitter(sampler->next1D());

			eRec.ref = isect.p;
			Color3f L_sampled_incoming = random_emitter->sample(eRec, sampler->next2D());
			pdf *= random_emitter->pdf(eRec);

			// Evalute the rendering equation
			Color3f partial_eval = bsdf->eval(BSDFQueryRecord(isect.toLocal(-traced_ray.d), isect.toLocal(eRec.wi), ESolidAngle)) * L_sampled_incoming * fmaxf(isect.shFrame.n.dot(eRec.wi), 0.0f) / pdf;
			if (pdf != 0.0f && !partial_eval.isZero())
			{
				float V = scene->rayIntersect(Ray3f(isect.p, eRec.wi, Epsilon, (1.0f - Epsilon) * eRec.dist)) ? 0.0f : 1.0f;
				L += throughput * partial_eval * V;
			}
			*/

			// Sample a reflection ray
			BSDFQueryRecord bRec(isect.toLocal(-traced_ray.d));
			Color3f f = bsdf->sample(bRec, sampler->next2D());
			Vector3f reflected_dir = isect.toWorld(bRec.wo);
			
			float cos_theta = fabsf(isect.geoFrame.n.dot(reflected_dir));
			throughput *= f * cos_theta;

			// Check if we've fa
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