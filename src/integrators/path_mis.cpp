#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathIntegratorMis : public Integrator
{
public:

	enum class DirectSamplingStrategy
	{
		SAMPLE_ALL_LIGHTS,
		SAMPLE_ONE_LIGHT
	};

	PathIntegratorMis(const PropertyList& props)
	{
		m_rrStart = props.getInteger("rrStart", 5);
		m_maxDepth = props.getInteger("maxDepth", -1);
		
		// Which strategy to use for our NEE scheme
		std::string strategy = props.getString("strategy", "sample_one_light");
		if (strategy == "sample_all_lights")
			m_strategy = DirectSamplingStrategy::SAMPLE_ALL_LIGHTS;
		else m_strategy = DirectSamplingStrategy::SAMPLE_ONE_LIGHT;

	}

	~PathIntegratorMis() {}

	// Estimate Direct Lighting using MIS
	// return appropriately weighted terms.
	Color3f LiDirect(const Scene* scene, Sampler* sampler, const Ray3f& ray, const Intersection& isect) const
	{
		Color3f L_ems(0.0f), L_mats(0.0f);		
		const BSDF* bsdf = isect.mesh->getBSDF();

		// Choose a light
		float pdf = 1.0f / scene->getLights().size();
		const Emitter* random_emitter = scene->getRandomEmitter(sampler->next1D());

		// Emitter Sampling
		// Perform only if not a delta bsdf		
		if(!bsdf->isDelta())
		{
			EmitterQueryRecord eRec;
			float pdf_e, pdf_m;
			eRec.ref = isect.p;

			Color3f Li = random_emitter->sample(eRec, sampler->next2D());
			pdf_e = eRec.pdf;
			
			BSDFQueryRecord bRec(isect.toLocal(-ray.d), isect.toLocal(eRec.wi), ESolidAngle);
			Color3f f = bsdf->eval(bRec);
			pdf_m = bsdf->pdf(bRec);
			if (pdf_e != 0.0f)
			{
				float mis = pdf_e / (pdf_m + pdf_e);

				// Compute lighting
				L_ems = f * Li * fabsf(isect.shFrame.n.dot(eRec.wi)) / pdf;

				// Compute shadow ray only when 
				if (L_ems.isValid() && !L_ems.isZero())
				{
					// Trace a shadow ray only now
					float V = scene->rayIntersect(Ray3f(isect.p, eRec.wi, Epsilon, (1.0f - Epsilon) * eRec.dist)) ? 0.0f : 1.0f;
					L_ems *= V;
				}
				else
					L_ems = Color3f(0.0f);


				if (!random_emitter->isDelta())
				{
					// The BSDF has no way of generating a direction that would hit this light
					// Hence multiply by MIS value only when 
					L_ems *= mis;
				}				
			}			
		}
				
		// BSDF sampling
		// If there were only one light and that light was a delta light, we can't do this kind of sampling.
		if (!(scene->getLights().size() == 1 && scene->getLights()[0]->isDelta()))
		{
			
			BSDFQueryRecord bRec(isect.toLocal(-ray.d));
			Color3f f = bsdf->sample(bRec, sampler->next2D());
			float pdf_m = bRec.pdf;

			if (!f.isZero() && pdf_m != 0.0f && !isnan(pdf_m))
			{
				Intersection light_isect;
				if (scene->rayIntersect(Ray3f(isect.p, isect.toWorld(bRec.wo), Epsilon, INFINITY), light_isect))
				{
					// check if a light soruce
					if (light_isect.mesh->isEmitter())
					{
						const Emitter* light = light_isect.mesh->getEmitter();

						EmitterQueryRecord eRec;
						eRec.ref = isect.p;
						eRec.wi = isect.toWorld(bRec.wo);
						eRec.n = light_isect.shFrame.n;
						eRec.emitter = light;
						eRec.p = light_isect.p;
						eRec.dist = light_isect.t;

						Color3f Li = light->eval(eRec);
						float mis = 1.0f;

						// If the BSDF was delta, the light could have never generated this reverse direction
						if (!bsdf->isDelta())
						{
							float pdf_e = light->pdf(eRec);
							mis = pdf_m / (pdf_m + pdf_e);
						}

						L_mats += (f * Li * fabsf(Frame::cosTheta(bRec.wo)));// *mis;
					}
				}
			}
		}


		return L_ems + L_mats;
	}


	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
	{
		Intersection isect;

		Color3f L(0.0f);
		Color3f throughput(1.0f);
		int depth = 0;
		Ray3f traced_ray = ray;
		bool wasLastBounceSpecular = false;

		while (depth < m_maxDepth || m_maxDepth == -1)
		{
			// Check if ray misses the scene
			if (!scene->rayIntersect(traced_ray, isect))
			{
				L += throughput * scene->getBackground(traced_ray);
				break;
			}

			// Check if direct emitter intersection
			// return this only if direct hit or previous hit was from a specular surface
			// because we couldnt have sampled it using a light sampling strategy.
			if (isect.mesh->isEmitter() && (depth == 0 || !wasLastBounceSpecular))
			{
				EmitterQueryRecord eRec;
				eRec.ref = traced_ray.o;
				eRec.wi = traced_ray.d;
				eRec.n = isect.shFrame.n;
				L += throughput * isect.mesh->getEmitter()->eval(eRec);

				// Assume for now we dont bounce off the light sources.				
			}

			const BSDF* bsdf = isect.mesh->getBSDF();			

			// NEE
			L += throughput * LiDirect(scene, sampler, ray, isect);
						
			// Sample a reflection ray
			BSDFQueryRecord bRec(isect.toLocal(-traced_ray.d));
			Color3f f = bsdf->sample(bRec, sampler->next2D());
			Vector3f reflected_dir = isect.toWorld(bRec.wo);

			throughput *= f * fabsf(Frame::cosTheta(bRec.wo));

			// Check if specular bounce
			wasLastBounceSpecular = bsdf->isDelta();

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
		return tfm::format("PathIntegratorMis[\nrrStart = %d\n]", m_rrStart);
	}

private:
	int m_rrStart;				// from which bounce should russian roulette start.
	int m_maxDepth;				// Fixed length cutoff
	DirectSamplingStrategy m_strategy;

};

NORI_REGISTER_CLASS(PathIntegratorMis, "path_mis")
NORI_NAMESPACE_END