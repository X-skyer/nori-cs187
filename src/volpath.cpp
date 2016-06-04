#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/medium.h>
#include <nori/phase.h>

NORI_NAMESPACE_BEGIN

class Volpath : public Integrator
{
public:

	Volpath(const PropertyList& props)
	{
		m_maxDepth = props.getInteger("maxDepth", -1);
		m_rrStart = props.getInteger("rrStart", 5);
	}

	Color3f LiAttenuatedDirect(const Scene* scene, Sampler* sampler, const Ray3f& ray, const Intersection& isect) const
	{
		Color3f L_ems(0.0f), L_mats(0.0f);
		const BSDF* bsdf = isect.mesh->getBSDF();

		// Choose a light
		float pdf = 1.0f / scene->getLights().size();
		const Emitter* random_emitter = scene->getRandomEmitter(sampler->next1D());

		// Emitter Sampling
		// Perform only if not a delta bsdf		
		if (!bsdf->isDelta())
		{
			EmitterQueryRecord eRec;
			float pdf_e, pdf_m;
			eRec.ref = isect.p;

			Color3f Li = random_emitter->sample(eRec, sampler->next2D(), sampler->next1D());
			pdf_e = eRec.pdf;

			BSDFQueryRecord bRec(isect.toLocal(-ray.d), isect.toLocal(eRec.wi), ESolidAngle);
			Color3f f = bsdf->eval(bRec);
			pdf_m = bsdf->pdf(bRec);
			if (pdf_e != 0.0f)
			{
				float mis = pdf_e / (pdf_m + pdf_e);

				// Compute lighting
				L_ems = f * Li * fabsf(isect.shFrame.n.dot(eRec.wi));

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
		// If the chosen light was a delta light, we can't expect the BSDF to sample a direction that will hit the light
		if (!random_emitter->isDelta())
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
					if (light_isect.mesh->isEmitter() && light_isect.mesh->getEmitter() == random_emitter)
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

						L_mats = (f * Li * fabsf(Frame::cosTheta(bRec.wo)));
						L_mats *= mis;
					}
				}
			}
		}

		// Divide by the pdf of choosing the random light
		return (L_ems + L_mats) / pdf;
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
	{
		Color3f L(0.0f);
		int depth = 0;
		Intersection isect;
		MediumSamplingRecord mRec;
				
		Ray3f traced_ray = ray;
		Color3f throughput(1.0f);
		bool wasLastBounceSpecular = true;

		scene->rayIntersect(traced_ray, isect);
		const Medium* m = scene->getSceneMedium();
		depth++;
		
		while (depth < m_maxDepth || m_maxDepth == -1)
		{
			// check for medium interaction or surface interaction
			if(m->sample_distance(traced_ray, mRec, sampler->next2D()))
			{

				throughput *= mRec.transmittance / mRec.pdf_success;
				
				// handle medium event
												
				// handle single bounce of light
				const Emitter* e = scene->getRandomEmitter(sampler->next1D());
				EmitterQueryRecord eRec;
				eRec.ref = mRec.p;
				
				Color3f L_sampled = e->sample(eRec, sampler->next2D(), sampler->next1D());

				// compute phase function
				const PhaseFunction* pfun = mRec.phase_function;
				PhaseFunctionSamplingRecord pRec(mRec, mRec.p, ray.d, eRec.wi);

				float pfun_value = pfun->eval(pRec);

				// Compute transmittance between sampled point and light point
				Color3f transmittance(1.0f);

				// check for visibility
				if (!scene->rayIntersect(Ray3f(mRec.p, eRec.wi, Epsilon, eRec.dist * (1.0f - Epsilon))))
				{
					L += throughput * mRec.m_sigmaS * pfun_value * L_sampled * transmittance;
				}				

				// MIS with phase function

				// sample the phase function for next bounce
				float phase_value = pfun->sample(pRec, sampler->next2D());

				traced_ray = Ray3f(mRec.p, pRec.wo);
				scene->rayIntersect(traced_ray, isect);
				depth++;
			}
			else
			{
				// Surface interaction
				throughput *= mRec.transmittance / mRec.pdf_failure;

				// Check if light source intersected
				if (isect.mesh->isEmitter())
				{
					const Emitter* emitter = isect.mesh->getEmitter();
					EmitterQueryRecord eRec;
					eRec.ref = traced_ray.o;
					eRec.wi = traced_ray.d;
					eRec.n = isect.shFrame.n;
					L += throughput * emitter->eval(eRec);
				}

				const BSDF* bsdf = isect.mesh->getBSDF();

				// NEE
				Color3f Li = LiAttenuatedDirect(scene, sampler, traced_ray, isect);
				Color3f debug = throughput * Li;
				L += throughput * Li;

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
					float rrcond = throughput.getLuminance();
					if (sampler->next1D() > rrcond)
						break;
					else throughput *= rrcond;
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
		}

		return L;
	}

	std::string toString() const
	{
		return tfm::format("Volpath Integrator["
			"Max Depth : %d\n"
			"Russian Roulette Start : %d\n]",
			m_maxDepth,
			m_rrStart
		);
	}

private:
	int m_maxDepth;
	int m_rrStart;
};

NORI_REGISTER_CLASS(Volpath, "volpath");
NORI_NAMESPACE_END