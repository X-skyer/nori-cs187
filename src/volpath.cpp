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

	Color3f LiAttenuatedDirect(const Scene* scene, Sampler* sampler, const Ray3f& ray, const Intersection& isect) const
	{
		Color3f L_ems(0.0f), L_mats(0.0f);
		const BSDF* bsdf = isect.mesh->getBSDF();
		float distance_to_light;
		const Medium* m = scene->getSceneMedium();
		Color3f transmittance(1.0);

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
				distance_to_light = eRec.dist;

				// Compute shadow ray only when 
				if (L_ems.isValid() && !L_ems.isZero())
				{
					// Trace a shadow ray only now
					Ray3f shadow_ray(isect.p, eRec.wi, Epsilon, (1.0f - Epsilon) * eRec.dist);
					float V = scene->rayIntersect(shadow_ray) ? 0.0f : 1.0f;
					transmittance = m->eval_transmittance(shadow_ray);

					L_ems *= V * transmittance;
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
				Ray3f next_bounce_ray(isect.p, isect.toWorld(bRec.wo), Epsilon, INFINITY);
				if (scene->rayIntersect(next_bounce_ray, light_isect))
				{
					// check if a light soruce
					if (light_isect.mesh->isEmitter() && light_isect.mesh->getEmitter() == random_emitter)
					{
						next_bounce_ray.maxt = light_isect.t;
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

						transmittance = m->eval_transmittance(next_bounce_ray);
						L_mats = transmittance * (f * Li * fabsf(Frame::cosTheta(bRec.wo)));
						L_mats *= mis;
					}
				}
			}
		}

		// Divide by the pdf of choosing the random light
		return ((L_ems + L_mats) / pdf);
	}

	Color3f LmSingleScatter(const Scene* scene, Sampler* sampler, const Ray3f& ray, const MediumSamplingRecord& mRec, const Point3f& scatter_point) const
	{
		float pdf = 1.0f;

		const Medium* m = scene->getSceneMedium();

		// Compute direct lighting to emitter
		const Emitter* emitter = scene->getRandomEmitter(sampler->next1D());
		float light_pdf = 1.0f / scene->getLights().size();

		Color3f Lm_emit(0.f), Lm_phase(0.0f);
		// Choose a point on the light source
		{
			EmitterQueryRecord eRec;
			float pdf_e;
			eRec.ref = scatter_point;

			Color3f Li = emitter->sample(eRec, sampler->next2D(), sampler->next1D());
			pdf_e = eRec.pdf;

			float phase_fun = INV_FOURPI;
			float pdf_phase = INV_FOURPI;
			if (pdf_e != 0.0f)
			{
				float mis = pdf_e / (pdf_phase + pdf_e);

				// Compute lighting
				Color3f transmittance_term1 = m->eval_transmittance(Ray3f(ray.o, ray.d, ray.mint, mRec.t));
				Color3f transmittance_term2 = m->eval_transmittance(Ray3f(scatter_point, eRec.wi, 0, eRec.dist));
				Lm_emit = Li * transmittance_term1 * transmittance_term2 * phase_fun * m->get_sigmaS();

				// Compute shadow ray only when 
				if (Lm_emit.isValid() && !Lm_emit.isZero())
				{
					// Trace a shadow ray only now
					Ray3f shadow_ray(scatter_point, eRec.wi, Epsilon, (1.0f - Epsilon) * eRec.dist);
					float V = scene->rayIntersect(shadow_ray) ? 0.0f : 1.0f;
					Lm_emit *= V;
				}
				else
					Lm_emit = Color3f(0.0f);

				if (!emitter->isDelta())
				{
					// The phase function has no way of generating a direction that would hit this light
					// Hence multiply by MIS value only when 
					Lm_emit *= mis;
				}
			}
		}

		// importance sample phase function with MIS
		{
			auto pFun = m->m_phase_funtion;
			Vector3f random_direction = Warp::squareToUniformSphere(sampler->next2D());
			float pdf_phase = INV_FOURPI;
			float phase_fun = INV_FOURPI;

			Ray3f scattered_ray(scatter_point, random_direction, Epsilon, INFINITY);
			Intersection isect;
			if (scene->rayIntersect(scattered_ray, isect))
			{
				if (isect.mesh->isEmitter() && isect.mesh->getEmitter() == emitter)
				{
					const Emitter* light = isect.mesh->getEmitter();

					EmitterQueryRecord eRec;
					eRec.wi = scattered_ray.d;
					eRec.ref = scatter_point;
					eRec.dist = isect.t;
					eRec.p = isect.p;
					eRec.n = isect.shFrame.n;
					eRec.emitter = light;

					Color3f Li = light->eval(eRec);
					float pdf_emit = light->pdf(eRec);
					float mis = pdf_phase / (pdf_emit + pdf_phase);

					Color3f transmittance_term1 = m->eval_transmittance(Ray3f(ray.o, ray.d, ray.mint, mRec.t));
					Color3f transmittance_term2 = m->eval_transmittance(Ray3f(scatter_point, eRec.wi, 0, eRec.dist));
					Lm_phase = Li * transmittance_term1 * transmittance_term2 * phase_fun * m->get_sigmaS() / pdf_phase;
					Lm_phase *= mis;
				}
			}

		}

		float total_pdf = pdf * light_pdf;
		if (total_pdf == 0.0f) return 0.0f;
		return (Lm_emit + Lm_phase) / total_pdf;
	}

	Color3f Li(const Scene* scene, Sampler* sampler, const Ray3f& ray) const
	{
		Color3f L(0.0f);
		Color3f throughput(1.0f);
		
		auto m = scene->getSceneMedium();
		Intersection isect;
		Ray3f traced_ray = ray;
		int depth = 0;

		MediumSamplingRecord mRec;

		// compute intersection and update the valid range for the first ray
		if (scene->rayIntersect(traced_ray, isect))
			traced_ray.maxt = isect.t;
		else
			return L;			// for now dont trace these rays.
		
		while (depth < m_maxDepth || m_maxDepth == -1)
		{
			if (m->sample_distance(traced_ray, mRec, sampler->next2D()))
			{
				// Compute scattered point
				Point3f scattered_pt = traced_ray(mRec.t);

				// Compute scattering term
				// SingleScatter includes all transmittance terms for this scattering event.
				// So don't bother adding it here
				L += throughput * mRec.m_sigmaS * LmSingleScatter(scene, sampler, traced_ray, mRec, scattered_pt) / mRec.pdf_success;

				// update throughput
				// however throughput has to be updated correctly.
				throughput *= mRec.transmittance * mRec.m_sigmaS / mRec.pdf_success;

				// sample a next direction
				Vector3f scatter_dir;
				PhaseFunctionSamplingRecord pRec(mRec, scattered_pt, -ray.d, scatter_dir);
				float phase_fn = m->m_phase_funtion->sample(pRec, sampler->next2D());

				traced_ray = Ray3f(scattered_pt, pRec.wo);
				depth++;
			}
			else
			{
				// transmitted branch from surface
				throughput *= mRec.transmittance / mRec.pdf_failure;

				// check if emitter
				if (isect.mesh->isEmitter() && (depth == 0))
				{

					EmitterQueryRecord eRec;
					eRec.ref = traced_ray.o;
					eRec.wi = traced_ray.d;
					eRec.n = isect.shFrame.n;
					L += throughput * mRec.transmittance * isect.mesh->getEmitter()->eval(eRec);
				}

				const BSDF* bsdf = isect.mesh->getBSDF();

				// NEE
				Color3f Li = LiAttenuatedDirect(scene, sampler, traced_ray, isect);
				Color3f debug = throughput * mRec.transmittance * Li;
				L += throughput * Li;

				// Sample a reflection ray
				BSDFQueryRecord bRec(isect.toLocal(-traced_ray.d));
				Color3f f = bsdf->sample(bRec, sampler->next2D());
				Vector3f reflected_dir = isect.toWorld(bRec.wo);

				throughput *= f * fabsf(Frame::cosTheta(bRec.wo));

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

			// update the ray params
			if (scene->rayIntersect(traced_ray, isect))
				traced_ray.maxt = isect.t;			
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