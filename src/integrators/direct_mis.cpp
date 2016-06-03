#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectMisIntegrator : public Integrator
{
public:
	DirectMisIntegrator(const PropertyList& prop)
	{
		// Empty constructor for now.
	}

	~DirectMisIntegrator()
	{
		// Empty destructor for now.	
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
	{
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return scene->getBackground(ray);

		/* 
			Valid intersection found.
			Sample both emitter and bsdf strategies
		*/

		// If intersecting a direct light source, return radiance immediately.
		Color3f Le(0.0f);
		if (its.mesh->isEmitter())
		{
			EmitterQueryRecord eRec;
			eRec.ref = ray.o;
			eRec.wi = ray.d;
			eRec.n = its.geoFrame.n;
			const Emitter* e = its.mesh->getEmitter();
			Le += e->eval(eRec);
		}
		
		Color3f L_bsdf(0.0f), L_emitter(0.0f);
		
		{
			// BSDF sampling
			Color3f Ld(0.0f);
			const BSDF* bsdf = its.mesh->getBSDF();

			// Construct a BSDF query record
			BSDFQueryRecord bRec(its.toLocal(-ray.d));
			bRec.p = its.p;
			bRec.uv = its.uv;
			
			Color3f f = bsdf->sample(bRec, sampler->next2D(), sampler->next1D());
			float bpdf = bsdf->pdf(bRec);
			const Ray3f shadow_ray(its.p, its.toWorld(bRec.wo), Epsilon, INFINITY);
			Intersection s_isect;
			if (scene->rayIntersect(shadow_ray, s_isect))
			{
				// check if the intersected object was an emitter
				if (s_isect.mesh->isEmitter())
				{
					// Construct an emitter query record
					EmitterQueryRecord eRec;
					eRec.ref = shadow_ray.o;
					eRec.emitter = s_isect.mesh->getEmitter();
					eRec.wi = its.toWorld(bRec.wo);
					eRec.n = s_isect.geoFrame.n;
					eRec.p = s_isect.p;
					eRec.dist = (eRec.p - eRec.ref).norm();

					// Get the radiance along the intersected direction
					const Emitter* e = s_isect.mesh->getEmitter();
					float lpdf = e->pdf(eRec);
					Color3f Li = e->eval(eRec);

					// Compute the direct lighting equation.
					Color3f evalTerm = f * Li * fmaxf(its.shFrame.n.dot(eRec.wi), 0.0f);
					float mis = bpdf / (bpdf + lpdf);

					if (mis != 0.0f && lpdf != 0.0f && bpdf != 0.0f)
						evalTerm *= mis;

					L_bsdf += evalTerm;
				}
			}
			else
			{
				// Check if light is directional?
				const Emitter* bgEmitter = scene->getBackgroundEmitter();
				if (bgEmitter != nullptr)
				{
					// get the distant light as of now and return the radiance for the direction
					EmitterQueryRecord eRec;
					eRec.wi = shadow_ray.d;
					eRec.ref = shadow_ray.o;
					eRec.emitter = bgEmitter;

					Color3f Li = bgEmitter->eval(eRec);
					float lpdf = bgEmitter->pdf(eRec);
					//std::cout << "BSDF_lpdf : " << lpdf << std::endl;
					Color3f evalTerm = f * Li * fmaxf(its.shFrame.n.dot(eRec.wi), 0.0f);
					float mis = bpdf / (bpdf + lpdf);
					if (!isnan(mis) && mis != 0.0f && lpdf != 0.0f && bpdf != 0.0f)
						evalTerm *= mis;

					L_bsdf += evalTerm;
				}
			}
		}

		
		{
			// Emitter sampling
			const BSDF* bsdf = its.mesh->getBSDF();
			for (auto e : scene->getLights())
			{
				// Construct an Emitter query record
				EmitterQueryRecord eRec;
				eRec.ref = its.p;

				// Get the incoming radiance and create shadow ray.
				// Assume Li has the pdf included in it.
				Color3f Li = e->sample(eRec, sampler->next2D(), sampler->next1D());
				float lpdf = 1.0f; 
				float bpdf = 0.0f;
				lpdf = e->pdf(eRec);
				const Ray3f shadow_ray(its.p, eRec.wi, Epsilon, (1.0f - Epsilon) * eRec.dist);
				Intersection s_isect;
				if (!scene->rayIntersect(shadow_ray, s_isect))
				{
					// If unoccluded to the light source, compute the lighting term and add contributions.
					BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(eRec.wi), ESolidAngle);
					bRec.p = its.p;
					bRec.uv = its.uv;

					Color3f evalTerm = bsdf->eval(bRec) * Li * fmaxf(its.shFrame.n.dot(eRec.wi), 0.0f);
									
					// compute MIS term
					bpdf = bsdf->pdf(bRec);
					float mis = lpdf / (bpdf + lpdf);
					if(mis != 0.0f && lpdf != 0.0f && bpdf != 0.0f)
						evalTerm *= mis;
					L_emitter += evalTerm;
				}
			}
		}
		return Le + L_bsdf + L_emitter;
	}

	std::string toString() const
	{
		return "DirectMisIntegrator[]";
	}

private:

};


NORI_REGISTER_CLASS(DirectMisIntegrator, "direct_mis");
NORI_NAMESPACE_END