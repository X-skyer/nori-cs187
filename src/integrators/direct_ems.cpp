#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <fstream>

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
		// Else do light sampling.
		Color3f Ld(0.0f);

		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return scene->getBackground(ray);

		/* Intersection found */
		// Check for intersection if a direct light source.
		if (its.mesh->isEmitter())
		{
			EmitterQueryRecord eRec;
			eRec.p = its.p;
			eRec.ref = ray.o;
			eRec.wi = ray.d;
			eRec.n = its.shFrame.n;
			const Emitter* e = its.mesh->getEmitter();
			Ld += e->eval(eRec);
		}
		
		const BSDF* bsdf = its.mesh->getBSDF();
		for (auto e : scene->getLights())
		{
			// Construct an Emitter query record
			EmitterQueryRecord eRec;
			eRec.ref = its.p;
			
			// Get the incoming radiance and create shadow ray.
			// Assume Li has the pdf included in it.
			//std::ofstream park_sampled_pts("park_sampled.csv");

			Color3f Li = e->sample(eRec, sampler->next2D(), sampler->next1D());

			/*
			for (int i = 0; i < 512; i++)
			{
				Li = e->sample(eRec, sampler->next2D(), sampler->next1D());
				park_sampled_pts << eRec.p.x() << "," << eRec.p.y() << std::endl;
			}
			*/			
			
			// Compute BSDF contribution for chosen direction.
			BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(eRec.wi), ESolidAngle);
			bRec.p = its.p;
			bRec.uv = its.uv;

			Color3f evalTerm = bsdf->eval(bRec) * Li * fabsf(its.shFrame.n.dot(eRec.wi));
			
			if (!evalTerm.isZero() && evalTerm.isValid())
			{
				const Ray3f shadow_ray(its.p, eRec.wi, Epsilon, (1.0f - Epsilon) * eRec.dist);
				Intersection s_isect;
				if (!scene->rayIntersect(shadow_ray, s_isect))
				{
					Ld += evalTerm;
				}
			}
		}

		return Ld;
	}
	
	std::string toString() const
	{
		return "DirectEmsIntegrator[]";
	}

private:

};


NORI_REGISTER_CLASS(DirectEmsIntegrator, "direct_ems");
NORI_NAMESPACE_END