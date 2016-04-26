#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>

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
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
		{
			for (auto e : scene->getLights())
			{
				if (e->getEmitterType() == EmitterType::EMITTER_DISTANT_DISK)
				{
					// get the distant light as of now and return the radiance for the direction
					EmitterQueryRecord eRec;
					eRec.wi = ray.d;
					return e->eval(eRec);
				}
			}
			return Color3f(0.0f);
		}

		/* Intersection found */
		// Check for intersection if a direct light source.
		if (its.mesh->isEmitter())
		{
			EmitterQueryRecord eRec;
			eRec.ref = ray.o;
			eRec.wi = ray.d;
			const Emitter* e = its.mesh->getEmitter();
			return e->eval(eRec);
		}

		// Else do light sampling.
		Color3f Ld(0.0f);
		const BSDF* bsdf = its.mesh->getBSDF();
		for (auto e : scene->getLights())
		{
			// Construct an Emitter query record
			EmitterQueryRecord eRec;
			eRec.ref = its.p;
			
			// Get the incoming radiance and create shadow ray.
			// Assume Li has the pdf included in it.
			Color3f Li = e->sample(eRec, sampler->next2D());
			const Ray3f shadow_ray(its.p, eRec.wi, Epsilon, (1.0f - Epsilon) * eRec.dist);
			Intersection s_isect;
			if (!scene->rayIntersect(shadow_ray, s_isect))
			{
				// If unoccluded to the light source, compute the lighting term and add contributions.
				BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(eRec.wi), ESolidAngle);
				Color3f evalTerm = bsdf->eval(bRec) * Li * fmaxf(its.shFrame.n.dot(eRec.wi), 0.0f);
				if (!evalTerm.isValid())
				{
					std::cout << "Invalid term" << std::endl;
				}
				Ld += evalTerm;
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