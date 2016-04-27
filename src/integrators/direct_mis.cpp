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


		/* 
			Valid intersection found.
			Sample both emitter and bsdf strategies
		*/

		// If intersecting a direct light source, return radiance immediately.
		if (its.mesh->isEmitter())
		{
			EmitterQueryRecord eRec;
			eRec.ref = ray.o;
			eRec.wi = ray.d;
			const Emitter* e = its.mesh->getEmitter();
			return e->eval(eRec);
		}

		Color3f L_bsdf(0.0f), L_emitter(0.0f);
		float w_bsdf = 0.0f, w_emitter = 0.0f;				

		{
			// BSDF sampling

		}

		{
			// Emitter sampling
		}

		return L_bsdf * w_bsdf + L_emitter * w_emitter;
	}

	std::string toString() const
	{
		return "DirectMisIntegrator[]";
	}

private:

};


NORI_REGISTER_CLASS(DirectMisIntegrator, "direct_mis");
NORI_NAMESPACE_END