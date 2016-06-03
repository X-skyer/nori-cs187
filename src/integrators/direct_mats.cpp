#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectMatsIntegrator : public Integrator
{
public:
	DirectMatsIntegrator(const PropertyList& prop)
	{
		// Empty constructor for now.
	}

	~DirectMatsIntegrator()
	{
		// Empty destructor for now.	
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const
	{
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return scene->getBackground(ray);

		/* Intersection found */
		Color3f Ld(0.0f);

		// Check for intersection if a direct light source.
		if (its.mesh->isEmitter())
		{
			EmitterQueryRecord eRec;
			eRec.ref = ray.o;
			eRec.wi = ray.d;
			eRec.n = its.shFrame.n;
			const Emitter* e = its.mesh->getEmitter();
			Ld += e->eval(eRec);
		}

		// Do brdf sampling.		
		const BSDF* bsdf = its.mesh->getBSDF();
		
		// Construct a BSDF query record
		BSDFQueryRecord bRec(its.toLocal(-ray.d));
		bRec.p = its.p;
		bRec.uv = its.uv;
		
		Color3f f = bsdf->sample(bRec, sampler->next2D(), sampler->next1D());
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
				eRec.n = s_isect.shFrame.n;
				eRec.p = s_isect.p;
				eRec.dist = (eRec.p - eRec.ref).norm();

				// Get the radiance along the intersected direction
				const Emitter* e = s_isect.mesh->getEmitter();
				Color3f Li = e->eval(eRec);

				// Compute the direct lighting equation.
				Color3f evalTerm = f * Li * fmaxf(its.shFrame.n.dot(eRec.wi), 0.0f);
				Ld += evalTerm;
			}
		}
		else
		{
			// Check if light is directional?
			Color3f Li = scene->getBackground(shadow_ray);
			Ld += f * Li * fmaxf(its.shFrame.n.dot(shadow_ray.d), 0.0f);
		}

		return Ld;
	}

	std::string toString() const
	{
		return "DirectMatsIntegrator[]";
	}

private:

};


NORI_REGISTER_CLASS(DirectMatsIntegrator, "direct_mats");
NORI_NAMESPACE_END