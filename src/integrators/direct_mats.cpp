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

		// Else do brdf sampling.
		Color3f Ld(0.0f);
		const BSDF* bsdf = its.mesh->getBSDF();
		//for (auto e : scene->getLights())
		{
			// Construct a BSDF query record
			BSDFQueryRecord bRec(its.toLocal(-ray.d));
			bRec.p = its.p;
			
			Color3f f = bsdf->sample(bRec, sampler->next2D());
			float pdf_f = bsdf->pdf(bRec);
			const Ray3f shadow_ray(its.p, its.toWorld(bRec.wi), Epsilon, INFINITY);
			Intersection s_isect;
			if (scene->rayIntersect(shadow_ray, s_isect))
			{
				// check if the intersected object was an emitter
				if (s_isect.mesh->isEmitter())
				{
					// Construct an emitter query record
					EmitterQueryRecord eRec;
					eRec.wi = bRec.wi;

					// Get the radiance along the intersected direction
					const Emitter* e = s_isect.mesh->getEmitter();
					Color3f Li = e->eval(eRec);

					// Compute the direct lighting equation.
					Color3f evalTerm = f * Li * fmaxf(its.shFrame.n.dot(eRec.wi), 0.0f) / pdf_f;
					if (!evalTerm.isValid())
					{
						std::cout << "Invalid term" << std::endl;
					}
					Ld += evalTerm;
				}
			}
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