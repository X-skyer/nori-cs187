#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class DirectLightingIntegrator : public Integrator
{
public:
	DirectLightingIntegrator(const PropertyList& prop) { }

	~DirectLightingIntegrator() { }

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(0.0f);

		/* Intersection found */
		Color3f Ld(0.0f);
		const BSDF* bsdf = its.mesh->getBSDF();

		for (auto e : scene->getLights())
		{
			// Construct an Emitter query record
			EmitterQueryRecord eRec;
			eRec.ref = its.p;
			eRec.n = its.geoFrame.n;
			
			// Get the incoming radiance and create shadow ray.
			Color3f Li = e->sample(eRec, sampler->next2D(), sampler->next1D());
			const Ray3f shadow_ray(its.p, eRec.wi, Epsilon, (1.0f - Epsilon) * eRec.dist);
			Intersection s_isect;
			if (!scene->rayIntersect(shadow_ray, s_isect))
			{
				// If unoccluded to the light source, compute the lighting term and add contributions.
				BSDFQueryRecord bRec(its.toLocal(-ray.d), its.toLocal(eRec.wi), ESolidAngle);
				Ld += bsdf->eval(bRec) * Li * std::max(its.geoFrame.n.dot(eRec.wi), 0.0f);
			}
		}

		return Ld;
	}

	std::string toString() const {
		return "DirectLightingIntegrator[]";
	}

private:

};

NORI_REGISTER_CLASS(DirectLightingIntegrator, "direct");
NORI_NAMESPACE_END