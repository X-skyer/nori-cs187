#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AVIntegrator : public Integrator {
public:
	AVIntegrator(const PropertyList &props) {
		m_length = props.getFloat("length", 1.0f);
	}

	Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
		/* Find the surface that is visible in the requested direction */
		Intersection its;
		if (!scene->rayIntersect(ray, its))
			return Color3f(1.0f);


		/*
			Compose a random reflected ray from the intersection point in the upper hemisphere.
			The max length of this ray is provided by the input parameter to the integrator.
		*/

		Vector3f rand_dir = Warp::sampleUniformHemisphere(sampler, its.shFrame.n);
        Vector3f wrand_dir = its.toWorld(rand_dir);
        
		Ray3f refl_ray(its.p, wrand_dir, Epsilon, m_length);
		Intersection av_isect;


		/*
			Return Color3f(1.0f) if no intersection is found within the max distance.
			Else return Color3f(0.0f) if occluded.
		*/

		if (scene->rayIntersect(refl_ray, av_isect))
			return Color3f(0.0f);
		else return Color3f(1.0f);
	}

	std::string toString() const {
		return tfm::format("AVIntegrator[length=%i]", m_length);
	}

private:
	float m_length;						// maximum length of the average visibility ray shot from intersection point.
};

NORI_REGISTER_CLASS(AVIntegrator, "av");
NORI_NAMESPACE_END