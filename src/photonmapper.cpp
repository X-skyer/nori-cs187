/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/scene.h>
#include <nori/photon.h>

NORI_NAMESPACE_BEGIN

class PhotonMapper : public Integrator {
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;

    PhotonMapper(const PropertyList &props) {
        /* Lookup parameters */
        m_photonCount  = props.getInteger("photonCount", 1000);
        m_photonRadius = props.getFloat("photonRadius", 0.0f /* Default: automatic */);
		m_rrStart = props.getInteger("rrStart", 5);
		m_maxDepth = props.getInteger("maxDepth", -1);
    }

    virtual void preprocess(const Scene *scene) {
        cout << "Gathering " << m_photonCount << " photons .. ";
        cout.flush();

        /* Create a sample generator for the preprocess step */
        Sampler *sampler = static_cast<Sampler *>(
            NoriObjectFactory::createInstance("independent", PropertyList()));

        /* Allocate memory for the photon map */
        m_photonMap = std::unique_ptr<PhotonMap>(new PhotonMap());
        m_photonMap->reserve(m_photonCount);

		/* Estimate a default photon radius */
		if (m_photonRadius == 0)
			m_photonRadius = scene->getBoundingBox().getExtents().norm() / 500.0f;

		bool bAreaLightPresent = false;
		for (auto e : scene->getLights())
		{
			if (e->getEmitterType() == EmitterType::EMITTER_AREA)
			{
				bAreaLightPresent = true;
				break;
			}
		}

		if (!bAreaLightPresent)
			throw NoriException("No Area lights present to shoot photons");
	

        /* How to add a photon?
	 * m_photonMap->push_back(Photon(
	 *	Point3f(0, 0, 0),  // Position
	 *	Vector3f(0, 0, 1), // Direction
	 *	Color3f(1, 2, 3)   // Power
	 * ));
	 */

	// put your code to trace photons here
		int stored_photons = 0;
		int emitted_photons = 0;
		int n_lights = scene->getLights().size();

		while (stored_photons < m_photonCount)
		{
			// First choose a light
			const Emitter* emitter = scene->getRandomEmitter(sampler->next1D());
			
			if (emitter->getEmitterType() == EmitterType::EMITTER_AREA)
			{
				Ray3f photon_ray;
				Color3f photon_power = emitter->samplePhoton(photon_ray, sampler->next2D(), sampler->next2D(), sampler->next1D());
				if (!photon_power.isZero())
				{
					emitted_photons++;			// keep track of how many photons we shot to divide the contrib of all stored photons finally

					// start the trace of the photon
					Intersection isect;
					int depth = 0;
					while ((depth < m_maxDepth || m_maxDepth == -1) && stored_photons < m_photonCount)
					{
						// Check for intersection
						// If not intersection, don't do anything
						if (!scene->rayIntersect(photon_ray, isect))
						{
							// This photon was waste. So we are not going to use it at all.
							// So negate it.
							//if(depth == 0)
							//	emitted_photons--;

							// else break
							break;
						}

						const BSDF* bsdf = isect.mesh->getBSDF();
						
						if (bsdf->isDiffuse())
						{
							// Store photon
							m_photonMap->push_back(Photon(isect.p, photon_ray.d, photon_power));
							stored_photons++;
						}

						// Now sample next direction
						BSDFQueryRecord bRec(isect.toLocal(-photon_ray.d));
						Color3f f = bsdf->sample(bRec, sampler->next2D());
						Vector3f reflected_dir = isect.toWorld(bRec.wo);

						//update the photon power
						// pdf is included in the f term
						Color3f incoming_power = photon_power;
						photon_power *= f * fabsf(Frame::cosTheta(bRec.wo));

						// Check for zero bsdf
						if (f.isZero())
							break;

						// Check for russian roulette
						if (depth > m_rrStart)
						{
							float p = 1.0f - photon_power.getLuminance() / incoming_power.getLuminance();
							if (sampler->next1D() < p)
								break;
							else photon_power /= (1.0f - p);
						}
						if (depth > m_maxDepth && m_maxDepth != -1)
						{
							break;
						}

						// Generate the next bounce direction
						photon_ray = Ray3f(isect.p, reflected_dir, Epsilon, INFINITY);
						depth++;
					}
				}
			}
		}

		// Divide all photons by the total emitted
		m_photonMap->scale(emitted_photons);

		/* Build the photon map */
        m_photonMap->build();
    }

    virtual Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const {
    	
	/* How to find photons?
	 * std::vector<uint32_t> results;
	 * m_photonMap->search(Point3f(0, 0, 0), // lookup position
         *                     m_photonRadius,   // search radius
         *                     results);
         *
	 * for (uint32_t i : results) {
         *    const Photon &photon = (*m_photonMap)[i];
	 *    cout << "Found photon!" << endl;
	 *    cout << " Position  : " << photon.getPosition().toString() << endl;
	 *    cout << " Power     : " << photon.getPower().toString() << endl;
	 *    cout << " Direction : " << photon.getDirection().toString() << endl;
	 * }
	 */

	// put your code for path tracing with photon gathering here

		Color3f L(0.0f);
		Color3f throughput(1.0f);
		int depth = 0;
		Intersection isect;
		bool bLastBounceWasSpecular = false;
		Ray3f traced_ray = _ray;

		while (depth < m_maxDepth || m_maxDepth == -1)
		{
			if (!scene->rayIntersect(traced_ray, isect))
				return scene->getBackground(traced_ray);

			if (isect.mesh->isEmitter())
			{
				EmitterQueryRecord eRec;
				eRec.ref = traced_ray.o;
				eRec.wi = traced_ray.d;
				eRec.n = isect.shFrame.n;
				L += throughput * isect.mesh->getEmitter()->eval(eRec);

				// Assume for now we dont bounce off the light sources.
				break;
			}

			const BSDF* bsdf = isect.mesh->getBSDF();

			// Compute indirect contribution only from diffuse surfaces
			if (bsdf->isDiffuse())
			{
				std::vector<uint32_t> results;
				m_photonMap->search(isect.p, m_photonRadius, results);
				float area = M_PI * square(m_photonRadius);

				// The uint32_t makes the size() - 1 wrap around. Subtle bug.
				if (results.size() > 0)
				{
					for (uint32_t i = 0; i < results.size() - 1; i++)
					{
						const Photon &photon = (*m_photonMap)[results[i]];

						// Compute the integral equation
						BSDFQueryRecord bRec(isect.toLocal(-traced_ray.d), isect.toLocal(-photon.getDirection()), ESolidAngle);

						L += throughput * bsdf->eval(bRec) * photon.getPower() / area;
					}
				}
				// Immediately stop recursion
				break;
			}

			// Sample a reflection ray
			BSDFQueryRecord bRec(isect.toLocal(-traced_ray.d));
			Color3f f = bsdf->sample(bRec, sampler->next2D());
			Vector3f reflected_dir = isect.toWorld(bRec.wo);

			float cos_theta = fabsf(Frame::cosTheta(bRec.wo));
			throughput *= f * cos_theta;

			// Check if we've fa
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
    
		return L;
    }

    virtual std::string toString() const {
        return tfm::format(
            "PhotonMapper[\n"
            "  photonCount = %i,\n"
            "  photonRadius = %f\n"
            "]",
            m_photonCount,
            m_photonRadius
        );
    }
private:
    int m_photonCount;
    float m_photonRadius;
	int m_rrStart;
	int m_maxDepth;
    std::unique_ptr<PhotonMap> m_photonMap;
};

NORI_REGISTER_CLASS(PhotonMapper, "photonmapper");
NORI_NAMESPACE_END
