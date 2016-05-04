#include <nori/bsdf.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

class PathIntegratorMis : public Integrator
{
public:

private:
	int m_rrStart;				// from which bounce should russian roulette start.

};

NORI_NAMESPACE_END