#include <nori/medium.h>
#include <nori/phase.h>

NORI_NAMESPACE_BEGIN

Medium::~Medium()
{
	// THe parent is responsible for deleting the child
	delete m_phase_funtion;
}

void Medium::addChild(NoriObject* obj)
{
	if (obj->getClassType() == EPhaseFunction)
	{
		m_phase_funtion = static_cast<PhaseFunction*>(obj);
	}
}

NORI_NAMESPACE_END