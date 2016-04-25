#include "area.h"

NORI_NAMESPACE_BEGIN

AreaEmitter::AreaEmitter(const PropertyList& props)
{
	m_radiance = props.getColor("radiance");
}

AreaEmitter::~AreaEmitter()
{
	// Empty destructor
}




NORI_REGISTER_CLASS(AreaEmitter, "area");
NORI_NAMESPACE_END