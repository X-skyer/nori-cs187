#include <nori/pointlight.h>

NORI_NAMESPACE_BEGIN

PointLight::PointLight(const PropertyList& prop)
{
	m_power = prop.getColor("power");
	m_position = prop.getPoint3("position");
	m_type = EmitterType::EMITTER_POINT;
}

Color3f PointLight::sample(EmitterQueryRecord &lRec, const Point2f &sample, float optional_u) const
{
	// Explicit sampling of the point light.
	lRec.p = m_position;
	lRec.pdf = 1.0f;
	lRec.wi = (m_position - lRec.ref).normalized();
	lRec.dist = (m_position - lRec.ref).norm();
	lRec.emitter = this;
	return m_power * INV_FOURPI / (lRec.dist * lRec.dist);
}

float PointLight::pdf(const EmitterQueryRecord &lRec) const
{
	// We can never sample the pointlight through a random process without explicit connection.
	return 0.0f;
}

Color3f PointLight::eval(const EmitterQueryRecord &lRec) const
{
	return Color3f(0.0f);
}

std::string PointLight::toString() const
{
	return tfm::format("PointLight[Position={%f,%f,%f}; Power={%f,%f,%f}]", m_position.x(), m_position.y(), m_position.z(), m_power.x(), m_power.y(), m_power.z());
}

NORI_REGISTER_CLASS(PointLight, "point");
NORI_NAMESPACE_END