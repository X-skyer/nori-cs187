#include "distantdisk.h"
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

DistantDisk::DistantDisk(const PropertyList & props)
{
	m_radiance = props.getColor("radiance");
	m_thetaA = degToRad(props.getFloat("thetaA"));
	m_localToWorld = props.getTransform("toWorld", Transform());
	m_worldToLocal = m_localToWorld.getInverseMatrix();

	// Compute cosThetaMax that will be used for the hemispherical cap sampling
	m_cosThetaMax = cosf(m_thetaA);
}

DistantDisk::~DistantDisk()
{
}

Color3f DistantDisk::sample(EmitterQueryRecord & lRec, const Point2f & sample) const
{
	// Sample in local coordinate frame
	Vector3f sampled_dir = Warp::squareToUniformSphereCap(sample, m_cosThetaMax);
	float sampled_pdf = Warp::squareToUniformSphereCapPdf(sampled_dir, m_cosThetaMax);

	// Convert to world coordinate frame
	lRec.wi = -(m_localToWorld * sampled_dir);
	lRec.pdf = sampled_pdf;
	lRec.dist = INFINITY;
	lRec.emitter = this;
	lRec.n = m_localToWorld * Vector3f(0.0f, 0.0f, 1.0f);
	return m_radiance / sampled_pdf;
}

float DistantDisk::pdf(const EmitterQueryRecord & lRec) const
{
	Vector3f world_dir = -lRec.wi;
	Vector3f local_dir = m_worldToLocal * world_dir;
	return Warp::squareToUniformSphereCapPdf(local_dir, m_cosThetaMax);
}

Color3f DistantDisk::eval(const EmitterQueryRecord & lRec) const
{
	Vector3f world_dir = -lRec.wi;
	Vector3f local_dir = m_worldToLocal * world_dir;
	if (Frame::cosTheta(local_dir) >= m_cosThetaMax) return m_radiance;
	else return 0.0f;
}

std::string DistantDisk::toString() const
{
	return tfm::format("DistantDisk[Radiance={%f,%f,%f}\n thetaA:%f]", m_radiance.x(), m_radiance.y(), m_radiance.z(), m_thetaA);
}

NORI_REGISTER_CLASS(DistantDisk, "distantdisk")
NORI_NAMESPACE_END