#include "distantdisk.h"
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

DistantDisk::DistantDisk(const PropertyList & props)
{
	m_type = EmitterType::EMITTER_DISTANT_DISK;
	m_radiance = props.getColor("radiance");
	m_thetaA = degToRad(props.getFloat("thetaA"));
	m_localToWorld = props.getTransform("toWorld", Transform());
	m_worldToLocal = m_localToWorld.getInverseMatrix();

	// Compute cosThetaMax that will be used for the hemispherical cap sampling
	m_cosThetaMax = cosf(m_thetaA);
}

DistantDisk::~DistantDisk()
{
	// Empty destructor.
}

Color3f DistantDisk::sample(EmitterQueryRecord & lRec, const Point2f & sample, float optional_u) const
{
	// Sample in local coordinate frame
	Vector3f sampled_dir = Warp::squareToUniformSphereCap(sample, m_cosThetaMax);
	float sampled_pdf = Warp::squareToUniformSphereCapPdf(sampled_dir, m_cosThetaMax);
	
	// Convert to world coordinate frame
	lRec.wi = (m_localToWorld * sampled_dir);
	lRec.pdf = sampled_pdf;
	lRec.dist = INFINITY;
	lRec.emitter = this;
	lRec.n = m_localToWorld * Vector3f(0.0f, 0.0f, 1.0f);
	
	// Appropriately scale the radiance and return back.
	if (sampled_pdf > 0.0f) return m_radiance / sampled_pdf;
	else return 0.0f;
}

float DistantDisk::pdf(const EmitterQueryRecord & lRec) const
{
	// Compute the pdf of sampling the direction.
	Vector3f world_dir = lRec.wi;
	Vector3f local_dir = m_worldToLocal * world_dir;
	return Warp::squareToUniformSphereCapPdf(local_dir, m_cosThetaMax);
}

Color3f DistantDisk::eval(const EmitterQueryRecord & lRec) const
{
	// Return radiance only if it's within the accepted limits.
	Vector3f world_dir = lRec.wi;
	Vector3f local_dir = m_worldToLocal * world_dir;
	float cos_theta = Frame::cosTheta(local_dir);
	if (cos_theta >= m_cosThetaMax) return m_radiance;
	else return 0.0f;
}

std::string DistantDisk::toString() const
{
	return tfm::format("DistantDisk[Radiance={%f,%f,%f}\n thetaA:%f]", m_radiance.x(), m_radiance.y(), m_radiance.z(), m_thetaA);
}

NORI_REGISTER_CLASS(DistantDisk, "distantdisk")
NORI_NAMESPACE_END
