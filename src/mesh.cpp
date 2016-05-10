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

#include <nori/mesh.h>
#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

Mesh::Mesh() { 
	m_totalSurfaceArea = 0.0f;
}

Mesh::~Mesh() {
    delete m_bsdf;
    delete m_emitter;
}

void Mesh::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
    }

	// create the pdf
	for (uint32_t i = 0; i < getTriangleCount(); i++)
	{
		float _area = surfaceArea(i);
		m_pdfs.append(_area);
		m_totalSurfaceArea += _area;
	}
	m_pdfs.normalize();
}

void Mesh::samplePosition(const Point2f &sample, Point3f &p, Normal3f &n, float optional_u) const
{
	auto id = m_pdfs.sample(optional_u);
	uint32_t i0 = m_F(0, id), i1 = m_F(1, id), i2 = m_F(2, id);

	// barycentric sampling of triangle.
	float u1 = sqrtf(sample.x());
	float u = 1.0f - u1;
	float v = sample.y() * u1;

	const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);
	if (m_N.size() != 0)
	{
		const Normal3f n0 = m_N.col(i0), n1 = m_N.col(i1), n2 = m_N.col(i2);
		n = (1.0f - u - v) * n0 + u * n1 + v * n2;
	}
	else
	{
		n = (p1 - p0).cross(p2 - p0).normalized();
	}			

	p = (1.0f - u - v) * p0 + u * p1 + v * p2;
	
}

float Mesh::surfaceArea(uint32_t index) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    /* Find vectors for two edges sharing v[0] */
    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = ray.d.cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f)
        return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = ray.o - p0;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = ray.d.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    /* Ray intersects triangle -> compute t */
    t = edge2.dot(qvec) * inv_det;

    return t >= ray.mint && t <= ray.maxt;
}

bool Mesh::rayMeshIntersect(const Ray3f& ray, Intersection& isect) const
{
	float u, v, t;
	float min_t = INFINITY, min_u = 0.0f, min_v = 0.0f;
	int min_id = -1;
	bool intersected = false;
	for (auto i = 0; i < getTriangleCount(); i++)
	{
		if (rayIntersect(i, ray, u, v, t))
		{
			if (t < min_t)
			{
				min_id = i;
				min_t = t;
				min_u = u;
				min_v = v;
				intersected = true;
			}
		}
	}

	if (min_id != -1)
	{
		isect.p = ray(min_t);
		isect.t = min_t;
		isect.mesh = this;

		// compose the normal for a frame
		uint32_t i0 = m_F(0, min_id), i1 = m_F(1, min_id), i2 = m_F(2, min_id);
		const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);
		Normal3f surface_normal = (p1 - p0).cross(p2 - p0).normalized();
		isect.geoFrame = Frame(surface_normal);
		isect.shFrame = isect.geoFrame;
		isect.uv = Vector2f(min_u, min_v);
	}

	return intersected;
}

bool Mesh::rayMeshIntersectP(const Ray3f& ray) const
{
	float u, v, t;
	for (auto i = 0; i < getTriangleCount(); i++)
	{
		if (rayIntersect(i, ray, u, v, t))
			return true;
	}
	return false;
}

BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
    BoundingBox3f result(m_V.col(m_F(0, index)));
    result.expandBy(m_V.col(m_F(1, index)));
    result.expandBy(m_V.col(m_F(2, index)));
    return result;
}

Point3f Mesh::getCentroid(uint32_t index) const {
    return (1.0f / 3.0f) *
        (m_V.col(m_F(0, index)) +
         m_V.col(m_F(1, index)) +
         m_V.col(m_F(2, index)));
}

void Mesh::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EBSDF:
            if (m_bsdf)
                throw NoriException(
                    "Mesh: tried to register multiple BSDF instances!");
            m_bsdf = static_cast<BSDF *>(obj);
            break;

        case EEmitter: {
                Emitter *emitter = static_cast<Emitter *>(obj);
                if (m_emitter)
                    throw NoriException(
                        "Mesh: tried to register multiple Emitter instances!");
                m_emitter = emitter;
            }
            break;

        default:
            throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                classTypeName(obj->getClassType()));
    }
}

std::string Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  name = \"%s\",\n"
        "  vertexCount = %i,\n"
        "  triangleCount = %i,\n"
        "  bsdf = %s,\n"
        "  emitter = %s\n"
        "]",
        m_name,
        m_V.cols(),
        m_F.cols(),
        m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
        m_emitter ? indent(m_emitter->toString()) : std::string("null")
    );
}

float Mesh::pdf() const
{
	return 1.0f / m_totalSurfaceArea;
}

std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format(
        "Intersection[\n"
        "  p = %s,\n"
        "  t = %f,\n"
        "  uv = %s,\n"
        "  shFrame = %s,\n"
        "  geoFrame = %s,\n"
        "  mesh = %s\n"
        "]",
        p.toString(),
        t,
        uv.toString(),
        indent(shFrame.toString()),
        indent(geoFrame.toString()),
        mesh ? mesh->toString() : std::string("null")
    );
}

NORI_NAMESPACE_END
