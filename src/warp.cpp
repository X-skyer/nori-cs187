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

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Vector3f Warp::sampleUniformHemisphere(Sampler *sampler, const Normal3f &pole) {
    // Naive implementation using rejection sampling
    Vector3f v;
    do {
        v.x() = 1.f - 2.f * sampler->next1D();
        v.y() = 1.f - 2.f * sampler->next1D();
        v.z() = 1.f - 2.f * sampler->next1D();
    } while (v.squaredNorm() > 1.f);

    if (v.dot(pole) < 0.f)
        v = -v;
    v /= v.norm();

    return v;
}

/* Square */
Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

/* Disk */
Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float r = sqrtf(sample.x());
    float theta = 2 * M_PI * sample.y();
    float x = r * cos(theta);
    float y = r * sin(theta);
    return Point2f(x, y);
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    if(p.x() * p.x() + p.y() * p.y() < 1.0f) return INV_PI;
    else return 0.0f;
}

/* Sphere Cap */
Vector3f Warp::squareToUniformSphereCap(const Point2f &sample, float cosThetaMax) {
    throw NoriException("Warp::squareToUniformSphereCap() is not yet implemented!");
}

float Warp::squareToUniformSphereCapPdf(const Vector3f &v, float cosThetaMax) {
    throw NoriException("Warp::squareToUniformSphereCapPdf() is not yet implemented!");
}

/* Uniform Sphere */
Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float theta = acosf(2.0f * sample.x() - 1.0f);
    float phi = 2.0f * M_PI * sample.y();
    return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return INV_FOURPI;
}

/* Uniform Hemisphere */
Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float z = sample.x();
    float r = sqrtf(std::max(0.0f, 1.0f - z * z));
    float phi = 2 * M_PI * sample.y();
    float x = r * cos(phi);
    float y = r * sin(phi);
    return Vector3f(x, y, z);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    return INV_TWOPI;
}

/* CosineHemisphere */
Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    float theta = 0.5f * acosf(1.0f - 2 * sample.x());
    float phi = 2.0f * M_PI * sample.y();
    return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    float cos_theta = Frame::cosTheta(v);
    return cos_theta / M_PI;
}


/* Beckmann */
Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    throw NoriException("Warp::squareToBeckmann() is not yet implemented!");
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    throw NoriException("Warp::squareToBeckmannPdf() is not yet implemented!");
}

NORI_NAMESPACE_END
