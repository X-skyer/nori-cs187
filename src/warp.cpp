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
    if(p.x() * p.x() + p.y() * p.y() <= 1.0f) return INV_PI;
    else return 0.0f;
}

/* Uniform Sphere */
Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float theta = acosf(1.0f - 2 * sample.x());
    float phi = 2.0f * M_PI * sample.y();
    return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    return INV_FOURPI;
}


/* Uniform Hemisphere */
Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float theta = acosf(sample.x());
    float phi = 2 * M_PI * sample.y();
    return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    float cos_theta = Frame::cosTheta(v);
    if(cos_theta < 0.0f) return 0.0f;
    return 1.0f / (2.0f * M_PI);
}


/* Sphere Cap */
Vector3f Warp::squareToUniformSphereCap(const Point2f &sample, float cosThetaMax) {
    float theta = acosf(1.0f - (1.0f - cosThetaMax) * sample.x());
    float phi = 2.0f * M_PI * sample.y();
    return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

float Warp::squareToUniformSphereCapPdf(const Vector3f &v, float cosThetaMax) {
    float cos_theta = Frame::cosTheta(v);
    if(cos_theta < cosThetaMax) return 0.0f;
    return INV_TWOPI * (1.0f / (1.0f - cosThetaMax));
}


/* CosineHemisphere */
Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    float theta = 0.5f * acosf(1.0f - 2 * sample.x());
    float phi = 2.0f * M_PI * sample.y();
    return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    float cos_theta = Frame::cosTheta(v);
    if(cos_theta < 0.0f) return 0.0f;
    return cos_theta / M_PI;
}


/* Beckmann */
Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    float theta = atanf(sqrtf(-(alpha * alpha * log(sample.x()))));
    float phi = 2.0f * M_PI * sample.y();
    return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    float cos_theta = Frame::cosTheta(m);
    if(cos_theta <= 0.0f) return 0.0f;
    
    float theta = acosf(cos_theta);
    float cos_theta2 = cos_theta * cos_theta;
    float cos_theta4 = cos_theta2 * cos_theta2;
    float tan_theta = tan(theta);
    float tan_theta2 = tan_theta * tan_theta;
    float alpha2 = alpha * alpha;
    
    float D_theta = (expf(-tan_theta2/alpha2))/(M_PI * alpha2 * cos_theta4);
    return D_theta * cos_theta;
}

/* Ggx */
Vector3f Warp::squareToGgx(const Point2f& sample, float alpha_g)
{
	float theta = atanf(alpha_g * sqrtf(sample.x()) / sqrtf(1 - sample.x()));
	float phi = 2.0f * M_PI * sample.y();
	return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

float Warp::squareToGgxPdf(const Vector3f& m, float alpha_g)
{
	float cos_theta = Frame::cosTheta(m);
	if (cos_theta <= 0.0f) return 0.0f;
	
	float theta = acosf(cos_theta);
	float cos_theta2 = cos_theta * cos_theta;
	float cos_theta4 = cos_theta2 * cos_theta2;
	float tan_theta = tan(theta);
	float tan_theta2 = tan_theta * tan_theta;

	float alpha2 = alpha_g * alpha_g;
	float term = alpha2 + tan_theta2;
	float term2 = term * term;

	float D = alpha2 / (M_PI * cos_theta4 * term2);
	return D * cos_theta;
}

/* Phong */
Vector3f Warp::squareToPhong(const Point2f& sample, float alpha_p)
{
	float theta = acosf(std::pow(sample.x(), 1.0f / (alpha_p + 2.0f)));
	float phi = 2 * M_PI * sample.y();
	return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

float Warp::squareToPhongPdf(const Vector3f& m, float alpha_p)
{
	float cos_theta = Frame::cosTheta(m);
	if (cos_theta <= 0.0f) return 0.0f;

	return (alpha_p + 2.0f) * INV_TWOPI * std::pow(cos_theta, alpha_p) * cos_theta;
}

NORI_NAMESPACE_END
