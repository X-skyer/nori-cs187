#include <nori/distributions.h>

NORI_NAMESPACE_BEGIN

Distribution::Distribution(const std::string& dist, float alpha)
{
	if (dist == "beckmann")
		m_distribution = Distributions::BECKMANN;
	else if (dist == "ggx")
		m_distribution = Distributions::GGX;
	else if (dist == "phong")
		m_distribution = Distributions::PHONG;
	else
		m_distribution = Distributions::BECKMANN;

	m_alpha = alpha;
}

Vector3f Distribution::sample(const Vector2f& sample, float& pdf) const
{
	Vector3f ret;
	
	switch (m_distribution)
	{
	case Distributions::BECKMANN:
		ret = Warp::squareToBeckmann(sample, m_alpha);
		pdf = Warp::squareToBeckmannPdf(ret, m_alpha);
		break;
	case Distributions::GGX:
		ret = Warp::squareToGgx(sample, m_alpha);
		pdf = Warp::squareToGgxPdf(ret, m_alpha);
		break;
	case Distributions::PHONG:
		ret = Warp::squareToPhong(sample, m_alpha);
		pdf = Warp::squareToPhongPdf(ret, m_alpha);
		break;
	}
	return ret;
}

float Distribution::pdf(const Vector3f& m) const
{
	float ret = 0.0f;
	switch (m_distribution)
	{
	case Distributions::BECKMANN:
		ret = Warp::squareToBeckmannPdf(m, m_alpha);
	case Distributions::GGX:
		ret = Warp::squareToGgxPdf(m, m_alpha);
	case Distributions::PHONG:
		ret = Warp::squareToPhongPdf(m, m_alpha);
	}
	return ret;
}

float Distribution::G(const Vector3f& i, const Vector3f& o, const Vector3f& m) const
{
	float g = 0.0f;
	switch (m_distribution)
	{
	case Distributions::BECKMANN:
		g = smithBeckmannG1(i, m) * smithBeckmannG1(o, m);
	case Distributions::PHONG:
		g = smithPhongG1(i, m) * smithPhongG1(o, m);
	case Distributions::GGX:
		g = smithGgxG1(i, m) * smithGgxG1(o, m);
	}

	return g;
}

float Distribution::smithBeckmannG1(const Vector3f& v, const Vector3f& m) const
{
	float tanTheta = Frame::tanTheta(v);

	/* Perpendicular incidence -- no shadowing/masking */
	if (tanTheta == 0.0f)
		return 1.0f;

	/* Can't see the back side from the front and vice versa */
	if (m.dot(v) * Frame::cosTheta(v) <= 0)
		return 0.0f;

	float a = 1.0f / (m_alpha * tanTheta);
	if (a >= 1.6f)
		return 1.0f;
	float a2 = a * a;

	/* Use a fast and accurate (<0.35% rel. error) rational
	approximation to the shadowing-masking function */
	return (3.535f * a + 2.181f * a2)
		/ (1.0f + 2.276f * a + 2.577f * a2);
}

float Distribution::smithGgxG1(const Vector3f& v, const Vector3f& m) const
{
	float tanTheta = Frame::tanTheta(v);

	/* Perpendicular incidence -- no shadowing/masking */
	if (tanTheta == 0.0f)
		return 1.0f;

	/* Can't see the back side from the front and vice versa */
	if (m.dot(v) * Frame::cosTheta(v) <= 0)
		return 0.0f;

	float tan_theta = Frame::tanTheta(v);
	return 2.0f / (1.0f + sqrtf(1.0f + m_alpha * m_alpha * tan_theta * tan_theta));
}

float Distribution::smithPhongG1(const Vector3f& v, const Vector3f& m) const
{
	float tanTheta = Frame::tanTheta(v);

	/* Perpendicular incidence -- no shadowing/masking */
	if (tanTheta == 0.0f)
		return 1.0f;

	/* Can't see the back side from the front and vice versa */
	if (m.dot(v) * Frame::cosTheta(v) <= 0)
		return 0.0f;

	float a = sqrtf(0.5f * m_alpha + 1) / Frame::tanTheta(v);

	if (a >= 1.6f)
		return 1.0f;
	float a2 = a * a;

	/* Use a fast and accurate (<0.35% rel. error) rational
	approximation to the shadowing-masking function */
	return (3.535f * a + 2.181f * a2)
		/ (1.0f + 2.276f * a + 2.577f * a2);
}



NORI_NAMESPACE_END