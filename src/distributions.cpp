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

Vector3f Distribution::sample(const Vector2f& _sample, float& pdf) const
{
	Vector3f ret;
	
	switch (m_distribution)
	{
	case Distributions::BECKMANN:
		ret = Warp::squareToBeckmann(_sample, m_alpha);
		pdf = Warp::squareToBeckmannPdf(ret, m_alpha);
		break;
	case Distributions::GGX:
		ret = Warp::squareToGgx(_sample, m_alpha);
		pdf = Warp::squareToGgxPdf(ret, m_alpha);
		break;
	case Distributions::PHONG:
		ret = Warp::squareToPhong(_sample, m_alpha);
		pdf = Warp::squareToPhongPdf(ret, m_alpha);
		break;
	}
	return ret;
}

Vector3f Distribution::sample(const Vector2f& _sample, float& pdf, float _alpha) const
{
	Vector3f ret;

	switch (m_distribution)
	{
	case Distributions::BECKMANN:
		ret = Warp::squareToBeckmann(_sample, _alpha);
		pdf = Warp::squareToBeckmannPdf(ret, _alpha);
		break;
	case Distributions::GGX:
		ret = Warp::squareToGgx(_sample, _alpha);
		pdf = Warp::squareToGgxPdf(ret, _alpha);
		break;
	case Distributions::PHONG:
		ret = Warp::squareToPhong(_sample, _alpha);
		pdf = Warp::squareToPhongPdf(ret, _alpha);
		break;
	}
	return ret;
}

float Distribution::D(const Vector3f& m) const
{
	float cos_theta = Frame::cosTheta(m);
	if (cos_theta <= 0.0f) return 0.0f;

	switch (m_distribution)
	{
		case Distributions::BECKMANN:
		{
			float temp = Frame::tanTheta(m) / m_alpha, ct = Frame::cosTheta(m), ct2 = ct*ct;
			float val = std::exp(-temp*temp) / (M_PI * m_alpha * m_alpha * ct2 * ct2);
			return val;
		}
		case Distributions::GGX:
		{
			float theta = acosf(cos_theta);
			float cos_theta2 = cos_theta * cos_theta;
			float cos_theta4 = cos_theta2 * cos_theta2;
			float tan_theta = tan(theta);
			float tan_theta2 = tan_theta * tan_theta;

			float alpha2 = m_alpha * m_alpha;
			float term = alpha2 + tan_theta2;
			float term2 = term * term;

			return alpha2 / (M_PI * cos_theta4 * term2);
		}
		case Distributions::PHONG:
		{
			return (m_alpha + 2.0f) * INV_TWOPI * std::pow(Frame::cosTheta(m), m_alpha);
		}
	}
	return 0.0f;
}

float Distribution::D(const Vector3f& m, float alpha) const
{
	float cos_theta = Frame::cosTheta(m);
	if (cos_theta <= 0.0f) return 0.0f;

	switch (m_distribution)
	{
	case Distributions::BECKMANN:
	{
		float temp = Frame::tanTheta(m) / alpha, ct = Frame::cosTheta(m), ct2 = ct*ct;
		return std::exp(-temp*temp) / (M_PI * alpha * alpha * ct2 * ct2);
	}
	case Distributions::GGX:
	{
		float theta = acosf(cos_theta);
		float cos_theta2 = cos_theta * cos_theta;
		float cos_theta4 = cos_theta2 * cos_theta2;
		float tan_theta = tan(theta);
		float tan_theta2 = tan_theta * tan_theta;

		float alpha2 = alpha * alpha;
		float term = alpha2 + tan_theta2;
		float term2 = term * term;

		return alpha2 / (M_PI * cos_theta4 * term2);
	}
	case Distributions::PHONG:
	{
		return (alpha + 2.0f) * INV_TWOPI * std::pow(Frame::cosTheta(m), alpha);
	}
	}
	return 0.0f;
}

float Distribution::pdf(const Vector3f& m) const
{
	float ret = 0.0f;
	switch (m_distribution)
	{
	case Distributions::BECKMANN:
		ret = Warp::squareToBeckmannPdf(m, m_alpha);
		break;
	case Distributions::GGX:
		ret = Warp::squareToGgxPdf(m, m_alpha);
		break;
	case Distributions::PHONG:
		ret = Warp::squareToPhongPdf(m, m_alpha);
		break;
	}
	return ret;
}

float Distribution::pdf(const Vector3f& m, float alpha) const
{
	float ret = 0.0f;
	switch (m_distribution)
	{
	case Distributions::BECKMANN:
		ret = Warp::squareToBeckmannPdf(m, alpha);
		break;
	case Distributions::GGX:
		ret = Warp::squareToGgxPdf(m, alpha);
		break;
	case Distributions::PHONG:
		ret = Warp::squareToPhongPdf(m, alpha);
		break;
	}
	return ret;
}

float Distribution::G(const Vector3f& i, const Vector3f& o, const Vector3f& m) const
{
	float g = 0.0f, g1 = 0.0f, g2 = 0.0f;
	switch (m_distribution)
	{
	case Distributions::BECKMANN:
		g1 = smithBeckmannG1(i, m);
		g2 = smithBeckmannG1(o, m);
		g = g1 * g2;
		break;
	case Distributions::PHONG:
		g = smithPhongG1(i, m) * smithPhongG1(o, m);
		break;
	case Distributions::GGX:
		g = smithGgxG1(i, m) * smithGgxG1(o, m);
		break;
	}
	if (g < 0.0f)
	{
		std::cout << "Something in G is wrong" << std::endl;
	}
	return g;
}

float Distribution::G(const Vector3f& i, const Vector3f& o, const Vector3f& m, float alpha) const
{
	float g = 0.0f;
	float g1 = 0.0f, g2 = 0.0f;
	switch (m_distribution)
	{
	case Distributions::BECKMANN:
		g1 = smithBeckmannG1(i, m, alpha);
		g2 = smithBeckmannG1(o, m, alpha);
		g = g1 * g2;
		break;
	case Distributions::PHONG:
		g = smithPhongG1(i, m, alpha) * smithPhongG1(o, m, alpha);
		break;
	case Distributions::GGX:
		g = smithGgxG1(i, m, alpha) * smithGgxG1(o, m, alpha);
		break;
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

	if (m.dot(v) / Frame::cosTheta(v) <= 0.0f) return 0.0f;

	float a = 1.0f / (m_alpha * tanTheta);
	if (a >= 1.6f)
		return 1.0f;
	float a2 = a * a;

	if (isnan(a) || isinf(a))
	{
		std::cout << "Something is fishing" << std::endl;
	}

	/* Use a fast and accurate (<0.35% rel. error) rational
	approximation to the shadowing-masking function */
	float val = (3.535f * a + 2.181f * a2)
		/ (1.0f + 2.276f * a + 2.577f * a2);
	if (val < 0.0f)
	{
		std::cout << "Something is fishing inside deep G" << std::endl;
	}
	return val;
}

float Distribution::smithBeckmannG1(const Vector3f & v, const Vector3f & m, float alpha) const
{
	float tanTheta = Frame::tanTheta(v);

	/* Perpendicular incidence -- no shadowing/masking */
	if (tanTheta == 0.0f)
		return 1.0f;

	/* Can't see the back side from the front and vice versa */
	if (m.dot(v) * Frame::cosTheta(v) <= 0)
		return 0.0f;

	float a = 1.0f / (alpha * tanTheta);
	if (a >= 1.6f)
		return 1.0f;
	float a2 = a * a;

	if (isnan(a) || isinf(a))
	{
		std::cout << "Something is fishing" << std::endl;
	}

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

float Distribution::smithGgxG1(const Vector3f& v, const Vector3f& m, float alpha) const
{
	float tanTheta = Frame::tanTheta(v);

	/* Perpendicular incidence -- no shadowing/masking */
	if (tanTheta == 0.0f)
		return 1.0f;

	/* Can't see the back side from the front and vice versa */
	if (m.dot(v) * Frame::cosTheta(v) <= 0)
		return 0.0f;

	float tan_theta = Frame::tanTheta(v);
	return 2.0f / (1.0f + sqrtf(1.0f + alpha * alpha * tan_theta * tan_theta));
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

float Distribution::smithPhongG1(const Vector3f& v, const Vector3f& m, float alpha) const
{
	float tanTheta = Frame::tanTheta(v);

	/* Perpendicular incidence -- no shadowing/masking */
	if (tanTheta == 0.0f)
		return 1.0f;

	/* Can't see the back side from the front and vice versa */
	if (m.dot(v) * Frame::cosTheta(v) <= 0)
		return 0.0f;

	float a = sqrtf(0.5f * alpha + 1) / Frame::tanTheta(v);

	if (a >= 1.6f)
		return 1.0f;
	float a2 = a * a;

	/* Use a fast and accurate (<0.35% rel. error) rational
	approximation to the shadowing-masking function */
	return (3.535f * a + 2.181f * a2)
		/ (1.0f + 2.276f * a + 2.577f * a2);
}



NORI_NAMESPACE_END