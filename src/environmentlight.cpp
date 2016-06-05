#include <nori/emitter.h>
#include <nori/bitmap.h>
#include <nori/frame.h>
#include <nori/sample.h>
#include <cmath>

NORI_NAMESPACE_BEGIN

/**
	This class implements an environment light that is a far sphere from the center of the scene
*/

class EnvironmentLight : public Emitter
{
public:

	EnvironmentLight(const PropertyList& prop)
	{
		m_tex_filename = prop.getString("filename");
		m_localToWorld = prop.getTransform("toWorld", Transform());
		m_worldToLocal = m_localToWorld.getInverseMatrix();
		m_texture = Texture(m_tex_filename);		
		m_type = EmitterType::EMITTER_ENVIRONMENT;
	}

	~EnvironmentLight()
	{

	}

	// Setup the discrete pdfs so that it can be sampled
	virtual void activate()
	{
		Color3f *data = m_texture.data();
		// Compute luminance of each pixel
		float* m_luminance = new float[m_texture.get_width() * m_texture.get_height()];

		for (int y = 0; y < m_texture.get_height(); y++)
		{
			float sineTheta = sinf(M_PI * (y + 0.5f) / m_texture.get_height());
			for (int x = 0; x < m_texture.get_width(); x++)
			{
				int index = x + y * m_texture.get_width();
				float lum = data[index].getLuminance();
				m_luminance[index] = lum * sineTheta;
			}
		}

		m_pdf = std::unique_ptr<Distribution2D>(new Distribution2D(m_luminance, m_texture.get_width(), m_texture.get_height()));		
		delete[] m_luminance;
	}

	// Methods that have to be implemented for all types of emitters.
	Color3f sample(EmitterQueryRecord &lRec, const Point2f &sample, float optional_u) const
	{
		float pdf;
		Point2f sampled_pt = m_pdf->sample_continuous(sample, &pdf);

		float theta = sampled_pt.y() * M_PI;
		float phi = sampled_pt.x() * 2.0f * M_PI;
		float sin_theta = sin(theta);
		float cos_theta = cos(theta);
		float sin_phi = sin(phi);
		float cos_phi = cos(phi);
		
		lRec.emitter = this;
		lRec.dist = INFINITY;
		lRec.p = Vector3f(sampled_pt.x(), sampled_pt.y(), 0.0f);
		lRec.wi = m_localToWorld * Vector3f(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);			// probably need to transform
		lRec.pdf = pdf / (2.0f * M_PI * M_PI * sin_theta);

		if (lRec.pdf == 0.0f)
			return 0.0f;
		else
			return m_texture.getval(sampled_pt.x(), sampled_pt.y()) / lRec.pdf;
	}

	float pdf(const EmitterQueryRecord &lRec) const
	{
		Vector3f world_dir = lRec.wi;
		Vector3f local_dir = m_worldToLocal * world_dir;

		float theta = Frame::spherical_theta(local_dir);
		float phi = Frame::spherical_phi(local_dir);
		float sin_theta = sinf(theta);
		if (sin_theta == 0.0f) 
			return 0.0f;
		else
			return m_pdf->pdf(Point2f(phi * INV_TWOPI, theta * INV_PI)) / (2.0f * M_PI * M_PI * sin_theta);
	}

	Color3f eval(const EmitterQueryRecord &lRec) const
	{

		Vector3f world_dir = lRec.wi;
		Vector3f local_dir = m_worldToLocal * world_dir;

		// compute theta and phi
		float theta = acos(clamp(Frame::cosTheta(local_dir), -1.0f, 1.0f));
		float phi = atan2f(local_dir.x(), local_dir.y());

		// compute the corresponding pixel coords
		// we assume height ranges from 0 -> pi/2
		// we assume width ranges from 0 - 2pi
		float s = (M_PI + phi) * INV_TWOPI;
		float t = theta * INV_PI;

		return m_texture.getval(s, t);
	}

	std::string toString() const
	{
		return tfm::format("Environment Light : \n File : %s \n Dimensions : [%d,%d]", m_tex_filename, m_texture.get_width(), m_texture.get_height());
	}

private:
	Texture m_texture;
	Transform m_localToWorld;
	Transform m_worldToLocal;
	std::string m_tex_filename;
	std::unique_ptr<Distribution2D> m_pdf;
};

NORI_REGISTER_CLASS(EnvironmentLight, "environment");
NORI_NAMESPACE_END