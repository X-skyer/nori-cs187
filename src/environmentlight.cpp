#include <nori/emitter.h>
#include <nori/bitmap.h>
#include <nori/frame.h>
#include <nori/sample.h>

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
			for (int x = 0; x < m_texture.get_width(); x++)
			{
				int index = x + y * m_texture.get_width();
				float lum = data[index].getLuminance();
				float sineTheta = sinf(M_PI * (y + 0.5f) / m_texture.get_height());
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

		float theta = sampled_pt.x() * M_PI;
		float phi = sampled_pt.y() * 2.0f * M_PI;
		float sin_theta = sin(theta);
		float cos_theta = cos(theta);
		float sin_phi = sin(phi);
		float cos_phi = cos(phi);
		
		lRec.emitter = this;
		lRec.dist = INFINITY;
		lRec.p = Vector3f(sampled_pt.x(), sampled_pt.y(), 0.0f);
		lRec.wi = Vector3f(sin_theta * cos_phi, sin_theta * sin_phi, cos_theta);			// probably need to transform

		return m_texture.getval(sampled_pt.x(), sampled_pt.y());
	}

	float pdf(const EmitterQueryRecord &lRec) const
	{
		float theta = Frame::spherical_theta(lRec.wi);
		float phi = Frame::spherical_phi(lRec.wi);
		float sin_theta = sinf(theta);
		if (sin_theta == 0.0f) return 0.0f;
		return m_pdf->pdf(Point2f(phi * INV_TWOPI, theta * INV_PI)) / (2.0f * M_PI * M_PI * sin_theta);
	}

	Color3f eval(const EmitterQueryRecord &lRec) const
	{
		// compute theta and phi
		float theta = std::acos(Frame::cosTheta(lRec.wi));
		float phi = std::atan2f(lRec.wi.x(), lRec.wi.y());

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
	Transform m_transform;
	Texture m_texture;
	std::string m_tex_filename;
	std::unique_ptr<Distribution2D> m_pdf;
};

NORI_REGISTER_CLASS(EnvironmentLight, "environment");
NORI_NAMESPACE_END