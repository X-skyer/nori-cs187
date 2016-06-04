#ifndef __SAMPLE_H__
#define __SAMPLE_H__

#pragma once
#include <nori/common.h>
#include <nori/vector.h>
#include <nori/color.h>
#include <memory>

NORI_NAMESPACE_BEGIN

// We are following the book for our distribution as this is highly unclear as in how to get pdfs from 2d distributions
// for our env lights.
struct Distribution1D {

	Distribution1D(const float *f, int n) : func(f, f + n), cdf(n + 1)
	{
		// Compute integral of step function at $x_i$
		cdf[0] = 0;
		for (int i = 1; i < n + 1; ++i) cdf[i] = cdf[i - 1] + func[i - 1] / n;

		// Transform step function integral into CDF
		funcInt = cdf[n];
		if (funcInt == 0) {
			for (int i = 1; i < n + 1; ++i) cdf[i] = float(i) / float(n);
		}
		else {
			for (int i = 1; i < n + 1; ++i) cdf[i] /= funcInt;
		}
	}

	int count() const { return static_cast<int>(func.size()); }

	float sample_continuous(float u, float *pdf, int *off = nullptr) const
	{
		// Find surrounding CDF segments and _offset_
		std::vector<float>::const_iterator entry =
			std::lower_bound(cdf.begin(), cdf.end(), u);

		size_t index = (size_t)std::max((ptrdiff_t)0, entry - cdf.begin() - 1);
		int offset = std::min(index, cdf.size() - 2);

		if (off) *off = offset;

		// Compute offset along CDF segment
		float du = u - cdf[offset];

		if ((cdf[offset + 1] - cdf[offset]) > 0)
			du /= (cdf[offset + 1] - cdf[offset]);
		assert(!std::isnan(du));

		// Compute PDF for sampled offset
		if (pdf) *pdf = func[offset] / funcInt;

		return (offset + du) / count();
	}

	int sample_discrete(float u, float *pdf = nullptr, float *uRemapped = nullptr) const
	{
		// Find surrounding CDF segments and _offset_
		std::vector<float>::const_iterator entry =
			std::lower_bound(cdf.begin(), cdf.end(), u);

		size_t index = (size_t)std::max((ptrdiff_t)0, entry - cdf.begin() - 1);
		int offset = std::min(index, cdf.size() - 2);

		if (pdf) *pdf = func[offset] / (funcInt * count());
		if (uRemapped)
			*uRemapped = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
		if (uRemapped) assert(*uRemapped >= 0.f && *uRemapped <= 1.f);
		return offset;
	}

	// Return the discrete pdf
	float pdf(int index) const {
		assert(index >= 0 && index < count());
		return func[index] / (funcInt * count());
	}

	// Distribution1D Public Data
	std::vector<float> func, cdf;
	float funcInt;
};

// Create a discrete 2d pdf for sampling images
struct Distribution2D
{
public:

	Distribution2D(const float *data, int nu, int nv)
	{
		// copy the data
		pConditionalV.reserve(nv);
		for (int v = 0; v < nv; ++v) {
			// Compute conditional sampling distribution for the column with the row
			pConditionalV.emplace_back(new Distribution1D(&data[v * nu], nu));
		}
		// Compute marginal sampling distribution corresponding row
		std::vector<float> marginalFunc;
		marginalFunc.reserve(nv);
		for (int v = 0; v < nv; ++v)
			marginalFunc.push_back(pConditionalV[v]->funcInt);
		pMarginal.reset(new Distribution1D(&marginalFunc[0], nv));
	}

	Point2f sample_continuous(const Point2f &u, float *pdf) const
	{
		float pdfs[2];
		int v;
		float d1 = pMarginal->sample_continuous(u[1], &pdfs[1], &v);
		float d0 = pConditionalV[v]->sample_continuous(u[0], &pdfs[0]);
		*pdf = pdfs[0] * pdfs[1];
		return Point2f(d0, d1);
	}

	float pdf(const Point2f &p) const
	{
		int iu = clamp(int(p[0] * pConditionalV[0]->count()), 0, pConditionalV[0]->count() - 1);
		int iv = clamp(int(p[1] * pMarginal->count()), 0, pMarginal->count() - 1);
		return pConditionalV[iv]->func[iu] / pMarginal->funcInt;
	}

private:
	// Distribution2D Private Data
	std::vector<std::unique_ptr<Distribution1D>> pConditionalV;
	std::unique_ptr<Distribution1D> pMarginal;
};

NORI_NAMESPACE_END

#endif
