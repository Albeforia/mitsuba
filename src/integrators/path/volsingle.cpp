#include <mitsuba/render/scene.h>

MTS_NAMESPACE_BEGIN

class SingleScatteringPathTracer : public MonteCarloIntegrator
{
public:
	SingleScatteringPathTracer(const Properties &props) : MonteCarloIntegrator(props)
	{
		m_stepSize = props.getFloat("stepSize", 0.1f);
	}

	/// Unserialize from a binary data stream
	SingleScatteringPathTracer(Stream *stream, InstanceManager *manager)
		: MonteCarloIntegrator(stream, manager) {}

	Spectrum Li(const RayDifferential &r, RadianceQueryRecord &rRec) const
	{
		const Scene *scene = rRec.scene;
		RayDifferential ray(r);

		/* Perform the first ray intersection (or ignore if the
		   intersection has already been provided). */
		rRec.rayIntersect(ray);
		Intersection &its = rRec.its;

		// if (m_maxDepth == 1)
		// 	rRec.type &= RadianceQueryRecord::EEmittedRadiance;

		Spectrum Li(0.0f);
		Spectrum Tr(1.0f);
		/*
			(TODO)1. camera is in medium, surface is in medium
			(TODO)2. camera is in medium, surface is not in medium
			3. camera is not in medium, surface is in medium
		    4. camera is not in medium, surface is not in medium
		*/
		while (rRec.depth <= m_maxDepth || m_maxDepth < 0)
		{
			if (!its.isValid())
			{
				/* If no intersection could be found, possibly return
					   attenuated radiance from a background luminaire */
				if ((rRec.type & RadianceQueryRecord::EEmittedRadiance) && !m_hideEmitters)
				{
					Spectrum value = Tr * scene->evalEnvironment(ray);
					if (rRec.medium)
						value *= rRec.medium->evalTransmittance(ray);
					Li += value;
				}
				break;
			}

			/* Possibly include emitted radiance if requested */
			if (its.isEmitter() && (rRec.type & RadianceQueryRecord::EEmittedRadiance) && !m_hideEmitters)
				Li += Tr * its.Le(-ray.d);

			/* If 'strictNormals'=true, the geometric and shading normals
				   must classify the incident direction to the same side */
			if (m_strictNormals && dot(ray.d, its.geoFrame.n) * Frame::cosTheta(its.wi) >= 0)
				break;

			const BSDF *bsdf = its.getBSDF(ray);
			if (!(its.getBSDF()->getType() & BSDF::ENull))
			{
				Li += Tr * surfaceDirect(ray, bsdf, rRec);
				break;
			}
			else
			{
				/* Hit a medium boundary ? */
				if (its.isMediumTransition())
				{
					Assert(rRec.medium = its.getTargetMedium(ray.d));

					Ray ray2(its.p, ray.d, ray.time);
					ray2.mint = 4.0f * Epsilon;
					Li += rayMarchMedium(ray2, rRec, Tr);
				}
				else
				{
					SLog(EWarn, "Hit something that is neither surface nor medium boundary?");
					break;
				}
			}

			/* Note the use of rRec.depth is inconsistent with other integrators */
			rRec.depth++;
		}

		return Li;
	}

	Spectrum surfaceDirect(const RayDifferential &ray, const BSDF *bsdf,
						   RadianceQueryRecord &rRec) const
	{
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		Spectrum Li(0.0f);

		/* ==================================================================== */
		/*                          Emitter sampling                          */
		/* ==================================================================== */

		DirectSamplingRecord dRec(its);
		/* Only use direct illumination sampling when the surface's
		   BSDF has smooth (i.e. non-Dirac delta) component */
		if (bsdf->getType() & BSDF::ESmooth)
		{
			DirectSamplingRecord dRec(its);
			int maxInteractions = m_maxDepth - rRec.depth - 1;

			Spectrum value = scene->sampleAttenuatedEmitterDirect(
				dRec, its, rRec.medium, maxInteractions,
				rRec.nextSample2D(), rRec.sampler);

			if (!value.isZero())
			{
				const Emitter *emitter = static_cast<const Emitter *>(dRec.object);

				/* Allocate a record for querying the BSDF */
				BSDFSamplingRecord bRec(its, its.toLocal(dRec.d));

				/* Evaluate BSDF * cos(theta) */
				const Spectrum bsdfVal = bsdf->eval(bRec);

				if (!bsdfVal.isZero() && (!m_strictNormals || dot(its.geoFrame.n, dRec.d) * Frame::cosTheta(bRec.wo) > 0))
				{
					/* Calculate prob. of sampling that direction using BSDF sampling */
					Float bsdfPdf = emitter->isOnSurface() ? bsdf->pdf(bRec) : 0;

					/* Weight using the power heuristic */
					const Float weight = miWeight(dRec.pdf * 0.5f, bsdfPdf * 0.5f);

					Li += value * bsdfVal * weight;
				}
			}
		}

		/* ==================================================================== */
		/*                            BSDF sampling                             */
		/* ==================================================================== */

		/* Sample BSDF * cos(theta) and also request the local density */
		Intersection bsdfIts;
		Float bsdfPdf;

		BSDFSamplingRecord bRec(its, rRec.sampler, ERadiance);
		Spectrum bsdfVal = bsdf->sample(bRec, bsdfPdf, rRec.nextSample2D());
		if (bsdfVal.isZero())
			return Li;

		/* Prevent light leaks due to the use of shading normals */
		const Vector wo = its.toWorld(bRec.wo);
		Float woDotGeoN = dot(its.geoFrame.n, wo);
		if (m_strictNormals && woDotGeoN * Frame::cosTheta(bRec.wo) <= 0)
			return Li;

		/* Trace a ray in this direction */
		Ray bsdfRay(its.p, wo, ray.time);

		Spectrum value;
		if (scene->rayIntersect(bsdfRay, bsdfIts))
		{
			/* Intersected something - check if it was an emitter */
			if (!bsdfIts.isEmitter())
				return Li;

			value = bsdfIts.Le(-bsdfRay.d);
			dRec.setQuery(bsdfRay, bsdfIts);
		}
		else
		{
			/* Intersected nothing -- perhaps there is an environment map? */
			const Emitter *env = scene->getEnvironmentEmitter();

			if (!env || (m_hideEmitters && bRec.sampledType == BSDF::ENull))
				return Li;

			value = env->evalEnvironment(RayDifferential(bsdfRay));
			if (!env->fillDirectSamplingRecord(dRec, bsdfRay))
				return Li;
		}

		/* Compute the prob. of generating that direction using the
			   implemented direct illumination sampling technique */
		const Float lumPdf = (!(bRec.sampledType & BSDF::EDelta)) ? scene->pdfEmitterDirect(dRec) : 0;

		/* Weight using the power heuristic */
		const Float weight = miWeight(bsdfPdf * 0.5f, lumPdf * 0.5f);

		Li += value * bsdfVal * weight;

		return Li;
	}

	Spectrum rayMarchMedium(const Ray &ray, RadianceQueryRecord &rRec, Spectrum &Tr) const
	{
		const Scene *scene = rRec.scene;
		Intersection &its = rRec.its;
		Spectrum inscattering(0.0f);

		/* Trace ray inside the medium */
		if (!scene->rayIntersect(ray, its))
		{
			/* Numerical issue ? */
			Tr = Spectrum(1.0f);
			return inscattering;
		}

		MediumSamplingRecord mRec;
		const PhaseFunction *phase = rRec.medium->getPhaseFunction();

		Float mint = ray.mint;
		Float maxt = its.t;

		/* Calculate volume integration via ray marching */
		int nSamples = (int)std::ceil((maxt - mint) / m_stepSize);
		Float step = (maxt - mint) / nSamples;

		/*  The first sample is placed randomly in the first segment */
		Float tPrev = mint;
		Float tCurr = mint + rRec.nextSample1D() * step;

		for (int i = 0; i < nSamples; ++i, tCurr += step, ++rRec.depth)
		{
			Tr *= rRec.medium->evalTransmittance(Ray(ray, tPrev, tCurr), rRec.sampler);

			/* Possibly terminate ray marching if transmittance is small */
			if (Tr.max() < 1e-3)
			{
				const Float continueProb = 0.5f;
				if (rRec.nextSample1D() > continueProb)
				{
					Tr = Spectrum(0.0f);
					break;
				}
				Tr /= continueProb;
			}

			/* Calculate in-scattering */
			mRec.transmittance = Tr;
			mRec.medium = rRec.medium;
			mRec.time = ray.time;
			mRec.t = tCurr;
			mRec.p = ray(mRec.t);
			if (mRec.medium->isHomogeneous())
			{
				mRec.sigmaA = mRec.medium->getSigmaA();
				mRec.sigmaS = mRec.medium->getSigmaS();
			}
			else
			{
				/* code */
			}

			DirectSamplingRecord dRec(mRec.p, mRec.time);
			int maxInteractions = m_maxDepth - rRec.depth - 1;

			Spectrum value = scene->sampleAttenuatedEmitterDirect(
				dRec, rRec.medium, maxInteractions,
				rRec.nextSample2D(), rRec.sampler);

			if (!value.isZero())
				inscattering += Tr * mRec.sigmaS * value * phase->eval(PhaseFunctionSamplingRecord(mRec, -ray.d, dRec.d));

			tPrev = tCurr;
		}

		const Medium *oldMedium = rRec.medium;
		rRec.medium = its.getTargetMedium(ray.d);
		/* For now, connected media (i.e. shared boundary) are not allowed.
		   Therefore we should be either heading out of the medium or hitting a surface */
		Assert(rRec.medium == nullptr);
		if (its.isMediumTransition() && (its.getBSDF()->getType() & BSDF::ENull))
		{
			Ray ray2(its.p, ray.d, ray.time);
			ray2.mint = 4.0f * Epsilon;
			scene->rayIntersect(ray2, its);
		}
		else
		{
			/* This is a point on surface but still in the medium */
			rRec.medium = oldMedium;
		}

		return inscattering * step;
	}

	inline Float miWeight(Float pdfA, Float pdfB) const
	{
		pdfA *= pdfA;
		pdfB *= pdfB;
		return pdfA / (pdfA + pdfB);
	}

	void serialize(Stream *stream, InstanceManager *manager) const
	{
		MonteCarloIntegrator::serialize(stream, manager);
	}

	std::string toString() const
	{
		std::ostringstream oss;
		oss << "SingleScatteringPathTracer[" << endl
			<< "  maxDepth = " << m_maxDepth << "," << endl
			<< "  rrDepth = " << m_rrDepth << "," << endl
			<< "  strictNormals = " << m_strictNormals << endl
			<< "  stepSize = " << m_stepSize << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()

private:
	Float m_stepSize;
};

MTS_IMPLEMENT_CLASS_S(SingleScatteringPathTracer, false, MonteCarloIntegrator)
MTS_EXPORT_PLUGIN(SingleScatteringPathTracer, "Single scattering path tracer");

MTS_NAMESPACE_END