#if !defined(__MICROFACET_DISCRETE_H)
#define __MICROFACET_DISCRETE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/properties.h>
#include <boost/algorithm/string.hpp>
#include <queue>

MTS_NAMESPACE_BEGIN

class DiscreteMicrofacetDistribution
{
  public:
    /// Supported distribution types
    enum EType
    {
        /// Beckmann distribution derived from Gaussian random surfaces
        EBeckmann = 0
    };

    /**
	 * Create an isotropic microfacet distribution of the specified type
	 *
	 * \param type
	 *     The desired type of microfacet distribution
	 * \param alpha
	 *     The surface roughness
	 */
    inline DiscreteMicrofacetDistribution(EType type, Float alpha, uint32_t totalFacets, bool sampleVisible = true)
        : m_type(type), m_alphaU(alpha), m_alphaV(alpha), m_totalFacets(totalFacets), m_sampleVisible(sampleVisible),
          m_exponentU(0.0f), m_exponentV(0.0f)
    {
        m_alphaU = std::max(m_alphaU, (Float)1e-4f);
        m_alphaV = std::max(m_alphaV, (Float)1e-4f);
    }

    /**
	 * Create an anisotropic microfacet distribution of the specified type
	 *
	 * \param type
	 *     The desired type of microfacet distribution
	 * \param alphaU
	 *     The surface roughness in the tangent direction
	 * \param alphaV
	 *     The surface roughness in the bitangent direction
	 */
    inline DiscreteMicrofacetDistribution(EType type, Float alphaU, Float alphaV, uint32_t totalFacets, bool sampleVisible = true)
        : m_type(type), m_alphaU(alphaU), m_alphaV(alphaV), m_totalFacets(totalFacets), m_sampleVisible(sampleVisible),
          m_exponentU(0.0f), m_exponentV(0.0f)
    {
        m_alphaU = std::max(m_alphaU, (Float)1e-4f);
        m_alphaV = std::max(m_alphaV, (Float)1e-4f);
    }

    /**
	 * \brief Create a microfacet distribution from a Property data
	 * structure
	 */
    DiscreteMicrofacetDistribution(const Properties &props, EType type = EBeckmann,
                                   Float alphaU = 0.1f, Float alphaV = 0.1f, bool sampleVisible = true)
        : m_type(type), m_alphaU(alphaU), m_alphaV(alphaV), m_exponentU(0.0f),
          m_exponentV(0.0f)
    {

        if (props.hasProperty("distribution"))
        {
            std::string distr = boost::to_lower_copy(props.getString("distribution"));
            if (distr == "beckmann")
                m_type = EBeckmann;
            else
                SLog(EError, "Specified an invalid distribution \"%s\", currently only "
                             "support \"beckmann\"!",
                     distr.c_str());
        }

        if (props.hasProperty("alpha"))
        {
            m_alphaU = m_alphaV = props.getFloat("alpha");
            if (props.hasProperty("alphaU") || props.hasProperty("alphaV"))
                SLog(EError, "Microfacet model: please specify either 'alpha' or 'alphaU'/'alphaV'.");
        }
        else if (props.hasProperty("alphaU") || props.hasProperty("alphaV"))
        {
            if (!props.hasProperty("alphaU") || !props.hasProperty("alphaV"))
                SLog(EError, "Microfacet model: both 'alphaU' and 'alphaV' must be specified.");
            if (props.hasProperty("alpha"))
                SLog(EError, "Microfacet model: please specify either 'alpha' or 'alphaU'/'alphaV'.");
            m_alphaU = props.getFloat("alphaU");
            m_alphaV = props.getFloat("alphaV");
        }

        if (m_alphaU == 0 || m_alphaV == 0)
        {
            SLog(EWarn, "Cannot create a microfacet distribution with alphaU/alphaV=0 (clamped to 0.0001)."
                        "Please use the corresponding smooth reflectance model to get zero roughness.");
        }

        m_alphaU = std::max(m_alphaU, (Float)1e-4f);
        m_alphaV = std::max(m_alphaV, (Float)1e-4f);

        if (props.hasProperty("totalFacets"))
        {
            m_totalFacets = props.getInteger("totalFacets");
        }
        else
        {
            SLog(EError, "Discrete microfacet model: please specify 'totalFacets'.");
        }

        m_sampleVisible = props.getBoolean("sampleVisible", sampleVisible);
    }

    /// Return the distribution type
    inline EType getType() const { return m_type; }

    /// Return the roughness (isotropic case)
    inline Float getAlpha() const { return m_alphaU; }

    /// Return the roughness along the tangent direction
    inline Float getAlphaU() const { return m_alphaU; }

    /// Return the roughness along the bitangent direction
    inline Float getAlphaV() const { return m_alphaV; }

    /// Return the Phong exponent (isotropic case)
    inline Float getExponent() const { return m_exponentU; }

    /// Return the Phong exponent along the tangent direction
    inline Float getExponentU() const { return m_exponentU; }

    /// Return the Phong exponent along the bitangent direction
    inline Float getExponentV() const { return m_exponentV; }

    /// Return the total number of facets
    inline uint32_t getTotalFacets() const { return m_totalFacets; }

    /// Return whether or not only visible normals are sampled?
    inline bool getSampleVisible() const { return m_sampleVisible; }

    /// Is this an anisotropic microfacet distribution?
    inline bool isAnisotropic() const { return m_alphaU != m_alphaV; }

    /// Is this an anisotropic microfacet distribution?
    inline bool isIsotropic() const { return m_alphaU == m_alphaV; }

    /// Scale the roughness values by some constant
    inline void scaleAlpha(Float value)
    {
        m_alphaU *= value;
        m_alphaV *= value;
    }

    /**
	 * \brief Evaluate the discrete microfacet distribution function
	 */
    inline Float eval(const BSDFSamplingRecord &bRec) const
    {
        /* Calculate the reflection half-vector */
        Vector m = normalize(bRec.wo + bRec.wi);

        if (Frame::cosTheta(m) <= 0)
            return 0.0f;

        // pixel footprint in texture space
        const auto &its = bRec.its;
        Point2 center = its.uv;
        Point2 extentU(its.dudx, its.dvdx);
        Point2 extentV(its.dudy, its.dvdy);

        Float cosTheta2 = Frame::cosTheta2(m);
        Float beckmannExponent = ((m.x * m.x) / (m_alphaU * m_alphaU) + (m.y * m.y) / (m_alphaV * m_alphaV)) / cosTheta2;

        Float result;
        switch (m_type)
        {
        case EBeckmann:
        {
            /* Beckmann distribution function for Gaussian random surfaces - [Walter 2005] evaluation */
            result = math::fastexp(-beckmannExponent) /
                     (M_PI * m_alphaU * m_alphaV * cosTheta2 * cosTheta2);
        }
        break;

        default:
            SLog(EError, "Invalid distribution type!");
            return -1;
        }

        /* Prevent potential numerical issues in other stages of the model */
        if (result * Frame::cosTheta(m) < 1e-20f)
            result = 0;

        return result;
    }

    /**
	 * \brief Evaluate the smooth microfacet distribution function
	 *
	 * \param m
	 *     The microfacet normal
	 */
    inline Float eval(const Vector &m) const
    {
        if (Frame::cosTheta(m) <= 0)
            return 0.0f;

        Float cosTheta2 = Frame::cosTheta2(m);
        Float beckmannExponent = ((m.x * m.x) / (m_alphaU * m_alphaU) + (m.y * m.y) / (m_alphaV * m_alphaV)) / cosTheta2;

        Float result;
        switch (m_type)
        {
        case EBeckmann:
        {
            /* Beckmann distribution function for Gaussian random surfaces - [Walter 2005] evaluation */
            result = math::fastexp(-beckmannExponent) /
                     (M_PI * m_alphaU * m_alphaV * cosTheta2 * cosTheta2);
        }
        break;

        default:
            SLog(EError, "Invalid distribution type!");
            return -1;
        }

        /* Prevent potential numerical issues in other stages of the model */
        if (result * Frame::cosTheta(m) < 1e-20f)
            result = 0;

        return result;
    }

    /**
	 * \brief Wrapper function which calls \ref sampleAll() or \ref sampleVisible()
	 * depending on the parameters of this class
	 */
    inline Normal sample(const Vector &wi, const Point2 &sample, Float &pdf) const
    {
        Normal m;
        if (m_sampleVisible)
        {
            m = sampleVisible(wi, sample);
            pdf = pdfVisible(wi, m);
        }
        else
        {
            m = sampleAll(sample, pdf);
        }
        return m;
    }

    /**
	 * \brief Wrapper function which calls \ref sampleAll() or \ref sampleVisible()
	 * depending on the parameters of this class
	 */
    inline Normal sample(const Vector &wi, const Point2 &sample) const
    {
        Normal m;
        if (m_sampleVisible)
        {
            m = sampleVisible(wi, sample);
        }
        else
        {
            Float pdf;
            m = sampleAll(sample, pdf);
        }
        return m;
    }

    /**
	 * \brief Wrapper function which calls \ref pdfAll() or \ref pdfVisible()
	 * depending on the parameters of this class
	 */
    inline Float pdf(const Vector &wi, const Vector &m) const
    {
        if (m_sampleVisible)
            return pdfVisible(wi, m);
        else
            return pdfAll(m);
    }

    /**
	 * \brief Draw a sample from the microfacet normal distribution
	 * (including *all* normals) and return the associated
	 * probability density
	 *
	 * \param sample
	 *    A uniformly distributed 2D sample
	 * \param pdf
	 *    The probability density wrt. solid angles
	 */
    inline Normal sampleAll(const Point2 &sample, Float &pdf) const
    {
        /* The azimuthal component is always selected
		   uniformly regardless of the distribution */
        Float cosThetaM = 0.0f;
        Float sinPhiM, cosPhiM;
        Float alphaSqr;

        switch (m_type)
        {
        case EBeckmann:
        {
            /* Beckmann distribution function for Gaussian random surfaces */
            if (isIsotropic())
            {
                /* Sample phi component (isotropic case) */
                math::sincos((2.0f * M_PI) * sample.y, &sinPhiM, &cosPhiM);

                alphaSqr = m_alphaU * m_alphaU;
            }
            else
            {
                /* Sample phi component (anisotropic case) */
                Float phiM = std::atan(m_alphaV / m_alphaU *
                                       std::tan(M_PI + 2 * M_PI * sample.y)) +
                             M_PI * std::floor(2 * sample.y + 0.5f);
                math::sincos(phiM, &sinPhiM, &cosPhiM);

                Float cosSc = cosPhiM / m_alphaU, sinSc = sinPhiM / m_alphaV;
                alphaSqr = 1.0f / (cosSc * cosSc + sinSc * sinSc);
            }

            /* Sample theta component */
            Float tanThetaMSqr = alphaSqr * -math::fastlog(1.0f - sample.x);
            cosThetaM = 1.0f / std::sqrt(1.0f + tanThetaMSqr);

            /* Compute probability density of the sampled position */
            pdf = (1.0f - sample.x) / (M_PI * m_alphaU * m_alphaV * cosThetaM * cosThetaM * cosThetaM);
        }
        break;

        default:
            SLog(EError, "Invalid distribution type!");
            pdf = -1;
            return Vector(-1);
        }

        /* Prevent potential numerical issues in other stages of the model */
        if (pdf < 1e-20f)
            pdf = 0;

        Float sinThetaM = std::sqrt(
            std::max((Float)0, 1 - cosThetaM * cosThetaM));

        return Vector(
            sinThetaM * cosPhiM,
            sinThetaM * sinPhiM,
            cosThetaM);
    }

    /**
	 * \brief Returns the density function associated with
	 * the \ref sampleAll() function.
	 *
	 * \param m
	 *     The microfacet normal
	 */
    inline Float pdfAll(const Vector &m) const
    {
        /* PDF is just D(m) * cos(theta_M) */
        return eval(m) * Frame::cosTheta(m);
    }

    /**
	 * \brief Draw a sample from the distribution of visible normals
	 * and return the associated probability density
	 *
	 * \param _wi
	 *    A reference direction that defines the set of visible normals
	 * \param sample
	 *    A uniformly distributed 2D sample
	 * \param pdf
	 *    The probability density wrt. solid angles
	 */
    inline Normal sampleVisible(const Vector &_wi, const Point2 &sample) const
    {
        /* Step 1: stretch wi */
        Vector wi = normalize(Vector(
            m_alphaU * _wi.x,
            m_alphaV * _wi.y,
            _wi.z));

        /* Get polar coordinates */
        Float theta = 0, phi = 0;
        if (wi.z < (Float)0.99999)
        {
            theta = std::acos(wi.z);
            phi = std::atan2(wi.y, wi.x);
        }
        Float sinPhi, cosPhi;
        math::sincos(phi, &sinPhi, &cosPhi);

        /* Step 2: simulate P22_{wi}(slope.x, slope.y, 1, 1) */
        Vector2 slope = sampleVisible11(theta, sample);

        /* Step 3: rotate */
        slope = Vector2(
            cosPhi * slope.x - sinPhi * slope.y,
            sinPhi * slope.x + cosPhi * slope.y);

        /* Step 4: unstretch */
        slope.x *= m_alphaU;
        slope.y *= m_alphaV;

        /* Step 5: compute normal */
        Float normalization = (Float)1 / std::sqrt(slope.x * slope.x + slope.y * slope.y + (Float)1.0);

        return Normal(
            -slope.x * normalization,
            -slope.y * normalization,
            normalization);
    }

    /// Implements the probability density of the function \ref sampleVisible()
    Float pdfVisible(const Vector &wi, const Vector &m) const
    {
        if (Frame::cosTheta(wi) == 0)
            return 0.0f;
        return smithG1(wi, m) * absDot(wi, m) * eval(m) / std::abs(Frame::cosTheta(wi));
    }

    /**
	 * \brief Smith's shadowing-masking function G1 for each
	 * of the supported microfacet distributions
	 *
	 * \param v
	 *     An arbitrary direction
	 * \param m
	 *     The microfacet normal
	 */
    Float smithG1(const Vector &v, const Vector &m) const
    {
        /* Ensure consistent orientation (can't see the back
		   of the microfacet from the front and vice versa) */
        if (dot(v, m) * Frame::cosTheta(v) <= 0)
            return 0.0f;

        /* Perpendicular incidence -- no shadowing/masking */
        Float tanTheta = std::abs(Frame::tanTheta(v));
        if (tanTheta == 0.0f)
            return 1.0f;

        Float alpha = projectRoughness(v);
        switch (m_type)
        {
        case EBeckmann:
        {
            Float a = 1.0f / (alpha * tanTheta);
            if (a >= 1.6f)
                return 1.0f;

            /* Use a fast and accurate (<0.35% rel. error) rational
					   approximation to the shadowing-masking function */
            Float aSqr = a * a;
            return (3.535f * a + 2.181f * aSqr) / (1.0f + 2.276f * a + 2.577f * aSqr);
        }
        break;

        default:
            SLog(EError, "Invalid distribution type!");
            return -1.0f;
        }
    }

    /**
	 * \brief Separable shadow-masking function based on Smith's
	 * one-dimensional masking model
	 */
    inline Float G(const Vector &wi, const Vector &wo, const Vector &m) const
    {
        return smithG1(wi, m) * smithG1(wo, m);
    }

    /// Return a string representation of the name of a distribution
    inline static std::string distributionName(EType type)
    {
        switch (type)
        {
        case EBeckmann:
            return "beckmann";
            break;
        default:
            return "invalid";
            break;
        }
    }

    /// Return a string representation of the contents of this instance
    std::string toString() const
    {
        return formatString("DiscreteMicrofacetDistribution[type=\"%s\", alphaU=%f, alphaV=%f]",
                            distributionName(m_type).c_str(), m_alphaU, m_alphaV);
    }

  protected:
    /// Compute the effective roughness projected on direction \c v
    inline Float projectRoughness(const Vector &v) const
    {
        Float invSinTheta2 = 1 / Frame::sinTheta2(v);

        if (isIsotropic() || invSinTheta2 <= 0)
            return m_alphaU;

        Float cosPhi2 = v.x * v.x * invSinTheta2;
        Float sinPhi2 = v.y * v.y * invSinTheta2;

        return std::sqrt(cosPhi2 * m_alphaU * m_alphaU + sinPhi2 * m_alphaV * m_alphaV);
    }

    /// Compute the interpolated roughness for the Phong model
    inline Float interpolatePhongExponent(const Vector &v) const
    {
        const Float sinTheta2 = Frame::sinTheta2(v);

        if (isIsotropic() || sinTheta2 <= RCPOVERFLOW)
            return m_exponentU;

        Float invSinTheta2 = 1 / sinTheta2;
        Float cosPhi2 = v.x * v.x * invSinTheta2;
        Float sinPhi2 = v.y * v.y * invSinTheta2;

        return m_exponentU * cosPhi2 + m_exponentV * sinPhi2;
    }

    /**
	 * \brief Visible normal sampling code for the alpha=1 case
	 *
	 * Source: supplemental material of "Importance Sampling
	 * Microfacet-Based BSDFs using the Distribution of Visible Normals"
	 */
    Vector2 sampleVisible11(Float thetaI, Point2 sample) const
    {
        const Float SQRT_PI_INV = 1 / std::sqrt(M_PI);
        Vector2 slope;

        switch (m_type)
        {
        case EBeckmann:
        {
            /* Special case (normal incidence) */
            if (thetaI < 1e-4f)
            {
                Float sinPhi, cosPhi;
                Float r = std::sqrt(-math::fastlog(1.0f - sample.x));
                math::sincos(2 * M_PI * sample.y, &sinPhi, &cosPhi);
                return Vector2(r * cosPhi, r * sinPhi);
            }

            /* The original inversion routine from the paper contained
					   discontinuities, which causes issues for QMC integration
					   and techniques like Kelemen-style MLT. The following code
					   performs a numerical inversion with better behavior */
            Float tanThetaI = std::tan(thetaI);
            Float cotThetaI = 1 / tanThetaI;

            /* Search interval -- everything is parameterized
					   in the erf() domain */
            Float a = -1, c = math::erf(cotThetaI);
            Float sample_x = std::max(sample.x, (Float)1e-6f);

            /* Start with a good initial guess */
            //Float b = (1-sample_x) * a + sample_x * c;

            /* We can do better (inverse of an approximation computed in Mathematica) */
            Float fit = 1 + thetaI * (-0.876f + thetaI * (0.4265f - 0.0594f * thetaI));
            Float b = c - (1 + c) * std::pow(1 - sample_x, fit);

            /* Normalization factor for the CDF */
            Float normalization = 1 / (1 + c + SQRT_PI_INV * tanThetaI * std::exp(-cotThetaI * cotThetaI));

            int it = 0;
            while (++it < 10)
            {
                /* Bisection criterion -- the oddly-looking
						   boolean expression are intentional to check
						   for NaNs at little additional cost */
                if (!(b >= a && b <= c))
                    b = 0.5f * (a + c);

                /* Evaluate the CDF and its derivative
						   (i.e. the density function) */
                Float invErf = math::erfinv(b);
                Float value = normalization * (1 + b + SQRT_PI_INV * tanThetaI * std::exp(-invErf * invErf)) - sample_x;
                Float derivative = normalization * (1 - invErf * tanThetaI);

                if (std::abs(value) < 1e-5f)
                    break;

                /* Update bisection intervals */
                if (value > 0)
                    c = b;
                else
                    a = b;

                b -= value / derivative;
            }

            /* Now convert back into a slope value */
            slope.x = math::erfinv(b);

            /* Simulate Y component */
            slope.y = math::erfinv(2.0f * std::max(sample.y, (Float)1e-6f) - 1.0f);
        };
        break;

        default:
            SLog(EError, "Invalid distribution type!");
            return Vector2(-1);
        };
        return slope;
    }

    /// Helper routine: convert from Beckmann-style roughness values to Phong exponents (Walter et al.)
    void computePhongExponent()
    {
        m_exponentU = std::max(2.0f / (m_alphaU * m_alphaU) - 2.0f, (Float)0.0f);
        m_exponentV = std::max(2.0f / (m_alphaV * m_alphaV) - 2.0f, (Float)0.0f);
    }

    /// Helper routine: sample the azimuthal part of the first quadrant of the A&S distribution
    void sampleFirstQuadrant(Float u1, Float &phi, Float &exponent) const
    {
        Float cosPhi, sinPhi;
        phi = std::atan(
            std::sqrt((m_exponentU + 2.0f) / (m_exponentV + 2.0f)) *
            std::tan(M_PI * u1 * 0.5f));
        math::sincos(phi, &sinPhi, &cosPhi);
        /* Return the interpolated roughness */
        exponent = m_exponentU * cosPhi * cosPhi + m_exponentV * sinPhi * sinPhi;
    }

    // [Algorithm 1]
    // for now we do the spatial-directional split in lock-step
    uint32_t countParticles() const
    {
        // std::queue<Node> queue;

        return 0;
    }

  protected:
    EType m_type;
    Float m_alphaU, m_alphaV;
    bool m_sampleVisible;
    Float m_exponentU, m_exponentV;
    uint32_t m_totalFacets;
};

MTS_NAMESPACE_END

#endif /* __MICROFACET_DISCRETE_H */
