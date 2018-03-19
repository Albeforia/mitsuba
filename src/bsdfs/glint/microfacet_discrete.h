#if !defined(__MICROFACET_DISCRETE_H)
#define __MICROFACET_DISCRETE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/properties.h>
#include <boost/algorithm/string.hpp>
#include <queue>
#include "spherical_conic_section.h"
#include "multinomial.h"

#define inAABB(p, min, max) ((min).x <= (p).x && (p).x <= (max).x && (min).y <= (p).y && (p).y <= (max).y)
#define inHalfPlane(p, p1, p2) (((p).x - (p2).x) * ((p1).y - (p2).y) - ((p1).x - (p2).x) * ((p).y - (p2).y) <= 0)
#define inTriangle(p, p1, p2, p3) (inHalfPlane(p, p1, p2) && inHalfPlane(p, p2, p3) && inHalfPlane(p, p3, p1))
#define aabbArea(aabb) (((aabb).max - (aabb).min).x * ((aabb).max - (aabb).min).y)

// NOTE average query area(i.e. pixel footprint) varies according to resolution and sample count
// the value below is pre-computed at the resolution of 683x512 with sample count 1
#define AVG_QUERY_AREA 0.0000016662

#define AVG_QUERY_SOLID_ANGLE 0.0131010329

MTS_NAMESPACE_BEGIN

using Parallelogram = std::array<Vector2, 4>;

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
    inline Float eval(const Vector &m,
                      const Parallelogram &pixel, const SphericalConicSection &scs,
                      const std::unordered_map<std::string, Float> &integrations,
                      size_t sampleCount) const
    {
        if (Frame::cosTheta(m) <= 0)
            return 0.0f;

        Float result = countParticles(pixel, scs, integrations, sampleCount) / static_cast<float>(m_totalFacets);

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

    // [Algorithm 1]
    uint32_t countParticles(const Parallelogram &pixel, const SphericalConicSection &scs,
                            const std::unordered_map<std::string, Float> &integrations,
                            size_t sampleCount) const
    {
        static unsigned int seed = 0;

        uint32_t count = 0;
        std::queue<Node> queue;
        Float p0 = integrations.at("00"),
              p1 = integrations.at("01"),
              p2 = integrations.at("02"),
              p3 = integrations.at("03");
        auto sum = p0 + p1 + p2 + p3;
        std::array<float, 4> pv{p0 / sum, p1 / sum, p2 / sum, p3 / sum};
        auto root_counts = multinomial(m_totalFacets, pv, ++seed);
        queue.emplace(Point2(0, 0), Point2(1, 1), Vector(0, 0, 1), Vector(1, 0, 0), Vector(0, 1, 0), "00", root_counts[0]);
        queue.emplace(Point2(0, 0), Point2(1, 1), Vector(0, 0, 1), Vector(-1, 0, 0), Vector(0, 1, 0), "01", root_counts[1]);
        queue.emplace(Point2(0, 0), Point2(1, 1), Vector(0, 0, 1), Vector(1, 0, 0), Vector(0, -1, 0), "02", root_counts[2]);
        queue.emplace(Point2(0, 0), Point2(1, 1), Vector(0, 0, 1), Vector(-1, 0, 0), Vector(0, -1, 0), "03", root_counts[3]);
        while (!queue.empty())
        {
            const auto &curr = queue.front();
            auto overlapSpatial = curr.overlap(pixel);
            auto overlapDirectional = curr.overlap(scs);
            if (overlapSpatial == 0 || overlapDirectional == 0 || curr.m_count == 0)
            {
            }
            else if (overlapSpatial == 2 && overlapDirectional == 2)
            {
                count += curr.m_count;
            }
            else
            {
                if ((AVG_QUERY_AREA / sampleCount) / aabbArea(curr.m_spatial) < AVG_QUERY_SOLID_ANGLE / curr.m_directional.excess())
                {
                    std::array<float, 4> pv{0.25f, 0.25f, 0.25f, 0.25f};
                    auto counts = multinomial(curr.m_count, pv, ++seed);
                    auto min = curr.m_spatial.min;
                    auto max = curr.m_spatial.max;
                    auto center = (min + max) * 0.5f;
                    queue.emplace(min, center, curr.m_directional[0], curr.m_directional[1], curr.m_directional[2], curr.m_triID, counts[0]);
                    queue.emplace(Point2(center.x, min.y), Point2(max.x, center.y), curr.m_directional[0], curr.m_directional[1], curr.m_directional[2], curr.m_triID, counts[1]);
                    queue.emplace(center, max, curr.m_directional[0], curr.m_directional[1], curr.m_directional[2], curr.m_triID, counts[2]);
                    queue.emplace(Point2(min.x, center.y), Point2(center.x, max.y), curr.m_directional[0], curr.m_directional[1], curr.m_directional[2], curr.m_triID, counts[3]);
                }
                else
                {
                    Float p0 = integrations.count(curr.m_triID + "0") == 0 ? 0 : integrations.at(curr.m_triID + "0");
                    Float p1 = integrations.count(curr.m_triID + "1") == 0 ? 0 : integrations.at(curr.m_triID + "1");
                    Float p2 = integrations.count(curr.m_triID + "2") == 0 ? 0 : integrations.at(curr.m_triID + "2");
                    Float p3 = integrations.count(curr.m_triID + "3") == 0 ? 0 : integrations.at(curr.m_triID + "3");
                    auto sum = p0 + p1 + p2 + p3;
                    std::array<float, 4> pv{p0 / sum, p1 / sum, p2 / sum, p3 / sum};
                    auto counts = multinomial(curr.m_count, pv, ++seed);
                    auto children = curr.m_directional.split();
                    queue.emplace(curr.m_spatial.min, curr.m_spatial.max, children[0][0], children[0][1], children[0][2], curr.m_triID + "0", counts[0]);
                    queue.emplace(curr.m_spatial.min, curr.m_spatial.max, children[1][0], children[1][1], children[1][2], curr.m_triID + "1", counts[1]);
                    queue.emplace(curr.m_spatial.min, curr.m_spatial.max, children[2][0], children[2][1], children[2][2], curr.m_triID + "2", counts[2]);
                    queue.emplace(curr.m_spatial.min, curr.m_spatial.max, children[3][0], children[3][1], children[3][2], curr.m_triID + "3", counts[3]);
                }
            }
            queue.pop();
        }
        return count;
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

  protected:
    EType m_type;
    Float m_alphaU, m_alphaV;
    bool m_sampleVisible;
    Float m_exponentU, m_exponentV;
    uint32_t m_totalFacets;

    struct Node
    {
        AABB2 m_spatial;
        SphericalTriangle m_directional;
        std::string m_triID;
        uint32_t m_count;

        Node(){};

        Node(AABB2::PointType min, AABB2::PointType max, Vector a, Vector b, Vector c,
             const std::string &triID, uint32_t count = 0)
            : m_spatial(min, max), m_directional(a, b, c), m_triID(triID), m_count(count) {}

        // this test include the case that aabb contains the parallelogram
        bool intersect(const Parallelogram &paral) const
        {
            auto center = (m_spatial.min + m_spatial.max) * 0.5f;
            auto extent = (m_spatial.max - m_spatial.min) * 0.5f;

            // intersect four edges of parallelogram against the aabb
            int result = 0;
            for (int i = 0; i < 4; i++)
            {
                auto p = paral[i];
                auto delta = paral[(i + 1) % 4] - p;
                auto scaleX = 1.0f / delta.x, scaleY = 1.0f / delta.y;
                int signX = scaleX < 0 ? -1 : 1, signY = scaleY < 0 ? -1 : 1;
                auto nearTimeX = (center.x - signX * (extent.x) - p.x) * scaleX,
                     nearTimeY = (center.y - signY * (extent.y) - p.y) * scaleY,
                     farTimeX = (center.x + signX * (extent.x) - p.x) * scaleX,
                     farTimeY = (center.y + signY * (extent.y) - p.y) * scaleY;
                if (nearTimeX > farTimeY || nearTimeY > farTimeX)
                {
                    continue;
                }
                auto nearTime = std::max(nearTimeX, nearTimeY),
                     farTime = std::min(farTimeX, farTimeY);
                if (nearTime >= 1 || farTime <= 0)
                {
                    continue;
                }
                result++;
                break;
            }
            // SLog(EInfo, (std::string("intersect: ") + std::to_string(result)).c_str());
            return result != 0;
        }

        int overlap(const Parallelogram &paral) const
        {
            // test if the parallelogram contains the aabb
            auto center = (m_spatial.min + m_spatial.max) * 0.5f;
            auto extent = (m_spatial.max - m_spatial.min) * 0.5f;
            AABB2::PointType aabb[]{
                {center.x + extent.x, center.y + extent.y},
                {center.x + extent.x, center.y - extent.y},
                {center.x - extent.x, center.y + extent.y},
                {center.x - extent.x, center.y - extent.y}};
            // vertices should be in clockwise order
            if ((inTriangle(aabb[0], paral[0], paral[2], paral[1]) || inTriangle(aabb[0], paral[0], paral[3], paral[2])) &&
                (inTriangle(aabb[1], paral[0], paral[2], paral[1]) || inTriangle(aabb[1], paral[0], paral[3], paral[2])) &&
                (inTriangle(aabb[2], paral[0], paral[2], paral[1]) || inTriangle(aabb[2], paral[0], paral[3], paral[2])) &&
                (inTriangle(aabb[3], paral[0], paral[2], paral[1]) || inTriangle(aabb[3], paral[0], paral[3], paral[2])))
            {
                // SLog(EInfo, "overlap: aabb is inside the parallelogram");
                return 2;
            }

            if (intersect(paral))
            {
                return 1;
            }

            return 0;
        }

        bool intersect(const SphericalConicSection &scs) const
        {
            return scs.intersect(m_directional);
        }

        int overlap(const SphericalConicSection &scs) const
        {
            return scs.overlap(m_directional);
        }
    };
};

MTS_NAMESPACE_END

#endif /* __MICROFACET_DISCRETE_H */
