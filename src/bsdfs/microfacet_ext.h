/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#if !defined(__MICROFACET_EXT_H)
#define __MICROFACET_EXT_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/frame.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/vmf.h>
#include <boost/algorithm/string.hpp>

MTS_NAMESPACE_BEGIN

/**
 * \brief Implementation of the Beckman and GGX / Trowbridge-Reitz microfacet
 * distributions and various useful sampling routines
 *
 * Based on the papers
 *
 *   "Microfacet Models for Refraction through Rough Surfaces"
 *    by Bruce Walter, Stephen R. Marschner, Hongsong Li, and Kenneth E. Torrance
 *
 * and
 *
 *   "Importance Sampling Microfacet-Based BSDFs using the Distribution of Visible Normals"
 *    by Eric Heitz and Eugene D'Eon
 *
 *  The visible normal sampling code was provided by Eric Heitz and Eugene D'Eon.
 */
class LayeredMicrofacetDistribution {
public:

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
	inline LayeredMicrofacetDistribution(Float kappa1, Float kappa2, Float eta)
		: m_kappa1(kappa1), m_kappa2(kappa2), m_eta(eta) {}

	/// Return the roughness along the tangent direction
	inline Float getKappa1() const { return m_kappa1; }

	/// Return the roughness along the bitangent direction
	inline Float getKappa2() const { return m_kappa2; }

    /// Return the refractive index of the layer
	inline Float getEta() const { return m_eta; }

	/// Is this an anisotropic microfacet distribution?
	inline bool isAnisotropic() const { return false; }

	/// Is this an anisotropic microfacet distribution?
	inline bool isIsotropic() const { return true; }

	inline Float eval0(const Vector &m) const {
		if (Frame::cosTheta(m) <= 0)
			return 0.0f;

        VonMisesFisherDistr vMF(m_kappa1);
		return vMF.eval(Frame::cosTheta(m));
	}

    inline Float eval1(const Vector &m, const Vector &wi) const {
		if (Frame::cosTheta(m) <= 0)
			return 0.0f;

        Vector t1 = refract(wi, Normal(0, 0, 1), m_eta);
        Vector t2 = reflect(-t1, Normal(0, 0, 1));

        Float cosThetaI = Frame::cosTheta(wi);
        Float cosThetaT1 = Frame::cosTheta(t1);
        Float cosThetaT2 = Frame::cosTheta(t2);
        Float temp1 = m_eta * cosThetaT2 + cosThetaI;
        Float jacobian1 = std::abs(cosThetaI) / (temp1 * temp1);
        Float jacobian2 = std::abs(cosThetaI)
                        / ( m_eta*m_eta * 4*std::abs(cosThetaT2)*std::abs(cosThetaT1) );
        Float jacobian3 = jacobian1;

        VonMisesFisherDistr vMF;
        Float kappa_temp = vMF.convolve(jacobian1*m_kappa1, jacobian2*m_kappa2);
        Float kappa = vMF.convolve(kappa_temp, jacobian3 * m_kappa1) *
                      (4 * std::abs(cosThetaI));
        vMF.setKappa(kappa);

		return vMF.eval(Frame::cosTheta(m));
	}

	inline Normal sample0(const Vector &wi, const Point2 &sample, Float &pdf) const {
        VonMisesFisherDistr vMF(m_kappa1);
        Normal m = Normal(vMF.sample(sample));

        pdf = vMF.eval(Frame::cosTheta(m));
        /* Prevent potential numerical issues in other stages of the model */
		if (pdf < 1e-20f)
			pdf = 0;

        return m;
	}

	inline Normal sample1(const Vector &wi, const Point2 &sample, Float &pdf) const {
        Vector t1 = refract(wi, Normal(0, 0, 1), m_eta);
        Vector t2 = reflect(-t1, Normal(0, 0, 1));

        Float cosThetaI = Frame::cosTheta(wi);
        Float cosThetaT1 = Frame::cosTheta(t1);
        Float cosThetaT2 = Frame::cosTheta(t2);
        Float temp1 = m_eta * cosThetaT2 + cosThetaI;
        Float jacobian1 = std::abs(cosThetaI) / (temp1 * temp1);
        Float jacobian2 = std::abs(cosThetaI)
                        / ( m_eta*m_eta * 4*std::abs(cosThetaT2)*std::abs(cosThetaT1) );
        Float jacobian3 = jacobian1;

        VonMisesFisherDistr vMF;
        Float kappa_temp = vMF.convolve(jacobian1*m_kappa1, jacobian2*m_kappa2);
        Float kappa = vMF.convolve(kappa_temp, jacobian3 * m_kappa1) *
                      (4 * std::abs(cosThetaI));
        vMF.setKappa(kappa);

        Normal m = Normal(vMF.sample(sample));

        pdf = vMF.eval(Frame::cosTheta(m));
        /* Prevent potential numerical issues in other stages of the model */
		if (pdf < 1e-20f)
			pdf = 0;

        return m;
	}

	inline Float pdf0(const Vector &m) const {
        return eval0(m);
	}

	inline Float pdf1(const Vector &m, const Vector &wi) const {
        return eval1(m, wi);
	}

	Float G0(const Vector &wi, const Vector &wo, const Vector &m, Float kappa) const {
        Float cosThetaO = dot(wo, m);
        Float cosThetaI = dot(wi, m);

        if (cosThetaO <= 0 || cosThetaI <= 0)
            return 0.0f;

        Float a3 = A3(kappa);
        Float a = 0.25f * (a3+1) * (a3+1),
              b = std::pow(a3, 0.3333333f);

        Float t1 = a * std::cos(b*std::acos(cosThetaO));
        Float t2 = a * std::cos(b*std::acos(cosThetaI));

        Float g1 = a3*cosThetaO / t1;
        Float g2 = a3*cosThetaI / t2;

        return g1 * g2;
	}

	inline Float G1(const Vector &wi, const Vector &wo, const Vector &m) const {
		Vector t1 = refract(wi, Normal(0, 0, 1), m_eta);
        Vector t2 = reflect(-t1, Normal(0, 0, 1));

        return G0(wi, -t1, m, m_kappa1) * G0(-t1, t2, m, m_kappa2) *
               G0(t2, wo, m, m_kappa1);
	}

	/// Return a string representation of the contents of this instance
	std::string toString() const {
		return formatString("LayeredMicrofacetDistribution[kappa1=%f, kappa2=%f, eta=%f]",
			m_kappa1, m_kappa2, m_eta);
	}

protected:
    Float A3(Float kappa) const {
        return 1 / std::tanh(kappa) - 1 / kappa;
    }

	Float m_kappa1, m_kappa2;
    Float m_eta;
};

MTS_NAMESPACE_END

#endif /* __MICROFACET_EXT_H */
