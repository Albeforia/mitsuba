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

#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include "microfacet.h"
#include "microfacet_ext.h"
#include "ior.h"

MTS_NAMESPACE_BEGIN

class LayeredConductor : public BSDF {
public:
	LayeredConductor(const Properties &props) : BSDF(props) {
		ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

		m_specularReflectance = new ConstantSpectrumTexture(
			props.getSpectrum("specularReflectance", Spectrum(1.0f)));

		std::string materialName = props.getString("material", "Cu");

		Spectrum intEta, intK;
		if (boost::to_lower_copy(materialName) == "none") {
			intEta = Spectrum(0.0f);
			intK = Spectrum(1.0f);
		} else {
			intEta.fromContinuousSpectrum(InterpolatedSpectrum(
				fResolver->resolve("data/ior/" + materialName + ".eta.spd")));
			intK.fromContinuousSpectrum(InterpolatedSpectrum(
				fResolver->resolve("data/ior/" + materialName + ".k.spd")));
		}

        m_layer_eta = props.getFloat("layerEta", 1.5f);

		// Float extEta = lookupIOR(props, "extEta", "air");
        Float extEta = m_layer_eta;

		m_eta = props.getSpectrum("eta", intEta) / extEta;
		m_k   = props.getSpectrum("k", intK) / extEta;

        m_kappa1 = props.getFloat("kappa1", 100.0f);
        m_kappa2 = props.getFloat("kappa2", 100.0f);
	}

	LayeredConductor(Stream *stream, InstanceManager *manager)
	 : BSDF(stream, manager) {
		m_kappa1 = stream->readFloat();
		m_kappa2 = stream->readFloat();
		m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
		m_eta = Spectrum(stream);
		m_k = Spectrum(stream);
		m_layer_eta = stream->readFloat();

		configure();
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		BSDF::serialize(stream, manager);

        stream->writeFloat(m_kappa1);
        stream->writeFloat(m_kappa2);
		manager->serialize(stream, m_specularReflectance.get());
		m_eta.serialize(stream);
		m_k.serialize(stream);
        stream->writeFloat(m_layer_eta);
	}

	void configure() {
		m_components.clear();
		m_components.push_back(EGlossyReflection | EFrontSide);

		/* Verify the input parameters and fix them if necessary */
		m_specularReflectance = ensureEnergyConservation(
			m_specularReflectance, "specularReflectance", 1.0f);

		m_usesRayDifferentials = m_specularReflectance->usesRayDifferentials();

		BSDF::configure();
	}

	/// Helper function: reflect \c wi with respect to a given surface normal
	inline Vector reflect(const Vector &wi, const Normal &m) const {
		return 2 * dot(wi, m) * Vector(m) - wi;
	}

	Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		/* Stop if this component was not requested */
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo+bRec.wi);

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		LayeredMicrofacetDistribution distr(
            m_kappa1, m_kappa2, m_layer_eta
		);

        /* surface reflection part */
        Spectrum value0(0.0f);

        const Float D0 = distr.eval0(H);
        Float F0 = fresnelDielectricExt(Frame::cosTheta(bRec.wi), m_layer_eta);
		if (D0 > 0) {

            const Float G0 = distr.G0(bRec.wi, bRec.wo, H, m_kappa1);

            value0 = D0 * G0 / (4.0f * Frame::cosTheta(bRec.wi)) * Spectrum(F0);
        }

        /* one bounce part */
        Spectrum value1(0.0f);

        const Float D1 = distr.eval1(H, bRec.wi);
		if (D1 > 0) {
            Vector t1 = refract(bRec.wi, Normal(0, 0, 1), m_layer_eta);
            Vector t2 = reflect(-t1, Normal(0, 0, 1));

            const Spectrum F1 = (1.0f - F0) *
                                ( fresnelConductorExact(Frame::cosTheta(-t1), m_eta, m_k) *
                                m_specularReflectance->eval(bRec.its) ) *
                                (1.0f - fresnelDielectricExt(Frame::cosTheta(-t2), m_layer_eta));

            const Float G1 = distr.G1(bRec.wi, bRec.wo, H);

            value1 = D1 * G1 / (4.0f * Frame::cosTheta(bRec.wi)) * F1;
        }

		return value0 + value1;
	}

	Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const {
		if (measure != ESolidAngle ||
			Frame::cosTheta(bRec.wi) <= 0 ||
			Frame::cosTheta(bRec.wo) <= 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return 0.0f;

		/* Calculate the reflection half-vector */
		Vector H = normalize(bRec.wo+bRec.wi);

		/* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
		LayeredMicrofacetDistribution distr(
            m_kappa1, m_kappa2, m_layer_eta
		);

        /* surface reflection part */

        Float pdf0 = distr.pdf0(H);

		Float F0 = fresnelDielectricExt(Frame::cosTheta(bRec.wi), m_layer_eta);

        pdf0 *= F0;
        pdf0 *= 1.0f / (4.0f * dot(bRec.wo, H));    // Jacobian

        /* one bounce part */

        Float pdf1 = distr.pdf1(H, bRec.wi);

        pdf1 *= 1.0f - F0;
        pdf1 *= 1.0f / (4.0f * dot(bRec.wo, H));    // Jacobian

		return pdf0 + pdf1;
	}

	Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

        LayeredMicrofacetDistribution distr(
            m_kappa1, m_kappa2, m_layer_eta
		);

        Float cosThetaT;
		Float F0 = fresnelDielectricExt(Frame::cosTheta(bRec.wi), cosThetaT, m_layer_eta);
        Float pdf = 1.0f;

        if (bRec.sampler->next1D() > F0) {
            // subsurface
            pdf *= 1.0f - F0;

            Float pdf1;
            const Normal m = distr.sample1(bRec.wi, sample, pdf1);
            if (pdf1 == 0)
			    return Spectrum(0.0f);
            pdf *= pdf1;

            bRec.wo = reflect(bRec.wi, m);
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);

            pdf *= 1.0f / (4.0f * dot(bRec.wo, m));    // Jacobian

            Vector t1 = refract(bRec.wi, Normal(0, 0, 1), m_layer_eta);
            Vector t2 = reflect(-t1, Normal(0, 0, 1));

            const Spectrum F1 = (1.0f - F0) *
                                ( fresnelConductorExact(Frame::cosTheta(-t1), m_eta, m_k) *
                                m_specularReflectance->eval(bRec.its) ) *
                                (1.0f - fresnelDielectricExt(Frame::cosTheta(-t2), m_layer_eta));

            const Float G1 = distr.G1(bRec.wi, bRec.wo, m);
            const Spectrum value1 = distr.eval1(m, bRec.wi) * G1 / (4.0f * Frame::cosTheta(bRec.wi)) * F1;

            return value1 / pdf;
        } else {
            // surface
            pdf *= F0;

            Float pdf0;
            const Normal m = distr.sample0(bRec.wi, sample, pdf0);
            if (pdf0 == 0)
			    return Spectrum(0.0f);
            pdf *= pdf0;

            bRec.wo = reflect(bRec.wi, m);
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);

            pdf *= 1.0f / (4.0f * dot(bRec.wo, m));    // Jacobian

            const Float G0 = distr.G0(bRec.wi, bRec.wo, m, m_kappa1);
            const Float value0 = distr.eval0(m) * G0 / (4.0f * Frame::cosTheta(bRec.wi)) * F0;

            return Spectrum(value0) / pdf;
        }
	}

	Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const {
		if (Frame::cosTheta(bRec.wi) < 0 ||
			((bRec.component != -1 && bRec.component != 0) ||
			!(bRec.typeMask & EGlossyReflection)))
			return Spectrum(0.0f);

        LayeredMicrofacetDistribution distr(
            m_kappa1, m_kappa2, m_layer_eta
		);

        Float cosThetaT;
		Float F0 = fresnelDielectricExt(Frame::cosTheta(bRec.wi), cosThetaT, m_layer_eta);
        pdf = 1.0f;

        if (bRec.sampler->next1D() > F0) {
            // subsurface
            pdf *= 1.0f - F0;

            Float pdf1;
            const Normal m = distr.sample1(bRec.wi, sample, pdf1);
            if (pdf1 == 0)
			    return Spectrum(0.0f);
            pdf *= pdf1;

            bRec.wo = reflect(bRec.wi, m);
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);

            pdf *= 1.0f / (4.0f * dot(bRec.wo, m));    // Jacobian

            Vector t1 = refract(bRec.wi, Normal(0, 0, 1), m_layer_eta);
            Vector t2 = reflect(-t1, Normal(0, 0, 1));

            const Spectrum F1 = (1.0f - F0) *
                                ( fresnelConductorExact(Frame::cosTheta(-t1), m_eta, m_k) *
                                m_specularReflectance->eval(bRec.its) ) *
                                (1.0f - fresnelDielectricExt(Frame::cosTheta(-t2), m_layer_eta));

            const Float G1 = distr.G1(bRec.wi, bRec.wo, m);
            const Spectrum value1 = distr.eval1(m, bRec.wi) * G1 / (4.0f * Frame::cosTheta(bRec.wi)) * F1;

            return value1 / pdf;
        } else {
            // surface
            pdf *= F0;

            Float pdf0;
            const Normal m = distr.sample0(bRec.wi, sample, pdf0);
            if (pdf0 == 0)
			    return Spectrum(0.0f);
            pdf *= pdf0;

            bRec.wo = reflect(bRec.wi, m);
			bRec.eta = 1.0f;
			bRec.sampledComponent = 0;
			bRec.sampledType = EGlossyReflection;

			/* Side check */
			if (Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0)
				return Spectrum(0.0f);

            pdf *= 1.0f / (4.0f * dot(bRec.wo, m));    // Jacobian

            const Float G0 = distr.G0(bRec.wi, bRec.wo, m, m_kappa1);
            const Float value0 = distr.eval0(m) * G0 / (4.0f * Frame::cosTheta(bRec.wi)) * F0;

            return Spectrum(value0) / pdf;
        }
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(Texture))) {
			if (name == "specularReflectance")
				m_specularReflectance = static_cast<Texture *>(child);
			else
				BSDF::addChild(name, child);
		} else {
			BSDF::addChild(name, child);
		}
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "LayeredConductor[" << endl
			<< "  id = \"" << getID() << "\"," << endl
			<< "  kappa1 = " << m_kappa1 << "," << endl
			<< "  kappa2 = " << m_kappa2 << "," << endl
			<< "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
			<< "  eta = " << m_eta.toString() << "," << endl
			<< "  k = " << m_k.toString() << endl
			<< "]";
		return oss.str();
	}

	Shader *createShader(Renderer *renderer) const;

	MTS_DECLARE_CLASS()
private:

	ref<Texture> m_specularReflectance;
	Float m_kappa1, m_kappa2;
	Spectrum m_eta, m_k;
    Float m_layer_eta;
};

/**
 * GLSL port of the rough conductor shader. This version is much more
 * approximate -- it only supports the Ashikhmin-Shirley distribution,
 * does everything in RGB, and it uses the Schlick approximation to the
 * Fresnel reflectance of conductors. When the roughness is lower than
 * \alpha < 0.2, the shader clamps it to 0.2 so that it will still perform
 * reasonably well in a VPL-based preview.
 */
class LayeredConductorShader : public Shader {
public:
	LayeredConductorShader(Renderer *renderer, const Texture *specularReflectance,
			const Spectrum &eta, const Spectrum &k) : Shader(renderer, EBSDFShader),
			m_specularReflectance(specularReflectance) {
		m_alphaU = new ConstantFloatTexture(0.1f);
		m_alphaV = new ConstantFloatTexture(0.1f);

		m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
		m_alphaUShader = renderer->registerShaderForResource(m_alphaU.get());
		m_alphaVShader = renderer->registerShaderForResource(m_alphaV.get());

		/* Compute the reflectance at perpendicular incidence */
		m_R0 = fresnelConductorExact(1.0f, eta, k);
	}

	bool isComplete() const {
		return m_specularReflectanceShader.get() != NULL &&
			   m_alphaUShader.get() != NULL &&
			   m_alphaVShader.get() != NULL;
	}

	void putDependencies(std::vector<Shader *> &deps) {
		deps.push_back(m_specularReflectanceShader.get());
		deps.push_back(m_alphaUShader.get());
		deps.push_back(m_alphaVShader.get());
	}

	void cleanup(Renderer *renderer) {
		renderer->unregisterShaderForResource(m_specularReflectance.get());
		renderer->unregisterShaderForResource(m_alphaU.get());
		renderer->unregisterShaderForResource(m_alphaV.get());
	}

	void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const {
		parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
	}

	void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const {
		program->setParameter(parameterIDs[0], m_R0);
	}

	void generateCode(std::ostringstream &oss,
			const std::string &evalName,
			const std::vector<std::string> &depNames) const {
		oss << "uniform vec3 " << evalName << "_R0;" << endl
			<< endl
			<< "float " << evalName << "_D(vec3 m, float alphaU, float alphaV) {" << endl
			<< "    float ct = cosTheta(m), ds = 1-ct*ct;" << endl
			<< "    if (ds <= 0.0)" << endl
			<< "        return 0.0f;" << endl
			<< "    alphaU = 2 / (alphaU * alphaU) - 2;" << endl
			<< "    alphaV = 2 / (alphaV * alphaV) - 2;" << endl
			<< "    float exponent = (alphaU*m.x*m.x + alphaV*m.y*m.y)/ds;" << endl
			<< "    return sqrt((alphaU+2) * (alphaV+2)) * 0.15915 * pow(ct, exponent);" << endl
			<< "}" << endl
			<< endl
			<< "float " << evalName << "_G(vec3 m, vec3 wi, vec3 wo) {" << endl
			<< "    if ((dot(wi, m) * cosTheta(wi)) <= 0 || " << endl
			<< "        (dot(wo, m) * cosTheta(wo)) <= 0)" << endl
			<< "        return 0.0;" << endl
			<< "    float nDotM = cosTheta(m);" << endl
			<< "    return min(1.0, min(" << endl
			<< "        abs(2 * nDotM * cosTheta(wo) / dot(wo, m))," << endl
			<< "        abs(2 * nDotM * cosTheta(wi) / dot(wi, m))));" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_schlick(float ct) {" << endl
			<< "    float ctSqr = ct*ct, ct5 = ctSqr*ctSqr*ct;" << endl
			<< "    return " << evalName << "_R0 + (vec3(1.0) - " << evalName << "_R0) * ct5;" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "   if (cosTheta(wi) <= 0 || cosTheta(wo) <= 0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "   vec3 H = normalize(wi + wo);" << endl
			<< "   vec3 reflectance = " << depNames[0] << "(uv);" << endl
			<< "   float alphaU = max(0.2, " << depNames[1] << "(uv).r);" << endl
			<< "   float alphaV = max(0.2, " << depNames[2] << "(uv).r);" << endl
			<< "   float D = " << evalName << "_D(H, alphaU, alphaV)" << ";" << endl
			<< "   float G = " << evalName << "_G(H, wi, wo);" << endl
			<< "   vec3 F = " << evalName << "_schlick(1-dot(wi, H));" << endl
			<< "   return reflectance * F * (D * G / (4*cosTheta(wi)));" << endl
			<< "}" << endl
			<< endl
			<< "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
			<< "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
			<< "    	return vec3(0.0);" << endl
			<< "    return " << evalName << "_R0 * inv_pi * inv_pi * cosTheta(wo);"<< endl
			<< "}" << endl;
	}
	MTS_DECLARE_CLASS()
private:
	ref<const Texture> m_specularReflectance;
	ref<const Texture> m_alphaU;
	ref<const Texture> m_alphaV;
	ref<Shader> m_specularReflectanceShader;
	ref<Shader> m_alphaUShader;
	ref<Shader> m_alphaVShader;
	Spectrum m_R0;
};

Shader *LayeredConductor::createShader(Renderer *renderer) const {
	return new LayeredConductorShader(renderer,
		m_specularReflectance.get(), m_eta, m_k);
}

MTS_IMPLEMENT_CLASS(LayeredConductorShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(LayeredConductor, false, BSDF)
MTS_EXPORT_PLUGIN(LayeredConductor, "Layered conductor BRDF");
MTS_NAMESPACE_END
