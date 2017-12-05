#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/hw/basicshader.h>
#include <unordered_map>
#include "microfacet_discrete.h"
#include "../ior.h"

#define triangleArea(p0, p1, p2) (std::abs((p0).x * ((p1).y - (p2).y) + (p1).x * ((p2).y - (p0).y) + (p2).x * ((p0).y - (p1).y)) / 2.0f)
#define GAMMA_RADIUS 0.0873f

MTS_NAMESPACE_BEGIN

class Glittery : public BSDF
{
  public:
    Glittery(const Properties &props) : BSDF(props)
    {
        ref<FileResolver> fResolver = Thread::getThread()->getFileResolver();

        m_specularReflectance = new ConstantSpectrumTexture(
            props.getSpectrum("specularReflectance", Spectrum(1.0f)));

        std::string materialName = props.getString("material", "Cu");

        Spectrum intEta, intK;
        if (boost::to_lower_copy(materialName) == "none")
        {
            intEta = Spectrum(0.0f);
            intK = Spectrum(1.0f);
        }
        else
        {
            intEta.fromContinuousSpectrum(InterpolatedSpectrum(
                fResolver->resolve("data/ior/" + materialName + ".eta.spd")));
            intK.fromContinuousSpectrum(InterpolatedSpectrum(
                fResolver->resolve("data/ior/" + materialName + ".k.spd")));
        }

        Float extEta = lookupIOR(props, "extEta", "air");

        m_eta = props.getSpectrum("eta", intEta) / extEta;
        m_k = props.getSpectrum("k", intK) / extEta;

        DiscreteMicrofacetDistribution distr(props);
        m_type = distr.getType();
        m_totalFacets = distr.getTotalFacets();
        m_sampleVisible = distr.getSampleVisible();

        m_alphaU = new ConstantFloatTexture(distr.getAlphaU());
        if (distr.getAlphaU() == distr.getAlphaV())
            m_alphaV = m_alphaU;
        else
            m_alphaV = new ConstantFloatTexture(distr.getAlphaV());

        // four spherical triangles corresponding to quadrants of the hemisphere
        SphericalTriangle tri0(Vector(0, 0, 1), Vector(1, 0, 1e-3), Vector(0, 1, 1e-3));
        SphericalTriangle tri1(Vector(0, 0, 1), Vector(-1, 0, 1e-3), Vector(0, 1, 1e-3));
        SphericalTriangle tri2(Vector(0, 0, 1), Vector(1, 0, 1e-3), Vector(0, -1, 1e-3));
        SphericalTriangle tri3(Vector(0, 0, 1), Vector(-1, 0, 1e-3), Vector(0, -1, 1e-3));
        // "0" is the id of root i.e. the hemisphere
        integrate(distr, tri0, "0", 0);
        integrate(distr, tri1, "0", 1);
        integrate(distr, tri2, "0", 2);
        integrate(distr, tri3, "0", 3);
        SLog(EInfo, "Integration table entries: %d", integrations.size());
    }

    Glittery(Stream *stream, InstanceManager *manager)
        : BSDF(stream, manager)
    {
        m_type = (DiscreteMicrofacetDistribution::EType)stream->readUInt();
        m_sampleVisible = stream->readBool();
        m_alphaU = static_cast<Texture *>(manager->getInstance(stream));
        m_alphaV = static_cast<Texture *>(manager->getInstance(stream));
        m_specularReflectance = static_cast<Texture *>(manager->getInstance(stream));
        m_eta = Spectrum(stream);
        m_k = Spectrum(stream);

        configure();
    }

    void serialize(Stream *stream, InstanceManager *manager) const
    {
        BSDF::serialize(stream, manager);

        stream->writeUInt((uint32_t)m_type);
        stream->writeBool(m_sampleVisible);
        manager->serialize(stream, m_alphaU.get());
        manager->serialize(stream, m_alphaV.get());
        manager->serialize(stream, m_specularReflectance.get());
        m_eta.serialize(stream);
        m_k.serialize(stream);
    }

    void configure()
    {
        unsigned int extraFlags = 0;
        if (m_alphaU != m_alphaV)
            extraFlags |= EAnisotropic;

        if (!m_alphaU->isConstant() || !m_alphaV->isConstant() ||
            !m_specularReflectance->isConstant())
            extraFlags |= ESpatiallyVarying;

        m_components.clear();
        m_components.push_back(EGlossyReflection | EFrontSide | extraFlags);

        /* Verify the input parameters and fix them if necessary */
        m_specularReflectance = ensureEnergyConservation(
            m_specularReflectance, "specularReflectance", 1.0f);

        /*
        m_usesRayDifferentials =
            m_alphaU->usesRayDifferentials() ||
            m_alphaV->usesRayDifferentials() ||
            m_specularReflectance->usesRayDifferentials();
        */
        m_usesRayDifferentials = true;

        BSDF::configure();
    }

    /// Helper function: reflect \c wi with respect to a given surface normal
    inline Vector reflect(const Vector &wi, const Normal &m) const
    {
        return 2 * dot(wi, m) * Vector(m) - wi;
    }

    Spectrum eval(const BSDFSamplingRecord &bRec, EMeasure measure) const
    {
        /* Stop if this component was not requested */
        if (measure != ESolidAngle ||
            Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 ||
            ((bRec.component != -1 && bRec.component != 0) ||
             !(bRec.typeMask & EGlossyReflection)))
            return Spectrum(0.0f);

        /* Calculate the reflection half-vector */
        Vector H = normalize(bRec.wo + bRec.wi);

        /* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
        DiscreteMicrofacetDistribution distr(
            m_type,
            m_alphaU->eval(bRec.its).average(),
            m_alphaV->eval(bRec.its).average(),
            m_totalFacets,
            m_sampleVisible);

        /* Evaluate the microfacet normal distribution */
        Float pixelArea;
        Float D = evaluateD(distr, bRec, pixelArea);

        if (D == 0)
            return Spectrum(0.0f);

        /* Fresnel factor */
        const Spectrum F = fresnelConductorExact(dot(bRec.wi, H), m_eta, m_k) *
                           m_specularReflectance->eval(bRec.its);

        /* Smith's shadow-masking function */
        const Float G = distr.G(bRec.wi, bRec.wo, H);

        /* Calculate the total amount of reflection */
        Float model = dot(bRec.wi, H) * D * G /
                      (pixelArea * (M_PI * (1 - cosf(GAMMA_RADIUS))) * Frame::cosTheta(bRec.wi));

        return F * model;
    }

    Float pdf(const BSDFSamplingRecord &bRec, EMeasure measure) const
    {
        if (measure != ESolidAngle ||
            Frame::cosTheta(bRec.wi) <= 0 ||
            Frame::cosTheta(bRec.wo) <= 0 ||
            ((bRec.component != -1 && bRec.component != 0) ||
             !(bRec.typeMask & EGlossyReflection)))
            return 0.0f;

        /* Calculate the reflection half-vector */
        Vector H = normalize(bRec.wo + bRec.wi);

        /* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
        DiscreteMicrofacetDistribution distr(
            m_type,
            m_alphaU->eval(bRec.its).average(),
            m_alphaV->eval(bRec.its).average(),
            m_totalFacets,
            m_sampleVisible);

        if (m_sampleVisible)
            return distr.eval(H) * distr.smithG1(bRec.wi, H) / (4.0f * Frame::cosTheta(bRec.wi));
        else
            return distr.pdf(bRec.wi, H) / (4 * absDot(bRec.wo, H));
    }

    Spectrum sample(BSDFSamplingRecord &bRec, const Point2 &sample) const
    {
        if (Frame::cosTheta(bRec.wi) < 0 ||
            ((bRec.component != -1 && bRec.component != 0) ||
             !(bRec.typeMask & EGlossyReflection)))
            return Spectrum(0.0f);

        /* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
        DiscreteMicrofacetDistribution distr(
            m_type,
            m_alphaU->eval(bRec.its).average(),
            m_alphaV->eval(bRec.its).average(),
            m_totalFacets,
            m_sampleVisible);

        /* Sample M, the microfacet normal */
        Float pdf;
        Normal m = distr.sample(bRec.wi, sample, pdf);

        if (pdf == 0)
            return Spectrum(0.0f);

        /* Perfect specular reflection based on the microfacet normal */
        bRec.wo = reflect(bRec.wi, m);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EGlossyReflection;

        /* Side check */
        if (Frame::cosTheta(bRec.wo) <= 0)
            return Spectrum(0.0f);

        Spectrum F = fresnelConductorExact(dot(bRec.wi, m),
                                           m_eta, m_k) *
                     m_specularReflectance->eval(bRec.its);

        Float pixelArea;
        Float D = evaluateD(distr, bRec, pixelArea);

        Float weight;
        if (m_sampleVisible)
        {
            weight = distr.smithG1(bRec.wo, m);
        }
        else
        {
            auto iDotm = dot(bRec.wi, m);
            weight = D * distr.G(bRec.wi, bRec.wo, m) * iDotm * iDotm /
                     (pdf * pixelArea * (M_PI * (1 - cosf(GAMMA_RADIUS))) * Frame::cosTheta(bRec.wi));
        }

        return F * weight;
    }

    Spectrum sample(BSDFSamplingRecord &bRec, Float &pdf, const Point2 &sample) const
    {
        if (Frame::cosTheta(bRec.wi) < 0 ||
            ((bRec.component != -1 && bRec.component != 0) ||
             !(bRec.typeMask & EGlossyReflection)))
            return Spectrum(0.0f);

        /* Construct the microfacet distribution matching the
		   roughness values at the current surface position. */
        DiscreteMicrofacetDistribution distr(
            m_type,
            m_alphaU->eval(bRec.its).average(),
            m_alphaV->eval(bRec.its).average(),
            m_totalFacets,
            m_sampleVisible);

        /* Sample M, the microfacet normal */
        Normal m = distr.sample(bRec.wi, sample, pdf);

        if (pdf == 0)
            return Spectrum(0.0f);

        /* Perfect specular reflection based on the microfacet normal */
        bRec.wo = reflect(bRec.wi, m);
        bRec.eta = 1.0f;
        bRec.sampledComponent = 0;
        bRec.sampledType = EGlossyReflection;

        /* Side check */
        if (Frame::cosTheta(bRec.wo) <= 0)
            return Spectrum(0.0f);

        Spectrum F = fresnelConductorExact(dot(bRec.wi, m),
                                           m_eta, m_k) *
                     m_specularReflectance->eval(bRec.its);

        Float pixelArea;
        Float D = evaluateD(distr, bRec, pixelArea);

        Float weight;
        if (m_sampleVisible)
        {
            weight = distr.smithG1(bRec.wo, m);
        }
        else
        {
            auto iDotm = dot(bRec.wi, m);
            weight = D * distr.G(bRec.wi, bRec.wo, m) * iDotm * iDotm /
                     (pdf * pixelArea * (M_PI * (1 - cosf(GAMMA_RADIUS))) * Frame::cosTheta(bRec.wi));
        }

        /* Jacobian of the half-direction mapping */
        pdf /= 4.0f * dot(bRec.wo, m);

        return F * weight;
    }

    void addChild(const std::string &name, ConfigurableObject *child)
    {
        if (child->getClass()->derivesFrom(MTS_CLASS(Texture)))
        {
            if (name == "alpha")
                m_alphaU = m_alphaV = static_cast<Texture *>(child);
            else if (name == "alphaU")
                m_alphaU = static_cast<Texture *>(child);
            else if (name == "alphaV")
                m_alphaV = static_cast<Texture *>(child);
            else if (name == "specularReflectance")
                m_specularReflectance = static_cast<Texture *>(child);
            else
                BSDF::addChild(name, child);
        }
        else
        {
            BSDF::addChild(name, child);
        }
    }

    Float getRoughness(const Intersection &its, int component) const
    {
        return 0.5f * (m_alphaU->eval(its).average() + m_alphaV->eval(its).average());
    }

    std::string toString() const
    {
        std::ostringstream oss;
        oss << "Glittery[" << endl
            << "  id = \"" << getID() << "\"," << endl
            << "  distribution = " << DiscreteMicrofacetDistribution::distributionName(m_type) << "," << endl
            << "  sampleVisible = " << m_sampleVisible << "," << endl
            << "  alphaU = " << indent(m_alphaU->toString()) << "," << endl
            << "  alphaV = " << indent(m_alphaV->toString()) << "," << endl
            << "  specularReflectance = " << indent(m_specularReflectance->toString()) << "," << endl
            << "  eta = " << m_eta.toString() << "," << endl
            << "  k = " << m_k.toString() << endl
            << "]";
        return oss.str();
    }

    Shader *createShader(Renderer *renderer) const;

    MTS_DECLARE_CLASS()
  private:
    DiscreteMicrofacetDistribution::EType m_type;
    ref<Texture> m_specularReflectance;
    ref<Texture> m_alphaU, m_alphaV;
    uint32_t m_totalFacets;
    bool m_sampleVisible;
    Spectrum m_eta, m_k;

    // <spherical triangle id, pre-computed integration value>
    std::unordered_map<std::string, Float> integrations;

    // [Algorithm 2]
    Float integrate(const DiscreteMicrofacetDistribution &distr, const SphericalTriangle &tri,
                    const std::string &parentID, int splitNumber)
    {
        auto id = parentID + std::to_string(splitNumber);

        Float excess = tri.excess();
        Float rule1 = excess * distr.eval(tri.center());
        Float rule2 = excess * (distr.eval(tri[0]) + distr.eval(tri[1]) + distr.eval(tri[2])) / 3;
        Float error = std::abs(rule1 - rule2);
        if (error < 1e-5 || error / rule2 < 1e-5)
        {
            // we also store the result for 'uniform' triangles
            // NOTE this will take up much more memory
            integrations[id] = rule2;
            return rule2;
        }
        else
        {
            Float rule3 = 0;
            auto children = tri.split();
            for (int i = 0; i < 4; i++)
            {
                rule3 += integrate(distr, children[i], id, i);
            }
            integrations[id] = rule3;
            return rule3;
        }
    }

    // get four points of a (2D)parallelogram in counterclockwise order
    Parallelogram extentsToPoint(Vector2 center, Vector2 extentU, Vector2 extentV) const
    {
        Parallelogram points;
        points[1] = center;
        points[3] = center + extentU + extentV;
        auto winding = cross(Vector3(extentU.x, extentU.y, 0), Vector3(extentV.x, extentV.y, 0));
        if (winding.z < 0)
        {
            points[0] = center + extentU;
            points[2] = center + extentV;
        }
        else
        {
            points[0] = center + extentV;
            points[2] = center + extentU;
        }
        return points;
    }

    // evaluate the microfacet normal distribution
    Float evaluateD(const DiscreteMicrofacetDistribution &distr, const BSDFSamplingRecord &bRec,
                    Float &pixelArea) const
    {
        Float D;

        /* Calculate the reflection half-vector */
        Vector H = normalize(bRec.wo + bRec.wi);

        // pixel footprint in texture space
        const auto &its = bRec.its;
        Vector2 center(its.uv.x, its.uv.y);
        Vector2 extentU(its.dudx, its.dvdx);
        Vector2 extentV(its.dudy, its.dvdy);
        auto pixel = extentsToPoint(center, extentU, extentV);
        pixelArea = 2 * triangleArea(pixel[0], pixel[2], pixel[1]);
        // no differentials are specified, reverts back to the smooth case
        if (extentU.isZero() || extentV.isZero())
        {
            D = distr.eval(H);
        }
        else
        {
            size_t sampleCount = bRec.sampler == nullptr ? 1 : bRec.sampler->getSampleCount();
            SphericalConicSection scs(bRec.wi, bRec.wo, GAMMA_RADIUS);
            D = distr.eval(H, pixel, scs, integrations, sampleCount);
        }

        return D;
    }

    std::string to_string(const Parallelogram &paral) const
    {
        std::stringstream ss;
        ss << "[" << paral[0].x << "," << paral[0].y << "] ";
        ss << "[" << paral[1].x << "," << paral[1].y << "] ";
        ss << "[" << paral[2].x << "," << paral[2].y << "] ";
        ss << "[" << paral[3].x << "," << paral[3].y << "]";
        return ss.str();
    }
};

/**
 * GLSL port of the rough conductor shader. This version is much more
 * approximate -- it only supports the Ashikhmin-Shirley distribution,
 * does everything in RGB, and it uses the Schlick approximation to the
 * Fresnel reflectance of conductors. When the roughness is lower than
 * \alpha < 0.2, the shader clamps it to 0.2 so that it will still perform
 * reasonably well in a VPL-based preview.
 */
class GlitteryShader : public Shader
{
  public:
    GlitteryShader(Renderer *renderer, const Texture *specularReflectance,
                   const Texture *alphaU, const Texture *alphaV, const Spectrum &eta,
                   const Spectrum &k) : Shader(renderer, EBSDFShader),
                                        m_specularReflectance(specularReflectance), m_alphaU(alphaU), m_alphaV(alphaV)
    {
        m_specularReflectanceShader = renderer->registerShaderForResource(m_specularReflectance.get());
        m_alphaUShader = renderer->registerShaderForResource(m_alphaU.get());
        m_alphaVShader = renderer->registerShaderForResource(m_alphaV.get());

        /* Compute the reflectance at perpendicular incidence */
        m_R0 = fresnelConductorExact(1.0f, eta, k);
    }

    bool isComplete() const
    {
        return m_specularReflectanceShader.get() != NULL &&
               m_alphaUShader.get() != NULL &&
               m_alphaVShader.get() != NULL;
    }

    void putDependencies(std::vector<Shader *> &deps)
    {
        deps.push_back(m_specularReflectanceShader.get());
        deps.push_back(m_alphaUShader.get());
        deps.push_back(m_alphaVShader.get());
    }

    void cleanup(Renderer *renderer)
    {
        renderer->unregisterShaderForResource(m_specularReflectance.get());
        renderer->unregisterShaderForResource(m_alphaU.get());
        renderer->unregisterShaderForResource(m_alphaV.get());
    }

    void resolve(const GPUProgram *program, const std::string &evalName, std::vector<int> &parameterIDs) const
    {
        parameterIDs.push_back(program->getParameterID(evalName + "_R0", false));
    }

    void bind(GPUProgram *program, const std::vector<int> &parameterIDs, int &textureUnitOffset) const
    {
        program->setParameter(parameterIDs[0], m_R0);
    }

    void generateCode(std::ostringstream &oss,
                      const std::string &evalName,
                      const std::vector<std::string> &depNames) const
    {
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
            << "   float D = " << evalName << "_D(H, alphaU, alphaV)"
            << ";" << endl
            << "   float G = " << evalName << "_G(H, wi, wo);" << endl
            << "   vec3 F = " << evalName << "_schlick(1-dot(wi, H));" << endl
            << "   return reflectance * F * (D * G / (4*cosTheta(wi)));" << endl
            << "}" << endl
            << endl
            << "vec3 " << evalName << "_diffuse(vec2 uv, vec3 wi, vec3 wo) {" << endl
            << "    if (cosTheta(wi) < 0.0 || cosTheta(wo) < 0.0)" << endl
            << "    	return vec3(0.0);" << endl
            << "    return " << evalName << "_R0 * inv_pi * inv_pi * cosTheta(wo);" << endl
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

Shader *Glittery::createShader(Renderer *renderer) const
{
    return new GlitteryShader(renderer,
                              m_specularReflectance.get(), m_alphaU.get(), m_alphaV.get(), m_eta, m_k);
}

MTS_IMPLEMENT_CLASS(GlitteryShader, false, Shader)
MTS_IMPLEMENT_CLASS_S(Glittery, false, BSDF)
MTS_EXPORT_PLUGIN(Glittery, "Glittery BRDF");
MTS_NAMESPACE_END
