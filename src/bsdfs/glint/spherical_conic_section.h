#if !defined(__SPHERICAL_CONIC_SECTION_H)
#define __SPHERICAL_CONIC_SECTION_H

#include <mitsuba/mitsuba.h>
#include "spherical_triangle.h"

MTS_NAMESPACE_BEGIN

struct SphericalConicSection
{
    Vector m_wi, m_wo;
    Float m_gamma;
    Matrix3x3 C;

    SphericalConicSection(Vector wi, Vector wo, Float gamma) : m_wi(wi), m_wo(wo), m_gamma(gamma)
    {
        auto x = normalize(cross(m_wi, m_wo)), y = normalize(m_wi - m_wo), z = normalize(m_wi + m_wo);
        Float lambda1 = (dot(m_wi, m_wo) + cosf(m_gamma)) / (1 - cosf(m_gamma)),
              lambda2 = 1 / (tanf(m_gamma / 2) * tanf(m_gamma / 2));
        Matrix3x3 Q(x, y, z), Qt,
            A(lambda1, 0.0f, 0.0f,
              0.0f, lambda2, 0.0f,
              0.0f, 0.0f, -1.0f);
        Q.transpose(Qt);
        C = Q;
        C *= A;
        C *= Qt;
    }

    // mT*C*m <= 0 assuming 'm' is normalized
    bool isInside(const Vector &m) const
    {
        return dot(C.preMult(m), m) <= 0;
    }
};

MTS_NAMESPACE_END

#endif /* __SPHERICAL_CONIC_SECTION_H */
