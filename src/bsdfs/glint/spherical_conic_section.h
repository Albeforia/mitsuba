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

    bool intersect(const SphericalTriangle &tri) const
    {
        // test if the three spherical arc intersects the spherical conic section
        Vector c, d;
        Float pa, pb, pc, axis;
        for (int i = 0; i < 3; i++)
        {
            c = tri[i] + tri[(i + 1) % 3];
            d = tri[i] - tri[(i + 1) % 3];
            // quadratic form: a*t^2 + b*t + c = 0
            pa = dot(C.preMult(d), d);
            pb = 2 * dot(C.preMult(d), c);
            pc = dot(C.preMult(c), c);
            // true iff there is a root in [-1,1]
            if ((pa - pb + pc) * (pa + pb + pc) < 0)
            {
                return true;
            }
            else
            {
                axis = -pb / (2 * pa);
                if (axis < -1.0f || axis > 2.0f || pa * pc - pb * pb * 0.25f >= 0)
                {
                    continue;
                }
                else
                {
                    return true;
                }
            }
        }
        return false;
    }

    bool contain(const SphericalTriangle &tri) const
    {
        return isInside(tri[0]) && isInside(tri[1]) && isInside(tri[2]);
    }

    int overlap(const SphericalTriangle &tri) const
    {
        if (contain(tri))
        {
            return 2;
        }
        else if (intersect(tri))
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
};

MTS_NAMESPACE_END

#endif /* __SPHERICAL_CONIC_SECTION_H */
