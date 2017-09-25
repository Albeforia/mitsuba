#if !defined(__SPHERICAL_TRIANGLE_H)
#define __SPHERICAL_TRIANGLE_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

// Spherical triangle on unit sphere
struct SphericalTriangle
{
    Vector m_vertices[3];

    SphericalTriangle() {}

    SphericalTriangle(Vector a, Vector b, Vector c)
    {
        m_vertices[0] = normalize(a);
        m_vertices[1] = normalize(b);
        m_vertices[2] = normalize(c);
    }

    Vector &operator[](int index)
    {
        assert(0 <= index && index < 3);
        return m_vertices[index];
    }

    const Vector &operator[](int index) const
    {
        assert(0 <= index && index < 3);
        return m_vertices[index];
    }

    Vector center() const
    {
        return (m_vertices[0] + m_vertices[1] + m_vertices[2]) / 3;
    }

    // https://codegolf.stackexchange.com/questions/63870/spherical-excess-of-a-triangle
    Float excess() const
    {
        Float cosa = dot(m_vertices[1], m_vertices[2]),
              cosb = dot(m_vertices[0], m_vertices[2]),
              cosc = dot(m_vertices[0], m_vertices[1]);

        Float A = acosf((cosa - cosb * cosc) / sqrt((1 - cosb * cosb) * (1 - cosc * cosc))),
              B = acosf((cosb - cosa * cosc) / sqrt((1 - cosa * cosa) * (1 - cosc * cosc))),
              C = acosf((cosc - cosa * cosb) / sqrt((1 - cosa * cosa) * (1 - cosb * cosb)));

        return A + B + C - M_PI;
    }

    void split(std::vector<SphericalTriangle> &children) const
    {
        assert(children.size() >= 4);
        Vector c0 = normalize(m_vertices[1] + m_vertices[2]);
        Vector c1 = normalize(m_vertices[0] + m_vertices[2]);
        Vector c2 = normalize(m_vertices[0] + m_vertices[1]);
        children[0][0] = m_vertices[0];
        children[0][1] = c1;
        children[0][2] = c2;
        children[1][0] = m_vertices[1];
        children[1][1] = c0;
        children[1][2] = c2;
        children[2][0] = m_vertices[2];
        children[2][1] = c0;
        children[2][2] = c1;
        children[3][0] = c0;
        children[3][1] = c1;
        children[3][2] = c2;
    }
};

MTS_NAMESPACE_END

#endif /* __SPHERICAL_TRIANGLE_H */
