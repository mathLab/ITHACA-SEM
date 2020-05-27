////////////////////////////////////////////////////////////////////////////////
//
//  File: TransfiniteSurface.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: OCC subclases that define transfinite surfaces of 4 curves.
//
////////////////////////////////////////////////////////////////////////////////

#include <Geom_BoundedCurve.hxx>
#include <Geom_BoundedSurface.hxx>
#include <Standard_Type.hxx>
#include <Standard_Version.hxx>
#include <gp_GTrsf2d.hxx>

#include <boost/core/ignore_unused.hpp>
#include <vector>

#if OCC_VERSION_HEX < 0x070000
  #define DEFINE_STANDARD_RTTIEXT(C1,C2) DEFINE_STANDARD_RTTI(C1)
  #define DEFINE_STANDARD_RTTI_INLINE(C1,C2) DEFINE_STANDARD_RTTI(C1)
#endif

class Geom_TransfiniteCurve;
DEFINE_STANDARD_HANDLE(Geom_TransfiniteCurve, Geom_BoundedCurve)
class Geom_TransfiniteSurface;
DEFINE_STANDARD_HANDLE(Geom_TransfiniteSurface, Geom_BoundedSurface)

/**
 * @brief A class to describe an isoline of a transfinite surface.
 *
 * @see Geom_TransfiniteSurface
 */
class Geom_TransfiniteCurve : public Geom_BoundedCurve
{
    /// If m_uvDir is true, then the isoline is defined for a constant value of
    /// U.
    bool m_uvDir;
    /// m_val stores the value of the isoline.
    double m_val;
    /// Vector of edges that define the surface.
    std::vector<Handle(Geom_Curve)> m_edges;
    /// fwd[i] is true if m_edges[u]->Value(0) == m_verts[i].
    std::vector<bool> m_fwd;
    /// Stores pairs that denote the start and end parameters of curve `i`.
    std::vector<std::pair<double, double>> m_clims;
    /// Coordinates of the patch \f$ P_i \f$.
    std::vector<gp_Pnt> m_verts;

    /**
     * @brief Helper function that maps @p val which lies in the interval \f$
     * [0,1] \f$ to the appropriate parameter of curve @p edge, accounting for
     * orientation.
     */
    double Map(int edge, double val) const
    {
        if (m_fwd[edge])
        {
            return m_clims[edge].first + val * (
                m_clims[edge].second - m_clims[edge].first);
        }
        else
        {
            return m_clims[edge].second + val * (
                m_clims[edge].first - m_clims[edge].second);
        }
    }

    DEFINE_STANDARD_RTTIEXT(Geom_TransfiniteCurve, Geom_BoundedCurve)

public:

    /**
     * @brief Construct a isoline of a transfinite surface in either the \f$ u
     * \f$ or \f$ v \f$ directions.
     *
     * @param uvDir    If true, the isoline is for a constant \f$ u \f$
     * @param val      The value of \f$ u \f$ or \f$ v \f$ for the isoline.
     * @param edges    Vector of edges that define the transfinite surface.
     * @param fwd      Orientations of the edges.
     * @param clims    Lower and upper bounds on the edge parametrisation.
     * @param verts    Outer vertices of the transfinite patch.
     */
    Geom_TransfiniteCurve(bool                                   uvDir,
                          double                                 val,
                          std::vector<Handle(Geom_Curve)>        edges,
                          std::vector<bool>                      fwd,
                          std::vector<std::pair<double, double>> clims,
                          std::vector<gp_Pnt>                    verts)
        : m_uvDir(uvDir), m_val(val), m_fwd(fwd), m_clims(clims), m_verts(verts)
    {
        // Take a copy, just to be on the safe side...
        m_edges.resize(edges.size());
        for (int i = 0; i < edges.size(); ++i)
        {
            m_edges[i] = Handle(Geom_Curve)::DownCast(edges[i]->Copy());
        }
    }

    /**
     * @brief Returns the endpoint of the isoline; i.e. if #m_uvDir is true,
     * then given a value #m_val \f$ a \f$ this will be defined by \f$
     * \mathbf{c}_2(1-a) \f$.
     */
    virtual gp_Pnt EndPoint() const override
    {
        return m_uvDir ? m_edges[2]->Value(Map(0, 1-m_val))
            : m_edges[1]->Value(Map(0, m_val));
    }

    /**
     * @brief Returns the starting point of the isoline; i.e. if #m_uvDir is
     * true, then given a value #m_val \f$ a \f$ this will be defined by \f$
     * \mathbf{c}_0(a) \f$.
     */
    virtual gp_Pnt StartPoint() const override
    {
        return m_uvDir ? m_edges[0]->Value(Map(0, m_val))
            : m_edges[3]->Value(Map(0, 1-m_val));
    }

    /**
     * @brief Reverses the direction of the curve. Unimplemented in this class.
     */
    virtual void Reverse() override
    {
        abort();
    }

    /**
     * @brief Returns the reversed parameter for this curve. Unimplemented in
     * this class.
     */
    virtual Standard_Real ReversedParameter (
        const Standard_Real U) const override
    {
        boost::ignore_unused(U);
        abort();
        return 0.0;
    }

    /**
     * @brief Returns the transformed parameter for this curve. Unimplemented in
     * this class.
     */
    virtual Standard_Real TransformedParameter (
        const Standard_Real U, const gp_Trsf& T) const override
    {
        boost::ignore_unused(U, T);
        abort();
        return 0.0;
    }

    /**
     * @brief Returns a parametric transformation for this curve. Unimplemented
     * in this class.
     */
    virtual Standard_Real ParametricTransformation (
        const gp_Trsf& T) const override
    {
        boost::ignore_unused(T);
        abort();
        return 0.0;
    }

    /**
     * @brief Defines the lower range of the parametrisation. Always `0.0` for
     * this class.
     */
    virtual Standard_Real FirstParameter() const override
    {
        return 0.0;
    }

    /**
     * @brief Defines the upper range of the parametrisation. Always `1.0` for
     * this class.
     */
    virtual Standard_Real LastParameter() const override
    {
        return 1.0;
    }

    /**
     * @brief Determines whether this curve is closed; always true for this
     * class.
     */
    virtual Standard_Boolean IsClosed() const override
    {
        return true;
    }

    /**
     * @brief Determines whether this curve is periodic; always false for this
     * class.
     */
    virtual Standard_Boolean IsPeriodic() const override
    {
        return false;
    }

    /**
     * @brief Determines the period of this curve, which is not defined for this
     * class.
     */
    virtual Standard_Real Period() const override
    {
        abort();
        return 0.0;
    }

    /**
     * @brief Returns the continuity of this curve; we only define up to
     * second-order derivatives, so return C^2 continuity.
     */
    virtual GeomAbs_Shape Continuity() const override
    {
        return GeomAbs_C2;
    }

    /**
     * @brief Determines whether the curve is \f$ C^N \f$ continuously
     * differentiable for a given order \f$ N \f$; we only define up to
     * second-order derivatives, so return true only if \f$ 0\leq N\leq 2 \f$.
     */
    virtual Standard_Boolean IsCN (const Standard_Integer N) const override
    {
        if (N > 2 || N < 0)
        {
            return false;
        }
        return true;
    }

    /**
     * @brief Compute a point on the curve @p P given the parametric coordinate
     * @p U.
     *
     * This evaluates the expression from Geom_TransfiniteSurface::D0.
     */
    virtual void D0 (const Standard_Real U1, gp_Pnt& P) const override
    {
        double U = m_uvDir ? m_val : U1;
        double V = m_uvDir ? U1 : m_val;

        gp_XYZ c0 = m_edges[0]->Value(Map(0, U)).XYZ();
        gp_XYZ c1 = m_edges[1]->Value(Map(1, V)).XYZ();
        gp_XYZ c2 = m_edges[2]->Value(Map(2, 1-U)).XYZ();
        gp_XYZ c3 = m_edges[3]->Value(Map(3, 1-V)).XYZ();
        gp_XYZ v0 = m_verts[0].XYZ(), v1 = m_verts[1].XYZ();
        gp_XYZ v2 = m_verts[2].XYZ(), v3 = m_verts[3].XYZ();

        gp_XYZ pnt = (1-U) * c3 + U * c1 + (1-V) * c0 + V * c2 - (
            (1-U) * (1-V) * v0 + U * (1-V) * v1 + U * V * v2 + (1-U) * V * v3);
        P = gp_Pnt(pnt);
    }

    /**
     * @brief Compute a point on the curve @p P and the first order derivative
     * @p V1 given the parametric coordinate @p U.
     *
     * This evaluates the expression from Geom_TransfiniteSurface::D1.
     */
    virtual void D1 (
        const Standard_Real U1, gp_Pnt& P, gp_Vec& V1) const override
    {
        double U = m_uvDir ? m_val : U1;
        double V = m_uvDir ? U1 : m_val;

        // Compute first-order derivatives from OCC curves.
        gp_Pnt t0, t1, t2, t3;
        gp_Vec d0, d1, d2, d3;

        m_edges[0]->D1(Map(0, U), t0, d0);
        m_edges[1]->D1(Map(1, V), t1, d1);
        m_edges[2]->D1(Map(2, 1-U), t2, d2);
        m_edges[3]->D1(Map(3, 1-V), t3, d3);

        // Convert everything to a gp_XYZ so that we can actually do stuff with
        // this.
        gp_XYZ c0 = t0.XYZ(), c1 = t1.XYZ(), c2 = t2.XYZ(), c3 = t3.XYZ();
        gp_XYZ c0p = d0.XYZ(), c1p = d1.XYZ(), c2p = d2.XYZ(), c3p = d3.XYZ();
        gp_XYZ v0 = m_verts[0].XYZ(), v1 = m_verts[1].XYZ();
        gp_XYZ v2 = m_verts[2].XYZ(), v3 = m_verts[3].XYZ();

        gp_XYZ pnt = (1-U) * c3 + U * c1 + (1-V) * c0 + V * c2 - (
            (1-U) * (1-V) * v0 + U * (1-V) * v1 + U * V * v2 + (1-U) * V * v3);
        P = gp_Pnt(pnt);

        // Multiply by chain rule factors.
        c0p *= (m_clims[0].second - m_clims[0].first);
        c1p *= (m_clims[1].second - m_clims[1].first);
        c2p *= (m_clims[2].second - m_clims[2].first);
        c3p *= (m_clims[3].second - m_clims[3].first);
        c0p *= m_fwd[0] ? 1.0 : -1.0;
        c1p *= m_fwd[1] ? 1.0 : -1.0;
        c2p *= m_fwd[2] ? -1.0 : 1.0; // These two edges contain the (1-u) and
        c3p *= m_fwd[3] ? -1.0 : 1.0; // (1-v) terms, so negate here.

        // Analytic derivatives of the transfinite surface.
        if (!m_uvDir)
        {
            gp_XYZ du = -1.0 * c3 + c1 + (1-V) * c0p + V * c2p - (
                -(1-V) * v0 + (1-V) * v1 + V * v2 - V * v3);
            V1 = gp_Vec(du);
        }
        else
        {
            gp_XYZ dv = (1-U) * c3p + U * c1p - c0 + c2 - (
                -(1-U) * v0 - U * v1 + U * v2 + (1-U) * v3);
            V1 = gp_Vec(dv);
        }
    }

    /**
     * @brief Compute a point on the curve @p P, the first order derivative @p
     * V1 and the second order derivative @p V2 given the parametric coordinate
     * @p U.
     *
     * This evaluates the expression from Geom_TransfiniteSurface::D2.
     */
    virtual void D2 (
        const Standard_Real U1, gp_Pnt& P, gp_Vec& V1,
        gp_Vec& V2) const override
    {
        double U = m_uvDir ? m_val : U1;
        double V = m_uvDir ? U1 : m_val;

        gp_Pnt t0, t1, t2, t3;
        gp_Vec d0, d1, d2, d3, dd0, dd1, dd2, dd3;

        m_edges[0]->D2(Map(0, U), t0, d0, dd0);
        m_edges[1]->D2(Map(1, V), t1, d1, dd1);
        m_edges[2]->D2(Map(2, 1-U), t2, d2, dd2);
        m_edges[3]->D2(Map(3, 1-V), t3, d3, dd3);

        gp_XYZ c0 = t0.XYZ(), c1 = t1.XYZ(), c2 = t2.XYZ(), c3 = t3.XYZ();
        gp_XYZ c0p = d0.XYZ(), c1p = d1.XYZ(), c2p = d2.XYZ(), c3p = d3.XYZ();
        gp_XYZ c0pp = dd0.XYZ(), c1pp = dd1.XYZ(), c2pp = dd2.XYZ();
        gp_XYZ c3pp = dd3.XYZ();
        gp_XYZ v0 = m_verts[0].XYZ(), v1 = m_verts[1].XYZ();
        gp_XYZ v2 = m_verts[2].XYZ(), v3 = m_verts[3].XYZ();

        gp_XYZ pnt = (1-U) * c3 + U * c1 + (1-V) * c0 + V * c2 - (
            (1-U) * (1-V) * v0 + U * (1-V) * v1 + U * V * v2 + (1-U) * V * v3);
        P = gp_Pnt(pnt);

        // Multiply by chain rule factors.
        c0p *= (m_clims[0].second - m_clims[0].first);
        c1p *= (m_clims[1].second - m_clims[1].first);
        c2p *= (m_clims[2].second - m_clims[2].first);
        c3p *= (m_clims[3].second - m_clims[3].first);
        c0p *= m_fwd[0] ? 1.0 : -1.0;
        c1p *= m_fwd[1] ? 1.0 : -1.0;
        c2p *= m_fwd[2] ? -1.0 : 1.0; // These two edges are flipped
        c3p *= m_fwd[3] ? -1.0 : 1.0; // so negate here.

        // Second order derivative chain rule factors.
        c0pp *= (m_clims[0].second - m_clims[0].first) *
                (m_clims[0].second - m_clims[0].first);
        c1pp *= (m_clims[1].second - m_clims[1].first) *
                (m_clims[1].second - m_clims[1].first);
        c2pp *= (m_clims[2].second - m_clims[2].first) *
                (m_clims[2].second - m_clims[2].first);
        c3pp *= (m_clims[3].second - m_clims[3].first) *
                (m_clims[3].second - m_clims[3].first);

        // Analytic derivatives of the transfinite surface.
        if (!m_uvDir)
        {
            gp_XYZ du = -1.0 * c3 + c1 + (1-V) * c0p + V * c2p - (
                -(1-V) * v0 + (1-V) * v1 + V * v2 - V * v3);
            gp_XYZ duu = (1-V) * c0pp + V * c2pp;
            V1 = gp_Vec(du);
            V2 = gp_Vec(duu);
        }
        else
        {
            gp_XYZ dv = (1-U) * c3p + U * c1p - c0 + c2 - (
                -(1-U) * v0 - U * v1 + U * v2 + (1-U) * v3);
            gp_XYZ dvv = (1-U) * c3pp + U * c1pp;
            V1 = gp_Vec(dv);
            V2 = gp_Vec(dvv);
        }
    }

    /**
     * @brief Compute a point on the curve @p P, the first order derivative @p
     * V1, the second order derivative @p V2 and third-order derivative @p V3
     * given the parametric coordinate @p U.
     *
     * Since we're lazy and D3 isn't called in the mesh generation pipeline,
     * this isn't implemented in this class (even though it is well-defined).
     */
    virtual void D3 (const Standard_Real U, gp_Pnt& P, gp_Vec& V1, gp_Vec& V2,
                     gp_Vec& V3) const override
    {
        boost::ignore_unused(U, P, V1, V2, V3);
        abort();
    }

    /**
     * @brief Compute the @p N -th order derivative of the curve given the
     * parametric coordinate @p U.
     *
     * Since we're lazy and DN isn't called in the mesh generation pipeline,
     * this isn't implemented in this class (even though it is well-defined).
     */
    virtual gp_Vec DN (
        const Standard_Real U, const Standard_Integer N) const override
    {
        boost::ignore_unused(U, N);
        abort();
        return gp_Vec();
    }

    /**
     * @brief Returns a copy of this geometry.
     */
    virtual Handle(Geom_Geometry) Copy() const override
    {
        Handle(Geom_TransfiniteCurve) tmp = new Geom_TransfiniteCurve(
            m_uvDir, m_val, m_edges, m_fwd, m_clims, m_verts);
        return tmp;
    }

    /**
     * @brief Transform this shape according to the transformation @p
     * T. Unimplemented in this class.
     */
    virtual void Transform (const gp_Trsf& T) override
    {
        abort();

        // Transform our constituent parts.
        for (int i = 0; i < m_edges.size(); ++i)
        {
            m_verts[i].Transform(T);
            m_edges[i]->Transform(T);
            m_clims[i].first =
                m_edges[i]->TransformedParameter(m_clims[i].first, T);
            m_clims[i].second =
                m_edges[i]->TransformedParameter(m_clims[i].second, T);
        }
    }
};

/**
 * @brief A class that defines a transfinite surface defined by four curves.
 *
 * Transfinite surfaces are defined using either three or four curves, and blend
 * the curvature of each in order to define a surface filling. These were
 * introduced by the paper of Gordon & Hall (see Int. J. Num. Meth. Eng., 7 (4)
 * pp 461-477, 1973).  For a surface with parameters \f$ (u,v)\in [0,1]^2 \f$
 * comprised of four curves orientated in a counter-clockwise manner \f$
 * \mathbf{c}_0(u) \f$, \f$ \mathbf{c}_1(v) \f$, \f$ \mathbf{c}_2(u) \f$ and \f$
 * \mathbf{c}_3(v) \f$, with corresponding vertices \f$ P_i \f$ such that \f$
 * \mathbf{c}_i(0) = P_i \f$, we define the surface parametrisation as
 *
 * \f{align*}{
 * \mathbf{S}(u, v) &= (1-u)\mathbf{c}_3(1-v) + u\mathbf{c}_1(v) +
 * (1-v)\mathbf{c}_0(u) + v \mathbf{c}_2(1-u) \\
 * & - [ (1-u)(1-v) P_0 + u(1-v) P_1 + uv P_2 + (1-u)v P_3 ]
 * \f}
 *
 * An extension of this approach allows one to define transfinite surfaces with
 * three curves, but this is not presently supported.
 */
class Geom_TransfiniteSurface : public Geom_BoundedSurface
{
    /// Vector of edges that define the surface.
    std::vector<Handle(Geom_Curve)> m_edges;
    /// fwd[i] is true if m_edges[u]->Value(0) == m_verts[i].
    std::vector<bool> m_fwd;
    /// Stores pairs that denote the start and end parameters of curve `i`.
    std::vector<std::pair<double, double>> m_clims;
    /// Coordinates of the patch \f$ P_i \f$.
    std::vector<gp_Pnt> m_verts;

    /**
     * @brief Helper function that maps @p val which lies in the interval \f$
     * [0,1] \f$ to the appropriate parameter of curve @p edge, accounting for
     * orientation.
     */
    double Map(int edge, double val) const
    {
        if (m_fwd[edge])
        {
            return m_clims[edge].first + val * (
                m_clims[edge].second - m_clims[edge].first);
        }
        else
        {
            return m_clims[edge].second + val * (
                m_clims[edge].first - m_clims[edge].second);
        }
    }

    DEFINE_STANDARD_RTTIEXT(Geom_TransfiniteSurface, Geom_BoundedSurface)

public:
    /**
     * @brief Construct a transfinite surface.
     *
     * @param edges    Vector of edges that define the transfinite surface.
     * @param fwd      Orientations of the edges.
     * @param clims    Lower and upper bounds on the edge parametrisation.
     * @param verts    Outer vertices of the transfinite patch.
     */
    Geom_TransfiniteSurface(std::vector<Handle(Geom_Curve)>        edges,
                            std::vector<bool>                      fwd,
                            std::vector<std::pair<double, double>> clims,
                            std::vector<gp_Pnt>                    verts)
        : m_fwd(fwd), m_clims(clims), m_verts(verts)
    {
        // Take a copy, just to be on the safe side...
        m_edges.resize(edges.size());
        for (int i = 0; i < edges.size(); ++i)
        {
            m_edges[i] = Handle(Geom_Curve)::DownCast(edges[i]->Copy());
        }
    }

    /**
     * @brief Reverses the direction of U. Unimplemented in this class.
     */
    virtual void UReverse() override
    {
        abort();
    }

    virtual Standard_Real UReversedParameter(
        const Standard_Real U) const override
    {
        boost::ignore_unused(U);
        abort();
        return 0.0;
    }

    /**
     * @brief Reverses the direction of V. Unimplemented in this class.
     */
    virtual void VReverse() override
    {
        abort();
    }

    virtual Standard_Real VReversedParameter(
        const Standard_Real V) const override
    {
        boost::ignore_unused(V);
        abort();
        return 0.0;
    }

    /**
     * @brief Transform this shape according to the transformation @p
     * T. Unimplemented in this class.
     */
    virtual void Transform(const gp_Trsf& T) override
    {
        abort();

        // Transform our constituent parts.
        for (int i = 0; i < m_edges.size(); ++i)
        {
            m_verts[i].Transform(T);
            m_edges[i]->Transform(T);
            m_clims[i].first =
                m_edges[i]->TransformedParameter(m_clims[i].first, T);
            m_clims[i].second =
                m_edges[i]->TransformedParameter(m_clims[i].second, T);
        }
    }

    virtual Handle(Geom_Geometry) Copy() const override
    {
        Handle(Geom_TransfiniteSurface) tmp = new Geom_TransfiniteSurface(
            m_edges, m_fwd, m_clims, m_verts);
        return tmp;
    }

    virtual void TransformParameters (
        Standard_Real& U, Standard_Real& V, const gp_Trsf& T) const override
    {
        boost::ignore_unused(U, V, T);
        abort();
    }

    virtual gp_GTrsf2d ParametricTransformation (
        const gp_Trsf& T) const override
    {
        boost::ignore_unused(T);
        abort();
        return gp_GTrsf2d();
    }

    /**
     * @brief Returns the bounds of the surface. For our surface we always have
     * that \f$ U_{\min{}} = V_{\min{}} = 0 \f$ and \f$ U_{\max{}} = V_{\max{}}
     * = 1 \f$.
     */
    virtual void Bounds (Standard_Real& U1, Standard_Real& U2,
                         Standard_Real& V1, Standard_Real& V2) const override
    {
        U1 = V1 = 0.0;
        U2 = V2 = 1.0;
    }

    /**
     * @brief Returns whether the surface is closed in the \f$ u\f$ direction,
     * which is always true for this surface.
     */
    virtual Standard_Boolean IsUClosed() const override
    {
        return true;
    }

    /**
     * @brief Returns whether the surface is closed in the \f$ v\f$ direction,
     * which is always true for this surface.
     */
    virtual Standard_Boolean IsVClosed() const override
    {
        return true;
    }

    /**
     * @brief Returns whether the surface is periodic in the \f$ u\f$ direction,
     * which is always false for this surface.
     */
    virtual Standard_Boolean IsUPeriodic() const override
    {
        return false;
    }

    /**
     * @brief Returns the period of this surface in the \f$ u\f$ direction,
     * which is not defined for this surface.
     */
    virtual Standard_Real UPeriod() const override
    {
        abort();
        return 0.0;
    }

    /**
     * @brief Returns whether the surface is periodic in the \f$ v\f$ direction,
     * which is always false for this surface.
     */
    virtual Standard_Boolean IsVPeriodic() const override
    {
        return false;
    }

    /**
     * @brief Returns the period of this surface in the \f$ v\f$ direction,
     * which is not defined for this surface.
     */
    virtual Standard_Real VPeriod() const override
    {
        abort();
        return 0.0;
    }

    /**
     * @brief Construct an isoline in the \f$ v\f$ direction of the surface
     * according to a fixed parameter @p U.
     */
    virtual Handle(Geom_Curve) UIso (const Standard_Real U) const override
    {
        Handle(Geom_Curve) c = new Geom_TransfiniteCurve(
            true, U, m_edges, m_fwd, m_clims, m_verts);
        return c;
    }

    /**
     * @brief Construct an isoline in the \f$ u\f$ direction of the surface
     * according to a fixed parameter @p V.
     */
    virtual Handle(Geom_Curve) VIso (const Standard_Real V) const override
    {
        Handle(Geom_Curve) c = new Geom_TransfiniteCurve(
            false, V, m_edges, m_fwd, m_clims, m_verts);
        return c;
    }

    /**
     * @brief Returns the continuity of this curve; we only define up to
     * second-order derivatives, so return C^2 continuity.
     */
    virtual GeomAbs_Shape Continuity() const override
    {
        return GeomAbs_C2;
    }

    /**
     * @brief Determines whether the curve is \f$ C^N \f$ continuously
     * differentiable for a given order \f$ N \f$ in the \f$ u\f$ direction; we
     * only define up to second-order derivatives, so return true only if \f$
     * 0\leq N\leq 2 \f$.
     */
    virtual Standard_Boolean IsCNu (const Standard_Integer N) const override
    {
        if (N > 2 || N < 0)
        {
            return false;
        }
        return true;
    }

    /**
     * @brief Determines whether the curve is \f$ C^N \f$ continuously
     * differentiable for a given order \f$ N \f$ in the \f$ v\f$ direction; we
     * only define up to second-order derivatives, so return true only if \f$
     * 0\leq N\leq 2 \f$.
     */
    virtual Standard_Boolean IsCNv (const Standard_Integer N) const override
    {
        if (N > 2 || N < 0)
        {
            return false;
        }
        return true;
    }

    /**
     * @brief Compute a point on the surface @p P given the parametric
     * coordinates @p U and @p V.
     *
     * This evaluates the expression:
     *
     * \f{align*}{
     * \mathbf{S}(u, v) &= (1-u)\mathbf{c}_3(1-v) + u\mathbf{c}_1(v) +
     * (1-v)\mathbf{c}_0(u) + v \mathbf{c}_2(1-u) \\
     * & - [ (1-u)(1-v) P_0 + u(1-v) P_1 + uv P_2 + (1-u)v P_3 ]
     * \f}
     */
    virtual void D0 (
        const Standard_Real U, const Standard_Real V, gp_Pnt& P) const override
    {
        gp_XYZ c0 = m_edges[0]->Value(Map(0, U)).XYZ();
        gp_XYZ c1 = m_edges[1]->Value(Map(1, V)).XYZ();
        gp_XYZ c2 = m_edges[2]->Value(Map(2, 1-U)).XYZ();
        gp_XYZ c3 = m_edges[3]->Value(Map(3, 1-V)).XYZ();
        gp_XYZ v0 = m_verts[0].XYZ(), v1 = m_verts[1].XYZ();
        gp_XYZ v2 = m_verts[2].XYZ(), v3 = m_verts[3].XYZ();

        gp_XYZ pnt = (1-U) * c3 + U * c1 + (1-V) * c0 + V * c2 - (
            (1-U) * (1-V) * v0 + U * (1-V) * v1 + U * V * v2 + (1-U) * V * v3);
        P = gp_Pnt(pnt);
    }

    /**
     * @brief Compute a point on the surface @p P and the first order
     * derivatives @p D1U and @p D1V given the parametric coordinates @p U and
     * @p V.
     *
     * As well as the evaluation of @p P from #D0, this function evaluates the
     * first-order derivatives. Note that since the parameterisations for each
     * curve \f$ i \f$ are of the form \f$ u^* = a_i + b_i u \f$, where \f$ a_i
     * \f$ and \f$ b_i \f$ are the parameter limits defined via #m_clims, we
     * also need to apply the chain rule to account for these extra
     * factors. Finally, if a curve has reversed orientation, the derivative
     * must also be negated.
     *
     * The expression for first order derivatives in the context of the
     * expression in #D0 (which does not account for these terms) is therefore:
     *
     * \f{align*}{
     * \frac{\partial\mathbf{S}}{\partial u} &= -\mathbf{c}_3(1-v) +
     * \mathbf{c}_1(v) + (1-v)\mathbf{c}'_0(u) - v\mathbf{c}_2(1-u) \\
     * & - [ -(1-v)P_0 + (1-v)P_1 + vP_2 - vP_3 ] \\
     * \frac{\partial\mathbf{S}}{\partial v} = & u\mathbf{c}_3'(1-v) +
     * u\mathbf{c}_1'(v) -\mathbf{c}_0(u) +\mathbf{c}_2(1-u) \\
     * & - [ -(1-u)P_0 - u P_1 + u P_2 + (1-u) P_3 ]
     * \f}
     *
     * where \f$ \mathbf{c}'_i \f$ denotes the derivative of curve \f$ i \f$,
     * which can be obtained from the D1 routines of each curve.
     */
    virtual void D1 (const Standard_Real U, const Standard_Real V, gp_Pnt& P,
                     gp_Vec& D1U, gp_Vec& D1V) const override
    {
        // Compute first-order derivatives from OCC curves.
        gp_Pnt t0, t1, t2, t3;
        gp_Vec d0, d1, d2, d3;

        m_edges[0]->D1(Map(0, U), t0, d0);
        m_edges[1]->D1(Map(1, V), t1, d1);
        m_edges[2]->D1(Map(2, 1-U), t2, d2);
        m_edges[3]->D1(Map(3, 1-V), t3, d3);

        // Convert everything to a gp_XYZ so that we can actually do stuff with
        // this.
        gp_XYZ c0 = t0.XYZ(), c1 = t1.XYZ(), c2 = t2.XYZ(), c3 = t3.XYZ();
        gp_XYZ c0p = d0.XYZ(), c1p = d1.XYZ(), c2p = d2.XYZ(), c3p = d3.XYZ();
        gp_XYZ v0 = m_verts[0].XYZ(), v1 = m_verts[1].XYZ();
        gp_XYZ v2 = m_verts[2].XYZ(), v3 = m_verts[3].XYZ();

        gp_XYZ pnt = (1-U) * c3 + U * c1 + (1-V) * c0 + V * c2 - (
            (1-U) * (1-V) * v0 + U * (1-V) * v1 + U * V * v2 + (1-U) * V * v3);
        P = gp_Pnt(pnt);

        // Multiply by chain rule factors.
        c0p *= (m_clims[0].second - m_clims[0].first);
        c1p *= (m_clims[1].second - m_clims[1].first);
        c2p *= (m_clims[2].second - m_clims[2].first);
        c3p *= (m_clims[3].second - m_clims[3].first);
        c0p *= m_fwd[0] ? 1.0 : -1.0;
        c1p *= m_fwd[1] ? 1.0 : -1.0;
        c2p *= m_fwd[2] ? -1.0 : 1.0; // These two edges contain the (1-u) and
        c3p *= m_fwd[3] ? -1.0 : 1.0; // (1-v) terms, so negate here.

        // Analytic derivatives of the transfinite surface.
        gp_XYZ du = -1.0 * c3 + c1 + (1-V) * c0p + V * c2p - (
            -(1-V) * v0 + (1-V) * v1 + V * v2 - V * v3);
        gp_XYZ dv = (1-U) * c3p + U * c1p - c0 + c2 - (
            -(1-U) * v0 - U * v1 + U * v2 + (1-U) * v3);

        D1U = gp_Vec(du);
        D1V = gp_Vec(dv);
    }

    /**
     * @brief Compute a point on the surface @p P and the first order
     * derivatives @p D1U and @p D1V given the parametric coordinates @p U and
     * @p V.
     *
     * As well as the evaluation of @p P, @p D1U and @p D1V from #D1, this
     * function evaluates the second-order derivatives.
     *
     * The expression for sdcond order derivatives in the context of the
     * expression in #D0 (which does not account for the chain rule terms
     * defined in #D1) is therefore:
     *
     * \f{align*}{
     * \frac{\partial^2\mathbf{S}}{\partial v^2} &=
     * (1-u) \mathbf{c}_4''(1-v) + u \mathbf{c}_2''(v) \\
     * \frac{\partial^2\mathbf{S}}{\partial v^2} &=
     * (1-v) \mathbf{c}_1''(u) + v \mathbf{c}_3''(1-u) \\
     * \frac{\partial^2\mathbf{S}}{\partial u\partial v} &= -\mathbf{c}_3'(1-v)
     * + \mathbf{c}_1'(v) - \mathbf{c}_0'(u) + \mathbf{c}_2'(1-u) \\
     * & - [ P_0 - P_1 + P_2 - P_3 ]
     * \f}
     *
     * where \f$ \mathbf{c}'_i \f$ denotes the derivative of curve \f$ i \f$,
     * which can be obtained from the D1 routines of each curve.
     */
    virtual void D2 (const Standard_Real U, const Standard_Real V, gp_Pnt& P,
                     gp_Vec& D1U, gp_Vec& D1V, gp_Vec& D2U, gp_Vec& D2V,
                     gp_Vec& D2UV) const override
    {
        gp_Pnt t0, t1, t2, t3;
        gp_Vec d0, d1, d2, d3, dd0, dd1, dd2, dd3;

        m_edges[0]->D2(Map(0, U), t0, d0, dd0);
        m_edges[1]->D2(Map(1, V), t1, d1, dd1);
        m_edges[2]->D2(Map(2, 1-U), t2, d2, dd2);
        m_edges[3]->D2(Map(3, 1-V), t3, d3, dd3);

        gp_XYZ c0 = t0.XYZ(), c1 = t1.XYZ(), c2 = t2.XYZ(), c3 = t3.XYZ();
        gp_XYZ c0p = d0.XYZ(), c1p = d1.XYZ(), c2p = d2.XYZ(), c3p = d3.XYZ();
        gp_XYZ c0pp = dd0.XYZ(), c1pp = dd1.XYZ(), c2pp = dd2.XYZ();
        gp_XYZ c3pp = dd3.XYZ();
        gp_XYZ v0 = m_verts[0].XYZ(), v1 = m_verts[1].XYZ();
        gp_XYZ v2 = m_verts[2].XYZ(), v3 = m_verts[3].XYZ();

        gp_XYZ pnt = (1-U) * c3 + U * c1 + (1-V) * c0 + V * c2 - (
            (1-U) * (1-V) * v0 + U * (1-V) * v1 + U * V * v2 + (1-U) * V * v3);
        P = gp_Pnt(pnt);

        // Multiply by chain rule factors.
        c0p *= (m_clims[0].second - m_clims[0].first);
        c1p *= (m_clims[1].second - m_clims[1].first);
        c2p *= (m_clims[2].second - m_clims[2].first);
        c3p *= (m_clims[3].second - m_clims[3].first);
        c0p *= m_fwd[0] ? 1.0 : -1.0;
        c1p *= m_fwd[1] ? 1.0 : -1.0;
        c2p *= m_fwd[2] ? -1.0 : 1.0; // These two edges are flipped
        c3p *= m_fwd[3] ? -1.0 : 1.0; // so negate here.

        // Second order derivative chain rule factors.
        c0pp *= (m_clims[0].second - m_clims[0].first) *
                (m_clims[0].second - m_clims[0].first);
        c1pp *= (m_clims[1].second - m_clims[1].first) *
                (m_clims[1].second - m_clims[1].first);
        c2pp *= (m_clims[2].second - m_clims[2].first) *
                (m_clims[2].second - m_clims[2].first);
        c3pp *= (m_clims[3].second - m_clims[3].first) *
                (m_clims[3].second - m_clims[3].first);

        // Analytic derivatives of the transfinite surface.
        gp_XYZ du = -1.0 * c3 + c1 + (1-V) * c0p + V * c2p - (
            -(1-V) * v0 + (1-V) * v1 + V * v2 - V * v3);
        gp_XYZ dv = (1-U) * c3p + U * c1p - c0 + c2 - (
            -(1-U) * v0 - U * v1 + U * v2 + (1-U) * v3);
        gp_XYZ duu = (1-V) * c0pp + V * c2pp;
        gp_XYZ dvv = (1-U) * c3pp + U * c1pp;
        gp_XYZ duv = -1.0 * c3p + c1p - c0p + c2p - (v0 - v1 + v2 - v3);

        D1U = gp_Vec(du);
        D1V = gp_Vec(dv);
        D2U = gp_Vec(duu);
        D2V = gp_Vec(dvv);
        D2UV = gp_Vec(duv);
    }

    /**
     * @brief Compute a point on the curve @p P, alongside its first-, second-
     * and third-order derivatives.
     */
    virtual void D3 (const Standard_Real U, const Standard_Real V, gp_Pnt& P,
                     gp_Vec& D1U, gp_Vec& D1V, gp_Vec& D2U, gp_Vec& D2V,
                     gp_Vec& D2UV, gp_Vec& D3U, gp_Vec& D3V, gp_Vec& D3UUV,
                     gp_Vec& D3UVV) const override
    {
        gp_Pnt t0, t1, t2, t3;
        gp_Vec d0, d1, d2, d3, dd0, dd1, dd2, dd3;
        gp_Vec ddd0, ddd1, ddd2, ddd3;

        m_edges[0]->D3(Map(0, U), t0, d0, dd0, ddd0);
        m_edges[1]->D3(Map(1, V), t1, d1, dd1, ddd1);
        m_edges[2]->D3(Map(2, 1-U), t2, d2, dd2, ddd2);
        m_edges[3]->D3(Map(3, 1-V), t3, d3, dd3, ddd3);

        gp_XYZ c0 = t0.XYZ(), c1 = t1.XYZ(), c2 = t2.XYZ(), c3 = t3.XYZ();
        gp_XYZ c0p = d0.XYZ(), c1p = d1.XYZ(), c2p = d2.XYZ(), c3p = d3.XYZ();
        gp_XYZ c0pp = dd0.XYZ(), c1pp = dd1.XYZ(), c2pp = dd2.XYZ();
        gp_XYZ c3pp = dd3.XYZ();
        gp_XYZ c0ppp = ddd0.XYZ(), c1ppp = ddd1.XYZ(), c2ppp = ddd2.XYZ();
        gp_XYZ c3ppp = ddd3.XYZ();
        gp_XYZ v0 = m_verts[0].XYZ(), v1 = m_verts[1].XYZ();
        gp_XYZ v2 = m_verts[2].XYZ(), v3 = m_verts[3].XYZ();

        gp_XYZ pnt = (1-U) * c3 + U * c1 + (1-V) * c0 + V * c2 - (
            (1-U) * (1-V) * v0 + U * (1-V) * v1 + U * V * v2 + (1-U) * V * v3);
        P = gp_Pnt(pnt);

        // Multiply by chain rule factors.
        c0p *= (m_clims[0].second - m_clims[0].first);
        c1p *= (m_clims[1].second - m_clims[1].first);
        c2p *= (m_clims[2].second - m_clims[2].first);
        c3p *= (m_clims[3].second - m_clims[3].first);
        c0p *= m_fwd[0] ? 1.0 : -1.0;
        c1p *= m_fwd[1] ? 1.0 : -1.0;
        c2p *= m_fwd[2] ? -1.0 : 1.0; // These two edges are flipped
        c3p *= m_fwd[3] ? -1.0 : 1.0; // so negate here.

        // Second order derivative chain rule factors.
        c0pp *= (m_clims[0].second - m_clims[0].first) *
                (m_clims[0].second - m_clims[0].first);
        c1pp *= (m_clims[1].second - m_clims[1].first) *
                (m_clims[1].second - m_clims[1].first);
        c2pp *= (m_clims[2].second - m_clims[2].first) *
                (m_clims[2].second - m_clims[2].first);
        c3pp *= (m_clims[3].second - m_clims[3].first) *
                (m_clims[3].second - m_clims[3].first);
        c0ppp *= pow(m_clims[0].second - m_clims[0].first, 3);
        c1ppp *= pow(m_clims[1].second - m_clims[1].first, 3);
        c2ppp *= pow(m_clims[2].second - m_clims[2].first, 3);
        c3ppp *= pow(m_clims[3].second - m_clims[3].first, 3);
        c0ppp *= m_fwd[0] ? 1.0 : -1.0;
        c1ppp *= m_fwd[1] ? 1.0 : -1.0;
        c2ppp *= m_fwd[2] ? -1.0 : 1.0; // These two edges are flipped
        c3ppp *= m_fwd[3] ? -1.0 : 1.0; // so negate here.

        // Analytic derivatives of the transfinite surface.
        gp_XYZ du = -1.0 * c3 + c1 + (1-V) * c0p + V * c2p - (
            -(1-V) * v0 + (1-V) * v1 + V * v2 - V * v3);
        gp_XYZ dv = (1-U) * c3p + U * c1p - c0 + c2 - (
            -(1-U) * v0 - U * v1 + U * v2 + (1-U) * v3);
        gp_XYZ duu = (1-V) * c0pp + V * c2pp;
        gp_XYZ dvv = (1-U) * c3pp + U * c1pp;
        gp_XYZ duv = -1.0 * c3p + c1p - c0p + c2p - (v0 - v1 + v2 - v3);

        gp_XYZ duuu = (1-V) * c0ppp + V * c2ppp;
        gp_XYZ dvvv = (1-U) * c3ppp + U * c1ppp;
        gp_XYZ duuv = -1.0 * c0pp + c2pp;
        gp_XYZ duvv = -1.0 * c3pp + c1pp;

        D1U = gp_Vec(du);
        D1V = gp_Vec(dv);
        D2U = gp_Vec(duu);
        D2V = gp_Vec(dvv);
        D2UV = gp_Vec(duv);
        D3U = gp_Vec(duuu);
        D3V = gp_Vec(dvvv);
        D3UUV = gp_Vec(duuv);
        D3UVV = gp_Vec(duvv);
    }

    /**
     * @brief Compute the @p Nu and @p Nu order derivative of the curve given
     * the parametric coordinates @p U and @p V.
     *
     * Since we're lazy and DN isn't called in the mesh generation pipeline,
     * this isn't implemented in this class (even though it is well-defined).
     */
    virtual gp_Vec DN (const Standard_Real U,
                       const Standard_Real V,
                       const Standard_Integer Nu,
                       const Standard_Integer Nv) const override
    {
        boost::ignore_unused(U, V, Nu, Nv);
        abort();
        return gp_Vec();
    }
};
