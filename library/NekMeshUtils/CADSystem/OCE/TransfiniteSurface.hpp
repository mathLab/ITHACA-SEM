////////////////////////////////////////////////////////////////////////////////
//
//  File: TransfiniteSurface.cpp
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
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <Geom_BoundedCurve.hxx>
#include <Geom_BoundedSurface.hxx>
#include <gp_GTrsf2d.hxx>

class Geom_TransfiniteCurve : public Geom_BoundedCurve
{
    bool m_uvDir;
    double m_val;
    std::vector<Handle(Geom_Curve)> m_edges;
    std::vector<bool> m_fwd;
    std::vector<std::pair<double, double>> m_clims;
    std::vector<gp_Pnt> m_verts;

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

public:
    Geom_TransfiniteCurve(bool                                   uvDir,
                          double                                 val,
                          std::vector<Handle(Geom_Curve)>        edges,
                          std::vector<bool>                      fwd,
                          std::vector<std::pair<double, double>> clims,
                          std::vector<gp_Pnt>                    verts)
        : m_uvDir(uvDir), m_val(val), m_edges(edges), m_fwd(fwd),
          m_clims(clims), m_verts(verts)
    {
    }

    virtual gp_Pnt EndPoint() const
    {
        return m_uvDir ? m_edges[2]->Value(Map(0, m_val))
            : m_edges[1]->Value(Map(0, m_val));
    }

    virtual gp_Pnt StartPoint() const
    {
        return m_uvDir ? m_edges[0]->Value(Map(0, m_val))
            : m_edges[3]->Value(Map(0, m_val));
    }

    virtual void Reverse()
    {
        abort();
    }

    virtual Standard_Real ReversedParameter (const Standard_Real U) const
    {
        boost::ignore_unused(U);
        abort();
        return 0.0;
    }

    virtual Standard_Real TransformedParameter (const Standard_Real U, const gp_Trsf& T) const
    {
        boost::ignore_unused(U, T);
        abort();
        return 0.0;
    }

    virtual Standard_Real ParametricTransformation (const gp_Trsf& T) const
    {
        boost::ignore_unused(T);
        abort();
        return 0.0;
    }

    virtual Standard_Real FirstParameter() const
    {
        return 0.0;
    }

    virtual Standard_Real LastParameter() const
    {
        return 1.0;
    }

    virtual Standard_Boolean IsClosed() const
    {
        return true;
    }

    virtual Standard_Boolean IsPeriodic() const
    {
        return false;
    }

    virtual Standard_Real Period() const
    {
        abort();
        return 0.0;
    }

    virtual GeomAbs_Shape Continuity() const
    {
        return GeomAbs_C1;
    }

    virtual Standard_Boolean IsCN (const Standard_Integer N) const
    {
        if (N > 2 || N < 0)
        {
            return false;
        }
        return true;
    }

    virtual void D0 (const Standard_Real U1, gp_Pnt& P) const
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

    virtual void D1 (const Standard_Real U1, gp_Pnt& P, gp_Vec& V1) const
    {
        boost::ignore_unused(U1, P, V1);
        // double U = m_uvDir ? m_val : U1;
        // double V = m_uvDir ? U1 : m_val;

        // // Calculate derivatives
        // gp_Pnt u0, u1, u2, u3;
        // gp_Vec u0p, u1p, u2p, u3p;
        // m_edges[0]->D1(Map(0, U), u0, u0p);
        // m_edges[1]->D1(Map(1, V), u1, u1p);
        // m_edges[2]->D1(Map(2, 1-U), u2, u2p);
        // m_edges[3]->D1(Map(3, 1-V), u3, u3p);

        // gp_XYZ c0 = u0.XYZ(), c1 = u1.XYZ(), c2 = u2.XYZ(), c3 = u3.XYZ();
        // gp_XYZ c0p = u0p.XYZ(), c1p = u1p.XYZ(), c2p = u2p.XYZ(), c3p = u3p.XYZ();

        // gp_XYZ v0 = m_verts[0].XYZ(), v1 = m_verts[1].XYZ();
        // gp_XYZ v2 = m_verts[2].XYZ(), v3 = m_verts[3].XYZ();

        // gp_XYZ pnt = (1-U) * c3 + U * c1 + (1-V) * c0 + V * c2 - (
        //     (1-U)*(1-V) * v0 + U*(1-V) * v1 + U*V * v2 + (1-U)*V * v3);
        // P = gp_Pnt(pnt);

        // if (m_uvDir)
        // {
        //     gp_XYZ dv = (1-U) * c3p + U * c1p - c0 + c2 - (
        //         -(1-U) * v0 - U * v1 + U * v2 + (1-U) * v3);
        //     V1 = gp_Vec(dv);
        // }
        // else
        // {
        //     gp_XYZ du = -1.0 * c3 + c1 + (1-V)*c0p + V*c2p - (
        //         -(1-V) * v0 + (1-V) * v1 + V * v2 - V * v3);
        //     V1 = gp_Vec(du);
        // }
        abort();
    }

    virtual void D2 (const Standard_Real U1, gp_Pnt& P, gp_Vec& V1, gp_Vec& V2) const
    {
        boost::ignore_unused(U1, P, V1, V2);
        // double U = m_uvDir ? m_val : U1;
        // double V = m_uvDir ? U1 : m_val;

        // // Calculate derivatives
        // gp_Pnt u0, u1, u2, u3;
        // gp_Vec u0p, u1p, u2p, u3p;
        // gp_Vec u0pp, u1pp, u2pp, u3pp;
        // m_edges[0]->D2(Map(0, U), u0, u0p, u0pp);
        // m_edges[1]->D2(Map(1, V), u1, u1p, u1pp);
        // m_edges[2]->D2(Map(2, 1-U), u2, u2p, u2pp);
        // m_edges[3]->D2(Map(3, 1-V), u3, u3p, u3pp);

        // gp_XYZ c0 = u0.XYZ(), c1 = u1.XYZ(), c2 = u2.XYZ(), c3 = u3.XYZ();
        // gp_XYZ c0p = u0p.XYZ(), c1p = u1p.XYZ(), c2p = u2p.XYZ(), c3p = u3p.XYZ();
        // gp_XYZ c0pp = u0pp.XYZ(), c1pp = u1pp.XYZ(), c2pp = u2pp.XYZ(), c3pp = u3pp.XYZ();

        // gp_XYZ v0 = m_verts[0].XYZ(), v1 = m_verts[1].XYZ();
        // gp_XYZ v2 = m_verts[2].XYZ(), v3 = m_verts[3].XYZ();

        // gp_XYZ pnt = (1-U) * c3 + U * c1 + (1-V) * c0 + V * c2 - (
        //     (1-U) * (1-V) * v0 + U * (1-V) * v1 + U * V * v2 + (1-U) * V * v3);
        // P = gp_Pnt(pnt);

        // if (m_uvDir)
        // {
        //     gp_XYZ dv = (1-U) * c3p + U * c1p - c0 + c2 - (
        //         -(1-U) * v0 - U * v1 + U * v2 + (1-U) * v3);
        //     V1 = gp_Vec(dv);

        //     gp_XYZ d2v = (1-U) * c3pp + U * c1pp;
        //     V2 = gp_Vec(d2v);
        // }
        // else
        // {
        //     gp_XYZ du = -1.0 * c3 + c1 + (1-V)*c0p + V*c2p - (
        //         -(1-V) * v0 + (1-V) * v1 + V * v2 - V * v3);
        //     V1 = gp_Vec(du);

        //     gp_XYZ d2u = (1-V) * c0pp + V * c2pp;
        //     V2 = gp_Vec(d2u);
        // }
        abort();
    }

    virtual void D3 (const Standard_Real U, gp_Pnt& P, gp_Vec& V1, gp_Vec& V2, gp_Vec& V3) const
    {
        boost::ignore_unused(U, P, V1, V2, V3);
        abort();
    }

    virtual gp_Vec DN (const Standard_Real U, const Standard_Integer N) const
    {
        boost::ignore_unused(U, N);
        abort();
        return gp_Vec();
    }

    virtual Handle(Geom_Geometry) Copy() const
    {
        Handle(Geom_TransfiniteCurve) tmp = new Geom_TransfiniteCurve(
            m_uvDir, m_val, m_edges, m_fwd, m_clims, m_verts);
        return tmp;
    }

    virtual void Transform (const gp_Trsf& T)
    {
        // Transform our constituent parts.
        abort();
        m_verts[0].Transform(T);
        m_verts[1].Transform(T);
        m_verts[2].Transform(T);
        m_verts[3].Transform(T);
        m_edges[0]->Transform(T);
        m_edges[1]->Transform(T);
        m_edges[2]->Transform(T);
        m_edges[3]->Transform(T);
    }
};

class Geom_TransfiniteSurface : public Geom_BoundedSurface
{
    std::vector<Handle(Geom_Curve)> m_edges;
    std::vector<bool> m_fwd;
    std::vector<std::pair<double, double>> m_clims;
    std::vector<gp_Pnt> m_verts;

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

public:
    Geom_TransfiniteSurface(std::vector<Handle(Geom_Curve)>        edges,
                            std::vector<bool>                      fwd,
                            std::vector<std::pair<double, double>> clims,
                            std::vector<gp_Pnt>                    verts)
        : m_edges(edges), m_fwd(fwd), m_clims(clims), m_verts(verts)
    {
    }

    virtual void UReverse()
    {
        abort();
    }
    Handle(Geom_Surface) UReversed() const
    {
        abort();
        return Handle(Geom_Surface)();
    }
    virtual Standard_Real UReversedParameter(const Standard_Real U) const
    {
        boost::ignore_unused(U);
        abort();
        return 0.0;
    }
    virtual void VReverse()
    {
        abort();
    }
    Handle(Geom_Surface) VReversed() const
    {
        abort();
        return Handle(Geom_Surface)();
    }
    virtual Standard_Real VReversedParameter(const Standard_Real V) const
    {
        boost::ignore_unused(V);
        abort();
        return 0.0;
    }
    virtual void Transform(const gp_Trsf& T)
    {
        abort();
        // Transform our constituent parts.
        m_verts[0].Transform(T);
        m_verts[1].Transform(T);
        m_verts[2].Transform(T);
        m_verts[3].Transform(T);
        m_edges[0]->Transform(T);
        m_edges[1]->Transform(T);
        m_edges[2]->Transform(T);
        m_edges[3]->Transform(T);
    }
    virtual Handle(Geom_Geometry) Copy() const
    {
        Handle(Geom_TransfiniteSurface) tmp = new Geom_TransfiniteSurface(
            m_edges, m_fwd, m_clims, m_verts);
        return tmp;
    }
    virtual void TransformParameters (
        Standard_Real& U, Standard_Real& V, const gp_Trsf& T) const
    {
        boost::ignore_unused(U, V, T);
        abort();
    }
    virtual gp_GTrsf2d ParametricTransformation (const gp_Trsf& T) const
    {
        boost::ignore_unused(T);
        abort();
        return gp_GTrsf2d();
    }
    virtual void Bounds (Standard_Real& U1, Standard_Real& U2,
                         Standard_Real& V1, Standard_Real& V2) const
    {
        U1 = V1 = 0.0;
        U2 = V2 = 1.0;
    }
    virtual Standard_Boolean IsUClosed() const
    {
        return true;
    }
    virtual Standard_Boolean IsVClosed() const
    {
        return true;
    }
    virtual Standard_Boolean IsUPeriodic() const
    {
        return false;
    }
    virtual Standard_Real UPeriod() const
    {
        abort();
        return 0.0;
    }
    virtual Standard_Boolean IsVPeriodic() const
    {
        return false;
    }
    virtual Standard_Real VPeriod() const
    {
        abort();
        return 0.0;
    }
    virtual Handle(Geom_Curve) UIso (const Standard_Real U) const
    {
        Handle(Geom_Curve) c = new Geom_TransfiniteCurve(
            true, U, m_edges, m_fwd, m_clims, m_verts);
        return c;
    }
    virtual Handle(Geom_Curve) VIso (const Standard_Real V) const
    {
        Handle(Geom_Curve) c = new Geom_TransfiniteCurve(
            false, V, m_edges, m_fwd, m_clims, m_verts);
        return c;
    }
    virtual GeomAbs_Shape Continuity() const
    {
        return GeomAbs_C2;
    }
    virtual Standard_Boolean IsCNu (const Standard_Integer N) const
    {
        if (N > 2 || N < 0)
        {
            return false;
        }
        return true;
    }
    virtual Standard_Boolean IsCNv (const Standard_Integer N) const
    {
        if (N > 2 || N < 0)
        {
            return false;
        }
        return true;
    }
    virtual void D0 (const Standard_Real U, const Standard_Real V, gp_Pnt& P) const
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
    virtual void D1 (const Standard_Real U, const Standard_Real V, gp_Pnt& P, gp_Vec& D1U, gp_Vec& D1V) const
    {
        gp_Pnt t0, t1, t2, t3;
        gp_Vec d0, d1, d2, d3;

        m_edges[0]->D1(Map(0, U), t0, d0);
        m_edges[1]->D1(Map(1, V), t1, d1);
        m_edges[2]->D1(Map(2, 1-U), t2, d2);
        m_edges[3]->D1(Map(3, 1-V), t3, d3);

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
        c2p *= m_fwd[2] ? -1.0 : 1.0; // These two edges are flipped
        c3p *= m_fwd[3] ? -1.0 : 1.0; // so negate here.

        // Analytic derivatives of the transfinite surface.
        gp_XYZ du = -1.0 * c3 + c1 + (1-V) * c0p + V * c2p - (
            -(1-V) * v0 + (1-V) * v1 + V * v2 - V * v3);
        gp_XYZ dv = (1-U) * c3p + U * c1p - c0 + c2 - (
            -(1-U) * v0 - U * v1 + U * v2 + (1-U) * v3);

        D1U = gp_Vec(du);
        D1V = gp_Vec(dv);
    }
    virtual void D2 (const Standard_Real U, const Standard_Real V, gp_Pnt& P, gp_Vec& D1U, gp_Vec& D1V, gp_Vec& D2U, gp_Vec& D2V, gp_Vec& D2UV) const
    {
        gp_Pnt t0, t1, t2, t3;
        gp_Vec d0, d1, d2, d3, dd0, dd1, dd2, dd3;

        m_edges[0]->D2(Map(0, U), t0, d0, dd0);
        m_edges[1]->D2(Map(1, V), t1, d1, dd1);
        m_edges[2]->D2(Map(2, 1-U), t2, d2, dd2);
        m_edges[3]->D2(Map(3, 1-V), t3, d3, dd3);

        gp_XYZ c0 = t0.XYZ(), c1 = t1.XYZ(), c2 = t2.XYZ(), c3 = t3.XYZ();
        gp_XYZ c0p = d0.XYZ(), c1p = d1.XYZ(), c2p = d2.XYZ(), c3p = d3.XYZ();
        gp_XYZ c0pp = dd0.XYZ(), c1pp = dd1.XYZ(), c2pp = dd2.XYZ(), c3pp = dd3.XYZ();
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
    virtual void D3 (const Standard_Real U, const Standard_Real V, gp_Pnt& P, gp_Vec& D1U, gp_Vec& D1V, gp_Vec& D2U, gp_Vec& D2V, gp_Vec& D2UV, gp_Vec& D3U, gp_Vec& D3V, gp_Vec& D3UUV, gp_Vec& D3UVV) const
    {
        boost::ignore_unused(U, V, P, D1U, D1V, D2U, D2V, D2UV, D3U, D3V, D3UUV, D3UVV);
        abort();
    }
    virtual gp_Vec DN (const Standard_Real U, const Standard_Real V, const Standard_Integer Nu, const Standard_Integer Nv) const
    {
        boost::ignore_unused(U, V, Nu, Nv);
        abort();
        return gp_Vec();
    }
};
