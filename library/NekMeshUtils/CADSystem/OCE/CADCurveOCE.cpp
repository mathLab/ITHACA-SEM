////////////////////////////////////////////////////////////////////////////////
//
//  File: CADCurve.cpp
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
//  Description: CAD object curve methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/CADSystem/OCE/CADCurveOCE.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

std::string CADCurveOCE::key = GetCADCurveFactory().RegisterCreatorFunction(
    "oce", CADCurveOCE::create, "CADCurveOCE");

void CADCurveOCE::Initialise(int i, TopoDS_Shape in)
{
    m_occEdge = TopoDS::Edge(in);

    GProp_GProps System;
    BRepGProp::LinearProperties(m_occEdge, System);
    m_length = System.Mass() / 1000.0;

    m_b = Array<OneD, NekDouble>(2);
    m_c = BRep_Tool::Curve(TopoDS::Edge(in), m_b[0], m_b[1]);

    m_id = i;
}

NekDouble CADCurveOCE::tAtArcLength(NekDouble s)
{
    GeomAdaptor_Curve c(m_c);
    GCPnts_AbscissaPoint ap(c, s * 1000.0, m_b[0]);
    return ap.Parameter();
}

NekDouble CADCurveOCE::Length(NekDouble ti, NekDouble tf)
{
    Handle(Geom_Curve) NewCurve = new Geom_TrimmedCurve(m_c, ti, tf);
    TopoDS_Edge NewEdge         = BRepBuilderAPI_MakeEdge(NewCurve);
    GProp_GProps System;
    BRepGProp::LinearProperties(NewEdge, System);
    return System.Mass() / 1000.0;
}

NekDouble CADCurveOCE::loct(Array<OneD, NekDouble> xyz, NekDouble &t)
{
    t = 0.0;

    gp_Pnt loc(xyz[0] * 1000.0, xyz[1] * 1000.0, xyz[2] * 1000.0);

    ShapeAnalysis_Curve sac;
    gp_Pnt p;
    sac.Project(m_c, loc, Precision::Confusion(), p, t);

    return p.Distance(loc) / 1000.0;
}

Array<OneD, NekDouble> CADCurveOCE::P(NekDouble t)
{
    Array<OneD, NekDouble> location(3);
    gp_Pnt loc = m_c->Value(t);

    location[0] = loc.X() / 1000.0;
    location[1] = loc.Y() / 1000.0;
    location[2] = loc.Z() / 1000.0;

    return location;
}

void CADCurveOCE::P(NekDouble t, NekDouble &x, NekDouble &y, NekDouble &z)
{
    gp_Pnt loc = m_c->Value(t);

    x = loc.X() / 1000.0;
    y = loc.Y() / 1000.0;
    z = loc.Z() / 1000.0;
}

Array<OneD, NekDouble> CADCurveOCE::D2(NekDouble t)
{
    Array<OneD, NekDouble> out(9);
    gp_Pnt loc;
    gp_Vec d1, d2;
    m_c->D2(t, loc, d1, d2);

    out[0] = loc.X() / 1000.0;
    out[1] = loc.Y() / 1000.0;
    out[2] = loc.Z() / 1000.0;
    out[3] = d1.X() / 1000.0;
    out[4] = d1.Y() / 1000.0;
    out[5] = d1.Z() / 1000.0;
    out[6] = d2.X() / 1000.0;
    out[7] = d2.Y() / 1000.0;
    out[8] = d2.Z() / 1000.0;

    return out;
}

Array<OneD, NekDouble> CADCurveOCE::N(NekDouble t)
{
    GeomLProp_CLProps d(m_c, 2, Precision::Confusion());
    d.SetParameter(t + 1e-8);

    gp_Vec d2 = d.D2();
    if (d2.Magnitude() < 1e-8)
    {
        // no normal, stright line
        return Array<OneD, NekDouble>(3, 0.0);
    }

    gp_Dir n;
    d.Normal(n);

    Array<OneD, NekDouble> N(3);
    N[0] = n.X();
    N[1] = n.Y();
    N[2] = n.Z();

    return N;
}

NekDouble CADCurveOCE::Curvature(NekDouble t)
{
    GeomLProp_CLProps d(m_c, 2, Precision::Confusion());
    d.SetParameter(t);

    return d.Curvature() * 1000.0;
}

Array<OneD, NekDouble> CADCurveOCE::GetBounds()
{
    return m_b;
}

void CADCurveOCE::GetBounds(NekDouble &tmin, NekDouble &tmax)
{
    tmin = m_b[0];
    tmax = m_b[1];
}

Array<OneD, NekDouble> CADCurveOCE::GetMinMax()
{
    Array<OneD, NekDouble> locs(6);

    gp_Pnt start =
        BRep_Tool::Pnt(TopExp::FirstVertex(m_occEdge, Standard_True));
    gp_Pnt end = BRep_Tool::Pnt(TopExp::LastVertex(m_occEdge, Standard_True));

    locs[0] = start.X() / 1000.0;
    locs[1] = start.Y() / 1000.0;
    locs[2] = start.Z() / 1000.0;
    locs[3] = end.X() / 1000.0;
    locs[4] = end.Y() / 1000.0;
    locs[5] = end.Z() / 1000.0;

    return locs;
}
}
}
