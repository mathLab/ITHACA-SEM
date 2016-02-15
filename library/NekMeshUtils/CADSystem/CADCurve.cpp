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
//  License for the specific language governing rights and limitations under
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

#include <NekMeshUtils/CADSystem/CADCurve.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

CADCurve::CADCurve(int i, TopoDS_Shape in)
{
    gp_Trsf transform;
    gp_Pnt ori(0.0, 0.0, 0.0);
    transform.SetScale(ori, 1.0 / 1000.0);
    TopLoc_Location mv(transform);
    in.Move(mv);

    m_occEdge  = TopoDS::Edge(in);
    m_occCurve = BRepAdaptor_Curve(m_occEdge);

    GProp_GProps System;
    BRepGProp::LinearProperties(m_occEdge, System);
    m_length = System.Mass();

    m_id   = i;
    m_type = curve;
}

NekDouble CADCurve::tAtArcLength(NekDouble s)
{
    NekDouble dt =
        (m_occCurve.LastParameter() - m_occCurve.FirstParameter()) / (5000);
    NekDouble t = m_occCurve.FirstParameter();

    NekDouble len = 0.0;

    while (len <= s)
    {
        gp_Pnt P1, P2;
        gp_Vec drdt1, drdt2;

        m_occCurve.D1(t, P1, drdt1);
        t += dt;
        m_occCurve.D1(t, P2, drdt2);

        len += (drdt1.Magnitude() + drdt2.Magnitude()) / 2.0 * dt;
    }

    return t - dt;
}

NekDouble CADCurve::Length(NekDouble ti, NekDouble tf)
{
    Array<OneD, NekDouble> b = Bounds();
    Handle(Geom_Curve) m_c = BRep_Tool::Curve(m_occEdge, b[0], b[1]);
    Handle(Geom_Curve) NewCurve = new Geom_TrimmedCurve(m_c, ti, tf);
    TopoDS_Edge NewEdge = BRepBuilderAPI_MakeEdge(NewCurve);
    GProp_GProps System;
    BRepGProp::LinearProperties(NewEdge, System);
    return System.Mass() / 1000.0;
}

Array<OneD, NekDouble> CADCurve::P(NekDouble t)
{
    Array<OneD, NekDouble> location(3);
    gp_Pnt loc = m_occCurve.Value(t);

    location[0] = loc.X();
    location[1] = loc.Y();
    location[2] = loc.Z();

    return location;
}

Array<OneD, NekDouble> CADCurve::D2(NekDouble t)
{
    Array<OneD, NekDouble> out(9);
    gp_Pnt loc;
    gp_Vec d1, d2;
    m_occCurve.D2(t, loc, d1, d2);

    out[0] = loc.X();
    out[1] = loc.Y();
    out[2] = loc.Z();
    out[3] = d1.X();
    out[4] = d1.Y();
    out[5] = d1.Z();
    out[6] = d2.X();
    out[7] = d2.Y();
    out[8] = d2.Z();

    return out;
}

Array<OneD, NekDouble> CADCurve::Bounds()
{
    Array<OneD, NekDouble> t(2);
    t[0] = m_occCurve.FirstParameter();
    t[1] = m_occCurve.LastParameter();

    return t;
}

Array<OneD, NekDouble> CADCurve::GetMinMax()
{
    Array<OneD, NekDouble> locs(6);

    gp_Pnt start =
        BRep_Tool::Pnt(TopExp::FirstVertex(m_occEdge, Standard_True));
    gp_Pnt end = BRep_Tool::Pnt(TopExp::LastVertex(m_occEdge, Standard_True));

    locs[0] = start.X();
    locs[1] = start.Y();
    locs[2] = start.Z();
    locs[3] = end.X();
    locs[4] = end.Y();
    locs[5] = end.Z();

    return locs;
}
}
}
