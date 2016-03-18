////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSurf.cpp
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
//  Description: cad object surface methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <sstream>

#include <NekMeshUtils/CADSystem/CADSurf.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

CADSurf::CADSurf(int i, TopoDS_Shape in, vector<EdgeLoop> ein) : m_edges(ein)
{
    // this bit of code changes the units of the cad from mm opencascade
    // defualt to m
    gp_Trsf transform;
    gp_Pnt ori(0.0, 0.0, 0.0);
    transform.SetScale(ori, 1.0 / 1000.0);
    TopLoc_Location mv(transform);
    m_s = BRep_Tool::Surface(TopoDS::Face(in));
    in.Move(mv);
    m_occSurface    = BRepAdaptor_Surface(TopoDS::Face(in));
    m_correctNormal = true;
    m_hasTwoCurves  = false;
    m_id            = i;
    m_type          = surf;
}

Array<OneD, NekDouble> CADSurf::locuv(Array<OneD, NekDouble> p)
{
    // has to transfer back to mm
    gp_Pnt loc(p[0] * 1000.0, p[1] * 1000.0, p[2] * 1000.0);

    GeomAPI_ProjectPointOnSurf projection(
        loc, m_s, m_occSurface.FirstUParameter(), m_occSurface.LastUParameter(),
        m_occSurface.FirstVParameter(), m_occSurface.LastVParameter(),
        Extrema_ExtAlgo_Tree);

    Array<OneD, NekDouble> uvr(2);
    if (projection.NbPoints() == 0)
    {
        // alternative locuv methods
        ShapeAnalysis_Surface sas(m_s);
        sas.SetDomain(
            m_occSurface.FirstUParameter(), m_occSurface.LastUParameter(),
            m_occSurface.FirstVParameter(), m_occSurface.LastVParameter());

        gp_Pnt2d p2 = sas.ValueOfUV(loc, 1e-7);
        uvr[0]      = p2.X();
        uvr[1]      = p2.Y();

        gp_Pnt p3 = sas.Value(p2);
        ASSERTL0(p3.Distance(loc) < 1e-3, "large locuv distance sas");
    }
    else
    {
        Quantity_Parameter ui;
        Quantity_Parameter vi;

        projection.Parameters(1, ui, vi);

        uvr[0] = ui;
        uvr[1] = vi;

        if (projection.Distance(1) > 1.0)
        {
            stringstream ss;
            cerr << "large locuv distance " << projection.Distance(1) / 1000.0
                 << endl;
        }
    }

    // if the uv returned is slightly off the surface
    //(which ShapeAnalysis_Surface can do sometimes)
    if (uvr[0] < m_occSurface.FirstUParameter() ||
        uvr[0] > m_occSurface.LastUParameter() ||
        uvr[1] < m_occSurface.FirstVParameter() ||
        uvr[1] > m_occSurface.LastVParameter())
    {
        if (uvr[0] < m_occSurface.FirstUParameter() &&
            fabs(m_occSurface.FirstUParameter() - uvr[0]) < 1E-6)
        {
            uvr[0] = m_occSurface.FirstUParameter();
        }
        else if (uvr[0] > m_occSurface.LastUParameter() &&
                 fabs(m_occSurface.LastUParameter() - uvr[0]) < 1E-6)
        {
            uvr[0] = m_occSurface.LastUParameter();
        }
        else if (uvr[1] < m_occSurface.FirstVParameter() &&
                 fabs(m_occSurface.FirstVParameter() - uvr[1]) < 1E-6)
        {
            uvr[1] = m_occSurface.FirstVParameter();
        }
        else if (uvr[1] > m_occSurface.LastVParameter() &&
                 fabs(m_occSurface.LastVParameter() - uvr[1]) < 1E-6)
        {
            uvr[1] = m_occSurface.LastVParameter();
        }
        else
        {
            ASSERTL0(false, "Cannot correct locuv");
        }
    }

    return uvr;
}

NekDouble CADSurf::Curvature(Array<OneD, NekDouble> uv)
{
    Test(uv);

    Array<OneD, NekDouble> n = N(uv);

    // a zero normal occurs at a signularity, CurvaturePoint
    // cannot be sampled here
    if (n[0] == 0 && n[1] == 0 && n[2] == 0)
    {
        return 0.0;
    }

    Array<OneD, NekDouble> r = D2(uv);

    // metric and curvature tensors
    NekDouble E = r[3] * r[3] + r[4] * r[4] + r[5] * r[5];
    NekDouble F = r[3] * r[6] + r[4] * r[7] + r[5] * r[8];
    NekDouble G = r[6] * r[6] + r[7] * r[7] + r[8] * r[8];
    NekDouble e = n[0] * r[9] + n[1] * r[10] + n[2] * r[11];
    NekDouble f = n[0] * r[15] + n[1] * r[16] + n[2] * r[17];
    NekDouble g = n[0] * r[12] + n[1] * r[13] + n[2] * r[14];

    // if det is zero cannot invert matrix, R=0 so must skip
    if (E * G - F * F < 1E-30)
    {
        return 0.0;
    }

    NekDouble K, H;

    K = (e * g - f * f) / (E * G - F * F);
    H = 0.5 * (e * G - 2 * f * F + g * E) / (E * G - F * F);

    NekDouble kv[2];
    kv[0] = abs(H + sqrt(H * H - K));
    kv[1] = abs(H - sqrt(H * H - K));

    return kv[0] > kv[1] ? kv[0] : kv[1];
}

NekDouble CADSurf::DistanceTo(Array<OneD, NekDouble> p)
{
    gp_Pnt loc(p[0] * 1000.0, p[1] * 1000.0, p[2] * 1000.0);

    // alternative locuv methods
    ShapeAnalysis_Surface sas(m_s);
    sas.SetDomain(m_occSurface.FirstUParameter(), m_occSurface.LastUParameter(),
                  m_occSurface.FirstVParameter(),
                  m_occSurface.LastVParameter());

    gp_Pnt2d p2 = sas.ValueOfUV(loc, 1e-7);

    gp_Pnt p3 = sas.Value(p2);

    return p3.Distance(loc);
}

void CADSurf::ProjectTo(Array<OneD, NekDouble> &tp, Array<OneD, NekDouble> &uv)
{
    gp_Pnt loc(tp[0] * 1000.0, tp[1] * 1000.0, tp[2] * 1000.0);

    // alternative locuv methods
    ShapeAnalysis_Surface sas(m_s);
    sas.SetDomain(m_occSurface.FirstUParameter(), m_occSurface.LastUParameter(),
                  m_occSurface.FirstVParameter(),
                  m_occSurface.LastVParameter());

    gp_Pnt2d p2 = sas.ValueOfUV(loc, 1e-7);

    gp_Pnt p3 = sas.Value(p2);

    tp[0] = p3.X() / 1000.0;
    tp[1] = p3.Y() / 1000.0;
    tp[2] = p3.Z() / 1000.0;

    uv[0] = p2.X();
    uv[1] = p2.Y();
}

bool CADSurf::IsPlane()
{
    if (m_occSurface.GetType() == GeomAbs_Plane)
    {
        return true;
    }
    else
    {
        return false;
    }
}

Array<OneD, NekDouble> CADSurf::P(Array<OneD, NekDouble> uv)
{
    Test(uv);

    Array<OneD, NekDouble> location(3);
    gp_Pnt loc;
    loc         = m_occSurface.Value(uv[0], uv[1]);
    location[0] = loc.X();
    location[1] = loc.Y();
    location[2] = loc.Z();
    return location;
}

Array<OneD, NekDouble> CADSurf::N(Array<OneD, NekDouble> uv)
{
    Test(uv);

    Array<OneD, NekDouble> normal(3);
    gp_Pnt Loc;
    gp_Vec D1U, D1V;
    m_occSurface.D1(uv[0], uv[1], Loc, D1U, D1V);
    gp_Vec n = D1U.Crossed(D1V);

    if (!m_correctNormal)
    {
        n.Reverse();
    }

    if (n.X() == 0 && n.Y() == 0 && n.Z() == 0)
    {
        // Return bad normal
        normal[0] = 0.0;
        normal[1] = 0.0;
        normal[2] = 0.0;
    }
    else
    {
        n.Normalize();
        normal[0] = n.X();
        normal[1] = n.Y();
        normal[2] = n.Z();
    }

    return normal;
}

Array<OneD, NekDouble> CADSurf::D1(Array<OneD, NekDouble> uv)
{
    Test(uv);

    Array<OneD, NekDouble> r(9);
    gp_Pnt Loc;
    gp_Vec D1U, D1V;
    m_occSurface.D1(uv[0], uv[1], Loc, D1U, D1V);

    r[0] = Loc.X(); // x
    r[1] = Loc.Y(); // y
    r[2] = Loc.Z(); // z
    r[3] = D1U.X(); // dx/du
    r[4] = D1U.Y(); // dy/du
    r[5] = D1U.Z(); // dz/du
    r[6] = D1V.X(); // dx/dv
    r[7] = D1V.Y(); // dy/dv
    r[8] = D1V.Z(); // dz/dv

    return r;
}

Array<OneD, NekDouble> CADSurf::D2(Array<OneD, NekDouble> uv)
{
    Test(uv);

    Array<OneD, NekDouble> r(18);
    gp_Pnt Loc;
    gp_Vec D1U, D1V, D2U, D2V, D2UV;
    m_occSurface.D2(uv[0], uv[1], Loc, D1U, D1V, D2U, D2V, D2UV);

    r[0]  = Loc.X(); // x
    r[1]  = Loc.Y(); // y
    r[2]  = Loc.Z(); // z
    r[3]  = D1U.X(); // dx/dx
    r[4]  = D1U.Y(); // dy/dy
    r[5]  = D1U.Z(); // dz/dz
    r[6]  = D1V.X(); // dx/dx
    r[7]  = D1V.Y(); // dy/dy
    r[8]  = D1V.Z(); // dz/dz
    r[9]  = D2U.X(); // d2x/du2
    r[10] = D2U.Y(); // d2y/du2
    r[11] = D2U.Z(); // d2z/du2
    r[12] = D2V.X(); // d2x/dv2
    r[13] = D2V.Y(); // d2y/dv2
    r[14] = D2V.Z(); // d2z/dv2
    r[15] = D2UV.X(); // d2x/dudv
    r[16] = D2UV.Y(); // d2y/dudv
    r[17] = D2UV.Z(); // d2z/dudv

    return r;
}

void CADSurf::Test(Array<OneD, NekDouble> uv)
{
    stringstream error;

    error << "Point not within parameter plane: ";

    bool passed = true;

    if (uv[0] < m_occSurface.FirstUParameter())
    {
        if (fabs(uv[0] - m_occSurface.FirstUParameter()) > 1E-8)
        {
            error << "U(" << uv[0] << ") is less than Umin("
                  << m_occSurface.FirstUParameter() << ")";
            passed = false;
        }
    }
    else if (uv[0] > m_occSurface.LastUParameter())
    {
        if (fabs(uv[0] - m_occSurface.LastUParameter()) > 1E-8)
        {
            error << "U(" << uv[0] << ") is greater than Umax("
                  << m_occSurface.LastUParameter() << ")";
            passed = false;
        }
    }
    else if (uv[1] < m_occSurface.FirstVParameter())
    {
        if (fabs(uv[1] - m_occSurface.FirstVParameter()) > 1E-8)
        {
            error << "V(" << uv[1] << ") is less than Vmin("
                  << m_occSurface.FirstVParameter() << ")";
            passed = false;
        }
    }
    else if (uv[1] > m_occSurface.LastVParameter())
    {
        if (fabs(uv[1] - m_occSurface.LastVParameter()) > 1E-8)
        {
            error << "V(" << uv[1] << ") is greater than Vmax("
                  << m_occSurface.LastVParameter() << ")";
            passed = false;
        }
    }

    error << " On Surface: " << GetId();

    ASSERTL0(passed, error.str());
}
}
}
