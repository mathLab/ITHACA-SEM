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

#include <NekMeshUtils/CADSystem/OCE/CADSurfOCE.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

std::string CADSurfOCE::key = GetCADSurfFactory().RegisterCreatorFunction(
    "oce", CADSurfOCE::create, "CADSurfOCE");

void CADSurfOCE::Initialise(int i, TopoDS_Shape in)
{
    // this bit of code changes the units of the cad from mm opencascade
    // defualt to m

    m_s = BRep_Tool::Surface(TopoDS::Face(in));

    if (in.Orientation() == 1)
    {
        m_orientation = CADOrientation::eBackwards;
    }

    gp_Trsf transform;
    gp_Pnt ori(0.0, 0.0, 0.0);
    transform.SetScale(ori, 1.0 / 1000.0);
    TopLoc_Location mv(transform);

    in.Move(mv);
    m_occSurface = BRepAdaptor_Surface(TopoDS::Face(in));
    m_id         = i;

    m_bounds = Array<OneD, NekDouble>(4);
    BRepTools::UVBounds(TopoDS::Face(in), m_bounds[0], m_bounds[1], m_bounds[2],
                        m_bounds[3]);
    m_sas = new ShapeAnalysis_Surface(m_s);
    m_sas->SetDomain(m_bounds[0], m_bounds[1], m_bounds[2], m_bounds[3]);
}

Array<OneD, NekDouble> CADSurfOCE::GetBounds()
{
    return m_bounds;
}

Array<OneD, NekDouble> CADSurfOCE::locuv(Array<OneD, NekDouble> p)
{
    // has to transfer back to mm
    gp_Pnt loc(p[0] * 1000.0, p[1] * 1000.0, p[2] * 1000.0);

    Array<OneD, NekDouble> uvr(2);

    gp_Pnt2d p2 = m_sas->ValueOfUV(loc, 1e-3);
    uvr[0]      = p2.X();
    uvr[1]      = p2.Y();

    gp_Pnt p3 = m_sas->Value(p2);
    if (p3.Distance(loc) > 1.0)
    {
        cout << "large locuv distance " << p3.Distance(loc) << " " << m_id
             << endl;
    }

    // if the uv returned is slightly off the surface
    //(which ShapeAnalysis_Surface can do sometimes)
    if (uvr[0] < m_bounds[0] || uvr[0] > m_bounds[1] || uvr[1] < m_bounds[2] ||
        uvr[1] > m_bounds[3])
    {
        if (uvr[0] < m_bounds[0])
        {
            uvr[0] = m_bounds[0];
        }
        else if (uvr[0] > m_bounds[1])
        {
            uvr[0] = m_bounds[1];
        }
        else if (uvr[1] < m_bounds[2])
        {
            uvr[1] = m_bounds[2];
        }
        else if (uvr[1] > m_bounds[3])
        {
            uvr[1] = m_bounds[3];
        }
        else
        {
            ASSERTL0(false, "Cannot correct locuv");
        }
    }

    return uvr;
}

NekDouble CADSurfOCE::Curvature(Array<OneD, NekDouble> uv)
{
#if defined(NEKTAR_DEBUG)
    Test(uv);
#endif

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

NekDouble CADSurfOCE::DistanceTo(Array<OneD, NekDouble> p)
{
    gp_Pnt loc(p[0] * 1000.0, p[1] * 1000.0, p[2] * 1000.0);

    // alternative locuv methods
    ShapeAnalysis_Surface sas(m_s);
    sas.SetDomain(m_bounds[0], m_bounds[1], m_bounds[2], m_bounds[3]);

    gp_Pnt2d p2 = sas.ValueOfUV(loc, 1e-7);

    gp_Pnt p3 = sas.Value(p2);

    return p3.Distance(loc);
}

void CADSurfOCE::ProjectTo(Array<OneD, NekDouble> &tp,
                           Array<OneD, NekDouble> &uv)
{
    gp_Pnt loc(tp[0] * 1000.0, tp[1] * 1000.0, tp[2] * 1000.0);

    // alternative locuv methods
    ShapeAnalysis_Surface sas(m_s);
    sas.SetDomain(m_bounds[0], m_bounds[1], m_bounds[2], m_bounds[3]);

    gp_Pnt2d p2 = sas.ValueOfUV(loc, 1e-7);

    gp_Pnt p3 = sas.Value(p2);

    tp[0] = p3.X() / 1000.0;
    tp[1] = p3.Y() / 1000.0;
    tp[2] = p3.Z() / 1000.0;

    uv[0] = p2.X();
    uv[1] = p2.Y();
}

Array<OneD, NekDouble> CADSurfOCE::P(Array<OneD, NekDouble> uv)
{
#if defined(NEKTAR_DEBUG)
    Test(uv);
#endif

    Array<OneD, NekDouble> location(3);
    gp_Pnt loc;
    loc         = m_occSurface.Value(uv[0], uv[1]);
    location[0] = loc.X();
    location[1] = loc.Y();
    location[2] = loc.Z();
    return location;
}

Array<OneD, NekDouble> CADSurfOCE::N(Array<OneD, NekDouble> uv)
{
#if defined(NEKTAR_DEBUG)
    Test(uv);
#endif

    BRepLProp_SLProps slp(m_occSurface, 2, 1e-6);
    slp.SetParameters(uv[0], uv[1]);

    if (!slp.IsNormalDefined())
    {
        return Array<OneD, NekDouble>(3, 0.0);
    }

    gp_Dir d = slp.Normal();

    Array<OneD, NekDouble> normal(3);

    if (m_orientation == CADOrientation::eBackwards)
    {
        d.Reverse();
    }

    normal[0] = d.X();
    normal[1] = d.Y();
    normal[2] = d.Z();

    return normal;
}

Array<OneD, NekDouble> CADSurfOCE::D1(Array<OneD, NekDouble> uv)
{
#if defined(NEKTAR_DEBUG)
    Test(uv);
#endif

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

Array<OneD, NekDouble> CADSurfOCE::D2(Array<OneD, NekDouble> uv)
{
#if defined(NEKTAR_DEBUG)
    Test(uv);
#endif

    Array<OneD, NekDouble> r(18);
    gp_Pnt Loc;
    gp_Vec D1U, D1V, D2U, D2V, D2UV;
    m_occSurface.D2(uv[0], uv[1], Loc, D1U, D1V, D2U, D2V, D2UV);

    r[0]  = Loc.X();  // x
    r[1]  = Loc.Y();  // y
    r[2]  = Loc.Z();  // z
    r[3]  = D1U.X();  // dx/dx
    r[4]  = D1U.Y();  // dy/dy
    r[5]  = D1U.Z();  // dz/dz
    r[6]  = D1V.X();  // dx/dx
    r[7]  = D1V.Y();  // dy/dy
    r[8]  = D1V.Z();  // dz/dz
    r[9]  = D2U.X();  // d2x/du2
    r[10] = D2U.Y();  // d2y/du2
    r[11] = D2U.Z();  // d2z/du2
    r[12] = D2V.X();  // d2x/dv2
    r[13] = D2V.Y();  // d2y/dv2
    r[14] = D2V.Z();  // d2z/dv2
    r[15] = D2UV.X(); // d2x/dudv
    r[16] = D2UV.Y(); // d2y/dudv
    r[17] = D2UV.Z(); // d2z/dudv

    return r;
}

void CADSurfOCE::Test(Array<OneD, NekDouble> uv)
{
    stringstream error;

    error << "Point not within parameter plane: ";

    bool passed = true;

    if (uv[0] < m_bounds[0])
    {
        if (fabs(uv[0] - m_bounds[0]) > 1E-6)
        {
            error << "U(" << uv[0] << ") is less than Umin(" << m_bounds[0]
                  << ")";
            passed = false;
        }
    }
    else if (uv[0] > m_bounds[1])
    {
        if (fabs(uv[0] - m_bounds[1]) > 1E-6)
        {
            error << "U(" << uv[0] << ") is greater than Umax(" << m_bounds[1]
                  << ")";
            passed = false;
        }
    }
    else if (uv[1] < m_bounds[2])
    {
        if (fabs(uv[1] - m_bounds[2]) > 1E-6)
        {
            error << "V(" << uv[1] << ") is less than Vmin(" << m_bounds[2]
                  << ")";
            passed = false;
        }
    }
    else if (uv[1] > m_bounds[3])
    {
        if (fabs(uv[1] - m_bounds[3]) > 1E-6)
        {
            error << "V(" << uv[1] << ") is greater than Vmax(" << m_bounds[3]
                  << ")";
            passed = false;
        }
    }

    error << " On Surface: " << GetId();
    ASSERTL1(passed, "Warning: " + error.str());
}
}
}
