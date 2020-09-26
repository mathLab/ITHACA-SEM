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
#include <NekMeshUtils/CADSystem/OCE/TransfiniteSurface.h>
#include <GeomAPI_ProjectPointOnSurf.hxx>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

std::string CADSurfOCE::key = GetCADSurfFactory().RegisterCreatorFunction(
    "oce", CADSurfOCE::create, "CADSurfOCE");

void CADSurfOCE::Initialise(int i, TopoDS_Shape in)
{
    m_s = BRep_Tool::Surface(TopoDS::Face(in));

    // Test to see if this surface is a transfinite surface, since some OCC
    // functions will not work on our custom type.
    Handle(Geom_TransfiniteSurface) tf = Handle(Geom_TransfiniteSurface)::
        DownCast(m_s);
    m_isTransfiniteSurf = !tf.IsNull();

    if (in.Orientation() == TopAbs_REVERSED)
    {
        m_orientation = CADOrientation::eBackwards;
    }

    m_id = i;

    m_bounds = Array<OneD, NekDouble>(4);
    BRepTools::UVBounds(TopoDS::Face(in), m_bounds[0], m_bounds[1], m_bounds[2],
                        m_bounds[3]);
    m_sas = new ShapeAnalysis_Surface(m_s);
    m_sas->SetDomain(m_bounds[0], m_bounds[1], m_bounds[2], m_bounds[3]);

    m_shape = in;

    m_2Dclass = new BRepTopAdaptor_FClass2d(TopoDS::Face(m_shape), 1e-4);
}

Array<OneD, NekDouble> CADSurfOCE::GetBounds()
{
    return m_bounds;
}

void CADSurfOCE::GetBounds(NekDouble &umin, NekDouble &umax, NekDouble &vmin,
                           NekDouble &vmax)
{
    umin = m_bounds[0];
    umax = m_bounds[1];
    vmin = m_bounds[2];
    vmax = m_bounds[3];
}

bool CADSurfOCE::IsPlanar()
{
    if (m_sas->Adaptor3d()->GetType() == GeomAbs_Plane)
    {
        return true;
    }

    return false;
}

Array<OneD, NekDouble> CADSurfOCE::BoundingBox()
{
    BRepMesh_IncrementalMesh brmsh(m_shape, 0.005);

    Bnd_Box B;
    BRepBndLib::Add(m_shape, B);
    NekDouble e = sqrt(B.SquareExtent()) * 0.01;
    e           = min(e, 5e-3);
    B.Enlarge(e);
    Array<OneD, NekDouble> ret(6);
    B.Get(ret[0], ret[1], ret[2], ret[3], ret[4], ret[5]);
    ret[0] /= 1000.0;
    ret[1] /= 1000.0;
    ret[2] /= 1000.0;
    ret[3] /= 1000.0;
    ret[4] /= 1000.0;
    ret[5] /= 1000.0;
    return ret;
}

Array<OneD, NekDouble> CADSurfOCE::locuv(Array<OneD, NekDouble> p,
                                         NekDouble &dist)
{
    gp_Pnt loc(p[0] * 1000.0, p[1] * 1000.0, p[2] * 1000.0);
    Array<OneD, NekDouble> uv(2);

    if (!m_isTransfiniteSurf)
    {
        gp_Pnt2d p2 = m_sas->ValueOfUV(loc, Precision::Confusion());

        TopAbs_State s = m_2Dclass->Perform(p2);

        if (s == TopAbs_OUT)
        {
            BRepBuilderAPI_MakeVertex v(loc);
            BRepExtrema_DistShapeShape dss(
                BRepTools::OuterWire(TopoDS::Face(m_shape)), v.Shape());
            dss.Perform();
            gp_Pnt np = dss.PointOnShape1(1);
            p2        = m_sas->ValueOfUV(np, Precision::Confusion());
        }

        uv[0] = p2.X();
        uv[1] = p2.Y();

        gp_Pnt p3 = m_sas->Value(p2);

        dist = p3.Distance(loc) / 1000.0;
    }
    else
    {
        Array<OneD, NekDouble> out(3);
        GeomAPI_ProjectPointOnSurf proj(loc, m_s, Precision::Confusion());
        proj.Perform(loc);
        ASSERTL1(proj.NbPoints() > 0, "Unable to find a projection!");
        proj.LowerDistanceParameters(uv[0], uv[1]);
        dist = proj.LowerDistance();
    }

    return uv;
}

NekDouble CADSurfOCE::Curvature(Array<OneD, NekDouble> uv)
{
#if defined(NEKTAR_DEBUG)
    Test(uv);
#endif

    GeomLProp_SLProps d(m_s, 2, Precision::Confusion());
    d.SetParameters(uv[0], uv[1]);

    if (!d.IsCurvatureDefined())
    {
        return -1.0;
    }

    return d.MaxCurvature() * 1000.0;
}

Array<OneD, NekDouble> CADSurfOCE::P(Array<OneD, NekDouble> uv)
{
#if defined(NEKTAR_DEBUG)
    Test(uv);
#endif

    gp_Pnt loc = m_s->Value(uv[0], uv[1]);
    Array<OneD, NekDouble> location(3);
    location[0] = loc.X() / 1000.0;
    location[1] = loc.Y() / 1000.0;
    location[2] = loc.Z() / 1000.0;
    return location;
}

void CADSurfOCE::P(Array<OneD, NekDouble> uv, NekDouble &x, NekDouble &y,
                   NekDouble &z)
{
#if defined(NEKTAR_DEBUG)
    Test(uv);
#endif

    gp_Pnt loc = m_s->Value(uv[0], uv[1]);
    x          = loc.X() / 1000.0;
    y          = loc.Y() / 1000.0;
    z          = loc.Z() / 1000.0;
}

Array<OneD, NekDouble> CADSurfOCE::N(Array<OneD, NekDouble> uv)
{
#if defined(NEKTAR_DEBUG)
    Test(uv);
#endif

    GeomLProp_SLProps d(m_s, 2, Precision::Confusion());
    d.SetParameters(uv[0], uv[1]);

    Array<OneD, NekDouble> normal(3);

    if (!d.IsNormalDefined())
    {
        normal = Array<OneD, NekDouble>(3, 0.0);
        return normal;
    }

    gp_Dir n = d.Normal();

    if (m_orientation == CADOrientation::eBackwards)
    {
        n.Reverse();
    }

    normal[0] = n.X();
    normal[1] = n.Y();
    normal[2] = n.Z();

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
    m_s->D1(uv[0], uv[1], Loc, D1U, D1V);

    r[0] = Loc.X() / 1000.0; // x
    r[1] = Loc.Y() / 1000.0; // y
    r[2] = Loc.Z() / 1000.0; // z
    r[3] = D1U.X() / 1000.0; // dx/du
    r[4] = D1U.Y() / 1000.0; // dy/du
    r[5] = D1U.Z() / 1000.0; // dz/du
    r[6] = D1V.X() / 1000.0; // dx/dv
    r[7] = D1V.Y() / 1000.0; // dy/dv
    r[8] = D1V.Z() / 1000.0; // dz/dv

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
    m_s->D2(uv[0], uv[1], Loc, D1U, D1V, D2U, D2V, D2UV);

    r[0]  = Loc.X() / 1000.0;  // x
    r[1]  = Loc.Y() / 1000.0;  // y
    r[2]  = Loc.Z() / 1000.0;  // z
    r[3]  = D1U.X() / 1000.0;  // dx/dx
    r[4]  = D1U.Y() / 1000.0;  // dy/dy
    r[5]  = D1U.Z() / 1000.0;  // dz/dz
    r[6]  = D1V.X() / 1000.0;  // dx/dx
    r[7]  = D1V.Y() / 1000.0;  // dy/dy
    r[8]  = D1V.Z() / 1000.0;  // dz/dz
    r[9]  = D2U.X() / 1000.0;  // d2x/du2
    r[10] = D2U.Y() / 1000.0;  // d2y/du2
    r[11] = D2U.Z() / 1000.0;  // d2z/du2
    r[12] = D2V.X() / 1000.0;  // d2x/dv2
    r[13] = D2V.Y() / 1000.0;  // d2y/dv2
    r[14] = D2V.Z() / 1000.0;  // d2z/dv2
    r[15] = D2UV.X() / 1000.0; // d2x/dudv
    r[16] = D2UV.Y() / 1000.0; // d2y/dudv
    r[17] = D2UV.Z() / 1000.0; // d2z/dudv

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
    WARNINGL1(passed, "Warning: " + error.str());
    (void)passed; // suppress warning
}


} // namespace NekMeshUtils
} // namespace Nektar
