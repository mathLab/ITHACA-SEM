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

#include <LibUtilities/CADSystem/CADSurf.h>

using namespace std;

namespace Nektar {
namespace LibUtilities {

CADSurf::CADSurf(int i, TopoDS_Shape in,
                 vector<vector<pair<int,int> > > ein) : m_ID(i), m_edges(ein)
{
    // this bit of code changes the units of the cad from mm opencascade
    // defualt to m
    gp_Trsf transform;
    gp_Pnt ori(0.0, 0.0, 0.0);
    transform.SetScale(ori, 1.0 / 1000.0);
    TopLoc_Location mv(transform);
    m_s = BRep_Tool::Surface(TopoDS::Face(in));
    in.Move(mv);
    m_occSurface = BRepAdaptor_Surface(TopoDS::Face(in));
}

/**
 * @brief Performs a reverse look up to find u,v and x,y,z.
 *
 * @param p Array of xyz location
 * @return The parametric location of xyz on this surface
 */
Array<OneD, NekDouble> CADSurf::locuv(Array<OneD, NekDouble> p)
{
    //has to transfer back to mm
    gp_Pnt loc(p[0] * 1000.0, p[1] * 1000.0, p[2] * 1000.0);

    GeomAPI_ProjectPointOnSurf projection(loc, m_s,
                        m_occSurface.FirstUParameter(),
                        m_occSurface.LastUParameter(),
                        m_occSurface.FirstVParameter(),
                        m_occSurface.LastVParameter(),
                        Extrema_ExtAlgo_Tree);

    ASSERTL0(projection.NbPoints() > 0, "locuv failed");

    Quantity_Parameter ui;
    Quantity_Parameter vi;

    projection.Parameters(1,ui,vi);

    /// @todo create a check so that if the calculated uv is out of bounds

    Array<OneD, NekDouble> uvr(2);
    uvr[0] = ui;
    uvr[1] = vi;

    ASSERTL1(projection.Distance(1) < NekConstants::GeomTol,
                "large locuv distance");

    return uvr;
}

/**
 * @brief Get the x,y,z at parametric point u,v.
 *
 * @param uv Array of u and v parametric coords.
 * @return Array of xyz location.
 */
Array<OneD, NekDouble> CADSurf::P(Array<OneD, NekDouble> uv)
{
    /// @todo create bound checking
    Array<OneD, NekDouble> location(3);
    gp_Pnt loc;
    loc = m_occSurface.Value(uv[0], uv[1]);
    location[0] = loc.X();
    location[1] = loc.Y();
    location[2] = loc.Z();
    return location;
}

/**
 * @brief Get the normal vector at parametric point u,v.
 *
 * @param uv Array of u and v parametric coords.
 * @return Array of xyz components of normal vector.
 */
Array<OneD, NekDouble> CADSurf::N(Array<OneD, NekDouble> uv)
{
    Array<OneD, NekDouble> normal(3);
    gp_Pnt Loc;
    gp_Vec D1U, D1V;
    m_occSurface.D1(uv[0], uv[1], Loc, D1U, D1V);
    gp_Vec n = D1U.Crossed(D1V);

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

/**
 * @brief Get the set of first derivatives at parametric point u,v
 *
 * @param uv Array of u and v parametric coords.
 * @return Array of xyz copmonents of first derivatives.
 */

Array<OneD, NekDouble> CADSurf::D1(Array<OneD, NekDouble> uv)
{
    Array<OneD, NekDouble> r(9);
    gp_Pnt Loc;
    gp_Vec D1U, D1V;
    m_occSurface.D1(uv[0], uv[1], Loc, D1U, D1V);

    r[0] = Loc.X();  //x
    r[1] = Loc.Y();  //y
    r[2] = Loc.Z();  //z
    r[3] = D1U.X();  //dx/du
    r[4] = D1U.Y();  //dy/du
    r[5] = D1U.Z();  //dz/du
    r[6] = D1V.X();  //dx/dv
    r[7] = D1V.Y();  //dy/dv
    r[8] = D1V.Z();  //dz/dv

    return r;
}

/**
 * @brief Get the set of second derivatives at parametric point u,v
 *
 * @param uv array of u and v parametric coords
 * @return array of xyz copmonents of second derivatives
 */

Array<OneD, NekDouble> CADSurf::D2(Array<OneD, NekDouble> uv)
{
    Array<OneD, NekDouble> r(18);
    gp_Pnt Loc;
    gp_Vec D1U, D1V, D2U, D2V, D2UV;
    m_occSurface.D2(uv[0], uv[1], Loc, D1U, D1V, D2U, D2V, D2UV);

    r[0] = Loc.X();    //x
    r[1] = Loc.Y();    //y
    r[2] = Loc.Z();    //z
    r[3] = D1U.X();    //dx/dx
    r[4] = D1U.Y();    //dy/dy
    r[5] = D1U.Z();    //dz/dz
    r[6] = D1V.X();    //dx/dx
    r[7] = D1V.Y();    //dy/dy
    r[8] = D1V.Z();    //dz/dz
    r[9] = D2U.X();    //d2x/du2
    r[10] = D2U.Y();   //d2y/du2
    r[11] = D2U.Z();   //d2z/du2
    r[12] = D2V.X();   //d2x/dv2
    r[13] = D2V.Y();   //d2y/dv2
    r[14] = D2V.Z();   //d2z/dv2
    r[15] = D2UV.X();  //d2x/dudv
    r[16] = D2UV.Y();  //d2y/dudv
    r[17] = D2UV.Z();  //d2z/dudv

    return r;
}

}
}
