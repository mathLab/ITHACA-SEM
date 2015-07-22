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
//  Description: cad object surfcae methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/CADSystem/CADSurf.h>

using namespace std;
namespace Nektar{
namespace LibUtilities{

    CADSurf::CADSurf(int i, TopoDS_Shape in,
                     vector<vector<pair<int,int> > > ein) : ID(i), edges(ein)
    {
        //this bit of code changes the units of the cad from mm opencascade
        //defualt to m
        gp_Trsf transform;
        gp_Pnt ori(0.0, 0.0, 0.0);
        transform.SetScale(ori, 1.0 / 1000.0);
        TopLoc_Location mv(transform);
        s = BRep_Tool::Surface(TopoDS::Face(in));
        in.Move(mv);
        occSurface = BRepAdaptor_Surface(TopoDS::Face(in));

    }

    /**
     * @brief performs a reverse look up to find u,v and x,y,z.
     */

    void CADSurf::locuv(NekDouble &u, NekDouble &v,
                        Array<OneD, NekDouble> p)
    {
        //has to transfer back to mm
        gp_Pnt loc(p[0] * 1000.0, p[1] * 1000.0, p[2] * 1000.0);

        GeomAPI_ProjectPointOnSurf projection(loc, s,
                            occSurface.FirstUParameter(),
                            occSurface.LastUParameter(),
                            occSurface.FirstVParameter(),
                            occSurface.LastVParameter(),
                            Extrema_ExtAlgo_Tree);

        ASSERTL0(projection.NbPoints()>0, "locuv failed");

        Quantity_Parameter ui;
        Quantity_Parameter vi;

        projection.Parameters(1,ui,vi);

        //this is always true, stupid check
        ASSERTL0(abs(ui - occSurface.FirstUParameter()) > -1E-6 &&
                 abs(occSurface.LastUParameter() - ui ) > -1E-6 &&
                 abs(vi - occSurface.FirstVParameter()) > -1E-6 &&
                 abs(occSurface.LastVParameter() - vi ) > -1E-6,
                 "locuv exceeded bounds");

        u = ui;
        v = vi;

        if(!(projection.Distance(1) < 1E-3))
            cout << projection.Distance(1) << endl;

    }

    /**
     * @brief Get the x,y,z at parametric point u,v.
     */

    Array<OneD, NekDouble> CADSurf::P(NekDouble u, NekDouble v)
    {
        ASSERTL0(u-occSurface.FirstUParameter() > -1E-6 &&
                 occSurface.LastUParameter()-u  > -1E-6 &&
                 v-occSurface.FirstVParameter() > -1E-6 &&
                 occSurface.LastVParameter()-v  > -1E-6, "out of bounds");
        Array<OneD, NekDouble> out(3);
        gp_Pnt loc;
        loc = occSurface.Value(u,v);
        out[0] = loc.X();
        out[1] = loc.Y();
        out[2] = loc.Z();
        return out;
    }

    /**
     * @brief Get the normal vector at parametric point u,v.
     */

    Array<OneD, NekDouble> CADSurf::N(NekDouble u, NekDouble v)
    {
        Array<OneD, NekDouble> out(3);
        gp_Pnt Loc;
        gp_Vec D1U, D1V;
        occSurface.D1(u, v, Loc, D1U, D1V);
        gp_Vec n = D1U.Crossed(D1V);
        if(n.X() == 0 && n.Y() == 0 && n.Z() == 0)
        {
            //return bad normal
            out[0] = 0.0;
            out[1] = 0.0;
            out[2] = 0.0;
        }
        else
        {
            n.Normalize();
            out[0] = n.X();
            out[1] = n.Y();
            out[2] = n.Z();
        }

        return out;
    }

    /**
     * @brief Get the set of first derivatives at parametric point u,v..
     */

    Array<OneD, NekDouble> CADSurf::D1(NekDouble u, NekDouble v)
    {
        Array<OneD, NekDouble> out(9);
        gp_Pnt Loc;
        gp_Vec D1U, D1V;
        occSurface.D1(u, v, Loc, D1U, D1V);

        out[0] = Loc.X();
        out[1] = Loc.Y();
        out[2] = Loc.Z();
        out[3] = D1U.X();
        out[4] = D1U.Y();
        out[5] = D1U.Z();
        out[6] = D1V.X();
        out[7] = D1V.Y();
        out[8] = D1V.Z();

        return out;
    }

    /**
     * @brief Get the set of second derivatives at parametric point u,v..
     */

    Array<OneD, NekDouble> CADSurf::D2(NekDouble u, NekDouble v)
    {
        Array<OneD, NekDouble> out(18);
        gp_Pnt Loc;
        gp_Vec D1U, D1V, D2U, D2V, D2UV;
        occSurface.D2(u, v, Loc, D1U, D1V, D2U, D2V, D2UV);

        out[0] = Loc.X();
        out[1] = Loc.Y();
        out[2] = Loc.Z();
        out[3] = D1U.X();
        out[4] = D1U.Y();
        out[5] = D1U.Z();
        out[6] = D1V.X();
        out[7] = D1V.Y();
        out[8] = D1V.Z();
        out[9] = D2U.X();
        out[10] = D2U.Y();
        out[11] = D2U.Z();
        out[12] = D2V.X();
        out[13] = D2V.Y();
        out[14] = D2V.Z();
        out[15] = D2UV.X();
        out[16] = D2UV.Y();
        out[17] = D2UV.Z();

        return out;
    }
}
}
