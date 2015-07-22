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

#include <LibUtilities/CADSystem/CADCurve.h>

using namespace std;
namespace Nektar{
namespace LibUtilities{

    /**
     * @brief Calculates the parametric coordinate t and arclength location s.
     */

    NekDouble CADCurve::tAtArcLength(NekDouble s)
    {

        NekDouble dt = (occCurve.LastParameter() -
                        occCurve.FirstParameter()) / (5000);
        NekDouble t = occCurve.FirstParameter();

        NekDouble len = 0.0;

        while(len <= s)
        {
            gp_Pnt P1,P2;
            gp_Vec drdt1,drdt2;

            occCurve.D1(t,P1,drdt1);
            t += dt;
            occCurve.D1(t,P2,drdt2);

            len += (drdt1.Magnitude() + drdt2.Magnitude()) / 2.0 * dt;
        }

        return t - dt;

        //this really needs improving for accuracy
    }

    /**
     * @brief Calculates the arclength between the two paremetric points ti
     * and tf.
     */

    NekDouble CADCurve::Length(NekDouble ti, NekDouble tf)
    {

        NekDouble len = 0;
        NekDouble dt = (occCurve.LastParameter() -
                        occCurve.FirstParameter()) / (1000 - 1);
        NekDouble t = ti;

        while(t + dt <= tf)
        {
            gp_Pnt P1,P2;
            gp_Vec drdt1,drdt2;

            occCurve.D1(t,P1,drdt1);
            t += dt;
            occCurve.D1(t,P2,drdt2);

            len += (drdt1.Magnitude() + drdt2.Magnitude()) / 2.0 * dt;
        }

        return len;
    }

    /**
     * @brief Gets the location x,y,z in array out of the curve at point t.
     */

    void CADCurve::P(NekDouble t, Array<OneD, NekDouble> &out)
    {

        out = Array<OneD, NekDouble>(3);
        gp_Pnt loc = occCurve.Value(t);

        out[0] = loc.X();
        out[1] = loc.Y();
        out[2] = loc.Z();
    }

    /**
     * @brief returns the minimum and maximum parametric coords t of the curve.
     */

    void CADCurve::Bounds(Array<OneD, NekDouble> &out)
    {
        out = Array<OneD, NekDouble>(2);
        out[0] = occCurve.FirstParameter();
        out[1] = occCurve.LastParameter();
    }

    CADCurve::CADCurve(int i, TopoDS_Shape in) : ID(i)
    {
        gp_Trsf transform;
        gp_Pnt ori(0.0, 0.0, 0.0);
        transform.SetScale(ori, 1.0 / 1000.0);
        TopLoc_Location mv(transform);
        in.Move(mv);
        occCurve = BRepAdaptor_Curve(TopoDS::Edge(in));
    }

    /**
     * @brief gets opencascade point objects for the start and end of the curve.
     */

    //this should really be removed and something better put in its place.
    void CADCurve::GetMinMax(gp_Pnt &start, gp_Pnt &end)
    {
        start = occCurve.Value(occCurve.FirstParameter());
        end  = occCurve.Value(occCurve.LastParameter());
    }
}
}
