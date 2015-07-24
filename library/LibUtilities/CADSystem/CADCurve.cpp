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
 * @brief Calculates the parametric coordinate and arclength location
 * defined by \p s.
 *
 * @param s arclength location
 * @return calculated parametric coordinate
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

    ///@todo this really needs improving for accuracy
}

/**
 * @brief Calculates the arclength between the two paremetric points \p ti
 * and \p tf. \p ti must be less than \p tf
 *
 * @param ti first parametric coordinate
 * @param tf second parametric coordinate
 * @return arc length between \p ti and \p tf
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
 * @brief Gets the location x,y,z in array out of the curve at point \p t.
 *
 * @param t parametric coordinate
 * @return array of x,y,z
 */

Array<OneD, NekDouble> CADCurve::P(NekDouble t)
{

    Array<OneD, NekDouble> location(3);
    gp_Pnt loc = occCurve.Value(t);

    location[0] = loc.X();
    location[1] = loc.Y();
    location[2] = loc.Z();

    return location;
}

/**
 * @brief returns the minimum and maximum parametric coords t of the curve.
 *
 * @return array of two entries, min and max parametric coordinate
 */

Array<OneD, NekDouble> CADCurve::Bounds()
{
    Array<OneD, NekDouble> t(2);
    t[0] = occCurve.FirstParameter();
    t[1] = occCurve.LastParameter();

    return t;
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
 *
 * @return array with 6 entries of endpoints x1,y1,z1,x2,y2,z2
 */

Array<OneD, NekDouble> CADCurve::GetMinMax()
{
    Array<OneD, NekDouble> locs(6);
    gp_Pnt start = occCurve.Value(occCurve.FirstParameter());
    gp_Pnt end  = occCurve.Value(occCurve.LastParameter());
    locs[0]=start.X();
    locs[1]=start.Y();
    locs[2]=start.Z();
    locs[3]=end.X();
    locs[4]=end.Y();
    locs[5]=end.Z();

    return locs;

}

}
}
