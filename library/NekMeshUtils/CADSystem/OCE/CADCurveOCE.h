////////////////////////////////////////////////////////////////////////////////
//
//  File: CADCurve.h
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
//  Description: CAD object curve.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_CADSYSTEM_OCE_CADCURVEOCE
#define NEKMESHUTILS_CADSYSTEM_OCE_CADCURVEOCE

#include "../CADCurve.h"
#include "OpenCascade.h"

namespace Nektar
{
namespace NekMeshUtils
{

/**
 * @brief class for CAD curves.
 *
 * This class wraps the OpenCascade BRepAdaptor_Curve class for use with
 * Nektar++.
 */
class CADCurveOCE : public CADCurve
{
public:

    static CADCurveSharedPtr create()
    {
        return MemoryManager<CADCurveOCE>::AllocateSharedPtr();
    }

    static EngineKey key;

    CADCurveOCE(){};

    ~CADCurveOCE(){};

    Array<OneD, NekDouble> Bounds();

    NekDouble Length(NekDouble ti, NekDouble tf);

    Array<OneD, NekDouble> P(NekDouble t);

    Array<OneD, NekDouble> D2(NekDouble t);

    NekDouble tAtArcLength(NekDouble s);

    Array<OneD, NekDouble> GetMinMax();

    void Initialise(int i, TopoDS_Shape in)
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

private:
    /// OpenCascade object of the curve.
    BRepAdaptor_Curve m_occCurve;
    /// OpenCascade edge
    TopoDS_Edge m_occEdge;
};

}
}

#endif
