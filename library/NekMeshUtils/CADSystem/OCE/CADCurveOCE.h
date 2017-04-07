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

#include <NekMeshUtils/CADSystem/CADCurve.h>
#include <NekMeshUtils/CADSystem/OCE/OpenCascade.h>

namespace Nektar
{
namespace NekMeshUtils
{

class CADCurveOCE : public CADCurve
{
public:

    static CADCurveSharedPtr create()
    {
        return MemoryManager<CADCurveOCE>::AllocateSharedPtr();
    }

    static std::string key;

    CADCurveOCE()
    {
    }

    ~CADCurveOCE()
    {
    }

    virtual Array<OneD, NekDouble> GetBounds();
    virtual NekDouble Length(NekDouble ti, NekDouble tf);
    virtual Array<OneD, NekDouble> P(NekDouble t);
    virtual Array<OneD, NekDouble> D2(NekDouble t);
    virtual NekDouble tAtArcLength(NekDouble s);
    virtual Array<OneD, NekDouble> GetMinMax();
    virtual NekDouble loct(Array<OneD, NekDouble> xyz);
    virtual NekDouble Curvature(NekDouble t);
    virtual Array<OneD, NekDouble> NormalWRT(NekDouble t, int surf);
    virtual Array<OneD, NekDouble> N(NekDouble t);

    void Initialise(int i, TopoDS_Shape in)
    {
        gp_Trsf transform;
        gp_Pnt ori(0.0, 0.0, 0.0);
        transform.SetScale(ori, 1.0 / 1000.0);
        TopLoc_Location mv(transform);
        TopoDS_Shape cp = in;
        in.Move(mv);

        m_occEdge  = TopoDS::Edge(in);
        m_occCurve = BRepAdaptor_Curve(m_occEdge);

        GProp_GProps System;
        BRepGProp::LinearProperties(m_occEdge, System);
        m_length = System.Mass();

        Array<OneD, NekDouble> b = GetBounds();
        m_c = BRep_Tool::Curve(TopoDS::Edge(cp), b[0], b[1]);

        m_id   = i;
    }

private:
    /// OpenCascade object of the curve.
    BRepAdaptor_Curve m_occCurve;
    /// OpenCascade edge
    TopoDS_Edge m_occEdge;
    /// Alternate object used for reverse lookups
    Handle(Geom_Curve) m_c;
};

}
}

#endif
