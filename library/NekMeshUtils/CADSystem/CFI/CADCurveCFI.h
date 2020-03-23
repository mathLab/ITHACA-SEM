////////////////////////////////////////////////////////////////////////////////
//
//  File: CADCurveCFI.h
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
//  Description: CAD object curve.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_CADSYSTEM_CFI_CADCURVECFI
#define NEKMESHUTILS_CADSYSTEM_CFI_CADCURVECFI

#include "../CADCurve.h"

#ifndef NEK_CADFIXAPI_HXX
#define NEK_CADFIXAPI_HXX
#include "cadfixapi.hxx"
#endif

namespace Nektar
{
namespace NekMeshUtils
{
class CADCurveCFI : public CADCurve
{
public:
    static CADCurveSharedPtr create()
    {
        return MemoryManager<CADCurveCFI>::AllocateSharedPtr();
    }

    static std::string key;

    CADCurveCFI()
    {
    }

    ~CADCurveCFI()
    {
    }

    virtual Array<OneD, NekDouble> GetBounds();
    virtual void GetBounds(NekDouble &tmin, NekDouble &tmax);
    virtual NekDouble Length(NekDouble ti, NekDouble tf);
    virtual Array<OneD, NekDouble> P(NekDouble t);
    virtual void P(NekDouble t, NekDouble &x, NekDouble &y, NekDouble &z);
    virtual Array<OneD, NekDouble> D2(NekDouble t);
    virtual NekDouble Curvature(NekDouble t)
    {
        boost::ignore_unused(t);
        ASSERTL0(false, "Function: Curvature not implemented in CFI engine");
        return 0;
    }
    virtual Array<OneD, NekDouble> N(NekDouble t)
    {
        boost::ignore_unused(t);
        ASSERTL0(false, "Function: N not implemented in CFI engine");
        return Array<OneD, NekDouble>();
    }
    virtual NekDouble tAtArcLength(NekDouble s);
    virtual Array<OneD, NekDouble> GetMinMax();
    virtual NekDouble loct(Array<OneD, NekDouble> xyz, NekDouble &t);

    void Initialise(int i, cfi::Line *in, NekDouble s);

    cfi::Line *GetCfiPointer()
    {
        return m_cfiEdge;
    }

private:
    /// cfi object
    cfi::Line *m_cfiEdge;
    NekDouble m_scal;
};
}
}

#endif
