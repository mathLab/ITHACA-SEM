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

#ifndef NEKMESHUTILS_CADSYSTEM_CFI_CADCURVECFI
#define NEKMESHUTILS_CADSYSTEM_CFI_CADCURVECFI

#include "../CADCurve.h"
#include "CADSystemCFI.h"

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

    static EngineKey key;

    CADCurveCFI(){};

    ~CADCurveCFI(){};

    Array<OneD, NekDouble> Bounds();

    NekDouble Length(NekDouble ti, NekDouble tf);

    Array<OneD, NekDouble> P(NekDouble t);

    Array<OneD, NekDouble> D2(NekDouble t);

    NekDouble tAtArcLength(NekDouble s);

    Array<OneD, NekDouble> GetMinMax();

    void Initialise(int i, cfi::Line* in)
    {
        m_cfiEdge = in;
        m_length = m_cfiEdge->calcLength();

        m_id   = i;
        m_type = curve;
    }

private:
    ///cfi object
    cfi::Line* m_cfiEdge;
};

}
}

#endif
