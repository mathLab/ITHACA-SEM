////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSurfCFI.h
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
//  Description: CAD object surface.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_CADSYSTEM_CFI_CADSURFCFI
#define NEKMESHUTILS_CADSYSTEM_CFI_CADSURFCFI

#include "../CADSurf.h"

#ifndef NEK_CADFIXAPI_HXX
#define NEK_CADFIXAPI_HXX
#include "cadfixapi.hxx"
#endif
namespace Nektar
{
namespace NekMeshUtils
{

class CADSurfCFI : public CADSurf
{
public:
    static CADSurfSharedPtr create()
    {
        return MemoryManager<CADSurfCFI>::AllocateSharedPtr();
    }

    static std::string key;

    CADSurfCFI(){};

    ~CADSurfCFI(){};

    void Initialise(int i, cfi::Face *in, NekDouble s);
    void SetScaling(NekDouble i)
    {
        m_scal = i;
    }

    Array<OneD, NekDouble> GetBounds();

    void GetBounds(NekDouble &umin, NekDouble &umax, NekDouble &vmin,
                           NekDouble &vmax);

    Array<OneD, NekDouble> N(Array<OneD, NekDouble> uv);

    Array<OneD, NekDouble> D1(Array<OneD, NekDouble> uv);

    Array<OneD, NekDouble> D2(Array<OneD, NekDouble> uv);

    Array<OneD, NekDouble> P(Array<OneD, NekDouble> uv);

    void P(Array<OneD, NekDouble> uv, NekDouble &x, NekDouble &y, NekDouble &z);

    Array<OneD, NekDouble> locuv(Array<OneD, NekDouble> p, NekDouble &dist);

    NekDouble Curvature(Array<OneD, NekDouble> uv);

    Array<OneD, NekDouble> BoundingBox()
    {
        ASSERTL0(false, "Not implemented in CFI");
        return Array<OneD, NekDouble>();
    }

    bool IsPlanar()
    {
        ASSERTL0(false, "Not implemented in CFI");
        return false;
    }

    cfi::Face *GetCfiPointer()
    {
        return m_cfiSurface;
    }

private:
    /// Function which tests the the value of uv used is within the surface
    void Test(Array<OneD, NekDouble> uv)
    {
        boost::ignore_unused(uv);
    }
    /// CFI object for surface.
    cfi::Face *m_cfiSurface;
    NekDouble m_scal;
};
}
}

#endif
