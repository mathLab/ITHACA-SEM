////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSurf.h
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
//  Description: CAD object surface.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NekMeshUtils_CADSYSTEM_CFI_CADSURFCFI
#define NekMeshUtils_CADSYSTEM_CFI_CADSURFCFI

#include "../CADSurf.h"
#include "CADSystemCFI.h"

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

    CADSurfCFI() {};

    ~CADSurfCFI(){};

    void Initialise(int i, cfi::Face* in, std::vector<EdgeLoop> ein);
    void SetScaling(NekDouble i)
    {
        m_scal = i;
    }

    Array<OneD, NekDouble> GetBounds();

    Array<OneD, NekDouble> N    (Array<OneD, NekDouble> uv);

    Array<OneD, NekDouble> D1   (Array<OneD, NekDouble> uv);

    Array<OneD, NekDouble> D2   (Array<OneD, NekDouble> uv);

    Array<OneD, NekDouble> P    (Array<OneD, NekDouble> uv);

    Array<OneD, NekDouble> locuv(Array<OneD, NekDouble> p);

    NekDouble DistanceTo(Array<OneD, NekDouble> p);

    void ProjectTo(Array<OneD, NekDouble> &tp, Array<OneD, NekDouble> &uv);

    NekDouble Curvature(Array<OneD, NekDouble> uv);

private:
    /// Function which tests the the value of uv used is within the surface
    void Test(Array<OneD, NekDouble> uv){}
    /// CFI object for surface.
    cfi::Face* m_cfiSurface;
    NekDouble m_scal;
};

}
}

#endif
