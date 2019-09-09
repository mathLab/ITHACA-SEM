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

#ifndef NEKMESHUTILS_CADSYSTEM_OCE_CADVERTOCE
#define NEKMESHUTILS_CADSYSTEM_OCE_CADVERTOCE

#include <NekMeshUtils/CADSystem/CADVert.h>
#include <NekMeshUtils/CADSystem/OCE/OpenCascade.h>

namespace Nektar
{
namespace NekMeshUtils
{

class CADVertOCE : public CADVert
{
public:

    static CADVertSharedPtr create()
    {
        return MemoryManager<CADVertOCE>::AllocateSharedPtr();
    }

    static std::string key;

    /**
     * @brief Default constructor.
     */
    CADVertOCE()
    {
    }

    ~CADVertOCE()
    {
    }

    void Initialise(int i, TopoDS_Shape in);
    
    NekDouble DistanceTo(Array<OneD, NekDouble> l)
    {
        gp_Pnt lp(l[0],l[1],l[2]);
        return m_occVert.Distance(lp);
    }

private:
    /// OpenCascade object of the curve.
    gp_Pnt m_occVert;
};

}
}

#endif
