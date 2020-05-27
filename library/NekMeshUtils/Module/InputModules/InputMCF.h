////////////////////////////////////////////////////////////////////////////////
//
//  File: InputCAD.h
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
//  Description: Create mesh from CAD.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef UTILITIES_NEKMESH_INPUTMCF
#define UTILITIES_NEKMESH_INPUTMCF

#include <NekMeshUtils/Module/Module.h>

namespace Nektar
{
namespace Utilities
{

class InputMCF : public NekMeshUtils::InputModule
{
public:
    InputMCF(NekMeshUtils::MeshSharedPtr m);
    virtual ~InputMCF();
    virtual void Process();

    /// Creates an instance of this class
    static NekMeshUtils::ModuleSharedPtr create(NekMeshUtils::MeshSharedPtr m)
    {
        return MemoryManager<InputMCF>::AllocateSharedPtr(m);
    }
    /// %ModuleKey for class.
    static NekMeshUtils::ModuleKey className;

    void ParseFile(std::string nm);

private:
    std::string m_minDelta, m_maxDelta, m_eps, m_cadfile, m_order, m_blsurfs,
        m_blthick, m_blprog, m_bllayers, m_refinement, m_nacadomain, m_periodic,
        m_adjustment, m_spaceoutblthr, m_nospaceoutsurf, m_voidPts;

    bool m_makeBL, m_surfopti, m_varopti, m_refine, m_woct, m_2D, m_splitBL,
        m_naca, m_adjust, m_adjustall, m_smoothbl, m_manifold, m_spaceoutbl;
};
}
}

#endif
