///////////////////////////////////////////////////////////////////////////////
//
// File LEE.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Linearized Euler Equations
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_LEESOLVER_EQUATIONSYSTEMS_LEE_H
#define NEKTAR_SOLVERS_LEESOLVER_EQUATIONSYSTEMS_LEE_H

#include <AcousticSolver/EquationSystems/AcousticSystem.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{

class LEE : public AcousticSystem
{
public:
    friend class MemoryManager<LEE>;

    /// Creates an instance of this class
    static EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        EquationSystemSharedPtr p =
            MemoryManager<LEE>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }
    /// Name of class
    static std::string className;

    /// Destructor
    virtual ~LEE();

protected:
    /// Initialises UnsteadySystem class members.
    LEE(const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph);

    virtual void v_AddLinTerm(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    virtual void v_InitObject();

    virtual void v_GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);

private:
    virtual void v_RiemannInvariantBC(
        int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble>> &Fwd,
        Array<OneD, Array<OneD, NekDouble>> &BfFwd,
        Array<OneD, Array<OneD, NekDouble>> &physarray);
};
} // namespace Nektar

#endif
