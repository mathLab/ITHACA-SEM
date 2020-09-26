///////////////////////////////////////////////////////////////////////////////
//
// File APE.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2018 Kilian Lackhove
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
// Description: APE1/APE4 (Acoustic Perturbation Equations)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ACOUSTICSOLVER_EQUATIONSYSTEMS_APE_H
#define NEKTAR_SOLVERS_ACOUSTICSOLVER_EQUATIONSYSTEMS_APE_H

#include <AcousticSolver/EquationSystems/AcousticSystem.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{

class APE : public AcousticSystem
{
public:
    friend class MemoryManager<APE>;

    /// Creates an instance of this class
    static EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        EquationSystemSharedPtr p =
            MemoryManager<APE>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }
    /// Name of class
    static std::string className;

    /// Destructor
    virtual ~APE();

protected:
    /// Initialises UnsteadySystem class members.
    APE(const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph);

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
