///////////////////////////////////////////////////////////////////////////////
//
// File AcousticSystem.h
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
// Description: AcousticSystem base class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ACOUSTICSOLVER_EQUATIONSYSTEMS_ACOUSTICSYSTEM_H
#define NEKTAR_SOLVERS_ACOUSTICSOLVER_EQUATIONSYSTEMS_ACOUSTICSYSTEM_H

// Define variable to avoid deprecated warning in Boost 1.69.
#include <boost/version.hpp>
#if BOOST_VERSION >= 106900 && BOOST_VERSION < 107000
#define BOOST_ALLOW_DEPRECATED_HEADERS
#endif

#include <boost/core/ignore_unused.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <SolverUtils/Advection/Advection.h>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Core/Coupling.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include <SolverUtils/UnsteadySystem.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{

class AcousticSystem : public AdvectionSystem
{
public:
    friend class MemoryManager<AcousticSystem>;

    /// Destructor
    virtual ~AcousticSystem();

protected:
    /// indices of the fields
    int m_ip, m_irho, m_iu;
    /// we are dealing with a conservative formualtion
    bool m_conservative;
    SolverUtils::CouplingSharedPtr m_coupling;
    SolverUtils::AdvectionSharedPtr m_advection;
    std::vector<SolverUtils::ForcingSharedPtr> m_forcing;
    SolverUtils::RiemannSolverSharedPtr m_riemannSolver;
    Array<OneD, Array<OneD, NekDouble>> m_bfFwdBwd;
    Array<OneD, Array<OneD, NekDouble>> m_vecLocs;
    Array<OneD, Array<OneD, NekDouble>> m_bf;
    std::vector<std::string> m_bfNames;

    /// Initialises UnsteadySystem class members.
    AcousticSystem(const LibUtilities::SessionReaderSharedPtr &pSession,
                   const SpatialDomains::MeshGraphSharedPtr &pGraph);

    virtual void v_InitObject();

    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    virtual void v_GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux) = 0;

    virtual void v_AddLinTerm(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray)
    {
        boost::ignore_unused(inarray, outarray);
    }

    virtual bool v_PreIntegrate(int step);

    virtual void v_Output();

    virtual Array<OneD, NekDouble> v_GetMaxStdVelocity();

    virtual void v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
        std::vector<std::string> &variables);

    const Array<OneD, const Array<OneD, NekDouble>> &GetNormals();

    const Array<OneD, const Array<OneD, NekDouble>> &GetVecLocs();

    const Array<OneD, const Array<OneD, NekDouble>> &GetBasefieldFwdBwd();

private:
    std::map<int, boost::mt19937> m_rng;
    NekDouble m_whiteNoiseBC_lastUpdate;
    NekDouble m_whiteNoiseBC_p;

    NekDouble GetCFLEstimate();

    void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble>> &physarray,
                               NekDouble time);

    virtual void v_WallBC(int bcRegion, int cnt,
                          Array<OneD, Array<OneD, NekDouble>> &Fwd,
                          Array<OneD, Array<OneD, NekDouble>> &physarray);

    virtual void v_RiemannInvariantBC(
        int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble>> &Fwd,
        Array<OneD, Array<OneD, NekDouble>> &BfFwd,
        Array<OneD, Array<OneD, NekDouble>> &physarray) = 0;

    virtual void v_WhiteNoiseBC(int bcRegion, int cnt,
                                Array<OneD, Array<OneD, NekDouble>> &Fwd,
                                Array<OneD, Array<OneD, NekDouble>> &BfFwd,
                                Array<OneD, Array<OneD, NekDouble>> &physarray);

    void CopyBoundaryTrace(const Array<OneD, NekDouble> &Fwd,
                           Array<OneD, NekDouble> &Bwd);

    void UpdateBasefieldFwdBwd();
};
} // namespace Nektar

#endif
