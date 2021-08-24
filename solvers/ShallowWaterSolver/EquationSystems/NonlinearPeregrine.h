///////////////////////////////////////////////////////////////////////////////
//
// File NonlinearPeregrine.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
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
// Description: Nonlinear Boussinesq equations of Peregrine in conservative
//              variables
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_NONLINEARPEREGRINE_H
#define NEKTAR_SOLVERS_SHALLOWWATERSOLVER_EQUATIONSYSTEMS_NONLINEARPEREGRINE_H

#include <ShallowWaterSolver/EquationSystems/ShallowWaterSystem.h>

namespace Nektar
{

enum ProblemType
{
    eGeneral,          ///< No problem defined - Default Inital data
    eSolitaryWave,     ///< First order Laitone solitary wave
    SIZE_ProblemType   ///< Length of enum list
};

const char* const ProblemTypeMap[] = { "General", "SolitaryWave" };

/**
 *
 *
 **/

class NonlinearPeregrine: public ShallowWaterSystem
{
public:
    friend class MemoryManager<NonlinearPeregrine> ;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p = MemoryManager<
            NonlinearPeregrine>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }
    /// Name of class
    static std::string className;

    virtual ~NonlinearPeregrine();

    ///< problem type selector
    ProblemType m_problemType;

protected:
    StdRegions::ConstFactorMap m_factors;

    NonlinearPeregrine(const LibUtilities::SessionReaderSharedPtr& pSession,
                       const SpatialDomains::MeshGraphSharedPtr& pGraph);

    virtual void v_InitObject();

    /// Still water depth traces
    Array<OneD, NekDouble> m_dFwd;
    Array<OneD, NekDouble> m_dBwd;

    void DoOdeRhs(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble time);

    void DoOdeProjection(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble time);

    void GetFluxVector(
            const Array<OneD, const Array<OneD, NekDouble> > &physfield,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux);

    virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

    virtual void v_PrimitiveToConservative();

    virtual void v_ConservativeToPrimitive();

    virtual void v_SetInitialConditions(
            NekDouble initialtime = 0.0,
            bool dumpInitialConditions = true,
            const int domain = 0);

    const Array<OneD, NekDouble> &GetDepthFwd()
    {
        return m_dFwd;
    }
    const Array<OneD, NekDouble> &GetDepthBwd()
    {
        return m_dBwd;
    }

private:
    NekDouble m_const_depth;

    void NumericalFlux1D(
            Array<OneD, Array<OneD, NekDouble> > &physfield,
            Array<OneD, Array<OneD, NekDouble> > &numfluxX);

    void NumericalFlux2D(
            Array<OneD, Array<OneD, NekDouble> > &physfield,
            Array<OneD, Array<OneD, NekDouble> > &numfluxX,
            Array<OneD, Array<OneD, NekDouble> > &numfluxY);

    void LaitoneSolitaryWave(
            NekDouble amp,
            NekDouble d,
            NekDouble time,
            NekDouble x_offset);

    void SetBoundaryConditions(
            Array<OneD, Array<OneD, NekDouble> > &physarray,
            NekDouble time);

    void WallBoundary2D(
            int bcRegion,
            int cnt,
            Array<OneD, Array<OneD, NekDouble> > &Fwd,
            Array<OneD, Array<OneD, NekDouble> > &physarray);

    void WallBoundary(
            int bcRegion,
            int cnt,
            Array<OneD, Array<OneD, NekDouble> > &Fwd,
            Array<OneD, Array<OneD, NekDouble> > &physarray);

    void AddCoriolis(
            const Array<OneD, const Array<OneD, NekDouble> > &physarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray);

    void AddVariableDepth(
            const Array<OneD, const Array<OneD, NekDouble> > &physarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray);

    void ConservativeToPrimitive(
            const Array<OneD, const Array<OneD, NekDouble> >&physin,
            Array<OneD, Array<OneD, NekDouble> >&physout);

    void PrimitiveToConservative(
            const Array<OneD, const Array<OneD, NekDouble> >&physin,
            Array<OneD, Array<OneD, NekDouble> >&physout);

    void GetVelocityVector(
            const Array<OneD, Array<OneD, NekDouble> > &physfield,
            Array<OneD, Array<OneD, NekDouble> > &velocity);

    // Dispersive parts
    void WCESolve(Array<OneD, NekDouble> &fce, NekDouble lambda);

    void NumericalFluxForcing(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            Array<OneD, NekDouble> &numfluxX,
            Array<OneD, NekDouble> &numfluxY);

    void SetBoundaryConditionsForcing(
            Array<OneD, Array<OneD, NekDouble> > &inarray, NekDouble time);

    void WallBoundaryForcing(
            int bcRegion,
            int cnt,
            Array<OneD, Array<OneD, NekDouble> >&inarray);

    void SetBoundaryConditionsContVariables(
            Array<OneD, NekDouble> &inarray,
            NekDouble time);

    void WallBoundaryContVariables(
            int bcRegion,
            int cnt,
            Array<OneD, NekDouble>&inarray);

    void NumericalFluxConsVariables(
            Array<OneD, NekDouble> &physfield,
            Array<OneD, NekDouble> &outX,
            Array<OneD, NekDouble> &outY);

};

}

#endif

