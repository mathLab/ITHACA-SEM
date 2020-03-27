///////////////////////////////////////////////////////////////////////////////
//
// File MMFSWE.h
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
// Description: MMF advection solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFSWE_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFSWE_H

#include <ShallowWaterSolver/EquationSystems/ShallowWaterSystem.h>
#include <SolverUtils/MMFSystem.h>

namespace Nektar
{

enum TestType
{
    eTestPlane,
    eTestSteadyZonal,
    eTestUnsteadyZonal,
    eTestIsolatedMountain,
    eTestUnstableJet,
    eTestRossbyWave,
    SIZE_TestType ///< Length of enum list
};

const char *const TestTypeMap[] = {
    "TestPlane",         "TestSteadyZonal",
    "TestUnsteadyZonal", "TestIsolatedMountain",
    "TestUnstableJet",   "TestRossbyWave",
};

class MMFSWE : public SolverUtils::MMFSystem
{
public:
    friend class MemoryManager<MMFSWE>;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<MMFSWE>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }
    /// Name of class
    static std::string className;

    TestType m_TestType;

    /// Destructor
    virtual ~MMFSWE();

protected:
    /// Still water depth
    Array<OneD, NekDouble> m_depth;
    Array<OneD, Array<OneD, NekDouble>> m_Derivdepth;

    /// Coriolis force
    Array<OneD, NekDouble> m_coriolis;

    int m_AddCoriolis, m_AddRotation;

    NekDouble m_g;                                  // Gravity
    NekDouble m_alpha, m_u0, m_Omega, m_H0, m_Hvar; // TestSteadyZonal
    NekDouble m_k2;                                 // eTestUnsteadyZonal
    NekDouble m_hs0;                                // eTestIsolatedMountain
    NekDouble m_uthetamax, m_theta0, m_theta1, m_en, m_hbar; // eTestUnstableJet
    NekDouble m_angfreq, m_K;                                // eTestRossbyWave

    NekDouble m_Vorticity0, m_Mass0, m_Energy0, m_Enstrophy0;

    /// Indicates if variables are primitive or conservative
    bool m_primitive;

    /// Advection velocity
    Array<OneD, Array<OneD, NekDouble>> m_velocity;
    Array<OneD, NekDouble> m_traceVn;

    Array<OneD, Array<OneD, NekDouble>> m_veldotMF;
    Array<OneD, NekDouble> m_vellc;

    // Plane (used only for Discontinous projection
    //        with 3DHomogenoeus1D expansion)
    int m_planeNumber;

    /// Session reader
    MMFSWE(const LibUtilities::SessionReaderSharedPtr &pSession,
           const SpatialDomains::MeshGraphSharedPtr& pGraph);

    /// Compute the RHS
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    /// Compute the projection
    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void SteadyZonalFlow(unsigned int field, Array<OneD, NekDouble> &outfield);

    void UnsteadyZonalFlow(unsigned int field, const NekDouble time,
                           Array<OneD, NekDouble> &outfield);

    void IsolatedMountainFlow(unsigned int field, const NekDouble time,
                              Array<OneD, NekDouble> &outfield);

    void UnstableJetFlow(unsigned int field, const NekDouble time,
                         Array<OneD, NekDouble> &outfield);

    void RossbyWave(unsigned int field, Array<OneD, NekDouble> &outfield);

    NekDouble ComputeUnstableJetEta(const NekDouble theta);

    NekDouble ComputeUnstableJetuphi(const NekDouble theta);

    NekDouble ComputeMass(const Array<OneD, const NekDouble> &eta);

    NekDouble ComputeEnergy(const Array<OneD, const NekDouble> &eta,
                            const Array<OneD, const NekDouble> &u,
                            const Array<OneD, const NekDouble> &v);

    NekDouble ComputeEnstrophy(const Array<OneD, const NekDouble> &eta,
                               const Array<OneD, const NekDouble> &u,
                               const Array<OneD, const NekDouble> &v);

    void ComputeVorticity(const Array<OneD, const NekDouble> &u,
                          const Array<OneD, const NekDouble> &v,
                          Array<OneD, NekDouble> &Vorticity);

    /// Get the normal velocity
    void GetNormalVelocity(Array<OneD, NekDouble> &traceVn);

    void ComputeNablaCdotVelocity(Array<OneD, NekDouble> &vellc);

    void PrimitiveToConservative(
        const Array<OneD, const Array<OneD, NekDouble>> &physin,
        Array<OneD, Array<OneD, NekDouble>> &physout);

    void ConservativeToPrimitive(
        const Array<OneD, const Array<OneD, NekDouble>> &physin,
        Array<OneD, Array<OneD, NekDouble>> &physout);

    void PrimitiveToConservative();
    void ConservativeToPrimitive();

    void WeakDGSWEDirDeriv(const Array<OneD, Array<OneD, NekDouble>> &InField,
                           Array<OneD, Array<OneD, NekDouble>> &OutField);

    void AddDivForGradient(Array<OneD, Array<OneD, NekDouble>> &physarray,
                           Array<OneD, Array<OneD, NekDouble>> &outarray);

    void NumericalSWEFlux(Array<OneD, Array<OneD, NekDouble>> &physfield,
                          Array<OneD, Array<OneD, NekDouble>> &numfluxFwd,
                          Array<OneD, Array<OneD, NekDouble>> &numfluxBwd);

    void GetSWEFluxVector(
        const int i, const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &flux);

    void RiemannSolverHLLC(const int index, NekDouble hL, NekDouble uL,
                           NekDouble vL, NekDouble hR, NekDouble uR,
                           NekDouble vR, Array<OneD, NekDouble> &numfluxF,
                           Array<OneD, NekDouble> &numfluxB);

    void Computehhuhvflux(NekDouble hL, NekDouble uL, NekDouble vL,
                          NekDouble hR, NekDouble uR, NekDouble vR,
                          NekDouble hstar, NekDouble &hflux, NekDouble &huflux,
                          NekDouble &hvflux);

    void AverageFlux(const int index, NekDouble hL, NekDouble uL, NekDouble vL,
                     NekDouble hR, NekDouble uR, NekDouble vR,
                     Array<OneD, NekDouble> &numfluxF,
                     Array<OneD, NekDouble> &numfluxB);

    void LaxFriedrichFlux(const int index, NekDouble hL, NekDouble uL,
                          NekDouble vL, NekDouble hR, NekDouble uR,
                          NekDouble vR, Array<OneD, NekDouble> &numfluxF,
                          Array<OneD, NekDouble> &numfluxB);

    void RusanovFlux(const int index, NekDouble hL, NekDouble uL, NekDouble vL,
                     NekDouble hR, NekDouble uR, NekDouble vR,
                     Array<OneD, NekDouble> &numfluxF,
                     Array<OneD, NekDouble> &numfluxB);

    void ComputeMagAndDot(const int index, NekDouble &MageF1, NekDouble &MageF2,
                          NekDouble &MageB1, NekDouble &MageB2,
                          NekDouble &eF1_cdot_eB1, NekDouble &eF1_cdot_eB2,
                          NekDouble &eF2_cdot_eB1, NekDouble &eF2_cdot_eB2);

    void AddCoriolis(Array<OneD, Array<OneD, NekDouble>> &physarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray);

    void AddElevationEffect(Array<OneD, Array<OneD, NekDouble>> &physarray,
                            Array<OneD, Array<OneD, NekDouble>> &outarray);

    void AddRotation(Array<OneD, Array<OneD, NekDouble>> &physarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray);

    void Compute_demdt_cdot_ek(
        const int indm, const int indk,
        const Array<OneD, const Array<OneD, NekDouble>> &physarray,
        Array<OneD, NekDouble> &outarray);

    void TestVorticityComputation();

    void EvaluateWaterDepth(void);

    void EvaluateCoriolis(void);
    void EvaluateCoriolisForZonalFlow(Array<OneD, NekDouble> &outarray);
    void EvaluateStandardCoriolis(Array<OneD, NekDouble> &outarray);

    void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble>> &inarray,
                               NekDouble time);
    void WallBoundary2D(int bcRegion, int cnt,
                        Array<OneD, Array<OneD, NekDouble>> &physarray);

    /// Initialise the object
    virtual void v_InitObject();

    virtual void v_DoSolve();

    virtual void v_DoInitialise();

    /// Print Summary
    virtual void v_GenerateSummary(SolverUtils::SummaryList &s);

    virtual void v_SetInitialConditions(const NekDouble initialtime,
                                        bool dumpInitialConditions,
                                        const int domain);

    virtual void v_EvaluateExactSolution(unsigned int field,
                                         Array<OneD, NekDouble> &outfield,
                                         const NekDouble time);

    virtual NekDouble v_L2Error(unsigned int field,
                                const Array<OneD, NekDouble> &exactsoln,
                                bool Normalised);

    virtual NekDouble v_LinfError(unsigned int field,
                                  const Array<OneD, NekDouble> &exactsoln);

private:
    int m_RossbyDisturbance;
    int m_PurturbedJet;

    void TestSWE2Dproblem(const NekDouble time, unsigned int field,
                          Array<OneD, NekDouble> &outfield);

    void Checkpoint_Output_Cartesian(std::string outname);
};
}

#endif
