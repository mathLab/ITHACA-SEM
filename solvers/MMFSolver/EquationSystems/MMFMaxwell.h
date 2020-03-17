///////////////////////////////////////////////////////////////////////////////
//
// File MMFMaxwell.h
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
// Description: MMF Maxwell solve routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFMAXWELL_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFMAXWELL_H

#include <SolverUtils/MMFSystem.h>
#include <SolverUtils/UnsteadySystem.h>

enum CloakType
{
    eNoCloak,
    eOpticalCloak,
    eOpticalConstCloak,
    eOpticalDispersiveCloak,
    eMicroWaveCloak,
    SIZE_CloakType
};

const char *const CloakTypeMap[] = {
    "NoCloak",           "OpticalCloak",
    "OpticalConstCloak", "OpticalDispersiveCloak",
    "MicroWaveCloak",
};

enum SourceType
{
    eNoSource,
    ePointSource,
    ePlanarSource,
    SIZE_SourceType
};

const char *const SourceTypeMap[] = {
    "NoSource", "PointSource", "PlanarSource",
};

namespace Nektar
{
class MMFMaxwell : public SolverUtils::MMFSystem
{
public:
    friend class MemoryManager<MMFMaxwell>;

    CloakType m_CloakType;
    SourceType m_SourceType;
    bool m_DispersiveCloak;

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<MMFMaxwell>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }
    /// Name of class
    static std::string className;

    /// Initialise the object
    virtual void v_InitObject();

    virtual void v_DoSolve();

    /// Destructor
    virtual ~MMFMaxwell();

protected:
    int m_ElemtGroup0;
    int m_ElemtGroup1;
    int m_boundaryforSF;
    int m_PrintoutSurfaceCurrent;

    int m_AddPML;
    int m_PMLorder;

    int m_AddRotation;

    bool m_Cloaking;
    NekDouble m_CloakNlayer;
    NekDouble m_Cloakraddelta;

    NekDouble m_wp2Tol;
    Array<OneD, NekDouble> m_wp2;

    Array<OneD, NekDouble> m_SourceVector;
    NekDouble m_Psx, m_Psy, m_Psz;
    NekDouble m_PSduration, m_Gaussianradius;

    Array<OneD, Array<OneD, NekDouble>> m_CrossProductMF;

    /// Session reader
    MMFMaxwell(const LibUtilities::SessionReaderSharedPtr &pSession,
               const SpatialDomains::MeshGraphSharedPtr& pGraph);

    NekDouble m_freq;

    NekDouble m_n1, m_n2, m_n3;
    Array<OneD, NekDouble> m_varepsilon;
    Array<OneD, NekDouble> m_mu;

    int m_TestPML;
    int m_PMLelement, m_RecPML;
    NekDouble m_PMLthickness, m_PMLstart, m_PMLmaxsigma;
    Array<OneD, Array<OneD, NekDouble>> m_SigmaPML;

    int m_NoInc;

    Array<OneD, NekDouble> m_coriolis;

    /// Compute the RHS
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &inarray,
                  Array<OneD, Array<OneD, NekDouble>> &outarray,
                  const NekDouble time);

    /// Compute the projection
    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time);

    void AddGreenDerivCompensate(
        const Array<OneD, const Array<OneD, NekDouble>> &physarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray);

    void WeakDGMaxwellDirDeriv(
        const Array<OneD, const Array<OneD, NekDouble>> &InField,
        Array<OneD, Array<OneD, NekDouble>> &OutField,
        const NekDouble time = 0.0);

    NekDouble ComputeEnergyDensity(Array<OneD, Array<OneD, NekDouble>> &fields);

    Array<OneD, NekDouble> TestMaxwell1D(const NekDouble time,
                                         unsigned int field);
    Array<OneD, NekDouble> TestMaxwell2DPEC(
        const NekDouble time, unsigned int field,
        const SolverUtils::PolType Polarization);

    Array<OneD, NekDouble> TestMaxwell2DPMC(
        const NekDouble time, unsigned int field,
        const SolverUtils::PolType Polarization);

    Array<OneD, NekDouble> TestMaxwellSphere(const NekDouble time,
                                             const NekDouble omega,
                                             unsigned int field);

    void Printout_SurfaceCurrent(Array<OneD, Array<OneD, NekDouble>> &fields,
                                 const int time);

    Array<OneD, NekDouble> ComputeSurfaceCurrent(
        const int time,
        const Array<OneD, const Array<OneD, NekDouble>> &fields);

    void GenerateSigmaPML(const NekDouble PMLthickness,
                          const NekDouble PMLstart, const NekDouble PMLmaxsigma,
                          Array<OneD, Array<OneD, NekDouble>> &SigmaPML);

    void ComputeMaterialVector(Array<OneD, Array<OneD, NekDouble>> &epsvec,
                               Array<OneD, Array<OneD, NekDouble>> &muvec);

    void ComputeMaterialOpticalCloak(
        const Array<OneD, const NekDouble> &radvec,
        Array<OneD, Array<OneD, NekDouble>> &epsvec,
        Array<OneD, Array<OneD, NekDouble>> &muvec,
        const bool Dispersion = false);

    void ComputeMaterialMicroWaveCloak(
        const Array<OneD, const NekDouble> &radvec,
        Array<OneD, Array<OneD, NekDouble>> &epsvec,
        Array<OneD, Array<OneD, NekDouble>> &muvec);

    void Checkpoint_TotalFieldOutput(
        const int n, const NekDouble time,
        const Array<OneD, const Array<OneD, NekDouble>> &fieldphys);

    void Checkpoint_PlotOutput(
        const int n,
        const Array<OneD, const Array<OneD, NekDouble>> &fieldphys);

    void Checkpoint_TotPlotOutput(
        const int n, const NekDouble time,
        const Array<OneD, const Array<OneD, NekDouble>> &fieldphys);

    void Checkpoint_EDFluxOutput(
        const int n, const NekDouble time,
        const Array<OneD, const Array<OneD, NekDouble>> &fieldphys);

    void Checkpoint_EnergyOutput(
        const int n, const NekDouble time,
        const Array<OneD, const Array<OneD, NekDouble>> &fieldphys);

    Array<OneD, NekDouble> GaussianPulse(const NekDouble time,
                                         const NekDouble Psx,
                                         const NekDouble Psy,
                                         const NekDouble Psz,
                                         const NekDouble Gaussianradius);

    void AddPML(const Array<OneD, const Array<OneD, NekDouble>> &physarray,
                Array<OneD, Array<OneD, NekDouble>> &outarray);

    Array<OneD, NekDouble> EvaluateCoriolis();
    void AddCoriolis(Array<OneD, Array<OneD, NekDouble>> &physarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray);

    Array<OneD, NekDouble> ComputeRadCloak(const int Nlayer = 5);

    /// Print Summary
    virtual void v_GenerateSummary(SolverUtils::SummaryList &s);

    virtual void v_SetInitialConditions(const NekDouble initialtime,
                                        bool dumpInitialConditions,
                                        const int domain);

    virtual void v_EvaluateExactSolution(unsigned int field,
                                         Array<OneD, NekDouble> &outfield,
                                         const NekDouble time);
    void print_MMF(Array<OneD, Array<OneD, NekDouble>> &inarray);

private:
};
}

#endif
