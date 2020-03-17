///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionTerm.h
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
// Description: Driver class for the stability solver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DRIVERSTEADYSTATE_H
#define NEKTAR_SOLVERUTILS_DRIVERSTEADYSTATE_H

#include <SolverUtils/Driver.h>
#include <SolverUtils/DriverArnoldi.h>
#include <SolverUtils/DriverModifiedArnoldi.h>
#include <LibUtilities/BasicUtils/Timer.h>

namespace Nektar
{
namespace SolverUtils
{

class DriverSteadyState: public DriverModifiedArnoldi
{
public:
    friend class MemoryManager<DriverSteadyState>;

    /// Creates an instance of this class
    static DriverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
    {
        DriverSharedPtr p = MemoryManager<DriverSteadyState>
            ::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    ///Name of the class
    static std::string className;

    void PrintSummarySFD();

    void ConvergenceHistory(
            const Array<OneD, const Array<OneD, NekDouble> > &qBar1,
            const Array<OneD, const Array<OneD, NekDouble> > &q0,
                  NekDouble &MaxNormDiff_q_qBar,
                  NekDouble &MaxNormDiff_q1_q0);

    void ComputeSFD(
            const int i,
            const Array<OneD, const Array<OneD, NekDouble> > &q0,
            const Array<OneD, const Array<OneD, NekDouble> > &qBar0,
                  Array<OneD, Array<OneD, NekDouble> > &q1,
                  Array<OneD, Array<OneD, NekDouble> > &qBar1);

    void EvalEV_ScalarSFD(
            const NekDouble &X_input,
            const NekDouble &Delta_input,
            const std::complex<NekDouble> &alpha,
                  NekDouble &MaxEV);

    void GradientDescentMethod(
            const std::complex<NekDouble> &alpha,
                  NekDouble &X_output,
                  NekDouble &Delta_output);

    void ReadEVfile(
                  int &KrylovSubspaceDim,
                  NekDouble &growthEV,
                  NekDouble &frequencyEV);

    void ComputeOptimization();

    void SetSFDOperator(
            const NekDouble X_input,
            const NekDouble Delta_input);


protected:
    /// Constructor
    SOLVER_UTILS_EXPORT DriverSteadyState(
        const LibUtilities::SessionReaderSharedPtr pSession,
        const SpatialDomains::MeshGraphSharedPtr pGraph);

    /// Destructor
    SOLVER_UTILS_EXPORT virtual ~DriverSteadyState();

    /// Second-stage initialisation
    SOLVER_UTILS_EXPORT virtual void v_InitObject(std::ostream &out = std::cout);

    /// Virtual function for solve implementation.
    SOLVER_UTILS_EXPORT virtual void v_Execute(std::ostream &out = std::cout);

    static std::string driverLookupId;

private:
    int m_stepCounter;
    int m_Check;
    int m_Check_BaseFlow;
    NekDouble TOL;
    int m_infosteps;
    int m_checksteps;
    int NumVar_SFD;

    LibUtilities::Timer     timer;
    NekDouble cpuTime;
    NekDouble totalTime;
    NekDouble elapsed;
    std::ofstream m_file;

    ///Definition of the SFD parameters:
    NekDouble m_Delta;
    NekDouble m_X;
    NekDouble m_dt;

    ///Definition of the SFD operator
    NekDouble M11;
    NekDouble M12;
    NekDouble M21;
    NekDouble M22;

    NekDouble Diff_q_qBar;
    NekDouble Diff_q1_q0;

    ///For adaptive SFD method
    bool FlowPartiallyConverged;
    NekDouble AdaptiveTOL;
    NekDouble AdaptiveTime;
    int m_NonConvergingStepsCounter;
    NekDouble GrowthRateEV;
    NekDouble FrequencyEV;
};

}
} //end of namespace


#endif //NEKTAR_SOLVERUTILS_DRIVERSTEADYSTATE_H

