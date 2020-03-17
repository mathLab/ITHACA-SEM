///////////////////////////////////////////////////////////////////////////////
//
// File MMFDiffusion.h
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
// Description: MMFDiffusion
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFDIFFUSION_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFDIFFUSION_H

#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/MMFSystem.h>

//using namespace Nektar::SolverUtils;

namespace Nektar
{

  enum TestType
  {           
    eTestPlane,
    eTestCube,
    eTestLinearSphere,     
    eTestNonlinearSphere,
    eFHNStandard,
    eFHNRogers,
    eFHNAlievPanf,
    SIZE_TestType   ///< Length of enum list
  };
  
  const char* const TestTypeMap[] =
    {
      "TestPlane",
      "TestCube",
      "TestLinearSphere",
      "TestNonlinearSphere",
      "FHNStandard",
      "FHNRogers",
      "FHNAlievPanf",
    };

  enum InitWaveType
  {          
    eLeft,
    eBothEnds,
    eCenter,
    eLeftBottomCorner,
    ePoint,
    eSpiralDock,
    SIZE_InitWaveType   ///< Length of enum list
  };

  const char* const InitWaveTypeMap[] =
    {
      "Left",
      "BothEnd",
      "Center",
      "LeftBottomCorner",
      "Point",
      "SpiralDock",
    };
 
    /// A model for cardiac conduction.
    class MMFDiffusion : public SolverUtils::MMFSystem
    {
    public:
        friend class MemoryManager<MMFDiffusion>;
	
        /// Creates an instance of this class
	    static SolverUtils::EquationSystemSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const SpatialDomains::MeshGraphSharedPtr& pGraph)
        {
	  SolverUtils::EquationSystemSharedPtr p = MemoryManager<MMFDiffusion>::AllocateSharedPtr(pSession, pGraph);
            p->InitObject();
            return p;
        }

        /// Name of class
        static std::string className;

	TestType                                     m_TestType; 

        /// Desctructor
        virtual ~MMFDiffusion();

    protected:
        /// Constructor
        MMFDiffusion(const LibUtilities::SessionReaderSharedPtr& pSession,
                    const SpatialDomains::MeshGraphSharedPtr& pGraph);

	InitWaveType                                    m_InitWaveType;

	
        virtual void v_InitObject();

        /// Solve for the diffusion term.
        void DoImplicitSolve(
                const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                      Array<OneD, Array<OneD, NekDouble> >&outarray,
                      NekDouble time,
                      NekDouble lambda);

        /// Computes the reaction terms \f$f(u,v)\f$ and \f$g(u,v)\f$.
        void DoOdeRhs(
                const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                      Array<OneD,        Array<OneD, NekDouble> >&outarray,
                const NekDouble time);

	void TestPlaneProblem(const NekDouble time,
			      Array<OneD, NekDouble> &outfield);

	void TestCubeProblem(const NekDouble time,
			      Array<OneD, NekDouble> &outfield);

	void Morphogenesis(const NekDouble time,
			   unsigned int field,
			   Array<OneD, NekDouble> &outfield);

	Array<OneD, NekDouble> PlanePhiWave();

        /// Sets a custom initial condition.
        virtual void v_SetInitialConditions(NekDouble initialtime,
                                bool dumpInitialConditions,
                                const int domain);

        /// Prints a summary of the model parameters.
        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

	
	virtual void v_EvaluateExactSolution(unsigned int field,
					     Array<OneD, NekDouble> &outfield,
					     const NekDouble time);

	NekDouble m_InitPtx, m_InitPty, m_InitPtz;

    private:
        /// Variable diffusivity
        StdRegions::VarCoeffMap m_varcoeff;

	Array<OneD, NekDouble> m_epsilon;
	Array<OneD, NekDouble> m_epsu;
    };

}

#endif
