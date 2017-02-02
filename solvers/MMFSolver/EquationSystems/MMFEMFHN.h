///////////////////////////////////////////////////////////////////////////////
//
// File MMFEMFHN.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: MMF EMFHN solver routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFEMFHN_H
#define NEKTAR_SOLVERS_ADRSOLVER_EQUATIONSYSTEMS_MMFEMFHN_H

#include <SolverUtils/UnsteadySystem.h>
#include <SolverUtils/MMFSystem.h>

namespace Nektar
{
  enum TestType
  {          
    eTestEMFHN,
    eTestCoulombField,
    eFHNStandard,
    eRogers,
    eAlievPanf,
    eFentonKarma,
    SIZE_TestType   ///< Length of enum list
  };

  const char* const TestTypeMap[] =
    {
      "TestEMFHN",
      "TestCoulombField",
      "FHNStandard",
      "Rogers",
      "AlievPanf",
      "FentonKarma",
    };

  enum FentonKarmaType
  {
    eBR,
    eMBR,
    eMLR1,
    eGP,
    eCF1,
    eCF2a,
    eCF2b,
    eCF2c,
    eCF3a,
    eCF3b,
    SIZE_FentonKarmaType
  };
  
  const char* const FentonKarmaTypeMap[] =
    {
      "BR",
      "MBR",
      "MLR1",
      "GP",
      "CF1",
      "CF2a",
      "CF2b",
      "CF2c",
      "CF3a",
      "CF3b"
    };
  
  
  enum InitWaveType
  {          
    eLeft,
    eBothEnds,
    eCenter,
    eLeftTopCorner,
    ePoint,
    eSpiralDock,
    SIZE_InitWaveType   ///< Length of enum list
  };

  const char* const InitWaveTypeMap[] =
    {
      "Left",
      "BothEnd",
      "Center",
      "LeftTopCorner",
      "Point",
      "SpiralDock",
    };

    class MMFEMFHN : public SolverUtils::MMFSystem
    {
    public:
        friend class MemoryManager<MMFEMFHN>;

        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession) {
            SolverUtils::EquationSystemSharedPtr p
                = MemoryManager<MMFEMFHN>::AllocateSharedPtr(pSession);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;

        virtual void v_InitObject();
	virtual void v_DoSolve();

	TestType                                        m_TestType;
	InitWaveType                                    m_InitWaveType;
	FentonKarmaType                                 m_FentonKarmaType;

        /// Destructor
        virtual ~MMFEMFHN();
	
    protected:
	int m_E1, m_E2, m_D1, m_D2, m_Hz;
	int m_phi, m_psi, m_chi, m_Fs;
	int m_Gr, m_dFdt;
	
	NekDouble m_alpha;
	NekDouble m_Jac;
	NekDouble m_PoissonTau;
	NekDouble m_PoissonNeumann;
	
	int m_ElemtGroup1;
	int m_StimulusPeriod;
	NekDouble m_radiusofinit;

	// Fenton Karma Variables
	NekDouble m_Cm, m_V0;
	NekDouble m_u_fi, m_u_c, m_u_v, m_u_r, m_g_fi_max;
	NekDouble m_tau_d, m_tau_v1_minus, m_tau_v2_minus, m_tau_v_plus, m_tau_0;
	NekDouble m_tau_r, m_tau_si, m_tau_y_plus, m_tau_y_minus;
	NekDouble m_u_csi, m_k1, m_k2, m_kp;
	NekDouble m_tau_w_minus, m_tau_w_plus;
	
        /// Session reader
        MMFEMFHN(const LibUtilities::SessionReaderSharedPtr& pSession);

	StdRegions::VarCoeffMap m_varcoeff;

	Array<OneD, NekDouble> m_sigma;
	
	Array<OneD, Array<OneD, NekDouble> > m_SigmaBlock;

	Array<OneD, Array<OneD, NekDouble> > m_CrossProductMF;
	
	Array<OneD, NekDouble> m_varepsilon;
	Array<OneD, NekDouble> m_mu;

	Array<OneD, Array<OneD, NekDouble> > m_Vn;
	
	Array<OneD, Array<OneD, NekDouble> > m_ZimFwd;
	Array<OneD, Array<OneD, NekDouble> > m_ZimBwd;
	Array<OneD, Array<OneD, NekDouble> > m_YimFwd;
	Array<OneD, Array<OneD, NekDouble> > m_YimBwd;

	Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_ntimesMFFwd;
	Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_ntimesMFBwd;
	Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_ntimes_ntimesMFFwd;
	Array<OneD, Array<OneD, Array<OneD, NekDouble> > > m_ntimes_ntimesMFBwd;

	NekDouble m_InitPtx, m_InitPty, m_InitPtz;
	
        /// Compute the RHS
        void DoOdeRhs(
            const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                  Array<OneD,         Array<OneD, NekDouble> > &outarray,
            const NekDouble time);

        /// Compute the projection
        void DoOdeProjection(
            const Array<OneD,  const  Array<OneD, NekDouble> > &inarray,
                  Array<OneD,         Array<OneD, NekDouble> > &outarray,
            const NekDouble time);
	
	void WeakDGAdv(const int indx,
		       const int dir,
		       const Array<OneD, const NekDouble> &inarray,
		       Array<OneD, NekDouble> &outarray);

	Array<OneD, NekDouble> FlowTestEMFHN(const NekDouble time, const int var);

	void SolvePoisson(const Array<OneD, const NekDouble> &inarray,
			  const NekDouble time = 0.0);

	void SetBoundaryConditions(const int var, NekDouble time);

	void WeakDGMaxwellDirDeriv(const Array<OneD, const Array<OneD, NekDouble> >& InField,
				   Array<OneD, Array<OneD, NekDouble> >& OutField,
				   const NekDouble time);

	void DoImplicitSolve(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
			     Array<OneD, Array<OneD, NekDouble> >&outarray,
			     const NekDouble time,
			     const NekDouble lambda);

	void ComputeDRReaction(const Array<OneD, const NekDouble> &uarray,
			       const Array<OneD, const NekDouble> &varray,
			       const Array<OneD, const NekDouble> &warray,
			       Array<OneD, NekDouble> &uReaction,
			       Array<OneD, NekDouble> &vReaction,
			       Array<OneD, NekDouble> &wReaction);

	void ComputedFdt(const Array<OneD, const NekDouble> &uarray,
			 const Array<OneD, const NekDouble> &varray,
			 const Array<OneD, const NekDouble> &warray,
			 const Array<OneD, const NekDouble> &dudtarray,
			 const Array<OneD, const NekDouble> &dvdtarray,
			 const Array<OneD, const NekDouble> &dwdtarray,
			 Array<OneD, NekDouble> &outarray);

	void AddGreenDerivCompensate(const Array<OneD, const Array<OneD, NekDouble> > &physarray,
				     Array<OneD, Array<OneD, NekDouble> > &outarray);
	
	void ComputeVarcoeffDiffusion();

	void PlaneWaveForPhi(Array<OneD, NekDouble> &outarray);

	void ComputeEnergy(const Array<OneD, const Array<OneD, NekDouble> > &fields, 
			   Array<OneD, Array<OneD, NekDouble> > &energy);

	void ComputeAPD(const Array<OneD, const NekDouble> &inarray,
			Array<OneD, NekDouble> &outarray);

	void SetUpParametersFentonKarma();

	void Checkpoint_XYZOutput(const int n,
				  const Array<OneD, const Array<OneD, NekDouble> > &fieldphys);

	void Checkpoint_EnergyOutput(const int n,
				     const Array<OneD, const Array<OneD, NekDouble> > &fieldphys);


        /// Print Summary
        virtual void v_GenerateSummary(SolverUtils::SummaryList& s);

	virtual void v_SetInitialConditions(const NekDouble initialtime,
					    bool dumpInitialConditions,
					    const int domain);

	virtual void v_EvaluateExactSolution(unsigned int field,
					     Array<OneD, NekDouble> &outfield,
					     const NekDouble time);

    private:

    };
}

#endif
