///////////////////////////////////////////////////////////////////////////////
//
// File PulseWaveSystem.h
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
// Description: Generic timestepping for PulseWaveSolver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_PULSEWAVESOLVER_EQUATIONSYSTEMS_PULSEWAVESYSTEM_H
#define NEKTAR_SOLVERS_PULSEWAVESOLVER_EQUATIONSYSTEMS_PULSEWAVESYSTEM_H

#include <SolverUtils/UnsteadySystem.h>
#include <time.h>

using namespace Nektar::SolverUtils;

namespace Nektar
{
  
	enum UpwindTypePulse
	{           
		eNotSetPulse,            ///< flux not defined
		eUpwindPulse,			 ///< simple upwinding scheme
		SIZE_UpwindTypePulse     ///< Length of enum list
	};
	
	const char* const UpwindTypeMapPulse[] =
    {
		"NoSetPulse",
		"UpwindPulse",
	};
	
	
    /// Base class for unsteady solvers.
    class PulseWaveSystem : public UnsteadySystem
    {
    public:
        /// Destructor
        virtual ~PulseWaveSystem();
		
    protected:
		
		Array<OneD, MultiRegions::ExpListSharedPtr>     m_vessels;
		
		Array<OneD, MultiRegions::ExpListSharedPtr>     m_outfields;

		int												m_domainsize;
		
		int												m_omega;

		UpwindTypePulse                                 m_upwindTypePulse;
		
		NekDouble                                       m_rho;
		
		NekDouble                                       m_pext;
		
		NekDouble                                       m_nue;
		
		NekDouble                                       m_h0;
		
		NekDouble m_C;
		NekDouble m_RT;
		NekDouble m_pout;
		
		Array<OneD, Array<OneD, NekDouble> >			m_A_0global;
		
		Array<OneD, NekDouble>							m_A_0;

		Array<OneD, Array<OneD, NekDouble> >			m_betaglobal;
		
		Array<OneD, NekDouble>							m_beta;

        /// Initialises PulseWaveSystem class members.
        PulseWaveSystem(const LibUtilities::SessionReaderSharedPtr& m_session);
		
		virtual void v_InitObject();
		
		/// Sets up initial conditions.
        virtual void v_DoInitialise();
	
		/// Solves an unsteady problem.
        virtual void v_DoSolve();
		
		/// Links the subdomains
        void LinkSubdomains(Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &fields);
		
		/// Riemann Problem for Bifurcation
		void BifurcationRiemann1_to_2(Array<OneD, NekDouble> &Au, Array<OneD, NekDouble> &uu,
									  Array<OneD, NekDouble> &beta, Array<OneD, NekDouble> &A_0);
		
		/// Riemann Problem for Junction
		void JunctionRiemann(Array<OneD, NekDouble> &Au, Array<OneD, NekDouble> &uu,
							 Array<OneD, NekDouble> &beta, Array<OneD, NekDouble> &A_0);
		
		/// Riemann Problem for Merging Flow
		void MergingRiemann2_to_1(Array<OneD, NekDouble> &Au, Array<OneD, NekDouble> &uu,
								  Array<OneD, NekDouble> &beta, Array<OneD, NekDouble> &A_0);
		
		/// Get the cross sectional area from the inputfile 
		void StaticArea(void);
		
		/// Get the material parameters from the inputfile 
		void MaterialProperties(void);
		
		// Ouptut field information
        virtual void v_Output(void);
		
		/// Prepares the multidomain output
        void PrepareMultidomainOutput(void);
		
		/// Compute the L2 error between fields and a given exact solution.
        NekDouble v_L2Error(unsigned int field,
							const Array<OneD,NekDouble> &exactsoln = NullNekDouble1DArray,
							bool Normalised = false);

        /// Compute the L_inf error between fields and a given exact solution.
        NekDouble v_LinfError(unsigned int field,
							  const Array<OneD,NekDouble> &exactsoln = NullNekDouble1DArray);
		
    private:
		
				
    };
}

#endif
