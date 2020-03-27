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
    
    struct InterfacePoint
    {
        InterfacePoint(const int vid, 
                       const int domain, 
                       const int elmt, 
                       const int elmtVert, 
                       const int traceId,
                       const int bcpos):
        m_vid(vid),
            m_domain(domain),
            m_elmt(elmt),
            m_elmtVert(elmtVert),
            m_traceId(traceId),
            m_bcPosition(bcpos)
        {
        };
        int m_vid;        // Global Vid of interface point 
        int m_domain;     // domain interface point belongs to
        int m_elmt;       // element id of vertex 
        int m_elmtVert;   // vertex id within local element 
        int m_traceId;    // Element id  within the trace 
        int m_bcPosition; // Position of boundary condition in region
    };

    typedef std::shared_ptr<InterfacePoint> InterfacePointShPtr;

    /// Base class for unsteady solvers.
    class PulseWaveSystem : public UnsteadySystem
    {
    public:
        /// Destructor
        virtual ~PulseWaveSystem();

        int  GetNdomains()
        {
            return m_nDomains;
        }

        Array<OneD, MultiRegions::ExpListSharedPtr> UpdateVessels(void)
        {
            return m_vessels;
        }

        void CalcCharacteristicVariables(int omega);

    protected:
        Array<OneD, MultiRegions::ExpListSharedPtr>     m_vessels;
        int				                m_nDomains; 
        int                                             m_currentDomain;
        int                                             m_nVariables;
        UpwindTypePulse                                 m_upwindTypePulse;
		
        Array<OneD, int>                                m_fieldPhysOffset;
        NekDouble                                       m_rho;
        NekDouble                                       m_pext;
	
        NekDouble m_C;
        NekDouble m_RT;
        NekDouble m_pout;
	
        Array<OneD, Array<OneD, NekDouble> >		m_A_0;
        Array<OneD, Array<OneD, NekDouble> >		m_A_0_trace;
        Array<OneD, Array<OneD, NekDouble> >		m_beta;
        Array<OneD, Array<OneD, NekDouble> >		m_beta_trace;
        Array<OneD, Array<OneD, NekDouble> >		m_trace_fwd_normal;


        std::vector<std::vector<InterfacePointShPtr> >  m_vesselJcts;
        std::vector<std::vector<InterfacePointShPtr> >  m_bifurcations;
        std::vector<std::vector<InterfacePointShPtr> >  m_mergingJcts;
    
        /// Initialises PulseWaveSystem class members.
        PulseWaveSystem(const LibUtilities::SessionReaderSharedPtr& pSession,
                        const SpatialDomains::MeshGraphSharedPtr& pGraph);
        
        virtual void v_InitObject();
	
        /// Sets up initial conditions.
        virtual void v_DoInitialise();
	
        /// Solves an unsteady problem.
        virtual void v_DoSolve();
	
        /// Links the subdomains
        void LinkSubdomains(Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  &fields);
	
        /// Riemann Problem for Bifurcation
        void BifurcationRiemann(Array<OneD, NekDouble> &Au, Array<OneD, NekDouble> &uu,
                                Array<OneD, NekDouble> &beta, Array<OneD, NekDouble> &A_0);
        
        /// Riemann Problem for Merging Flow
        void MergingRiemann(Array<OneD, NekDouble> &Au, Array<OneD, NekDouble> &uu,
                            Array<OneD, NekDouble> &beta, Array<OneD, NekDouble> &A_0);
        
        /// Riemann Problem for Junction
        void JunctionRiemann(Array<OneD, NekDouble> &Au, Array<OneD, NekDouble> &uu,
                             Array<OneD, NekDouble> &beta, Array<OneD, NekDouble> &A_0);
        
        // Ouptut field information
        virtual void v_Output(void);

        // Checkpoint field output
        void CheckPoint_Output(const int n);
	
        /// Compute the L2 error between fields and a given exact solution.
        NekDouble v_L2Error(unsigned int field,
                            const Array<OneD,NekDouble> &exactsoln = NullNekDouble1DArray,
                            bool Normalised = false);
        
        /// Compute the L_inf error between fields and a given exact solution.
        NekDouble v_LinfError(unsigned int field,
                              const Array<OneD,NekDouble> &exactsoln = NullNekDouble1DArray);
        
        /// Write input fields to the given filename.
        void WriteVessels(const std::string &outname);
        
        void EnforceInterfaceConditions(const Array<OneD, const Array<OneD, NekDouble> > &fields);
        
    private:
        void SetUpDomainInterfaces(void);
        void FillDataFromInterfacePoint(InterfacePointShPtr &I, 
                         const Array<OneD, const Array<OneD, NekDouble> >&field, 
                                        NekDouble &A, NekDouble &u,
                                        NekDouble &beta, NekDouble &A_0);

            
    };
        
    typedef std::shared_ptr<PulseWaveSystem> PulseWaveSystemSharedPtr;

       
}

#endif
