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
// Description: Base class for Navier-Stokes advection term
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADVECTIONTERM_H
#define NEKTAR_SOLVERS_ADVECTIONTERM_H

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/FFT/NektarFFT.h>  // for NektarFFTSharedPtr
#include <SpatialDomains/MeshGraph.h>   // for MeshGraphSharedPtr
#include <MultiRegions/ExpList.h>       // for ExpListSharedPtr


namespace Nektar
{
    class AdvectionTerm;

    /// A shared pointer to an EquationSystem object
    typedef boost::shared_ptr<AdvectionTerm> AdvectionTermSharedPtr;
    /// Datatype of the NekFactory used to instantiate classes derived from
    /// the EquationSystem class.
    typedef LibUtilities::NekFactory<
            std::string, AdvectionTerm,
            const LibUtilities::SessionReaderSharedPtr&,
            const SpatialDomains::MeshGraphSharedPtr&
        > AdvectionTermFactory;
    AdvectionTermFactory& GetAdvectionTermFactory();

    /// Base class for the development of solvers.
    class AdvectionTerm
    {
    public:
        /// Destructor
        virtual ~AdvectionTerm();
        
        inline void InitObject();
        
        /// Compute advection term
        void DoAdvection(Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
                         const int nConvectiveFields,
                         const Array<OneD, int>  &vel_loc,
                         const Array<OneD, const Array<OneD, NekDouble> > &pInarray,
                         Array<OneD, Array<OneD, NekDouble> > &pOutarray,
						 NekDouble m_time,
                         Array<OneD, NekDouble> &pWk = NullNekDouble1DArray);

        void DoAdvection(Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
                         const Array<OneD, const Array<OneD, NekDouble> > &Velocity,
                         const Array<OneD, const Array<OneD, NekDouble> > &pInarray,
                         Array<OneD, Array<OneD, NekDouble> > &pOutarray,
						 NekDouble m_time,
                         Array<OneD, NekDouble> &pWk = NullNekDouble1DArray);
	protected:
        LibUtilities::SessionReaderSharedPtr        m_session;
        /// Name of the session
        std::string m_sessionName;
        /// Pointer to mesh graph
        SpatialDomains::MeshGraphSharedPtr          m_graph;
		
        bool m_dealiasing;           ///< flag to determine if use dealising or not
        bool m_SingleMode;               ///< Flag to determine if use single mode or not
        bool m_HalfMode;                 ///< Flag to determine if use half mode or not
        
        MultiRegions::CoeffState m_CoeffState;

        /// Type of projection, i.e. Galerkin or DG.
        enum MultiRegions::ProjectionType m_projectionType;
        
        int m_spacedim;              ///< Spatial dimension (> expansion dim)
        int m_expdim;                ///< Dimension of the expansion
        int nvariables;              ///< Number of variables
        
        int m_nConvectiveFields;     /// Number of fields to be convected;
	
        //number of slices
        int                                             m_slices;
        //period length
        NekDouble										m_period;
        //interpolation vector
        Array<OneD, Array<OneD, NekDouble> >			m_interp;
        //auxiliary variables for time depedent base flows
        LibUtilities::NektarFFTSharedPtr				m_FFT;
        Array<OneD,NekDouble>							m_tmpIN;
        Array<OneD,NekDouble>							m_tmpOUT;
        bool											    m_useFFTW;
	
        /// Constructor
        AdvectionTerm(const LibUtilities::SessionReaderSharedPtr&        pSession,
                      const SpatialDomains::MeshGraphSharedPtr&          pGraph);
        
        virtual void v_InitObject();
        
        virtual void v_ComputeAdvectionTerm(Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
                                            const Array<OneD, Array<OneD, NekDouble> > &pV,
                                            const Array<OneD, const NekDouble> &pU,
                                            Array<OneD, NekDouble> &pOutarray,
                                            int pVelocityComponent,
                                            NekDouble m_time,
                                            Array<OneD, NekDouble> &pWk)
        {
            ASSERTL0(false,"This function is not defined in parent class");
        };
	
        int NoCaseStringCompare(const string & s1, const string& s2);
    };
    
    inline void AdvectionTerm::InitObject()
    {
        v_InitObject();
    }
} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H

