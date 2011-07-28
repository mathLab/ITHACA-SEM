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
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/HistoryPoints.h>
#include <SpatialDomains/SpatialData.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/DisContField2D.h>

#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <MultiRegions/DisContField3DHomogeneous2D.h>

#include <Auxiliary/EquationSystem.h>

namespace Nektar
{
    class AdvectionTerm;

    /// A shared pointer to an EquationSystem object
    typedef boost::shared_ptr<AdvectionTerm> AdvectionTermSharedPtr;
    /// Datatype of the NekFactory used to instantiate classes derived from
    /// the EquationSystem class.
    typedef LibUtilities::NekFactory<
            std::string, AdvectionTerm,
            LibUtilities::CommSharedPtr&,
            LibUtilities::SessionReaderSharedPtr&,
            SpatialDomains::MeshGraphSharedPtr&,
            SpatialDomains::BoundaryConditionsSharedPtr&
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
		inline void DoAdvection(
                               Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
                               const Array<OneD, const Array<OneD, NekDouble> > &pInarray,
                               Array<OneD, Array<OneD, NekDouble> > &pOutarray,
                               Array<OneD, NekDouble> &pWk = NullNekDouble1DArray);
		
	protected:
		LibUtilities::CommSharedPtr                 m_comm;
		/// Filename of session
		LibUtilities::SessionReaderSharedPtr        m_session;
		/// Name of the session
        std::string m_sessionName;
        /// Pointer to boundary conditions object.
        SpatialDomains::BoundaryConditionsSharedPtr m_boundaryConditions;
        /// Pointer to mesh graph
        SpatialDomains::MeshGraphSharedPtr          m_graph;

        /// Type of projection, i.e. Galerkin or DG.
        enum MultiRegions::ProjectionType m_projectionType;

        int m_spacedim;              ///< Spatial dimension (> expansion dim)
        int m_expdim;                ///< Dimension of the expansion
		int nvariables;              ///< Number of variables

        /// Constructor
        AdvectionTerm(
                LibUtilities::CommSharedPtr&                 pComm,
                LibUtilities::SessionReaderSharedPtr&        pSession,
                SpatialDomains::MeshGraphSharedPtr&          pGraph,
                SpatialDomains::BoundaryConditionsSharedPtr& pBoundaryConditions);

        virtual void v_InitObject();

        //Virtual function for the evaluation of the advective term
        virtual void v_DoAdvection(
                                   Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
                                   const Array<OneD, const Array<OneD, NekDouble> > &pInarray,
                                   Array<OneD, Array<OneD, NekDouble> > &pOutarray,
                                   Array<OneD, NekDouble> &pWk = NullNekDouble1DArray);
		
		
		int NoCaseStringCompare(const string & s1, const string& s2) ;

	};
    
    inline void AdvectionTerm::InitObject()
    {
        v_InitObject();
    }

    inline void AdvectionTerm::DoAdvection(
                           Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
                           const Array<OneD, const Array<OneD, NekDouble> > &pInarray,
                           Array<OneD, Array<OneD, NekDouble> > &pOutarray,
                           Array<OneD, NekDouble> &pWk)
    {
        v_DoAdvection(pFields, pInarray, pOutarray, pWk);
    }
	
} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H

