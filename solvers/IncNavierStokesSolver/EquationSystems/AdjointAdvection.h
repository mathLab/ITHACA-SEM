///////////////////////////////////////////////////////////////////////////////
//
// File LinearisedAdvection.h
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
// Description: TBA
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_ADJOINTADVECTION_H
#define NEKTAR_SOLVERS_ADJOINTADVECTION_H

#include <SpatialDomains/MeshComponents.h>
#include <SpatialDomains/HistoryPoints.h>
#include <SpatialDomains/SpatialData.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/DisContField2D.h>

#include <IncNavierStokesSolver/EquationSystems/AdvectionTerm.h>

//#define TIMING

//#ifdef TIMING
//#include <time.h>
//#include <sys/time.h>
//#endif


namespace Nektar
{     


    class AdjointAdvection: public AdvectionTerm
    {
    public:           

        /**
         * Default constructor. 
         * 
         */ 
        AdjointAdvection();


        /**
         * Constructor.
         * \param 
         * 
         */

        AdjointAdvection(
                LibUtilities::CommSharedPtr                 pComm,
                LibUtilities::SessionReaderSharedPtr        pSession,
                SpatialDomains::MeshGraphSharedPtr          pGraph,
                SpatialDomains::BoundaryConditionsSharedPtr pBoundaryConditions);

		virtual ~AdjointAdvection();
     
		//Virtual function for the evaluation of the advective terms
		virtual void v_DoAdvection(
								   Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
								   const Array<OneD, const Array<OneD, NekDouble> > &pInarray,
								   Array<OneD, Array<OneD, NekDouble> > &pOutarray,
								   Array<OneD, NekDouble> &pWk);
		
	protected:
		//Storage of the base flow
		Array<OneD, MultiRegions::ExpListSharedPtr>     m_base;
 		int                                             m_nConvectiveFields;
		Array<OneD, int>                                m_velocity;

        void SetUpBaseFields(SpatialDomains::MeshGraphSharedPtr &mesh);

        /// Import Base flow
        void ImportFldBase(std::string pInfile,
                SpatialDomains::MeshGraphSharedPtr pGraph,
                SpatialDomains::BoundaryConditionsSharedPtr &pBoundaryConditions);

        //Function for the evaluation of the adjoint advective terms
        void ComputeAdvectionTerm(
                         Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
                         int pVelocityComponent,
                         const Array<OneD, Array<OneD, NekDouble> > &pVelocity,
                         Array<OneD, NekDouble> &pOutarray);

    };
    
    
} //end of namespace

#endif //NEKTAR_SOLVERS_INCNAVIERSTOKES_H
