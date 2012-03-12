///////////////////////////////////////////////////////////////////////////////
//
// File AdjointAdvection.h
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

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/DisContField2D.h>

#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <MultiRegions/DisContField3DHomogeneous2D.h>

#include <IncNavierStokesSolver/AdvectionTerms/AdvectionTerm.h>

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
        friend class MemoryManager<AdjointAdvection>;

        /// Creates an instance of this class
        static AdvectionTermSharedPtr create(
                                const LibUtilities::SessionReaderSharedPtr& pSession,
                                const SpatialDomains::MeshGraphSharedPtr& pGraph) {
            AdvectionTermSharedPtr p = MemoryManager<AdjointAdvection>::AllocateSharedPtr(pSession, pGraph);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;

	protected:
        //Storage of the base flow
        Array<OneD, MultiRegions::ExpListSharedPtr>     m_base;

        AdjointAdvection(
                const LibUtilities::SessionReaderSharedPtr&        pSession,
                const SpatialDomains::MeshGraphSharedPtr&          pGraph);

        virtual ~AdjointAdvection();

        virtual void v_InitObject();

        void SetUpBaseFields(SpatialDomains::MeshGraphSharedPtr &mesh);
		void UpdateBase(const NekDouble m_slices,
						Array<OneD, const NekDouble> &inarray,
						Array<OneD, NekDouble> &outarray,
						const NekDouble m_time,
						const NekDouble m_period);
		

        /// Import Base flow
        void ImportFldBase(std::string pInfile,
                SpatialDomains::MeshGraphSharedPtr pGraph,int cnt);
		void ImportFldBase(std::string pInfile,
						   SpatialDomains::MeshGraphSharedPtr pGraph);
		
		/// Write field data to the given filename.
        void WriteFldBase(std::string &outname);
		
        /// Write input fields to the given filename.
        void WriteFldBase(
						  std::string &outname,
						  MultiRegions::ExpListSharedPtr &field,
						  Array<OneD, Array<OneD, NekDouble> > &fieldcoeffs,
						  Array<OneD, std::string> &variables);
		


    private:
        //Function for the evaluation of the Adjoint advective terms
        virtual void v_ComputeAdvectionTerm(
                         Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
                         const Array<OneD, Array<OneD, NekDouble> > &pV,
                         const Array<OneD, const NekDouble> &pU,
                         Array<OneD, NekDouble> &pOutarray,
                         int pVelocityComponent,
						 NekDouble m_time,
                         Array<OneD, NekDouble> &pWk);
    };
    
    
} //end of namespace

#endif //NEKTAR_SOLVERS_INCNAVIERSTOKES_H
