///////////////////////////////////////////////////////////////////////////////
//
// File NoAdvection.h
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
// Description: NoAdvection
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_NOADVECTION_H
#define NEKTAR_SOLVERS_NOADVECTION_H

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

namespace Nektar
{     
    class NoAdvection: public AdvectionTerm
    {
    public:
        friend class MemoryManager<NoAdvection>;

        /// Creates an instance of this class
        static AdvectionTermSharedPtr create(
                                const LibUtilities::SessionReaderSharedPtr& pSession,
                                const SpatialDomains::MeshGraphSharedPtr& pGraph) {
            AdvectionTermSharedPtr p = MemoryManager<NoAdvection>::AllocateSharedPtr(pSession, pGraph);
            p->InitObject();
            return p;
        }
        /// Name of class
        static std::string className;
        
	protected:
        
        NoAdvection(const LibUtilities::SessionReaderSharedPtr&        pSession,
					const SpatialDomains::MeshGraphSharedPtr&          pGraph);


        virtual ~NoAdvection();

	private:

        //Function for the evaluation of the linearised advective terms
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
