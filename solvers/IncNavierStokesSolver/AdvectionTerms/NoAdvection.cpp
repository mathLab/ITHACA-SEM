///////////////////////////////////////////////////////////////////////////////
//
// File NoAdvection.cpp
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
// Description: Evaluation of the Navier Stokes advective term
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/AdvectionTerms/NoAdvection.h>
#include <cstdio>
#include <cstdlib>

namespace Nektar
{
    string NoAdvection::className  = GetAdvectionTermFactory().RegisterCreatorFunction("NoAdvection", NoAdvection::create);
    
    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */

    NoAdvection::NoAdvection(const LibUtilities::SessionReaderSharedPtr&        pSession,
							 const SpatialDomains::MeshGraphSharedPtr&          pGraph):
	AdvectionTerm(pSession, pGraph)
	
    {
        
    }
    
    NoAdvection::~NoAdvection()
    {
    }
    
    //Advection function
    
    
    //Evaluation of the advective terms
    void NoAdvection::v_ComputeAdvectionTerm(
            Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const Array<OneD, Array<OneD, NekDouble> > &pV,
            const Array<OneD, const NekDouble> &pU,
            Array<OneD, NekDouble> &pOutarray,
            int pVelocityComponent,
            NekDouble m_time,
            Array<OneD, NekDouble> &pWk)
    {
        int nPointsTot = pFields[0]->GetNpoints();
		Vmath::Zero(nPointsTot,pOutarray,1);
	}

} //end of namespace

