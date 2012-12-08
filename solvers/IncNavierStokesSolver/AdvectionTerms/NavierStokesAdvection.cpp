///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesAdvection.cpp
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

#include <IncNavierStokesSolver/AdvectionTerms/NavierStokesAdvection.h>

namespace Nektar
{
    string NavierStokesAdvection::className  = GetAdvectionTermFactory().RegisterCreatorFunction("Convective", NavierStokesAdvection::create);
    string NavierStokesAdvection::className2 = GetAdvectionTermFactory().RegisterCreatorFunction("NonConservative", NavierStokesAdvection::create);
    
    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */

    NavierStokesAdvection::NavierStokesAdvection(
            const LibUtilities::SessionReaderSharedPtr&        pSession,
            const SpatialDomains::MeshGraphSharedPtr&          pGraph):
        AdvectionTerm(pSession, pGraph)
	
    {
        
    }
    
    NavierStokesAdvection::~NavierStokesAdvection()
    {
    }
    
    //Advection function
    
    
    //Evaluation of the advective terms
    void NavierStokesAdvection::v_ComputeAdvectionTerm(
            Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const Array<OneD, Array<OneD, NekDouble> > &pV,
            const Array<OneD, const NekDouble> &pU,
            Array<OneD, NekDouble> &pOutarray,
            int pVelocityComponent,
            NekDouble m_time,
            Array<OneD, NekDouble> &pWk)
    {
        // use dimension of Velocity vector to dictate dimension of operation
        int ndim       = pV.num_elements();
        
        // ToDo: here we should add a check that V has right dimension
	
        int nPointsTot = pFields[0]->GetNpoints();
        Array<OneD, NekDouble> grad0,grad1,grad2;
		
        grad0 = Array<OneD, NekDouble> (nPointsTot);
		
        // Evaluate V\cdot Grad(u)
        switch(ndim)
        {
        case 1:
            pFields[0]->PhysDeriv(pU,grad0);
            Vmath::Vmul(nPointsTot,grad0,1,pV[0],1,pOutarray,1);
            break;
        case 2:
            grad1 = Array<OneD, NekDouble> (nPointsTot);
            pFields[0]->PhysDeriv(pU,grad0,grad1);
            Vmath::Vmul (nPointsTot,grad0,1,pV[0],1,pOutarray,1);
            Vmath::Vvtvp(nPointsTot,grad1,1,pV[1],1,pOutarray,1,pOutarray,1);
            break;	 
        case 3:
            grad1 = Array<OneD, NekDouble> (nPointsTot);
            grad2 = Array<OneD, NekDouble> (nPointsTot);
            pFields[0]->PhysDeriv(pU,grad0,grad1,grad2);
            
            if(m_dealiasing == true && pFields[0]->GetWaveSpace() == false)
            {
                pFields[0]->DealiasedProd(pV[0],grad0,grad0,m_CoeffState);
                pFields[0]->DealiasedProd(pV[1],grad1,grad1,m_CoeffState);
                pFields[0]->DealiasedProd(pV[2],grad2,grad2,m_CoeffState);
                Vmath::Vadd(nPointsTot,grad0,1,grad1,1,pOutarray,1);
                Vmath::Vadd(nPointsTot,grad2,1,pOutarray,1,pOutarray,1);
            }
            else if(pFields[0]->GetWaveSpace() == true && m_dealiasing == false)
            {
                //vector reused to avoid even more memory requirements
				//names may be misleading
                pFields[0]->HomogeneousBwdTrans(grad0,pOutarray);
                Vmath::Vmul(nPointsTot,pOutarray,1,pV[0],1,pOutarray,1);
		
                pFields[0]->HomogeneousBwdTrans(grad1,grad0);
                Vmath::Vvtvp(nPointsTot,grad0,1,pV[1],1,pOutarray,1,pOutarray,1);
		
                pFields[0]->HomogeneousBwdTrans(grad2,grad1);
                Vmath::Vvtvp(nPointsTot,grad1,1,pV[2],1,pOutarray,1,grad0,1);
		
                pFields[0]->HomogeneousFwdTrans(grad0,pOutarray);
            }
            else if(pFields[0]->GetWaveSpace() == false && m_dealiasing == false) 
            {
                Vmath::Vmul(nPointsTot,grad0,1,pV[0],1,pOutarray,1);
                Vmath::Vvtvp(nPointsTot,grad1,1,pV[1],1,pOutarray,1,pOutarray,1);
                Vmath::Vvtvp(nPointsTot,grad2,1,pV[2],1,pOutarray,1,pOutarray,1);
            }
            else if(pFields[0]->GetWaveSpace() == true && m_dealiasing == true) 
            {
                pFields[0]->HomogeneousBwdTrans(grad0, pOutarray);
                pFields[0]->DealiasedProd(pV[0], pOutarray, grad0, 
                                          m_CoeffState);

                pFields[0]->HomogeneousBwdTrans(grad1,pOutarray);
                pFields[0]->DealiasedProd(pV[1], pOutarray, grad1,
                                          m_CoeffState);

                pFields[0]->HomogeneousBwdTrans(grad2,pOutarray);
                pFields[0]->DealiasedProd(pV[2], pOutarray, grad2,
                                          m_CoeffState);

                Vmath::Vadd(nPointsTot, grad0, 1, grad1, 1, grad0, 1);
                Vmath::Vadd(nPointsTot, grad0, 1, grad2, 1, grad0, 1);

                pFields[0]->HomogeneousFwdTrans(grad0,pOutarray);
            }
            else 
            {
                ASSERTL0(false, "Advection term calculation not implented or "
                                "possible with the current problem set up");
            }
            break;
        default:
            ASSERTL0(false,"dimension unknown");
        }
    }

} //end of namespace

