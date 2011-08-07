///////////////////////////////////////////////////////////////////////////////
//
// File AdjointAdvection.cpp
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
// Description: Evaluation of the adjoint advective term
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/AdvectionTerms/AdjointAdvection.h>
#include <cstdio>
#include <cstdlib>

namespace Nektar
{
    string AdjointAdvection::className = GetAdvectionTermFactory().RegisterCreatorFunction("Adjoint", AdjointAdvection::create);


    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */
    
    AdjointAdvection::AdjointAdvection(
            LibUtilities::SessionReaderSharedPtr        pSession,
            SpatialDomains::MeshGraphSharedPtr          pGraph):
        LinearisedAdvection(pSession, pGraph)
    {
        SetUpBaseFields(pGraph);
        ImportFldBase(pSession->GetFilename().substr(0,pSession->GetFilename().find_last_of('.')) + ".bse",pGraph);
    }
    
    
    AdjointAdvection::~AdjointAdvection()
    {
    }
    

    //Evaluation of the advective terms
    void AdjointAdvection::v_ComputeAdvectionTerm(
            Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const Array<OneD, Array<OneD, NekDouble> > &pVelocity,
            const Array<OneD, const NekDouble> &pU,
            Array<OneD, NekDouble> &pOutarray,
            int pVelocityComponent,
            Array<OneD, NekDouble> &pWk)
    {
        int ndim       = m_nConvectiveFields;
        int nPointsTot = pFields[0]->GetNpoints();
        Array<OneD, NekDouble> grad0,grad1,grad2;
	
        //Evaluation of the gradiend of each component of the base flow
        //\nabla U
        Array<OneD, NekDouble> grad_base_u0,grad_base_u1,grad_base_u2;
        // \nabla V
        Array<OneD, NekDouble> grad_base_v0,grad_base_v1,grad_base_v2;
        // \nabla W
        Array<OneD, NekDouble> grad_base_w0,grad_base_w1,grad_base_w2;
	
        
        grad0 = Array<OneD, NekDouble> (ndim*nPointsTot, 0.0);
        grad_base_u0 = Array<OneD, NekDouble> (ndim*nPointsTot, 0.0);
        grad_base_v0 = Array<OneD, NekDouble> (ndim*nPointsTot, 0.0);
        grad_base_w0 = Array<OneD, NekDouble> (ndim*nPointsTot, 0.0);
	
        //Evaluate the adjoint advection term
        switch(ndim) 
        {
            // 1D
        case 1:
            pFields[0]->PhysDeriv(pVelocity[pVelocityComponent],grad0);
            pFields[0]->PhysDeriv(m_base[0]->GetPhys(),grad_base_u0);
            //Evaluate  U du'/dx
            Vmath::Vmul(nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
            //Evaluate U du'/dx+ u' dU/dx
            Vmath::Vvtvp(nPointsTot,grad_base_u0,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
            break;
            
            //2D
        case 2:
            
            grad1 = grad0 + nPointsTot;
            grad_base_u1 = grad_base_u0 + nPointsTot;
            grad_base_v1 = grad_base_v0 +nPointsTot;
            
            pFields[0]->PhysDeriv(pVelocity[pVelocityComponent],grad0,grad1);
            
            //Derivates of the base flow
            pFields[0]-> PhysDeriv(m_base[0]->GetPhys(), grad_base_u0, grad_base_u1);
            pFields[0]-> PhysDeriv(m_base[1]->GetPhys(), grad_base_v0, grad_base_v1);
            
            //Since the components of the velocity are passed one by one, it is necessary to distinguish which
            //term is consider
            switch (pVelocityComponent)
            {
                //x-equation
            case 0:
                // Evaluate U du'/dx
                Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
                //Evaluate U du'/dx+ V du'/dy
                Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
                //Evaluate - (U du'/dx+ V du'/dy)
                Vmath::Neg(nPointsTot,pOutarray,1);
                //Evaluate -(U du'/dx+ V du'/dy)+u' dU/dx
                Vmath::Vvtvp(nPointsTot,grad_base_u0,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
                        //Evaluate -(U du'/dx+ V du'/dy) +u' dU/dx +v' dV/dx
                Vmath::Vvtvp(nPointsTot,grad_base_v0,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
                break;
		
                //y-equation
            case 1:
                // Evaluate U dv'/dx
                Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
                //Evaluate U dv'/dx+ V dv'/dy
                Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
                //Evaluate -(U dv'/dx+ V dv'/dy)
                Vmath::Neg(nPointsTot,pOutarray,1);
                //Evaluate (U dv'/dx+ V dv'/dy)+u' dU/dy
                Vmath::Vvtvp(nPointsTot,grad_base_u1,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
                //Evaluate (U dv'/dx+ V dv'/dy +u' dv/dx)+v' dV/dy
                Vmath::Vvtvp(nPointsTot,grad_base_v1,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
                break;
            }
            break;
            
            
            //3D
        case 3:
            
            grad1 = grad0 + nPointsTot;
            grad2 = grad1 + nPointsTot;
            grad_base_u1 = grad_base_u0 + nPointsTot;
            grad_base_v1 = grad_base_v0 +nPointsTot;
            grad_base_u2 = grad_base_u1 +nPointsTot;
            grad_base_v2 = grad_base_v1 +nPointsTot;
            
            pFields[0]->PhysDeriv(pVelocity[pVelocityComponent], grad0, grad1, grad2);
            
            pFields[0]->PhysDeriv(m_base[0]->GetPhys(), grad_base_u0, grad_base_u1,grad_base_u2);
            pFields[0]->PhysDeriv(m_base[1]->GetPhys(), grad_base_v0, grad_base_v1,grad_base_v2);
            pFields[0]->PhysDeriv(m_base[2]->GetPhys(), grad_base_w0, grad_base_w1, grad_base_w2);
            
            switch (pVelocityComponent)
            {
					//x-equation	
				case 0:
					//Evaluate U du'/dx
					Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
					//Evaluate U du'/dx+ V du'/dy
					Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
					//Evaluate U du'/dx+ V du'/dy+W du'/dz
					Vmath::Vvtvp(nPointsTot,grad2,1,m_base[2]->GetPhys(),1,pOutarray,1,pOutarray,1);
					//Evaluate -(U du'/dx+ V du'/dy+W du'/dz)
					Vmath::Neg(nPointsTot,pOutarray,1);
					//Evaluate -(U du'/dx+ V du'/dy+W du'/dz)+u' dU/dx
					Vmath::Vvtvp(nPointsTot,grad_base_u0,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
					//Evaluate -(U du'/dx+ V du'/dy+W du'/dz)+u'dU/dx+ v' dV/dx
					Vmath::Vvtvp(nPointsTot,grad_base_v0,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
					//Evaluate -(U du'/dx+ V du'/dy+W du'/dz)+u'dU/dx+ v' dV/dx+ w' dW/dz
					Vmath::Vvtvp(nPointsTot,grad_base_w0,1,pVelocity[2],1,pOutarray,1,pOutarray,1);
					break;
					//y-equation	
				case 1:
					//Evaluate U dv'/dx
					Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
					//Evaluate U dv'/dx+ V dv'/dy
					Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
					//Evaluate U dv'/dx+ V dv'/dy+W dv'/dz
					Vmath::Vvtvp(nPointsTot,grad2,1,m_base[2]->GetPhys(),1,pOutarray,1,pOutarray,1);
					//Evaluate -(U dv'/dx+ V dv'/dy+W dv'/dz)
					Vmath::Neg(nPointsTot,pOutarray,1);
					//Evaluate  -(U dv'/dx+ V dv'/dy+W dv'/dz)+u' dU/dy
					Vmath::Vvtvp(nPointsTot,grad_base_u1,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
					//Evaluate  -(U dv'/dx+ V dv'/dy+W dv'/dz)+u' dU/dy +v' dV/dy
					Vmath::Vvtvp(nPointsTot,grad_base_v1,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
					//Evaluate  -(U dv'/dx+ V dv'/dy+W dv'/dz)+u' dU/dy +v' dV/dy+ w' dW/dy
					Vmath::Vvtvp(nPointsTot,grad_base_w1,1,pVelocity[2],1,pOutarray,1,pOutarray,1);
					break;
					
					//z-equation	
				case 2:
					//Evaluate U dw'/dx
					Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
					//Evaluate U dw'/dx+ V dw'/dx
					Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
					//Evaluate U dw'/dx+ V dw'/dx+ W dw'/dz
					Vmath::Vvtvp(nPointsTot,grad2,1,m_base[2]->GetPhys(),1,pOutarray,1,pOutarray,1);
					//Evaluate -(U dw'/dx+ V dw'/dx+ W dw'/dz)
					Vmath::Neg(nPointsTot,pOutarray,1);
					//Evaluate -(U dw'/dx+ V dw'/dx+ W dw'/dz)+u' dU/dz
					Vmath::Vvtvp(nPointsTot,grad_base_u2,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
					//Evaluate -(U dw'/dx+ V dw'/dx+ W dw'/dz)+u' dU/dz+v'dV/dz
					Vmath::Vvtvp(nPointsTot,grad_base_v2,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
					//Evaluate -(U dw'/dx+ V dw'/dx+ W dw'/dz)+u' dU/dz+v'dV/dz + w' dW/dz
					Vmath::Vvtvp(nPointsTot,grad_base_w2,1,pVelocity[2],1,pOutarray,1,pOutarray,1);
					break;
            }
            break;
            
        default:
            ASSERTL0(false,"dimension unknown");
        }
    }	
	
} //end of namespace

