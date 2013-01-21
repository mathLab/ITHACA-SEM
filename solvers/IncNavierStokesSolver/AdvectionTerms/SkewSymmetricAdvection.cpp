///////////////////////////////////////////////////////////////////////////////
//
// File SkewSymmetricAdvection.cpp
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

#include <IncNavierStokesSolver/AdvectionTerms/SkewSymmetricAdvection.h>

namespace Nektar
{
    string SkewSymmetricAdvection::className  = GetAdvectionTermFactory().RegisterCreatorFunction("SkewSymmetric", SkewSymmetricAdvection::create);
    
    
    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */

    SkewSymmetricAdvection::SkewSymmetricAdvection(
            const LibUtilities::SessionReaderSharedPtr&        pSession,
            const SpatialDomains::MeshGraphSharedPtr&          pGraph):
        AdvectionTerm(pSession, pGraph)
	
    {
        
    }
    
    SkewSymmetricAdvection::~SkewSymmetricAdvection()
    {
    }
    
    //Advection function
    
    
    //Evaluation of the advective terms
    void SkewSymmetricAdvection::v_ComputeAdvectionTerm(
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
        Array<OneD, NekDouble> gradV0,gradV1,gradV2, tmp, Up;
		
        gradV0   = Array<OneD, NekDouble> (nPointsTot);
		tmp = Array<OneD, NekDouble> (nPointsTot);
		
        // Evaluate V\cdot Grad(u)
        switch(ndim)
        {
        case 1:
			pFields[0]->PhysDeriv(pU,gradV0);
			Vmath::Vmul(nPointsTot,gradV0,1,pV[0],1,pOutarray,1);
			Vmath::Vmul(nPointsTot,pU,1,pV[0],1,gradV0,1);
			pFields[0]->PhysDeriv(gradV0,tmp);
			Vmath::Vadd(nPointsTot,tmp,1,pOutarray,1,pOutarray,1);
			Vmath::Smul(nPointsTot,0.5,pOutarray,1,pOutarray,1);
            break;
        case 2:
			gradV1 = Array<OneD, NekDouble> (nPointsTot);
			pFields[0]->PhysDeriv(pU,gradV0,gradV1);
			Vmath::Vmul (nPointsTot,gradV0,1,pV[0],1,pOutarray,1);
			Vmath::Vvtvp(nPointsTot,gradV1,1,pV[1],1,pOutarray,1,pOutarray,1);
			Vmath::Vmul(nPointsTot,pU,1,pV[0],1,gradV0,1);
			Vmath::Vmul(nPointsTot,pU,1,pV[1],1,gradV1,1);
			pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],gradV0,tmp);
			Vmath::Vadd(nPointsTot,tmp,1,pOutarray,1,pOutarray,1);
			pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],gradV1,tmp);
			Vmath::Vadd(nPointsTot,tmp,1,pOutarray,1,pOutarray,1);
			Vmath::Smul(nPointsTot,0.5,pOutarray,1,pOutarray,1);
            break;	 
        case 3:
			gradV1 = Array<OneD, NekDouble> (nPointsTot);
			gradV2 = Array<OneD, NekDouble> (nPointsTot);
			
			pFields[0]->PhysDeriv(pU,gradV0,gradV1,gradV2);
			
			//pOutarray = 1/2(u*du/dx + v*du/dy + w*du/dz + duu/dx + duv/dy + duw/dz)
				
			if(m_dealiasing == true && pFields[0]->GetWaveSpace() == false) 
			{
				pFields[0]->DealiasedProd(pV[0],gradV0,gradV0,m_CoeffState);
				pFields[0]->DealiasedProd(pV[1],gradV1,gradV1,m_CoeffState);
				pFields[0]->DealiasedProd(pV[2],gradV2,gradV2,m_CoeffState);
				Vmath::Vadd(nPointsTot,gradV0,1,gradV1,1,pOutarray,1);
				Vmath::Vadd(nPointsTot,gradV2,1,pOutarray,1,pOutarray,1);
				pFields[0]->DealiasedProd(pU,pV[0],gradV0,m_CoeffState);
				pFields[0]->DealiasedProd(pU,pV[1],gradV1,m_CoeffState);
				pFields[0]->DealiasedProd(pU,pV[2],gradV2,m_CoeffState);
				pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],gradV0,tmp);
				Vmath::Vadd(nPointsTot,tmp,1,pOutarray,1,pOutarray,1);
				pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],gradV1,tmp);
				Vmath::Vadd(nPointsTot,tmp,1,pOutarray,1,pOutarray,1);
				pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],gradV2,tmp);
				Vmath::Vadd(nPointsTot,tmp,1,pOutarray,1,pOutarray,1);
				Vmath::Smul(nPointsTot,0.5,pOutarray,1,pOutarray,1);
			}
			else if(pFields[0]->GetWaveSpace() == true && m_dealiasing == false)
			{
				Up = Array<OneD, NekDouble> (nPointsTot);
				//vector reused to avoid even more memory requirements
				//names may be misleading
				pFields[0]->HomogeneousBwdTrans(gradV0,tmp);
				Vmath::Vmul(nPointsTot,tmp,1,pV[0],1,pOutarray,1); // + u*du/dx
				pFields[0]->HomogeneousBwdTrans(gradV1,tmp);
				Vmath::Vvtvp(nPointsTot,tmp,1,pV[1],1,pOutarray,1,pOutarray,1);// + v*du/dy
				pFields[0]->HomogeneousBwdTrans(gradV2,tmp);
				Vmath::Vvtvp(nPointsTot,tmp,1,pV[2],1,pOutarray,1,pOutarray,1);// + w*du/dz
				
				pFields[0]->HomogeneousBwdTrans(pU,Up);
				Vmath::Vmul(nPointsTot,Up,1,pV[0],1,gradV0,1);
				Vmath::Vmul(nPointsTot,Up,1,pV[1],1,gradV1,1);
				Vmath::Vmul(nPointsTot,Up,1,pV[2],1,gradV2,1);
				
				pFields[0]->SetWaveSpace(false);
				pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],gradV0,tmp);//duu/dx
				Vmath::Vadd(nPointsTot,tmp,1,pOutarray,1,pOutarray,1);
				pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],gradV1,tmp);//duv/dy
				Vmath::Vadd(nPointsTot,tmp,1,pOutarray,1,pOutarray,1);
				pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],gradV2,tmp);//duw/dz
				Vmath::Vadd(nPointsTot,tmp,1,pOutarray,1,pOutarray,1);
				pFields[0]->SetWaveSpace(true);
				
				Vmath::Smul(nPointsTot,0.5,pOutarray,1,tmp,1);
				pFields[0]->HomogeneousFwdTrans(tmp,pOutarray);
			}
			else if(pFields[0]->GetWaveSpace() == false && m_dealiasing == false) 
			{
				Vmath::Vmul(nPointsTot,gradV0,1,pV[0],1,pOutarray,1);
				Vmath::Vvtvp(nPointsTot,gradV1,1,pV[1],1,pOutarray,1,pOutarray,1);
				Vmath::Vvtvp(nPointsTot,gradV2,1,pV[2],1,pOutarray,1,pOutarray,1);
				Vmath::Vmul(nPointsTot,pU,1,pV[0],1,gradV0,1);
				Vmath::Vmul(nPointsTot,pU,1,pV[1],1,gradV1,1);
				Vmath::Vmul(nPointsTot,pU,1,pV[2],1,gradV2,1);
				pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],gradV0,tmp);
				Vmath::Vadd(nPointsTot,tmp,1,pOutarray,1,pOutarray,1);
				pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],gradV1,tmp);
				Vmath::Vadd(nPointsTot,tmp,1,pOutarray,1,pOutarray,1);
				pFields[0]->PhysDeriv(MultiRegions::DirCartesianMap[2],gradV2,tmp);
				Vmath::Vadd(nPointsTot,tmp,1,pOutarray,1,pOutarray,1);
				Vmath::Smul(nPointsTot,0.5,pOutarray,1,pOutarray,1);
			}
			else 
			{
                ASSERTL0(false, "Dealiasing is not allowed in combination "
                                "with the Skew-Symmetric advection form for "
                                "efficiency reasons.");	
			}
            break;
        default:
            ASSERTL0(false,"dimension unknown");
        }
    }

} //end of namespace

