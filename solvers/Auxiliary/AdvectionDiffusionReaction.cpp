///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionDiffusionReaction.cpp
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
// Description: Basic Advection Diffusion Reaction definition for 
// multiple fields in a 2D domain with boundary conditions
//
///////////////////////////////////////////////////////////////////////////////

#include <../solvers/Auxiliary/AdvectionDiffusionReaction.h>
#include <cstdio>
#include <cstdlib>

namespace Nektar
{
    /**
     * Basic construnctor
     */
    AdvectionDiffusionReaction::AdvectionDiffusionReaction(void):
        m_nvariables(0),
        m_fields(0)
    {     
    }
    
    /**
     * Constructor. Creates ... of #DisContField2D fields
     *
     * \param 
     * \param
     */
    AdvectionDiffusionReaction::AdvectionDiffusionReaction(SpatialDomains::MeshGraph2D &graph2D,
                                                           SpatialDomains::BoundaryConditions &bcs,
                                                           int variables):
      m_fields(variables)
    {
        m_nvariables = variables;
           
        for(int i = 0 ; i < m_fields.num_elements(); i++)
	{
            m_fields[i] = MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(graph2D,bcs,i);
	}
    }
    
    
    void AdvectionDiffusionReaction::SetInitialConditions(SpatialDomains::BoundaryConditions &bcs, int initialtime)
    {
        int nq = m_fields[0]->GetPointsTot();
      
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
      
        // get the coordinates of the quad points
        // (assuming all fields have the same discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);
      
        for(int i = 0 ; i < m_fields.num_elements(); i++)
	{
            SpatialDomains::ConstInitialConditionShPtr ifunc = bcs.GetInitialCondition(i);
            for(int j = 0; j < nq; j++)
	    {
                (m_fields[i]->UpdatePhys())[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],initialtime);
	    }
            m_fields[i]->SetPhysState(true);
	}
    }
    
    void AdvectionDiffusionReaction::GetExactSolutions(SpatialDomains::BoundaryConditions &bcs)
    {
        int nq = m_fields[0]->GetPointsTot();
      
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
      
        // get the coordinates of the quad points
        m_fields[0]->GetCoords(x0,x1,x2);
      
        for(int i = 0 ; i < m_fields.num_elements(); i++)
	{
            SpatialDomains::ConstExactSolutionShPtr ifunc = bcs.GetExactSolution(i);
            for(int j = 0; j < nq; j++)
	    {
                (m_fields[i]->UpdatePhys())[j] = ifunc->Evaluate(x0[j],x1[j],x2[j]);
	    }
            m_fields[i]->SetPhysState(true);
	}
    }
    
    
    void AdvectionDiffusionReaction::AdvectionOperation(const Array<OneD,const NekDouble>& a, 
                                                        const Array<OneD,const NekDouble>& b)
    {
        // not general enough!
      
        int nq = m_fields[0]->GetPointsTot();
      

        for(int i = 0 ; i < m_fields.num_elements(); i++)
	{
            Array<OneD, NekDouble> fluxvector0(nq);
            Array<OneD, NekDouble> fluxvector1(nq);
            Array<OneD, NekDouble> dummy(nq);
	  
            // compute the fluxvector
            Vmath::Vmul(nq,a,1,m_fields[i]->GetPhys(),1,fluxvector0,1);
            Vmath::Vmul(nq,b,1,m_fields[i]->GetPhys(),1,fluxvector1,1);
	  
            // differentiate
            m_fields[i]->ExpList::PhysDeriv(fluxvector0,fluxvector0,NullNekDouble1DArray);
            m_fields[i]->ExpList::PhysDeriv(fluxvector1,dummy,fluxvector1);
	  
            // add the derivatives
            Vmath::Vadd(nq,fluxvector0,1,fluxvector1,1,fluxvector0,1);
	  
            // negate due to move to rhs 
            Vmath::Smul(nq,-1.0,fluxvector0,1,m_fields[i]->UpdatePhys(),1);
	  
            m_fields[i]->SetPhysState(true);
	}
      
    }
    
    void AdvectionDiffusionReaction::FwdTrans(void)
    {  
        int i;
      
        for(i = 0 ; i < m_fields.num_elements(); i++)
	{
            m_fields[i]->FwdTrans(*m_fields[i]);
	}
      
    }
    
    void AdvectionDiffusionReaction::BwdTrans(void)
    {  
      int i;
      
      for(i = 0 ; i < m_fields.num_elements(); i++)
	{
            m_fields[i]->BwdTrans(*m_fields[i]);
	}
      
    }
    
    
    void AdvectionDiffusionReaction::BwdTrans(int field_no)
    {  
      m_fields[field_no]->BwdTrans(*m_fields[field_no]);
    }

    void AdvectionDiffusionReaction::BwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                              Array<OneD, NekDouble> &outarray, int field_no)
    {
        //for(int i = 0 ; i < m_fields.num_elements(); i++)
        {
            m_fields[0]->BwdTrans(inarray,outarray);
        }
    }
    
    void AdvectionDiffusionReaction::BwdTrans(Array<OneD, Array<OneD, NekDouble> > &outarray)
    {
        for(int i = 0 ; i < m_fields.num_elements(); i++)
	{
            m_fields[0]->BwdTrans(m_fields[i]->GetCoeffs(),outarray[i]);
	}
    }
    
    
  
    void AdvectionDiffusionReaction::GetCoords(Array<OneD, NekDouble>& x0,
                                               Array<OneD, NekDouble>& x1,
                                               Array<OneD, NekDouble>& x2)
    {
        m_fields[0]->GetCoords(x0,x1,x2);
    }

     
    void AdvectionDiffusionReaction::SetPhys(Array<OneD, Array<OneD, NekDouble> >&inarray, int field_no)
    {
        int i;
      
        // check that in array large enough
      
      
        // check that field_no no larger than m_fields.num_elements()
	
	
        if (field_no != -1)
	  {
            m_fields[field_no]->SetPhys(inarray[field_no]);
	  }
        else
	  {
            for (i = 0; i < m_fields.num_elements(); ++i)
	      {
                m_fields[i]->SetPhys(inarray[i]);
	      }
	  }
    }
  
    void AdvectionDiffusionReaction::SetPhys(Array<OneD, NekDouble> &inarray, int field_no)
    {
      int i;
	  
      // check that field_no no larger than m_fields.num_elements()
      
      
      if (field_no != -1)
	{
	  m_fields[field_no]->SetPhys(inarray);
	}
    }
  
    const Array<OneD, const NekDouble> &AdvectionDiffusionReaction::GetPhys(int field_no)
    {
        return m_fields[field_no]->GetPhys();
    }
  
  
    void AdvectionDiffusionReaction::GetPhys(Array<OneD, Array<OneD, NekDouble> >&outarray)
    {
      // check size
      
      for (int i = 0; i < m_fields.num_elements(); ++i)
	{
	  Vmath::Vcopy(GetPointsTot(),m_fields[i]->GetPhys(),1,outarray[i],1);
	  //outarray[i] = m_fields[i]->GetPhys();
	}
    }
  
    void AdvectionDiffusionReaction::WriteToFile(std::ofstream &out, OutputFormat format, int field_no)
    {
        m_fields[field_no]->WriteToFile(out,format);
    }

    NekDouble AdvectionDiffusionReaction::L2(const MultiRegions::ExpList2D &In, int field_no)
    {
        return m_fields[field_no]->L2(In);
    }

    
    void AdvectionDiffusionReaction::ExtractTracePhys(Array<OneD, NekDouble> &out, int field_no)
    {
        m_fields[field_no]->ExtractTracePhys(out);
    }   


    /**
     * Fwd 
     * Bwd
     * field_no gives the 
     */
    void AdvectionDiffusionReaction::GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd, 
                                                        Array<OneD,NekDouble> &Bwd,
                                                        int field_no)
    {
        m_fields[field_no]->GetFwdBwdTracePhys(Fwd,Bwd);
    }
    
    void AdvectionDiffusionReaction::GetTraceNormals(Array<OneD, Array<OneD, NekDouble> > &Normals)
    {        
        Array<OneD,NekDouble> normals; 
      
        m_fields[0]->GetTrace()->GetNormals(Normals);
    }
    
    void AdvectionDiffusionReaction::UpwindTrace(Array<OneD,Array<OneD, const NekDouble> > &Vec, 
                                                 Array<OneD, const NekDouble> &Fwd, 
                                                 Array<OneD, const NekDouble> &Bwd, 
                                                 Array<OneD, NekDouble> &Upwind)
    {
        m_fields[0]->GetTrace()->Upwind(Vec, Fwd, Bwd, Upwind);
    }
    
    void AdvectionDiffusionReaction::IProductWRTDerivBase(const int dir, const Array< OneD, const NekDouble > &in,
                                                          Array< OneD, NekDouble > &out, int field_no)
    {
        m_fields[field_no]->IProductWRTDerivBase(dir,in,out);
    }
    
    void AdvectionDiffusionReaction::IProductWRTBase(const Array< OneD, const NekDouble > &in,
                                                     Array< OneD, NekDouble > &out, int field_no)
    {
        m_fields[field_no]->IProductWRTBase(in,out);
    }

    void AdvectionDiffusionReaction::MultiplyByElmtInvMass(const Array< OneD, const NekDouble > &in,
                                                           Array< OneD, NekDouble > &out, int field_no)
    {
        m_fields[field_no]->MultiplyByElmtInvMass(in,out);
    }
    
    
    
    void AdvectionDiffusionReaction::AddTraceIntegral(Array<OneD, const NekDouble> &Fx, 
                                                      Array<OneD, const NekDouble> &Fy, 
                                                      Array<OneD, NekDouble> &outarray,
                                                      int field_no)
    {
        m_fields[field_no]->AddTraceIntegral(Fx,Fy,outarray);
    }
    
    Array<OneD, NekDouble> &AdvectionDiffusionReaction::UpdateCoeffs(int field_no)
    {
        return m_fields[field_no]->UpdateCoeffs();
    }
    
    const Array<OneD, const NekDouble> &AdvectionDiffusionReaction::GetCoeffs(int field_no)
    {
        return m_fields[field_no]->GetCoeffs();
    }

    void AdvectionDiffusionReaction::SetCoeffs(const Array<OneD, const NekDouble > &coeffs, int field_no)
    {
      for(int j = 0; j < m_fields[field_no]->GetNcoeffs(); j++)
	{
	  (m_fields[field_no]->UpdateCoeffs())[j] = coeffs[j];
	}
    }
  

} //end of namespace

/**
* $Log: AdvectionDiffusionReaction.cpp,v $
* Revision 1.1  2008/08/22 09:48:23  pvos
* Added Claes' AdvectionDiffusionReaction, ShallowWater and Euler solver
*
**/
