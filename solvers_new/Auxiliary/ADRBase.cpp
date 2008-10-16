///////////////////////////////////////////////////////////////////////////////
//
// File ADRBase.cpp
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
// Description: Base definitions definiton for
// AdvectionReactiondiffusion, Euler and ShallowWater classes.
//
///////////////////////////////////////////////////////////////////////////////

#include <../solvers_new/Auxiliary/ADRBase.h>
#include <cstdio>
#include <cstdlib>

#include <string>

namespace Nektar
{
    /**
     * Basic construnctor
     */
    ADRBase::ADRBase(void):
        m_fields(0)
    {     
    }
    
    /**
     * Constructor. Creates ... of #DisContField2D fields
     *
     * \param 
     * \param
     */
    ADRBase::ADRBase(string &fileNameString, bool UseInputFileForProjectionType,
                     bool UseContinuousField)
    {
        SpatialDomains::MeshGraph graph; 

        // Both the geometry and the expansion information should be read
        SpatialDomains::MeshGraphSharedPtr mesh = graph.Read(fileNameString);

        // Also read and store the boundary conditions
        // SpatialDomains::BoundaryConditions boundaryConds(&mesh); 
        // boundaryConds.Read(fileNameString);        
        SpatialDomains::MeshGraph *meshptr = mesh.get();
        m_boundaryConditions = MemoryManager<SpatialDomains::BoundaryConditions>::AllocateSharedPtr(meshptr);
        m_boundaryConditions->Read(fileNameString);
        
        // set space dimension for use in class
        m_spacedim = mesh->GetSpaceDimension();
        // save the input file name for output details. 
        m_sessionName = fileNameString;
        m_sessionName = m_sessionName.substr(0,m_sessionName.find_first_of(".")); // Pull out ending
        

        // Options to determine type of projection from file or
        // directly from constructor
        if(UseInputFileForProjectionType == true)
        {
            m_projectionType = eGalerkin;
            
            if(m_boundaryConditions->CheckForParameter("DisContinuous") == true)
            {
                if((int) m_boundaryConditions->GetParameter("DisContinuous") != 0)
                {
                    m_projectionType = eDiscontinuousGalerkin;
                }
            }
        }
        else 
        {
            if(UseContinuousField == true)
            {
                m_projectionType == eGalerkin;
                
            }
            else
            {
                m_projectionType == eDiscontinuousGalerkin;
            }
        }


        SetADRBase(mesh,m_boundaryConditions->GetNumVariables());
    }
    
    void ADRBase::SetADRBase(SpatialDomains::MeshGraphSharedPtr &mesh,
                             int nvariables)
    {
        int i;

        m_fields = Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);
        m_spacedim = mesh->GetSpaceDimension();

        // May want to attach this information in future. 
        int expdim  = mesh->GetMeshDimension();

        if(m_projectionType == eGalerkin)
        {
            switch(expdim)
            {
            case 1:
                {
                    SpatialDomains::MeshGraph1DSharedPtr mesh1D;

                    if(!(mesh1D = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph1D>(mesh)))
                    {
                        ASSERTL0(false,"Dynamics cast failed");
                    }

                    for(i = 0 ; i < m_fields.num_elements(); i++)
                    {
                        m_fields[i] = MemoryManager<MultiRegions::ContField1D>::AllocateSharedPtr(*mesh1D,*m_boundaryConditions);
                    }
                }
                break;
            case 2:
                {
                    SpatialDomains::MeshGraph2DSharedPtr mesh2D;
                    
                    if(!(mesh2D = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(mesh)))
                    {
                        ASSERTL0(false,"Dynamics cast failed");
                    }

                    for(i = 0 ; i < m_fields.num_elements(); i++)
                    {
                        m_fields[i] = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(*mesh2D,*m_boundaryConditions);
                    }
                    break;
                }
            case 3:
                ASSERTL0(false,"3 D not set up");
            default:
                ASSERTL0(false,"Expansion dimension not recognised");
                break;
            }
        }
        else // Discontinuous Field
        {
          
            switch(expdim)
            {
            case 1:
                {
                    SpatialDomains::MeshGraph1DSharedPtr mesh1D;
                    
                    if(!(mesh1D = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph1D>(mesh)))
                    {
                        ASSERTL0(false,"Dynamics cast failed");
                    }
                    
                    for(i = 0 ; i < m_fields.num_elements(); i++)
                    {
                        m_fields[i] = MemoryManager<MultiRegions::DisContField1D>::AllocateSharedPtr(*mesh1D,*m_boundaryConditions,i);
                    }
                }
                break;
            case 2:
                {
                    SpatialDomains::MeshGraph2DSharedPtr mesh2D;
                    
                    if(!(mesh2D = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(mesh)))
                    {
                        ASSERTL0(false,"Dynamics cast failed");
                    }
                    
                    for(i = 0 ; i < m_fields.num_elements(); i++)
                    {
                        m_fields[i] = MemoryManager<MultiRegions::DisContField2D>::AllocateSharedPtr(*mesh2D,*m_boundaryConditions,i);
                    }
                }
                break;
            case 3:
                ASSERTL0(false,"3 D not set up");
            default:
                ASSERTL0(false,"Expansion dimension not recognised");
                break;
            }

            // Set up Normals. 
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
            
            for(i = 0; i < m_spacedim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (m_fields[0]->GetTrace()->GetNpoints());
            }
            
            m_fields[0]->GetTrace()->GetNormals(m_traceNormals);
        }
        
        // Set Default Parameter

        if(m_boundaryConditions->CheckForParameter("Time") == true)
        {
            m_time  = m_boundaryConditions->GetParameter("Time");
        }
        else
        {
            m_time  = 0.0;
        }
        m_timestep   = m_boundaryConditions->GetParameter("TimeStep");
        m_steps      = m_boundaryConditions->GetParameter("Steps");
        
        if(m_boundaryConditions->CheckForParameter("CheckSteps") == true)
        {
            m_checksteps = m_boundaryConditions->GetParameter("CheckSteps");
        }
        else
        {
            m_checksteps = m_steps;
        }
    }

    void ADRBase::SetInitialConditions(NekDouble initialtime)
    {
        int nq = m_fields[0]->GetPointsTot();
      
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
      
        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);
      
        for(int i = 0 ; i < m_fields.num_elements(); i++)
	{
            SpatialDomains::ConstInitialConditionShPtr ifunc = m_boundaryConditions->GetInitialCondition(i);
            for(int j = 0; j < nq; j++)
	    {
                (m_fields[i]->UpdatePhys())[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],initialtime);
	    }
            m_fields[i]->SetPhysState(true);

            m_fields[i]->FwdTrans(*(m_fields[i]));
	}
    }
    
    void ADRBase::EvaluateExactSolution(int field, Array<OneD, NekDouble> &outfield)
    {
        int nq = m_fields[field]->GetPointsTot();
      
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
      
        // get the coordinates of the quad points
        m_fields[field]->GetCoords(x0,x1,x2);
      
        SpatialDomains::ConstExactSolutionShPtr ifunc = m_boundaryConditions->GetExactSolution(field);
        for(int j = 0; j < nq; j++)
        {
            outfield[j] = ifunc->Evaluate(x0[j],x1[j],x2[j]);
        }
    }
    
    void ADRBase::EvaluateUserDefinedEqn(Array<OneD, Array<OneD, NekDouble> > &outfield)
    {
        int nq = m_fields[0]->GetPointsTot();
        
        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);
      
        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);
      
        for(int i = 0 ; i < m_fields.num_elements(); i++)
	{
            SpatialDomains::ConstUserDefinedEqnShPtr ifunc = m_boundaryConditions->GetUserDefinedEqn(i);
            for(int j = 0; j < nq; j++)
	    {
                outfield[i][j] = ifunc->Evaluate(x0[j],x1[j],x2[j]);
	    }
	}
        
    }

    NekDouble ADRBase::L2Error(int field)
    {
        Array<OneD, NekDouble> exactsoln(m_fields[field]->GetPointsTot());
        
        EvaluateExactSolution(field,exactsoln);
        
        return m_fields[field]->L2(exactsoln);
    }


    //-------------------------------------------------------------
    // Compute weak Green form of advection terms (without boundary
    // integral, i.e (\grad \phi \cdot F) where for example F = uV
    //
    // Note: Assuming all fields are of the same expansion and order
    // so that we can use the parameters of m_fields[0].
    // ------------------------------------------------------------

    
    void ADRBase::WeakAdvectionGreensDivergenceForm(const Array<OneD, Array<OneD, NekDouble> > &F, Array<OneD, NekDouble> &outarray)
    {
        // use dimension of Velocity vector to dictate dimension of operation
        int ndim    = F.num_elements();
        int nCoeffs = m_fields[0]->GetNcoeffs();

        Array<OneD, NekDouble> iprod(nCoeffs);
        Vmath::Zero(nCoeffs,outarray,1);
        
        for (int i = 0; i < ndim; ++i)
        {
            m_fields[0]->IProductWRTDerivBase(i,F[i],iprod);
            Vmath::Vadd(nCoeffs,iprod,1,outarray,1,outarray,1);
        }
       
    }

    //-------------------------------------------------------------
    // Calculate Inner product of the divergence advection form
    // .....(\phi, Div \cdot F) where for example F = uV
    // -------------------------------------------------------------
    
    void ADRBase::WeakAdvectionDivergenceForm(const Array<OneD, Array<OneD, NekDouble> > &F, Array<OneD, NekDouble> &outarray)
    {
        // use dimension of Velocity vector to dictate dimension of operation
        int ndim       = F.num_elements();
        int nPointsTot = m_fields[0]->GetPointsTot();
        Array<OneD, NekDouble> tmp(nPointsTot);
        Array<OneD, NekDouble> div(nPointsTot,0.0);
        
        // Evaluate the divergence 
        for(int i = 0; i < ndim; ++i)
        {
            m_fields[0]->PhysDeriv(i,F[i],tmp);
            Vmath::Vadd(nPointsTot,tmp,1,div,1,div,1);
        }

        m_fields[0]->IProductWRTBase(div,outarray);
    }

    //-------------------------------------------------------------
    // Calculate Inner product of the divergence advection form
    // ..... (\phi, V\cdot Grad(u))
    // -------------------------------------------------------------

    void ADRBase::WeakAdvectionNonConservativeForm(const Array<OneD, Array<OneD, NekDouble> > &V, const Array<OneD, const NekDouble> &u, Array<OneD, NekDouble> &outarray)
    {
        // use dimension of Velocity vector to dictate dimension of operation
        int ndim       = V.num_elements();
        int nPointsTot = m_fields[0]->GetPointsTot();
        Array<OneD, NekDouble> tmp(nPointsTot);
        Array<OneD, NekDouble> grad0(ndim*nPointsTot,0.0),grad1,grad2;

        // Evaluate V\cdot Grad(u)
        switch(ndim)
        {
        case 1:
            m_fields[0]->PhysDeriv(u,grad0);
            Vmath::Vmul(nPointsTot,grad0,1,V[0],1,tmp,1);
            break;
        case 2:
            grad1 = grad0 + nPointsTot;
            m_fields[0]->PhysDeriv(u,grad0,grad1);
            Vmath::Vmul (nPointsTot,grad0,1,V[0],1,tmp,1);
            Vmath::Vvtvp(nPointsTot,grad1,1,V[1],1,tmp,1,tmp,1);
            break;
        case 3:
            grad1 = grad0 + nPointsTot;
            grad2 = grad1 + nPointsTot;
            m_fields[0]->PhysDeriv(u,grad0,grad1,grad2);
            Vmath::Vmul (nPointsTot,grad0,1,V[0],1,tmp,1);
            Vmath::Vvtvp(nPointsTot,grad1,1,V[1],1,tmp,1,tmp,1);
            Vmath::Vvtvp(nPointsTot,grad2,1,V[2],1,tmp,1,tmp,1);
            break;
        default:
            ASSERTL0(false,"dimension unknown");
        }

        m_fields[0]->IProductWRTBase(tmp,outarray);
    }
                                       
    //-------------------------------------------------------------
    // Calculate weak DG advection in the form 
    //  <\phi, \hat{F}\cdot n> - (\grad \phi \cdot F)
    // -------------------------------------------------------------
    void ADRBase::WeakDGAdvection(Array<OneD, Array<OneD, NekDouble> >& OutField )
    {
        int i;
        int nVelDim = m_spacedim;
        int nPointsTot = GetPointsTot();
        int ncoeffs    = GetNcoeffs();
        int nTracePointsTot = GetTracePointsTot();
        int nvariables = m_fields.num_elements();

        Array<OneD, Array<OneD, NekDouble> > flux    (nVelDim);
        Array<OneD, Array<OneD, NekDouble> > numflux (nvariables);
        
        for(i = 0; i < nVelDim; ++i)
        {
            flux[i]    = Array<OneD, NekDouble>(nPointsTot);
        }

        for(i = 0; i < nvariables; ++i)
        {
            numflux[i] = Array<OneD, NekDouble>(nTracePointsTot);
        }
        
        for(i = 0; i < nvariables; ++i)
        {
            // Get the ith component of the  flux vector in (physical space)
            GetFluxVector(i, flux);
            
            // Calcuate the i^th value of (\grad_i \phi, F)
            WeakAdvectionGreensDivergenceForm(flux,OutField[i]);
        }
        
        // Evaluate numerical flux in physical space which may in
        // general couple all component of vectors
        NumericalFlux(numflux);
        
        // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i] 
        for(i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(ncoeffs,OutField[i],1);
            m_fields[i]->AddTraceIntegral(numflux[i],OutField[i]);
        }
    }

    
    void ADRBase::FwdTrans(const ADRBase &In)
    {  
        int i;
      
        for(i = 0 ; i < m_fields.num_elements(); i++)
	{
            m_fields[i]->FwdTrans(*m_fields[i]);
	}
      
    }
    
    void ADRBase::BwdTrans(void)
    {  
        int i;
      
        for(i = 0 ; i < m_fields.num_elements(); i++)
	{
            m_fields[i]->BwdTrans(*m_fields[i]);
	}
      
    }
    
    void ADRBase::BwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                              Array<OneD, NekDouble> &outarray, int field_no)
    {
        //for(int i = 0 ; i < m_fields.num_elements(); i++)
        {
            m_fields[0]->BwdTrans(inarray,outarray);
        }
    }
    
    void ADRBase::BwdTrans(Array<OneD, Array<OneD, NekDouble> > &outarray)
    {
        for(int i = 0 ; i < m_fields.num_elements(); i++)
	{
            m_fields[0]->BwdTrans(m_fields[i]->GetCoeffs(),outarray[i]);
	}
    }
    
    
  
    void ADRBase::GetCoords(Array<OneD, NekDouble>& x0,
                                               Array<OneD, NekDouble>& x1,
                                               Array<OneD, NekDouble>& x2)
    {
        m_fields[0]->GetCoords(x0,x1,x2);
    }

     
    void ADRBase::SetPhys(Array<OneD, Array<OneD, NekDouble> >&inarray, int field_no)
    {
        int i;
      
        // check that in array large enough
      
      
        // check that field_no no larger than m_fields.num_elements()
      
	  
        if (field_no != -1)
	{
            m_fields[field_no]->SetPhys(inarray[0]);
	}
        else
	{
            for (i = 0; i < m_fields.num_elements(); ++i)
	    {
                m_fields[i]->SetPhys(inarray[i]);
	    }
	}
      
    }

    void ADRBase::SetPhys(Array<OneD, NekDouble> &inarray, int field_no)
    {
        int i;
      
        // check that field_no no larger than m_fields.num_elements()
	
	  
        if (field_no != -1)
	{
            m_fields[field_no]->SetPhys(inarray);
	}
    }
    
    const Array<OneD, const NekDouble> &ADRBase::GetPhys(int field_no)
    {
        return m_fields[field_no]->GetPhys();
    }

  
    void ADRBase::GetPhys(Array<OneD, Array<OneD, NekDouble> >&outarray)
    {
        // check size
      
        for (int i = 0; i < m_fields.num_elements(); ++i)
	{
            outarray[i] = m_fields[i]->GetPhys();
	}

    }
    
    void ADRBase::WriteToFile(std::ofstream &out, OutputFormat format, int field_no)
    {
        m_fields[field_no]->WriteToFile(out,format);
    }

    NekDouble ADRBase::L2(const MultiRegions::ExpList2D &In, int field_no)
    {
        return m_fields[field_no]->L2(In);
    }

    
    void ADRBase::ExtractTracePhys(Array<OneD, NekDouble> &out, int field_no)
    {
        m_fields[field_no]->ExtractTracePhys(out);
    }   


    /**
     * Fwd 
     * Bwd
     * field_no gives the 
     */
    void ADRBase::GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd, 
                                                        Array<OneD,NekDouble> &Bwd,
                                                        int field_no)
    {
        m_fields[field_no]->GetFwdBwdTracePhys(Fwd,Bwd);
    }
    
    void ADRBase::GetTraceNormals(Array<OneD, Array<OneD, NekDouble> > &Normals)
    {        
        Array<OneD,NekDouble> normals; 
      
        m_fields[0]->GetTrace()->GetNormals(Normals);
    }
    
    void ADRBase::UpwindTrace(const Array<OneD,Array<OneD,NekDouble> > &Vel, 
                              const Array<OneD, const NekDouble> &Fwd, 
                              const Array<OneD, const NekDouble> &Bwd, 
                              Array<OneD, NekDouble> &Upwind)
    {
        m_fields[0]->GetTrace()->Upwind(Vel, Fwd, Bwd, Upwind);
    }
    
    void ADRBase::IProductWRTDerivBase(const int dir, const Array< OneD, const NekDouble > &in,
                                                          Array< OneD, NekDouble > &out, int field_no)
    {
        m_fields[field_no]->IProductWRTDerivBase(dir,in,out);
    }
    
    void ADRBase::IProductWRTBase(const Array< OneD, const NekDouble > &in,
                                                     Array< OneD, NekDouble > &out, int field_no)
    {
        m_fields[field_no]->IProductWRTBase(in,out);
    }

    void ADRBase::MultiplyByElmtInvMass(const Array< OneD, const NekDouble > &in,
                                                           Array< OneD, NekDouble > &out, int field_no)
    {
        m_fields[field_no]->MultiplyByElmtInvMass(in,out);
    }
    
    
    
    void ADRBase::AddTraceIntegral(Array<OneD, const NekDouble> &Fx, 
                                   Array<OneD, const NekDouble> &Fy, 
                                   Array<OneD, NekDouble> &outarray,
                                   int field_no)
    {
        m_fields[field_no]->AddTraceIntegral(Fx,Fy,outarray);
    }
    
    Array<OneD, NekDouble> &ADRBase::UpdateCoeffs(int field_no)
    {
        return m_fields[field_no]->UpdateCoeffs();
    }
    
    const Array<OneD, const NekDouble> &ADRBase::GetCoeffs(int field_no)
    {
        return m_fields[field_no]->GetCoeffs();
    }

    void ADRBase::Output(void)
    {
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            std::string outname = m_sessionName +"_" + m_boundaryConditions->GetVariable(i) +  ".dat";
            ofstream outfile(outname.c_str());
            m_fields[i]->BwdTrans(*(m_fields[i]));
            m_fields[i]->WriteToFile(outfile,eTecplot);
        }
    }

    void ADRBase::Checkpoint_Output(const int n)
    {
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            char chkout[16] = "";
            sprintf(chkout, "%d", n);
            std::string outname = m_sessionName +"_" + m_boundaryConditions->GetVariable(i) + "_" + chkout + ".chk";
            ofstream outfile(outname.c_str());
            m_fields[i]->BwdTrans(*(m_fields[i]));
            m_fields[i]->WriteToFile(outfile,eTecplot);
        }
    }

} //end of namespace

/**
* $Log:$
**/
