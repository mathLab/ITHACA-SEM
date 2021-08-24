///////////////////////////////////////////////////////////////////////////////
//
// File NonlinearSWE.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Nonlinear Shallow water equations in conservative variables
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <boost/algorithm/string.hpp>

#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>
#include <ShallowWaterSolver/EquationSystems/NonlinearSWE.h>

using namespace std;

namespace Nektar
{
  string NonlinearSWE::className =
     SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
	  "NonlinearSWE", NonlinearSWE::create,
          "Nonlinear shallow water equation in conservative variables.");

  NonlinearSWE::NonlinearSWE(
      const LibUtilities::SessionReaderSharedPtr& pSession,
      const SpatialDomains::MeshGraphSharedPtr& pGraph)
      : ShallowWaterSystem(pSession, pGraph)
  {
  }

  void NonlinearSWE::v_InitObject()
  {
      ShallowWaterSystem::v_InitObject();

    if (m_explicitAdvection)
      {
	m_ode.DefineOdeRhs     (&NonlinearSWE::DoOdeRhs,        this);
	m_ode.DefineProjection (&NonlinearSWE::DoOdeProjection, this);
      }
    else
      {
	ASSERTL0(false, "Implicit SWE not set up.");
      }

    // Type of advection class to be used
    switch(m_projectionType)
      {
	// Continuous field
      case MultiRegions::eGalerkin:
	{
	  //  Do nothing
	  break;
	}
	// Discontinuous field
      case MultiRegions::eDiscontinuous:
	{
	  string advName;
	  string diffName;
	  string riemName;

	  //---------------------------------------------------------------
	  // Setting up advection and diffusion operators
	  // NB: diffusion not set up for SWE at the moment
	  //     but kept here for future use ...
	  m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
	  // m_session->LoadSolverInfo("DiffusionType", diffName, "LDGEddy");
	  m_advection = SolverUtils::GetAdvectionFactory()
	    .CreateInstance(advName, advName);
	  // m_diffusion = SolverUtils::GetDiffusionFactory()
	  //                             .CreateInstance(diffName, diffName);

	  m_advection->SetFluxVector(&NonlinearSWE::GetFluxVector, this);
	  // m_diffusion->SetFluxVectorNS(&ShallowWaterSystem::
	  //                                  GetEddyViscosityFluxVector, this);

	  // Setting up Riemann solver for advection operator
	  m_session->LoadSolverInfo("UpwindType", riemName, "Average");
	  m_riemannSolver = SolverUtils::GetRiemannSolverFactory()
	    .CreateInstance(riemName, m_session);

	  // Setting up upwind solver for diffusion operator
	  // m_riemannSolverLDG = SolverUtils::GetRiemannSolverFactory()
	  //                                 .CreateInstance("UpwindLDG");

	  // Setting up parameters for advection operator Riemann solver
	  m_riemannSolver->SetParam (
                                     "gravity",
                                     &NonlinearSWE::GetGravity, this);
          m_riemannSolver->SetAuxVec(
                                     "vecLocs",
                                     &NonlinearSWE::GetVecLocs, this);
	  m_riemannSolver->SetVector(
				     "N",
				     &NonlinearSWE::GetNormals, this);
	  m_riemannSolver->SetScalar(
				     "depth",
				     &NonlinearSWE::GetDepth, this);

	  // Setting up parameters for diffusion operator Riemann solver
	  // m_riemannSolverLDG->AddParam (
	  //                     "gravity",
	  //                     &NonlinearSWE::GetGravity,   this);
      // m_riemannSolverLDG->SetAuxVec(
      //                     "vecLocs",
      //                     &NonlinearSWE::GetVecLocs,  this);
	  // m_riemannSolverLDG->AddVector(
	  //                     "N",
	  //                     &NonlinearSWE::GetNormals, this);

	  // Concluding initialisation of advection / diffusion operators
	  m_advection->SetRiemannSolver   (m_riemannSolver);
	  //m_diffusion->SetRiemannSolver   (m_riemannSolverLDG);
	  m_advection->InitObject         (m_session, m_fields);
	  //m_diffusion->InitObject         (m_session, m_fields);
	  break;
	}
      default:
	{
	  ASSERTL0(false, "Unsupported projection type.");
	  break;
	}
      }


  }

  NonlinearSWE::~NonlinearSWE()
  {

  }

  // physarray contains the conservative variables
  void NonlinearSWE::AddCoriolis(const Array<OneD, const Array<OneD, NekDouble> > &physarray,
			      Array<OneD, Array<OneD, NekDouble> > &outarray)
  {

    int ncoeffs    = GetNcoeffs();
    int nq         = GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> mod(ncoeffs);

    switch(m_projectionType)
      {
      case MultiRegions::eDiscontinuous:
	{
	  // add to hu equation
	  Vmath::Vmul(nq,m_coriolis,1,physarray[2],1,tmp,1);
	  m_fields[0]->IProductWRTBase(tmp,mod);
	  m_fields[0]->MultiplyByElmtInvMass(mod,mod);
	  m_fields[0]->BwdTrans(mod,tmp);
	  Vmath::Vadd(nq,tmp,1,outarray[1],1,outarray[1],1);

	  // add to hv equation
	  Vmath::Vmul(nq,m_coriolis,1,physarray[1],1,tmp,1);
	  Vmath::Neg(nq,tmp,1);
	  m_fields[0]->IProductWRTBase(tmp,mod);
	  m_fields[0]->MultiplyByElmtInvMass(mod,mod);
	  m_fields[0]->BwdTrans(mod,tmp);
	  Vmath::Vadd(nq,tmp,1,outarray[2],1,outarray[2],1);
	}
	break;
      case MultiRegions::eGalerkin:
      case MultiRegions::eMixed_CG_Discontinuous:
	{
	  // add to hu equation
	  Vmath::Vmul(nq,m_coriolis,1,physarray[2],1,tmp,1);
	  Vmath::Vadd(nq,tmp,1,outarray[1],1,outarray[1],1);

	  // add to hv equation
	  Vmath::Vmul(nq,m_coriolis,1,physarray[1],1,tmp,1);
	  Vmath::Neg(nq,tmp,1);
	  Vmath::Vadd(nq,tmp,1,outarray[2],1,outarray[2],1);
	}
	break;
      default:
	ASSERTL0(false,"Unknown projection scheme for the NonlinearSWE");
	break;
      }


  }


  // physarray contains the conservative variables
  void NonlinearSWE::AddVariableDepth(const Array<OneD, const Array<OneD, NekDouble> > &physarray,
				      Array<OneD, Array<OneD, NekDouble> > &outarray)
  {

    int ncoeffs    = GetNcoeffs();
    int nq         = GetTotPoints();

    Array<OneD, NekDouble> tmp(nq);
    Array<OneD, NekDouble> mod(ncoeffs);

    switch(m_projectionType)
      {
      case MultiRegions::eDiscontinuous:
	{
	  for (int i = 0; i < m_spacedim; ++i)
	    {
	      Vmath::Vmul(nq,m_bottomSlope[i],1,physarray[0],1,tmp,1);
	      Vmath::Smul(nq,m_g,tmp,1,tmp,1);
	      m_fields[0]->IProductWRTBase(tmp,mod);
	      m_fields[0]->MultiplyByElmtInvMass(mod,mod);
	      m_fields[0]->BwdTrans(mod,tmp);
	      Vmath::Vadd(nq,tmp,1,outarray[i+1],1,outarray[i+1],1);
	    }
	}
	break;
      case MultiRegions::eGalerkin:
      case MultiRegions::eMixed_CG_Discontinuous:
	{
	 for (int i = 0; i < m_spacedim; ++i)
	    {
	      Vmath::Vmul(nq,m_bottomSlope[i],1,physarray[0],1,tmp,1);
	      Vmath::Smul(nq,m_g,tmp,1,tmp,1);
	      Vmath::Vadd(nq,tmp,1,outarray[i+1],1,outarray[i+1],1);
	    }
	}
	break;
      default:
	ASSERTL0(false,"Unknown projection scheme for the NonlinearSWE");
	break;
      }


  }

  void NonlinearSWE::DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
			            Array<OneD,       Array<OneD, NekDouble> >&outarray,
			      const NekDouble time)
  {
    int i, j;
    int ndim       = m_spacedim;
    int nvariables = inarray.size();
    int nq         = GetTotPoints();


    switch(m_projectionType)
      {
      case MultiRegions::eDiscontinuous:
	{

	  //-------------------------------------------------------
	  // Compute the DG advection including the numerical flux
	  // by using SolverUtils/Advection
	  // Input and output in physical space
	  Array<OneD, Array<OneD, NekDouble> > advVel;

	  m_advection->Advect(nvariables, m_fields, advVel, inarray,
	                      outarray, time);
	  //-------------------------------------------------------


	  //-------------------------------------------------------
	  // negate the outarray since moving terms to the rhs
	  for(i = 0; i < nvariables; ++i)
	    {
	      Vmath::Neg(nq,outarray[i],1);
	    }
	  //-------------------------------------------------------


	  //-------------------------------------------------
	  // Add "source terms"
	  // Input and output in physical space

	  // Coriolis forcing
	  if (m_coriolis.size() != 0)
	    {
	      AddCoriolis(inarray,outarray);
	    }

	  // Variable Depth
	  if (m_constantDepth != true)
	    {
	      AddVariableDepth(inarray,outarray);
	    }
	  //-------------------------------------------------

	}
	break;
      case MultiRegions::eGalerkin:
      case MultiRegions::eMixed_CG_Discontinuous:
	{

	  //-------------------------------------------------------
	  // Compute the fluxvector in physical space
	  Array<OneD, Array<OneD, Array<OneD, NekDouble> > >
	    fluxvector(nvariables);

	  for (i = 0; i < nvariables; ++i)
            {
	      fluxvector[i] = Array<OneD, Array<OneD, NekDouble> >(ndim);
	      for(j = 0; j < ndim; ++j)
                {
		  fluxvector[i][j] = Array<OneD, NekDouble>(nq);
                }
            }

	  NonlinearSWE::GetFluxVector(inarray, fluxvector);
	  //-------------------------------------------------------



 	  //-------------------------------------------------------
	  // Take the derivative of the flux terms
	  // and negate the outarray since moving terms to the rhs
	  Array<OneD,NekDouble> tmp(nq);
	  Array<OneD, NekDouble>tmp1(nq);

	  for(i = 0; i < nvariables; ++i)
	    {
	      m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[0],fluxvector[i][0],tmp);
	      m_fields[i]->PhysDeriv(MultiRegions::DirCartesianMap[1],fluxvector[i][1],tmp1);
	      Vmath::Vadd(nq,tmp,1,tmp1,1,outarray[i],1);
	      Vmath::Neg(nq,outarray[i],1);
	    }


	  //-------------------------------------------------
	  // Add "source terms"
	  // Input and output in physical space

	  // Coriolis forcing
	  if (m_coriolis.size() != 0)
	    {
	      AddCoriolis(inarray,outarray);
	    }

	  // Variable Depth
	  if (m_constantDepth != true)
	    {
	      AddVariableDepth(inarray,outarray);
	    }
	  //-------------------------------------------------
	}
	break;
      default:
	ASSERTL0(false,"Unknown projection scheme for the NonlinearSWE");
	break;
      }
  }


  void NonlinearSWE::DoOdeProjection(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
				           Array<OneD,       Array<OneD, NekDouble> >&outarray,
				     const NekDouble time)
  {
    int i;
    int nvariables = inarray.size();


    switch(m_projectionType)
      {
      case MultiRegions::eDiscontinuous:
	{

	  // Just copy over array
	  int npoints = GetNpoints();

	  for(i = 0; i < nvariables; ++i)
	    {
	      Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
	    }
	  SetBoundaryConditions(outarray, time);
	  break;
	}
      case MultiRegions::eGalerkin:
      case MultiRegions::eMixed_CG_Discontinuous:
	{

	  EquationSystem::SetBoundaryConditions(time);
	  Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs(),0.0);

	  for(i = 0; i < nvariables; ++i)
          {
              m_fields[i]->FwdTrans(inarray[i],coeffs);
	      m_fields[i]->BwdTrans_IterPerExp(coeffs,outarray[i]);
          }
	  break;
	}
      default:
	ASSERTL0(false,"Unknown projection scheme");
	break;
      }
  }


   //----------------------------------------------------
  void NonlinearSWE::SetBoundaryConditions(
    Array<OneD, Array<OneD, NekDouble> > &inarray,
    NekDouble time)
  {
      std::string varName;
      int nvariables = m_fields.size();
      int cnt = 0;
      int nTracePts  = GetTraceTotPoints();

      // Extract trace for boundaries. Needs to be done on all processors to avoid
      // deadlock.
      Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
      for (int i = 0; i < nvariables; ++i)
      {
          Fwd[i] = Array<OneD, NekDouble>(nTracePts);
          m_fields[i]->ExtractTracePhys(inarray[i], Fwd[i]);
      }

      // Loop over Boundary Regions
      for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
      {

          if (m_fields[0]->GetBndConditions()[n]->GetBoundaryConditionType()
              == SpatialDomains::ePeriodic)
          {
              continue;
          }

          // Wall Boundary Condition
          if (boost::iequals(m_fields[0]->GetBndConditions()[n]->GetUserDefined(),"Wall"))
          {
              WallBoundary2D(n, cnt, Fwd, inarray);
          }

          // Time Dependent Boundary Condition (specified in meshfile)
          if (m_fields[0]->GetBndConditions()[n]->IsTimeDependent())
          {
              for (int i = 0; i < nvariables; ++i)
              {
                  varName = m_session->GetVariable(i);
                  m_fields[i]->EvaluateBoundaryConditions(time, varName);
              }
          }
          cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
      }
  }

  //----------------------------------------------------
  /**
     * @brief Wall boundary condition.
     */
    void NonlinearSWE::WallBoundary(
        int                                   bcRegion,
        int                                   cnt,
        Array<OneD, Array<OneD, NekDouble> > &Fwd,
        Array<OneD, Array<OneD, NekDouble> > &physarray)
    {
        int i;
        int nvariables = physarray.size();

        // Adjust the physical values of the trace to take
        // user defined boundaries into account
        int e, id1, id2, npts;

        for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]
                 ->GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
            id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetPhys_Offset(e);
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(
                        m_fields[0]->GetTraceMap()->
                                    GetBndCondIDToGlobalTraceID(cnt+e));

            // For 2D/3D, define: v* = v - 2(v.n)n
            Array<OneD, NekDouble> tmp(npts, 0.0);

            // Calculate (v.n)
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(npts,
                             &Fwd[1+i][id2], 1,
                             &m_traceNormals[i][id2], 1,
                             &tmp[0], 1,
                             &tmp[0], 1);
            }

            // Calculate 2.0(v.n)
            Vmath::Smul(npts, -2.0, &tmp[0], 1, &tmp[0], 1);

            // Calculate v* = v - 2.0(v.n)n
            for (i = 0; i < m_spacedim; ++i)
            {
                Vmath::Vvtvp(npts,
                             &tmp[0], 1,
                             &m_traceNormals[i][id2], 1,
                             &Fwd[1+i][id2], 1,
                             &Fwd[1+i][id2], 1);
            }

            // copy boundary adjusted values into the boundary expansion
            for (i = 0; i < nvariables; ++i)
            {
                Vmath::Vcopy(npts, &Fwd[i][id2], 1,
                             &(m_fields[i]->GetBndCondExpansions()[bcRegion]->
                             UpdatePhys())[id1], 1);
            }
        }
    }


  void NonlinearSWE::WallBoundary2D(int bcRegion, int cnt, Array<OneD, Array<OneD, NekDouble> > &Fwd, Array<OneD, Array<OneD, NekDouble> > &physarray)
  {

    int i;
    int nvariables = physarray.size();

    // Adjust the physical values of the trace to take
    // user defined boundaries into account
    int e, id1, id2, npts;

    for(e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize(); ++e)
      {
	npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExp(e)->GetNumPoints(0);
	id1  = m_fields[0]->GetBndCondExpansions()[bcRegion]->GetPhys_Offset(e) ;
	id2  = m_fields[0]->GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt+e));

	switch(m_expdim)
	  {
	  case 1:
	    {
	      // negate the forward flux
	      Vmath::Neg(npts,&Fwd[1][id2],1);
	    }
	    break;
	  case 2:
	    {
	      Array<OneD, NekDouble> tmp_n(npts);
	      Array<OneD, NekDouble> tmp_t(npts);

	      Vmath::Vmul(npts,&Fwd[1][id2],1,&m_traceNormals[0][id2],1,&tmp_n[0],1);
	      Vmath::Vvtvp(npts,&Fwd[2][id2],1,&m_traceNormals[1][id2],1,&tmp_n[0],1,&tmp_n[0],1);

	      Vmath::Vmul(npts,&Fwd[1][id2],1,&m_traceNormals[1][id2],1,&tmp_t[0],1);
	      Vmath::Vvtvm(npts,&Fwd[2][id2],1,&m_traceNormals[0][id2],1,&tmp_t[0],1,&tmp_t[0],1);

	      // negate the normal flux
	      Vmath::Neg(npts,tmp_n,1);

	      // rotate back to Cartesian
	      Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[1][id2],1,&Fwd[1][id2],1);
	      Vmath::Vvtvm(npts,&tmp_n[0],1,&m_traceNormals[0][id2],1,&Fwd[1][id2],1,&Fwd[1][id2],1);

	      Vmath::Vmul(npts,&tmp_t[0],1,&m_traceNormals[0][id2],1,&Fwd[2][id2],1);
	      Vmath::Vvtvp(npts,&tmp_n[0],1,&m_traceNormals[1][id2],1,&Fwd[2][id2],1,&Fwd[2][id2],1);
	    }
	    break;
	  case 3:
	    ASSERTL0(false,"3D not implemented for Shallow Water Equations");
	    break;
	  default:
	    ASSERTL0(false,"Illegal expansion dimension");
	  }



	// copy boundary adjusted values into the boundary expansion
	for (i = 0; i < nvariables; ++i)
	  {
	    Vmath::Vcopy(npts,&Fwd[i][id2], 1,&(m_fields[i]->GetBndCondExpansions()[bcRegion]->UpdatePhys())[id1],1);
	  }
      }
  }


  // Physfield in conservative Form
 void NonlinearSWE::GetFluxVector(
     const Array<OneD, const Array<OneD, NekDouble> > &physfield,
           Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
  {
    int i, j;
    int nq = m_fields[0]->GetTotPoints();

    NekDouble g = m_g;
    Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);

    // Flux vector for the mass equation
    for (i = 0; i < m_spacedim; ++i)
      {
	velocity[i] = Array<OneD, NekDouble>(nq);
	Vmath::Vcopy(nq, physfield[i+1], 1, flux[0][i], 1);
      }

     GetVelocityVector(physfield, velocity);

     // Put (0.5 g h h) in tmp
     Array<OneD, NekDouble> tmp(nq);
     Vmath::Vmul(nq, physfield[0], 1, physfield[0], 1, tmp, 1);
     Vmath::Smul(nq, 0.5*g, tmp, 1, tmp, 1);

     // Flux vector for the momentum equations
     for (i = 0; i < m_spacedim; ++i)
       {
	 for (j = 0; j < m_spacedim; ++j)
	   {
	     Vmath::Vmul(nq, velocity[j], 1, physfield[i+1], 1,
			 flux[i+1][j], 1);
	   }

	 // Add (0.5 g h h) to appropriate field
	 Vmath::Vadd(nq, flux[i+1][i], 1, tmp, 1, flux[i+1][i], 1);
       }

  }

  void NonlinearSWE::ConservativeToPrimitive(const Array<OneD, const Array<OneD, NekDouble> >&physin,
					      Array<OneD,       Array<OneD, NekDouble> >&physout)
  {
    int nq = GetTotPoints();

    if(physin.get() == physout.get())
      {
	// copy indata and work with tmp array
	Array<OneD, Array<OneD, NekDouble> >tmp(3);
	for (int i = 0; i < 3; ++i)
	  {
	    // deep copy
	    tmp[i] = Array<OneD, NekDouble>(nq);
	    Vmath::Vcopy(nq,physin[i],1,tmp[i],1);
	  }

	// \eta = h - d
	Vmath::Vsub(nq,tmp[0],1,m_depth,1,physout[0],1);

	// u = hu/h
	Vmath::Vdiv(nq,tmp[1],1,tmp[0],1,physout[1],1);

	// v = hv/ v
	Vmath::Vdiv(nq,tmp[2],1,tmp[0],1,physout[2],1);
      }
    else
      {
	// \eta = h - d
	Vmath::Vsub(nq,physin[0],1,m_depth,1,physout[0],1);

	// u = hu/h
	Vmath::Vdiv(nq,physin[1],1,physin[0],1,physout[1],1);

	// v = hv/ v
	Vmath::Vdiv(nq,physin[2],1,physin[0],1,physout[2],1);
      }
  }


   void NonlinearSWE::v_ConservativeToPrimitive( )
  {
    int nq = GetTotPoints();

    // u = hu/h
    Vmath::Vdiv(nq,m_fields[1]->GetPhys(),1,m_fields[0]->GetPhys(),1,m_fields[1]->UpdatePhys(),1);

    // v = hv/ v
    Vmath::Vdiv(nq,m_fields[2]->GetPhys(),1,m_fields[0]->GetPhys(),1,m_fields[2]->UpdatePhys(),1);

    // \eta = h - d
    Vmath::Vsub(nq,m_fields[0]->GetPhys(),1,m_depth,1,m_fields[0]->UpdatePhys(),1);
  }

  void NonlinearSWE::PrimitiveToConservative(const Array<OneD, const Array<OneD, NekDouble> >&physin,
					     Array<OneD,       Array<OneD, NekDouble> >&physout)
  {

    int nq = GetTotPoints();

    if(physin.get() == physout.get())
      {
	// copy indata and work with tmp array
	Array<OneD, Array<OneD, NekDouble> >tmp(3);
	for (int i = 0; i < 3; ++i)
	  {
	    // deep copy
	    tmp[i] = Array<OneD, NekDouble>(nq);
	    Vmath::Vcopy(nq,physin[i],1,tmp[i],1);
	  }

	// h = \eta + d
	Vmath::Vadd(nq,tmp[0],1,m_depth,1,physout[0],1);

	// hu = h * u
	Vmath::Vmul(nq,physout[0],1,tmp[1],1,physout[1],1);

	// hv = h * v
	Vmath::Vmul(nq,physout[0],1,tmp[2],1,physout[2],1);

      }
    else
      {
	// h = \eta + d
	Vmath::Vadd(nq,physin[0],1,m_depth,1,physout[0],1);

	// hu = h * u
	Vmath::Vmul(nq,physout[0],1,physin[1],1,physout[1],1);

	// hv = h * v
	Vmath::Vmul(nq,physout[0],1,physin[2],1,physout[2],1);

      }

  }

  void NonlinearSWE::v_PrimitiveToConservative( )
  {
    int nq = GetTotPoints();

    // h = \eta + d
    Vmath::Vadd(nq,m_fields[0]->GetPhys(),1,m_depth,1,m_fields[0]->UpdatePhys(),1);

    // hu = h * u
    Vmath::Vmul(nq,m_fields[0]->GetPhys(),1,m_fields[1]->GetPhys(),1,m_fields[1]->UpdatePhys(),1);

    // hv = h * v
    Vmath::Vmul(nq,m_fields[0]->GetPhys(),1,m_fields[2]->GetPhys(),1,m_fields[2]->UpdatePhys(),1);
  }


 /**
     * @brief Compute the velocity field \f$ \mathbf{v} \f$ given the momentum
     * \f$ h\mathbf{v} \f$.
     *
     * @param physfield  Momentum field.
     * @param velocity   Velocity field.
     */
    void NonlinearSWE::GetVelocityVector(
        const Array<OneD, Array<OneD, NekDouble> > &physfield,
              Array<OneD, Array<OneD, NekDouble> > &velocity)
    {
        const int npts = physfield[0].size();

        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Vdiv(npts, physfield[1+i], 1, physfield[0], 1,
                              velocity[i],    1);
        }
    }


    void NonlinearSWE::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        ShallowWaterSystem::v_GenerateSummary(s);
        SolverUtils::AddSummaryItem(s, "Variables", "h  should be in field[0]");
        SolverUtils::AddSummaryItem(s, "",          "hu should be in field[1]");
        SolverUtils::AddSummaryItem(s, "",          "hv should be in field[2]");
    }

} //end of namespace

