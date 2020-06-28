///////////////////////////////////////////////////////////////////////////////
//
// File ShallowWaterSystem.cpp
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
// Description: Generic timestepping for shallow water solvers
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <ShallowWaterSolver/EquationSystems/ShallowWaterSystem.h>

using namespace std;

namespace Nektar
{
    /**
     * @class ShallowWaterSystem
     *
     * Provides the underlying timestepping framework for shallow water flow solvers
     * including the general timestepping routines. This class is not intended
     * to be directly instantiated, but rather is a base class on which to
     * define shallow water solvers, e.g. SWE, Boussinesq, linear and nonlinear versions.
     *
     * For details on implementing unsteady solvers see
     * \ref sectionADRSolverModuleImplementation
     */

    /**
     * Processes SolverInfo parameters from the session file and sets up
     * timestepping-specific code.
     * @param   pSession        Session object to read parameters from.
     */

   string ShallowWaterSystem::className =
        SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
            "ShallowWaterSystem",
            ShallowWaterSystem::create,
            "Auxiliary functions for the shallow water system.");


    ShallowWaterSystem::ShallowWaterSystem(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph)
    {
    }

    void ShallowWaterSystem::v_InitObject()
    {
        UnsteadySystem::v_InitObject();

       	// if discontinuous Galerkin determine numerical flux to use
	if (m_projectionType == MultiRegions::eDiscontinuous)
	  {
	    ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
		     "No UPWINDTYPE defined in session.");
	  }

	 // Set up locations of velocity vector.
	 m_vecLocs = Array<OneD, Array<OneD, NekDouble> >(1);
	 m_vecLocs[0] = Array<OneD, NekDouble>(m_spacedim);
	 for (int i = 0; i < m_spacedim; ++i)
	 {
             m_vecLocs[0][i] = 1 + i;
	 }

 	// Load generic input parameters
        m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);

	// Load acceleration of gravity
	m_session->LoadParameter("Gravity", m_g, 9.81);

	// input/output in primitive variables
	m_primitive = true;

       	EvaluateWaterDepth();

	m_constantDepth = true;
	NekDouble depth = m_depth[0];
	for (int i = 0; i < GetTotPoints(); ++i)
        {
	    if (m_depth[i] != depth)
            {
		m_constantDepth = false;
		break;
            }
        }
        
	// Compute the bottom slopes
	int nq = GetTotPoints();
	if (m_constantDepth != true)
        {
	    m_bottomSlope = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
	    for (int i = 0; i < m_spacedim; ++i)
            {
		m_bottomSlope[i] = Array<OneD, NekDouble>(nq);
		m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i],
                                       m_depth,m_bottomSlope[i]);
		Vmath::Neg(nq,m_bottomSlope[i],1);
            }
        }
        
	EvaluateCoriolis();
    }


    /**
     *
     */
    ShallowWaterSystem::~ShallowWaterSystem()
    {
    }


    void ShallowWaterSystem::v_GenerateSummary(SolverUtils::SummaryList& s)
    {
        UnsteadySystem::v_GenerateSummary(s);
	if (m_constantDepth == true)
	  {
	    SolverUtils::AddSummaryItem(s, "Depth", "constant");
	  }
	else
	  {
	    SolverUtils::AddSummaryItem(s, "Depth", "variable");
	  }


    }

   void ShallowWaterSystem::v_PrimitiveToConservative()
  {
    ASSERTL0(false, "This function is not implemented for this equation.");
  }

  void ShallowWaterSystem::v_ConservativeToPrimitive()
  {
    ASSERTL0(false, "This function is not implemented for this equation.");
  }

  void ShallowWaterSystem::EvaluateWaterDepth(void)
  {
    GetFunction("WaterDepth")->Evaluate("d", m_depth);
  }


  void ShallowWaterSystem::EvaluateCoriolis(void)
  {
    GetFunction("Coriolis")->Evaluate("f", m_coriolis);
  }

  void ShallowWaterSystem::CopyBoundaryTrace(const Array<OneD, NekDouble>&Fwd, Array<OneD, NekDouble>&Bwd)
  {

    int cnt = 0;
    // loop over Boundary Regions
    for(int bcRegion = 0; bcRegion < m_fields[0]->GetBndConditions().size(); ++bcRegion)
      {
        if (m_fields[0]->GetBndConditions()[bcRegion]->GetBoundaryConditionType()
            == SpatialDomains::ePeriodic)
        {
            continue;
        }

	// Copy the forward trace of the field to the backward trace
        int e, id2, npts;

        for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]
                 ->GetExpSize(); ++e)
        {
            npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
                GetExp(e)->GetTotPoints();
            id2  = m_fields[0]->GetTrace()->GetPhys_Offset(
                        m_fields[0]->GetTraceMap()->
                                    GetBndCondIDToGlobalTraceID(cnt+e));

	    Vmath::Vcopy(npts, &Fwd[id2], 1, &Bwd[id2], 1);
	}

	cnt +=m_fields[0]->GetBndCondExpansions()[bcRegion]->GetExpSize();
      }
  }

}
