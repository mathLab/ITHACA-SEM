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

#include <Auxiliary/ADRBase.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <string>

namespace Nektar
{
    /**
     * @class ADRBase
     * It is basically a class handling vector valued fields
     */


    /**
     * Basic construnctor
     */
    ADRBase::ADRBase(void):
        m_fields(0)
    {
    }


    /**
     * Sets up the base class. Creates a MeshGraph and loads the geometry,
     * boundary conditions, problem dimension and session name from the session
     * file.
     * @param   fileNameString                  Session file name.
     * @param   UseInputFileForProjectionType   Default: false.
     * @param   UseContinuousField              Default: false.
     * @param   globoptfile                     Global optimisations.
     */
    ADRBase::ADRBase(const string &fileNameString, bool UseInputFileForProjectionType,
                     bool UseContinuousField, string &globoptfile)
    {
        SpatialDomains::MeshGraph graph;
        m_filename = fileNameString;

        // Read the geometry and the expansion information
        m_graph = graph.Read(m_filename);

        // Also read and store the boundary conditions
        SpatialDomains::MeshGraph *meshptr = m_graph.get();
        m_boundaryConditions = MemoryManager<SpatialDomains::BoundaryConditions>
                                        ::AllocateSharedPtr(meshptr);
        m_boundaryConditions->Read(m_filename);

        // Read and store history point data
        m_historyPoints = MemoryManager<SpatialDomains::History>
                                        ::AllocateSharedPtr(meshptr);
        m_historyPoints->Read(m_filename);

        // Set space dimension for use in class
        m_spacedim = m_graph->GetSpaceDimension();

        // Save the basename of input file name for output details.
        m_sessionName = fileNameString;
        m_sessionName = m_sessionName.substr(0,
                                m_sessionName.find_last_of("."));

        // Options to determine type of projection from file or
        // directly from constructor
        if(UseInputFileForProjectionType == true)
        {
            if(m_boundaryConditions->SolverInfoExists("PROJECTION"))
            {
                std::string ProjectStr
                        = m_boundaryConditions->GetSolverInfo("PROJECTION");

                if((ProjectStr == "Continuous")||(ProjectStr == "Galerkin")||
                   (ProjectStr == "CONTINUOUS")||(ProjectStr == "GALERKIN"))
                {
                    m_projectionType = eGalerkin;
                }
                else if(ProjectStr == "DisContinuous")
                {
                    m_projectionType = eDiscontinuousGalerkin;
                }
                else
                {
                    ASSERTL0(false,"PROJECTION value not recognised");
                }
            }
            else
            {
                cerr << "Projection type not specified in SOLVERINFO,"
                        "defaulting to continuous Galerkin" << endl;
                m_projectionType = eGalerkin;
            }
        }
        else
        {
            if(UseContinuousField == true)
            {
                m_projectionType = eGalerkin;
            }
            else
            {
                m_projectionType = eDiscontinuousGalerkin;
            }
        }

        SetADRBase(m_graph,m_boundaryConditions->GetNumVariables(),globoptfile);
    }


    /**
     * Set up the dependent variable fields using MultiRegions::ContField1D,
     * MultiRegions::ContField2D, MultiRegions::DisContField1D or
     * MultiRegions::DisContField2D classes. For 2D set up the normals.
     * Also set default values for the total time, time step, number of steps
     * and the frequency of checkpoints.
     * @param   mesh            Session
     * @param   nvariables      Number of dependent variables.
     * @param   globoptfile     Global optimisations.
     */
    void ADRBase::SetADRBase(SpatialDomains::MeshGraphSharedPtr &mesh,
                             int nvariables,  string &globoptfile)
    {
        int i;

        m_fields   = Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);
        m_spacedim = mesh->GetSpaceDimension();
        m_expdim   = mesh->GetMeshDimension();

        // Continuous Galerkin projection
        if(m_projectionType == eGalerkin)
        {
            switch(m_expdim)
            {
                case 1:
                {
                    SpatialDomains::MeshGraph1DSharedPtr mesh1D;

                    if( !(mesh1D = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph1D>(mesh)) )
                    {
                        ASSERTL0(false,"Dynamics cast failed");
                    }

                    for(i = 0 ; i < m_fields.num_elements(); i++)
                    {
                        m_fields[i] = MemoryManager<MultiRegions::ContField1D>
                                            ::AllocateSharedPtr(*mesh1D,
                                                *m_boundaryConditions,i);
                    }
                    break;
                }
                case 2:
                {
                    SpatialDomains::MeshGraph2DSharedPtr mesh2D;

                    if(!(mesh2D = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph2D>(mesh)))
                    {
                        ASSERTL0(false,"Dynamics cast failed");
                    }

                    i = 0;
                    MultiRegions::ContField2DSharedPtr firstfield =
                            MemoryManager<MultiRegions::ContField2D>
                                    ::AllocateSharedPtr(*mesh2D,
                                        *m_boundaryConditions,i);

                    firstfield->ReadGlobalOptimizationParameters(m_filename);

                    m_fields[0] = firstfield;
                    for(i = 1 ; i < m_fields.num_elements(); i++)
                    {
                        if(mesh2D->SameExpansions(m_boundaryConditions->GetVariable(0),
                                                 m_boundaryConditions->GetVariable(i)))
                        {
                            m_fields[i] = MemoryManager<MultiRegions::ContField2D>
                                ::AllocateSharedPtr(*firstfield,
                                                    *mesh2D,*m_boundaryConditions,i);
                        }
                        else
                        {
                            m_fields[i] = MemoryManager<MultiRegions::ContField2D>
                                ::AllocateSharedPtr(*mesh2D,*m_boundaryConditions,i);
                        }
                    }
                    break;
                }
                case 3:
                {
                    SpatialDomains::MeshGraph3DSharedPtr mesh3D;

                    if(!(mesh3D = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph3D>(mesh)))
                    {
                        ASSERTL0(false,"Dynamics cast failed");
                    }

                    i = 0;
                    MultiRegions::ContField3DSharedPtr firstfield =
                            MemoryManager<MultiRegions::ContField3D>
                                    ::AllocateSharedPtr(*mesh3D,
                                        *m_boundaryConditions,i);

                    firstfield->ReadGlobalOptimizationParameters(m_filename);

                    m_fields[0] = firstfield;
                    for(i = 1 ; i < m_fields.num_elements(); i++)
                    {
                        m_fields[i] = MemoryManager<MultiRegions::ContField3D>
                                        ::AllocateSharedPtr(*firstfield,
                                            *mesh3D,*m_boundaryConditions,i);
                    }
                    break;
                }
                default:
                    ASSERTL0(false,"Expansion dimension not recognised");
                    break;
            }
        }
        else // Discontinuous Field
        {
            switch(m_expdim)
            {
                case 1:
                {
                    SpatialDomains::MeshGraph1DSharedPtr mesh1D;

                    if(!(mesh1D = boost::dynamic_pointer_cast<SpatialDomains
                                        ::MeshGraph1D>(mesh)))
                    {
                        ASSERTL0(false,"Dynamics cast failed");
                    }

                    for(i = 0 ; i < m_fields.num_elements(); i++)
                    {
                        m_fields[i] = MemoryManager<MultiRegions
                                ::DisContField1D>::AllocateSharedPtr(*mesh1D,
                                                *m_boundaryConditions,i);
                    }
                    break;
                }
                case 2:
                {
                    SpatialDomains::MeshGraph2DSharedPtr mesh2D;

                    if(!(mesh2D = boost::dynamic_pointer_cast<SpatialDomains
                                        ::MeshGraph2D>(mesh)))
                    {
                        ASSERTL0(false,"Dynamics cast failed");
                    }

                    for(i = 0 ; i < m_fields.num_elements(); i++)
                    {
                        m_fields[i] = MemoryManager<MultiRegions
                                ::DisContField2D>::AllocateSharedPtr(*mesh2D,
                                                *m_boundaryConditions,i);
                    }
                    break;
                }
                case 3:
                    ASSERTL0(false,"3 D not set up");
                default:
                    ASSERTL0(false,"Expansion dimension not recognised");
                    break;
            }

            // Set up Normals.
            switch(m_expdim)
            {
                case 1:
                    // no need??...
                    break;
                case 2:
                {
                    m_traceNormals
                            = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);

                    for(i = 0; i < m_spacedim; ++i)
                    {
                        m_traceNormals[i] = Array<OneD, NekDouble> (m_fields[0]
                                                    ->GetTrace()->GetNpoints());
                    }

                    m_fields[0]->GetTrace()->GetNormals(m_traceNormals);
                    break;
                }
                case 3:
                    ASSERTL0(false,"3 D not set up");
                    break;
                default:
                    ASSERTL0(false,"Expansion dimension not recognised");
                    break;
            }
        }
        
        NekOptimize::LoadElementalOptimizationParameters(m_filename);

        // Set Default Parameter
        if(m_boundaryConditions->CheckForParameter("Time") == true)
        {
            m_time = m_boundaryConditions->GetParameter("Time");
        }
        else
        {
            m_time = 0.0;
        }

        if(m_boundaryConditions->CheckForParameter("TimeStep") == true)
        {
            m_timestep = m_boundaryConditions->GetParameter("TimeStep");
        }

        if(m_boundaryConditions->CheckForParameter("NumSteps") == true)
        {
            m_steps = m_boundaryConditions->GetParameter("NumSteps");
        }
        else
        {
            m_steps  = 0;
        }

        if(m_boundaryConditions->CheckForParameter("IO_CheckSteps") == true)
        {
            m_checksteps = m_boundaryConditions->GetParameter("IO_CheckSteps");
        }
        else
        {
            m_checksteps = m_steps;
        }

        if(m_boundaryConditions->CheckForParameter("FinTime") == true)
        {
            m_fintime = m_boundaryConditions->GetParameter("FinTime");
        }
        else
        {
            m_fintime = 0;
        }
		
		if(m_boundaryConditions->CheckForParameter("NumQuadPointsError") == true)
        {
             m_NumQuadPointsError = m_boundaryConditions->GetParameter("NumQuadPointsError");
        }
        else
        {
             m_NumQuadPointsError = 0;
        }

        // Read in spatial data
        int nq = m_fields[0]->GetNpoints();
        m_spatialParameters = MemoryManager<SpatialDomains::SpatialParameters>
          ::AllocateSharedPtr(nq);
        m_spatialParameters->Read(m_filename);

        Array<OneD, NekDouble> x(nq), y(nq), z(nq);
        m_fields[0]->GetCoords(x,y,z);
        m_spatialParameters->EvaluateParameters(x,y,z);

        ScanForHistoryPoints();
    }


    /**
     * Zero the physical fields.
     */
    void ADRBase::ZeroPhysFields(void)
    {
        for(int i = 0; i < m_fields.num_elements(); i++)
        {
            Vmath::Zero(m_fields[i]->GetNpoints(),m_fields[i]->UpdatePhys(),1);
        }
    }

    /**
     * FwdTrans the m_fields members
     */
    void ADRBase::FwdTransFields(void)
    {
        for(int i = 0; i < m_fields.num_elements(); i++)
        {
            m_fields[i]->FwdTrans(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
            m_fields[i]->SetPhysState(false);
        }
    }

    /**
     * Set the physical fields based on a restart file, or a function
     * describing the initial condition given in the session.
     * @param   initialtime             Time at which to evaluate the function.
     * @param   dumpInitialConditions   Write the initial condition to file?
     */
     void ADRBase::v_SetInitialConditions(NekDouble initialtime,
                                          bool dumpInitialConditions)
    {
        std::string restartstr = "RESTART";

        cout << "Initial Conditions:" << endl;

        // Check for restart file.
        if(m_boundaryConditions->FoundInitialCondition(restartstr))
        {
            SpatialDomains::ConstInitialConditionShPtr ifunc
                    = m_boundaryConditions->GetInitialCondition(restartstr);

            std::string restartfile = ifunc->GetEquation();
            cout << "\tRestart file: "<< restartfile << endl;
            ImportFld(restartfile);
        }
        else
        {
            int nq = m_fields[0]->GetNpoints();

            Array<OneD,NekDouble> x0(nq);
            Array<OneD,NekDouble> x1(nq);
            Array<OneD,NekDouble> x2(nq);

            // get the coordinates (assuming all fields have the same
            // discretisation)
            m_fields[0]->GetCoords(x0,x1,x2);

            for(int i = 0 ; i < m_fields.num_elements(); i++)
            {
                SpatialDomains::ConstInitialConditionShPtr ifunc
                        = m_boundaryConditions->GetInitialCondition(i);
                for(int j = 0; j < nq; j++)
                {
                    (m_fields[i]->UpdatePhys())[j]
                            = ifunc->Evaluate(x0[j],x1[j],x2[j],initialtime);
                }
                m_fields[i]->SetPhysState(true);
                m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->GetPhys(),
                                                 m_fields[i]->UpdateCoeffs());
                cout << "\tField "<< m_boundaryConditions->GetVariable(i)
                     <<": " << ifunc->GetEquation() << endl;
            }
        }
        if(dumpInitialConditions)
        {
            // dump initial conditions to file
            std::string outname = m_sessionName + "_initial.chk";
            WriteFld(outname);
        }
    }


    /**
     * Populates a forcing function for each of the dependent variables using
     * the expression provided by the BoundaryConditions object.
     * @param   force           Array of fields to assign forcing.
     */
    void ADRBase::SetPhysForcingFunctions(
                        Array<OneD, MultiRegions::ExpListSharedPtr> &force)
    {
        int nq = m_fields[0]->GetNpoints();

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        // get the coordinates (assuming all fields have the same
        // discretisation)
        force[0]->GetCoords(x0,x1,x2);

        for(int i = 0 ; i < m_fields.num_elements(); i++)
        {
            SpatialDomains::ConstForcingFunctionShPtr ffunc
                    = m_boundaryConditions->GetForcingFunction(i);

            for(int j = 0; j < nq; j++)
            {
                (force[i]->UpdatePhys())[j]
                        = ffunc->Evaluate(x0[j],x1[j],x2[j]);
            }
            force[i]->SetPhysState(true);
        }
    }

    /**
     * Evaluates the exact solution provided in the session for a given
     * dependent variable.
     * @param   field               The index of the field to evaluate.
     * @param   outfield            Storage for exact solution.
     * @param   time                The time at which to evaluate the solution.
     */
  void ADRBase::EvaluateExactSolution(int field, Array<OneD, NekDouble> &exactsoln, const NekDouble time)
    { 
      v_EvaluateExactSolution(field,exactsoln,time);
    }

    /**
     * Evaluates a user-defined expression for all dependent variables.
     * @param   outfield            Array of fields to assign values.
     */
    void ADRBase::EvaluateUserDefinedEqn(
                            Array<OneD, Array<OneD, NekDouble> > &outfield)
    {
        int nq = m_fields[0]->GetNpoints();

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        // get the coordinates (assuming all fields have the same
        // discretisation)
        m_fields[0]->GetCoords(x0,x1,x2);

        for(int i = 0 ; i < m_fields.num_elements(); i++)
        {
            SpatialDomains::ConstUserDefinedEqnShPtr ifunc
                            = m_boundaryConditions->GetUserDefinedEqn(i);
            for(int j = 0; j < nq; j++)
            {
                outfield[i][j] = ifunc->Evaluate(x0[j],x1[j],x2[j]);
            }
        }
    }


    /**
     * Compute the error in the L2-norm
     * @param   field           The field to compare.
     * @param   exactsoln       The exact solution to compare with.
     * @param   Normalised      Normalise L2-error.
     * @returns                 Error in the L2-norm.
     */
    NekDouble ADRBase::L2Error(int field,
                               const Array<OneD, NekDouble> &exactsoln,
                               bool Normalised)
    {
		NekDouble L2error = -1.0;
		
		if(m_NumQuadPointsError == 0)
		{
			if(m_fields[field]->GetPhysState() == false)
			{
				m_fields[field]->BwdTrans(m_fields[field]->GetCoeffs(),
                                       m_fields[field]->UpdatePhys());
			}

			if(exactsoln.num_elements())
			{
				L2error = m_fields[field]->L2(exactsoln);
			}
			else
			{
				Array<OneD, NekDouble> exactsoln(m_fields[field]->GetNpoints());

				EvaluateExactSolution(field,exactsoln,m_time);

				L2error = m_fields[field]->L2(exactsoln);
			}

			if(Normalised == true)
			{
				Array<OneD, NekDouble> one(m_fields[field]->GetNpoints(),1.0);

				NekDouble Vol = m_fields[field]->PhysIntegral(one);

				L2error = sqrt(L2error*L2error/Vol);
			}
		}
		else 
		{
			Array<OneD,NekDouble> L2INF(2);
			L2INF = ErrorExtraPoints(field);
			L2error = L2INF[0];
			
		}


        return L2error;
    }


    /**
     * Compute the error in the L_inf-norm
     * @param   field           The field to compare.
     * @param   exactsoln       The exact solution to compare with.
     * @returns                 Error in the L_inft-norm.
     */
    NekDouble ADRBase::LinfError(int field,
                                 const Array<OneD, NekDouble> &exactsoln)
    {
		NekDouble Linferror = -1.0;
		
		if(m_NumQuadPointsError == 0)
		{
			if(m_fields[field]->GetPhysState() == false)
			{
				m_fields[field]->BwdTrans(m_fields[field]->GetCoeffs(),
                                      m_fields[field]->UpdatePhys());
			}

			if(exactsoln.num_elements())
			{
				Linferror = m_fields[field]->Linf(exactsoln);
			}
			else
			{
				Array<OneD, NekDouble> exactsoln(m_fields[field]->GetNpoints());

				EvaluateExactSolution(field,exactsoln,m_time);

				Linferror = m_fields[field]->Linf(exactsoln);
			}
		}
		else 
		{
			Array<OneD,NekDouble> L2INF(2);
			L2INF = ErrorExtraPoints(field);
			Linferror = L2INF[1];
		}
		
		return Linferror;
    }
	
	/**
     * Compute the error in the L2-norm, L-inf for a larger number of Quadrature Points
     * @param   field              The field to compare.
     * @returns                    Error in the L2-norm and L-inf norm. 
     */
    Array<OneD,NekDouble> ADRBase::ErrorExtraPoints(int field)
    {
		SpatialDomains::MeshGraph2DSharedPtr mesh2D;
		mesh2D = boost::dynamic_pointer_cast<SpatialDomains::MeshGraph2D>(m_graph);
		
		int NumModes = GetNumExpModes();
		
		Array<OneD,NekDouble> L2INF(2);
		
		const LibUtilities::PointsKey PkeyT1(m_NumQuadPointsError,LibUtilities::eGaussLobattoLegendre);
		const LibUtilities::PointsKey PkeyT2(m_NumQuadPointsError,LibUtilities::eGaussRadauMAlpha1Beta0);
		const LibUtilities::PointsKey PkeyQ1(m_NumQuadPointsError,LibUtilities::eGaussLobattoLegendre);
		const LibUtilities::PointsKey PkeyQ2(m_NumQuadPointsError,LibUtilities::eGaussLobattoLegendre);
		const LibUtilities::BasisKey  BkeyT1(LibUtilities::eModified_A,NumModes,PkeyT1);
		const LibUtilities::BasisKey  BkeyT2(LibUtilities::eModified_B,NumModes,PkeyT2);
		const LibUtilities::BasisKey  BkeyQ1(LibUtilities::eModified_A,NumModes,PkeyQ1);
		const LibUtilities::BasisKey  BkeyQ2(LibUtilities::eModified_A,NumModes,PkeyQ2);
        
		MultiRegions::ExpList2DSharedPtr ErrorExp = 
        MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(BkeyT1,BkeyT2,BkeyQ1,BkeyQ2,*mesh2D);
		
		int ErrorCoordim = ErrorExp->GetCoordim(0);
		int ErrorNq      = ErrorExp->GetTotPoints();
		
		Array<OneD,NekDouble> ErrorXc0(ErrorNq,0.0);
		Array<OneD,NekDouble> ErrorXc1(ErrorNq,0.0);
		Array<OneD,NekDouble> ErrorXc2(ErrorNq,0.0);
		
		switch(ErrorCoordim)
		{
			case 1:
				ErrorExp->GetCoords(ErrorXc0);
				break;
			case 2:
				ErrorExp->GetCoords(ErrorXc0,ErrorXc1);
				break;
			case 3:
				ErrorExp->GetCoords(ErrorXc0,ErrorXc1,ErrorXc2);
				break;
		}
		
		SpatialDomains::ConstExactSolutionShPtr exSol = m_boundaryConditions->GetExactSolution(field);
		// evaluate exact solution 
		Array<OneD,NekDouble> ErrorSol(ErrorNq);
		for(int i = 0; i < ErrorNq; ++i)
		{
			ErrorSol[i] = exSol->Evaluate(ErrorXc0[i],ErrorXc1[i],ErrorXc2[i],m_time);
		}
		
		// calcualte spectral/hp approximation on the quad points of this new
		// expansion basis
		ErrorExp->BwdTrans_IterPerExp(m_fields[field]->GetCoeffs(),ErrorExp->UpdatePhys());
		
		L2INF[0]    = ErrorExp->L2  (ErrorSol);
		L2INF[1]    = ErrorExp->Linf(ErrorSol);
		
        return L2INF;
    }
	
	


    /**
     * Computes the weak Green form of advection terms (without boundary
     * integral), i.e. \f$ (\nabla \phi \cdot F) \f$ where for example
     * \f$ F=uV \f$.
     * @param   F           Fields.
     * @param   outarray    Storage for result.
     *
     * \note Assuming all fields are of the same expansion and order so that we
     * can use the parameters of m_fields[0].
     */
    void ADRBase::WeakAdvectionGreensDivergenceForm(
                const Array<OneD, Array<OneD, NekDouble> > &F,
                Array<OneD, NekDouble> &outarray)
    {
        // use dimension of Velocity vector to dictate dimension of operation
        int ndim    = F.num_elements();
        int nCoeffs = m_fields[0]->GetNcoeffs();

        Array<OneD, NekDouble> iprod(nCoeffs);
        Vmath::Zero(nCoeffs, outarray, 1);

        for (int i = 0; i < ndim; ++i)
        {
            m_fields[0]->IProductWRTDerivBase(i, F[i], iprod);
            Vmath::Vadd(nCoeffs, iprod, 1, outarray, 1, outarray, 1);
        }
    }


    /**
     * Calculate Inner product of the divergence advection form
     * \f$(\phi, \nabla \cdot F)\f$, where for example \f$ F = uV \f$.
     * @param   F           Fields.
     * @param   outarray    Storage for result.
     */
    void ADRBase::WeakAdvectionDivergenceForm(
                const Array<OneD, Array<OneD, NekDouble> > &F,
                Array<OneD, NekDouble> &outarray)
    {
        // use dimension of Velocity vector to dictate dimension of operation
        int ndim       = F.num_elements();
        int nPointsTot = m_fields[0]->GetNpoints();
        Array<OneD, NekDouble> tmp(nPointsTot);
        Array<OneD, NekDouble> div(nPointsTot, 0.0);

        // Evaluate the divergence
        for(int i = 0; i < ndim; ++i)
        {
            m_fields[0]->PhysDeriv(i,F[i],tmp);
            Vmath::Vadd(nPointsTot, tmp, 1, div, 1, div, 1);
        }

        m_fields[0]->IProductWRTBase(div, outarray);
    }


    /**
     * Calculate Inner product of the divergence advection form
     * \f$ (\phi, V\cdot \nabla u) \f$
     * @param   V           Fields.
     * @param   u           Fields.
     * @param   outarray    Storage for result.
     */
    void ADRBase::WeakAdvectionNonConservativeForm(
                const Array<OneD, Array<OneD, NekDouble> > &V,
                const Array<OneD, const NekDouble> &u,
                Array<OneD, NekDouble> &outarray)
    {
        // use dimension of Velocity vector to dictate dimension of operation
        int ndim       = V.num_elements();

        int nPointsTot = m_fields[0]->GetNpoints();
        Array<OneD, NekDouble> tmp(nPointsTot);
        Array<OneD, NekDouble> wk(ndim * nPointsTot, 0.0);

        AdvectionNonConservativeForm(V, u, tmp, wk);

        m_fields[0]->IProductWRTBase_IterPerExp(tmp, outarray);
    }


    /**
     * Calculate the inner product \f$ V\cdot \nabla u \f$
     * @param   V           Fields.
     * @param   u           Fields.
     * @param   outarray    Storage for result.
     * @param   wk          Workspace.
     */
    void ADRBase::AdvectionNonConservativeForm(
                const Array<OneD, Array<OneD, NekDouble> > &V,
                const Array<OneD, const NekDouble> &u,
                Array<OneD, NekDouble> &outarray,
                Array<OneD, NekDouble> &wk)
    {
        // use dimension of Velocity vector to dictate dimension of operation
        int ndim       = V.num_elements();
        //int ndim = m_expdim;

        // ToDo: here we should add a check that V has right dimension

        int nPointsTot = m_fields[0]->GetNpoints();
        Array<OneD, NekDouble> grad0,grad1,grad2;

        // check to see if wk space is defined
        if (wk.num_elements())
        {
            grad0 = wk;
        }
        else
        {
            grad0 = Array<OneD, NekDouble> (ndim*nPointsTot);
        }

        // Evaluate V\cdot Grad(u)
        switch(ndim)
        {
        case 1:
            m_fields[0]->PhysDeriv(u,grad0);
            Vmath::Vmul(nPointsTot,grad0,1,V[0],1,outarray,1);
            break;
        case 2:
            grad1 = grad0 + nPointsTot;
            m_fields[0]->PhysDeriv(u,grad0,grad1);
            Vmath::Vmul (nPointsTot,grad0,1,V[0],1,outarray,1);
            Vmath::Vvtvp(nPointsTot,grad1,1,V[1],1,outarray,1,outarray,1);
            break;
        case 3:
            grad1 = grad0 + nPointsTot;
            grad2 = grad1 + nPointsTot;
            m_fields[0]->PhysDeriv(u,grad0,grad1,grad2);
            Vmath::Vmul (nPointsTot,grad0,1,V[0],1,outarray,1);
            Vmath::Vvtvp(nPointsTot,grad1,1,V[1],1,outarray,1,outarray,1);
            Vmath::Vvtvp(nPointsTot,grad2,1,V[2],1,outarray,1,outarray,1);
            break;

        default:
            ASSERTL0(false,"dimension unknown");
        }
    }

    /*
     * Calculate weak DG advection in the form
     * \f$ \langle\phi, \hat{F}\cdot n\rangle - (\nabla \phi \cdot F) \f$
     * @param   InField         Fields.
     * @param   OutField        Storage for result.
     * @param   NumericalFluxIncludesNormal     Default: true.
     * @param   InFieldIsPhysSpace              Default: false.
     * @param   nvariables      Number of fields.
     */
    void ADRBase::WeakDGAdvection(
                const Array<OneD, Array<OneD, NekDouble> >& InField,
                Array<OneD, Array<OneD, NekDouble> >& OutField,
                bool NumericalFluxIncludesNormal,
                bool InFieldIsInPhysSpace,
                int nvariables)
    {
        int i;
        int nVelDim         = m_spacedim;
        int nPointsTot      = GetNpoints();
        int ncoeffs         = GetNcoeffs();
        int nTracePointsTot = GetTraceNpoints();

        if (!nvariables)
        {
            nvariables      = m_fields.num_elements();
        }

        Array<OneD, Array<OneD, NekDouble> > fluxvector(nVelDim);
        Array<OneD, Array<OneD, NekDouble> > physfield (nvariables);

        for(i = 0; i < nVelDim; ++i)
        {
            fluxvector[i]    = Array<OneD, NekDouble>(nPointsTot);
        }

        // Get the variables in physical space
        // already in physical space
        if(InFieldIsInPhysSpace == true)
        {
            for(i = 0; i < nvariables; ++i)
            {
                physfield[i] = InField[i];
            }
        }
        // otherwise do a backward transformation
        else
        {
            for(i = 0; i < nvariables; ++i)
            {
                // Could make this point to m_fields[i]->UpdatePhys();
                physfield[i] = Array<OneD, NekDouble>(nPointsTot);
                m_fields[i]->BwdTrans(InField[i],physfield[i]);
            }
        }

        // Get the advection part (without numerical flux)
        for(i = 0; i < nvariables; ++i)
        {
            // Get the ith component of the  flux vector in (physical space)
            GetFluxVector(i, physfield, fluxvector);

            // Calculate the i^th value of (\grad_i \phi, F)
            WeakAdvectionGreensDivergenceForm(fluxvector,OutField[i]);
        }

        // Get the numerical flux and add to the modal coeffs
        // if the NumericalFluxs function already includes the
        // normal in the output
        if (NumericalFluxIncludesNormal == true)
        {
            Array<OneD, Array<OneD, NekDouble> > numflux   (nvariables);

            for(i = 0; i < nvariables; ++i)
            {
                numflux[i]   = Array<OneD, NekDouble>(nTracePointsTot);
            }

            // Evaluate numerical flux in physical space which may in
            // general couple all component of vectors
            NumericalFlux(physfield, numflux);

            // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
            for(i = 0; i < nvariables; ++i)
            {
                Vmath::Neg(ncoeffs,OutField[i],1);
                m_fields[i]->AddTraceIntegral(numflux[i],OutField[i]);
                m_fields[i]->SetPhysState(false);
            }
        }
        // if the NumericalFlux function does not include the
        // normal in the output
        else
        {
            Array<OneD, Array<OneD, NekDouble> > numfluxX   (nvariables);
            Array<OneD, Array<OneD, NekDouble> > numfluxY   (nvariables);

            for(i = 0; i < nvariables; ++i)
            {
                numfluxX[i]   = Array<OneD, NekDouble>(nTracePointsTot);
                numfluxY[i]   = Array<OneD, NekDouble>(nTracePointsTot);
            }

            // Evaluate numerical flux in physical space which may in
            // general couple all component of vectors
            NumericalFlux(physfield, numfluxX, numfluxY);

            // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
            for(i = 0; i < nvariables; ++i)
            {
                Vmath::Neg(ncoeffs,OutField[i],1);
                m_fields[i]->AddTraceIntegral(numfluxX[i], numfluxY[i],
                                              OutField[i]);
                m_fields[i]->SetPhysState(false);
            }
        }
    }


    /**
     * Calculate weak DG Diffusion in the LDG form
     * \f$ \langle\psi, \hat{u}\cdot n\rangle
     * - \langle\nabla\psi \cdot u\rangle
     *  \langle\phi, \hat{q}\cdot n\rangle - (\nabla \phi \cdot q) \rangle \f$
     */
    void ADRBase::WeakDGDiffusion(
            const Array<OneD, Array<OneD, NekDouble> >& InField,
            Array<OneD, Array<OneD, NekDouble> >& OutField,
            bool NumericalFluxIncludesNormal,
            bool InFieldIsInPhysSpace)
    {
        int i,j,k;
        int nPointsTot      = GetNpoints();
        int ncoeffs         = GetNcoeffs();
        int nTracePointsTot = GetTraceNpoints();
        int nvariables = m_fields.num_elements();
        int nqvar = 2;

        Array<OneD, NekDouble>  qcoeffs (ncoeffs);
        Array<OneD, NekDouble>  temp (ncoeffs);

        Array<OneD, Array<OneD, NekDouble> > fluxvector (m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > ufield (nvariables);

        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  flux   (nqvar);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > >  qfield  (nqvar);

        for(j = 0; j < nqvar; ++j)
        {
            qfield[j] = Array<OneD, Array<OneD, NekDouble> >(nqvar);
            flux[j] = Array<OneD, Array<OneD, NekDouble> >(nqvar);

            for(i = 0; i< nvariables; ++i)
            {
                ufield[i] = Array<OneD, NekDouble>(nPointsTot,0.0);

                qfield[j][i]  = Array<OneD, NekDouble>(nPointsTot,0.0);
                flux[j][i] = Array<OneD, NekDouble>(nTracePointsTot,0.0);
            }
        }

        for(k = 0; k < m_spacedim; ++k)
        {
            fluxvector[k] = Array<OneD, NekDouble>(nPointsTot,0.0);
        }

        // Get the variables in physical space
        // already in physical space
        if(InFieldIsInPhysSpace == true)
        {
            for(i = 0; i < nvariables; ++i)
            {
                ufield[i] = InField[i];
            }
        }
        // otherwise do a backward transformation
        else
        {
            for(i = 0; i < nvariables; ++i)
            {
                // Could make this point to m_fields[i]->UpdatePhys();
                ufield[i] = Array<OneD, NekDouble>(nPointsTot);
                m_fields[i]->BwdTrans(InField[i],ufield[i]);
            }
        }

        // ##########################################################
        // Compute q_{\eta} and q_{\xi} from su
        // Obtain Numerical Fluxes
        // ##########################################################
        NumFluxforScalar(ufield, flux);

        for(j = 0; j < nqvar; ++j)
        {
            for(i = 0; i < nvariables; ++i)
            {
                // Get the ith component of the  flux vector in (physical space)
                // fluxvector = m_tanbasis * u, where m_tanbasis = 2 by
                // m_spacedim by nPointsTot
                if(m_tanbasis.num_elements())
                {
                    for (k = 0; k < m_spacedim; ++k)
                    {
                        Vmath::Vmul(nPointsTot, m_tanbasis[j][k], 1, ufield[i],
                                    1, fluxvector[k], 1);
                    }
                }
                else
                {
                    GetFluxVector(i, j, ufield, fluxvector);
                }

                // Calculate the i^th value of (\grad_i \phi, F)
                WeakAdvectionGreensDivergenceForm(fluxvector, qcoeffs);

                Vmath::Neg(ncoeffs,qcoeffs,1);
                m_fields[i]->AddTraceIntegral(flux[j][i], qcoeffs);
                m_fields[i]->SetPhysState(false);

                // Add weighted mass matrix = M ( \nabla \cdot Tanbasis )
                if(m_gradtan.num_elements())
                {
                    MultiRegions::GlobalMatrixKey key(StdRegions::eMass,
                                                        m_gradtan[j]);
                    m_fields[i]->MultiRegions::ExpList::GeneralMatrixOp(key,
                                                        InField[i], temp);
                    Vmath::Svtvp(ncoeffs, -1.0, temp, 1, qcoeffs, 1,
                                                        qcoeffs, 1);
                }

                 //Multiply by the inverse of mass matrix
                m_fields[i]->MultiplyByElmtInvMass(qcoeffs, qcoeffs);

                // Back to physical space
                m_fields[i]->BwdTrans(qcoeffs, qfield[j][i]);
            }
        }


        // ##########################################################
        //   Compute u from q_{\eta} and q_{\xi}
        // ##########################################################

        // Obtain Numerical Fluxes
        NumFluxforVector(ufield, qfield, flux[0]);

        for (i = 0; i < nvariables; ++i)
        {
            // L = L(tan_eta) q_eta + L(tan_xi) q_xi
            OutField[i] = Array<OneD, NekDouble>(ncoeffs, 0.0);
            temp = Array<OneD, NekDouble>(ncoeffs, 0.0);

            if(m_tanbasis.num_elements())
            {
                for(j = 0; j < nqvar; ++j)
                {
                    for (k = 0; k < m_spacedim; ++k)
                    {
                        Vmath::Vmul(nPointsTot, m_tanbasis[j][k], 1,
                                    qfield[j][i], 1, fluxvector[k], 1);
                    }

                    WeakAdvectionGreensDivergenceForm(fluxvector, temp);
                    Vmath::Vadd(ncoeffs, temp, 1, OutField[i], 1,
                                                    OutField[i], 1);
                }
            }
            else
            {
                for (k = 0; k < m_spacedim; ++k)
                {
                    Vmath::Vcopy(nPointsTot, qfield[k][i], 1, fluxvector[k], 1);
                }

                WeakAdvectionGreensDivergenceForm(fluxvector, OutField[i]);
            }

            // Evaulate  <\phi, \hat{F}\cdot n> - OutField[i]
            Vmath::Neg(ncoeffs,OutField[i],1);
            m_fields[i]->AddTraceIntegral(flux[0][i], OutField[i]);
            m_fields[i]->SetPhysState(false);
        }
    }


    /**
     * Write the field data to file. The file is named according to the session
     * name with the extension .fld appended.
     */
    void ADRBase::Output()
    {
        std::string outname = m_sessionName + ".fld";
        WriteFld(outname); 
    }


    /**
     * Write the field data to file. The file is named according to the session
     * name with the extension .fld appended.
     */
    void ADRBase::Output(MultiRegions::ExpListSharedPtr &field, Array< OneD, Array<OneD, NekDouble> > &fieldcoeffs, Array<OneD, std::string> &variables)
    {
        std::string outname = m_sessionName + ".fld";
        WriteFld(outname,field, fieldcoeffs, variables);
    }

    /**
     * Write the n-th checkpoint file.
     * @param   n           The index of the checkpoint file.
     */
    void ADRBase::Checkpoint_Output(const int n)
    {
        char chkout[16] = "";
        sprintf(chkout, "%d", n);
        std::string outname = m_sessionName +"_" + chkout + ".chk";
        WriteFld(outname);
    }

    /**
     * Write the n-th checkpoint file.
     * @param   n           The index of the checkpoint file.
     */
    void ADRBase::Checkpoint_Output(const int n, MultiRegions::ExpListSharedPtr &field, Array< OneD, Array<OneD, NekDouble> > &fieldcoeffs, Array<OneD, std::string> &variables)
    {
        char chkout[16] = "";
        sprintf(chkout, "%d", n);
        std::string outname = m_sessionName +"_" + chkout + ".chk";
        WriteFld(outname, field, fieldcoeffs, variables);
    }


    /**
     * Writes the field data to a file with the given filename.
     * @param   outname     Filename to write to.
     */
    void ADRBase::WriteFld(std::string &outname)
    {
        Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(m_fields.num_elements());
        Array<OneD, std::string>  variables(m_fields.num_elements());
        
        for(int i = 0; i < m_fields.num_elements(); ++i)
        {
            if (m_fields[i]->GetPhysState()==true)
            {	
                m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->GetPhys(),m_fields[i]->UpdateCoeffs());
            }
            fieldcoeffs[i] = m_fields[i]->UpdateCoeffs();
            variables[i] = m_boundaryConditions->GetVariable(i);
        }
        
        WriteFld(outname, m_fields[0], fieldcoeffs, variables);
    }


    /**
     * Writes the field data to a file with the given filename.
     * @param   outname     Filename to write to.
     * @param   field       ExpList on which data is based
     * @param fieldcoeffs   An array of array of expansion coefficients
     * @param  variables    An array of variable names
     */
    void ADRBase::WriteFld(std::string &outname, MultiRegions::ExpListSharedPtr &field, Array<OneD, Array<OneD, NekDouble> > &fieldcoeffs, Array<OneD, std::string> &variables)
    {
        
    	std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
            = field->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        // copy Data into FieldData and set variable
        for(int j = 0; j < fieldcoeffs.num_elements(); ++j)
        {
            for(int i = 0; i < FieldDef.size(); ++i)
            {
                // Could do a search here to find correct variable
                FieldDef[i]->m_fields.push_back(variables[j]);
                field->AppendFieldData(FieldDef[i], FieldData[i], fieldcoeffs[j]);
            }
        }

        m_graph->Write(outname,FieldDef,FieldData);
    }


    /**
     * Import field from infile and load into \a m_fields. This routine will
     * also perform a \a BwdTrans to ensure data is in both the physical and
     * coefficient storage.
     * @param   infile          Filename to read.
     */
    void ADRBase::ImportFld(std::string &infile)
    {
        std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef;
        std::vector<std::vector<NekDouble> > FieldData;

        m_graph->Import(infile,FieldDef,FieldData);

        // copy FieldData into m_fields
        for(int j = 0; j < m_fields.num_elements(); ++j)
        {
            for(int i = 0; i < FieldDef.size(); ++i)
            {
                bool flag = FieldDef[i]->m_fields[j]
                                    == m_boundaryConditions->GetVariable(j);
                ASSERTL1(flag, (std::string("Order of ") + infile
                            + std::string(" data and that defined in "
                                    "m_boundaryconditions differs")).c_str());

                m_fields[j]->ExtractDataToCoeffs(FieldDef[i], FieldData[i],
                                                 FieldDef[i]->m_fields[j]);
            }
            m_fields[j]->BwdTrans(m_fields[j]->GetCoeffs(),
                                  m_fields[j]->UpdatePhys());
        }
    }


    /**
     * Write data to file in Tecplot format?
     * @param   n           Checkpoint index.
     * @param   name        Additional name (appended to session name).
     * @param   inarray     Field data to write out.
     * @param   IsInPhysicalSpace   Indicates if field data is in phys space.
     */
    void ADRBase::Array_Output(const int n, std::string name,
                               const Array<OneD, const NekDouble>&inarray,
                               bool IsInPhysicalSpace)
    {
        int nq = m_fields[0]->GetTotPoints();

        Array<OneD, NekDouble> tmp(nq);

        // save values
        Vmath::Vcopy(nq, m_fields[0]->GetPhys(), 1, tmp, 1);

        // put inarray in m_phys
        if (IsInPhysicalSpace == false)
        {
            m_fields[0]->BwdTrans(inarray,(m_fields[0]->UpdatePhys()));
        }
        else
        {
            Vmath::Vcopy(nq,inarray,1,(m_fields[0]->UpdatePhys()),1);
        }

        char chkout[16] = "";
        sprintf(chkout, "%d", n);
        std::string outname = m_sessionName +"_" + name + "_" + chkout + ".chk";
        ofstream outfile(outname.c_str());
        m_fields[0]->WriteToFile(outfile,eTecplot);

        // copy back the original values
        Vmath::Vcopy(nq,tmp,1,m_fields[0]->UpdatePhys(),1);
    }

  /**
   * Write data to file in Tecplot format
   * @param   n                   Checkpoint index.
   * @param   name                Additional name (appended to session name).
   * @param   IsInPhysicalSpace   Indicates if field data is in phys space.
   */
  void ADRBase::WriteTecplotFile(const int n, std::string name, bool IsInPhysicalSpace)
  {
    int nq = m_fields[0]->GetTotPoints();

    std::string var = "";
    for(int j = 0; j < m_fields.num_elements(); ++j)
      {
	var = var + ", " + m_boundaryConditions->GetVariable(j);
      }

    char chkout[16] = "";
    sprintf(chkout, "%d", n);
    std::string outname = m_sessionName + "_" + name + "_" + chkout + ".dat";
    ofstream outfile(outname.c_str());

    // put inarray in m_phys
    if (IsInPhysicalSpace == false)
      {
	for(int i = 0; i < m_fields.num_elements(); ++i)
	  {
	    m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),m_fields[i]->UpdatePhys());
	  }
      }

    m_fields[0]->WriteTecplotHeader(outfile,var);

    for(int i = 0; i < m_fields[0]->GetExpSize(); ++i)
      {
	m_fields[0]->WriteTecplotZone(outfile,i);
	for(int j = 0; j < m_fields.num_elements(); ++j)
    	  {
    	    m_fields[j]->WriteTecplotField(outfile,i);
    	  }
      }
  }

    /**
     *
     */
    void ADRBase::ScanForHistoryPoints()
    {
        m_historyList.clear();
        Array<OneD, NekDouble> gloCoord(3,0.0);
        for (int i = 0; i < m_historyPoints->GetNumHistoryPoints(); ++i) {
            SpatialDomains::VertexComponentSharedPtr vtx = m_historyPoints->GetHistoryPoint(i);
            vtx->GetCoords(gloCoord[0], gloCoord[1], gloCoord[2]);

            int eId = m_fields[0]->GetExpIndex(gloCoord);

            m_historyList.push_back(
                std::pair<SpatialDomains::VertexComponentSharedPtr, int>(vtx, eId));
        }
    }

    /**
     * @todo Efficiency improvement required. PutPhysInToElmtExp needs to be
     * called only once per field.
     */
    void ADRBase::WriteHistoryData (std::ostream &out)
    {
        int numPoints = m_historyList.size();
        int numFields = m_fields.num_elements();

        vector<NekDouble> data(numPoints*numFields);
        int k;
        Array<OneD, NekDouble> gloCoord(3,0.0);
        std::list<pair<SpatialDomains::VertexComponentSharedPtr, int> >::iterator x;

        // Pull out data values field by field
        for (int j = 0; j < m_fields.num_elements(); ++j)
        {
            m_fields[j]->BwdTrans(m_fields[j]->GetCoeffs(),m_fields[j]->UpdatePhys());
            m_fields[j]->PutPhysInToElmtExp();
            for (k = 0, x = m_historyList.begin(); x != m_historyList.end(); ++x, ++k)
            {
                (*x).first->GetCoords(gloCoord[0], gloCoord[1], gloCoord[2]);
                data[k*numFields+j] = m_fields[j]->GetExp((*x).second)->PhysEvaluate(gloCoord);
            }
        }

        // Write data values point by point
        for (k = 0, x = m_historyList.begin(); x != m_historyList.end(); ++x, ++k)
        {
            (*x).first->GetCoords(gloCoord[0], gloCoord[1], gloCoord[2]);
            out.width(8);
            out << m_time;
            out.width(8);
            out << gloCoord[0];
            out.width(8);
            out << gloCoord[1];
            out.width(8);
            out << gloCoord[2];
            for (int j = 0; j < numFields; ++j)
            {
                //m_fields[j]->PutPhysInToElmtExp();
                out.width(14);
                //out << m_fields[j]->GetExp((*x).second)->PhysEvaluate(gloCoord);
                out << data[k*numFields+j];
            }
            out << endl;
        }
    }


    /**
     * Write out a summary of the session and timestepping to the given output
     * stream.
     * @param   out         Output stream to write data to.
     */
    void ADRBase::Summary(std::ostream &out)
    {
        SessionSummary(out);
        TimeParamSummary(out);
    }


    /**
     * Write out a summary of the session data.
     * @param   out         Output stream to write data to.
     */
    void ADRBase::SessionSummary(std::ostream &out)
    {

        out << "\tSession Name    : " << m_sessionName << endl;
        out << "\tExpansion Dim.  : " << m_expdim << endl;
        out << "\tSpatial   Dim.  : " << m_spacedim << endl;
        out << "\tMax Exp. Order  : " << m_fields[0]->EvalBasisNumModesMax()
                                      << endl;
        if(m_projectionType == eGalerkin)
        {
            out << "\tProjection Type : Galerkin" <<endl;
        }
        else
        {
            out << "\tProjection Type : Discontinuous Galerkin" <<endl;
        }
    }

    /**
     * Write out a summary of the time parameters.
     * @param   out         Output stream to write to.
     */
    void ADRBase::TimeParamSummary(std::ostream &out)
    {
        out << "\tTime Step       : " << m_timestep << endl;
        out << "\tNo. of Steps    : " << m_steps << endl;
        out << "\tCheckpoints     : " << m_checksteps << " steps" << endl;
//        out << "\tInformation     : " << m_infosteps << " steps" << endl;
    }


    /**
     * Performs a case-insensitive string comparison (from web).
     * @param   s1          First string to compare.
     * @param   s2          Second string to compare.
     * @returns             0 if the strings match.
     */
    int ADRBase::NoCaseStringCompare(const string & s1, const string& s2)
    {
        //if (s1.size() < s2.size()) return -1;
        //if (s1.size() > s2.size()) return 1;

        string::const_iterator it1=s1.begin();
        string::const_iterator it2=s2.begin();

        //stop when either string's end has been reached
        while ( (it1!=s1.end()) && (it2!=s2.end()) )
        {
            if(::toupper(*it1) != ::toupper(*it2)) //letters differ?
            {
                // return -1 to indicate smaller than, 1 otherwise
                return (::toupper(*it1)  < ::toupper(*it2)) ? -1 : 1;
            }

            //proceed to the next character in each string
            ++it1;
            ++it2;
        }

        size_t size1=s1.size();
        size_t size2=s2.size();// cache lengths

        //return -1,0 or 1 according to strings' lengths
        if (size1==size2)
        {
            return 0;
        }

        return (size1 < size2) ? -1 : 1;
    }

    void ADRBase::LoadParameter(std::string name, int &var, int def)
    {
        if(m_boundaryConditions->CheckForParameter(name) == true)
        {
            var = m_boundaryConditions->GetParameter(name);
        }
        else
        {
            var  = def;
        }
    }

    void ADRBase::LoadParameter(std::string name, NekDouble &var, NekDouble def)
    {
        if(m_boundaryConditions->CheckForParameter(name) == true)
        {
            var = m_boundaryConditions->GetParameter(name);
        }
        else
        {
            var  = def;
        }
    }

    /**
     * Evaluates the exact solution provided in the session for a given
     * dependent variable.
     * @param   field               The index of the field to evaluate.
     * @param   outfield            Storage for exact solution.
     * @param   time                The time at which to evaluate the solution.
     */
    void ADRBase::v_EvaluateExactSolution(int field,
					  Array<OneD, NekDouble> &outfield,
					  const NekDouble time)
    {
        int nq = m_fields[field]->GetNpoints();
        bool Readit = true;

        Array<OneD,NekDouble> x0(nq);
        Array<OneD,NekDouble> x1(nq);
        Array<OneD,NekDouble> x2(nq);

        // get the coordinates of the quad points
        m_fields[field]->GetCoords(x0,x1,x2);

        SpatialDomains::ConstExactSolutionShPtr ifunc
                        = m_boundaryConditions->GetExactSolution(field);
        for(int j = 0; j < nq; j++)
        {
            outfield[j] = ifunc->Evaluate(x0[j],x1[j],x2[j],time);
        }
    }

} //end of namespace

/**
* $Log: ADRBase.cpp,v $
* Revision 1.27  2010/02/26 13:52:47  cantwell
* Tested and fixed where necessary Hex/Tet projection and differentiation in
*   StdRegions, and LocalRegions for regular and deformed (where applicable).
* Added SpatialData and SpatialParameters classes for managing spatiall-varying
*   data.
* Added TimingGeneralMatrixOp3D for timing operations on 3D geometries along
*   with some associated input meshes.
* Added 3D std and loc projection demos for tet and hex.
* Added 3D std and loc regression tests for tet and hex.
* Fixed bugs in regression tests in relation to reading OK files.
* Extended Elemental and Global optimisation parameters for 3D expansions.
* Added GNUPlot output format option.
* Updated ADR2DManifoldSolver to use spatially varying data.
* Added Barkley model to ADR2DManifoldSolver.
* Added 3D support to FldToVtk and XmlToVtk.
* Renamed History.{h,cpp} to HistoryPoints.{h,cpp}
*
* Revision 1.26  2010/02/02 13:53:26  cantwell
* Moved reading in of history data to separate SpatialDomains class.
* Updated AlievPanfilov demo to move history specification.
* Replaced FindNektar line and NEKTAR_BIN_DIR def in regressionTests
* CMakeLists.txt as this is required to locate the regression test execs.
*
* Revision 1.25  2010/01/27 15:55:57  cantwell
* Fixed incorrect ordering of history point data.
* Fixed parsing of session name when session filename contains multiple
*   full-stops.
* Removed extra empty composite entries in XML generation from Gmsh.
*
* Revision 1.24  2010/01/27 13:19:13  cantwell
* Added functions to write history/probe data during timestepping.
*
* Revision 1.23  2010/01/26 17:43:08  cantwell
* Updated CMakeLists.txt to build FitzHughNagumoSolver
* Added Aliev-Panfilov model to ADR2DManifoldSolver
*
* Revision 1.22  2009/12/14 17:59:33  cbiotto
* Adding writing tecplot file
*
* Revision 1.21  2009/12/09 12:37:12  cbiotto
* Update for regression test
*
* Revision 1.20  2009/11/02 19:15:43  cantwell
* Moved ContField1D to inherit from DisContField1D.
* Moved ContField3D to inherit from DisContField3D.
* Incorporated GenExpList1D functionality into ExpList1D.
* Tidied up and added documentation to various classes.
* Moved Namespace documentation and introductions to separate files along with
* doxygen configuration.
* Added option to use system ZLIB library instead of libboost_zlib on UNIX.
* Added extra search paths to FindMetis.cmake and FindNektar++.cmake.
* Updated Linux compiling instructions.
* Updated regDemo to use Helmholtz2D-g when built as debug.
*
* Revision 1.19  2009/10/07 16:52:59  cbiotto
* Updating write functions
*
* Revision 1.18  2009/09/07 11:21:58  sherwin
* Updates related to Navier-Stokes solver
*
* Revision 1.17  2009/08/14 09:30:01  cbiotto
* Add WriteVar function
*
* Revision 1.16  2009/07/29 09:19:42  sehunchun
* Generalization of WeakDGDiffusion
*
* Revision 1.15  2009/07/23 05:23:21  sehunchun
* WeakDiffusion operator is updated
*
* Revision 1.14  2009/07/11 23:39:23  sehunchun
* Move uncommon functions to each solvers
*
* Revision 1.13  2009/07/09 21:29:13  sehunchun
* Add SetUpSurfaceNormal function..
*
* Revision 1.12  2009/07/02 15:57:36  sehunchun
* "ReadBoundaryCondition" options with extenstion to 2D geometry imbedded in 3D
*
* Revision 1.11  2009/07/01 21:55:00  sehunchun
* Changes of WeakDiffusion according to updates
*
* Revision 1.10  2009/06/11 01:54:08  claes
* Added Inviscid Burger
*
* Revision 1.9  2009/04/29 20:45:55  sherwin
* Update for new definition of enum
*
* Revision 1.8  2009/04/27 21:37:14  sherwin
* Updated to dump .fld and .chk file in compressed coefficient format
*
* Revision 1.7  2009/03/10 23:37:14  claes
* Updated the ShallowWaterSolver to work with the general timestepping scheme
*
* Revision 1.6  2009/03/04 14:17:38  pvos
* Removed all methods that take and Expansion as argument
*
* Revision 1.5  2009/02/28 22:00:38  sehunchun
*  Explicit Diffusion solver is added
*
* Revision 1.4  2009/02/03 14:33:44  pvos
* Modifications for solvers with time-dependent dirichlet BC's
*
* Revision 1.3  2009/02/02 16:10:16  claes
* Update to make SWE, Euler and Boussinesq solvers up to date with the time integrator scheme. Linear and classical Boussinsq solver working
*
* Revision 1.2  2009/01/27 12:07:18  pvos
* Modifications to make cont. Galerkin Advection solver working
*
* Revision 1.1  2009/01/13 10:59:32  pvos
* added new solvers file
*
* Revision 1.8  2009/01/10 23:50:32  sherwin
* Update ContField1D/2D to use different variable in constructors
*
* Revision 1.7  2009/01/06 21:11:03  sherwin
* Updates for Virtual ExpList calls
*
* Revision 1.6  2008/11/17 08:10:07  claes
* Removed functions that were no longer used after the solver library was restructured
*
* Revision 1.5  2008/11/02 22:39:27  sherwin
* Updated naming convention
*
* Revision 1.4  2008/10/31 10:50:10  pvos
* Restructured directory and CMakeFiles
*
* Revision 1.3  2008/10/29 22:51:07  sherwin
* Updates for const correctness and ODEforcing
*
* Revision 1.2  2008/10/19 15:59:20  sherwin
* Added Summary method
*
* Revision 1.1  2008/10/16 15:25:45  sherwin
* Working verion of restructured AdvectionDiffusionReactionSolver
*
**/
