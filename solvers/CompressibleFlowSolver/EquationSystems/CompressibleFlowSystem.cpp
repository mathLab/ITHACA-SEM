///////////////////////////////////////////////////////////////////////////////
//
// File CompressibleFlowSystem.cpp
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
// Description: Compressible flow system base class with auxiliary functions
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystem.h>

using namespace std;

namespace Nektar
{
    CompressibleFlowSystem::CompressibleFlowSystem(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          AdvectionSystem(pSession, pGraph)
    {
    }

    /**
     * @brief Initialization object for CompressibleFlowSystem class.
     */
    void CompressibleFlowSystem::v_InitObject()
    {
        AdvectionSystem::v_InitObject();

        m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
                    m_session, m_spacedim);

        ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
                 "No UPWINDTYPE defined in session.");

        // Do not forwards transform initial condition
        m_homoInitialFwd = false;

        // Set up locations of velocity vector.
        m_vecLocs = Array<OneD, Array<OneD, NekDouble> >(1);
        m_vecLocs[0] = Array<OneD, NekDouble>(m_spacedim);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_vecLocs[0][i] = 1 + i;
        }

        // Loading parameters from session file
        InitialiseParameters();

        // Setting up advection and diffusion operators
        InitAdvection();

        // Create artificial diffusion
        if (m_shockCaptureType != "Off")
        {
            if (m_shockCaptureType == "Physical")
            {
                int nPts = m_fields[0]->GetTotPoints();
                m_muav = Array<OneD, NekDouble>(nPts, 0.0);

                int nTracePts = m_fields[0]->GetTrace()->GetTotPoints();
                m_muavTrace = Array<OneD, NekDouble> (nTracePts,0.0);
            }
            else
            {
                m_artificialDiffusion = GetArtificialDiffusionFactory()
                                        .CreateInstance(m_shockCaptureType,
                                                        m_session,
                                                        m_fields,
                                                        m_spacedim);
            }
        }

        // Forcing terms for the sponge region
        m_forcing = SolverUtils::Forcing::Load(m_session, shared_from_this(),
                                        m_fields, m_fields.size());

        // User-defined boundary conditions
        int cnt = 0;
        for (int n = 0; n < m_fields[0]->GetBndConditions().size(); ++n)
        {
            std::string type =
                m_fields[0]->GetBndConditions()[n]->GetUserDefined();

            if (m_fields[0]->GetBndConditions()[n]->GetBoundaryConditionType()
                == SpatialDomains::ePeriodic)
            {
                continue;
            }

            if (!type.empty())
            {
                m_bndConds.push_back(GetCFSBndCondFactory().CreateInstance(
                        type,
                        m_session,
                        m_fields,
                        m_traceNormals,
                        m_spacedim,
                        n,
                        cnt));
            }
            cnt += m_fields[0]->GetBndCondExpansions()[n]->GetExpSize();
        }

        if (m_explicitAdvection)
        {
            m_ode.DefineOdeRhs    (&CompressibleFlowSystem::DoOdeRhs, this);
            m_ode.DefineProjection(&CompressibleFlowSystem::DoOdeProjection, this);
        }
        else
        {
            m_ode.DefineOdeRhs    (&CompressibleFlowSystem::DoOdeRhs, this);
            m_ode.DefineProjection(&CompressibleFlowSystem::DoOdeProjection, this);
            // m_ode.DefineImplicitSolve    (&CompressibleFlowSystem::DoImplicitSolve, this);
            m_ode.DefineImplicitSolve    (&CompressibleFlowSystem::DoImplicitSolve_phy2coeff, this);

            //TODO: NekLinSysIterative as a member to avoid repeted initialization
            LibUtilities::CommSharedPtr v_Comm  = m_fields[0]->GetComm()->GetRowComm();
            m_linsol    = MemoryManager<NekLinSysIterative>::AllocateSharedPtr(m_session,v_Comm); 

            m_maxLinItePerNewton = m_linsol->GetMaxLinIte()*m_MaxNonlinIte+m_MaxNonlinIte;
            // m_linsol    =   NekLinSysIterative(m_session,v_Comm);
            // m_LinSysOprtors.DefineMatrixMultiply(&CompressibleFlowSystem::MatrixMultiply_MatrixFree_coeff, this);
            m_LinSysOprtors.DefineMatrixMultiply(&CompressibleFlowSystem::MatrixMultiply_MatrixFree_coeff_dualtimestep, this);
            m_LinSysOprtors.DefinePrecond(&CompressibleFlowSystem::preconditioner_BlkSOR_coeff, this);
            // m_LinSysOprtors.DefinePrecond(&CompressibleFlowSystem::preconditioner_BlkDiag, this);
            // m_LinSysOprtors.DefinePrecond(&CompressibleFlowSystem::preconditioner, this);
            m_linsol->setLinSysOperators(m_LinSysOprtors);

            if (boost::iequals(m_session->GetSolverInfo("PRECONDITIONER"),
                               "IncompleteLU"))
            {
                // m_LinSysOprtors.DefinePrecond(&CompressibleFlowSystem::preconditioner_BlkILU_coeff, this);
                m_PrecMatStorage    =   eSparse;

                ASSERTL0(false,"IncompleteLU preconditioner not finished yet");

                // DNekSmvBsrMat::SparseStorageSharedPtr sparseStorage =
                //             MemoryManager<DNekSmvBsrMat::StorageType>::
                //                     AllocateSharedPtr(
                //                         brows, bcols, block_size, bcoMat, matStorage );

                // // Create sparse matrix
                // m_smvbsrmatrix = MemoryManager<DNekSmvBsrMat>::
                //                         AllocateSharedPtr( sparseStorage );

                // matBytes = m_smvbsrmatrix->GetMemoryFootprint();

            }
            else
            {
                int nvariables  =   m_fields.num_elements();
                m_LinSysOprtors.DefinePrecond(&CompressibleFlowSystem::preconditioner_BlkSOR_coeff, this);
                m_PrecMatStorage    =   eDiagonal;
                m_session->LoadParameter("nPadding",     m_nPadding      ,    4);
                
                int ntmp=0;
                m_session->LoadParameter("PrecondMatDataSingle",                 ntmp      ,    1);
                m_flagPrecMatVarsSingle             = true;
                if(0==ntmp)
                {
                    m_flagPrecMatVarsSingle = false;
                }
                if(m_DebugNumJacBSOR)
                {
                    m_flagPrecMatVarsSingle = false;
                }

                m_session->LoadParameter("flagPrecondCacheOptmis",                 ntmp      ,    1);
                m_flagPrecondCacheOptmis             = true;
                if(0==ntmp)
                {
                    m_flagPrecondCacheOptmis = false;
                }

                // cout << " flagPrecondCacheOptmis= "<<m_flagPrecondCacheOptmis<<endl;
                
                if(m_flagPrecMatVarsSingle)
                {
                    m_PrecMatVarsSingle = Array<OneD, Array<OneD, SNekBlkMatSharedPtr> >(nvariables);
                    for(int i = 0; i < nvariables; i++)
                    {
                        m_PrecMatVarsSingle[i] =  Array<OneD, SNekBlkMatSharedPtr> (nvariables);
                    }
                    AllocatePrecondBlkDiag_coeff(m_PrecMatVarsSingle);

                    if(m_flagPrecondCacheOptmis)
                    {
                        int nelmts  = m_fields[0]->GetNumElmts();
                        int nelmtcoef;
                        Array<OneD, unsigned int > nelmtmatdim(nelmts);
                        for(int i = 0; i < nelmts; i++)
                        {
                            nelmtcoef   =   m_fields[0]->GetExp(i)->GetNcoeffs();
                            nelmtmatdim[i]  =   nelmtcoef*nvariables;
                        }
                        AllocateNekBlkMatDig(m_PrecMatSingle,nelmtmatdim,nelmtmatdim);
                    }
                }
                else
                {
                    m_PrecMatVars = Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >(nvariables);
                    for(int i = 0; i < nvariables; i++)
                    {
                        m_PrecMatVars[i] =  Array<OneD, DNekBlkMatSharedPtr> (nvariables);
                    }
                    AllocatePrecondBlkDiag_coeff(m_PrecMatVars);
                    if(m_flagPrecondCacheOptmis)
                    {
                        int nelmts  = m_fields[0]->GetNumElmts();
                        int nelmtcoef;
                        Array<OneD, unsigned int > nelmtmatdim(nelmts);
                        for(int i = 0; i < nelmts; i++)
                        {
                            nelmtcoef   =   m_fields[0]->GetExp(i)->GetNcoeffs();
                            nelmtmatdim[i]  =   nelmtcoef*nvariables;
                        }
                        AllocateNekBlkMatDig(m_PrecMat,nelmtmatdim,nelmtmatdim);
                    }
                }
            }

            int nvariables  =   m_fields.num_elements();
            Array<OneD, Array<OneD, Array<OneD, int > > >   map;
            bool flag;
            const MultiRegions::LocTraceToTraceMapSharedPtr locTraceToTraceMap = m_fields[0]->GetlocTraceToTraceMap();
            m_fields[0]->CalcuTracephysToLeftRightExpphysMap(flag,map);
            locTraceToTraceMap->SetTracephysToLeftRightExpphysMap(map);
            locTraceToTraceMap->SetflagTracephysToLeftRightExpphysMap(flag);

            locTraceToTraceMap->CalcuLocTracephysToTraceIDMap(m_fields[0]->GetTrace(),m_spacedim);
            for(int i=1;i<nvariables;i++)
            {
                m_fields[i]->GetlocTraceToTraceMap()->SetTracephysToLeftRightExpphysMap(map);
                m_fields1[i]->GetlocTraceToTraceMap()->SetflagTracephysToLeftRightExpphysMap(flag);
                m_fields[i]->GetlocTraceToTraceMap()->SetLocTracephysToTraceIDMap(
                    locTraceToTraceMap->GetLocTracephysToTraceIDMap()    );
            }

#else

            ASSERTL0(false, "Implicit CFS not set up.");
        }

        SetBoundaryConditionsBwdWeight();

        string advName;
        m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");
        m_session->MatchSolverInfo(
            "useUnifiedWeakIntegration", "True",
            m_useUnifiedWeakIntegration, false);
        if(m_useUnifiedWeakIntegration)
        {
            if(advName=="WeakDG" && m_shockCaptureType!="Smooth"&&(eNotHomogeneous == m_HomogeneousType))
            {

            }
            else
            {
                m_useUnifiedWeakIntegration=false;
                if(m_session->DefinesCmdLineArgument("verbose"))
                {
                    WARNINGL0(false, "useUnifiedWeakIntegration not coded for these parameters of Advection");
                }
            }
        }
    }

    /**
     * @brief Destructor for CompressibleFlowSystem class.
     */
    CompressibleFlowSystem::~CompressibleFlowSystem()
    {

    }

    /**
     * @brief Load CFS parameters from the session file
     */
    void CompressibleFlowSystem::InitialiseParameters()
    {
        // Get gamma parameter from session file.
        m_session->LoadParameter("Gamma", m_gamma, 1.4);

        // Shock capture
        m_session->LoadSolverInfo("ShockCaptureType",
                                  m_shockCaptureType, "Off");

        // Load parameters for exponential filtering
        m_session->MatchSolverInfo("ExponentialFiltering","True",
                                   m_useFiltering, false);
        if (m_useFiltering)
        {
            m_session->LoadParameter ("FilterAlpha", m_filterAlpha, 36);
            m_session->LoadParameter ("FilterExponent", m_filterExponent, 16);
            m_session->LoadParameter ("FilterCutoff", m_filterCutoff, 0);
        }

        // Load CFL for local time-stepping (for steady state)
        m_session->MatchSolverInfo("LocalTimeStep","True",
                                   m_useLocalTimeStep, false);
        if (m_useLocalTimeStep)
        {
            ASSERTL0(m_cflSafetyFactor != 0,
                    "Local time stepping requires CFL parameter.");
        }

        m_session->LoadParameter ("JFEps", m_JFEps, 5.0E-8);

        int ntmp;
        m_session->LoadParameter("DEBUG_ADVECTION_JAC_MAT",     ntmp      ,    1);
        m_DEBUG_ADVECTION_JAC_MAT             = true;
        if(0==ntmp)
        {
            m_DEBUG_ADVECTION_JAC_MAT = false;
        }
        m_session->LoadParameter("DEBUG_VISCOUS_JAC_MAT",                 ntmp      ,    1);
        m_DEBUG_VISCOUS_JAC_MAT             = true;
        if(0==ntmp)
        {
            m_DEBUG_VISCOUS_JAC_MAT = false;
        }
        m_session->LoadParameter("DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT",     ntmp      ,    0);
        m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT = false;
        if(1==ntmp)
        {
            m_DEBUG_VISCOUS_TRACE_DERIV_JAC_MAT = true;
        }
    }

    /**
     * @brief Create advection and diffusion objects for CFS
     */
    void CompressibleFlowSystem::InitAdvection()
    {
        // Check if projection type is correct
        ASSERTL0(m_projectionType == MultiRegions::eDiscontinuous,
                "Unsupported projection type.");

        string advName, riemName;
        m_session->LoadSolverInfo("AdvectionType", advName, "WeakDG");

        m_advObject = SolverUtils::GetAdvectionFactory()
                                    .CreateInstance(advName, advName);

        if (m_specHP_dealiasing)
        {
            m_advObject->SetFluxVector(&CompressibleFlowSystem::
                                       GetFluxVectorDeAlias, this);
        }
        else
        {
            m_advObject->SetFluxVector  (&CompressibleFlowSystem::
                                          GetFluxVector, this);
        }

        // Setting up Riemann solver for advection operator
        m_session->LoadSolverInfo("UpwindType", riemName, "Average");

        SolverUtils::RiemannSolverSharedPtr riemannSolver;
        riemannSolver = SolverUtils::GetRiemannSolverFactory()
                                    .CreateInstance(riemName, m_session);

        // Setting up parameters for advection operator Riemann solver
        riemannSolver->SetParam (
            "gamma",   &CompressibleFlowSystem::GetGamma,   this);
        riemannSolver->SetAuxVec(
            "vecLocs", &CompressibleFlowSystem::GetVecLocs, this);
        riemannSolver->SetVector(
            "N",       &CompressibleFlowSystem::GetNormals, this);

        // Concluding initialisation of advection / diffusion operators
        m_advObject->SetRiemannSolver   (riemannSolver);
        m_advObject->InitObject         (m_session, m_fields);
    }

    /**
     * @brief Compute the right-hand side.
     */
    void CompressibleFlowSystem::DoOdeRhs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int nvariables = inarray.size();
        int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        m_BndEvaluateTime   = time;

        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble> > Fwd    (nvariables);
        Array<OneD, Array<OneD, NekDouble> > Bwd    (nvariables);

        if (m_HomogeneousType == eHomogeneous1D)
        {
            Fwd = NullNekDoubleArrayofArray;
            Bwd = NullNekDoubleArrayofArray;
        }
        else
        {
            for (int i = 0; i < nvariables; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }

        //Only test solver use the reduced coddes
        if(m_useUnifiedWeakIntegration)
        {
            int nDim        = m_spacedim;
            int nCoeffs     = GetNcoeffs();
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>> VolumeFlux1(nvariables);
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>> VolumeFlux2(nDim);
            Array<OneD, Array<OneD, NekDouble>> TraceFlux1(nvariables);
            Array<OneD, Array<OneD, NekDouble>> TraceFlux2(nvariables);

            for (int i = 0; i < nvariables; ++i)
            {
                VolumeFlux1[i] = Array<OneD, Array<OneD, NekDouble>>(nDim);
                for (int j= 0; j < nDim; ++j)
                {
                    VolumeFlux1[i][j] = Array<OneD, NekDouble>(npoints,0.0);
                }
            }
            for (int j = 0; j < nDim; ++j)
            {
                VolumeFlux2[j] = Array<OneD, Array<OneD, NekDouble>>(nvariables);
                for (int i = 0; i < nvariables; ++i)
                {
                    VolumeFlux2[j][i] = Array<OneD, NekDouble>(npoints,0.0);
                }
            }
            for (int i = 0; i < nvariables; ++i)
            {
                TraceFlux1[i]  = Array<OneD, NekDouble>(nTracePts, 0.0);
                TraceFlux2[i]  = Array<OneD, NekDouble>(nTracePts, 0.0);
            }

            Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);
            m_advObject->AdvectVolumeFlux(nvariables, m_fields, advVel, inarray, VolumeFlux1, time);
            m_advObject->AdvectTraceFlux(nvariables, m_fields, advVel, inarray,TraceFlux1, time, Fwd, Bwd);
            v_DoDiffusionFlux(inarray, VolumeFlux2, TraceFlux2, Fwd, Bwd);

            // Add Trace and Volume Integral together
            Array<OneD, NekDouble> tmp(nCoeffs, 0.0);
            for (int i = 0; i < nvariables; ++i)
            {
                // Add Volume integral
                for (int j = 0; j < nDim; ++j)
                {
                    // Advection term needs to be negative
                    // Add Advection and Diffusion part
                    Vmath::Vsub(npoints, &VolumeFlux2[j][i][0], 1,
                                &VolumeFlux1[i][j][0], 1, &VolumeFlux1[i][j][0], 1);
                }
                m_fields[i]->IProductWRTDerivBase(VolumeFlux1[i], tmp);

                Vmath::Neg(nCoeffs, &tmp[0], 1);
                // Add Trace integral
                // Advection term needs to be negative
                // Add Advection and Diffusion part
                Vmath::Vsub(nTracePts, &TraceFlux2[i][0], 1, &TraceFlux1[i][0], 1,
                            &TraceFlux1[i][0], 1);
                m_fields[i]->AddTraceIntegral(TraceFlux1[i], tmp);
                m_fields[i]->MultiplyByElmtInvMass(tmp, tmp);
                m_fields[i]->BwdTrans(tmp, outarray[i]);
            }

            // AddDiffusionSymmFluxToPhys(inarray, VolumeFlux2, outarray, Fwd, Bwd);
        }
        else
        {
            //Oringinal CompressibleFlowSolver
            // Calculate advection
            DoAdvection(inarray, outarray, time, Fwd, Bwd);

            // Negate results
            for (int i = 0; i < nvariables; ++i)
            {
                Vmath::Neg(npoints, outarray[i], 1);
            }

            // Add diffusion terms
            DoDiffusion(inarray, outarray, Fwd, Bwd);

        }

        // Add forcing terms
        for (auto &x : m_forcing)
        {
            x->Apply(m_fields, inarray, outarray, time);
        }

        if (m_useLocalTimeStep)
        {
            int nElements = m_fields[0]->GetExpSize();
            int nq, offset;
            NekDouble fac;
            Array<OneD, NekDouble> tmp;

            Array<OneD, NekDouble> tstep (nElements, 0.0);
            GetElmtTimeStep(inarray, tstep);

            // Loop over elements
            for(int n = 0; n < nElements; ++n)
            {
                nq     = m_fields[0]->GetExp(n)->GetTotPoints();
                offset = m_fields[0]->GetPhys_Offset(n);
                fac    = tstep[n] / m_timestep;
                for(int i = 0; i < nvariables; ++i)
                {
                    Vmath::Smul(nq, fac, outarray[i] + offset, 1,
                                         tmp = outarray[i] + offset, 1);
                }
            }
        }
    }

    /**
     * @brief Compute the projection and call the method for imposing the
     * boundary conditions in case of discontinuous projection.
     */
    void CompressibleFlowSystem::DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();

        switch(m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                int npoints = GetNpoints();

                for(i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                    if(m_useFiltering)
                    {
                        m_fields[i]->ExponentialFilter(outarray[i],
                            m_filterAlpha, m_filterExponent, m_filterCutoff);
                    }
                }
                SetBoundaryConditions(outarray, time);
                break;
            }
            case MultiRegions::eGalerkin:
            case MultiRegions::eMixed_CG_Discontinuous:
            {
                ASSERTL0(false, "No Continuous Galerkin for full compressible "
                                "Navier-Stokes equations");
                break;
            }
            default:
                ASSERTL0(false, "Unknown projection scheme");
                break;
        }
    }


    void CompressibleFlowSystem::preconditioner(
                                                 const Array<OneD, NekDouble> &inarray,
                                                 Array<OneD, NekDouble >&out)
    {
        int ntotal     = inarray.num_elements();
        Vmath::Vcopy(ntotal,inarray,1,out,1);
        return;
    }

    
    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::preconditioner_BlkDiag(
        const Array<OneD, NekDouble>                                &inarray,
        Array<OneD, NekDouble >                                     &outarray,
        const Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> >     &PrecMatVars,
        const DataType                                              &tmpDataType)
    {
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        Array<OneD, Array<OneD, DataType > >Sinarray(nvariables);
        Array<OneD, NekVector<DataType> >tmpVect(nvariables);
        Array<OneD, DataType > Soutarray(npoints);
        NekVector<DataType> outVect(npoints,Soutarray,eWrapper);

        for(int m = 0; m < nvariables; m++)
        {
            int moffset = m*npoints;
            Sinarray[m] = Array<OneD, DataType > (npoints);
            for(int i=0;i<npoints;i++)
            {
                Sinarray[m][i]  =  DataType(inarray[moffset+i]);
            }
            tmpVect[m] =  NekVector<DataType> (npoints,Sinarray[m],eWrapper);
        }

        Vmath::Fill(outarray.num_elements(),0.0,outarray,1);

        for(int m = 0; m < nvariables; m++)
        {
            Vmath::Zero(npoints,Soutarray,1);
            int moffset = m*npoints;
            for(int n = 0; n < nvariables; n++)
            {
                outVect += (*PrecMatVars[m][n])*tmpVect[n];
            }

            for(int i=0;i<npoints;i++)
            {
                outarray[moffset+i]  =  NekDouble(Soutarray[i]);
            }
        }
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::preconditioner_BlkDiag(
        const Array<OneD, NekDouble>                        &inarray,
        Array<OneD, NekDouble >                             &outarray,
        const TypeNekBlkMatSharedPtr                        &PrecMatVars,
        const DataType                                      &tmpDataType)
    {
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        unsigned int npointsVar = nvariables*npoints;
        Array<OneD, DataType >Sinarray(npointsVar);
        Array<OneD, DataType > Soutarray(npointsVar);
        NekVector<DataType> tmpVect(npointsVar,Sinarray,eWrapper);
        NekVector<DataType> outVect(npointsVar,Soutarray,eWrapper);

        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();

        for(int m = 0; m < nvariables; m++)
        {
            int nVarOffset = m*npoints;
            for(int ne=0;ne<ntotElmt;ne++)
            {
                int nCoefOffset = GetCoeff_Offset(ne);
                int nElmtCoef   = GetNcoeffs(ne);
                int inOffset    = nVarOffset+nCoefOffset;
                int outOffset   = nCoefOffset*nvariables+m*nElmtCoef;
                for(int i=0;i<nElmtCoef;i++)
                {
                    Sinarray[outOffset+i]  =  DataType(inarray[inOffset+i]);
                }
            }
        }

        outVect = (*PrecMatVars)*tmpVect;

        for(int m = 0; m < nvariables; m++)
        {
            int nVarOffset = m*npoints;
            for(int ne=0;ne<ntotElmt;ne++)
            {
                int nCoefOffset = GetCoeff_Offset(ne);
                int nElmtCoef   = GetNcoeffs(ne);
                int inOffset    = nVarOffset+nCoefOffset;
                int outOffset   = nCoefOffset*nvariables+m*nElmtCoef;
                for(int i=0;i<nElmtCoef;i++)
                {
                    outarray[inOffset+i]  =  NekDouble(Soutarray[outOffset+i]);
                }
            }
        }
    }

    void CompressibleFlowSystem::preconditioner_NumJac(
        const Array<OneD, NekDouble>                                                &inarray,
        Array<OneD, NekDouble >                                                     &outarray,
        const Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >                        &PrecMatVars,
        const Array<OneD, Array<OneD, NekDouble > >                                 &PrecMatVarsOffDiag)
    {
        const NekDouble SORParam        =   m_SORRelaxParam;
        const NekDouble OmSORParam      =   1.0-SORParam;

        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        unsigned int ntotpnt    = inarray.num_elements();
        
        ASSERTL0(nvariables*npoints==ntotpnt,"nvariables*npoints==ntotpnt not satisfied in preconditioner_BlkSOR");

        Array<OneD, NekDouble> rhs(ntotpnt);

        Array<OneD, NekDouble>  outN(ntotpnt);
        Array<OneD, NekDouble>  outTmp(ntotpnt);
        Vmath::Vcopy(ntotpnt,&inarray[0],1,&rhs[0],1);

        NekDouble tmpDouble = 0.0;
        preconditioner_BlkDiag(rhs,outarray,PrecMatVars,tmpDouble);

        int nSORTot   =   m_JFNKPrecondStep;
        for(int nsor = 0; nsor < nSORTot-1; nsor++)
        {
            Vmath::Smul(ntotpnt,OmSORParam,outarray,1,outN,1);
            
            MinusOffDiag2RhsNumJac(nvariables,npoints,rhs,outarray,PrecMatVarsOffDiag);

            preconditioner_BlkDiag(outarray,outTmp,PrecMatVars,tmpDouble);
            Vmath::Svtvp(ntotpnt,SORParam,outTmp,1,outN,1,outarray,1);
        }
    }

    void CompressibleFlowSystem::MinusOffDiag2RhsNumJac(
        const int                                                                   nvariables,
        const int                                                                   nCoeffs,
        const Array<OneD, NekDouble>                                                &rhs,
        Array<OneD, NekDouble>                                                      &outarray,
        const Array<OneD, Array<OneD, NekDouble > >                                 &PrecMatVarsOffDiag)
    {
        int nTotCoef = nvariables*nCoeffs;
        Array<OneD, NekDouble> tmp (nTotCoef,0.0);
        for(int i=0;i<nTotCoef;i++)
        {
            Vmath::Svtvp(nTotCoef,outarray[i],&PrecMatVarsOffDiag[i][0],1,&tmp[0],1,&tmp[0],1);
        }

        // for(int i=0;i<nTotCoef;i++)
        // {
        //     for(int j=0;j<nTotCoef;j++)
        //     {
        //         tmp[j] += outarray[i]*PrecMatVarsOffDiag[j][i];
        //     }
        // }

        Vmath::Vsub(nTotCoef,&rhs[0],1,&tmp[0],1,&outarray[0],1);
    }

    void CompressibleFlowSystem::preconditioner_BlkSOR_coeff(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&outarray,
            const bool                   &flag)
    {
        int nSORTot   =   m_JFNKPrecondStep;
    /**
     * @brief Compute the advection terms for the right-hand side
     */
    void CompressibleFlowSystem::DoAdvection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time,
        const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >       &pBwd)
    {
        int nvariables = inarray.size();
        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);

        m_advObject->Advect(nvariables, m_fields, advVel, inarray,
                            outarray, time, pFwd, pBwd);
    }

    /**
     * @brief Add the diffusions terms to the right-hand side
     */
    void CompressibleFlowSystem::DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd)
    {
        v_DoDiffusion(inarray, outarray, pFwd, pBwd);

        if (m_shockCaptureType != "Off" && m_shockCaptureType != "Physical")
        {
            m_artificialDiffusion->DoArtificialDiffusion(inarray, outarray);
        }
    }

    /**
     * @brief Add the diffusions terms to the right-hand side
     * Similar to DoDiffusion() but with outarray in coefficient space
     */
    void CompressibleFlowSystem::DoDiffusion_coeff(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd)
    {
        v_DoDiffusion_coeff(inarray, outarray, pFwd, pBwd);
    }

    void CompressibleFlowSystem::SetBoundaryConditions(
            Array<OneD, Array<OneD, NekDouble> >             &physarray,
            NekDouble                                         time)
    {
        int nTracePts  = GetTraceTotPoints();
        int nvariables = physarray.size();

        Array<OneD, Array<OneD, NekDouble> > Fwd(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            Fwd[i] = Array<OneD, NekDouble>(nTracePts);
            m_fields[i]->ExtractTracePhys(physarray[i], Fwd[i]);
        }

        if (m_bndConds.size())
        {
            // Loop over user-defined boundary conditions
            for (auto &x : m_bndConds)
            {
                x->Apply(Fwd, physarray, time);
            }
        }
    }
    /**
     * @brief Set up a weight on physical boundaries for boundary condition 
     * applications
     */
    void CompressibleFlowSystem::SetBoundaryConditionsBwdWeight()
    {
        if (m_bndConds.size())
        {
            // Loop over user-defined boundary conditions
            for (auto &x : m_bndConds)
            {
                x->ApplyBwdWeight();
            }
        }
    }

    /**
     * @brief Return the flux vector for the compressible Euler equations.
     *
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >  &physfield,
        TensorOfArray3D<NekDouble>                  &flux)
    {
        int i, j;
        int nq = physfield[0].size();
        int nVariables = m_fields.size();

        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);

        // Flux vector for the rho equation
        for (i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nq);
            Vmath::Vcopy(nq, physfield[i+1], 1, flux[0][i], 1);
        }

        m_varConv->GetVelocityVector(physfield, velocity);
        m_varConv->GetPressure(physfield, pressure);

        // Flux vector for the velocity fields
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vmul(nq, velocity[j], 1, physfield[i+1], 1,
                            flux[i+1][j], 1);
            }

            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux[i+1][i], 1, pressure, 1, flux[i+1][i], 1);
        }

        // Flux vector for energy.
        Vmath::Vadd(nq, physfield[m_spacedim+1], 1, pressure, 1,
                    pressure, 1);

        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, pressure, 1,
                        flux[m_spacedim+1][j], 1);
        }

        // For the smooth viscosity model
        if (nVariables == m_spacedim+3)
        {
            // Add a zero row for the advective fluxes
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Zero(nq, flux[m_spacedim+2][j], 1);
            }
        }
    }

    /**
     * @brief Return the flux vector for the compressible Euler equations
     * by using the de-aliasing technique.
     *
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >      &physfield,
        TensorOfArray3D<NekDouble>                      &flux)
    {
        int i, j;
        int nq = physfield[0].size();
        int nVariables = m_fields.size();

        // Factor to rescale 1d points in dealiasing
        NekDouble OneDptscale = 2;
        nq = m_fields[0]->Get1DScaledTotPoints(OneDptscale);

        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);

        Array<OneD, Array<OneD, NekDouble> > physfield_interp(nVariables);
        TensorOfArray3D<NekDouble> flux_interp(nVariables);

        for (i = 0; i < nVariables; ++ i)
        {
            physfield_interp[i] = Array<OneD, NekDouble>(nq);
            flux_interp[i] = Array<OneD, Array<OneD, NekDouble> >(m_spacedim);
            m_fields[0]->PhysInterp1DScaled(
                OneDptscale, physfield[i], physfield_interp[i]);

            for (j = 0; j < m_spacedim; ++j)
            {
                flux_interp[i][j] = Array<OneD, NekDouble>(nq);
            }
        }

        // Flux vector for the rho equation
        for (i = 0; i < m_spacedim; ++i)
        {
            velocity[i] = Array<OneD, NekDouble>(nq);

            // Galerkin project solution back to original space
            m_fields[0]->PhysGalerkinProjection1DScaled(
                OneDptscale, physfield_interp[i+1], flux[0][i]);
        }

        m_varConv->GetVelocityVector(physfield_interp, velocity);
        m_varConv->GetPressure      (physfield_interp, pressure);

        // Evaluation of flux vector for the velocity fields
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                Vmath::Vmul(nq, velocity[j], 1, physfield_interp[i+1], 1,
                            flux_interp[i+1][j], 1);
            }

            // Add pressure to appropriate field
            Vmath::Vadd(nq, flux_interp[i+1][i], 1, pressure,1,
                        flux_interp[i+1][i], 1);
        }

        // Galerkin project solution back to original space
        for (i = 0; i < m_spacedim; ++i)
        {
            for (j = 0; j < m_spacedim; ++j)
            {
                m_fields[0]->PhysGalerkinProjection1DScaled(
                    OneDptscale, flux_interp[i+1][j], flux[i+1][j]);
            }
        }

        // Evaluation of flux vector for energy
        Vmath::Vadd(nq, physfield_interp[m_spacedim+1], 1, pressure, 1,
                    pressure, 1);

        for (j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nq, velocity[j], 1, pressure, 1,
                        flux_interp[m_spacedim+1][j], 1);

            // Galerkin project solution back to original space
            m_fields[0]->PhysGalerkinProjection1DScaled(
                OneDptscale,
                flux_interp[m_spacedim+1][j],
                flux[m_spacedim+1][j]);
        }
    }

    /**
     * @brief Calculate the maximum timestep on each element
     *        subject to CFL restrictions.
     */
    void CompressibleFlowSystem::GetElmtTimeStep(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD, NekDouble> &tstep)
    {
        boost::ignore_unused(inarray);

        int nElements = m_fields[0]->GetExpSize();

        // Change value of m_timestep (in case it is set to zero)
        NekDouble tmp = m_timestep;
        m_timestep    = 1.0;

        Array<OneD, NekDouble> cfl(nElements);
        cfl = GetElmtCFLVals();

        // Factors to compute the time-step limit
        NekDouble alpha     = MaxTimeStepEstimator();

        // Loop over elements to compute the time-step limit for each element
        for (int n = 0; n < nElements; ++n)
        {
            tstep[n] = m_cflSafetyFactor * alpha / cfl[n];
        }

        // Restore value of m_timestep
        m_timestep = tmp;
    }

    /**
     * @brief Calculate the maximum timestep subject to CFL restrictions.
     */
    NekDouble CompressibleFlowSystem::v_GetTimeStep(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray)
    {
        int nElements = m_fields[0]->GetExpSize();
        Array<OneD, NekDouble> tstep (nElements, 0.0);

        GetElmtTimeStep(inarray, tstep);

        // Get the minimum time-step limit and return the time-step
        NekDouble TimeStep = Vmath::Vmin(nElements, tstep, 1);
        m_comm->AllReduce(TimeStep, LibUtilities::ReduceMin);

        NekDouble tmp = m_timestep;
        m_timestep    = TimeStep;

        Array<OneD, NekDouble> cflNonAcoustic(nElements,0.0);
        cflNonAcoustic = GetElmtCFLVals(false);

        // Get the minimum time-step limit and return the time-step
        NekDouble MaxcflNonAcoustic = Vmath::Vmax(nElements, cflNonAcoustic, 1);
        m_comm->AllReduce(MaxcflNonAcoustic, LibUtilities::ReduceMax);

        m_cflNonAcoustic = MaxcflNonAcoustic;
        m_timestep = tmp;

        return TimeStep;
    }

    /**
     * @brief Set up logic for residual calculation.
     */
    void CompressibleFlowSystem::v_SetInitialConditions(
        NekDouble initialtime,
        bool      dumpInitialConditions,
        const int domain)
    {
        boost::ignore_unused(domain);

        EquationSystem::v_SetInitialConditions(initialtime, false);

        // insert white noise in initial condition
        NekDouble Noise;
        int phystot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> noise(phystot);

        m_session->LoadParameter("Noise", Noise,0.0);
        int m_nConvectiveFields =  m_fields.size();

        if (Noise > 0.0)
        {
            int seed = - m_comm->GetRank()*m_nConvectiveFields;
            for (int i = 0; i < m_nConvectiveFields; i++)
            {
                Vmath::FillWhiteNoise(phystot, Noise, noise, 1,
                                      seed);
                --seed;
                Vmath::Vadd(phystot, m_fields[i]->GetPhys(), 1,
                            noise, 1, m_fields[i]->UpdatePhys(), 1);
                m_fields[i]->FwdTrans_IterPerExp(m_fields[i]->GetPhys(),
                                                 m_fields[i]->UpdateCoeffs());
            }
        }

        if (dumpInitialConditions && m_checksteps)
        {
            Checkpoint_Output(m_nchk);
            m_nchk++;
        }
    }

    /**
     * @brief Compute the advection velocity in the standard space
     * for each element of the expansion.
     */
    Array<OneD, NekDouble> CompressibleFlowSystem::v_GetMaxStdVelocity(const NekDouble SpeedSoundFactor)
    {
        int nTotQuadPoints = GetTotPoints();
        int n_element      = m_fields[0]->GetExpSize();
        int expdim         = m_fields[0]->GetGraph()->GetMeshDimension();
        int nfields        = m_fields.size();
        int offset;
        Array<OneD, NekDouble> tmp;

        Array<OneD, Array<OneD, NekDouble> > physfields(nfields);
        for (int i = 0; i < nfields; ++i)
        {
            physfields[i] = m_fields[i]->GetPhys();
        }

        Array<OneD, NekDouble> stdV(n_element, 0.0);

        // Getting the velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > velocity   (m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(m_spacedim);
        Array<OneD, Array<OneD, NekDouble> > stdSoundSpeed(m_spacedim);
        Array<OneD, NekDouble>               soundspeed (nTotQuadPoints);
        LibUtilities::PointsKeyVector        ptsKeys;

        for (int i = 0; i < m_spacedim; ++i)
        {
            velocity   [i]   = Array<OneD, NekDouble>(nTotQuadPoints);
            stdVelocity[i]   = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
            stdSoundSpeed[i] = Array<OneD, NekDouble>(nTotQuadPoints, 0.0);
        }

        m_varConv->GetVelocityVector(physfields, velocity);
        m_varConv->GetSoundSpeed    (physfields, soundspeed);

        for (int el = 0; el < n_element; ++el)
        {
            ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();
            offset  = m_fields[0]->GetPhys_Offset(el);
            int nq = m_fields[0]->GetExp(el)->GetTotPoints();

            const SpatialDomains::GeomFactorsSharedPtr metricInfo =
                m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo();
            const Array<TwoD, const NekDouble> &gmat =
                m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()
                                                  ->GetDerivFactors(ptsKeys);

            // Convert to standard element
            //    consider soundspeed in all directions
            //    (this might overestimate the cfl)
            if (metricInfo->GetGtype() == SpatialDomains::eDeformed)
            {
                // d xi/ dx = gmat = 1/J * d x/d xi
                for (int i = 0; i < expdim; ++i)
                {
                    Vmath::Vmul(nq, gmat[i], 1,
                                    velocity[0] + offset, 1,
                                    tmp = stdVelocity[i] + offset, 1);
                    Vmath::Vmul(nq, gmat[i], 1,
                                    soundspeed + offset, 1,
                                    tmp = stdSoundSpeed[i] + offset, 1);
                    for (int j = 1; j < expdim; ++j)
                    {
                        Vmath::Vvtvp(nq, gmat[expdim*j+i], 1,
                                         velocity[j] + offset, 1,
                                         stdVelocity[i] + offset, 1,
                                         tmp = stdVelocity[i] + offset, 1);
                        Vmath::Vvtvp(nq, gmat[expdim*j+i], 1,
                                         soundspeed + offset, 1,
                                         stdSoundSpeed[i] + offset, 1,
                                         tmp = stdSoundSpeed[i] + offset, 1);
                    }
                }
            }
            else
            {
                for (int i = 0; i < expdim; ++i)
                {
                    Vmath::Smul(nq, gmat[i][0],
                                    velocity[0] + offset, 1,
                                    tmp = stdVelocity[i] + offset, 1);
                    Vmath::Smul(nq, gmat[i][0],
                                    soundspeed + offset, 1,
                                    tmp = stdSoundSpeed[i] + offset, 1);
                    for (int j = 1; j < expdim; ++j)
                    {
                        Vmath::Svtvp(nq, gmat[expdim*j+i][0],
                                         velocity[j] + offset, 1,
                                         stdVelocity[i] + offset, 1,
                                         tmp = stdVelocity[i] + offset, 1);
                        Vmath::Svtvp(nq, gmat[expdim*j+i][0],
                                         soundspeed + offset, 1,
                                         stdSoundSpeed[i] + offset, 1,
                                         tmp = stdSoundSpeed[i] + offset, 1);
                    }
                }
            }

            NekDouble vel;
            for (int i = 0; i < nq; ++i)
            {
                NekDouble pntVelocity = 0.0;
                for (int j = 0; j < expdim; ++j)
                {
                    // Add sound speed
                    vel = std::abs(stdVelocity[j][offset + i]) +
                          SpeedSoundFactor * std::abs(stdSoundSpeed[j][offset + i]);
                    pntVelocity += vel * vel;
                }
                pntVelocity = sqrt(pntVelocity);
                if (pntVelocity > stdV[el])
                {
                    stdV[el] = pntVelocity;
                }
            }
        }

        return stdV;
    }

    /**
     * @brief Set the denominator to compute the time step when a cfl
     * control is employed. This function is no longer used but is still
     * here for being utilised in the future.
     *
     * @param n   Order of expansion element by element.
     */
    NekDouble CompressibleFlowSystem::GetStabilityLimit(int n)
    {
        ASSERTL0(n <= 20, "Illegal modes dimension for CFL calculation "
                          "(P has to be less then 20)");

        NekDouble CFLDG[21] = {  2.0000,   6.0000,  11.8424,  19.1569,
                                27.8419,  37.8247,  49.0518,  61.4815,
                                75.0797,  89.8181, 105.6700, 122.6200,
                               140.6400, 159.7300, 179.8500, 201.0100,
                               223.1800, 246.3600, 270.5300, 295.6900,
                               321.8300}; //CFLDG 1D [0-20]
        NekDouble CFL = 0.0;

        if (m_projectionType == MultiRegions::eDiscontinuous)
        {
            CFL = CFLDG[n];
        }
        else
        {
            ASSERTL0(false, "Continuous Galerkin stability coefficients "
                            "not introduced yet.");
        }

        return CFL;
    }

    /**
     * @brief Compute the vector of denominators to compute the time step
     * when a cfl control is employed. This function is no longer used but
     * is still here for being utilised in the future.
     *
     * @param ExpOrder   Order of expansion element by element.
     */
    Array<OneD, NekDouble> CompressibleFlowSystem::GetStabilityLimitVector(
        const Array<OneD,int> &ExpOrder)
    {
        int i;
        Array<OneD,NekDouble> returnval(m_fields[0]->GetExpSize(), 0.0);
        for (i =0; i<m_fields[0]->GetExpSize(); i++)
        {
            returnval[i] = GetStabilityLimit(ExpOrder[i]);
        }
        return returnval;
    }

    void CompressibleFlowSystem::v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble> > &fieldcoeffs,
        std::vector<std::string>             &variables)
    {
        bool extraFields;
        m_session->MatchSolverInfo("OutputExtraFields","True",
                                   extraFields, true);
        if (extraFields)
        {
            const int nPhys   = m_fields[0]->GetNpoints();
            const int nCoeffs = m_fields[0]->GetNcoeffs();
            Array<OneD, Array<OneD, NekDouble> > tmp(m_fields.size());

            for (int i = 0; i < m_fields.size(); ++i)
            {
                tmp[i] = m_fields[i]->GetPhys();
            }

            Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);
            Array<OneD, Array<OneD, NekDouble> > velFwd  (m_spacedim);
            for (int i = 0; i < m_spacedim; ++i)
            {
                velocity[i] = Array<OneD, NekDouble> (nPhys);
                velFwd[i]   = Array<OneD, NekDouble> (nCoeffs);
            }

            Array<OneD, NekDouble> pressure(nPhys), temperature(nPhys);
            Array<OneD, NekDouble> entropy(nPhys);
            Array<OneD, NekDouble> soundspeed(nPhys), mach(nPhys);
            Array<OneD, NekDouble> sensor(nPhys), SensorKappa(nPhys);

            m_varConv->GetVelocityVector(tmp, velocity);
            m_varConv->GetPressure  (tmp, pressure);
            m_varConv->GetTemperature(tmp, temperature);
            m_varConv->GetEntropy   (tmp, entropy);
            m_varConv->GetSoundSpeed(tmp, soundspeed);
            m_varConv->GetMach      (tmp, soundspeed, mach);

            Array<OneD, Array<OneD, NekDouble> > velocities(m_spacedim);
            for (int i=0;i<m_spacedim;i++)
            {
                velocities[i] = Array<OneD, NekDouble> (nPhys);
            }
            m_varConv->GetVelocityVector(tmp,velocities);

            int sensorOffset;
            m_session->LoadParameter ("SensorOffset", sensorOffset, 1);
            m_varConv->GetSensor (m_fields[0], tmp, sensor, SensorKappa,
                                    sensorOffset);

            Array<OneD, NekDouble> pFwd(nCoeffs), TFwd(nCoeffs);
            Array<OneD, NekDouble> sFwd(nCoeffs);
            Array<OneD, NekDouble> aFwd(nCoeffs), mFwd(nCoeffs);
            Array<OneD, NekDouble> sensFwd(nCoeffs);

            string velNames[3] = {"u", "v", "w"};
            for (int i = 0; i < m_spacedim; ++i)
            {
                m_fields[0]->FwdTrans_IterPerExp(velocity[i], velFwd[i]);
                variables.push_back(velNames[i]);
                fieldcoeffs.push_back(velFwd[i]);
            }

            m_fields[0]->FwdTrans_IterPerExp(pressure,   pFwd);
            m_fields[0]->FwdTrans_IterPerExp(temperature,TFwd);
            m_fields[0]->FwdTrans_IterPerExp(entropy,    sFwd);
            m_fields[0]->FwdTrans_IterPerExp(soundspeed, aFwd);
            m_fields[0]->FwdTrans_IterPerExp(mach,       mFwd);
            m_fields[0]->FwdTrans_IterPerExp(sensor,     sensFwd);

            variables.push_back  ("p");
            variables.push_back  ("T");
            variables.push_back  ("s");
            variables.push_back  ("a");
            variables.push_back  ("Mach");
            variables.push_back  ("Sensor");
            fieldcoeffs.push_back(pFwd);
            fieldcoeffs.push_back(TFwd);
            fieldcoeffs.push_back(sFwd);
            fieldcoeffs.push_back(aFwd);
            fieldcoeffs.push_back(mFwd);
            fieldcoeffs.push_back(sensFwd);

            Array<OneD, NekDouble> uFwd(nCoeffs);
            m_fields[0]->FwdTrans_IterPerExp(velocities[0],uFwd);
            variables.push_back  ("u");
            fieldcoeffs.push_back(uFwd);

            if(m_spacedim>1)
            {
                Array<OneD, NekDouble> vFwd(nCoeffs);
                variables.push_back  ("v");
                m_fields[0]->FwdTrans_IterPerExp(velocities[1],vFwd);
                fieldcoeffs.push_back(vFwd);
            }
            if(m_spacedim>2)
            {
                Array<OneD, NekDouble> wFwd(nCoeffs);
                variables.push_back  ("w");
                m_fields[0]->FwdTrans_IterPerExp(velocities[2],wFwd);
                fieldcoeffs.push_back(wFwd);
            }

            if (m_artificialDiffusion)
            {
                // Get min h/p
                m_artificialDiffusion->SetElmtHP(GetElmtMinHP());
                // reuse pressure
                Array<OneD, NekDouble> sensorFwd(nCoeffs);
                m_artificialDiffusion->GetArtificialViscosity(tmp, pressure);
                m_fields[0]->FwdTrans_IterPerExp(pressure,   sensorFwd);

                variables.push_back  ("ArtificialVisc");
                fieldcoeffs.push_back(sensorFwd);
            }
        }
    }

    void CompressibleFlowSystem::v_GetViscousSymmtrFluxConservVar(
            const int                                                       nConvectiveFields,
            const int                                                       nSpaceDim,
            const Array<OneD, Array<OneD, NekDouble> >                      &inaverg,
            const Array<OneD, Array<OneD, NekDouble > >                     &inarray,
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >              &outarray,
            Array< OneD, int >                                              &nonZeroIndex,
            const Array<OneD, Array<OneD, NekDouble> >                      &normals)
    {
        density = physfield[0];
    }
    
    void CompressibleFlowSystem::v_SteadyStateResidual(
                int                         step, 
                Array<OneD, NekDouble>      &L2)
    {
        const int nPoints = GetTotPoints();
        const int nFields = m_fields.num_elements();
        Array<OneD, Array<OneD, NekDouble> > rhs (nFields);
        Array<OneD, Array<OneD, NekDouble> > inarray (nFields);
        for (int i = 0; i < nFields; ++i)
        {
            rhs[i] =   Array<OneD, NekDouble> (nPoints,0.0);
            inarray[i] =   m_fields[i]->UpdatePhys();
        }

        DoOdeRhs(inarray,rhs,m_time);

        // Holds L2 errors.
        Array<OneD, NekDouble> tmp;
        Array<OneD, NekDouble> RHSL2    (nFields);
        Array<OneD, NekDouble> residual(nFields);

        for (int i = 0; i < nFields; ++i)
        {
            tmp = rhs[i];

            Vmath::Vmul(nPoints, tmp, 1, tmp, 1, tmp, 1);
            residual[i] = Vmath::Vsum(nPoints, tmp, 1);
        }

        m_comm->AllReduce(residual , LibUtilities::ReduceSum);

        NekDouble onPoints = 1.0/NekDouble(nPoints);
        for (int i = 0; i < nFields; ++i)
        {
            L2[i] = sqrt(residual[i]*onPoints);
        }
    }

    void CompressibleFlowSystem::v_GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr                            &explist,
            const Array<OneD, const Array<OneD, NekDouble> >                &normals,
            const int                                                       nDervDir,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > > &ElmtJacArray,
            const int                                                       nfluxDir)
    {
        ASSERTL0(false, "v_GetFluxDerivJacDirctn not coded");
    }

    void CompressibleFlowSystem::v_GetFluxDerivJacDirctnElmt(
            const int                                                       nConvectiveFields,
            const int                                                       nElmtPnt,
            const int                                                       nDervDir,
            const Array<OneD, Array<OneD, NekDouble> >                      &locVars,
            const Array<OneD, NekDouble>                                    &locmu,
            const Array<OneD, Array<OneD, NekDouble> >                      &locnormal,
            DNekMatSharedPtr                                                &wspMat,
            Array<OneD, Array<OneD, NekDouble> >                            &PntJacArray)
    {
        ASSERTL0(false, "v_GetFluxDerivJacDirctn not coded");
    }
    
    void CompressibleFlowSystem::v_GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr                            &explist,
            const Array<OneD, const Array<OneD, NekDouble> >                &normals,
            const int                                                       nDervDir,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
                  Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac)
    {
        m_varConv->GetVelocityVector(physfield, velocity);
    }

    // void CompressibleFlowSystem::v_GetFluxDerivJacDirctn(
    //         const MultiRegions::ExpListSharedPtr                            &explist,
    //         const int                                                       nFluxDir,
    //         const int                                                       nDervDir,
    //         const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
    //               Array<OneD, Array<OneD, DNekMatSharedPtr> >               &ElmtJac)
    // {
    //     ASSERTL0(false, "v_GetFluxDerivJacDirctn not coded");
    // }

    void CompressibleFlowSystem::v_GetDiffusionFluxJacPoint(
            const Array<OneD, NekDouble>                        &conservVar, 
            const Array<OneD, const Array<OneD, NekDouble> >    &conseDeriv, 
            const NekDouble                                     mu,
            const NekDouble                                     DmuDT,
            const Array<OneD, NekDouble>                        &normals,
                  DNekMatSharedPtr                              &fluxJac)
    {
        boost::ignore_unused(inarray, outarray, pFwd, pBwd);
        if (m_shockCaptureType != "Off")
        {
            m_artificialDiffusion->DoArtificialDiffusion(inarray, outarray);
        }
    }

    void CompressibleFlowSystem::v_MinusDiffusionFluxJacDirctn(
            const int                                                       nDirctn,
            const Array<OneD, const Array<OneD, NekDouble> >                &inarray,
            const Array<OneD, const Array<OneD, Array<OneD, NekDouble>> >   &qfields,
            Array<OneD, Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > > > &ElmtJacArray)
    {
        ASSERTL0(false, "not coded");
    }

    void CompressibleFlowSystem::v_MinusDiffusionFluxJacDirctnElmt(
            const int                                                       nConvectiveFields,
            const int                                                       nElmtPnt,
            const Array<OneD, Array<OneD, NekDouble> >                      &locVars,
            const Array<OneD, Array<OneD,  Array<OneD, NekDouble> > >       &locDerv,
            const Array<OneD, NekDouble>                                    &locmu,
            const Array<OneD, NekDouble>                                    &locDmuDT,
            const Array<OneD, NekDouble>                                    &normals,
            DNekMatSharedPtr                                                &wspMat,
            Array<OneD, Array<OneD, NekDouble> >                            &PntJacArray)
    {
        boost::ignore_unused(inarray, outarray, pFwd, pBwd);
        // Do nothing by default
    }


/**
 * @brief Compute an estimate of minimum h/p
 * for each element of the expansion.
 */
Array<OneD, NekDouble>  CompressibleFlowSystem::GetElmtMinHP(void)
{
    int nElements               = m_fields[0]->GetExpSize();
    Array<OneD, NekDouble> hOverP(nElements, 1.0);

    // Determine h/p scaling
    Array<OneD, int> pOrderElmt = m_fields[0]->EvalBasisNumModesMaxPerExp();
    for (int e = 0; e < nElements; e++)
    {
        NekDouble h = 1.0e+10;
        switch(m_expdim)
        {
            case 3:
            {
                LocalRegions::Expansion3DSharedPtr exp3D;
                exp3D = m_fields[0]->GetExp(e)->as<LocalRegions::Expansion3D>();
                for(int i = 0; i < exp3D->GetNtraces(); ++i)
                {
                    h = min(h, exp3D->GetGeom3D()->GetEdge(i)->GetVertex(0)->
                        dist(*(exp3D->GetGeom3D()->GetEdge(i)->GetVertex(1))));
                }
            break;
            }

            case 2:
            {
                LocalRegions::Expansion2DSharedPtr exp2D;
                exp2D = m_fields[0]->GetExp(e)->as<LocalRegions::Expansion2D>();
                for(int i = 0; i < exp2D->GetNtraces(); ++i)
                {
                    h = min(h, exp2D->GetGeom2D()->GetEdge(i)->GetVertex(0)->
                        dist(*(exp2D->GetGeom2D()->GetEdge(i)->GetVertex(1))));
                }
            break;
            }
            case 1:
            {
                LocalRegions::Expansion1DSharedPtr exp1D;
                exp1D = m_fields[0]->GetExp(e)->as<LocalRegions::Expansion1D>();

                h = min(h, exp1D->GetGeom1D()->GetVertex(0)->
                    dist(*(exp1D->GetGeom1D()->GetVertex(1))));

            break;
            }
            default:
            {
                ASSERTL0(false,"Dimension out of bound.")
            }
        }

        // Determine h/p scaling
        hOverP[e] = h/max(pOrderElmt[e]-1,1);

    }
    return hOverP;
}
}
