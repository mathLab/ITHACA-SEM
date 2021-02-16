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

#include <LibUtilities/BasicUtils/Timer.h>

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

        for (int i = 0; i < m_fields.size(); i++)
        {
            // Use BwdTrans to make sure initial condition is in solution space
            m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                  m_fields[i]->UpdatePhys());
        }

        m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
                    m_session, m_spacedim);

        ASSERTL0(m_session->DefinesSolverInfo("UPWINDTYPE"),
                 "No UPWINDTYPE defined in session.");

        // Do not forwards transform initial condition
        m_homoInitialFwd = false;

        // Set up locations of velocity vector.
        m_vecLocs = TensorOfArray2D<NekDouble>(1);
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
            m_ode.DefineProjection(
                &CompressibleFlowSystem::DoOdeProjection, this);
        }
        else
        {
            m_ode.DefineOdeRhs    (&CompressibleFlowSystem::DoOdeRhs, this);
            m_ode.DefineProjection(&CompressibleFlowSystem::
                                   DoOdeProjection, this);
            m_ode.DefineImplicitSolve(&CompressibleFlowSystem::
                                      DoImplicitSolvePhysToCoeff, this);
            InitialiseNonlinSysSolver();

            int nvariables  =   m_fields.size();
            Array<OneD, Array<OneD, Array<OneD, int > > >   map;
            const MultiRegions::LocTraceToTraceMapSharedPtr locTraceToTraceMap =
                m_fields[0]->GetLocTraceToTraceMap();

            locTraceToTraceMap->CalcLocTracePhysToTraceIDMap(
                m_fields[0]->GetTrace(),m_spacedim);
            for(int i=1;i<nvariables;i++)
            {
                m_fields[i]->GetLocTraceToTraceMap()
                    ->SetLocTracephysToTraceIDMap(
                    locTraceToTraceMap->GetLocTracephysToTraceIDMap());
            }
        }

        SetBoundaryConditionsBwdWeight();
    }

    void CompressibleFlowSystem::InitialiseNonlinSysSolver()
    {
        std::string SolverType = "Newton";
        if (m_session->DefinesSolverInfo("NonlinSysIterSovler"))
        {
            SolverType = m_session->GetSolverInfo("NonlinSysIterSovler");
        }
        ASSERTL0(LibUtilities::GetNekNonlinSysFactory().
            ModuleExists(SolverType), "NekNonlinSys '" + SolverType + 
            "' is not defined.\n");
        int ntotal = m_fields[0]->GetNcoeffs() * m_fields.size();

        LibUtilities::NekSysKey key = LibUtilities::NekSysKey();

        key.m_NonlinIterTolRelativeL2   = 1.0E-3;
        key.m_LinSysRelativeTolInNonlin = 5.0E-2;
        key.m_NekNonlinSysMaxIterations = 10;
        key.m_NekLinSysMaxIterations    = 30;
        key.m_LinSysMaxStorage          = 30;

        m_nonlinsol = LibUtilities::GetNekNonlinSysFactory().CreateInstance(
            SolverType, m_session, m_comm, ntotal, key);

        m_NekSysOp.DefineNekSysResEval(&CompressibleFlowSystem::
            NonlinSysEvaluatorCoeff1D, this);
        m_NekSysOp.DefineNekSysLhsEval(&CompressibleFlowSystem::
            MatrixMultiplyMatrixFreeCoeff, this);

        m_session->LoadParameter("PreconMatFreezNumb",     
            m_PreconMatFreezNumb    , 1);
        m_session->LoadParameter("NewtonAbsoluteIteTol", 
            m_NewtonAbsoluteIteTol, 1.0E-12);
        m_session->LoadParameter("JFNKTimeAccurate",     
            m_JFNKTimeAccurate, 1);
        m_session->LoadParameter("JFNKPrecondStep",      
            m_JFNKPrecondStep, 7);
        m_session->LoadParameter("SORRelaxParam",        
            m_SORRelaxParam, 1.0);

        // when no time accuracy needed
        if(m_JFNKTimeAccurate < 1)
        {
            m_NewtonAbsoluteIteTol = 1.0E-10;
        }

        if (boost::iequals(m_session->GetSolverInfo("PRECONDITIONER"),
                               "IncompleteLU"))
        {
            m_PreconMatStorage    =   eSparse;

            NEKERROR(ErrorUtil::efatal,
                "IncompleteLU preconditioner not finished yet");

        }
        else
        {
            int nvariables  =   m_fields.size();
            m_NekSysOp.DefineNekSysPrecond(
                    &CompressibleFlowSystem::PrecBlkSORCoeff, this);
            m_PreconMatStorage    =   eDiagonal;
            m_session->LoadParameter("nPadding", m_nPadding, 4);
    
            m_PreconMatVarsSingle = 
                TensorOfArray2D<SNekBlkMatSharedPtr>(nvariables);
            for(int i = 0; i < nvariables; i++)
            {
                m_PreconMatVarsSingle[i] =  
                    Array<OneD, SNekBlkMatSharedPtr> (nvariables);
            }
            AllocatePrecBlkDiagCoeff(m_PreconMatVarsSingle);

            int nelmts  = m_fields[0]->GetNumElmts();
            int nelmtcoef;
            Array<OneD, unsigned int > nelmtmatdim(nelmts);
            for(int i = 0; i < nelmts; i++)
            {
                nelmtcoef   =   m_fields[0]->GetExp(i)->GetNcoeffs();
                nelmtmatdim[i]  =   nelmtcoef*nvariables;
            }
            AllocateNekBlkMatDig(m_PreconMatSingle,nelmtmatdim,nelmtmatdim);
        }

        m_nonlinsol->SetSysOperators(m_NekSysOp);
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
        m_session->LoadParameter ("JacobiFreeEps", m_JacobiFreeEps, 5.0E-8);

        int ntmp;
        m_session->LoadParameter("AdvectionJacFlag", ntmp, 1);
        m_AdvectionJacFlag             = true;
        if(0==ntmp)
        {
            m_AdvectionJacFlag = false;
        }
        m_session->LoadParameter("ViscousJacFlag", ntmp, 1);
        m_ViscousJacFlag             = true;
        if(0==ntmp)
        {
            m_ViscousJacFlag = false;
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
              TensorOfArray2D<NekDouble> &outarray,
        const NekDouble                                   time)
    {
        int nvariables = inarray.size();
        int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        m_BndEvaluateTime   = time;

        // Store forwards/backwards space along trace space
        TensorOfArray2D<NekDouble> Fwd    (nvariables);
        TensorOfArray2D<NekDouble> Bwd    (nvariables);

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

        // Calculate advection
        LibUtilities::Timer timer;
        timer.Start();
                DoAdvection(inarray, outarray, time, Fwd, Bwd);
        timer.Stop();
        timer.AccumulateRegion("DoAdvection");

        // Negate results
        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(npoints, outarray[i], 1);
        }

        // Add diffusion terms
timer.Start();
        DoDiffusion(inarray, outarray, Fwd, Bwd);
timer.Stop();
timer.AccumulateRegion("DoDiffusion");

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
            for (int n = 0; n < nElements; ++n)
            {
                nq     = m_fields[0]->GetExp(n)->GetTotPoints();
                offset = m_fields[0]->GetPhys_Offset(n);
                fac    = tstep[n] / m_timestep;
                for (int i = 0; i < nvariables; ++i)
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
              TensorOfArray2D<NekDouble> &outarray,
        const NekDouble                                   time)
    {
        int nvariables = inarray.size();

        switch(m_projectionType)
        {
            case MultiRegions::eDiscontinuous:
            {
                // Just copy over array
                int npoints = GetNpoints();

                for (int i = 0; i < nvariables; ++i)
                {
                    Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
                    if (m_useFiltering)
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
                NEKERROR(ErrorUtil::efatal, "No Continuous Galerkin for full "
                    "compressible Navier-Stokes equations");
                break;
            }
            default:
                NEKERROR(ErrorUtil::efatal, "Unknown projection scheme");
                break;
        }
    }

    void CompressibleFlowSystem::Prec(
        const Array<OneD, NekDouble> &inarray,
        Array<OneD, NekDouble >&out)
    {
        int ntotal = inarray.size();
        Vmath::Vcopy(ntotal, inarray, 1, out, 1);
        return;
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::PrecBlkDiag(
        const Array<OneD, NekDouble>  &inarray,
        Array<OneD, NekDouble >       &outarray,
        const TypeNekBlkMatSharedPtr  &PreconMatVars,
        const DataType                &tmpDataType)
    {
        boost::ignore_unused(tmpDataType);

        unsigned int nvariables = m_fields.size();
        unsigned int npoints    = GetNcoeffs();
        unsigned int npointsVar = nvariables*npoints;
        Array<OneD, DataType >Sinarray(npointsVar);
        Array<OneD, DataType > Soutarray(npointsVar);
        NekVector<DataType> tmpVect(npointsVar,Sinarray,eWrapper);
        NekVector<DataType> outVect(npointsVar,Soutarray,eWrapper);

        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    
            m_fields[0]->GetExp();
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

        outVect = (*PreconMatVars)*tmpVect;

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

    void CompressibleFlowSystem::PrecBlkSORCoeff(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble >&outarray,
            const bool                   &flag)
    {

        boost::ignore_unused(flag);
        if (m_updatePreconMatFlag)
        {
            CalcPreconMat(m_solutionPhys, m_BndEvaluateTime, m_TimeIntegLambda);

            m_TimeIntegLambdaPrcMat = m_TimeIntegLambda;

            m_CalcPreconMatNumbSteps += m_CalcPreconMatCounter;
            m_CalcPreconMatNumbSteps = m_CalcPreconMatNumbSteps/2;
            m_CalcPreconMatCounter   = 1;

            m_CalcPreconMatFlag = false;

            m_updatePreconMatFlag = false;

            if (m_verbose&&m_root)
            {
                cout << "     ## CalcuPreconMat " << endl;
            }
        }


        int nSORTot   =   m_JFNKPrecondStep;
        if (0==nSORTot)
        {
            Prec(inarray,outarray);
        }
        else
        {
            const NekDouble SORParam        =   m_SORRelaxParam;
            const NekDouble OmSORParam      =   1.0-SORParam;

            unsigned int nvariables = m_fields.size();
            unsigned int npoints    = GetNcoeffs();
            unsigned int ntotpnt    = inarray.size();
            
            ASSERTL0(nvariables*npoints==ntotpnt,
                "nvariables*npoints!=ntotpnt in PrecBlkSORCoeff");

            Array<OneD, NekDouble> rhs(ntotpnt);

            Array<OneD, NekDouble>  outN(ntotpnt);
            Array<OneD, NekDouble>  outTmp(ntotpnt);
            TensorOfArray2D<NekDouble>rhs2d(nvariables);
            TensorOfArray2D<NekDouble>out_2d(nvariables);
            TensorOfArray2D<NekDouble>outTmp_2d(nvariables);
            for(int m = 0; m < nvariables; m++)
            {
                int moffset     = m*npoints;
                rhs2d[m]        = rhs       + moffset;
                out_2d[m]       = outarray  + moffset;
                outTmp_2d[m]    = outTmp    + moffset;
                m_fields[m]->MultiplyByMassMatrix(inarray+moffset,rhs2d[m]);
            }

            int nphysic    = GetNpoints();
            int nTracePts  = GetTraceTotPoints();
            TensorOfArray3D<NekDouble> qfield(m_spacedim);
            for(int i = 0; i< m_spacedim; i++)
            {
                qfield[i]   =   
                    TensorOfArray2D<NekDouble>(nvariables);
                for(int j = 0; j< nvariables; j++)
                {
                    qfield[i][j]   =   Array<OneD, NekDouble>(nphysic);
                }
            }
            int ntmpTrace = 4+2*m_spacedim;
            TensorOfArray3D<NekDouble> tmpTrace(ntmpTrace);
            for(int i = 0; i< ntmpTrace; i++)
            {
                tmpTrace[i]   =   
                    TensorOfArray2D<NekDouble>(nvariables);
                for(int j = 0; j< nvariables; j++)
                {
                    tmpTrace[i][j]   =   Array<OneD, NekDouble>(nTracePts);
                }
            }
            TensorOfArray2D<NekDouble> FwdFluxDeriv(nvariables);
            TensorOfArray2D<NekDouble> BwdFluxDeriv(nvariables);
            for(int j = 0; j< nvariables; j++)
            {
                FwdFluxDeriv[j]   =   Array<OneD, NekDouble>(nTracePts);
                BwdFluxDeriv[j]   =   Array<OneD, NekDouble>(nTracePts);
            }

            bool flagUpdateDervFlux = false;

            const int nwspTraceDataType = nvariables+1;
            NekSingle tmpSingle;
            TensorOfArray2D<NekSingle> wspTraceDataType(nwspTraceDataType);
            for(int m=0;m<nwspTraceDataType;m++)
            {
                wspTraceDataType[m] =   Array<OneD, NekSingle>(nTracePts);
            }

            PrecBlkDiag(rhs,outarray,m_PreconMatSingle,tmpSingle);

            for(int nsor = 0; nsor < nSORTot-1; nsor++)
            {
                Vmath::Smul(ntotpnt,OmSORParam,outarray,1,outN,1);
                
                MinusOffDiag2Rhs(nvariables,npoints,rhs2d,out_2d,
                    flagUpdateDervFlux,FwdFluxDeriv,BwdFluxDeriv,qfield,
                    tmpTrace,wspTraceDataType,m_TraceJacArraySingle, 
                    m_TraceJacDerivArraySingle, m_TraceJacDerivSignSingle,
                    m_TraceIPSymJacArraySingle);

                PrecBlkDiag(outarray,outTmp,m_PreconMatSingle,
                    tmpSingle);
                Vmath::Svtvp(ntotpnt,SORParam,outTmp,1,outN,1,outarray,1);
            }
        }
    }

    void CompressibleFlowSystem::CalcPreconMat(
        const TensorOfArray2D<NekDouble>    &inpnts,
        const NekDouble                     time,
        const NekDouble                     lambda)
    {
        boost::ignore_unused(time, lambda);
        int nvariables = inpnts.size();

        int nphspnt = inpnts[0].size();
        TensorOfArray2D<NekDouble> intmp(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            intmp[i]    =   Array<OneD, NekDouble>(nphspnt,0.0);
        }

        DoOdeProjection(inpnts,intmp,m_BndEvaluateTime);
        GetPrecNSBlkDiagCoeff(intmp, m_PreconMatVarsSingle, 
            m_PreconMatSingle, m_TraceJacSingle, m_TraceJacDerivSingle,
            m_TraceJacDerivSignSingle, m_TraceJacArraySingle, 
            m_TraceJacDerivArraySingle, m_TraceIPSymJacArraySingle,
            m_StdSMatDataDBB,m_StdSMatDataDBDB);
            
    
        m_CalcPreconMatFlag = false;
        m_TimeIntegLambdaPrcMat = m_TimeIntegLambda;

        // to free the storage
        for(int i = 0; i < nvariables; i++)
        {
            intmp[i]    =   NullNekDouble1DArray;
        }
    }

    template<typename DataType>
    void CompressibleFlowSystem::MinusOffDiag2Rhs(
        const int                                       nvariables,
        const int                                       nCoeffs,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        TensorOfArray2D<NekDouble>                      &outarray,
        bool                                            flagUpdateDervFlux,
        TensorOfArray2D<NekDouble>                      &FwdFluxDeriv,
        TensorOfArray2D<NekDouble>                      &BwdFluxDeriv,
        TensorOfArray3D<NekDouble>                      &qfield,
        TensorOfArray3D<NekDouble>                      &wspTrace,
        TensorOfArray2D<DataType>                       &wspTraceDataType,
        const TensorOfArray4D<DataType>                 &TraceJacArray,
        const TensorOfArray4D<DataType>                 &TraceJacDerivArray,
        const TensorOfArray2D<DataType>                 &TraceJacDerivSign,
        const TensorOfArray5D<DataType>                 &TraceIPSymJacArray)
    {
        boost::ignore_unused(flagUpdateDervFlux, qfield, TraceJacDerivArray, 
            TraceJacDerivSign, FwdFluxDeriv, BwdFluxDeriv, TraceIPSymJacArray);

        int nTracePts  = GetTraceTotPoints();
        int npoints    = GetNpoints();
        int nDim       = m_spacedim;

        TensorOfArray2D<NekDouble> outpnts(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            outpnts[i]  =  Array<OneD, NekDouble> (npoints,0.0);
            m_fields[i]->BwdTrans(outarray[i],outpnts[i]);
        }

        // Store forwards/backwards space along trace space
        TensorOfArray2D<NekDouble> Fwd    ;
        TensorOfArray2D<NekDouble> Bwd    ;
        TensorOfArray2D<NekDouble> FwdFlux;
        TensorOfArray2D<NekDouble> BwdFlux;
        TensorOfArray3D<NekDouble>    numDerivBwd(nDim);
        TensorOfArray3D<NekDouble>    numDerivFwd(nDim);
        int indexwspTrace = 0;
        Fwd     =   wspTrace[indexwspTrace], indexwspTrace++;
        Bwd     =   wspTrace[indexwspTrace], indexwspTrace++;
        FwdFlux =   wspTrace[indexwspTrace], indexwspTrace++;
        BwdFlux =   wspTrace[indexwspTrace], indexwspTrace++;
       
        for(int i = 0; i < nvariables; ++i)
        {
            m_fields[i]->GetFwdBwdTracePhys(outpnts[i], Fwd[i], Bwd[i]);
        }

        int indexwspTraceDataType = 0;
        TensorOfArray2D<DataType> Fwdarray (nvariables);
        for(int m = 0; m < nvariables; ++m)
        {
            Fwdarray[m] = wspTraceDataType[indexwspTraceDataType], indexwspTraceDataType++;
        }
        Array<OneD, DataType> Fwdreslt;
        Fwdreslt = wspTraceDataType[indexwspTraceDataType], indexwspTraceDataType++;

        for(int m = 0; m < nvariables; ++m)
        {
            for(int i = 0; i < nTracePts; ++i)
            {
                Fwdarray[m][i] =  DataType( Fwd[m][i] );
            }
        }
        for(int m = 0; m < nvariables; ++m)
        {
            Vmath::Zero(nTracePts, Fwdreslt,1);
            for(int n = 0; n < nvariables; ++n)
            {
                Vmath::Vvtvp(nTracePts, TraceJacArray[0][m][n], 1, 
                    Fwdarray[n], 1, Fwdreslt, 1, Fwdreslt, 1);
            }

            for(int i = 0; i < nTracePts; ++i)
            {
                FwdFlux[m][i] =  NekDouble( Fwdreslt[i] );
            }
        }

        for(int m = 0; m < nvariables; ++m)
        {
            for(int i = 0; i < nTracePts; ++i)
            {
                Fwdarray[m][i] =  DataType( Bwd[m][i] );
            }
        }
        for(int m = 0; m < nvariables; ++m)
        {
            Vmath::Zero(nTracePts, Fwdreslt,1);
            for(int n = 0; n < nvariables; ++n)
            {
                Vmath::Vvtvp(nTracePts,TraceJacArray[1][m][n],1,Fwdarray[n],1,
                    Fwdreslt,1,Fwdreslt,1);
            }
            for(int i = 0; i < nTracePts; ++i)
            {
                BwdFlux[m][i] =  NekDouble( Fwdreslt[i] );
            }
        }

        for(int i = 0; i < nvariables; ++i)
        {
            Vmath::Fill(nCoeffs,0.0,outarray[i],1);
            m_fields[i]->AddTraceIntegralToOffDiag(FwdFlux[i],BwdFlux[i], 
                outarray[i]);
        }

        for(int i = 0; i < nvariables; ++i)
        {
            Vmath::Svtvp(nCoeffs,-m_TimeIntegLambda,outarray[i],1,inarray[i],1,
                outarray[i],1);
        }
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::ElmtVarInvMtrx(
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        TypeNekBlkMatSharedPtr                  &gmtVar,
        const DataType                          &tmpDataType)
    {
        boost::ignore_unused(tmpDataType);

        int n1d = gmtxarray.size();
        int n2d = gmtxarray[0].size();
        int nConvectiveFields = n1d;

        ASSERTL0(n1d==n2d,"ElmtVarInvMtrx requires n1d==n2d");

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;

        gmtxarray[0][0]->GetBlockSizes(rowSizes,colSizes);
        int ntotElmt  = rowSizes.size();
        int nElmtCoef   =    rowSizes[0]-1;
        int nElmtCoef0  =    -1;
        int blocksize = -1;

        Array<OneD, unsigned int> tmprow(1);
        TypeNekBlkMatSharedPtr tmpGmtx;

        Array<OneD, DataType>    GMatData,ElmtMatData;
        Array<OneD, DataType>    tmpArray1,tmpArray2;

        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            int nrows = gmtxarray[0][0]->GetBlock(nelmt,nelmt)->GetRows();
            int ncols = gmtxarray[0][0]->GetBlock(nelmt,nelmt)->GetColumns();
            ASSERTL0(nrows==ncols,"ElmtVarInvMtrx requires nrows==ncols");

            nElmtCoef            = nrows;

            if (nElmtCoef0!=nElmtCoef)
            {
                nElmtCoef0 = nElmtCoef;
                int nElmtCoefVAr = nElmtCoef0*nConvectiveFields;
                blocksize = nElmtCoefVAr*nElmtCoefVAr;
                tmprow[0] = nElmtCoefVAr;
                AllocateNekBlkMatDig(tmpGmtx,tmprow,tmprow);
                GMatData = tmpGmtx->GetBlock(0,0)->GetPtr();
            }

            for(int n = 0; n < nConvectiveFields; n++)
            {
                for(int m = 0; m < nConvectiveFields; m++)
                {
                    ElmtMatData = gmtxarray[m][n]->
                        GetBlock(nelmt,nelmt)->GetPtr();

                    for(int ncl = 0; ncl < nElmtCoef; ncl++)
                    {
                        int Goffset = (n*nElmtCoef+ncl)*nConvectiveFields*
                            nElmtCoef+m*nElmtCoef;
                        int Eoffset = ncl*nElmtCoef;

                        Vmath::Vcopy(nElmtCoef,
                            tmpArray1 = ElmtMatData+Eoffset,1, 
                            tmpArray2 = GMatData+Goffset,1);
                    }
                }
            }

            tmpGmtx->GetBlock(0,0)->Invert();

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    ElmtMatData = 
                        gmtxarray[m][n]->GetBlock(nelmt,nelmt)->GetPtr();

                    for(int ncl = 0; ncl < nElmtCoef; ncl++)
                    {
                        int Goffset = (n*nElmtCoef+ncl)*nConvectiveFields*
                            nElmtCoef+m*nElmtCoef;
                        int Eoffset = ncl*nElmtCoef;

                        Vmath::Vcopy(nElmtCoef, 
                            tmpArray1 = GMatData+Goffset,1,
                            tmpArray2 = ElmtMatData+Eoffset,1);
                    }
                }
            }
            ElmtMatData = gmtVar->GetBlock(nelmt,nelmt)->GetPtr();
            Vmath::Vcopy(blocksize, GMatData,1,ElmtMatData,1);
        }
        return;
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::AddMatNSBlkDiag_volume(
        const Array<OneD, const Array<OneD, NekDouble>>     &inarray,
        const Array<OneD, const TensorOfArray2D<NekDouble>> &qfield,
        TensorOfArray2D<TypeNekBlkMatSharedPtr>             &gmtxarray,
        TensorOfArray4D<DataType>                           &StdMatDataDBB,
        TensorOfArray5D<DataType>                           &StdMatDataDBDB)
    {
        if (StdMatDataDBB.size() == 0)
        {
            CalcVolJacStdMat(StdMatDataDBB,StdMatDataDBDB);
        }
        
        int nSpaceDim = m_graph->GetSpaceDimension();
        int nvariable = inarray.size();
        int npoints   = m_fields[0]->GetTotPoints();
        int nVar2     = nvariable*nvariable;
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    
            m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();

        Array<OneD, NekDouble > mu      (npoints, 0.0);
        Array<OneD, NekDouble > DmuDT   (npoints, 0.0);
        if(m_ViscousJacFlag)
        {
            CalcMuDmuDT(inarray,mu,DmuDT);
        }

        Array<OneD, NekDouble> normals;
        TensorOfArray2D<NekDouble> normal3D(3);
        for(int i = 0; i < 3; i++)
        {
            normal3D[i] = Array<OneD, NekDouble>(3,0.0);
        }
        normal3D[0][0] = 1.0;
        normal3D[1][1] = 1.0;
        normal3D[2][2] = 1.0;
        TensorOfArray2D<NekDouble> normalPnt(3);
        
        DNekMatSharedPtr wspMat     = 
            MemoryManager<DNekMat>::AllocateSharedPtr(nvariable,nvariable,0.0);
        DNekMatSharedPtr wspMatDrv  = 
            MemoryManager<DNekMat>::AllocateSharedPtr(nvariable-1,nvariable,0.0);

        Array<OneD, DataType> GmatxData;
        Array<OneD, DataType> MatData;

        Array<OneD, NekDouble> tmppnts;
        TensorOfArray3D<NekDouble> PntJacCons(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray3D<DataType>  PntJacConsStd(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray2D<NekDouble> ConsStdd(m_spacedim);
        TensorOfArray2D<NekDouble> ConsCurv(m_spacedim);
        TensorOfArray4D<NekDouble> PntJacDerv(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray4D<DataType>  PntJacDervStd(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray3D<NekDouble> DervStdd(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray3D<NekDouble> DervCurv(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        for(int ndir=0; ndir<m_spacedim;ndir++)
        {
            PntJacDerv[ndir]  =   TensorOfArray3D<NekDouble>(m_spacedim);
            PntJacDervStd[ndir] = TensorOfArray3D<DataType>(m_spacedim);
            DervStdd[ndir]    =   TensorOfArray2D<NekDouble>(m_spacedim);
            DervCurv[ndir]    =   TensorOfArray2D<NekDouble>(m_spacedim);
        }

        Array<OneD, NekDouble> locmu;
        Array<OneD, NekDouble> locDmuDT;
        TensorOfArray2D<NekDouble> locVars(nvariable);
        TensorOfArray3D<NekDouble> locDerv(m_spacedim);
        for(int ndir=0; ndir<m_spacedim;ndir++)
        {
            locDerv[ndir] = TensorOfArray2D<NekDouble>(nvariable);
        }

        int nElmtCoefOld = -1;
        for(int ne=0; ne<ntotElmt;ne++)
        {
            int nElmtCoef           = (*expvect)[ne]->GetNcoeffs();
            int nElmtCoef2          = nElmtCoef*nElmtCoef;
            int nElmtPnt            = (*expvect)[ne]->GetTotPoints();

            int nQuot = nElmtCoef2/m_nPadding;
            int nRemd = nElmtCoef2- nQuot*m_nPadding;
            int nQuotPlus=nQuot;
            if(nRemd>0)
            {
                nQuotPlus++;
            }
            int nElmtCoef2Paded = nQuotPlus*m_nPadding;

            if(nElmtPnt>PntJacCons[0].size()||nElmtCoef>nElmtCoefOld)
            {
                nElmtCoefOld = nElmtCoef;
                for(int ndir=0; ndir<3;ndir++)
                {
                    normalPnt[ndir]        = Array<OneD, NekDouble>(npoints,0.0);
                }
                tmppnts = Array<OneD, NekDouble>  (nElmtPnt);
                MatData = Array<OneD, DataType>  (nElmtCoef2Paded*nVar2);
                for(int ndir=0; ndir<m_spacedim;ndir++)
                {
                    ConsCurv[ndir] = Array<OneD, NekDouble> (nElmtPnt);
                    ConsStdd[ndir] = Array<OneD, NekDouble> (nElmtPnt);
                    PntJacCons[ndir] = TensorOfArray2D<NekDouble> (nElmtPnt);
                    PntJacConsStd[ndir] = TensorOfArray2D<DataType> (nElmtPnt);
                    for(int i=0; i<nElmtPnt;i++)
                    {
                        PntJacCons[ndir][i] = Array<OneD, NekDouble>(nVar2);
                        PntJacConsStd[ndir][i] = Array<OneD, DataType>(nVar2);
                    }
                    
                    for(int ndir1=0; ndir1<m_spacedim;ndir1++)
                    {
                        PntJacDerv[ndir][ndir1] = TensorOfArray2D<NekDouble> (nElmtPnt);
                        PntJacDervStd[ndir][ndir1] = TensorOfArray2D<DataType> (nElmtPnt);
                        DervStdd[ndir][ndir1] = Array<OneD, NekDouble> (nElmtPnt);
                        DervCurv[ndir][ndir1] = Array<OneD, NekDouble> (nElmtPnt);
                        for(int i=0; i<nElmtPnt;i++)
                        {
                            PntJacDerv[ndir][ndir1][i] = 
                                Array<OneD, NekDouble>(nVar2);
                            PntJacDervStd[ndir][ndir1][i] = 
                                Array<OneD, DataType>(nVar2);
                        }
                    }
                }
            }
    
            int noffset = GetPhys_Offset(ne);
            for(int j = 0; j < nvariable; j++)
            {   
                locVars[j] = inarray[j]+noffset;
            }
            
            if(m_AdvectionJacFlag)
            {
                for(int nfluxDir = 0; nfluxDir < nSpaceDim; nfluxDir++)
                {
                    normals =   normal3D[nfluxDir];
                    GetFluxVectorJacDirctnElmt(nvariable,nElmtPnt,locVars,
                        normals,wspMat,PntJacCons[nfluxDir]);
                }
            }

            if(m_ViscousJacFlag)
            {
                for(int j = 0; j < nSpaceDim; j++)
                {   
                    for(int k = 0; k < nvariable; k++)
                    {
                        locDerv[j][k] = qfield[j][k]+noffset;
                    }
                }
                locmu       =   mu      + noffset;
                locDmuDT    =   DmuDT   + noffset;
                for(int nfluxDir = 0; nfluxDir < nSpaceDim; nfluxDir++)
                {
                    normals =   normal3D[nfluxDir];
                    MinusDiffusionFluxJacDirctnElmt(nvariable,nElmtPnt,
                        locVars,locDerv,locmu,locDmuDT,normals,wspMatDrv,
                        PntJacCons[nfluxDir]);
                }
            }

            if(m_ViscousJacFlag)
            {
                locmu = mu + noffset;
                for(int nfluxDir = 0; nfluxDir < nSpaceDim; nfluxDir++)
                {
                    Vmath::Fill(npoints,1.0,normalPnt[nfluxDir],1);
                    for(int nDervDir = 0; nDervDir < nSpaceDim; nDervDir++)
                    {
                        GetFluxDerivJacDirctnElmt(nvariable,nElmtPnt,nDervDir,
                            locVars,locmu,normalPnt,wspMatDrv,
                            PntJacDerv[nfluxDir][nDervDir]);
                    }
                    Vmath::Fill(npoints,0.0,normalPnt[nfluxDir],1);
                }
            }

            for(int n=0; n<nvariable;n++)
            {
                for(int m=0; m<nvariable;m++)
                {
                    int nvarOffset = m+n*nvariable;
                    GmatxData = gmtxarray[m][n]->GetBlock(ne,ne)->GetPtr();

                    for(int ndStd0 =0;ndStd0<m_spacedim;ndStd0++)
                    {
                        Vmath::Zero(nElmtPnt,ConsStdd[ndStd0],1);
                    }
                    for(int ndir =0;ndir<m_spacedim;ndir++)
                    {
                        for(int i=0; i<nElmtPnt;i++)
                        {
                            tmppnts[i] =  PntJacCons[ndir][i][nvarOffset];
                        }
                        (*expvect)[ne]->ProjectVectorintoStandardExp(ndir,
                            tmppnts,ConsCurv);
                        for(int nd =0;nd<m_spacedim;nd++)
                        {
                            Vmath::Vadd(nElmtPnt,ConsCurv[nd],1,ConsStdd[nd],1,
                                ConsStdd[nd],1);
                        }
                    }

                    for(int ndir =0;ndir<m_spacedim;ndir++)
                    {
                        (*expvect)[ne]->MultiplyByQuadratureMetric(
                            ConsStdd[ndir],ConsStdd[ndir]); // weight with metric
                        for(int i=0; i<nElmtPnt;i++)
                        {
                            PntJacConsStd[ndir][i][nvarOffset] = 
                                DataType(ConsStdd[ndir][i]);
                        }
                    }
                }
            }

            if(m_ViscousJacFlag)
            {
                for(int m=0; m<nvariable;m++)
                {
                    for(int n=0; n<nvariable;n++)
                    {
                        int nvarOffset = m+n*nvariable;
                        for(int ndStd0 =0;ndStd0<m_spacedim;ndStd0++)
                        {
                            for(int ndStd1 =0;ndStd1<m_spacedim;ndStd1++)
                            {
                                Vmath::Zero(nElmtPnt,
                                    DervStdd[ndStd0][ndStd1],1);
                            }
                        }
                        for(int nd0 =0;nd0<m_spacedim;nd0++)
                        {
                            for(int nd1 =0;nd1<m_spacedim;nd1++)
                            {
                                for(int i=0; i<nElmtPnt;i++)
                                {
                                    tmppnts[i] =  
                                        PntJacDerv[nd0][nd1][i][nvarOffset];
                                }

                                (*expvect)[ne]->ProjectVectorintoStandardExp(
                                    nd0,tmppnts,ConsCurv);
                                for(int nd =0;nd<m_spacedim;nd++)
                                {
                                    (*expvect)[ne]->
                                        ProjectVectorintoStandardExp(nd1,
                                            ConsCurv[nd],DervCurv[nd]);
                                }

                                for(int ndStd0 =0;ndStd0<m_spacedim;ndStd0++)
                                {
                                    for(int ndStd1 =0;ndStd1<m_spacedim;ndStd1++)
                                    {
                                        Vmath::Vadd(nElmtPnt,
                                            DervCurv[ndStd0][ndStd1],1,
                                            DervStdd[ndStd0][ndStd1],1,
                                            DervStdd[ndStd0][ndStd1],1);
                                    }
                                }
                            }
                        }
                        for(int nd0 =0;nd0<m_spacedim;nd0++)
                        {
                            for(int nd1 =0;nd1<m_spacedim;nd1++)
                            {
                                (*expvect)[ne]->
                                    MultiplyByQuadratureMetric(
                                        DervStdd[nd0][nd1],
                                        DervStdd[nd0][nd1]); // weight with metric
                                for(int i=0; i<nElmtPnt;i++)
                                {
                                    PntJacDervStd[nd0][nd1][i][nvarOffset] = 
                                        -DataType(DervStdd[nd0][nd1][i]);
                                }
                            }
                        }
                    }
                }
            }
            
            Vmath::Zero(nElmtCoef2Paded*nVar2,MatData,1);
            DataType one = 1.0;
            for(int ndir =0;ndir<m_spacedim;ndir++)
            {
                for(int i=0;i<nElmtPnt;i++)
                {
                    Blas::Ger (nElmtCoef2Paded,nVar2,one,
                                &StdMatDataDBB[ne][ndir][i][0],1,
                                &PntJacConsStd[ndir][i][0],1,
                                &MatData[0],nElmtCoef2Paded);
                }
            }

            if(m_ViscousJacFlag)
            {
                for(int nd0 =0;nd0<m_spacedim;nd0++)
                {
                    for(int nd1 =0;nd1<m_spacedim;nd1++)
                    {
                        for(int i=0;i<nElmtPnt;i++)
                        {
                            Blas::Ger (nElmtCoef2Paded,nVar2,one,
                                        &StdMatDataDBDB[ne][nd0][nd1][i][0],1,
                                        &PntJacDervStd[nd0][nd1][i][0],1,
                                        &MatData[0],nElmtCoef2Paded);
                        }
                    }
                }
            }


            Array<OneD, DataType> tmpA;

            for(int n=0; n<nvariable;n++)
            {
                for(int m=0; m<nvariable;m++)
                {
                    int nvarOffset = m+n*nvariable;
                    GmatxData = gmtxarray[m][n]->GetBlock(ne,ne)->GetPtr();
                    Vmath::Vcopy(nElmtCoef2, 
                        tmpA = MatData + nvarOffset*nElmtCoef2Paded,1,
                        GmatxData,1);
                }
            }
        }
    }

    template<typename DataType>
    void CompressibleFlowSystem::CalcVolJacStdMat(
        TensorOfArray4D<DataType> &StdMatDataDBB,
        TensorOfArray5D<DataType> &StdMatDataDBDB)
    {
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    
            m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();

        StdMatDataDBB = TensorOfArray4D<DataType> (ntotElmt);
        StdMatDataDBDB  = TensorOfArray5D<DataType> (ntotElmt);

        vector<DNekMatSharedPtr> VectStdDerivBase0;
        vector< TensorOfArray3D<DataType> > VectStdDerivBase_Base;
        vector< TensorOfArray4D<DataType> > VectStdDervBase_DervBase;
        DNekMatSharedPtr MatStdDerivBase0;
        Array<OneD, DNekMatSharedPtr> ArrayStdMat(m_spacedim);
        TensorOfArray2D<NekDouble>    ArrayStdMatData(m_spacedim);
        for(int ne=0; ne<ntotElmt;ne++)
        {
            StdRegions::StdExpansionSharedPtr stdExp;
            stdExp = (*expvect)[ne]->GetStdExp();
            StdRegions::StdMatrixKey  matkey(StdRegions::eDerivBase0,
                                stdExp->DetShapeType(), *stdExp);
            MatStdDerivBase0      =   stdExp->GetStdMatrix(matkey);

            int ntotStdExp = VectStdDerivBase0.size();
            int nfoundStdExp = -1;
            for(int i=0;i<ntotStdExp;i++) 
            {
                if((*VectStdDerivBase0[i])==(*MatStdDerivBase0))
                {
                    nfoundStdExp = i;
                }
            }
            if(nfoundStdExp>=0)
            {
                StdMatDataDBB[ne] = VectStdDerivBase_Base[nfoundStdExp];
                StdMatDataDBDB[ne] = VectStdDervBase_DervBase[nfoundStdExp];
            }
            else
            {
                int nElmtCoef           = (*expvect)[ne]->GetNcoeffs();
                int nElmtCoef2          = nElmtCoef*nElmtCoef;
                int nElmtPnt            = (*expvect)[ne]->GetTotPoints();

                int nQuot = nElmtCoef2/m_nPadding;
                int nRemd = nElmtCoef2- nQuot*m_nPadding;
                int nQuotPlus=nQuot;
                if(nRemd>0)
                {
                    nQuotPlus++;
                }
                int nPaded = nQuotPlus*m_nPadding;

                ArrayStdMat[0] = MatStdDerivBase0;
                if(m_spacedim>1)
                {
                    StdRegions::StdMatrixKey  matkey(StdRegions::eDerivBase1,
                                    stdExp->DetShapeType(), *stdExp);
                    ArrayStdMat[1]  =   stdExp->GetStdMatrix(matkey);
                    
                    if(m_spacedim>2)
                    {
                        StdRegions::StdMatrixKey  matkey(
                            StdRegions::eDerivBase2, 
                            stdExp->DetShapeType(), *stdExp);
                        ArrayStdMat[2]  =   stdExp->GetStdMatrix(matkey);
                    }
                }
                for(int nd0=0;nd0<m_spacedim;nd0++)
                {
                    ArrayStdMatData[nd0] =  ArrayStdMat[nd0]->GetPtr();
                }

                StdRegions::StdMatrixKey  matkey(StdRegions::eBwdMat,
                                        stdExp->DetShapeType(), *stdExp);
                DNekMatSharedPtr BwdMat =  stdExp->GetStdMatrix(matkey);
                Array<OneD, NekDouble> BwdMatData = BwdMat->GetPtr();

                TensorOfArray3D<DataType> tmpStdDBB (m_spacedim);
                TensorOfArray4D<DataType> tmpStdDBDB(m_spacedim);

                for(int nd0=0;nd0<m_spacedim;nd0++)
                {
                    tmpStdDBB[nd0]  = TensorOfArray2D<DataType> (nElmtPnt);
                    for(int i=0;i<nElmtPnt;i++)
                    {
                        tmpStdDBB[nd0][i] = Array<OneD, DataType> (nPaded,0.0);
                        for(int nc1=0;nc1<nElmtCoef;nc1++)
                        {
                            int noffset = nc1*nElmtCoef;
                            for(int nc0=0;nc0<nElmtCoef;nc0++)
                            {
                                tmpStdDBB[nd0][i][nc0+noffset] = 
                                    DataType (ArrayStdMatData[nd0]
                                        [i*nElmtCoef+nc0]*
                                        BwdMatData[i*nElmtCoef+nc1]);
                            }
                        }
                    }

                    tmpStdDBDB[nd0] = TensorOfArray3D<DataType> (m_spacedim);
                    for(int nd1=0;nd1<m_spacedim;nd1++)
                    {
                        tmpStdDBDB[nd0][nd1] = 
                            TensorOfArray2D<DataType> (nElmtPnt);
                        for(int i=0;i<nElmtPnt;i++)
                        {
                            tmpStdDBDB[nd0][nd1][i] = 
                                Array<OneD, DataType> (nPaded,0.0);
                            for(int nc1=0;nc1<nElmtCoef;nc1++)
                            {
                                int noffset = nc1*nElmtCoef;
                                for(int nc0=0;nc0<nElmtCoef;nc0++)
                                {
                                    tmpStdDBDB[nd0][nd1][i][nc0+noffset] = 
                                        DataType(ArrayStdMatData[nd0]
                                            [i*nElmtCoef+nc0]*
                                            ArrayStdMatData[nd1]
                                            [i*nElmtCoef+nc1]);
                                }
                            }
                        }
                    }
                }
                VectStdDerivBase0.push_back(MatStdDerivBase0);
                VectStdDerivBase_Base.push_back(tmpStdDBB);
                VectStdDervBase_DervBase.push_back(tmpStdDBDB);

                StdMatDataDBB[ne]  = tmpStdDBB;
                StdMatDataDBDB[ne] = tmpStdDBDB;
            }
        }
    }


    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::AddMatNSBlkDiag_boundary(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        TensorOfArray3D<NekDouble>                       &qfield,
        TensorOfArray2D<TypeNekBlkMatSharedPtr>          &gmtxarray,
        Array<OneD, TypeNekBlkMatSharedPtr >             &TraceJac,
        Array<OneD, TypeNekBlkMatSharedPtr >             &TraceJacDeriv,
        TensorOfArray2D<DataType>                        &TraceJacDerivSign,
        TensorOfArray5D<DataType>                        &TraceIPSymJacArray)
    {
        int nvariables = inarray.size();
        GetTraceJac(inarray,qfield,TraceJac,TraceJacDeriv,TraceJacDerivSign,
            TraceIPSymJacArray);
        
        Array<OneD, TypeNekBlkMatSharedPtr > tmpJac;
        TensorOfArray2D<DataType>  tmpSign;

        m_advObject->AddTraceJacToMat(nvariables,m_spacedim,m_fields, 
            TraceJac,gmtxarray,tmpJac,tmpSign);
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::GetTraceJac(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        TensorOfArray3D<NekDouble>                       &qfield,
        Array<OneD, TypeNekBlkMatSharedPtr >             &TraceJac,
        Array<OneD, TypeNekBlkMatSharedPtr >             &TraceJacDeriv,
        TensorOfArray2D<DataType>                        &TraceJacDerivSign,
        TensorOfArray5D<DataType>                        &TraceIPSymJacArray)
    {
        boost::ignore_unused(TraceJacDeriv, TraceJacDerivSign);

        int nvariables = inarray.size();
        int nTracePts  = GetTraceTotPoints();

        // Store forwards/backwards space along trace space
        TensorOfArray2D<NekDouble> Fwd (nvariables);
        TensorOfArray2D<NekDouble> Bwd (nvariables);

        TypeNekBlkMatSharedPtr FJac,BJac;
        Array<OneD, unsigned int> n_blks1(nTracePts, nvariables);

        if(TraceJac.size()>0)
        {
            FJac = TraceJac[0];
            BJac = TraceJac[1];
        }
        else
        {
            TraceJac = Array<OneD, TypeNekBlkMatSharedPtr>(2);

            AllocateNekBlkMatDig(FJac,n_blks1, n_blks1);
            AllocateNekBlkMatDig(BJac,n_blks1, n_blks1);
        }
        
        if (m_HomogeneousType == eHomogeneous1D)
        {
            Fwd = NullNekDoubleArrayofArray;
            Bwd = NullNekDoubleArrayofArray;
        }
        else
        {
            for(int i = 0; i < nvariables; ++i)
            {
                Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
                Bwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }

        TensorOfArray2D<NekDouble> AdvVel(m_spacedim);

        NumCalRiemFluxJac(nvariables, m_fields, AdvVel, inarray,qfield,
            m_BndEvaluateTime, Fwd, Bwd, FJac, BJac,TraceIPSymJacArray);

        TraceJac[0] = FJac;
        TraceJac[1] = BJac;
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::CalcVisFluxDerivJac(
        const int                           nConvectiveFields,
        const TensorOfArray2D<NekDouble>    &inarray,
        const TensorOfArray2D<NekDouble>    &Fwd,
        const TensorOfArray2D<NekDouble>    &Bwd,
        TypeNekBlkMatSharedPtr              &BJac,
        DataType                            &tmpDataType)
    {
        boost::ignore_unused(inarray, tmpDataType);
        MultiRegions::ExpListSharedPtr tracelist = 
            m_fields[0]->GetTrace();
        std::shared_ptr<LocalRegions::ExpansionVector> traceExp= 
            tracelist->GetExp();
        int ntotTrac            = (*traceExp).size();
        int noffset,pntoffset;

        int nTracePts  = tracelist->GetTotPoints();

        Array<OneD, NekDouble> tracBwdWeightAver(nTracePts,0.0);
        Array<OneD, NekDouble> tracBwdWeightJump(nTracePts,0.0);
        m_fields[0]->GetBwdWeight(tracBwdWeightAver,tracBwdWeightJump);
        tracBwdWeightAver = NullNekDouble1DArray;
        Vmath::Ssub(nTracePts,2.0,tracBwdWeightJump,1,tracBwdWeightJump,1);
        Vmath::Smul(nTracePts,0.5,tracBwdWeightJump,1,tracBwdWeightJump,1);
        
        TensorOfArray2D<NekDouble>    solution_jump(nConvectiveFields);
        TensorOfArray2D<NekDouble>    solution_Aver(nConvectiveFields);
        for (int i = 0; i < nConvectiveFields; ++i)
        {
            solution_jump[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
            solution_Aver[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
        }

        m_diffusion->ConsVarAveJump(nConvectiveFields,nTracePts,Fwd,Bwd,
            solution_Aver,solution_jump);

        const Array<OneD, const Array<OneD, NekDouble> > tracenormals =   
            m_diffusion->GetTraceNormal();

        TensorOfArray2D<DNekMatSharedPtr> ElmtJac;
        ElmtJac = TensorOfArray2D<DNekMatSharedPtr> (ntotTrac);
        for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
        {
            int nTracPnt = (*traceExp)[nelmt]->GetTotPoints();
            ElmtJac[nelmt] = Array<OneD, DNekMatSharedPtr>(nTracPnt);
            for(int npnt = 0; npnt < nTracPnt; npnt++)
            {
                ElmtJac[nelmt][npnt] = MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nConvectiveFields, nConvectiveFields);
            }
        }
        
        DNekMatSharedPtr tmpMatinn;
        Array<OneD, NekDouble > tmpMatinnData;
        Array<OneD, DataType > tmpMatoutData;
        for(int nDervDir = 0; nDervDir < m_spacedim; nDervDir++)
        {
            GetFluxDerivJacDirctn(tracelist,tracenormals,nDervDir,
                solution_Aver,ElmtJac);

            for(int  ntrace = 0; ntrace < ntotTrac; ntrace++)
            {
                int nTracPnt    = (*traceExp)[ntrace]->GetTotPoints();
                noffset         = tracelist->GetPhys_Offset(ntrace);
                for(int npnt = 0; npnt < nTracPnt; npnt++)
                {
                    pntoffset = noffset+npnt;
                    NekDouble weight = tracBwdWeightJump[pntoffset];

                    tmpMatinn       = ElmtJac[ntrace][npnt];
                    tmpMatinnData   = tmpMatinn->GetPtr(); 
                    tmpMatoutData = 
                        BJac->GetBlock(pntoffset,pntoffset)->GetPtr(); 
                    Vmath::Smul(nConvectiveFields*nConvectiveFields,weight,
                        tmpMatinnData,1,tmpMatinnData,1);
                    for(int j = 0; j< nConvectiveFields;j++)
                    {
                        for(int i=0;i<nConvectiveFields;i++)
                        {
                            tmpMatoutData[(j*m_spacedim+nDervDir)*
                                nConvectiveFields+i] =
                                DataType( tmpMatinnData[j*nConvectiveFields+i] );
                        }
                    }

                }
            }
        }
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::NumCalRiemFluxJac(
        const int                                         nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const TensorOfArray2D<NekDouble>                  &AdvVel,
        const TensorOfArray2D<NekDouble>                  &inarray,
        TensorOfArray3D<NekDouble>                        &qfield,
        const NekDouble                                   &time,
        const TensorOfArray2D<NekDouble>                  &Fwd,
        const TensorOfArray2D<NekDouble>                  &Bwd,
        TypeNekBlkMatSharedPtr                            &FJac,
        TypeNekBlkMatSharedPtr                            &BJac,
        TensorOfArray5D<DataType>                         &TraceIPSymJacArray)
    {
        boost::ignore_unused(TraceIPSymJacArray);

        const NekDouble     PenaltyFactor2  =   0.0;
        int nvariables  = nConvectiveFields;
        int npoints     = GetNpoints();
        // int nPts        = npoints;
        int nTracePts   = GetTraceTotPoints();
        int nDim        = m_spacedim;

        Array<OneD, int > nonZeroIndex;

        TensorOfArray2D<NekDouble>  tmpinarry(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            tmpinarry[i]    =    Array<OneD, NekDouble>(npoints,0.0);
            Vmath::Vcopy(npoints, inarray[i],1,tmpinarry[i],1);
        }

        // DmuDT of artificial diffusion is neglected
        // TODO: to consider the Jacobian of AV seperately
        Array<OneD, NekDouble> muvar        =   NullNekDouble1DArray;
        Array<OneD, NekDouble> MuVarTrace   =   NullNekDouble1DArray;
        if (m_shockCaptureType != "Off" && m_shockCaptureType != "Physical")
        {
            MuVarTrace  =   Array<OneD, NekDouble>(nTracePts, 0.0);
            muvar       =   Array<OneD, NekDouble>(npoints, 0.0);
            m_diffusion->GetAVmu(fields,inarray,muvar,MuVarTrace);
            muvar       =   NullNekDouble1DArray;
        }

        TensorOfArray2D<NekDouble> numflux(nvariables);
        for(int i = 0; i < nvariables; ++i)
        {
            numflux[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        }

        const MultiRegions::AssemblyMapDGSharedPtr  TraceMap = 
            fields[0]->GetTraceMap();
        TensorOfArray3D<NekDouble>    qBwd(nDim);
        TensorOfArray3D<NekDouble>    qFwd(nDim);
        if(m_ViscousJacFlag)
        {
            for (int nd = 0; nd < nDim; ++nd)
            {
                qBwd[nd] = TensorOfArray2D<NekDouble> (nConvectiveFields);
                qFwd[nd] = TensorOfArray2D<NekDouble> (nConvectiveFields);
                for (int i = 0; i < nConvectiveFields; ++i)
                {
                    qBwd[nd][i] = Array<OneD, NekDouble>(nTracePts,0.0);
                    qFwd[nd][i] = Array<OneD, NekDouble>(nTracePts,0.0);

                    fields[i]->GetFwdBwdTracePhys(qfield[nd][i], qFwd[nd][i], 
                        qBwd[nd][i], true, true, false);
                    TraceMap->GetAssemblyCommDG()->PerformExchange(qFwd[nd][i], 
                        qBwd[nd][i]);
                }
            }
        }

        CalcTraceNumericalFlux(nConvectiveFields,nDim,npoints,nTracePts,
            PenaltyFactor2, fields,AdvVel,inarray,time,qfield,Fwd,Bwd,
            qFwd,qBwd,MuVarTrace,nonZeroIndex,numflux);

        int nFields = nvariables;
        TensorOfArray2D<NekDouble>  plusFwd(nFields),plusBwd(nFields);
        TensorOfArray2D<NekDouble>  Jacvect(nFields);
        TensorOfArray2D<NekDouble>  FwdBnd(nFields);
        TensorOfArray2D<NekDouble>  plusflux(nFields);
        for(int i = 0; i < nFields; i++)
        {
            Jacvect[i]  = Array<OneD, NekDouble>(nTracePts,0.0);
            plusFwd[i]  = Array<OneD, NekDouble>(nTracePts,0.0);
            plusBwd[i]  = Array<OneD, NekDouble>(nTracePts,0.0);
            plusflux[i] = Array<OneD, NekDouble>(nTracePts,0.0);
            FwdBnd[i]   = Array<OneD, NekDouble>(nTracePts,0.0);
        }


        for(int i = 0; i < nFields; i++)
        {
            Vmath::Vcopy(nTracePts, Fwd[i],1,plusFwd[i],1);
            Vmath::Vcopy(nTracePts, Bwd[i],1,plusBwd[i],1);
        }

        NekDouble eps   =   1.0E-6;
        
        Array<OneD, DataType> tmpMatData;
        // Fwd Jacobian
        for(int i = 0; i < nFields; i++)
        {
            NekDouble epsvar = eps*m_magnitdEstimat[i];
            NekDouble oepsvar   =   1.0/epsvar;
            Vmath::Sadd(nTracePts,epsvar,Fwd[i],1,plusFwd[i],1);

            if (m_bndConds.size())
            {
                for(int i = 0; i < nFields; i++)
                {
                    Vmath::Vcopy(nTracePts, plusFwd[i],1,FwdBnd[i],1);
                }
                // Loop over user-defined boundary conditions
                for (auto &x : m_bndConds)
                {
                    x->Apply(FwdBnd, tmpinarry, time);
                }
            }

            for(int j = 0; j < nFields; j++)
            {
                m_fields[j]->FillBwdWithBoundCond(plusFwd[j], plusBwd[j]);
            }

            CalcTraceNumericalFlux(nConvectiveFields,nDim,npoints,nTracePts,
                PenaltyFactor2,fields,AdvVel,inarray,time,qfield,
                plusFwd,plusBwd,qFwd,qBwd,MuVarTrace,nonZeroIndex,plusflux);

            for (int n = 0; n < nFields; n++)
            {
                Vmath::Vsub(nTracePts,plusflux[n],1,numflux[n],1,Jacvect[n],1);
                Vmath::Smul(nTracePts, oepsvar ,Jacvect[n],1,Jacvect[n],1);
            }
            for(int j = 0; j < nTracePts; j++)
            {
                tmpMatData  =   FJac->GetBlock(j,j)->GetPtr();
                for (int n = 0; n < nFields; n++)
                {
                    tmpMatData[n+i*nFields] = DataType(Jacvect[n][j]);
                }
            }

            Vmath::Vcopy(nTracePts, Fwd[i],1,plusFwd[i],1);
        }

        // Reset the boundary conditions
        if (m_bndConds.size())
        {
            for(int i = 0; i < nFields; i++)
            {
                Vmath::Vcopy(nTracePts, Fwd[i],1,FwdBnd[i],1);
            }
            // Loop over user-defined boundary conditions
            for (auto &x : m_bndConds)
            {
                x->Apply(FwdBnd, tmpinarry, time);
            }
        }

        for(int i = 0; i < nFields; i++)
        {
            Vmath::Vcopy(nTracePts, Bwd[i],1,plusBwd[i],1);
        }

        for(int i = 0; i < nFields; i++)
        {
            NekDouble epsvar    = eps*m_magnitdEstimat[i];
            NekDouble oepsvar   =   1.0/epsvar;

            Vmath::Sadd(nTracePts,epsvar,Bwd[i],1,plusBwd[i],1);

            for(int j = 0; j < nFields; j++)
            {
                m_fields[j]->FillBwdWithBoundCond(Fwd[j], plusBwd[j]);
            }

            CalcTraceNumericalFlux(nConvectiveFields,nDim,npoints,nTracePts,
                PenaltyFactor2,fields,AdvVel,inarray,time,qfield,Fwd,
                plusBwd,qFwd,qBwd,MuVarTrace,nonZeroIndex,plusflux);

            for (int n = 0; n < nFields; n++)
            {
                Vmath::Vsub(nTracePts,plusflux[n],1,numflux[n],1,Jacvect[n],1);
                Vmath::Smul(nTracePts, oepsvar,Jacvect[n],1,Jacvect[n],1);
            }
            for(int j = 0; j < nTracePts; j++)
            {
                tmpMatData  =   BJac->GetBlock(j,j)->GetPtr();
                for (int n = 0; n < nFields; n++)
                {
                    tmpMatData[n+i*nFields] = DataType(Jacvect[n][j]);
                }
            }

            Vmath::Vcopy(nTracePts, Bwd[i],1,plusBwd[i],1);
        }
    }

    void CompressibleFlowSystem::CalcTraceNumericalFlux(
            const int                                         nConvectiveFields,
            const int                                         nDim,
            const int                                         nPts,
            const int                                         nTracePts,
            const NekDouble                                   PenaltyFactor2,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const TensorOfArray2D<NekDouble>                  &AdvVel,
            const TensorOfArray2D<NekDouble>                  &inarray,
            const NekDouble                                   time,
            TensorOfArray3D<NekDouble>                        &qfield,
            const TensorOfArray2D<NekDouble>                  &vFwd,
            const TensorOfArray2D<NekDouble>                  &vBwd,
            const Array<OneD, const TensorOfArray2D<NekDouble>> &qFwd,
            const Array<OneD, const TensorOfArray2D<NekDouble>> &qBwd,
            const Array<OneD, NekDouble >                       &MuVarTrace,
            Array<OneD, int >                                   &nonZeroIndex,
            TensorOfArray2D<NekDouble>                          &traceflux)
    {
        boost::ignore_unused(nDim, nPts, PenaltyFactor2, time, qFwd, qBwd, 
            MuVarTrace);

        if (m_AdvectionJacFlag)
        {
            m_advObject->AdvectTraceFlux(nConvectiveFields, m_fields, AdvVel,
                inarray, traceflux, m_BndEvaluateTime,vFwd, vBwd);
        }
        else
        {
            for (int i = 0; i < nConvectiveFields; i++)
            {
                traceflux[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
            }
        }
        
        if(m_ViscousJacFlag)
        {
            Array<OneD, Array<OneD, NekDouble > > visflux(nConvectiveFields);
            for(int i = 0; i < nConvectiveFields; i++)
            {
                visflux[i]  =    Array<OneD, NekDouble>(nTracePts,0.0);
            }

            m_diffusion->DiffuseTraceFlux(fields, inarray, qfield, 
                NullTensorOfArray3DDouble, visflux, vFwd, vBwd, nonZeroIndex);
            for(int i = 0; i < nConvectiveFields; i++)
            {
                Vmath::Vsub(nTracePts,traceflux[i],1,visflux[i],1,
                    traceflux[i],1);
            }
        }
    }

    void CompressibleFlowSystem::CalcTraceIPSymmFlux(
            const int                                         nConvectiveFields,
            const int                                         nTracePts,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const TensorOfArray2D<NekDouble>                  &inarray,
            const NekDouble                                   time,
            const Array<OneD, const TensorOfArray2D<NekDouble>> &qfield,
            const TensorOfArray2D<NekDouble>                    &vFwd,
            const TensorOfArray2D<NekDouble>                    &vBwd,
            const Array<OneD, NekDouble >                       &MuVarTrace,
            Array<OneD, int >                                   &nonZeroIndex,
            TensorOfArray3D<NekDouble>                          &traceflux)
    {
        boost::ignore_unused(time, MuVarTrace);

        if(m_ViscousJacFlag)
        {
            TensorOfArray3D<NekDouble>  VolumeFlux;
            m_diffusion->DiffuseTraceSymmFlux(nConvectiveFields,fields,
                inarray,qfield,VolumeFlux,traceflux,vFwd,vBwd,nonZeroIndex);
            for(int nd = 0; nd < m_spacedim; nd++)
            {
                for(int i = 0; i < nConvectiveFields; i++)
                {
                    Vmath::Ssub(nTracePts,0.0,traceflux[nd][i],1,
                        traceflux[nd][i],1);
                }
            }
            
        }
    }

    template<typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::AllocatePrecBlkDiagCoeff(
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        const int                                          &nscale )
    {

        int nvars = m_fields.size();
        int nelmts  = m_fields[0]->GetNumElmts();
        int nelmtcoef;
        Array<OneD, unsigned int > nelmtmatdim(nelmts);
        for(int i = 0; i < nelmts; i++)
        {
            nelmtcoef   =   m_fields[0]->GetExp(i)->GetNcoeffs();
            nelmtmatdim[i]  =   nelmtcoef*nscale;
        }

        for(int i = 0; i < nvars; i++)
        {
            for(int j = 0; j < nvars; j++)
            {
                AllocateNekBlkMatDig(gmtxarray[i][j],nelmtmatdim,nelmtmatdim);
            }
        }
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::GetPrecNSBlkDiagCoeff(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        TensorOfArray2D<TypeNekBlkMatSharedPtr>          &gmtxarray,
        TypeNekBlkMatSharedPtr                           &gmtVar,
        Array<OneD, TypeNekBlkMatSharedPtr >             &TraceJac,
        Array<OneD, TypeNekBlkMatSharedPtr >             &TraceJacDeriv,
        TensorOfArray2D<DataType>                        &TraceJacDerivSign,
        TensorOfArray4D<DataType>                        &TraceJacArray,
        TensorOfArray4D<DataType>                        &TraceJacDerivArray,
        TensorOfArray5D<DataType>                        &TraceIPSymJacArray,
        TensorOfArray4D<DataType>                        &StdMatDataDBB,
        TensorOfArray5D<DataType>                        &StdMatDataDBDB)
    {
        TensorOfArray3D<NekDouble> qfield;

        if(m_ViscousJacFlag)
        {
            CalcPhysDeriv(inarray,qfield);
        }

        DataType zero =0.0;
        Fill2DArrayOfBlkDiagonalMat(gmtxarray,zero);

        AddMatNSBlkDiag_volume(inarray,qfield,gmtxarray,StdMatDataDBB,
            StdMatDataDBDB);

        AddMatNSBlkDiag_boundary(inarray,qfield,gmtxarray,TraceJac,
            TraceJacDeriv,TraceJacDerivSign,TraceIPSymJacArray);

        MultiplyElmtInvMass_PlusSource(gmtxarray,m_TimeIntegLambda,zero);

        ElmtVarInvMtrx(gmtxarray,gmtVar,zero);

        TransTraceJacMatToArray(TraceJac,TraceJacDeriv,TraceJacArray, 
            TraceJacDerivArray);
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::MultiplyElmtInvMass_PlusSource(
        TensorOfArray2D<TypeNekBlkMatSharedPtr> &gmtxarray,
        const NekDouble                         dtlamda,
        const DataType                          tmpDataType)
    {
        boost::ignore_unused(tmpDataType);

        MultiRegions::ExpListSharedPtr explist = m_fields[0];
        std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
        int ntotElmt            = (*pexp).size();
        int nConvectiveFields = m_fields.size();

        NekDouble Negdtlamda    =   -dtlamda;

        Array<OneD, NekDouble> pseudotimefactor(ntotElmt,0.0);
        if(m_cflLocTimestep>0.0)
        {
            NekDouble Neglamda      =   Negdtlamda/m_timestep;
            Vmath::Svtvp(ntotElmt,Neglamda,m_locTimeStep,1,
                pseudotimefactor,1,pseudotimefactor,1);
        }
        else
        {
            Vmath::Fill(ntotElmt,Negdtlamda,pseudotimefactor,1);
        }

        Array<OneD, DataType>    GMatData;
        for(int m = 0; m < nConvectiveFields; m++)
        {
            for(int n = 0; n < nConvectiveFields; n++)
            {
                for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                {
                    GMatData = gmtxarray[m][n]->GetBlock(nelmt,nelmt)->GetPtr();
                    DataType factor = DataType(pseudotimefactor[nelmt]);

                    Vmath::Smul(GMatData.size(),factor,GMatData,1,GMatData,1);
                }
            }
        }

        DNekMatSharedPtr MassMat;
        Array<OneD,NekDouble> BwdMatData,MassMatData,tmp;
        Array<OneD,NekDouble> tmp2;
        Array<OneD,DataType> MassMatDataDataType;
        LibUtilities::ShapeType ElmtTypePrevious = LibUtilities::eNoShapeType;

        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            int nelmtcoef  = GetNcoeffs(nelmt);
            int nelmtpnts  = GetTotPoints(nelmt);
            LibUtilities::ShapeType ElmtTypeNow =   
                explist->GetExp(nelmt)->DetShapeType();

            if (tmp.size()!=nelmtcoef||(ElmtTypeNow!=ElmtTypePrevious)) 
            {
                StdRegions::StdExpansionSharedPtr stdExp;
                stdExp = explist->GetExp(nelmt)->GetStdExp();
                StdRegions::StdMatrixKey  matkey(StdRegions::eBwdTrans,
                                    stdExp->DetShapeType(), *stdExp);
                    
                DNekMatSharedPtr BwdMat =  stdExp->GetStdMatrix(matkey);
                BwdMatData = BwdMat->GetPtr();

                if(nelmtcoef!=tmp.size())
                {
                    tmp = Array<OneD,NekDouble> (nelmtcoef,0.0);
                    MassMat  =  MemoryManager<DNekMat>
                        ::AllocateSharedPtr(nelmtcoef, nelmtcoef, 0.0);
                    MassMatData = MassMat->GetPtr();
                    MassMatDataDataType = 
                        Array<OneD, DataType> (nelmtcoef*nelmtcoef);
                }

                ElmtTypePrevious    = ElmtTypeNow;
            }
            
            for(int np=0; np<nelmtcoef;np++)
            {
                explist->GetExp(nelmt)->IProductWRTBase(
                    BwdMatData+np*nelmtpnts,tmp);
                Vmath::Vcopy(nelmtcoef,tmp,1,tmp2 = MassMatData+np*nelmtcoef,1);
            }
            for(int i=0;i<MassMatData.size();i++)
            {
                MassMatDataDataType[i]    =   DataType(MassMatData[i]);
            }

            for(int m = 0; m < nConvectiveFields; m++)
            {
                GMatData = gmtxarray[m][m]->GetBlock(nelmt,nelmt)->GetPtr();
                Vmath::Vadd(MassMatData.size(),MassMatDataDataType,1,
                    GMatData,1,GMatData,1);
            }
        }
        return;
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::TransTraceJacMatToArray(
        const Array<OneD, TypeNekBlkMatSharedPtr > &TraceJac,
        const Array<OneD, TypeNekBlkMatSharedPtr > &TraceJacDeriv,
        TensorOfArray4D<DataType>                  &TraceJacArray,
        TensorOfArray4D<DataType>                  &TraceJacDerivArray)
    {
        boost::ignore_unused(TraceJacArray, TraceJacDeriv, TraceJacDerivArray);

        int nFwdBwd,nDiagBlks,nvar0Jac,nvar1Jac;

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;
        nFwdBwd = TraceJac.size();
        TraceJac[0]->GetBlockSizes(rowSizes,colSizes);
        nDiagBlks   = rowSizes.size();
        nvar0Jac    = rowSizes[1] - rowSizes[0];
        nvar1Jac    = colSizes[1] - colSizes[0];

        if(0==TraceJacArray.size())
        {
            TraceJacArray = TensorOfArray4D<DataType> (nFwdBwd);
            for(int nlr=0;nlr<nFwdBwd;nlr++)
            {
                TraceJacArray[nlr] = TensorOfArray3D<DataType> (nvar0Jac);
                for(int m=0;m<nvar0Jac;m++)
                {
                    TraceJacArray[nlr][m] = TensorOfArray2D<DataType>(nvar1Jac);
                    for(int n=0;n<nvar1Jac;n++)
                    {
                        TraceJacArray[nlr][m][n] = 
                            Array<OneD,DataType >(nDiagBlks);
                    }
                }
            }
        }

        for(int nlr=0;nlr<nFwdBwd;nlr++)
        {
            const TypeNekBlkMatSharedPtr tmpMat = TraceJac[nlr];
            TensorOfArray3D<DataType> tmpaa = TraceJacArray[nlr];
            TranSamesizeBlkDiagMatIntoArray(tmpMat,tmpaa);
        }

        return;
    }

    void CompressibleFlowSystem::MatrixMultiplyMatrixFreeCoeffCentral(
            const  Array<OneD, NekDouble> &inarray,
                Array<OneD, NekDouble >&out)
    {
        NekDouble eps = m_JacobiFreeEps;
        unsigned int ntotalGlobal     = inarray.size();
        NekDouble magninarray = Vmath::Dot(ntotalGlobal,inarray,inarray);
        m_comm->AllReduce(magninarray, Nektar::LibUtilities::ReduceSum);
        // eps *= magnitdEstimatMax;
        eps *= sqrt( (sqrt(m_inArrayNorm) + 1.0)/magninarray);
        NekDouble oeps = 1.0/eps;
        unsigned int nvariables = m_fields.size();
        unsigned int npoints    = GetNcoeffs();
        Array<OneD, NekDouble > tmp;
        TensorOfArray2D<NekDouble> solplus(nvariables);
        TensorOfArray2D<NekDouble> resplus(nvariables);
        TensorOfArray2D<NekDouble> resminus(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            solplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resminus[i] =  Array<OneD, NekDouble>(npoints,0.0);
        }

        for (int i = 0; i < nvariables; i++)
        {
            tmp = inarray + i*npoints;
        }
        NonlinSysEvaluatorCoeff(solplus,resplus);

        for (int i = 0; i < nvariables; i++)
        {
            tmp = inarray + i*npoints;
        }
        NonlinSysEvaluatorCoeff(solplus,resminus);

        for (int i = 0; i < nvariables; i++)
        {
            tmp = out + i*npoints;
            Vmath::Vsub(npoints,resplus[i],1,resminus[i],1,tmp,1);
            Vmath::Smul(npoints, 0.5*oeps ,tmp,1,tmp,1);
        }
       
        return;
    }

    void CompressibleFlowSystem::MatrixMultiplyMatrixFreeCoeffDualtimestep(
        const  Array<OneD, NekDouble> &inarray,
            Array<OneD, NekDouble >&out,
        const  bool                   &controlFlag)
    {
        if(controlFlag)
        {
            MatrixMultiplyMatrixFreeCoeffCentral(inarray,out);
        }
        else
        {
            MatrixMultiplyMatrixFreeCoeff(inarray,out);
        }

        if(m_cflLocTimestep>0.0)
        {
            int nElements = m_fields[0]->GetExpSize();
            int nelmtcoeffs,offset,varoffset;
            NekDouble fac;
            Array<OneD, NekDouble> tmp;
            Array<OneD, NekDouble> tmpout;
            unsigned int nvariables = m_fields.size();
            unsigned int ntotcoeffs = GetNcoeffs();

            Array<OneD, NekDouble> pseudotimefactor(nElements,0.0);
            NekDouble   otimestep   =   1.0/m_timestep;
            Vmath::Smul(nElements,otimestep,m_locTimeStep,1,pseudotimefactor,1);

            Vmath::Vsub(inarray.size(),out,1,inarray,1,out,1);
            for(int i = 0; i < nvariables; ++i)
            {
                varoffset = i*ntotcoeffs;
                tmpout = out + varoffset;
                // Loop over elements
                for(int n = 0; n < nElements; ++n)
                {
                    nelmtcoeffs = m_fields[0]->GetExp(n)->GetNcoeffs();
                    offset      = m_fields[0]->GetCoeff_Offset(n);
                    fac         = pseudotimefactor[n];
                    tmp         = tmpout + offset;
                    Vmath::Smul(nelmtcoeffs,fac,tmp,1,tmp,1);
                }
            }
            Vmath::Vadd(inarray.size(),out,1,inarray,1,out,1);
        }
        return;
    }

    /**
     * @brief Compute the right-hand side.
     */
    void CompressibleFlowSystem::DoOdeRhsCoeff(
        const TensorOfArray2D<NekDouble>    &inarray,
        TensorOfArray2D<NekDouble>          &outarray,
        const NekDouble                     time)
    {
        int nvariables = inarray.size();
        int nTracePts  = GetTraceTotPoints();
        int ncoeffs    = GetNcoeffs();
        // Store forwards/backwards space along trace space
        TensorOfArray2D<NekDouble> Fwd    (nvariables);
        TensorOfArray2D<NekDouble> Bwd    (nvariables);

        if (m_HomogeneousType == eHomogeneous1D)
        {
            Fwd = NullNekDoubleArrayofArray;
            Bwd = NullNekDoubleArrayofArray;
        }
        else
        {
            for (int i = 0; i < nvariables; ++i)
            {
                Fwd[i]     = TensorOfArray1D<NekDouble> (nTracePts, 0.0);
                Bwd[i]     = TensorOfArray1D<NekDouble> (nTracePts, 0.0);
                m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }
 
         // Calculate advection
        DoAdvectionCoeff(inarray, outarray, time, Fwd, Bwd);

        // Negate results
        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(ncoeffs, outarray[i], 1);
        }

        DoDiffusionCoeff(inarray, outarray, Fwd, Bwd);

        // Add forcing terms
        for (auto &x : m_forcing)
        {
            x->ApplyCoeff(m_fields, inarray, outarray, time);
        }

        if (m_useLocalTimeStep)
        {
            int nElements = m_fields[0]->GetExpSize();
            int nq, offset;
            NekDouble fac;
            TensorOfArray1D<NekDouble> tmp;

            TensorOfArray1D<NekDouble> tstep(nElements, 0.0);
            GetElmtTimeStep(inarray, tstep);

            // Loop over elements
            for (int n = 0; n < nElements; ++n)
            {
                nq     = m_fields[0]->GetExp(n)->GetNcoeffs();
                offset = m_fields[0]->GetCoeff_Offset(n);
                fac    = tstep[n] / m_timestep;
                for (int i = 0; i < nvariables; ++i)
                {
                    Vmath::Smul(nq, fac, outarray[i] + offset, 1,
                                tmp = outarray[i] + offset, 1);
                }
            }
        }
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::TranSamesizeBlkDiagMatIntoArray(
        const TypeNekBlkMatSharedPtr                        &BlkMat,
        TensorOfArray3D<DataType>       &MatArray)
    {
        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;
        BlkMat->GetBlockSizes(rowSizes,colSizes);
        int nDiagBlks   = rowSizes.size();
        int nvar0       = rowSizes[1] - rowSizes[0];
        int nvar1       = colSizes[1] - colSizes[0];

        Array<OneD, DataType>    ElmtMatData;

        for(int i=0;i<nDiagBlks;i++)
        {
            ElmtMatData = BlkMat->GetBlock(i,i)->GetPtr();
            for(int n=0;n<nvar1;n++)
            {
                int noffset = n*nvar0;
                for(int m=0;m<nvar0;m++)
                {
                    MatArray[m][n][i]   =   ElmtMatData[m+noffset];
                }
            }
        }
    }
    
    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::Fill2DArrayOfBlkDiagonalMat(
        TensorOfArray2D<TypeNekBlkMatSharedPtr>   &gmtxarray,
        const DataType                                      valu)
    {

        int n1d = gmtxarray.size();

        for(int n1 = 0; n1 < n1d; ++n1)
        {
            Fill1DArrayOfBlkDiagonalMat(gmtxarray[n1],valu);
        }
    }

    
    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CompressibleFlowSystem::Fill1DArrayOfBlkDiagonalMat( 
        Array<OneD, TypeNekBlkMatSharedPtr >    &gmtxarray,
        const DataType                          valu)
    {
        int n1d = gmtxarray.size();

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;

        Array<OneD, DataType > loc_mat_arr;


        for(int n1 = 0; n1 < n1d; ++n1)
        {
            gmtxarray[n1]->GetBlockSizes(rowSizes,colSizes);
            int nelmts  = rowSizes.size();

            for(int i = 0; i < nelmts; ++i)
            {
                loc_mat_arr = gmtxarray[n1]->GetBlock(i,i)->GetPtr();

                int nrows = gmtxarray[n1]->GetBlock(i,i)->GetRows();
                int ncols = gmtxarray[n1]->GetBlock(i,i)->GetColumns();

                Vmath::Fill(nrows*ncols,valu,loc_mat_arr,1);
            }
        }

    }

    void CompressibleFlowSystem::GetFluxVectorJacDirctnElmt(
        const int                           nConvectiveFields,
        const int                           nElmtPnt,
        const TensorOfArray2D<NekDouble>    &locVars,
        const Array<OneD, NekDouble>        &normals,
        DNekMatSharedPtr                    &wspMat,
        TensorOfArray2D<NekDouble>          &PntJacArray)
    {
        Array<OneD, NekDouble> wspMatData = wspMat->GetPtr();

        int matsize = nConvectiveFields*nConvectiveFields;

        Array<OneD, NekDouble> pointVar(nConvectiveFields);

        for(int npnt = 0; npnt < nElmtPnt; npnt++)
        {
            for(int j = 0; j < nConvectiveFields; j++)
            {
                pointVar[j] = locVars[j][npnt];
            }

            GetFluxVectorJacPoint(nConvectiveFields,pointVar,normals,wspMat);

            Vmath::Vcopy(matsize, wspMatData,1,PntJacArray[npnt],1);
        }
        return ;
    }

    void CompressibleFlowSystem::GetFluxVectorJacPoint(
        const int                                   nConvectiveFields,
        const Array<OneD, NekDouble>                &conservVar, 
        const Array<OneD, NekDouble>                &normals, 
        DNekMatSharedPtr                            &fluxJac)
    {
        int nvariables      = conservVar.size();
        const int nvariables3D    = 5;
        int expDim          = m_spacedim;

        NekDouble fsw,efix_StegerWarming;
        efix_StegerWarming = 0.0;
        fsw = 0.0; // exact flux Jacobian if fsw=0.0
        if (nvariables > expDim+2)
        {
            NEKERROR(ErrorUtil::efatal,"nvariables > expDim+2 case not coded")
        }

        Array<OneD, NekDouble> fluxJacData;
        ;
        fluxJacData = fluxJac->GetPtr();

        if(nConvectiveFields==nvariables3D)
        {
            PointFluxJacobian_pn(conservVar,normals,fluxJac,
                efix_StegerWarming,fsw);
        }
        else
        {
            DNekMatSharedPtr PointFJac3D = MemoryManager<DNekMat>
                ::AllocateSharedPtr(nvariables3D, nvariables3D,0.0);
            
            Array<OneD, NekDouble> PointFJac3DData;
            PointFJac3DData = PointFJac3D->GetPtr();

            Array<OneD, NekDouble> PointFwd(nvariables3D,0.0);

            Array<OneD, unsigned int> index(nvariables);

            index[nvariables-1] = 4;
            for(int i=0;i<nvariables-1;i++)
            {
                index[i] = i;
            }

            int nj=0;
            int nk=0;
            for(int j=0; j< nvariables; j++)
            {
                nj = index[j];
                PointFwd[nj] = conservVar[j];
            }
            
            PointFluxJacobian_pn(PointFwd,normals,PointFJac3D,
                efix_StegerWarming,fsw);

            for(int j=0; j< nvariables; j++)
            {
                nj = index[j];
                for(int k=0; k< nvariables; k++)
                {
                    nk = index[k];
                    fluxJacData[j+k*nConvectiveFields] = 
                        PointFJac3DData[nj+nk*nvariables3D]; 
                }
            }
        }
    }

    // Currently duplacate in compressibleFlowSys
    // if fsw=+-1 calculate the steger-Warming flux vector splitting flux Jacobian
    // if fsw=0   calculate the Jacobian of the exact flux
    // efix is the numerical flux entropy fix parameter
    void CompressibleFlowSystem::PointFluxJacobian_pn(
            const Array<OneD, NekDouble> &Fwd,
            const Array<OneD, NekDouble> &normals,
                  DNekMatSharedPtr       &FJac,
            const NekDouble efix,   const NekDouble fsw)
    {
        Array<OneD, NekDouble> FJacData = FJac->GetPtr();
        const int nvariables3D    = 5;

        NekDouble ro,vx,vy,vz,ps,gama,ae ;
        NekDouble a,a2,h,h0,v2,vn,eps,eps2;
        NekDouble nx,ny,nz;
        NekDouble sn,osn,nxa,nya,nza,vna;
        NekDouble l1,l4,l5,al1,al4,al5,x1,x2,x3,y1;
        NekDouble c1,d1,c2,d2,c3,d3,c4,d4,c5,d5;
        NekDouble sml_ssf= 1.0E-12;

        NekDouble fExactorSplt = 2.0-abs(fsw); // if fsw=+-1 calculate 

        NekDouble   rhoL  = Fwd[0];
        NekDouble   rhouL = Fwd[1];
        NekDouble   rhovL = Fwd[2];
        NekDouble   rhowL = Fwd[3];
        NekDouble   EL    = Fwd[4];

        ro = rhoL;
        vx = rhouL / rhoL;
        vy = rhovL / rhoL;
        vz = rhowL / rhoL;

        // Internal energy (per unit mass)
        NekDouble eL =
                (EL - 0.5 * (rhouL * vx + rhovL * vy + rhowL * vz)) / rhoL;
        // TODO:
        // ps = m_eos->GetPressure(rhoL, eL);
        // gama = m_eos->GetGamma();
        ps      = m_varConv->Geteos()->GetPressure(rhoL, eL);
        gama    = m_gamma;

        ae = gama - 1.0;
        v2 = vx*vx + vy*vy + vz*vz;
        a2 = gama*ps/ro;
        h = a2/ae;

        h0 = h + 0.5*v2;
        a = sqrt(a2);

        nx = normals[0];
        ny = normals[1];
        nz = normals[2];
        vn = nx*vx + ny*vy + nz*vz;
        sn = std::max(sqrt(nx*nx + ny*ny + nz*nz),sml_ssf);
        osn = 1.0/sn;

        nxa = nx * osn;
        nya = ny * osn;
        nza = nz * osn;
        vna = vn * osn;
        l1 = vn;
        l4 = vn + sn*a;
        l5 = vn - sn*a;

        eps = efix*sn;
        eps2 = eps*eps;

        al1 = sqrt(l1*l1 + eps2);
        al4 = sqrt(l4*l4 + eps2);
        al5 = sqrt(l5*l5 + eps2);

        l1 = 0.5*(fExactorSplt*l1 + fsw*al1);
        l4 = 0.5*(fExactorSplt*l4 + fsw*al4);
        l5 = 0.5*(fExactorSplt*l5 + fsw*al5);

        x1 = 0.5*(l4 + l5);
        x2 = 0.5*(l4 - l5);
        x3 = x1 - l1;
        y1 = 0.5*v2;
        c1 = ae*x3/a2;
        d1 = x2/a;

        int nVar0 = 0;
        int nVar1 = nvariables3D;
        int nVar2 = 2*nvariables3D;
        int nVar3 = 3*nvariables3D;
        int nVar4 = 4*nvariables3D;
        FJacData[     nVar0] = c1*y1 - d1*vna + l1;
        FJacData[     nVar1] = -c1*vx + d1*nxa;
        FJacData[     nVar2] = -c1*vy + d1*nya;
        FJacData[     nVar3] = -c1*vz + d1*nza;
        FJacData[     nVar4] = c1;
        c2 = c1*vx + d1*nxa*ae;
        d2 = x3*nxa + d1*vx;
        FJacData[ 1 + nVar0] = c2*y1 - d2*vna;
        FJacData[ 1 + nVar1] = -c2*vx + d2*nxa + l1;
        FJacData[ 1 + nVar2] = -c2*vy + d2*nya;
        FJacData[ 1 + nVar3] = -c2*vz + d2*nza;
        FJacData[ 1 + nVar4] = c2;
        c3 = c1*vy + d1*nya*ae;
        d3 = x3*nya + d1*vy;
        FJacData[ 2 + nVar0] = c3*y1 - d3*vna;
        FJacData[ 2 + nVar1] = -c3*vx + d3*nxa;
        FJacData[ 2 + nVar2] = -c3*vy + d3*nya + l1;
        FJacData[ 2 + nVar3] = -c3*vz + d3*nza;
        FJacData[ 2 + nVar4] = c3;
        c4 = c1*vz + d1*nza*ae;
        d4 = x3*nza + d1*vz;
        FJacData[ 3 + nVar0] = c4*y1 - d4*vna;
        FJacData[ 3 + nVar1] = -c4*vx + d4*nxa;
        FJacData[ 3 + nVar2] = -c4*vy + d4*nya;
        FJacData[ 3 + nVar3] = -c4*vz + d4*nza + l1;
        FJacData[ 3 + nVar4] = c4;
        c5 = c1*h0 + d1*vna*ae;
        d5 = x3*vna + d1*h0;
        FJacData[ 4 + nVar0] = c5*y1 - d5*vna;
        FJacData[ 4 + nVar1] = -c5*vx + d5*nxa;
        FJacData[ 4 + nVar2] = -c5*vy + d5*nya;
        FJacData[ 4 + nVar3] = -c5*vz + d5*nza;
        FJacData[ 4 + nVar4] = c5 + l1;
    }

    /**
     * @brief Compute the advection terms for the right-hand side
     */
    void CompressibleFlowSystem::DoAdvection(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              TensorOfArray2D<NekDouble> &outarray,
        const NekDouble                                   time,
        const TensorOfArray2D<NekDouble>       &pFwd,
        const TensorOfArray2D<NekDouble>       &pBwd)
    {
        int nvariables = inarray.size();
        TensorOfArray2D<NekDouble> advVel(m_spacedim);

        m_advObject->Advect(nvariables, m_fields, advVel, inarray,
                            outarray, time, pFwd, pBwd);
    }

    /**
     * @brief Add the diffusions terms to the right-hand side
     */
    void CompressibleFlowSystem::DoDiffusion(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              TensorOfArray2D<NekDouble> &outarray,
            const TensorOfArray2D<NekDouble>   &pFwd,
            const TensorOfArray2D<NekDouble>   &pBwd)
    {
        v_DoDiffusion(inarray, outarray, pFwd, pBwd);

        if (m_shockCaptureType != "Off" && m_shockCaptureType != "Physical")
        {
            m_artificialDiffusion->DoArtificialDiffusion(inarray, outarray);
        }
    }

    void CompressibleFlowSystem::NonlinSysEvaluatorCoeff1D(
        const TensorOfArray1D<NekDouble>    &inarray,
        TensorOfArray1D<NekDouble>          &out,
        const bool                          &flag)
    {
        const TensorOfArray1D<NekDouble> refsource 
            = m_nonlinsol->GetRefSourceVec();
        NonlinSysEvaluatorCoeff(inarray, out, flag, refsource);
    }

    void CompressibleFlowSystem::NonlinSysEvaluatorCoeff(
        const TensorOfArray1D<NekDouble>    &inarray,
        TensorOfArray1D<NekDouble>          &out,
        const bool                          &flag,
        const TensorOfArray1D<NekDouble>    &source)
    {
        boost::ignore_unused(flag);
        unsigned int nvariables     = m_fields.size();
        unsigned int npoints        = m_fields[0]->GetNcoeffs();
        TensorOfArray2D<NekDouble> in2D(nvariables);
        TensorOfArray2D<NekDouble> out2D(nvariables);
        TensorOfArray2D<NekDouble> source2D(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            int offset = i * npoints;
            in2D[i]    = inarray + offset;
            out2D[i]   = out + offset;
            source2D[i]   = source + offset;
        }
        NonlinSysEvaluatorCoeff(in2D, out2D, source2D);
    }

    void CompressibleFlowSystem::NonlinSysEvaluatorCoeff(
        const TensorOfArray2D<NekDouble> &inarray,
        TensorOfArray2D<NekDouble>       &out,
        const TensorOfArray2D<NekDouble> &source)
    {
        unsigned int nvariable  = inarray.size();
        unsigned int ncoeffs    = inarray[nvariable - 1].size();
        unsigned int npoints    = m_fields[0]->GetNpoints();

        TensorOfArray2D<NekDouble> inpnts(nvariable);

        for (int i = 0; i < nvariable; ++i)
        {
            inpnts[i] = TensorOfArray1D<NekDouble>(npoints, 0.0);
            m_fields[i]->BwdTrans(inarray[i], inpnts[i]);
        }
        
        DoOdeProjection(inpnts, inpnts, m_BndEvaluateTime);
        DoOdeRhsCoeff(inpnts, out, m_BndEvaluateTime);

        for (int i = 0; i < nvariable; ++i)
        {
            Vmath::Svtvp(ncoeffs, -m_TimeIntegLambda, out[i], 1,
                         inarray[i], 1, out[i], 1);
        }

        if (NullTensorOfArray2DDouble != source)
        {
            for (int i = 0; i < nvariable; ++i)
            {
                Vmath::Vsub(ncoeffs, out[i], 1,
                            source[i], 1, out[i], 1);
            }
        }
        return;
    }

    /**
     * @brief Compute the advection terms for the right-hand side
     */
    void CompressibleFlowSystem::DoAdvectionCoeff(
        const TensorOfArray2D<NekDouble>    &inarray,
        TensorOfArray2D<NekDouble>          &outarray,
        const NekDouble                     time,
        const TensorOfArray2D<NekDouble>    &pFwd,
        const TensorOfArray2D<NekDouble>    &pBwd)
    {
        int nvariables = inarray.size();
        TensorOfArray2D<NekDouble> advVel(m_spacedim);

        m_advObject->AdvectCoeffs(nvariables, m_fields, advVel, inarray,
                                  outarray, time, pFwd, pBwd);
    }

    /**
     * @brief Add the diffusions terms to the right-hand side
     * Similar to DoDiffusion() but with outarray in coefficient space
     */
    void CompressibleFlowSystem::DoDiffusionCoeff(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              TensorOfArray2D<NekDouble> &outarray,
            const TensorOfArray2D<NekDouble>   &pFwd,
            const TensorOfArray2D<NekDouble>   &pBwd)
    {
        v_DoDiffusionCoeff(inarray, outarray, pFwd, pBwd);
    }

    void CompressibleFlowSystem::DoImplicitSolvePhysToCoeff(
        const TensorOfArray2D<NekDouble>    &inpnts,
        TensorOfArray2D<NekDouble>          &outpnt,
        const NekDouble                     time,
        const NekDouble                     lambda)
    {
        unsigned int nvariables  = inpnts.size();
        unsigned int ncoeffs     = m_fields[0]->GetNcoeffs();
        unsigned int ntotal      = nvariables * ncoeffs;

        TensorOfArray1D<NekDouble>  inarray(ntotal);
        TensorOfArray1D<NekDouble>  out(ntotal);
        TensorOfArray1D<NekDouble>  tmpArray;

        for (int i = 0; i < nvariables; ++i)
        {
            int noffset = i * ncoeffs;
            tmpArray = inarray + noffset;
            m_fields[i]->FwdTrans(inpnts[i], tmpArray);
        }

        DoImplicitSolveCoeff(inpnts, inarray, out, time, lambda);

        for (int i = 0; i < nvariables; ++i)
        {
            int noffset = i * ncoeffs;
            tmpArray = out + noffset;
            m_fields[i]->BwdTrans(tmpArray, outpnt[i]);
        }
    }

    void CompressibleFlowSystem::DoImplicitSolveCoeff(
            const TensorOfArray2D<NekDouble>    &inpnts,
            const TensorOfArray1D<NekDouble>    &inarray,
            TensorOfArray1D<NekDouble>          &out,
            const NekDouble                     time,
            const NekDouble                     lambda)
    {
        boost::ignore_unused(inpnts);

        m_TimeIntegLambda   = lambda;
        m_BndEvaluateTime   = time;
        m_solutionPhys      = inpnts;
        unsigned int ntotal = inarray.size();

        if (m_inArrayNorm < 0.0)
        {
            CalcRefValues(inarray);
        }
        
        NekDouble tol2 = m_inArrayNorm
                        *m_NewtonAbsoluteIteTol * m_NewtonAbsoluteIteTol;

        m_nonlinsol->v_SetupNekNonlinSystem(ntotal, inarray, inarray, 0);
        m_updatePreconMatFlag = UpdatePreconMatCheck(m_nonlinsol->GetRefResidual());

        m_TotNewtonIts +=  m_nonlinsol->SolveSystem(ntotal,inarray,
                                                    out, 0, tol2);

        m_TotLinIts += m_nonlinsol->GetNtotLinSysIts();

        m_TotImpStages++;
        m_StagesPerStep++;
    }

    bool CompressibleFlowSystem::UpdatePreconMatCheck(
        const TensorOfArray1D<NekDouble>  &res)
    {
        boost::ignore_unused(res);

        bool flag = false;

        if(0 < m_JFNKPrecondStep)
        {
            if (m_cflLocTimestep > 0.0)
            {
                flag = true;
            }
            else
            {
                if((m_CalcPreconMatFlag) ||
                    (m_TimeIntegLambdaPrcMat!=m_TimeIntegLambda))
                {
                    flag = true;
                    m_CalcPreconMatNumbSteps = m_CalcPreconMatCounter;
                }
            }
        }
        return flag;
    }

    void CompressibleFlowSystem::CalcRefValues(
            const TensorOfArray1D<NekDouble>    &inarray)
    {
        unsigned int nvariables         = m_fields.size();
        unsigned int ntotal             = inarray.size();
        unsigned int npoints            = ntotal/nvariables;

        unsigned int ntotalGlobal       = ntotal;
        m_comm->AllReduce(ntotalGlobal, Nektar::LibUtilities::ReduceSum);
        unsigned int ntotalDOF          = ntotalGlobal / nvariables;
        NekDouble ototalDOF             = 1.0 / ntotalDOF;

        m_inArrayNorm = 0.0;
        m_magnitdEstimat = Array<OneD, NekDouble>  (nvariables, 0.0);

        for (int i = 0; i < nvariables; ++i)
        {
            int offset = i * npoints;
            m_magnitdEstimat[i] = Vmath::Dot(npoints, inarray + offset,
                                            inarray + offset);
        }
        m_comm->AllReduce(m_magnitdEstimat, Nektar::LibUtilities::ReduceSum);

        for (int i = 0; i < nvariables; ++i)
        {
            m_inArrayNorm += m_magnitdEstimat[i];
        }

        for (int i = 2; i < nvariables - 1; ++i)
        {
            m_magnitdEstimat[1]   +=   m_magnitdEstimat[i] ;
        }
        for (int i = 2; i < nvariables - 1; ++i)
        {
            m_magnitdEstimat[i]   =   m_magnitdEstimat[1] ;
        }

        for (int i = 0; i < nvariables; ++i)
        {
            m_magnitdEstimat[i] = sqrt(m_magnitdEstimat[i] * ototalDOF);
        }
        if (m_root && m_verbose)
        {
            for (int i = 0; i < nvariables; ++i)
            {
                cout << "m_magnitdEstimat[" << i << "]    = "
                     << m_magnitdEstimat[i] << endl;
            }
            cout << "m_inArrayNorm    = " << m_inArrayNorm << endl;
        }
    }

    void CompressibleFlowSystem::MatrixMultiplyMatrixFreeCoeff(
            const TensorOfArray1D<NekDouble>    &inarray,
            TensorOfArray1D<NekDouble>          &out,
            const bool                          &flag)
    {
        boost::ignore_unused(flag);
        const Array<OneD, const NekDouble> refsol =
            m_nonlinsol->GetRefSolution();
        const Array<OneD, const NekDouble> refres =
            m_nonlinsol->GetRefResidual();
        const Array<OneD, const NekDouble> refsource =
            m_nonlinsol->GetRefSourceVec();

        NekDouble eps = m_JacobiFreeEps;
        unsigned int ntotalGlobal     = inarray.size();
        NekDouble magninarray = Vmath::Dot(ntotalGlobal,inarray,inarray);
        m_comm->AllReduce(magninarray, Nektar::LibUtilities::ReduceSum);
        eps *= sqrt( (sqrt(m_inArrayNorm) + 1.0)/magninarray);

        NekDouble oeps = 1.0 / eps;
        unsigned int ntotal     = inarray.size();
        TensorOfArray1D<NekDouble> solplus{ntotal, 0.0};
        TensorOfArray1D<NekDouble> resplus{ntotal, 0.0};

        Vmath::Svtvp(ntotal, eps, inarray, 1, refsol, 1, solplus, 1);

        NonlinSysEvaluatorCoeff(solplus, resplus, flag, refsource);
        
        Vmath::Vsub(ntotal, resplus, 1, refres, 1, out,1);
        Vmath::Smul(ntotal, oeps, out, 1, out, 1);
       
        return;
    }

    void CompressibleFlowSystem::SetBoundaryConditions(
            TensorOfArray2D<NekDouble>             &physarray,
            NekDouble                                         time)
    {
        int nTracePts  = GetTraceTotPoints();
        int nvariables = physarray.size();

        TensorOfArray2D<NekDouble> Fwd(nvariables);
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
        const TensorOfArray2D<NekDouble>  &physfield,
        TensorOfArray3D<NekDouble>                  &flux)
    {
        int i, j;
        int nq = physfield[0].size();
        int nVariables = m_fields.size();

        Array<OneD, NekDouble> pressure(nq);
        TensorOfArray2D<NekDouble> velocity(m_spacedim);

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
        const TensorOfArray2D<NekDouble>      &physfield,
        TensorOfArray3D<NekDouble>                      &flux)
    {
        int i, j;
        int nq = physfield[0].size();
        int nVariables = m_fields.size();

        // Factor to rescale 1d points in dealiasing
        NekDouble OneDptscale = 2;
        nq = m_fields[0]->Get1DScaledTotPoints(OneDptscale);

        Array<OneD, NekDouble> pressure(nq);
        TensorOfArray2D<NekDouble> velocity(m_spacedim);

        TensorOfArray2D<NekDouble> physfield_interp(nVariables);
        TensorOfArray3D<NekDouble> flux_interp(nVariables);

        for (i = 0; i < nVariables; ++ i)
        {
            physfield_interp[i] = Array<OneD, NekDouble>(nq);
            flux_interp[i] = TensorOfArray2D<NekDouble>(m_spacedim);
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
    Array<OneD, NekDouble> CompressibleFlowSystem::v_GetMaxStdVelocity(
        const NekDouble SpeedSoundFactor)
    {
        int nTotQuadPoints = GetTotPoints();
        int n_element      = m_fields[0]->GetExpSize();
        int expdim         = m_fields[0]->GetGraph()->GetMeshDimension();
        int nfields        = m_fields.size();
        int offset;
        Array<OneD, NekDouble> tmp;

        TensorOfArray2D<NekDouble> physfields(nfields);
        for (int i = 0; i < nfields; ++i)
        {
            physfields[i] = m_fields[i]->GetPhys();
        }

        Array<OneD, NekDouble> stdV(n_element, 0.0);

        // Getting the velocity vector on the 2D normal space
        TensorOfArray2D<NekDouble> velocity   (m_spacedim);
        TensorOfArray2D<NekDouble> stdVelocity(m_spacedim);
        TensorOfArray2D<NekDouble> stdSoundSpeed(m_spacedim);
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
                          SpeedSoundFactor * 
                          std::abs(stdSoundSpeed[j][offset + i]);
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
            NEKERROR(ErrorUtil::efatal, "Continuous Galerkin stability "
                "coefficients not introduced yet.");
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
            TensorOfArray2D<NekDouble> tmp(m_fields.size());

            for (int i = 0; i < m_fields.size(); ++i)
            {
                tmp[i] = m_fields[i]->GetPhys();
            }

            TensorOfArray2D<NekDouble> velocity(m_spacedim);
            TensorOfArray2D<NekDouble> velFwd  (m_spacedim);
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

            TensorOfArray2D<NekDouble> velocities(m_spacedim);
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

    /**
     *
     */
    void CompressibleFlowSystem::GetPressure(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD, NekDouble>                     &pressure)
    {
        m_varConv->GetPressure(physfield, pressure);
    }

    /**
     *
     */
    void CompressibleFlowSystem::GetDensity(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              Array<OneD, NekDouble>                     &density)
    {
        density = physfield[0];
    }

    /**
     *
     */
    void CompressibleFlowSystem::GetVelocity(
        const Array<OneD, const Array<OneD, NekDouble> > &physfield,
              TensorOfArray2D<NekDouble>       &velocity)
    {
        m_varConv->GetVelocityVector(physfield, velocity);
    }

    void CompressibleFlowSystem::v_SteadyStateResidual(
                int                         step, 
                Array<OneD, NekDouble>      &L2)
    {
        boost::ignore_unused(step);
        const int nPoints = GetTotPoints();
        const int nFields = m_fields.size();
        TensorOfArray2D<NekDouble> rhs (nFields);
        TensorOfArray2D<NekDouble> inarray (nFields);
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
            const MultiRegions::ExpListSharedPtr             &explist,
            const Array<OneD, const Array<OneD, NekDouble> > &normals,
            const int                                        nDervDir,
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
            TensorOfArray5D<NekDouble>                       &ElmtJacArray,
            const int                                        nfluxDir)
    {
        boost::ignore_unused(explist, normals, nDervDir, inarray, ElmtJacArray,
            nfluxDir);
        NEKERROR(ErrorUtil::efatal, "v_GetFluxDerivJacDirctn not coded");
    }

    void CompressibleFlowSystem::v_GetFluxDerivJacDirctnElmt(
            const int                        nConvectiveFields,
            const int                        nElmtPnt,
            const int                        nDervDir,
            const TensorOfArray2D<NekDouble> &locVars,
            const Array<OneD, NekDouble>     &locmu,
            const TensorOfArray2D<NekDouble> &locnormal,
            DNekMatSharedPtr                 &wspMat,
            TensorOfArray2D<NekDouble>       &PntJacArray)
    {
        boost::ignore_unused(nConvectiveFields, nElmtPnt, nDervDir, locVars,
            locmu, locnormal, wspMat, PntJacArray);
        NEKERROR(ErrorUtil::efatal, "v_GetFluxDerivJacDirctn not coded");
    }
    
    void CompressibleFlowSystem::v_GetFluxDerivJacDirctn(
            const MultiRegions::ExpListSharedPtr            &explist,
            const Array<OneD, const Array<OneD, NekDouble>> &normals,
            const int                                       nDervDir,
            const Array<OneD, const Array<OneD, NekDouble>> &inarray,
            TensorOfArray2D<DNekMatSharedPtr>               &ElmtJac)
    {
        boost::ignore_unused(explist, normals, nDervDir, inarray, ElmtJac);
    }

    void CompressibleFlowSystem::v_GetDiffusionFluxJacPoint(
            const Array<OneD, NekDouble>                        &conservVar, 
            const Array<OneD, const Array<OneD, NekDouble> >    &conseDeriv, 
            const NekDouble                                     mu,
            const NekDouble                                     DmuDT,
            const Array<OneD, NekDouble>                        &normals,
                  DNekMatSharedPtr                              &fluxJac)
    {
        boost::ignore_unused(conservVar, conseDeriv, mu, DmuDT, normals, 
            fluxJac);
   
    }

    void CompressibleFlowSystem::v_MinusDiffusionFluxJacDirctnElmt(
            const int                           nConvectiveFields,
            const int                           nElmtPnt,
            const TensorOfArray2D<NekDouble>    &locVars,
            const TensorOfArray3D<NekDouble>    &locDerv,
            const Array<OneD, NekDouble>        &locmu,
            const Array<OneD, NekDouble>        &locDmuDT,
            const Array<OneD, NekDouble>        &normals,
            DNekMatSharedPtr                    &wspMat,
            TensorOfArray2D<NekDouble>          &PntJacArray)
    {
        boost::ignore_unused(nConvectiveFields, nElmtPnt, locVars, locDerv,
                locmu, locDmuDT, normals, wspMat, PntJacArray);
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
                for (int i = 0; i < exp3D->GetNtraces(); ++i)
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
                for (int i = 0; i < exp2D->GetNtraces(); ++i)
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
                NEKERROR(ErrorUtil::efatal,"Dimension out of bound.")
            }
        }

        // Determine h/p scaling
        hOverP[e] = h/max(pOrderElmt[e]-1,1);

    }
    return hOverP;
}
}
