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
// Description: Compressible flow system base class with auxiliary functions
//
///////////////////////////////////////////////////////////////////////////////

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

        for (int i = 0; i < m_fields.num_elements(); i++)
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
            m_artificialDiffusion = GetArtificialDiffusionFactory()
                                    .CreateInstance(m_shockCaptureType,
                                                    m_session,
                                                    m_fields,
                                                    m_spacedim);
        }

        // Forcing terms for the sponge region
        m_forcing = SolverUtils::Forcing::Load(m_session, m_fields,
                                               m_fields.num_elements());

        // User-defined boundary conditions
        int cnt = 0;
        for (int n = 0; n < m_fields[0]->GetBndConditions().num_elements(); ++n)
        {
            std::string type =
                m_fields[0]->GetBndConditions()[n]->GetUserDefined();

            if (m_fields[0]->GetBndConditions()[n]->GetBoundaryConditionType()
                == SpatialDomains::ePeriodic)
            {
                continue;
            }

            if(!type.empty())
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
            
            
            //ASSERTL0(false, "Implicit CFS not set up.");
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
        if(m_useFiltering)
        {
            m_session->LoadParameter ("FilterAlpha", m_filterAlpha, 36);
            m_session->LoadParameter ("FilterExponent", m_filterExponent, 16);
            m_session->LoadParameter ("FilterCutoff", m_filterCutoff, 0);
        }

        // Load CFL for local time-stepping (for steady state)
        m_session->MatchSolverInfo("LocalTimeStep","True",
                                   m_useLocalTimeStep, false);
        if(m_useLocalTimeStep)
        {
            ASSERTL0(m_cflSafetyFactor != 0,
                    "Local time stepping requires CFL parameter.");
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
        int i;
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

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
            for(i = 0; i < nvariables; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }
        
        // Calculate advection
        DoAdvection(inarray, outarray, time, Fwd, Bwd);

        // Negate results
        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(npoints, outarray[i], 1);
        }

        // Add diffusion terms
        DoDiffusion(inarray, outarray, Fwd, Bwd);

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
                for(i = 0; i < nvariables; ++i)
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

    void CompressibleFlowSystem::preconditioner_BlkDiag(
                                                 const Array<OneD, NekDouble> &inarray,
                                                 Array<OneD, NekDouble >&outarray)
    {
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        Array<OneD, NekVector<NekDouble> > tmparray(nvariables);


        PointerWrapper pwrapp = eWrapper;
        if(inarray.get() == outarray.get())
        {
            pwrapp = eCopy;
        }
      
        for(int m = 0; m < nvariables; m++)
        {
            int moffset = m*npoints;
            tmparray[m] =  NekVector<NekDouble> (npoints,inarray+moffset,pwrapp);
        }

        Vmath::Fill(outarray.num_elements(),0.0,outarray,1);

        for(int m = 0; m < nvariables; m++)
        {
            int moffset = m*npoints;
            NekVector<NekDouble> out(npoints,outarray+moffset,eWrapper);
            for(int n = 0; n < nvariables; n++)
            {
                out += (*m_PrecMatVars[m][n])*tmparray[n];
            }
        }
    }

    void CompressibleFlowSystem::preconditioner_BlkSOR_coeff(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble> &outarray)
    {
        int nSORTot   =   m_JFNKPrecondStep;
        const NekDouble SORParam        =   1.0;
        const NekDouble OmSORParam      =   1.0-SORParam;
        
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        unsigned int ntotpnt    = inarray.num_elements();

        ASSERTL0(nvariables*npoints==ntotpnt,"nvariables*npoints==ntotpnt not satisfied in preconditioner_BlkSOR");


        Array<OneD, NekDouble> rhs;
        if(inarray.get() == outarray.get())
        {
            rhs =   Array<OneD, NekDouble>(ntotpnt,0.0);
            Vmath::Vcopy(ntotpnt,inarray,1,rhs,1);
        }
        else
        {
            rhs =   inarray;
        }

        Array<OneD, NekDouble>  outN(ntotpnt,0.0);
        Array<OneD, Array<OneD, NekDouble> >rhs2d(nvariables),outN2d(nvariables);
        for(int m = 0; m < nvariables; m++)
        {
            int moffset = m*npoints;
            rhs2d[m]    =   rhs     + moffset;
            outN2d[m]   =   outN    + moffset;
        }

        Vmath::Fill(ntotpnt,0.0,outN,1);

        preconditioner_BlkDiag(inarray,outN);

        PointerWrapper pwrapp = eWrapper;
        Array<OneD, NekVector<NekDouble> > tmparray(nvariables);
        for(int m = 0; m < nvariables; m++)
        {
            int moffset = m*npoints;
            tmparray[m] =  NekVector<NekDouble> (npoints,outN+moffset,pwrapp);
        }

        for(int nsor = 0; nsor < nSORTot-1; nsor++)
        {
            Vmath::Smul(ntotpnt,OmSORParam,outN,1,outarray,1);
            
            MinusOffDiag2Rhs(nvariables,npoints,rhs2d,outN2d);

            for(int m = 0; m < nvariables; m++)
            {
                int moffset = m*npoints;
                NekVector<NekDouble> out(npoints,outarray+moffset,eWrapper);
                for(int n = 0; n < nvariables; n++)
                {
                    int noffset = n*npoints;
                    tmparray[n] =  NekVector<NekDouble> (npoints,outN+noffset,pwrapp);
                    out += (*m_PrecMatVars[m][n])*(tmparray[n]); //SORParam
                }
            }
            Vmath::Vcopy(ntotpnt,outarray,1,outN,1);
        }
    }

    void CompressibleFlowSystem::MinusOffDiag2Rhs(
            const int nvariables,
            const int nCoeffs,
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray)
    {
        int nTracePts  = GetTraceTotPoints();
        int npoints    = GetNpoints();

        Array<OneD, Array<OneD, NekDouble> > outpnts(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            outpnts[i]  =  Array<OneD, NekDouble> (npoints,0.0);
            m_fields[i]->BwdTrans(outarray[i],outpnts[i]);
        }

        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble> > Fwd        (nvariables);
        Array<OneD, Array<OneD, NekDouble> > Bwd        (nvariables);
        Array<OneD, Array<OneD, NekDouble> > FwdFlux  (nvariables);
        Array<OneD, Array<OneD, NekDouble> > BwdFlux  (nvariables);

        // for(int i = 0; i < nvariables; ++i)
        // {
        //     Fwd[i]          = Array<OneD, NekDouble>(nTracePts, -100000000000.0);
        //     Bwd[i]          = Array<OneD, NekDouble>(nTracePts, -100000000000.0);
        //     FwdFlux[i]      = Array<OneD, NekDouble>(nTracePts, 0.0);
        //     BwdFlux[i]      = Array<OneD, NekDouble>(nTracePts, 0.0);
        //     m_fields[i]->GetFwdBwdTracePhys_singlethread(outpnts[i], Fwd[i], Bwd[i]);
        //     // m_fields[i]->GetFwdBwdTracePhys(outpnts[i], Fwd[i], Bwd[i]);
        // }

        // for(int j=0;j<nTracePts;j++)
        // {
        //     if(Bwd[0][j]<-1000000.0)
        //     {
        //         for(int i = 0; i < nvariables; ++i)
        //         {
        //             Bwd[i][j] = Fwd[i][j];
        //         }

        //     }
        // }


        for(int i = 0; i < nvariables; ++i)
        {
            Fwd[i]          = Array<OneD, NekDouble>(nTracePts, 0.0);
            Bwd[i]          = Array<OneD, NekDouble>(nTracePts, 0.0);
            FwdFlux[i]      = Array<OneD, NekDouble>(nTracePts, 0.0);
            BwdFlux[i]      = Array<OneD, NekDouble>(nTracePts, 0.0);
            // m_fields[i]->GetFwdBwdTracePhys_singlethread(outpnts[i], Fwd[i], Bwd[i]);
            m_fields[i]->GetFwdBwdTracePhys(outpnts[i], Fwd[i], Bwd[i]);
        }

            // m_fields[i]->GetFwdBwdTracePhys(outpnts[i], Fwd[i], Bwd[i]);
        
        NekVector<NekDouble> VFwd(nvariables),VBwd(nvariables);
        NekVector<NekDouble> VFlux(nvariables);
        DNekMatSharedPtr    PJacFwd,PJacBwd;

        for(int n = 0; n < nTracePts; ++n)
        {
            for(int i = 0; i < nvariables; ++i)
            {
                VFwd[i] =   Fwd[i][n];
                VBwd[i] =   Bwd[i][n];
            }   
            PJacFwd =  m_TraceJac[0]->GetBlock(n,n); 
            PJacBwd =  m_TraceJac[1]->GetBlock(n,n); 

            VFlux   =   (*PJacFwd)*VFwd;

            for(int i = 0; i < nvariables; ++i)
            {
                FwdFlux[i][n] =   VFlux[i];
            }  

            VFlux   =   (*PJacBwd)*VBwd;

            for(int i = 0; i < nvariables; ++i)
            {
                BwdFlux[i][n] =   VFlux[i];
            }  
        }



        // // Store forwards/backwards space along trace space
        // Array<OneD, Array<OneD, NekDouble> > Fwd0        (nvariables);
        // Array<OneD, Array<OneD, NekDouble> > Bwd0        (nvariables);
        // Array<OneD, Array<OneD, NekDouble> > FwdFlux0  (nvariables);
        // Array<OneD, Array<OneD, NekDouble> > BwdFlux0  (nvariables);

        // for(int i = 0; i < nvariables; ++i)
        // {
        //     Fwd0[i]          = Array<OneD, NekDouble>(nTracePts, 0.0);
        //     Bwd0[i]          = Array<OneD, NekDouble>(nTracePts, 0.0);
        //     FwdFlux0[i]      = Array<OneD, NekDouble>(nTracePts, 0.0);
        //     BwdFlux0[i]      = Array<OneD, NekDouble>(nTracePts, 0.0);
        //     m_fields[i]->GetFwdBwdTracePhys_singlethread(outpnts[i], Fwd0[i], Bwd0[i]);
        //     // m_fields[i]->GetFwdBwdTracePhys(outpnts[i], Fwd[i], Bwd[i]);
        // }

        // for(int n = 0; n < nTracePts; ++n)
        // {
        //     for(int i = 0; i < nvariables; ++i)
        //     {
        //         VFwd[i] =   Fwd0[i][n];
        //         VBwd[i] =   Bwd0[i][n];
        //     }   
        //     PJacFwd =  m_TraceJac[0]->GetBlock(n,n); 
        //     PJacBwd =  m_TraceJac[1]->GetBlock(n,n); 

        //     VFlux   =   (*PJacFwd)*VFwd;

        //     for(int i = 0; i < nvariables; ++i)
        //     {
        //         FwdFlux0[i][n] =   VFlux[i];
        //     }  

        //     VFlux   =   (*PJacBwd)*VBwd;

        //     for(int i = 0; i < nvariables; ++i)
        //     {
        //         BwdFlux0[i][n] =   VFlux[i];
        //     }  
        // }

        // if(0==m_session->GetComm()->GetRank())
        // {
        //     int nwidthcolm = 20;
        //     for(int n = 0; n < nTracePts; ++n)
        //     {
        //         cout <<endl<<" /*///////*/*/*/*********************/*/*/*/**/*/*/**///**/*//******/*= "<<endl;
        //         cout <<std::right <<std::scientific<<std::setw(nwidthcolm)<<std::setprecision(nwidthcolm-8)
        //             <<" nTracePts= " <<n<<endl;
        //         for(int i = 0; i < nvariables; ++i)
        //         {
        //             cout    << "delt    = "<<abs(Fwd[i][n]      - Fwd0[i][n]    )   << "    Fwd    = "<<Fwd[i][n]       << "    Fwd0    = "<<Fwd0[i][n]     <<endl
        //                     << "delt    = "<<abs(Bwd[i][n]      - Bwd0[i][n]    )   << "    Bwd    = "<<Bwd[i][n]       << "    Bwd0    = "<<Bwd0[i][n]     <<endl
        //                     << "delt    = "<<abs(FwdFlux[i][n]  - FwdFlux0[i][n])   << "    FwdFlux= "<<FwdFlux[i][n]   << "    FwdFlux0= "<<FwdFlux0[i][n] <<endl
        //                     << "delt    = "<<abs(BwdFlux[i][n]  - BwdFlux0[i][n])   << "    BwdFlux= "<<BwdFlux[i][n]   << "    BwdFlux0= "<<BwdFlux0[i][n] <<endl;
        //         }
        //     }
        // }

        // int tmpcin;
        // cin >> tmpcin;

        // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]
        for(int i = 0; i < nvariables; ++i)
        {
            Vmath::Fill(nCoeffs,0.0,outarray[i],1);
            // Vmath::Neg                                  (nCoeffs, outarray[i], 1);
            m_fields[i]->AddTraceIntegral2OffDiag       (FwdFlux[i],BwdFlux[i], outarray[i]);
            m_fields[i]->MultiplyByElmtInvMass          (outarray[i], outarray[i]);
        }

        for(int i = 0; i < nvariables; ++i)
        {
            Vmath::Svtvp(nCoeffs,-m_TimeIntegLambda,outarray[i],1,inarray[i],1,outarray[i],1);
        }
    }
    
    void CompressibleFlowSystem::ElmtVarInvMtrx(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray)
    {
        int n1d = gmtxarray.num_elements();
        int n2d = gmtxarray[0].num_elements();
        int nConvectiveFields = n1d;

        ASSERTL0(n1d==n2d,"ElmtVarInvMtrx requires n1d==n2d");

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;

        gmtxarray[0][0]->GetBlockSizes(rowSizes,colSizes);
        int ntotElmt  = rowSizes.num_elements();
        int nElmtCoef   =    rowSizes[0]-1;
        int nElmtCoef0  =    nElmtCoef;

        int nElmtCoefVAr = nElmtCoef0*nConvectiveFields;
        DNekMatSharedPtr    tmpGmtx = MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nElmtCoefVAr, nElmtCoefVAr,0.0);

        DNekMatSharedPtr        ElmtMat;

        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            ElmtMat =   gmtxarray[0][0]->GetBlock(nelmt,nelmt);

            int nrows = ElmtMat->GetRows();
            int ncols = ElmtMat->GetColumns();
            ASSERTL0(nrows==ncols,"ElmtVarInvMtrx requires nrows==ncols");

            nElmtCoef            = nrows;
            
            if (nElmtCoef0!=nElmtCoef) 
            {
                nElmtCoef0 = nElmtCoef;
                nElmtCoefVAr = nElmtCoef0*nConvectiveFields;
                tmpGmtx = MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nElmtCoefVAr, nElmtCoefVAr,0.0);
            }

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    ElmtMat =   gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                    for(int nrw = 0; nrw < nElmtCoef; nrw++)
                    {
                        int nrwvar = m*nElmtCoef+nrw;
                        for(int ncl = 0; ncl < nElmtCoef; ncl++)
                        {
                            int nclvar = n*nElmtCoef+ncl;
                            (*tmpGmtx)(nrwvar,nclvar)=(*ElmtMat)(nrw,ncl);
                        }
                    }
                }
            }

            tmpGmtx->Invert();

            for(int m = 0; m < nConvectiveFields; m++)
            {
                for(int n = 0; n < nConvectiveFields; n++)
                {
                    ElmtMat =   gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                    for(int nrw = 0; nrw < nElmtCoef; nrw++)
                    {
                        int nrwvar = m*nElmtCoef+nrw;
                        for(int ncl = 0; ncl < nElmtCoef; ncl++)
                        {
                            int nclvar = n*nElmtCoef+ncl;
                            (*ElmtMat)(nrw,ncl) =   (*tmpGmtx)(nrwvar,nclvar);
                        }
                    }
                }
            }
        }
        return;
    }

    void CompressibleFlowSystem::AddMatNSBlkDiag_volume(
                                        const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray)
    {
        int nSpaceDim = m_graph->GetSpaceDimension();
        int nvariable = inarray.num_elements();



        Array<OneD, Array<OneD, DNekMatSharedPtr> > ElmtJac;

        
        for(int nDirctn = 0; nDirctn < nSpaceDim; nDirctn++)
        {
            ElmtJac = GetFluxVectorJacDirctn(nDirctn,inarray);

            m_advObject->AddVolumJac2Mat(nvariable,m_fields,ElmtJac,nDirctn,gmtxarray);
        }
    }


    void CompressibleFlowSystem::AddMatNSBlkDiag_boundary(const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,
                                        Array<OneD, DNekBlkMatSharedPtr > &TraceJac)
    {
        int nvariables = inarray.num_elements();
        GetTraceJac(inarray,TraceJac);
        m_advObject->AddTraceJac2Mat(nvariables,m_fields, TraceJac,gmtxarray);
    }

    void CompressibleFlowSystem::GetTraceJac(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, DNekBlkMatSharedPtr > &TraceJac)
    {
        int nvariables = inarray.num_elements();
        // int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble> > Fwd    (nvariables);
        Array<OneD, Array<OneD, NekDouble> > Bwd    (nvariables);

        Array<OneD, unsigned int> n_blks(nTracePts);
        for(int i=0;i<nTracePts;i++)
        {
            n_blks[i]    = nvariables;
        }
        DNekBlkMatSharedPtr FJac = MemoryManager<DNekBlkMat>
            ::AllocateSharedPtr(n_blks, n_blks, eDIAGONAL);
        DNekBlkMatSharedPtr BJac = MemoryManager<DNekBlkMat>
            ::AllocateSharedPtr(n_blks, n_blks, eDIAGONAL);

        if (m_HomogeneousType == eHomogeneous1D)
        {
            Fwd = NullNekDoubleArrayofArray;
            Bwd = NullNekDoubleArrayofArray;
        }
        else
        {
            for(int i = 0; i < nvariables; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }

        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);

        // Riemann flux Jacobian considering implicit boundary conditions
        NumCalRiemFluxJac(nvariables, m_fields, advVel, inarray,
                            m_BndEvaluateTime, Fwd, Bwd, FJac, BJac);

        // Riemann flux Jacobian without considering implicit boundary conditions
        // m_advObject->NumCalRiemFluxJac(nvariables, m_fields, advVel, inarray,
        //                     Fwd, Bwd, FJac, BJac);

        TraceJac    =   Array<OneD, DNekBlkMatSharedPtr>(2);
        TraceJac[0] = FJac;
        TraceJac[1] = BJac;
    }


    void CompressibleFlowSystem::NumCalRiemFluxJac(
        const int                                          nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &AdvVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const NekDouble                                   &time,
        const Array<OneD, Array<OneD, NekDouble> >        &Fwd,
        const Array<OneD, Array<OneD, NekDouble> >        &Bwd,
        DNekBlkMatSharedPtr &FJac,
        DNekBlkMatSharedPtr &BJac)
    {
        int nvariables = inarray.num_elements();
        int npoints    = GetNpoints();
        int nTracePts  = GetTraceTotPoints();

        Array<OneD,       Array<OneD, NekDouble> >  tmpinarry(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            tmpinarry[i]    =    Array<OneD, NekDouble>(npoints,0.0);
            Vmath::Vcopy(npoints, inarray[i],1,tmpinarry[i],1);
        }

        // Allocate the Jacobian matrix
        DNekMatSharedPtr tmpMat;
        for(int i=0;i<nTracePts;i++)
        {
            tmpMat = MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nvariables, nvariables);
            FJac->SetBlock(i,i,tmpMat);

            tmpMat = MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nvariables, nvariables);
            BJac->SetBlock(i,i,tmpMat);
        }


        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);

        Array<OneD, Array<OneD, NekDouble> > numflux(nvariables);
        for(int i = 0; i < nvariables; ++i)
        {
            numflux[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        }

        m_advObject->NumericalFlux(nvariables, m_fields, advVel, inarray, numflux,
                            m_BndEvaluateTime,Fwd, Bwd);
        
        NekDouble eps   =   1.0E-6;

        int nFields = nvariables;
        int nPts    = nTracePts;
        
             
        // Allocate temporary variables
        Array<OneD,       Array<OneD, NekDouble> >  plusFwd(nFields),plusBwd(nFields);
        Array<OneD,       Array<OneD, NekDouble> >  plusflux(nFields),Jacvect(nFields);
        Array<OneD,       Array<OneD, NekDouble> >  FwdBnd(nFields);
        for(int i = 0; i < nFields; i++)
        {
            Jacvect[i]      =    Array<OneD, NekDouble>(nPts,0.0);
            plusFwd[i]      =    Array<OneD, NekDouble>(nPts,0.0);
            plusBwd[i]      =    Array<OneD, NekDouble>(nPts,0.0);
            plusflux[i]     =    Array<OneD, NekDouble>(nPts,0.0);
            FwdBnd[i]       =    Array<OneD, NekDouble>(nPts,0.0);
        }


        for(int i = 0; i < nFields; i++)
        {
            Vmath::Vcopy(nPts, Fwd[i],1,plusFwd[i],1);
            Vmath::Vcopy(nPts, Bwd[i],1,plusBwd[i],1);
        }

        // Fwd Jacobian
        for(int i = 0; i < nFields; i++)
        {
            NekDouble epsvar = eps*m_magnitdEstimat[i];
            NekDouble oepsvar   =   1.0/epsvar;
            Vmath::Sadd(nPts,epsvar,Fwd[i],1,plusFwd[i],1);
            
            if (m_bndConds.size())
            {
                for(int i = 0; i < nFields; i++)
                {
                    Vmath::Vcopy(nPts, plusFwd[i],1,FwdBnd[i],1);
                }
                // Loop over user-defined boundary conditions
                for (auto &x : m_bndConds)
                {
                    x->Apply(FwdBnd, tmpinarry, time);
                }
            }

            for(int j = 0; j < nFields; j++)
            {
                m_fields[j]->FillBwdWITHBound(plusFwd[j], plusBwd[j]);
            }

            m_advObject->NumericalFlux(nvariables, m_fields, advVel, inarray, plusflux,
                            m_BndEvaluateTime,plusFwd, plusBwd);

            for (int n = 0; n < nFields; n++)
            {
                Vmath::Vsub(nPts,&plusflux[n][0],1,&numflux[n][0],1,&Jacvect[n][0],1);
                Vmath::Smul(nPts, oepsvar ,&Jacvect[n][0],1,&Jacvect[n][0],1);
            }
            for(int j = 0; j < nPts; j++)
            {
                tmpMat  =   FJac->GetBlock(j,j);
                for (int n = 0; n < nFields; n++)
                {
                    (*tmpMat)(n,i) = Jacvect[n][j];
                }
            }
            
            Vmath::Vcopy(nPts, Fwd[i],1,plusFwd[i],1);
        }

        // Reset the boundary conditions
        if (m_bndConds.size())
        {
            for(int i = 0; i < nFields; i++)
            {
                Vmath::Vcopy(nPts, Fwd[i],1,FwdBnd[i],1);
            }
            // Loop over user-defined boundary conditions
            for (auto &x : m_bndConds)
            {
                x->Apply(FwdBnd, tmpinarry, time);
            }
        }

        for(int i = 0; i < nFields; i++)
        {
            Vmath::Vcopy(nPts, Bwd[i],1,plusBwd[i],1);
        }

        for(int i = 0; i < nFields; i++)
        {
            NekDouble epsvar = eps*m_magnitdEstimat[i];
            NekDouble oepsvar   =   1.0/epsvar;

            Vmath::Sadd(nPts,epsvar,Bwd[i],1,plusBwd[i],1);

            for(int j = 0; j < nFields; j++)
            {
                m_fields[j]->FillBwdWITHBound(Fwd[j], plusBwd[j]);
            }

            m_advObject->NumericalFlux(nvariables, m_fields, advVel, inarray, plusflux,
                            m_BndEvaluateTime,Fwd, plusBwd);

            for (int n = 0; n < nFields; n++)
            {
                Vmath::Vsub(nPts,&plusflux[n][0],1,&numflux[n][0],1,&Jacvect[n][0],1);
                Vmath::Smul(nPts, oepsvar ,&Jacvect[n][0],1,&Jacvect[n][0],1);
            }
            for(int j = 0; j < nPts; j++)
            {
                tmpMat  =   BJac->GetBlock(j,j);
                for (int n = 0; n < nFields; n++)
                {
                    (*tmpMat)(n,i) = Jacvect[n][j];
                }
            }
            
            Vmath::Vcopy(nPts, Bwd[i],1,plusBwd[i],1);
        }

        // int nwidthcolm =23;
        // for(int i = 0; i < nvariables; i++)
        // {
        //     for(int j = 0; j < Fwd[i].num_elements(); j++)
        //     {
        //         cout    << "plusFwd["<<i<<"]["<<j<<"]"
        //                 <<std::scientific<<std::setw(nwidthcolm)<<std::setprecision(nwidthcolm-8)
        //                 <<Fwd[i][j]<<endl;
        //     }

        //     for(int j = 0; j < Bwd[i].num_elements(); j++)
        //     {
        //         cout    << "plusBwd["<<i<<"]["<<j<<"]"
        //                 <<std::scientific<<std::setw(nwidthcolm)<<std::setprecision(nwidthcolm-8)
        //                 <<Bwd[i][j]<<endl;
        //     }
        // }
    }

#ifdef DEMO_IMPLICITSOLVER_JFNK_COEFF

    void CompressibleFlowSystem::DoImplicitSolve_phy2coeff(
                                                 const Array<OneD, const Array<OneD, NekDouble> >&inpnts,
                                                       Array<OneD,       Array<OneD, NekDouble> >&outpnt,
                                                 const NekDouble time,
                                                 const NekDouble lambda)
    {
        
        unsigned int nvariables  = inpnts.num_elements();
        unsigned int ncoeffs     = m_fields[0]->GetNcoeffs();

        Array<OneD, Array<OneD, NekDouble> > inarray(nvariables);
        Array<OneD, Array<OneD, NekDouble> > out(nvariables);
        
        for(int i = 0; i < nvariables; i++)
        {
            inarray[i]  =  Array<OneD, NekDouble> (ncoeffs,0.0);
            out[i]      =  Array<OneD, NekDouble> (ncoeffs,0.0);
            m_fields[i]->FwdTrans(inpnts[i],inarray[i]);
        }

        DoImplicitSolve_coeff(inpnts,inarray,out,time,lambda);

        for(int i = 0; i < nvariables; i++)
        {
            m_fields[i]->BwdTrans(out[i],outpnt[i]);
        }

    }
    
    void CompressibleFlowSystem::DoImplicitSolve_coeff(
                                                 const Array<OneD, const Array<OneD, NekDouble> >&inpnts,
                                                 const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                                       Array<OneD,       Array<OneD, NekDouble> >&out,
                                                 const NekDouble time,
                                                 const NekDouble lambda)
    {
        m_TimeIntegtSol_n   = inarray;
        m_TimeIntegtSol_k   = out;
        m_TimeIntegLambda   = lambda;
        m_BndEvaluateTime   = time;
        bool l_verbose      = m_session->DefinesCmdLineArgument("verbose");
        bool converged      = false;
        const unsigned int MaxNonlinIte =   100;
        unsigned int nvariables  = inarray.num_elements();
        unsigned int npoints     = inarray[0].num_elements();
        unsigned int ntotal      = nvariables*npoints;
        unsigned int nElements = m_fields[0]->GetExpSize();

        LibUtilities::CommSharedPtr v_Comm  = m_fields[0]->GetComm()->GetRowComm();

        bool l_root=false;
        if(0==v_Comm->GetRank())
        {
            l_root =true;
        }

        unsigned int ntotalGlobal     = ntotal;
        v_Comm->AllReduce(ntotalGlobal, Nektar::LibUtilities::ReduceSum);
        unsigned int ntotalDOF        = ntotalGlobal/nvariables;
        NekDouble ototalGlobal  = 1.0/ntotalGlobal; 
        NekDouble ototalDOF     = 1.0/ntotalDOF; 

        if(m_inArrayNorm<0.0)
        {
            m_inArrayNorm = 0.0;

            // estimate the magnitude of each flow variables 
            m_magnitdEstimat = Array<OneD, NekDouble>  (nvariables,0.0);

            for(int i = 0; i < nvariables; i++)
            {
                m_magnitdEstimat[i] = Vmath::Dot(npoints,inarray[i],inarray[i]);
            }
            v_Comm->AllReduce(m_magnitdEstimat, Nektar::LibUtilities::ReduceSum);

            
            // NekDouble ototpnts = 1.0/(ntotalGlobal/nvariables);

            for(int i = 0; i < nvariables; i++)
            {
                m_inArrayNorm += m_magnitdEstimat[i];
            }

            for(int i = 2; i < nvariables-1; i++)
            {
                m_magnitdEstimat[1]   +=   m_magnitdEstimat[i] ;
            }
            for(int i = 2; i < nvariables-1; i++)
            {
                m_magnitdEstimat[i]   =   m_magnitdEstimat[1] ;
            }
            
            for(int i = 0; i < nvariables; i++)
            {
                m_magnitdEstimat[i] = sqrt(m_magnitdEstimat[i]*ototalDOF);
            }
            if(0==m_session->GetComm()->GetRank()&&l_verbose)
            {
                for(int i = 0; i < nvariables; i++)
                {
                    cout << "m_magnitdEstimat["<<i<<"]    = "<<m_magnitdEstimat[i]<<endl;
                }
                cout << "m_inArrayNorm    = "<<m_inArrayNorm<<endl;
            }
        }


        NekDouble LinSysTol = 0.0;
        NekDouble tolrnc    = m_NewtonAbsoluteIteTol;
        NekDouble tol2      = m_inArrayNorm*tolrnc*tolrnc*ntotalDOF;
        NekDouble tol2Max   = sqrt(m_inArrayNorm*ototalDOF*tolrnc);
        // NekDouble ratioTol = m_NewtonRelativeIteTol;
        NekDouble ratioTol = tolrnc;
        //TODO: if steady state the dt should change according to resnorm/m_inArrayNorm
        // so that quadrature convergence rate may be achieved
        // at this situation && may be better in convergence critiria
        NekDouble tol2Ratio = ratioTol*ratioTol*ntotalDOF;

        m_PrecMatVars = Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            m_PrecMatVars[i] =  Array<OneD, DNekBlkMatSharedPtr> (nvariables);
        }
        AllocatePrecondBlkDiag_coeff(m_PrecMatVars);

        Array<OneD, NekDouble> NonlinSysRes_1D(ntotal,0.0),sol_k_1D(ntotal,0.0),dsol_1D(ntotal,0.0);
        Array<OneD,       Array<OneD, NekDouble> > NonlinSysRes(nvariables),dsol(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            int offset = i*npoints;
            NonlinSysRes[i] =  NonlinSysRes_1D + offset;
            dsol[i]         =  dsol_1D + offset;
        }

        m_SysEquatResid_k = NonlinSysRes;
        //sol_k = out;
        for(int i = 0; i < nvariables; i++)
        {
            Vmath::Vcopy(npoints,inarray[i],1,m_TimeIntegtSol_k[i],1);
        }

        Array<OneD, NekDouble> tstep (nElements, 0.0);
        NekDouble tmpcfl    = m_cflSafetyFactor;
        m_cflSafetyFactor   = m_cflLocTimestep;
        GetElmtTimeStep(inpnts, tstep);
        m_cflSafetyFactor   = tmpcfl;

        m_locTimeStep = tstep;


        // TODO: 
        // v_Comm is based on the explist used may be different(m_tracemap or m_locToGloMap) for diffrent
        // should generate GlobalLinSys first and get the v_Comm from GlobalLinSys. here just give it a m_Comm no parallel support yet!!
        //const std::weak_ptr<ExpList> 
        // std::shared_ptr<MultiRegions::ExpList> vExpList = m_fields[0]->GetSharedThisPtr();

        NekLinSysIterative linsol(m_session,v_Comm);
        m_LinSysOprtors.DefineMatrixMultiply(&CompressibleFlowSystem::MatrixMultiply_MatrixFree_coeff, this);
        m_LinSysOprtors.DefinePrecond(&CompressibleFlowSystem::preconditioner_BlkSOR_coeff, this);
        // m_LinSysOprtors.DefinePrecond(&CompressibleFlowSystem::preconditioner_BlkDiag, this);
        // m_LinSysOprtors.DefinePrecond(&CompressibleFlowSystem::preconditioner, this);
        linsol.setLinSysOperators(m_LinSysOprtors);

        // NonlinSysEvaluator_coeff(m_TimeIntegtSol_k,m_SysEquatResid_k);
        // DebugNumCalJac_coeff(m_PrecMatVars);
        // ElmtVarInvMtrx_coeff(m_PrecMatVars);

        // if(m_TotLinItePrecondMat>0)
        // {
            int nphspnt = inpnts[0].num_elements();
            Array<OneD, Array<OneD, NekDouble> > intmp(nvariables);
            for(int i = 0; i < nvariables; i++)
            {
                intmp[i]    =   Array<OneD, NekDouble>(nphspnt,0.0);
            }

            DoOdeProjection(inpnts,intmp,m_BndEvaluateTime);
            GetpreconditionerNSBlkDiag_coeff(intmp,m_PrecMatVars,m_TraceJac);
            // cout << "m_TotLinItePrecondMat  =   "<<m_TotLinItePrecondMat<<endl;
            // m_TotLinItePrecondMat = 0;

            // to free the storage
            for(int i = 0; i < nvariables; i++)
            {
                intmp[i]    =   Array<OneD, NekDouble>(1,0.0);
            }
        // }

        

        converged = false;
        int nwidthcolm = 8;
        int NtotDoOdeRHS = 0;
        NekDouble resnorm0 = 0.0;
        NekDouble resnorm  = 0.0;
        NekDouble resmaxm  = 0.0;
        NekDouble resratio = 1.0;
        for (int k = 0; k < MaxNonlinIte; k++)
        {
            NonlinSysEvaluator_coeff(m_TimeIntegtSol_k,m_SysEquatResid_k);
            // NonlinSysEvaluator_coeff_out(m_TimeIntegtSol_k,m_SysEquatResid_k);
            NtotDoOdeRHS++;
            // NonlinSysRes_1D and m_SysEquatResid_k share the same storage
            resnorm = Vmath::Dot(ntotal,NonlinSysRes_1D,NonlinSysRes_1D);

            v_Comm->AllReduce(resnorm, Nektar::LibUtilities::ReduceSum);

            if(0==k)
            {
                resnorm0 = resnorm;
                resratio = 1.0;
            }
            else
            {
                resratio = resnorm/resnorm0;
            }

            if (resnorm<tol2||resratio<tol2Ratio)
            {
                // this works when resnorm0<<m_inArrayNorm ie. close to converge
                if(resratio<0.1)
                {
                    resmaxm = 0.0;
                    for(int i=0;i<ntotal;i++)
                    {
                        resmaxm = max(resmaxm,abs(NonlinSysRes_1D[i]));
                    }
                    v_Comm->AllReduce(resmaxm, Nektar::LibUtilities::ReduceMax);
                    if(resmaxm<tol2Max)
                    {
                    // at least one Newton Iteration 
                    // or else the flow field will not update 
                    // if(k>0)
                    // {
                        /// TODO: m_root
                        if(l_verbose&&l_root)
                        {
                            cout <<right<<scientific<<setw(nwidthcolm)<<setprecision(nwidthcolm-6)
                                <<" * Newton-Its converged (RES=" 
                                << sqrt(resnorm)<<" ResPerDOF/Q="<< sqrt(resnorm/m_inArrayNorm*ototalDOF)
                                <<" ResMax/Q="<< resmaxm*sqrt(ntotalDOF/m_inArrayNorm)<<" ResPerDOF/(DtRHS): "<<sqrt(resratio*ototalDOF)
                                <<" with "<<setw(3)<<k<<" Non-Its"<<" and "<<setw(4)<<NtotDoOdeRHS<<" Lin-Its)"<<endl;
                        }
                        converged = true;
                        break;
                    // }
                    }
                }
            }

            //TODO: currently  NonlinSysRes is 2D array and SolveLinearSystem needs 1D array
            // LinSysTol = sqrt(0.01*sqrt(ratioTol)*resnorm);
            LinSysTol = 1.0E-5*sqrt(resnorm);
            // LinSysTol = 0.005*sqrt(resnorm)*(k+1);
            NtotDoOdeRHS  +=   linsol.SolveLinearSystem(ntotal,NonlinSysRes_1D,dsol_1D,0,LinSysTol);
            // cout << "NtotDoOdeRHS    = "<<NtotDoOdeRHS<<endl;

            // LinSysEPS  =   linsol.SolveLinearSystem(ntotal,NonlinSysRes_1D,dsol_1D,0);
            for(int i = 0; i < nvariables; i++)
            {
                Vmath::Vsub(npoints,m_TimeIntegtSol_k[i],1,dsol[i],1,m_TimeIntegtSol_k[i],1);
            }

        }

        m_TotLinItePerStep += NtotDoOdeRHS;

        /// TODO: disconnect these from other arrays to avoid memory cannot release.
        // m_TimeIntegtSol_k   =   nullptr;
        // m_TimeIntegtSol_n   =   nullptr;
        // m_SysEquatResid_k   =   nullptr;
        ASSERTL0((converged),"Nonlinear system solver not converge in CompressibleFlowSystem::DoImplicitSolve ");
    }

    void CompressibleFlowSystem::AllocatePrecondBlkDiag_coeff(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray)
    {

        int nvars = m_fields.num_elements();
        int nelmts  = m_fields[0]->GetNumElmts();
        int nrowsVars,ncolsVars;
        int nelmtcoef;
        DNekMatSharedPtr loc_matNvar;
        Array<OneD, unsigned int > nelmtmatdim(nelmts);
        for(int i = 0; i < nelmts; i++)
        {
            nelmtcoef   =   m_fields[0]->GetExp(i)->GetNcoeffs();
            nelmtmatdim[i]  =   nelmtcoef;
        }

        for(int i = 0; i < nvars; i++)
        {
            for(int j = 0; j < nvars; j++)
            {
                gmtxarray[i][j] = MemoryManager<DNekBlkMat>
                    ::AllocateSharedPtr(nelmtmatdim, nelmtmatdim, eDIAGONAL);
                for(int nelm = 0; nelm < nelmts; ++nelm)
                {
                    nrowsVars = nelmtmatdim[nelm];
                    ncolsVars = nrowsVars;
                    
                    loc_matNvar = MemoryManager<DNekMat>::AllocateSharedPtr(nrowsVars,ncolsVars,0.0);
                    gmtxarray[i][j]->SetBlock(nelm,nelm,loc_matNvar);
                }
                
            }
        }
    }

    void CompressibleFlowSystem::GetpreconditionerNSBlkDiag_coeff(
                                            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                                            Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,
                                            Array<OneD, DNekBlkMatSharedPtr > &TraceJac)
    {
        Fill2DArrayOfBlkDiagonalMat(gmtxarray,0.0);
        // Cout2DArrayBlkMat(gmtxarray);
        AddMatNSBlkDiag_volume(inarray,gmtxarray);

        AddMatNSBlkDiag_boundary(inarray,gmtxarray,TraceJac);
            
        MultiplyElmtInvMass_PlusSource(gmtxarray,m_TimeIntegLambda);
        ElmtVarInvMtrx(gmtxarray);
        // if(m_session->GetComm()->GetRank()==0)
        // {
        //     Cout2DArrayBlkMat(gmtxarray);
        // }
        // int tmpcin = 0;
        // cin >> tmpcin;
    }

    void CompressibleFlowSystem::MultiplyElmtInvMass_PlusSource(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,const NekDouble dtlamda)
    {
        MultiRegions::ExpListSharedPtr explist = m_fields[0];
            std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
        int ntotElmt            = (*pexp).size();
        int nElmtCoef;
        int nConvectiveFields = m_fields.num_elements();

        NekDouble Negdtlamda    =   -dtlamda;
        DNekMatSharedPtr        tmpGmtx,ElmtMat;

        Array<OneD, NekDouble> coefarray(explist->GetNcoeffs(),0.0);

        Array<OneD, int > elmtcoef(ntotElmt);
        int nElmtCoef0 =    (*pexp)[0]->GetNcoeffs();
        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            nElmtCoef           = (*pexp)[nelmt]->GetNcoeffs();
            ASSERTL0(nElmtCoef==nElmtCoef0,"nElmtCoef==nElmtCoef0");
            elmtcoef[nelmt]     =   nElmtCoef;
        }


        for(int m = 0; m < nConvectiveFields; m++)
        {
            for(int n = 0; n < nConvectiveFields; n++)
            {
                for(int ncl = 0; ncl < nElmtCoef0; ncl++)
                {
                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtCoef           = elmtcoef[nelmt];
                        tmpGmtx =   gmtxarray[m][n]->GetBlock(nelmt,nelmt);

                        int n_offset    =   explist->GetCoeff_Offset(nelmt);
                        
                        for(int nrw = 0; nrw < nElmtCoef; nrw++)
                        {
                            coefarray[n_offset+nrw]   =   (*tmpGmtx)(nrw,ncl);
                        }
                    }
                    explist->MultiplyByElmtInvMass(coefarray, coefarray);

                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        nElmtCoef            = elmtcoef[nelmt];
                        tmpGmtx =   gmtxarray[m][n]->GetBlock(nelmt,nelmt);
                        int n_offset    =   explist->GetCoeff_Offset(nelmt);
                        for(int nrw = 0; nrw < nElmtCoef; nrw++)
                        {
                            (*tmpGmtx)(nrw,ncl)   = Negdtlamda*coefarray[n_offset+nrw];
                        }
                    }
                }
            }
        }

        for(int m = 0; m < nConvectiveFields; m++)
        {
            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {

                NekDouble fac         = m_timestep/m_locTimeStep[nelmt];
                tmpGmtx =   gmtxarray[m][m]->GetBlock(nelmt,nelmt);
                nElmtCoef            = elmtcoef[nelmt];
                for(int ncl = 0; ncl < nElmtCoef; ncl++)
                {
                    (*tmpGmtx)(ncl,ncl)   += (1.0+fac);
                }
            }
        }

        return;
    }

    void CompressibleFlowSystem::DebugNumCalJac_coeff(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray)
    {
        MultiRegions::ExpListSharedPtr explist = m_fields[0];
            std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
        int nElmts    = (*pexp).size();
        int nvariables= gmtxarray.num_elements();

        Array<OneD, Array<OneD, DNekMatSharedPtr> >  ElmtPrecMatVars(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            ElmtPrecMatVars[i]  =   Array<OneD, DNekMatSharedPtr>(nvariables);
        }
        for(int i = 0; i < nElmts; i++)
        {
            DebugNumCalElmtJac_coeff(ElmtPrecMatVars,i);
            for(int m = 0; m < nvariables; m++)
            {
                for(int n = 0; n < nvariables; n++)
                {
                    gmtxarray[m][n]->SetBlock(i,i,ElmtPrecMatVars[m][n]);
                }
            }
        }
    }

    void CompressibleFlowSystem::DebugNumCalElmtJac_coeff(Array<OneD, Array<OneD, DNekMatSharedPtr> > &ElmtPrecMatVars ,const int nelmt)
    {
        MultiRegions::ExpListSharedPtr explist = m_fields[0];
            std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
        int nElmtcoef    = (*pexp)[nelmt]->GetNcoeffs();
        int nElmtOffset = explist->GetCoeff_Offset(nelmt);
        
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        unsigned int ntotpnt    = nvariables*npoints;
        Array<OneD, NekDouble > tmpinn_1d(ntotpnt,0.0);
        Array<OneD, NekDouble > tmpout_1d(ntotpnt,0.0);
        Array<OneD,       Array<OneD, NekDouble> > tmpinn(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > tmpout(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            int noffset = i*npoints;
            tmpinn[i] = tmpinn_1d+noffset;
            tmpout[i] = tmpout_1d+noffset;
        }


        for(int i = 0; i < nvariables; i++)
        {
            for(int j = 0; j < nvariables; j++)
            {
                ElmtPrecMatVars[i][j] =  MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nElmtcoef, nElmtcoef, 0.0);
            }
        }

        DNekMatSharedPtr    tmpStdMat;
        for (int i = 0; i < nvariables; i++)
        {
            for (int npnt = 0; npnt < nElmtcoef; npnt++)
            {
                tmpinn[i][nElmtOffset+npnt] = 1.0;
                MatrixMultiply_MatrixFree_coeff(tmpinn_1d,tmpout_1d);
                
                for (int j = 0; j < nvariables; j++)
                {
                    for (int npntf = 0; npntf < nElmtcoef; npntf++)
                    {
                        tmpStdMat = ElmtPrecMatVars[j][i];
                        // NekDouble tmp = tmpout[j][nElmtOffset+npntf];
                        (*tmpStdMat)(npntf,npnt) = tmpout[j][nElmtOffset+npntf];
                    }
                }
                tmpinn[i][nElmtOffset+npnt] = 0.0;
            }
        }
    }
    
    void CompressibleFlowSystem::MatrixMultiply_MatrixFree_coeff(
                                                 const  Array<OneD, NekDouble> &inarray,
                                                        Array<OneD, NekDouble >&out)
    {
        NekDouble eps = 1.0E-6;
        NekDouble oeps = 1.0/eps;
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int npoints    = m_TimeIntegtSol_n[0].num_elements();
        Array<OneD, NekDouble > tmp;
        Array<OneD,       Array<OneD, NekDouble> > solplus(nvariables);
        Array<OneD,       Array<OneD, NekDouble> > resplus(nvariables);
        for(int i = 0; i < nvariables; i++)
        {
            solplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
            resplus[i] =  Array<OneD, NekDouble>(npoints,0.0);
        }

        for (int i = 0; i < nvariables; i++)
        {
            tmp = inarray + i*npoints;
            Vmath::Svtvp(npoints,eps,tmp,1,m_TimeIntegtSol_k[i],1,solplus[i],1);
        }
        
        NonlinSysEvaluator_coeff(solplus,resplus);

        for (int i = 0; i < nvariables; i++)
        {
            tmp = out + i*npoints;
            Vmath::Vsub(npoints,&resplus[i][0],1,&m_SysEquatResid_k[i][0],1,&tmp[0],1);
            Vmath::Smul(npoints, oeps ,&tmp[0],1,&tmp[0],1);
        }
       
        return;
    }

    void CompressibleFlowSystem::MatrixMultiply_MatrixFree_coeff_dualtimestep(
                                                 const  Array<OneD, NekDouble> &inarray,
                                                        Array<OneD, NekDouble >&out)
    {
        

        MatrixMultiply_MatrixFree_coeff(inarray,out);

        int nElements = m_fields[0]->GetExpSize();
        int nelmtcoeffs,offset,varoffset;
        NekDouble fac;
        Array<OneD, NekDouble> tmp;
        Array<OneD, NekDouble> tmpinn;
        Array<OneD, NekDouble> tmpout;
        unsigned int nvariables = m_TimeIntegtSol_n.num_elements();
        unsigned int ntotcoeffs = m_TimeIntegtSol_n[0].num_elements();

        
        for(int i = 0; i < nvariables; ++i)
        {
            varoffset = i*ntotcoeffs;

            tmpinn = inarray + varoffset;
            tmpout = out + varoffset;
            // Loop over elements
            for(int n = 0; n < nElements; ++n)
            {
                nelmtcoeffs = m_fields[0]->GetExp(n)->GetNcoeffs();
                offset      = m_fields[0]->GetCoeff_Offset(n);
                fac         = m_timestep/m_locTimeStep[n];
                tmp         = tmpout + offset;
                Vmath::Svtvp(nelmtcoeffs,fac,tmpinn+offset,1,tmp,1,tmp,1);
            }
        }
        return;
    }



    void CompressibleFlowSystem::NonlinSysEvaluator_coeff(
                                                       Array<OneD, Array<OneD, NekDouble> > &inarray,
                                                       Array<OneD, Array<OneD, NekDouble> > &out)
    {
        Array<OneD, Array<OneD, NekDouble> > sol_n;
        sol_n                  = m_TimeIntegtSol_n;
        //inforc = m_TimeIntegForce;
        unsigned int nvariable  = inarray.num_elements();
        unsigned int ncoeffs    = inarray[0].num_elements();
        unsigned int npoints    = m_fields[0]->GetNpoints();

        Array<OneD, Array<OneD, NekDouble> > inpnts(nvariable);

        for(int i = 0; i < nvariable; i++)
        {
            inpnts[i]   =   Array<OneD, NekDouble>(npoints,0.0);
            m_fields[i]->BwdTrans(inarray[i], inpnts[i]);
        }

        DoOdeProjection(inpnts,inpnts,m_BndEvaluateTime);
        DoOdeRhs_coeff(inpnts,out,m_BndEvaluateTime);

        for (int i = 0; i < nvariable; i++)
        {
            Vmath::Svtvp(ncoeffs,m_TimeIntegLambda,out[i],1,sol_n[i],1,out[i],1);
            Vmath::Vsub(ncoeffs,inarray[i],1,out[i],1,out[i],1);
        }
        return;
    }


    void CompressibleFlowSystem::NonlinSysEvaluator_coeff_out(
                                                       Array<OneD, Array<OneD, NekDouble> > &inarray,
                                                       Array<OneD, Array<OneD, NekDouble> > &out)
    {
        Array<OneD, Array<OneD, NekDouble> > sol_n;
        sol_n                  = m_TimeIntegtSol_n;
        //inforc = m_TimeIntegForce;
        unsigned int nvariable  = inarray.num_elements();
        unsigned int ncoeffs    = inarray[0].num_elements();
        unsigned int npoints    = m_fields[0]->GetNpoints();

        Array<OneD, Array<OneD, NekDouble> > inpnts(nvariable);

        for(int i = 0; i < nvariable; i++)
        {
            inpnts[i]   =   Array<OneD, NekDouble>(npoints,0.0);
            m_fields[i]->BwdTrans(inarray[i], inpnts[i]);
        }

        DoOdeProjection(inpnts,inpnts,m_BndEvaluateTime);
        DoOdeRhs_coeff(inpnts,out,m_BndEvaluateTime);


        LibUtilities::CommSharedPtr v_Comm  = m_fields[0]->GetComm()->GetRowComm();
        NekDouble resnorm = 0.0;
        for(int i = 0; i < nvariable; i++)
        {
            resnorm += Vmath::Dot(ncoeffs,out[i],out[i]);
        }
        v_Comm->AllReduce(resnorm, Nektar::LibUtilities::ReduceSum);
        if(0==m_session->GetComm()->GetRank())
        {
            cout << "        m_TimeIntegLambda= "<<m_TimeIntegLambda<<" DtRhs^2= "<<m_TimeIntegLambda*m_TimeIntegLambda*resnorm<<endl;
        }

        for (int i = 0; i < nvariable; i++)
        {
            Vmath::Svtvp(ncoeffs,m_TimeIntegLambda,out[i],1,sol_n[i],1,out[i],1);
            Vmath::Vsub(ncoeffs,inarray[i],1,out[i],1,out[i],1);
        }
        return;
    }


    /**
     * @brief Compute the right-hand side.
     */
    void CompressibleFlowSystem::DoOdeRhs_coeff(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time)
    {
        int i;
        int nvariables = inarray.num_elements();
        int nTracePts  = GetTraceTotPoints();
        int ncoeffs    = GetNcoeffs();

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
            for(i = 0; i < nvariables; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePts, 0.0);
                m_fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }


        // Calculate advection
        DoAdvection_coeff(inarray, outarray, time, Fwd, Bwd);
        // Negate results
        
        for (i = 0; i < nvariables; ++i)
        {
            Vmath::Neg(ncoeffs, outarray[i], 1);
        }

        // Add diffusion t
        // DoDiffusion(inarray, outarray, Fwd, Bwd);
        DoDiffusion_coeff(inarray, outarray, Fwd, Bwd);

        // Add forcing terms
        // for (auto &x : m_forcing)
        // {
        //     x->Apply(m_fields, inarray, outarray, time);
        // }

        // if (m_useLocalTimeStep)
        // {
        //     int nElements = m_fields[0]->GetExpSize();
        //     int nq, offset;
        //     NekDouble fac;
        //     Array<OneD, NekDouble> tmp;

        //     Array<OneD, NekDouble> tstep (nElements, 0.0);
        //     GetElmtTimeStep(inarray, tstep);

        //     // Loop over elements
        //     for(int n = 0; n < nElements; ++n)
        //     {
        //         nq     = m_fields[0]->GetExp(n)->GetTotPoints();
        //         offset = m_fields[0]->GetPhys_Offset(n);
        //         fac    = tstep[n] / m_timestep;
        //         for(i = 0; i < nvariables; ++i)
        //         {
        //             Vmath::Smul(nq, fac, outarray[i] + offset, 1,
        //                                  tmp = outarray[i] + offset, 1);
        //         }
        //     }
        // }
    }

    /**
     * @brief Compute the advection terms for the right-hand side
     */
    void CompressibleFlowSystem::DoAdvection_coeff(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
        const NekDouble                                   time,
        const Array<OneD, Array<OneD, NekDouble> >       &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >       &pBwd)
    {
        int nvariables = inarray.num_elements();
        Array<OneD, Array<OneD, NekDouble> > advVel(m_spacedim);

        m_advObject->Advect_coeff(nvariables, m_fields, advVel, inarray,
                            outarray, time, pFwd, pBwd);
    }


#endif
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
        int nvariables = inarray.num_elements();
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

        if (m_shockCaptureType != "Off")
        {
            m_artificialDiffusion->DoArtificialDiffusion(inarray, outarray);
        }
    }


    /**
     * @brief Add the diffusions terms to the right-hand side
     */
    void CompressibleFlowSystem::DoDiffusion_coeff(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
              Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const Array<OneD, Array<OneD, NekDouble> >   &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >   &pBwd)
    {
        v_DoDiffusion_coeff(inarray, outarray, pFwd, pBwd);

        if (m_shockCaptureType != "Off")
        {
            // m_artificialDiffusion->DoArtificialDiffusion(inarray, outarray);
            m_artificialDiffusion->DoArtificialDiffusion_coeff(inarray, outarray);
        }
    }

    void CompressibleFlowSystem::SetBoundaryConditions(
            Array<OneD, Array<OneD, NekDouble> >             &physarray,
            NekDouble                                         time)
    {
        int nTracePts  = GetTraceTotPoints();
        int nvariables = physarray.num_elements();

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
     * @brief Return the flux vector for the compressible Euler equations.
     *
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        int i, j;
        int nq = physfield[0].num_elements();
        int nVariables = m_fields.num_elements();

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

    Array<OneD, Array<OneD, DNekMatSharedPtr> > CompressibleFlowSystem::GetFluxVectorJacDirctn(const int nDirctn,
                                                        const Array<OneD, const Array<OneD, NekDouble> >&inarray)
    {
        int nConvectiveFields   = inarray.num_elements();
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    m_fields[0]->GetExp();
        int ntotElmt            = (*expvect).size();
        LocalRegions::ExpansionSharedPtr pexp = (*expvect)[0];
        // int nElmtPnt            = pexp->GetTotPoints();
        // int nElmtCoef           = pexp->GetNcoeffs();
        // int nSpaceDim           = m_graph->GetSpaceDimension();  
        Array<OneD, NekDouble> normals;
        Array<OneD, Array<OneD, NekDouble> > normal3D(3);
        for(int i = 0; i < 3; i++)
        {
            normal3D[i] = Array<OneD, NekDouble>(3,0.0);
        }
        normal3D[0][0] = 1.0;
        normal3D[1][1] = 1.0;
        normal3D[2][2] = 1.0;
        normals =   normal3D[nDirctn];

        Array<OneD, NekDouble> pointVar(nConvectiveFields,0.0);

        Array<OneD, Array<OneD, DNekMatSharedPtr> > ElmtJac(ntotElmt);

        Array<OneD, Array<OneD, NekDouble> > locVars(nConvectiveFields);


        for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
        {
            // nElmtCoef           = (*expvect)[nelmt]->GetNcoeffs();
            int nElmtPnt            = (*expvect)[nelmt]->GetTotPoints();
    
            for(int j = 0; j < nConvectiveFields; j++)
            {   
                locVars[j] = inarray[j]+GetPhys_Offset(nelmt);
            }

            ElmtJac[nelmt] =   Array<OneD, DNekMatSharedPtr>(nElmtPnt);
            for(int npnt = 0; npnt < nElmtPnt; npnt++)
            {
                ElmtJac[nelmt][npnt] = MemoryManager<DNekMat>
                    ::AllocateSharedPtr(nConvectiveFields, nConvectiveFields);
                for(int j = 0; j < nConvectiveFields; j++)
                {
                    pointVar[j] = locVars[j][npnt];
                }

                GetFluxVectorJacPoint(pointVar,normals,ElmtJac[nelmt][npnt]);
            }
        }
        return ElmtJac;
    }


    void CompressibleFlowSystem::GetFluxVectorJacPoint(
            const Array<OneD, NekDouble>                &conservVar, 
            const Array<OneD, NekDouble>                &normals, 
                 DNekMatSharedPtr                       &fluxJac)
    {
        int nvariables      = conservVar.num_elements();
        int nvariables3D    = 5;
        int expDim          = m_spacedim;

        if (nvariables > expDim+2)
        {
            ASSERTL0(false,"nvariables > expDim+2 case not coded")
        }

        DNekMatSharedPtr PointFJac3D = MemoryManager<DNekMat>
            ::AllocateSharedPtr(nvariables3D, nvariables3D);
        
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
        NekDouble fsw,efix_StegerWarming;
        efix_StegerWarming = 0.0;

        fsw = 0.0; // exact flux Jacobian if fsw=0.0
        PointFluxJacobian_pn(PointFwd,normals,PointFJac3D,efix_StegerWarming,fsw);

        for(int j=0; j< nvariables; j++)
        {
            nj = index[j];
            for(int k=0; k< nvariables; k++)
            {
                nk = index[k];
                (*fluxJac)(j,k) = (*PointFJac3D)(nj,nk); 
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
            gama    = m_varConv->Geteos()->GetGamma();

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

            int nsf = 0;
            (*FJac)(nsf  ,nsf  ) = c1*y1 - d1*vna + l1;
            (*FJac)(nsf  ,nsf+1) = -c1*vx + d1*nxa;
            (*FJac)(nsf  ,nsf+2) = -c1*vy + d1*nya;
            (*FJac)(nsf  ,nsf+3) = -c1*vz + d1*nza;
            (*FJac)(nsf  ,nsf+4) = c1;
            c2 = c1*vx + d1*nxa*ae;
            d2 = x3*nxa + d1*vx;
            (*FJac)(nsf+1,nsf  ) = c2*y1 - d2*vna;
            (*FJac)(nsf+1,nsf+1) = -c2*vx + d2*nxa + l1;
            (*FJac)(nsf+1,nsf+2) = -c2*vy + d2*nya;
            (*FJac)(nsf+1,nsf+3) = -c2*vz + d2*nza;
            (*FJac)(nsf+1,nsf+4) = c2;
            c3 = c1*vy + d1*nya*ae;
            d3 = x3*nya + d1*vy;
            (*FJac)(nsf+2,nsf  ) = c3*y1 - d3*vna;
            (*FJac)(nsf+2,nsf+1) = -c3*vx + d3*nxa;
            (*FJac)(nsf+2,nsf+2) = -c3*vy + d3*nya + l1;
            (*FJac)(nsf+2,nsf+3) = -c3*vz + d3*nza;
            (*FJac)(nsf+2,nsf+4) = c3;
            c4 = c1*vz + d1*nza*ae;
            d4 = x3*nza + d1*vz;
            (*FJac)(nsf+3,nsf  ) = c4*y1 - d4*vna;
            (*FJac)(nsf+3,nsf+1) = -c4*vx + d4*nxa;
            (*FJac)(nsf+3,nsf+2) = -c4*vy + d4*nya;
            (*FJac)(nsf+3,nsf+3) = -c4*vz + d4*nza + l1;
            (*FJac)(nsf+3,nsf+4) = c4;
            c5 = c1*h0 + d1*vna*ae;
            d5 = x3*vna + d1*h0;
            (*FJac)(nsf+4,nsf  ) = c5*y1 - d5*vna;
            (*FJac)(nsf+4,nsf+1) = -c5*vx + d5*nxa;
            (*FJac)(nsf+4,nsf+2) = -c5*vy + d5*nya;
            (*FJac)(nsf+4,nsf+3) = -c5*vz + d5*nza;
            (*FJac)(nsf+4,nsf+4) = c5 + l1;

    }

    /**
     * @brief Return the flux vector for the compressible Euler equations
     * by using the de-aliasing technique.
     *
     * @param physfield   Fields.
     * @param flux        Resulting flux.
     */
    void CompressibleFlowSystem::GetFluxVectorDeAlias(
        const Array<OneD, Array<OneD, NekDouble> >               &physfield,
              Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &flux)
    {
        int i, j;
        int nq = physfield[0].num_elements();
        int nVariables = m_fields.num_elements();

        // Factor to rescale 1d points in dealiasing
        NekDouble OneDptscale = 2;
        nq = m_fields[0]->Get1DScaledTotPoints(OneDptscale);

        Array<OneD, NekDouble> pressure(nq);
        Array<OneD, Array<OneD, NekDouble> > velocity(m_spacedim);

        Array<OneD, Array<OneD, NekDouble> > physfield_interp(nVariables);
        Array<OneD, Array<OneD, Array<OneD, NekDouble> > > flux_interp(
                                                            nVariables);

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

        // Galerkin project solution back to origianl space
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

            // Galerkin project solution back to origianl space
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
        int n;
        int nElements = m_fields[0]->GetExpSize();

        // Change value of m_timestep (in case it is set to zero)
        NekDouble tmp = m_timestep;
        m_timestep    = 1.0;

        Array<OneD, NekDouble> cfl(nElements);
        cfl = GetElmtCFLVals();

        // Factors to compute the time-step limit
        NekDouble alpha     = MaxTimeStepEstimator();

        // Loop over elements to compute the time-step limit for each element
        for(n = 0; n < nElements; ++n)
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
        EquationSystem::v_SetInitialConditions(initialtime, false);

        // insert white noise in initial condition
        NekDouble Noise;
        int phystot = m_fields[0]->GetTotPoints();
        Array<OneD, NekDouble> noise(phystot);

        m_session->LoadParameter("Noise", Noise,0.0);
        int m_nConvectiveFields =  m_fields.num_elements();

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
    Array<OneD, NekDouble> CompressibleFlowSystem::v_GetMaxStdVelocity()
    {
        int nTotQuadPoints = GetTotPoints();
        int n_element      = m_fields[0]->GetExpSize();
        int expdim         = m_fields[0]->GetGraph()->GetMeshDimension();
        int nfields        = m_fields.num_elements();
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

        for(int el = 0; el < n_element; ++el)
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
            if(metricInfo->GetGtype() == SpatialDomains::eDeformed)
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
            Array<OneD, Array<OneD, NekDouble> > tmp(m_fields.num_elements());

            for (int i = 0; i < m_fields.num_elements(); ++i)
            {
                tmp[i] = m_fields[i]->GetPhys();
            }

            Array<OneD, NekDouble> pressure(nPhys), temperature(nPhys);
            Array<OneD, NekDouble> entropy(nPhys);
            Array<OneD, NekDouble> soundspeed(nPhys), mach(nPhys);
            Array<OneD, NekDouble> sensor(nPhys), SensorKappa(nPhys);

            m_varConv->GetPressure  (tmp, pressure);
            m_varConv->GetTemperature(tmp, temperature);
            m_varConv->GetEntropy   (tmp, entropy);
            m_varConv->GetSoundSpeed(tmp, soundspeed);
            m_varConv->GetMach      (tmp, soundspeed, mach);

            int sensorOffset;
            m_session->LoadParameter ("SensorOffset", sensorOffset, 1);
            m_varConv->GetSensor (m_fields[0], tmp, sensor, SensorKappa,
                                    sensorOffset);

            Array<OneD, NekDouble> pFwd(nCoeffs), TFwd(nCoeffs);
            Array<OneD, NekDouble> sFwd(nCoeffs);
            Array<OneD, NekDouble> aFwd(nCoeffs), mFwd(nCoeffs);
            Array<OneD, NekDouble> sensFwd(nCoeffs);

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

            if (m_artificialDiffusion)
            {
                Array<OneD, NekDouble> sensorFwd(nCoeffs);
                // reuse pressure
                m_artificialDiffusion->GetArtificialViscosity(tmp, pressure);
                m_fields[0]->FwdTrans_IterPerExp(pressure,   sensorFwd);

                variables.push_back  ("ArtificialVisc");
                fieldcoeffs.push_back(sensorFwd);
            }
        }
    }


    void CompressibleFlowSystem::Cout1DArrayBlkMat(Array<OneD, DNekBlkMatSharedPtr> &gmtxarray,const unsigned int nwidthcolm)
    {
        int nvar1 = gmtxarray.num_elements();

        
        for(int i = 0; i < nvar1; i++)
        {
            cout<<endl<<"£$£$£$£$£$£$££$£$£$$$£$££$$£$££$£$$££££$$£$£$£$£$£$£$££$£$$"<<endl<< "Cout2DArrayBlkMat i= "<<i<<endl;
            CoutBlkMat(gmtxarray[i],nwidthcolm);
        }
    }
    
    void CompressibleFlowSystem::Cout2DArrayBlkMat(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,const unsigned int nwidthcolm)
    {
        int nvar1 = gmtxarray.num_elements();
        int nvar2 = gmtxarray[0].num_elements();

        
        for(int i = 0; i < nvar1; i++)
        {
            for(int j = 0; j < nvar2; j++)
            {
                cout<<endl<<"£$£$£$£$£$£$££$£$£$$$£$££$$£$££$£$$££££$$£$£$£$£$£$£$££$£$$"<<endl<< "Cout2DArrayBlkMat i= "<<i<<" j=  "<<j<<endl;
                CoutBlkMat(gmtxarray[i][j],nwidthcolm);
            }
        }
    }

    void CompressibleFlowSystem::CoutBlkMat(DNekBlkMatSharedPtr &gmtx,const unsigned int nwidthcolm)
    {
        DNekMatSharedPtr    loc_matNvar;

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;
        gmtx->GetBlockSizes(rowSizes,colSizes);

        int nelmts  = rowSizes.num_elements();
        
        // int noffset = 0;
        for(int i = 0; i < nelmts; ++i)
        {
            loc_matNvar =   gmtx->GetBlock(i,i);
            std::cout   <<std::endl<<"*********************************"<<std::endl<<"element :   "<<i<<std::endl;
            CoutStandardMat(loc_matNvar,nwidthcolm);
        }
        return;
    }


    void CompressibleFlowSystem::CoutStandardMat(DNekMatSharedPtr &loc_matNvar,const unsigned int nwidthcolm)
    {
        int nrows = loc_matNvar->GetRows();
        int ncols = loc_matNvar->GetColumns();
        NekDouble tmp=0.0;
        std::cout   <<"ROW="<<std::setw(3)<<-1<<" ";
        for(int k = 0; k < ncols; k++)
        {
            std::cout   <<"   COL="<<std::setw(nwidthcolm-7)<<k;
        }
        std::cout   << endl;

        for(int j = 0; j < nrows; j++)
        {
            std::cout   <<"ROW="<<std::setw(3)<<j<<" ";
            for(int k = 0; k < ncols; k++)
            {
                tmp =   (*loc_matNvar)(j,k);
                std::cout   <<std::scientific<<std::setw(nwidthcolm)<<std::setprecision(nwidthcolm-8)<<tmp;
            }
            std::cout   << endl;
        }
    }

    void CompressibleFlowSystem::Fill2DArrayOfBlkDiagonalMat(Array<OneD, Array<OneD, DNekBlkMatSharedPtr> > &gmtxarray,const NekDouble valu)
    {
        
        int n1d = gmtxarray.num_elements();
        int n2d = gmtxarray[0].num_elements();

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;

        DNekMatSharedPtr    loc_matNvar;

        for(int n1 = 0; n1 < n1d; ++n1)
        {
            for(int n2 = 0; n2 < n2d; ++n2)
            {
                
                gmtxarray[n1][n2]->GetBlockSizes(rowSizes,colSizes);
                int nelmts  = rowSizes.num_elements();

                for(int i = 0; i < nelmts; ++i)
                {
                    loc_matNvar =   gmtxarray[n1][n2]->GetBlock(i,i);

                    int nrows = loc_matNvar->GetRows();
                    int ncols = loc_matNvar->GetColumns();

                    for(int j = 0; j < nrows; j++)
                    {
                        for(int k = 0; k < ncols; k++)
                        {
                            (*loc_matNvar)(j,k)=valu;
                        }
                    }
                }
            }
        }

    }

}
