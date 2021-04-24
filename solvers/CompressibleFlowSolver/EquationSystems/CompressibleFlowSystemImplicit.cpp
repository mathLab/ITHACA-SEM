///////////////////////////////////////////////////////////////////////////////
//
// File CFSImplicit.cpp
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

#include <CompressibleFlowSolver/EquationSystems/CompressibleFlowSystemImplicit.h>

#include <LibUtilities/BasicUtils/Timer.h>

using namespace std;

namespace Nektar
{
    CFSImplicit::CFSImplicit(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
        : UnsteadySystem(pSession, pGraph),
          CompressibleFlowSystem(pSession, pGraph)
    {
    }

    /**
     * @brief Initialization object for CFSImplicit class.
     */
    void CFSImplicit::v_InitObject()
    {
        CompressibleFlowSystem::v_InitObject();

        // initialise implicit parameters
        m_session->LoadParameter ("JacobiFreeEps", m_jacobiFreeEps, 5.0E-8);

        int ntmp;
        m_session->LoadParameter("AdvectionJacFlag", ntmp, 1);
        m_advectionJacFlag = (0==ntmp)? false:true;

        m_session->LoadParameter("ViscousJacFlag", ntmp, 1);
        m_viscousJacFlag   = (0==ntmp)? false:true;

        // initialise implicit functors
        m_ode.DefineImplicitSolve(&CFSImplicit::DoImplicitSolve, this);
        
        InitialiseNonlinSysSolver();

        m_flagImplicitSolver = true;
    }
    
    void CFSImplicit::InitialiseNonlinSysSolver()
    {
        std::string SolverType = "Newton";
        if (m_session->DefinesSolverInfo("NonlinSysIterSolver"))
        {
            SolverType = m_session->GetSolverInfo("NonlinSysIterSolver");
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

        LibUtilities::NekSysOperators       nekSysOp;
        nekSysOp.DefineNekSysResEval(&CFSImplicit::
            NonlinSysEvaluatorCoeff1D, this);
        nekSysOp.DefineNekSysLhsEval(&CFSImplicit::
            MatrixMultiplyMatrixFreeCoeff, this);
        nekSysOp.DefineNekSysPrecon(
                    &CFSImplicit::PreconCoeff, this);

        int nvariables  =   m_fields.size();
        const MultiRegions::LocTraceToTraceMapSharedPtr locTraceToTraceMap =
            m_fields[0]->GetLocTraceToTraceMap();

        locTraceToTraceMap->CalcLocTracePhysToTraceIDMap(
            m_fields[0]->GetTrace(),m_spacedim);
        for (int i = 1; i < nvariables; i++)
        {
            m_fields[i]->GetLocTraceToTraceMap()
                ->SetLocTracePhysToTraceIDMap(
                locTraceToTraceMap->GetLocTracephysToTraceIDMap());
        }
        
        m_preconCfs = GetPreconCfsOpFactory().CreateInstance(
            "PreconCfsBRJ", m_fields, m_session, m_comm);
        NekPreconCfsOperators tmpPreconOp;
        tmpPreconOp.DefineCalcPreconMatBRJCoeff(&CFSImplicit::
            CalcPreconMatBRJCoeff, this);
        m_preconCfs->SetOperators(tmpPreconOp);

        m_session->LoadParameter("NewtonAbsoluteIteTol", 
            m_newtonAbsoluteIteTol, 1.0E-12);

        m_session->LoadParameter("nPadding", m_nPadding, 4);

        m_nonlinsol->SetSysOperators(nekSysOp);
    }

    /**
     * @brief Destructor for CFSImplicit class.
     */
    CFSImplicit::~CFSImplicit()
    {
    }

    void CFSImplicit::NonlinSysEvaluatorCoeff1D(
        const Array<OneD, const NekDouble>  &inarray,
              Array<OneD,       NekDouble>  &out,
        const bool                          &flag)
    {
        const Array<OneD, const NekDouble> refsource 
            = m_nonlinsol->GetRefSourceVec();
        NonlinSysEvaluatorCoeff(inarray, out, flag, refsource);
    }

    void CFSImplicit::NonlinSysEvaluatorCoeff(
        const Array<OneD, const NekDouble>  &inarray,
              Array<OneD,       NekDouble>  &out,
        const bool                          &flag,
        const Array<OneD, const NekDouble>  &source)
    {
        boost::ignore_unused(flag);
        unsigned int nvariables     = m_fields.size();
        unsigned int npoints        = m_fields[0]->GetNcoeffs();
        Array<OneD, Array<OneD, NekDouble>> in2D(nvariables);
        Array<OneD, Array<OneD, NekDouble>> out2D(nvariables);
        Array<OneD, Array<OneD, NekDouble>> source2D(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            int offset = i * npoints;
            in2D[i]    = inarray + offset;
            out2D[i]   = out + offset;
            source2D[i]   = source + offset;
        }
        NonlinSysEvaluatorCoeff(in2D, out2D, source2D);
    }

    void CFSImplicit::NonlinSysEvaluatorCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              Array<OneD,       Array<OneD, NekDouble>> &out,
        const Array<OneD, const Array<OneD, NekDouble>> &source)
    {
        unsigned int nvariable  = inarray.size();
        unsigned int ncoeffs    = inarray[nvariable - 1].size();
        unsigned int npoints    = m_fields[0]->GetNpoints();

        Array<OneD, Array<OneD, NekDouble>> inpnts(nvariable);

        for (int i = 0; i < nvariable; ++i)
        {
            inpnts[i] = Array<OneD, NekDouble>(npoints, 0.0);
            m_fields[i]->BwdTrans(inarray[i], inpnts[i]);
        }
        
        DoOdeProjection(inpnts, inpnts, m_bndEvaluateTime);
        DoOdeRhsCoeff(inpnts, out, m_bndEvaluateTime);

        for (int i = 0; i < nvariable; ++i)
        {
            Vmath::Svtvp(ncoeffs, -m_TimeIntegLambda, out[i], 1,
                         inarray[i], 1, out[i], 1);
        }

        if (NullNekDoubleArrayOfArray != source)
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
     * @brief Compute the right-hand side.
     */
    void CFSImplicit::DoOdeRhsCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              Array<OneD,       Array<OneD, NekDouble>> &outarray,
        const NekDouble                                 time)
    {
        int nvariables = inarray.size();
        int nTracePts  = GetTraceTotPoints();
        int ncoeffs    = GetNcoeffs();
        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble>> Fwd    (nvariables);
        Array<OneD, Array<OneD, NekDouble>> Bwd    (nvariables);

        if (m_HomogeneousType == eHomogeneous1D)
        {
            Fwd = NullNekDoubleArrayOfArray;
            Bwd = NullNekDoubleArrayOfArray;
        }
        else
        {
            for (int i = 0; i < nvariables; ++i)
            {
                Fwd[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
                Bwd[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
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
            Array<OneD, NekDouble> tmp;

            Array<OneD, NekDouble> tstep(nElements, 0.0);
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

    /**
     * @brief Compute the advection terms for the right-hand side
     */
    void CFSImplicit::DoAdvectionCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              Array<OneD,       Array<OneD, NekDouble>> &outarray,
        const NekDouble                                 time,
        const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &pBwd)
    {
        int nvariables = inarray.size();
        Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);

        m_advObject->AdvectCoeffs(nvariables, m_fields, advVel, inarray,
                                  outarray, time, pFwd, pBwd);
    }

    void CFSImplicit::DoImplicitSolve(
        const Array<OneD, const Array<OneD, NekDouble>> &inpnts,
              Array<OneD,       Array<OneD, NekDouble>> &outpnt,
        const NekDouble                                 time,
        const NekDouble                                 lambda)
    {
        unsigned int nvariables  = inpnts.size();
        unsigned int ncoeffs     = m_fields[0]->GetNcoeffs();
        unsigned int ntotal      = nvariables * ncoeffs;

        Array<OneD, NekDouble>  inarray(ntotal);
        Array<OneD, NekDouble>  out(ntotal);
        Array<OneD, NekDouble>  tmpArray;

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

    void CFSImplicit::DoImplicitSolveCoeff(
            const Array<OneD, const Array<OneD, NekDouble>> &inpnts,
            const Array<OneD, const NekDouble>              &inarray,
                  Array<OneD,       NekDouble>              &out,
            const NekDouble                                 time,
            const NekDouble                                 lambda)
    {
        boost::ignore_unused(inpnts);

        m_TimeIntegLambda   = lambda;
        m_bndEvaluateTime   = time;
        m_solutionPhys      = inpnts;
        unsigned int ntotal = inarray.size();

        if (m_inArrayNorm < 0.0)
        {
            CalcRefValues(inarray);
        }
        
        NekDouble tol2 = m_inArrayNorm
                        *m_newtonAbsoluteIteTol * m_newtonAbsoluteIteTol;

        m_nonlinsol->v_SetupNekNonlinSystem(ntotal, inarray, inarray, 0);

        m_TotNewtonIts +=  m_nonlinsol->SolveSystem(ntotal,inarray,
                                                    out, 0, tol2);

        m_TotLinIts += m_nonlinsol->GetNtotLinSysIts();

        m_TotImpStages++;
        m_StagesPerStep++;
    }

    void CFSImplicit::CalcRefValues(
            const Array<OneD, const NekDouble> &inarray)
    {
        unsigned int nvariables         = m_fields.size();
        unsigned int ntotal             = inarray.size();
        unsigned int npoints            = ntotal/nvariables;

        unsigned int nTotalGlobal       = ntotal;
        m_comm->AllReduce(nTotalGlobal, Nektar::LibUtilities::ReduceSum);
        unsigned int nTotalDOF          = nTotalGlobal / nvariables;
        NekDouble invTotalDOF             = 1.0 / nTotalDOF;

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
            m_magnitdEstimat[i] = sqrt(m_magnitdEstimat[i] * invTotalDOF);
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

    void CFSImplicit::PreconCoeff(
            const Array<OneD, NekDouble> &inarray,
                  Array<OneD, NekDouble> &outarray,
            const bool                   &flag)
    {
        if (m_preconCfs->UpdatePreconMatCheck(NullNekDouble1DArray, 
            m_TimeIntegLambda) && m_flagUpdatePreconMat)
        {
            int nvariables = m_solutionPhys.size();

            int nphspnt = m_solutionPhys[nvariables - 1].size();
            Array<OneD, Array<OneD, NekDouble>> intmp(nvariables);
            for(int i = 0; i < nvariables; i++)
            {
                intmp[i]    =   Array<OneD, NekDouble>(nphspnt,0.0);
            }

            DoOdeProjection(m_solutionPhys,intmp,m_bndEvaluateTime);

            m_preconCfs->BuildPreconCfs(m_fields, intmp, m_bndEvaluateTime,
                m_TimeIntegLambda);
        }

        m_flagUpdatePreconMat = false;

        m_preconCfs->DoPreconCfs(m_fields, inarray, outarray, flag);
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CFSImplicit::AddMatNSBlkDiagVol(
        const Array<OneD, const Array<OneD, NekDouble>>         &inarray,
        const Array<OneD, const TensorOfArray2D<NekDouble>>     &qfield,
              Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>>  &gmtxarray,
              TensorOfArray4D<DataType>                         &StdMatDataDBB,
              TensorOfArray5D<DataType>                         &StdMatDataDBDB)
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
        int nTotElmt            = (*expvect).size();

        Array<OneD, NekDouble > mu      (npoints, 0.0);
        Array<OneD, NekDouble > DmuDT   (npoints, 0.0);
        if(m_viscousJacFlag)
        {
            CalcMuDmuDT(inarray,mu,DmuDT);
        }

        Array<OneD, NekDouble> normals;
        Array<OneD, Array<OneD, NekDouble>> normal3D(3);
        for(int i = 0; i < 3; i++)
        {
            normal3D[i] = Array<OneD, NekDouble>(3,0.0);
        }
        normal3D[0][0] = 1.0;
        normal3D[1][1] = 1.0;
        normal3D[2][2] = 1.0;
        Array<OneD, Array<OneD, NekDouble>> normalPnt(3);
        
        DNekMatSharedPtr wspMat     = 
            MemoryManager<DNekMat>::AllocateSharedPtr(nvariable,nvariable,0.0);
        DNekMatSharedPtr wspMatDrv  = 
            MemoryManager<DNekMat>::AllocateSharedPtr(nvariable-1,nvariable,0.0);

        Array<OneD, DataType> GmatxData;
        Array<OneD, DataType> MatData;

        Array<OneD, NekDouble> tmppnts;
        TensorOfArray3D<NekDouble> PntJacCons(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray3D<DataType>  PntJacConsStd(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        Array<OneD, Array<OneD, NekDouble>> ConsStdd(m_spacedim);
        Array<OneD, Array<OneD, NekDouble>> ConsCurv(m_spacedim);
        TensorOfArray4D<NekDouble> PntJacDerv(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray4D<DataType>  PntJacDervStd(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray3D<NekDouble> DervStdd(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        TensorOfArray3D<NekDouble> DervCurv(m_spacedim); // Nvar*Nvar*Ndir*Nelmt*Npnt
        for(int ndir=0; ndir<m_spacedim;ndir++)
        {
            PntJacDerv[ndir]  =   TensorOfArray3D<NekDouble>(m_spacedim);
            PntJacDervStd[ndir] = TensorOfArray3D<DataType>(m_spacedim);
            DervStdd[ndir]    =   Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
            DervCurv[ndir]    =   Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        }

        Array<OneD, NekDouble> locmu;
        Array<OneD, NekDouble> locDmuDT;
        Array<OneD, Array<OneD, NekDouble>> locVars(nvariable);
        TensorOfArray3D<NekDouble> locDerv(m_spacedim);
        for(int ndir=0; ndir<m_spacedim;ndir++)
        {
            locDerv[ndir] = Array<OneD, Array<OneD, NekDouble>>(nvariable);
        }

        int nElmtCoefOld = -1;
        for(int ne=0; ne<nTotElmt;ne++)
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
                    PntJacCons[ndir] = Array<OneD, Array<OneD, NekDouble>> (nElmtPnt);
                    PntJacConsStd[ndir] = Array<OneD, Array<OneD, DataType>> (nElmtPnt);
                    for(int i=0; i<nElmtPnt;i++)
                    {
                        PntJacCons[ndir][i] = Array<OneD, NekDouble>(nVar2);
                        PntJacConsStd[ndir][i] = Array<OneD, DataType>(nVar2);
                    }
                    
                    for(int ndir1=0; ndir1<m_spacedim;ndir1++)
                    {
                        PntJacDerv[ndir][ndir1] = Array<OneD, Array<OneD, NekDouble>> (nElmtPnt);
                        PntJacDervStd[ndir][ndir1] = Array<OneD, Array<OneD, DataType>> (nElmtPnt);
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
            
            if(m_advectionJacFlag)
            {
                for(int nFluxDir = 0; nFluxDir < nSpaceDim; nFluxDir++)
                {
                    normals =   normal3D[nFluxDir];
                    GetFluxVectorJacDirElmt(nvariable,nElmtPnt,locVars,
                        normals,wspMat,PntJacCons[nFluxDir]);
                }
            }

            if(m_viscousJacFlag)
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
                for(int nFluxDir = 0; nFluxDir < nSpaceDim; nFluxDir++)
                {
                    normals =   normal3D[nFluxDir];
                    MinusDiffusionFluxJacPoint(nvariable,nElmtPnt,
                        locVars,locDerv,locmu,locDmuDT,normals,wspMatDrv,
                        PntJacCons[nFluxDir]);
                }
            }

            if(m_viscousJacFlag)
            {
                locmu = mu + noffset;
                for(int nFluxDir = 0; nFluxDir < nSpaceDim; nFluxDir++)
                {
                    Vmath::Fill(npoints,1.0,normalPnt[nFluxDir],1);
                    for(int nDervDir = 0; nDervDir < nSpaceDim; nDervDir++)
                    {
                        GetFluxDerivJacDirctnElmt(nvariable,nElmtPnt,nDervDir,
                            locVars,locmu,normalPnt,wspMatDrv,
                            PntJacDerv[nFluxDir][nDervDir]);
                    }
                    Vmath::Fill(npoints,0.0,normalPnt[nFluxDir],1);
                }
            }

            for(int n=0; n<nvariable;n++)
            {
                for(int m=0; m<nvariable;m++)
                {
                    int nVarOffset = m+n*nvariable;
                    GmatxData = gmtxarray[m][n]->GetBlock(ne,ne)->GetPtr();

                    for(int ndStd0 =0;ndStd0<m_spacedim;ndStd0++)
                    {
                        Vmath::Zero(nElmtPnt,ConsStdd[ndStd0],1);
                    }
                    for(int ndir =0;ndir<m_spacedim;ndir++)
                    {
                        for(int i=0; i<nElmtPnt;i++)
                        {
                            tmppnts[i] =  PntJacCons[ndir][i][nVarOffset];
                        }
                        (*expvect)[ne]->AlignVectorToCollapsedDir(ndir,
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
                            PntJacConsStd[ndir][i][nVarOffset] = 
                                DataType(ConsStdd[ndir][i]);
                        }
                    }
                }
            }

            if(m_viscousJacFlag)
            {
                for(int m=0; m<nvariable;m++)
                {
                    for(int n=0; n<nvariable;n++)
                    {
                        int nVarOffset = m+n*nvariable;
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
                                        PntJacDerv[nd0][nd1][i][nVarOffset];
                                }

                                (*expvect)[ne]->AlignVectorToCollapsedDir(
                                    nd0,tmppnts,ConsCurv);
                                for(int nd =0;nd<m_spacedim;nd++)
                                {
                                    (*expvect)[ne]->
                                        AlignVectorToCollapsedDir(nd1,
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
                                    PntJacDervStd[nd0][nd1][i][nVarOffset] = 
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

            if(m_viscousJacFlag)
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
                    int nVarOffset = m+n*nvariable;
                    GmatxData = gmtxarray[m][n]->GetBlock(ne,ne)->GetPtr();
                    Vmath::Vcopy(nElmtCoef2, 
                        tmpA = MatData + nVarOffset*nElmtCoef2Paded,1,
                        GmatxData,1);
                }
            }
        }
    }

    template<typename DataType>
    void CFSImplicit::CalcVolJacStdMat(
        TensorOfArray4D<DataType> &StdMatDataDBB,
        TensorOfArray5D<DataType> &StdMatDataDBDB)
    {
        std::shared_ptr<LocalRegions::ExpansionVector> expvect =    
            m_fields[0]->GetExp();
        int nTotElmt            = (*expvect).size();

        StdMatDataDBB = TensorOfArray4D<DataType> (nTotElmt);
        StdMatDataDBDB  = TensorOfArray5D<DataType> (nTotElmt);

        vector<DNekMatSharedPtr> VectStdDerivBase0;
        vector< TensorOfArray3D<DataType> > VectStdDerivBase_Base;
        vector< TensorOfArray4D<DataType> > VectStdDervBase_DervBase;
        DNekMatSharedPtr MatStdDerivBase0;
        Array<OneD, DNekMatSharedPtr> ArrayStdMat(m_spacedim);
        Array<OneD, Array<OneD, NekDouble>>    ArrayStdMatData(m_spacedim);
        for (int ne = 0; ne < nTotElmt; ne++)
        {
            StdRegions::StdExpansionSharedPtr stdExp;
            stdExp = (*expvect)[ne]->GetStdExp();
            StdRegions::StdMatrixKey  matkey(StdRegions::eDerivBase0,
                                stdExp->DetShapeType(), *stdExp);
            MatStdDerivBase0      =   stdExp->GetStdMatrix(matkey);

            int nTotStdExp = VectStdDerivBase0.size();
            int nFoundStdExp = -1;
            for (int i = 0; i < nTotStdExp; i++) 
            {
                if ((*VectStdDerivBase0[i])==(*MatStdDerivBase0))
                {
                    nFoundStdExp = i;
                }
            }
            if (nFoundStdExp >= 0)
            {
                StdMatDataDBB[ne] = VectStdDerivBase_Base[nFoundStdExp];
                StdMatDataDBDB[ne] = VectStdDervBase_DervBase[nFoundStdExp];
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
                    tmpStdDBB[nd0]  = Array<OneD, Array<OneD, DataType>> (nElmtPnt);
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
                            Array<OneD, Array<OneD, DataType>> (nElmtPnt);
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
    void CFSImplicit::AddMatNSBlkDiagBnd(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              TensorOfArray3D<NekDouble>                &qfield,
              TensorOfArray2D<TypeNekBlkMatSharedPtr>   &gmtxarray,
              Array<OneD, TypeNekBlkMatSharedPtr>       &TraceJac,
              Array<OneD, TypeNekBlkMatSharedPtr>       &TraceJacDeriv,
              Array<OneD, Array<OneD, DataType>>        &TraceJacDerivSign,
              TensorOfArray5D<DataType>                 &TraceIPSymJacArray)
    {
        int nvariables = inarray.size();
        GetTraceJac(inarray,qfield,TraceJac,TraceJacDeriv,TraceJacDerivSign,
            TraceIPSymJacArray);
        
        Array<OneD, TypeNekBlkMatSharedPtr > tmpJac;
        Array<OneD, Array<OneD, DataType>>  tmpSign;

        m_advObject->AddTraceJacToMat(nvariables,m_spacedim,m_fields, 
            TraceJac,gmtxarray,tmpJac,tmpSign);
    }


    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CFSImplicit::ElmtVarInvMtrx(
              Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>>  &gmtxarray,
              TypeNekBlkMatSharedPtr                            &gmtVar,
        const DataType                                          &tmpDataType)
    {
        boost::ignore_unused(tmpDataType);

        int n1d = gmtxarray.size();
        int n2d = gmtxarray[0].size();
        int nConvectiveFields = n1d;

        ASSERTL0(n1d==n2d,"ElmtVarInvMtrx requires n1d==n2d");

        Array<OneD, unsigned int> rowSizes;
        Array<OneD, unsigned int> colSizes;

        gmtxarray[0][0]->GetBlockSizes(rowSizes,colSizes);
        int nTotElmt  = rowSizes.size();
        int nElmtCoef   =    rowSizes[0]-1;
        int nElmtCoef0  =    -1;
        int blocksize = -1;

        Array<OneD, unsigned int> tmprow(1);
        TypeNekBlkMatSharedPtr tmpGmtx;

        Array<OneD, DataType>    GMatData,ElmtMatData;
        Array<OneD, DataType>    tmpArray1,tmpArray2;

        for(int nelmt = 0; nelmt < nTotElmt; nelmt++)
        {
            int nrows = gmtxarray[0][0]->GetBlock(nelmt,nelmt)->GetRows();
            int ncols = gmtxarray[0][0]->GetBlock(nelmt,nelmt)->GetColumns();
            ASSERTL0(nrows==ncols,"ElmtVarInvMtrx requires nrows==ncols");

            nElmtCoef            = nrows;

            if (nElmtCoef0!=nElmtCoef)
            {
                nElmtCoef0 = nElmtCoef;
                int nElmtCoefVar = nElmtCoef0*nConvectiveFields;
                blocksize = nElmtCoefVar*nElmtCoefVar;
                tmprow[0] = nElmtCoefVar;
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
    void CFSImplicit::GetTraceJac(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              TensorOfArray3D<NekDouble>                &qfield,
              Array<OneD, TypeNekBlkMatSharedPtr >      &TraceJac,
              Array<OneD, TypeNekBlkMatSharedPtr >      &TraceJacDeriv,
              Array<OneD,       Array<OneD, DataType>>  &TraceJacDerivSign,
              TensorOfArray5D<DataType>                 &TraceIPSymJacArray)
    {
        boost::ignore_unused(TraceJacDeriv, TraceJacDerivSign);

        int nvariables = inarray.size();
        int nTracePts  = GetTraceTotPoints();

        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble>> Fwd (nvariables);
        Array<OneD, Array<OneD, NekDouble>> Bwd (nvariables);

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
            Fwd = NullNekDoubleArrayOfArray;
            Bwd = NullNekDoubleArrayOfArray;
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

        Array<OneD, Array<OneD, NekDouble>> AdvVel(m_spacedim);

        NumCalcRiemFluxJac(nvariables, m_fields, AdvVel, inarray,qfield,
            m_bndEvaluateTime, Fwd, Bwd, FJac, BJac,TraceIPSymJacArray);

        TraceJac[0] = FJac;
        TraceJac[1] = BJac;
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CFSImplicit::NumCalcRiemFluxJac(
        const int                                         nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, const Array<OneD, NekDouble>>   &AdvVel,
        const Array<OneD, const Array<OneD, NekDouble>>   &inarray,
              TensorOfArray3D<NekDouble>                  &qfield,
        const NekDouble                                   &time,
        const Array<OneD, const Array<OneD, NekDouble>>   &Fwd,
        const Array<OneD, const Array<OneD, NekDouble>>   &Bwd,
              TypeNekBlkMatSharedPtr                      &FJac,
              TypeNekBlkMatSharedPtr                      &BJac,
              TensorOfArray5D<DataType>                   &TraceIPSymJacArray)
    {
        boost::ignore_unused(TraceIPSymJacArray);

        const NekDouble     PenaltyFactor2  =   0.0;
        int nvariables  = nConvectiveFields;
        int npoints     = GetNpoints();
        // int nPts        = npoints;
        int nTracePts   = GetTraceTotPoints();
        int nDim        = m_spacedim;

        Array<OneD, int > nonZeroIndex;

        Array<OneD, Array<OneD, NekDouble>>  tmpinarry(nvariables);
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

        Array<OneD, Array<OneD, NekDouble>> numflux(nvariables);
        for(int i = 0; i < nvariables; ++i)
        {
            numflux[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        }

        const MultiRegions::AssemblyMapDGSharedPtr  TraceMap = 
            fields[0]->GetTraceMap();
        TensorOfArray3D<NekDouble>    qBwd(nDim);
        TensorOfArray3D<NekDouble>    qFwd(nDim);
        if(m_viscousJacFlag)
        {
            for (int nd = 0; nd < nDim; ++nd)
            {
                qBwd[nd] = Array<OneD, Array<OneD, NekDouble>> (nConvectiveFields);
                qFwd[nd] = Array<OneD, Array<OneD, NekDouble>> (nConvectiveFields);
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
        Array<OneD, Array<OneD, NekDouble>>  plusFwd(nFields),plusBwd(nFields);
        Array<OneD, Array<OneD, NekDouble>>  Jacvect(nFields);
        Array<OneD, Array<OneD, NekDouble>>  FwdBnd(nFields);
        Array<OneD, Array<OneD, NekDouble>>  plusflux(nFields);
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

    void CFSImplicit::CalcTraceNumericalFlux(
            const int                                         nConvectiveFields,
            const int                                         nDim,
            const int                                         nPts,
            const int                                         nTracePts,
            const NekDouble                                   PenaltyFactor2,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, const Array<OneD, NekDouble>>   &AdvVel,
            const Array<OneD, const Array<OneD, NekDouble>>   &inarray,
            const NekDouble                                   time,
                  TensorOfArray3D<NekDouble>                  &qfield,
            const Array<OneD, const Array<OneD, NekDouble>>   &vFwd,
            const Array<OneD, const Array<OneD, NekDouble>>   &vBwd,
            const Array<OneD, const TensorOfArray2D<NekDouble>> &qFwd,
            const Array<OneD, const TensorOfArray2D<NekDouble>> &qBwd,
            const Array<OneD, NekDouble >                     &MuVarTrace,
                  Array<OneD, int >                           &nonZeroIndex,
                  Array<OneD, Array<OneD,       NekDouble>>   &traceflux)
    {
        boost::ignore_unused(nDim, nPts, PenaltyFactor2, time, qFwd, qBwd, 
            MuVarTrace);

        if (m_advectionJacFlag)
        {
            m_advObject->AdvectTraceFlux(nConvectiveFields, m_fields, AdvVel,
                inarray, traceflux, m_bndEvaluateTime,vFwd, vBwd);
        }
        else
        {
            for (int i = 0; i < nConvectiveFields; i++)
            {
                traceflux[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
            }
        }
        
        if(m_viscousJacFlag)
        {
            Array<OneD, Array<OneD, NekDouble > > visflux(nConvectiveFields);
            for(int i = 0; i < nConvectiveFields; i++)
            {
                visflux[i]  =    Array<OneD, NekDouble>(nTracePts,0.0);
            }

            m_diffusion->DiffuseTraceFlux(fields, inarray, qfield, 
                NullNekDoubleTensorOfArray3D, visflux, vFwd, vBwd, nonZeroIndex);
            for(int i = 0; i < nConvectiveFields; i++)
            {
                Vmath::Vsub(nTracePts,traceflux[i],1,visflux[i],1,
                    traceflux[i],1);
            }
        }
    }

    void CFSImplicit::CalcPreconMatBRJCoeff(
        const Array<OneD, const Array<OneD, NekDouble>>  &inarray,
        Array<OneD, Array<OneD, SNekBlkMatSharedPtr>>    &gmtxarray,
        SNekBlkMatSharedPtr              &gmtVar,
        Array<OneD, SNekBlkMatSharedPtr> &TraceJac,
        Array<OneD, SNekBlkMatSharedPtr> &TraceJacDeriv,
        Array<OneD, Array<OneD, NekSingle>> &TraceJacDerivSign,
        TensorOfArray4D<NekSingle>          &TraceJacArray,
        TensorOfArray4D<NekSingle>          &TraceJacDerivArray,
        TensorOfArray5D<NekSingle>          &TraceIPSymJacArray)
    {
        TensorOfArray3D<NekDouble> qfield;

        if(m_viscousJacFlag)
        {
            CalcPhysDeriv(inarray,qfield);
        }

        NekSingle zero =0.0;
        Fill2DArrayOfBlkDiagonalMat(gmtxarray,zero);

        AddMatNSBlkDiagVol(inarray,qfield,gmtxarray,m_stdSMatDataDBB,
            m_stdSMatDataDBDB);

        AddMatNSBlkDiagBnd(inarray,qfield,gmtxarray,TraceJac,
            TraceJacDeriv,TraceJacDerivSign,TraceIPSymJacArray);

        MultiplyElmtInvMassPlusSource(gmtxarray,m_TimeIntegLambda,zero);

        ElmtVarInvMtrx(gmtxarray,gmtVar,zero);

        TransTraceJacMatToArray(TraceJac,TraceJacDeriv,TraceJacArray, 
            TraceJacDerivArray);
    }

    void CFSImplicit::v_MinusDiffusionFluxJacPoint(
        const int                                       nConvectiveFields,
        const int                                       nElmtPnt,
        const Array<OneD, const Array<OneD, NekDouble>> &locVars,
        const TensorOfArray3D<NekDouble>                &locDerv,
        const Array<OneD, NekDouble>                    &locmu,
        const Array<OneD, NekDouble>                    &locDmuDT,
        const Array<OneD, NekDouble>                    &normals,
              DNekMatSharedPtr                          &wspMat,
              Array<OneD,       Array<OneD, NekDouble>> &PntJacArray)
    {
        boost::ignore_unused(nConvectiveFields, nElmtPnt, locVars, locDerv,
                locmu, locDmuDT, normals, wspMat, PntJacArray);
        // Do nothing by default
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CFSImplicit::MultiplyElmtInvMassPlusSource(
              Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>>  &gmtxarray,
        const NekDouble                                         dtlamda,
        const DataType                                          tmpDataType)
    {
        boost::ignore_unused(tmpDataType);

        MultiRegions::ExpListSharedPtr explist = m_fields[0];
        std::shared_ptr<LocalRegions::ExpansionVector> pexp = explist->GetExp();
        int nTotElmt            = (*pexp).size();
        int nConvectiveFields = m_fields.size();

        NekDouble Negdtlamda    =   -dtlamda;

        Array<OneD, NekDouble> pseudotimefactor(nTotElmt,0.0);
        Vmath::Fill(nTotElmt,Negdtlamda,pseudotimefactor,1);

        Array<OneD, DataType>    GMatData;
        for(int m = 0; m < nConvectiveFields; m++)
        {
            for(int n = 0; n < nConvectiveFields; n++)
            {
                for(int nelmt = 0; nelmt < nTotElmt; nelmt++)
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

        for(int nelmt = 0; nelmt < nTotElmt; nelmt++)
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
    void CFSImplicit::TransTraceJacMatToArray(
        const Array<OneD, TypeNekBlkMatSharedPtr> &TraceJac,
        const Array<OneD, TypeNekBlkMatSharedPtr> &TraceJacDeriv,
              TensorOfArray4D<DataType>           &TraceJacArray,
              TensorOfArray4D<DataType>           &TraceJacDerivArray)
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
                    TraceJacArray[nlr][m] = Array<OneD, Array<OneD, DataType>>(nvar1Jac);
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

     void CFSImplicit::v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr            &explist,
        const Array<OneD, const Array<OneD, NekDouble>> &normals,
        const int                                       nDervDir,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              TensorOfArray5D<NekDouble>                &ElmtJacArray,
        const int                                       nFluxDir)
    {
        boost::ignore_unused(explist, normals, nDervDir, inarray, ElmtJacArray,
            nFluxDir);
        NEKERROR(ErrorUtil::efatal, "v_GetFluxDerivJacDirctn not coded");
    }

    void CFSImplicit::v_GetFluxDerivJacDirctnElmt(
        const int                                       nConvectiveFields,
        const int                                       nElmtPnt,
        const int                                       nDervDir,
        const Array<OneD, const Array<OneD, NekDouble>> &locVars,
        const Array<OneD, NekDouble>                    &locmu,
        const Array<OneD, const Array<OneD, NekDouble>> &locnormal,
              DNekMatSharedPtr                          &wspMat,
              Array<OneD,       Array<OneD, NekDouble>> &PntJacArray)
    {
        boost::ignore_unused(nConvectiveFields, nElmtPnt, nDervDir, locVars,
            locmu, locnormal, wspMat, PntJacArray);
        NEKERROR(ErrorUtil::efatal, "v_GetFluxDerivJacDirctn not coded");
    }
    
    void CFSImplicit::v_GetFluxDerivJacDirctn(
        const MultiRegions::ExpListSharedPtr            &explist,
        const Array<OneD, const Array<OneD, NekDouble>> &normals,
        const int                                       nDervDir,
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        TensorOfArray2D<DNekMatSharedPtr>               &ElmtJac)
    {
        boost::ignore_unused(explist, normals, nDervDir, inarray, ElmtJac);
    }

    void CFSImplicit::GetFluxVectorJacDirElmt(
        const int                                       nConvectiveFields,
        const int                                       nElmtPnt,
        const Array<OneD, const Array<OneD, NekDouble>> &locVars,
        const Array<OneD, NekDouble>                    &normals,
              DNekMatSharedPtr                          &wspMat,
              Array<OneD,       Array<OneD, NekDouble>> &PntJacArray)
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

    void CFSImplicit::GetFluxVectorJacPoint(
        const int                       nConvectiveFields,
        const Array<OneD, NekDouble>    &conservVar, 
        const Array<OneD, NekDouble>    &normals, 
              DNekMatSharedPtr          &fluxJac)
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
            PointFluxJacobianPoint(conservVar,normals,fluxJac,
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
            
            PointFluxJacobianPoint(PointFwd,normals,PointFJac3D,
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



    /**
     * @brief Add the diffusions terms to the right-hand side
     * Similar to DoDiffusion() but with outarray in coefficient space
     */
    void CFSImplicit::DoDiffusionCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
              Array<OneD,       Array<OneD, NekDouble>> &outarray,
        const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &pBwd)
    {
        v_DoDiffusionCoeff(inarray, outarray, pFwd, pBwd);
    }



    void CFSImplicit::MatrixMultiplyMatrixFreeCoeff(
            const Array<OneD, const NekDouble>  &inarray,
                  Array<OneD,       NekDouble>  &out,
            const bool                          &flag)
    {
        boost::ignore_unused(flag);
        const Array<OneD, const NekDouble> refsol =
            m_nonlinsol->GetRefSolution();
        const Array<OneD, const NekDouble> refres =
            m_nonlinsol->GetRefResidual();
        const Array<OneD, const NekDouble> refsource =
            m_nonlinsol->GetRefSourceVec();

        NekDouble eps = m_jacobiFreeEps;
        unsigned int nTotalGlobal     = inarray.size();
        NekDouble magninarray = Vmath::Dot(nTotalGlobal,inarray,inarray);
        m_comm->AllReduce(magninarray, Nektar::LibUtilities::ReduceSum);
        eps *= sqrt( (sqrt(m_inArrayNorm) + 1.0)/magninarray);

        NekDouble oeps = 1.0 / eps;
        unsigned int ntotal     = inarray.size();
        Array<OneD, NekDouble> solplus{ntotal, 0.0};
        Array<OneD, NekDouble> resplus{ntotal, 0.0};

        Vmath::Svtvp(ntotal, eps, inarray, 1, refsol, 1, solplus, 1);

        NonlinSysEvaluatorCoeff(solplus, resplus, flag, refsource);
        
        Vmath::Vsub(ntotal, resplus, 1, refres, 1, out,1);
        Vmath::Smul(ntotal, oeps, out, 1, out, 1);
       
        return;
    }

    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CFSImplicit::TranSamesizeBlkDiagMatIntoArray(
        const TypeNekBlkMatSharedPtr    &BlkMat,
              TensorOfArray3D<DataType> &MatArray)
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
    void CFSImplicit::Fill2DArrayOfBlkDiagonalMat(
              Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr>>  &gmtxarray,
        const DataType                                          valu)
    {

        int n1d = gmtxarray.size();

        for(int n1 = 0; n1 < n1d; ++n1)
        {
            Fill1DArrayOfBlkDiagonalMat(gmtxarray[n1],valu);
        }
    }

    
    template<typename DataType, typename TypeNekBlkMatSharedPtr>
    void CFSImplicit::Fill1DArrayOfBlkDiagonalMat( 
              Array<OneD, TypeNekBlkMatSharedPtr>   &gmtxarray,
        const DataType                              valu)
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


    // Currently duplacate in compressibleFlowSys
    // if fsw=+-1 calculate the steger-Warming flux vector splitting flux Jacobian
    // if fsw=0   calculate the Jacobian of the exact flux
    // efix is the numerical flux entropy fix parameter
    void CFSImplicit::PointFluxJacobianPoint(
            const Array<OneD, NekDouble> &Fwd,
            const Array<OneD, NekDouble> &normals,
                  DNekMatSharedPtr       &FJac,
            const NekDouble efix, const NekDouble fsw)
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

    bool CFSImplicit::UpdateTimeStepCheck()
    {
        bool flag = (m_time + m_timestep > m_fintime && m_fintime > 0.0) || 
                (m_checktime && m_time + m_timestep - m_lastCheckTime >=
                        m_checktime);
        if (m_explicitAdvection)
        {
            flag = true;
        }
        else
        {
            flag = flag || m_preconCfs->UpdatePreconMatCheck(NullNekDouble1DArray, 
                m_TimeIntegLambda);
        }
        return flag;
    }

}
