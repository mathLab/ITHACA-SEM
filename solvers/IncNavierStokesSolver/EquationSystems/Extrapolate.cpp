///////////////////////////////////////////////////////////////////////////////
//
// File: Extrapolate.cpp
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
// Description: Abstract base class for Extrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/Extrapolate.h>
#include <LibUtilities/Communication/Comm.h>

using namespace std;

namespace Nektar
{
    NekDouble Extrapolate::StifflyStable_Betaq_Coeffs[3][3] = {
        { 1.0,  0.0, 0.0},{ 2.0, -1.0, 0.0},{ 3.0, -3.0, 1.0}};
    NekDouble Extrapolate::StifflyStable_Alpha_Coeffs[3][3] = {
        { 1.0,  0.0, 0.0},{ 2.0, -0.5, 0.0},{ 3.0, -1.5, 1.0/3.0}};
    NekDouble Extrapolate::StifflyStable_Gamma0_Coeffs[3] = {
          1.0,  1.5, 11.0/6.0};

    ExtrapolateFactory& GetExtrapolateFactory()
    {
        static ExtrapolateFactory instance;
        return instance;
    }

    Extrapolate::Extrapolate(
        const LibUtilities::SessionReaderSharedPtr  pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
        MultiRegions::ExpListSharedPtr              pPressure,
        const Array<OneD, int>                      pVel,
        const SolverUtils::AdvectionSharedPtr       advObject)
        : m_session(pSession),
          m_fields(pFields),
          m_pressure(pPressure),
          m_velocity(pVel),
          m_advObject(advObject)
    {
        m_session->LoadParameter("TimeStep", m_timestep,   0.01);
        m_comm = m_session->GetComm();
    }

    Extrapolate::~Extrapolate()
    {
    }

    std::string Extrapolate::def =
        LibUtilities::SessionReader::RegisterDefaultSolverInfo(
            "StandardExtrapolate", "StandardExtrapolate");

    /**
     *
     */
    void Extrapolate::AddDuDt(void)
    {
        if(m_numHBCDof)
        {
            // Update velocity BF at n+1 (actually only needs doing if
            // velocity is time dependent on HBCs)
            IProductNormVelocityBCOnHBC(m_iprodnormvel[m_intSteps]);

            //Calculate acceleration term at level n based on previous steps
            AccelerationBDF(m_iprodnormvel);

            // Subtract acceleration term off m_pressureHBCs[nlevels-1]
            Vmath::Svtvp(m_numHBCDof, -1.0/m_timestep,
                         m_iprodnormvel[m_intSteps],  1,
                         m_pressureHBCs[m_intSteps-1], 1,
                         m_pressureHBCs[m_intSteps-1], 1);
        }
    }

    /**
     *
     */
    void Extrapolate::AddVelBC(void)
    {
        if(m_numHBCDof)
        {
            int order = std::min(m_pressureCalls,m_intSteps);

            // Update velocity BF at n+1 (actually only needs doing if
            // velocity is time dependent on HBCs)
            IProductNormVelocityBCOnHBC(m_iprodnormvel[0]);

            // Subtract acceleration term off m_pressureHBCs[nlevels-1]
            Vmath::Svtvp(m_numHBCDof,
                         -1.0*StifflyStable_Gamma0_Coeffs[order-1]/m_timestep,
                         m_iprodnormvel[0],  1,
                         m_pressureHBCs[m_intSteps-1], 1,
                         m_pressureHBCs[m_intSteps-1], 1);
        }
    }

    /**
     * Unified routine for calculation high-oder terms
     */
    void Extrapolate::v_CalcNeumannPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        const Array<OneD, const Array<OneD, NekDouble> >  &N,
        NekDouble kinvis)
    {
        int n, cnt;

        Array<OneD, NekDouble> Pvals;

        Array<OneD, Array<OneD, NekDouble> > Velocity(m_curl_dim);
        Array<OneD, Array<OneD, NekDouble> > Advection(m_bnd_dim);

        Array<OneD, Array<OneD, NekDouble> > BndValues(m_bnd_dim);
        Array<OneD, Array<OneD, NekDouble> > Q(m_curl_dim);

        MultiRegions::ExpListSharedPtr BndElmtExp;
        for(n = cnt = 0; n < m_PBndConds.size(); ++n)
        {
            // High order boundary condition;
            if((m_hbcType[n] == eHBCNeumann)||(m_hbcType[n] == eConvectiveOBC))
            {
                m_fields[0]->GetBndElmtExpansion(n, BndElmtExp, false);
                int nqb = m_PBndExp[n]->GetTotPoints();
                int nq  = BndElmtExp->GetTotPoints();

                for(int i = 0; i < m_bnd_dim; i++)
                {
                    BndValues[i] = Array<OneD, NekDouble> (nqb,0.0);
                }

                for(int i = 0; i < m_curl_dim; i++)
                {
                    Q[i] = Array<OneD, NekDouble> (nq,0.0);
                }

                // Obtaining fields on BndElmtExp
                for(int i = 0; i < m_curl_dim; i++)
                {
                    m_fields[0]->ExtractPhysToBndElmt(n, fields[i],Velocity[i]);
                }

                if(N.size()) // not required for some extrapolation
                {
                    for(int i = 0; i < m_bnd_dim; i++)
                    {
                        m_fields[0]->ExtractPhysToBndElmt(n, N[i], Advection[i]);
                    }
                }

                // CurlCurl
                BndElmtExp->CurlCurl(Velocity, Q);

                // Mounting advection component into the high-order condition
                for(int i = 0; i < m_bnd_dim; i++)
                {
                    MountHOPBCs(nq, kinvis,Q[i],Advection[i]);
                }

                Pvals = (m_pressureHBCs[m_intSteps-1]) + cnt;

                // Getting values on the boundary and filling the pressure bnd
                // expansion. Multiplication by the normal is required
                for(int i = 0; i < m_bnd_dim; i++)
                {
                    m_fields[0]->ExtractElmtToBndPhys(n, Q[i],BndValues[i]);
                }

                m_PBndExp[n]->NormVectorIProductWRTBase(BndValues, Pvals);

                // Get offset for next terms
                cnt += m_PBndExp[n]->GetNcoeffs();
            }
        }
    }

    // do nothing unless otherwise defined.
    void Extrapolate::v_CorrectPressureBCs( const Array<OneD, NekDouble>  &pressure)
    {
    }

    // do nothing unless otherwise defined.
    void  Extrapolate::v_AddNormVelOnOBC(const int nbcoeffs, const int nreg,
                                         Array<OneD, Array<OneD, NekDouble> > &u)
    {
    }

    void Extrapolate::CalcOutflowBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        NekDouble kinvis)
    {
        if(!m_houtflow.get())
        {
           return;
        }

        Array<OneD, Array<OneD, NekDouble> > Velocity(m_curl_dim);

        MultiRegions::ExpListSharedPtr BndElmtExp;
        int cnt    = 0;

        // Evaluate robin primitive coefficient here so they can be
        // updated whem m_int > 1 Currently not using this update
        // since we only using u^n at outflow instead of BDF rule.
        UpdateRobinPrimCoeff();

        for(int n  = 0; n < m_PBndConds.size(); ++n)
        {
            if((m_hbcType[n] == eOBC)||(m_hbcType[n] == eConvectiveOBC))
            {
                // Get expansion with element on this boundary
                m_fields[0]->GetBndElmtExpansion(n, BndElmtExp, false);
                int nqb = m_PBndExp[n]->GetTotPoints();
                int nq  = BndElmtExp->GetTotPoints();

                // Get velocity and extrapolate
                for(int i = 0; i < m_curl_dim; i++)
                {
                    m_fields[0]->ExtractPhysToBndElmt(n, fields[i],
                            m_houtflow->m_outflowVel[cnt][i][m_intSteps-1]);
                    ExtrapolateArray(m_houtflow->m_outflowVel[cnt][i]);
                    Velocity[i] = m_houtflow->m_outflowVel[cnt][i][m_intSteps-1];

                }

                // Homogeneous case needs conversion to physical space
                if ( m_fields[0]->GetWaveSpace())
                {
                    for(int i = 0; i < m_curl_dim; i++)
                    {
                        BndElmtExp->HomogeneousBwdTrans(Velocity[i],
                                                        Velocity[i]);
                    }
                    BndElmtExp->SetWaveSpace(false);
                }

                // Get normal vector
                Array<OneD, Array<OneD, NekDouble> > normals;
                m_fields[0]->GetBoundaryNormals(n, normals);

                // Calculate n.gradU.n, div(U)
                Array<OneD, NekDouble> nGradUn (nqb, 0.0);
                Array<OneD, NekDouble> divU (nqb, 0.0);
                Array<OneD, Array<OneD, NekDouble> > grad(m_curl_dim);
                Array<OneD, NekDouble> bndVal (nqb, 0.0);
                for( int i = 0; i < m_curl_dim; i++)
                {
                    grad[i] = Array<OneD, NekDouble> (nq, 0.0);
                }
                for( int i = 0; i < m_curl_dim; i++)
                {
                    if( m_curl_dim == 2)
                    {
                        BndElmtExp->PhysDeriv(Velocity[i], grad[0], grad[1]);
                    }
                    else
                    {
                        BndElmtExp->PhysDeriv(Velocity[i], grad[0], grad[1],
                                                                    grad[2]);
                    }

                    for( int j = 0; j < m_curl_dim; j++)
                    {
                        m_fields[0]->ExtractElmtToBndPhys(n, grad[j],bndVal);
                        // div(U) = gradU_ii
                        if ( i == j)
                        {
                            Vmath::Vadd(nqb , divU, 1, bndVal, 1, divU, 1);
                        }
                        // n.gradU.n = gradU_ij n_i n_j
                        Vmath::Vmul(nqb , normals[i], 1, bndVal, 1,
                                                       bndVal, 1);
                        Vmath::Vvtvp(nqb , normals[j], 1, bndVal, 1,
                                         nGradUn, 1, nGradUn, 1);
                    }
                }

                // Obtain u at the boundary
                Array<OneD, Array<OneD, NekDouble> > u (m_curl_dim);
                for( int i = 0; i < m_curl_dim; i++)
                {
                    u[i] = Array<OneD, NekDouble> (nqb, 0.0);
                    m_fields[0]->ExtractElmtToBndPhys(n, Velocity[i],u[i]);
                }

                // Calculate u.n and u^2
                Array<OneD, NekDouble> un (nqb, 0.0);
                Array<OneD, NekDouble> u2 (nqb, 0.0);
                for( int i = 0; i < m_curl_dim; i++)
                {
                    Vmath::Vvtvp(nqb, normals[i], 1, u[i], 1,
                                 un, 1, un, 1);
                    Vmath::Vvtvp(nqb, u[i], 1, u[i], 1,
                                 u2, 1, u2, 1);
                }

                // Calculate S_0(u.n) = 0.5*(1-tanh(u.n/*U0*delta))
                Array<OneD, NekDouble> S0 (nqb, 0.0);
                for( int i = 0; i < nqb; i++)
                {
                    S0[i] = 0.5*(1.0-tanh(un[i]/(m_houtflow->m_U0*m_houtflow->m_delta)));
                }

                // Calculate E(n,u) = ((theta+alpha2)*0.5*(u^2)n +
                //                    (1-theta+alpha1)*0.5*(n.u)u ) * S_0(u.n)
                NekDouble k1 = 0.5*(m_houtflow->m_obcTheta +
                                    m_houtflow->m_obcAlpha2);
                NekDouble k2 = 0.5*(1-m_houtflow->m_obcTheta +
                                    m_houtflow->m_obcAlpha1);

                Array<OneD, Array<OneD, NekDouble> > E (m_curl_dim);
                for( int i = 0; i < m_curl_dim; i++)
                {
                    E[i] = Array<OneD, NekDouble> (nqb, 0.0);

                    Vmath::Smul(nqb, k1, u2, 1, E[i], 1);
                    Vmath::Vmul(nqb, E[i], 1, normals[i], 1, E[i], 1);
                    // Use bndVal as a temporary storage
                    Vmath::Smul(nqb, k2, un, 1, bndVal, 1);
                    Vmath::Vvtvp(nqb, u[i], 1, bndVal, 1, E[i], 1, E[i], 1);
                    Vmath::Vmul(nqb, E[i], 1, S0, 1, E[i], 1);
                }

                // if non-zero forcing is provided we want to subtract
                // value if we want to interpret values as being the
                // desired pressure value. This is now precribed from
                // the velocity forcing to be consistent with the
                // paper except f_b = -f_b

                // Calculate (E(n,u) + f_b).n
                Array<OneD, NekDouble> En (nqb, 0.0);
                for( int i = 0; i < m_bnd_dim; i++)
                {
                    // Use bndVal as temporary
                    Vmath::Vsub(nqb,E[i],1,m_houtflow->
                                m_UBndExp[i][n]->GetPhys(),
                                1,  bndVal, 1);

                    Vmath::Vvtvp(nqb, normals[i], 1, bndVal, 1,
                                 En, 1, En, 1);

                }

                // Calculate pressure bc = kinvis*n.gradU.n - E.n + f_b.n
                Array<OneD, NekDouble> pbc (nqb, 0.0);
                Vmath::Svtvm( nqb, kinvis, nGradUn, 1, En, 1, pbc, 1);

                if(m_hbcType[n] == eOBC)
                {

                    if ( m_PBndExp[n]->GetWaveSpace())
                    {
                        m_PBndExp[n]->HomogeneousFwdTrans(pbc, bndVal);
                        m_PBndExp[n]->FwdTrans(bndVal,
                                               m_PBndExp[n]->UpdateCoeffs());
                    }
                    else
                    {
                        m_PBndExp[n]->FwdTrans(pbc,
                                               m_PBndExp[n]->UpdateCoeffs());
                    }
                }
                else if(m_hbcType[n] == eConvectiveOBC) // add outflow values to calculation from HBC
                {
                    int nbcoeffs = m_PBndExp[n]->GetNcoeffs();
                    Array<OneD, NekDouble> bndCoeffs (nbcoeffs, 0.0);
                    if ( m_PBndExp[n]->GetWaveSpace())
                    {
                        m_PBndExp[n]->HomogeneousFwdTrans(pbc, bndVal);
                        m_PBndExp[n]->IProductWRTBase(bndVal,bndCoeffs);
                    }
                    else
                    {
                        m_PBndExp[n]->IProductWRTBase(pbc,bndCoeffs);
                    }
                    // Note we have the negative of what is in the Dong paper in bndVal
                    Vmath::Svtvp(nbcoeffs,m_houtflow->m_pressurePrimCoeff[n],
                                 bndCoeffs, 1,m_PBndExp[n]->UpdateCoeffs(),1,
                                 m_PBndExp[n]->UpdateCoeffs(),1);

                    // evaluate u^n at outflow boundary for velocity BC
                    for( int i = 0; i < m_curl_dim; i++)
                    {
                        m_fields[0]->ExtractElmtToBndPhys(n,
                                        m_houtflow->
                                        m_outflowVel[cnt][i][0],
                                        m_houtflow->
                                        m_outflowVelBnd[cnt][i][m_intSteps-1]);

                        EvaluateBDFArray(m_houtflow->m_outflowVelBnd[cnt][i]);

                        // point u[i] to BDF evalauted value \hat{u}
                        u[i] = m_houtflow->m_outflowVelBnd[cnt][i]
                            [m_intSteps-1];
                    }

                    // Add normal velocity if weak pressure
                    // formulation. In this case there is an
                    // additional \int \hat{u}.n ds on the outflow
                    // boundary since we use the inner product wrt
                    // deriv of basis in pressure solve.
                    AddNormVelOnOBC(cnt, n, u);
                }

                // Calculate velocity boundary conditions
                if(m_hbcType[n] == eOBC)
                {
                    //    = (pbc n - kinvis divU n)
                    Vmath::Smul(nqb, kinvis, divU, 1, divU, 1);
                    Vmath::Vsub(nqb, pbc, 1, divU, 1, bndVal, 1);
                }
                else if (m_hbcType[n] == eConvectiveOBC)
                {
                    //    = (-kinvis divU n)
                    Vmath::Smul(nqb, -1.0*kinvis, divU, 1, bndVal, 1);

                    // pbc  needs to be added after pressure solve
                }

                for(int i = 0; i < m_curl_dim; ++i)
                {
                    // Reuse divU  -> En
                    Vmath::Vvtvp( nqb, normals[i], 1, bndVal, 1, E[i], 1,
                                  divU, 1);
                    // - f_b
                    Vmath::Vsub( nqb, divU, 1,
                                 m_houtflow->m_UBndExp[i][n]->GetPhys(),
                                 1, divU, 1);
                    // * 1/kinvis
                    Vmath::Smul(nqb, 1.0/kinvis, divU, 1, divU, 1);

                    if(m_hbcType[n] == eConvectiveOBC)
                    {
                        Vmath::Svtvp(nqb,m_houtflow->m_velocityPrimCoeff[i][n],
                                     u[i], 1,divU,1,divU,1);
                    }

                    if ( m_houtflow->m_UBndExp[i][n]->GetWaveSpace())
                    {
                        m_houtflow->m_UBndExp[i][n]->HomogeneousFwdTrans(divU,
                                                                         divU);
                    }

                    m_houtflow->m_UBndExp[i][n]->IProductWRTBase(divU,
                                 m_houtflow->m_UBndExp[i][n]->UpdateCoeffs());

                }

                // Get offset for next terms
                cnt++;
            }
        }
    }


    void Extrapolate::AddPressureToOutflowBCs(NekDouble kinvis)
    {
        if(!m_houtflow.get())
        {
            return;
        }


        for(int n  = 0; n < m_PBndConds.size(); ++n)
        {
            if(m_hbcType[n] == eConvectiveOBC)
            {
                int nqb = m_PBndExp[n]->GetTotPoints();
                int ncb = m_PBndExp[n]->GetNcoeffs();

                m_pressure->FillBndCondFromField(n);
                Array<OneD, NekDouble> pbc(nqb);

                m_PBndExp[n]->BwdTrans(m_PBndExp[n]->GetCoeffs(), pbc);

                if (m_PBndExp[n]->GetWaveSpace())
                {
                    m_PBndExp[n]->HomogeneousBwdTrans(pbc, pbc);
                }

                Array<OneD, NekDouble> wk(nqb);
                Array<OneD, NekDouble> wk1(ncb);

                // Get normal vector
                Array<OneD, Array<OneD, NekDouble> > normals;
                m_fields[0]->GetBoundaryNormals(n, normals);

                //  Add  1/kinvis * (pbc n )
                for(int i = 0; i < m_curl_dim; ++i)
                {
                    Vmath::Vmul(nqb, normals[i], 1, pbc, 1, wk, 1);

                    Vmath::Smul(nqb, 1.0/kinvis, wk, 1, wk, 1);

                    if (m_houtflow->m_UBndExp[i][n]->GetWaveSpace())
                    {
                        m_houtflow->m_UBndExp[i][n]->
                                    HomogeneousFwdTrans(wk, wk);
                    }
                    m_houtflow->m_UBndExp[i][n]->IProductWRTBase(wk,wk1);

                    Vmath::Vadd(ncb, wk1,1,
                                m_houtflow->m_UBndExp[i][n]->GetCoeffs(), 1,
                                m_houtflow->m_UBndExp[i][n]->UpdateCoeffs(),1);

                }
            }
        }
    }


    void Extrapolate::IProductNormVelocityOnHBC(
                       const Array<OneD, const Array<OneD, NekDouble> >  &Vel,
                       Array<OneD, NekDouble> &IProdVn)
    {
        int i,n,cnt;
        Array<OneD, NekDouble> IProdVnTmp;
        Array<OneD, Array<OneD, NekDouble> > velbc(m_bnd_dim);

        for(n = cnt = 0; n < m_PBndConds.size(); ++n)
        {
            // High order boundary condition;
            if(m_hbcType[n] == eHBCNeumann)
            {
                for(i = 0; i < m_bnd_dim; ++i)
                {
                    m_fields[0]->ExtractPhysToBnd(n, Vel[i], velbc[i]);
                }
                IProdVnTmp = IProdVn + cnt;
                m_PBndExp[n]->NormVectorIProductWRTBase(velbc, IProdVnTmp);
                cnt += m_PBndExp[n]->GetNcoeffs();
            }
            else if(m_hbcType[n] == eConvectiveOBC) // skip over conective OBC
            {
                cnt += m_PBndExp[n]->GetNcoeffs();
            }
        }
    }

    void Extrapolate::IProductNormVelocityBCOnHBC(Array<OneD, NekDouble> &IProdVn)
    {

        if(!m_HBCnumber)
        {
            return;
        }
        int i,n,cnt;
        Array<OneD, NekDouble> IProdVnTmp;
        Array<OneD, Array<OneD, NekDouble> > velbc(m_bnd_dim);
        Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> > VelBndExp(m_bnd_dim);
        for(i = 0; i < m_bnd_dim; ++i)
        {
            VelBndExp[i] = m_fields[m_velocity[i]]->GetBndCondExpansions();
        }

        for(n = cnt = 0; n < m_PBndConds.size(); ++n)
        {
            // High order boundary condition;
            if(m_hbcType[n] == eHBCNeumann)
            {
                for(i = 0; i < m_bnd_dim; ++i)
                {
                    velbc[i] = Array<OneD, NekDouble>
                        (VelBndExp[i][n]->GetTotPoints(), 0.0);
                    VelBndExp[i][n]->SetWaveSpace(
                            m_fields[m_velocity[i]]->GetWaveSpace());
                    VelBndExp[i][n]->BwdTrans(VelBndExp[i][n]->GetCoeffs(),
                                              velbc[i]);
                }
                IProdVnTmp = IProdVn + cnt;
                m_PBndExp[n]->NormVectorIProductWRTBase(velbc, IProdVnTmp);
                cnt += m_PBndExp[n]->GetNcoeffs();
            }
            else if(m_hbcType[n] == eConvectiveOBC)
            {
                // skip over convective OBC
                cnt += m_PBndExp[n]->GetNcoeffs();
            }
        }
    }

    /**
     * Function to roll time-level storages to the next step layout.
     * The stored data associated with the oldest time-level
     * (not required anymore) are moved to the top, where they will
     * be overwritten as the solution process progresses.
     */
    void Extrapolate::RollOver(Array<OneD, Array<OneD, NekDouble> > &input)
    {
        int  nlevels = input.size();

        Array<OneD, NekDouble> tmp;

        tmp = input[nlevels-1];

        for(int n = nlevels-1; n > 0; --n)
        {
            input[n] = input[n-1];
        }

        input[0] = tmp;
    }


    /**
     * Initialize HOBCs
     */
    void Extrapolate::GenerateHOPBCMap(const LibUtilities::SessionReaderSharedPtr& pSession)
    {
        m_PBndConds   = m_pressure->GetBndConditions();
        m_PBndExp     = m_pressure->GetBndCondExpansions();

        int cnt, n;

        // Storage array for high order pressure BCs
        m_pressureHBCs = Array<OneD, Array<OneD, NekDouble> > (m_intSteps);
        m_iprodnormvel = Array<OneD, Array<OneD, NekDouble> > (m_intSteps + 1);

        // Get useful values for HOBCs
        m_HBCnumber = 0;
        m_numHBCDof = 0;

        int outHBCnumber = 0;
        int numOutHBCPts = 0;

        m_hbcType = Array<OneD, HBCType> (m_PBndConds.size(),eNOHBC);
        for( n = 0; n < m_PBndConds.size(); ++n)
        {
            // High order boundary Neumann Condiiton
            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
            {
                m_hbcType[n] = eHBCNeumann;
                m_numHBCDof += m_PBndExp[n]->GetNcoeffs();
                m_HBCnumber += m_PBndExp[n]->GetExpSize();
            }

            // High order outflow convective condition
            if(m_PBndConds[n]->GetBoundaryConditionType() ==
               SpatialDomains::eRobin &&
               boost::iequals(m_PBndConds[n]->GetUserDefined(),
                              "HOutflow"))
            {
                m_hbcType[n]  = eConvectiveOBC;
                m_numHBCDof  += m_PBndExp[n]->GetNcoeffs();
                m_HBCnumber  += m_PBndExp[n]->GetExpSize();
                numOutHBCPts += m_PBndExp[n]->GetTotPoints();
                outHBCnumber++;
            }
            // High order outflow boundary condition;
            else if(boost::iequals(m_PBndConds[n]->GetUserDefined(),
                                   "HOutflow"))
            {
                m_hbcType[n] = eOBC;
                numOutHBCPts += m_PBndExp[n]->GetTotPoints();
                outHBCnumber++;
            }
        }

        m_iprodnormvel[0] = Array<OneD, NekDouble>(m_numHBCDof, 0.0);
        for(n = 0; n < m_intSteps; ++n)
        {
            m_pressureHBCs[n]   = Array<OneD, NekDouble>(m_numHBCDof, 0.0);
            m_iprodnormvel[n+1] = Array<OneD, NekDouble>(m_numHBCDof, 0.0);
        }

        m_pressureCalls = 0;

        switch(m_pressure->GetExpType())
        {
        case MultiRegions::e2D:
            {
                m_curl_dim = 2;
                m_bnd_dim  = 2;
            }
            break;
        case MultiRegions::e3DH1D:
            {
                m_curl_dim = 3;
                m_bnd_dim  = 2;
            }
            break;
        case MultiRegions::e3DH2D:
            {
                m_curl_dim = 3;
                m_bnd_dim  = 1;
            }
            break;
        case MultiRegions::e3D:
            {
                m_curl_dim = 3;
                m_bnd_dim  = 3;
            }
            break;
        default:
            ASSERTL0(0,"Dimension not supported");
            break;
        }

        // Initialise storage for outflow HOBCs
        if(numOutHBCPts > 0)
        {
            m_houtflow  = MemoryManager<HighOrderOutflow>::AllocateSharedPtr(numOutHBCPts,  outHBCnumber, m_curl_dim, pSession);

            MultiRegions::ExpListSharedPtr BndElmtExp;

            // set up boundary expansions link
            for (int i = 0; i < m_curl_dim; ++i)
            {
                m_houtflow->m_UBndExp[i] =
                    m_fields[m_velocity[i]]->GetBndCondExpansions();
            }

            for(n = 0, cnt = 0; n < m_PBndConds.size(); ++n)
            {
                if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"HOutflow"))
                {
                    m_houtflow->m_outflowVel[cnt] =
                        Array<OneD, Array<OneD,
                             Array<OneD, NekDouble> > > (m_curl_dim);

                    m_houtflow->m_outflowVelBnd[cnt] =
                        Array<OneD, Array<OneD,
                                 Array<OneD, NekDouble> > > (m_curl_dim);

                    m_fields[0]->GetBndElmtExpansion(n, BndElmtExp, false);
                    int nqb = m_PBndExp[n]->GetTotPoints();
                    int nq  = BndElmtExp->GetTotPoints();
                    for(int j = 0; j < m_curl_dim; ++j)
                    {
                        m_houtflow->m_outflowVel[cnt][j] =
                            Array<OneD, Array<OneD, NekDouble> > (m_intSteps);

                        m_houtflow->m_outflowVelBnd[cnt][j] =
                            Array<OneD, Array<OneD, NekDouble> > (m_intSteps);

                        for(int k = 0; k < m_intSteps; ++k)
                        {
                            m_houtflow->m_outflowVel[cnt][j][k] =
                                Array<OneD, NekDouble>(nq,0.0);
                            m_houtflow->m_outflowVelBnd[cnt][j][k] =
                                Array<OneD, NekDouble>(nqb,0.0);
                        }
                    }
                    cnt++;
                }

                // evaluate convective primitive coefficient if
                // convective OBCs are used
                if(m_hbcType[n] == eConvectiveOBC)
                {
                    // initialise convective members of
                    // HighOrderOutflow struct
                    if(m_houtflow->m_pressurePrimCoeff.size() == 0)
                    {
                        m_houtflow->m_pressurePrimCoeff =
                            Array<OneD, NekDouble>
                            (m_PBndConds.size(),0.0);
                        m_houtflow->m_velocityPrimCoeff =
                            Array<OneD, Array< OneD, NekDouble> >(m_curl_dim);


                        for(int i = 0; i < m_curl_dim; ++i)
                        {
                            m_houtflow->m_velocityPrimCoeff[i] =
                                Array<OneD, NekDouble> (m_PBndConds.
                                                        size(),0.0);
                        }
                    }

                    LibUtilities::Equation coeff =
                        std::static_pointer_cast<
                            SpatialDomains::RobinBoundaryCondition
                        >(m_PBndConds[n])->m_robinPrimitiveCoeff;

                    // checkout equation evaluation options!!
                    m_houtflow->m_pressurePrimCoeff[n] = coeff.Evaluate();

                    for (int i = 0; i < m_curl_dim; ++i)
                    {
                        Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                            UBndConds = m_fields[m_velocity[i]]->GetBndConditions();

                        LibUtilities::Equation coeff1 =
                            std::static_pointer_cast<
                                SpatialDomains::RobinBoundaryCondition
                            >(UBndConds[n])->m_robinPrimitiveCoeff;


                        m_houtflow->m_defVelPrimCoeff[i] = coeff1.GetExpression();

                        ASSERTL1(UBndConds[n]->GetBoundaryConditionType()
                                 == SpatialDomains::eRobin,"Require Velocity "
                                 "conditions to be of Robin type when pressure"
                                 "outflow is specticied as Robin Boundary type");


                        // checkout equation evaluation options!!
                        m_houtflow->m_velocityPrimCoeff[i][n] = coeff1.Evaluate();
                    }
                }
            }

        }
    }

    void Extrapolate::UpdateRobinPrimCoeff(void)
    {

        if((m_pressureCalls == 1) || (m_pressureCalls > m_intSteps))
        {
            return;
        }

        for(int n  = 0; n < m_PBndConds.size(); ++n)
        {
            // Get expansion with element on this boundary
            if(m_hbcType[n] == eConvectiveOBC)
            {
                for (int i = 0; i < m_curl_dim; ++i)
                {
                    Array<OneD, SpatialDomains::BoundaryConditionShPtr>
                        UBndConds = m_fields[m_velocity[i]]->UpdateBndConditions();

                    std::string primcoeff = m_houtflow->m_defVelPrimCoeff[i] + "*" +
                        boost::lexical_cast<std::string>(StifflyStable_Gamma0_Coeffs
                                                         [m_pressureCalls-1]);

                    SpatialDomains::RobinBCShPtr rcond =
                        std::dynamic_pointer_cast<
                            SpatialDomains::RobinBoundaryCondition >(UBndConds[n]);

                    SpatialDomains::BoundaryConditionShPtr bcond =
                        MemoryManager<SpatialDomains::RobinBoundaryCondition>::AllocateSharedPtr(
                                             m_session,rcond->m_robinFunction.GetExpression(),
                                             primcoeff,
                                             rcond->GetUserDefined(),
                                             rcond->m_filename);

                    UBndConds[n] = bcond;
                }

            }
        }
    }

    /**
     *
     */
    Array<OneD, NekDouble> Extrapolate::GetMaxStdVelocity(
        const Array<OneD, Array<OneD,NekDouble> > inarray)
    {
        // Checking if the problem is 2D
        ASSERTL0(m_curl_dim >= 2, "Method not implemented for 1D");

        int n_points_0      = m_fields[0]->GetExp(0)->GetTotPoints();
        int n_element       = m_fields[0]->GetExpSize();
        int nvel            = inarray.size();
        int cnt;

        NekDouble pntVelocity;

        // Getting the standard velocity vector
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(nvel);
        Array<OneD, NekDouble> tmp;
        Array<OneD, NekDouble> maxV(n_element, 0.0);
        LibUtilities::PointsKeyVector ptsKeys;

        for (int i = 0; i < nvel; ++i)
        {
            stdVelocity[i] = Array<OneD, NekDouble>(n_points_0);
        }

        cnt = 0.0;
        for (int el = 0; el < n_element; ++el)
        {
            int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
            ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();

            // reset local space
            if(n_points != n_points_0)
            {
                for (int j = 0; j < nvel; ++j)
                {
                    stdVelocity[j] = Array<OneD, NekDouble>(n_points, 0.0);
                }
                n_points_0 = n_points;
            }
            else
            {
                for (int j = 0; j < nvel; ++j)
                {
                    Vmath::Zero( n_points, stdVelocity[j], 1);
                }
            }

            Array<TwoD, const NekDouble> gmat =
                m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetDerivFactors(ptsKeys);

            if (m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetGtype()
                == SpatialDomains::eDeformed)
            {
                for(int j = 0; j < nvel; ++j)
                {
                    for(int k = 0; k < nvel; ++k)
                    {
                        Vmath::Vvtvp( n_points,  gmat[k*nvel + j], 1,
                                                 tmp = inarray[k] + cnt, 1,
                                                 stdVelocity[j], 1,
                                                 stdVelocity[j], 1);
                    }
                }
            }
            else
            {
                for(int j = 0; j < nvel; ++j)
                {
                    for(int k = 0; k < nvel; ++k)
                    {
                        Vmath::Svtvp( n_points,  gmat[k*nvel + j][0],
                                                 tmp = inarray[k] + cnt, 1,
                                                 stdVelocity[j], 1,
                                                 stdVelocity[j], 1);
                    }
                }
            }
            cnt += n_points;

            // Calculate total velocity in stdVelocity[0]
            Vmath::Vmul( n_points, stdVelocity[0], 1, stdVelocity[0], 1,
                                                      stdVelocity[0], 1);
            for(int k = 1; k < nvel; ++k)
            {
                Vmath::Vvtvp( n_points,  stdVelocity[k], 1,
                                         stdVelocity[k], 1,
                                         stdVelocity[0], 1,
                                         stdVelocity[0], 1);
            }
            pntVelocity = Vmath::Vmax( n_points, stdVelocity[0], 1);
            maxV[el] = sqrt(pntVelocity);
        }

        return maxV;
    }


    std::string Extrapolate::v_GetSubStepName(void)
    {
        return "";
    }

    /**
     *    At the start, the newest value is stored in array[nlevels-1]
     *        and the previous values in the first positions
     *    At the end, the extrapolated value is stored in array[nlevels-1]
     *        and the storage has been updated to included the new value
     */
    void Extrapolate::ExtrapolateArray(
            Array<OneD, Array<OneD, NekDouble> > &array)
    {
        int nint     = min(m_pressureCalls,m_intSteps);
        int nlevels  = array.size();
        int nPts     = array[0].size();

        // Update array
        RollOver(array);

        // Extrapolate to outarray
        Vmath::Smul(nPts, StifflyStable_Betaq_Coeffs[nint-1][nint-1],
                         array[nint-1],    1,
                         array[nlevels-1], 1);

        for(int n = 0; n < nint-1; ++n)
        {
            Vmath::Svtvp(nPts, StifflyStable_Betaq_Coeffs[nint-1][n],
                         array[n],1, array[nlevels-1],1,
                         array[nlevels-1],1);
        }
    }


    /**
     *    At the start, the newest value is stored in array[nlevels-1]
     *        and the previous values in the first positions
     *    At the end, the value of the bdf explicit part is stored in array[nlevels-1]
     *        and the storage has been updated to included the new value
     */
    void Extrapolate::EvaluateBDFArray(
            Array<OneD, Array<OneD, NekDouble> > &array)
    {
        int nint     = min(m_pressureCalls,m_intSteps);
        int nlevels  = array.size();
        int nPts     = array[0].size();

        // Update array
        RollOver(array);

        // Extrapolate to outarray
        Vmath::Smul(nPts, StifflyStable_Alpha_Coeffs[nint-1][nint-1],
                         array[nint-1],    1,
                         array[nlevels-1], 1);

        for(int n = 0; n < nint-1; ++n)
        {
            Vmath::Svtvp(nPts, StifflyStable_Alpha_Coeffs[nint-1][n],
                         array[n],1, array[nlevels-1],1,
                         array[nlevels-1],1);
        }
    }

    /**
     *    At the start, the newest value is stored in array[nlevels-1]
     *        and the previous values in the first positions
     *    At the end, the acceleration from BDF is stored in array[nlevels-1]
     *        and the storage has been updated to included the new value
     */
    void Extrapolate::AccelerationBDF(
            Array<OneD, Array<OneD, NekDouble> > &array)
    {
        int nlevels  = array.size();
        int nPts     = array[0].size();


        if(nPts)
        {
            // Update array
            RollOver(array);

            // Calculate acceleration using Backward Differentiation Formula
            Array<OneD, NekDouble> accelerationTerm (nPts, 0.0);
            if (m_pressureCalls > 2)
            {
                int acc_order = min(m_pressureCalls-2,m_intSteps);
                Vmath::Smul(nPts,
                            StifflyStable_Gamma0_Coeffs[acc_order-1],
                            array[0], 1,
                             accelerationTerm,  1);

                for(int i = 0; i < acc_order; i++)
                {
                    Vmath::Svtvp(nPts,
                                 -1*StifflyStable_Alpha_Coeffs[acc_order-1][i],
                                 array[i+1], 1,
                                 accelerationTerm,    1,
                                 accelerationTerm,    1);
                }
            }
            array[nlevels-1] = accelerationTerm;
        }
    }

    void Extrapolate::CopyPressureHBCsToPbndExp(void)
    {
        int n, cnt;
        for(cnt = n = 0; n < m_PBndConds.size(); ++n)
        {
            if((m_hbcType[n] == eHBCNeumann)||(m_hbcType[n] == eConvectiveOBC))
            {
                int nq = m_PBndExp[n]->GetNcoeffs();
                Vmath::Vcopy(nq, &(m_pressureHBCs[m_intSteps-1])[cnt],  1,
                                 &(m_PBndExp[n]->UpdateCoeffs()[0]), 1);
                cnt += nq;
            }
        }
    }
}
