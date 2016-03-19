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
// Description: Abstract base class for Extrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/Extrapolate.h>
#include <LibUtilities/Communication/Comm.h>

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
        typedef Loki::SingletonHolder<ExtrapolateFactory,
                                      Loki::CreateUsingNew,
                                      Loki::NoDestroy,
                                      Loki::SingleThreaded > Type;
        return Type::Instance();
    }

    Extrapolate::Extrapolate(
        const LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
        MultiRegions::ExpListSharedPtr pPressure,
        const Array<OneD, int> pVel,
        const SolverUtils::AdvectionSharedPtr advObject)
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
     * Function to extrapolate the new pressure boundary condition.
     * Based on the velocity field and on the advection term.
     * Acceleration term is also computed.
     * This routine is a general one for 2d and 3D application and it can be called
     * directly from velocity correction scheme. Specialisation on dimensionality is
     * redirected to the CalcPressureBCs method.
     */
    void Extrapolate::EvaluatePressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        const Array<OneD, const Array<OneD, NekDouble> >  &N,
        NekDouble kinvis)
    {
        m_pressureCalls++;
        if(m_HBCnumber > 0)
        {
            int  n,cnt;

            // Calculate Neumann BCs at current level
            CalcNeumannPressureBCs(fields, N, kinvis);

            //Calculate acceleration term at level n based on previous steps
            AccelerationBDF(m_acceleration);

            // Adding acceleration term to HOPBCs
            Vmath::Svtvp(m_numHBCDof, -1.0/m_timestep,
                         m_acceleration[m_intSteps],  1,
                         m_pressureHBCs[m_intSteps-1], 1,
                         m_pressureHBCs[m_intSteps-1], 1);

            // Extrapolate to n+1
            ExtrapolateArray(m_pressureHBCs);

            // Copy values of [dP/dn]^{n+1} in the pressure bcs storage.
            // m_pressureHBCS[nlevels-1] will be cancelled at next time step
            for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
            {
                if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
                {
                    int nq = m_PBndExp[n]->GetNcoeffs();
                    Vmath::Vcopy(nq, &(m_pressureHBCs[m_intSteps-1])[cnt],  1,
                                     &(m_PBndExp[n]->UpdateCoeffs()[0]), 1);
                    cnt += nq;
                }
            }

        }

        CalcOutflowBCs(fields, N, kinvis);
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
        Array<OneD, NekDouble> Uvals;

        Array<OneD, Array<OneD, NekDouble> > Velocity(m_curl_dim);
        Array<OneD, Array<OneD, NekDouble> > Advection(m_bnd_dim);

        Array<OneD, Array<OneD, NekDouble> > BndValues(m_bnd_dim);
        Array<OneD, Array<OneD, NekDouble> > Q(m_curl_dim);

        MultiRegions::ExpListSharedPtr BndElmtExp;
        for(n = cnt = 0; n < m_PBndConds.num_elements(); ++n)
        {
            // High order boundary condition;
            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
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
                    Q[i]         = Array<OneD, NekDouble> (nq,0.0);
                }

                // Obtaining fields on BndElmtExp
                for(int i = 0; i < m_curl_dim; i++)
                {
                    m_fields[0]->ExtractPhysToBndElmt(n, fields[i],Velocity[i]);
                }
                for(int i = 0; i < m_bnd_dim; i++)
                {
                    m_fields[0]->ExtractPhysToBndElmt(n, N[i],Advection[i]);
                }

                // CurlCurl
                BndElmtExp->CurlCurl(Velocity, Q);

                // Mounting advection component into the high-order condition
                for(int i = 0; i < m_bnd_dim; i++)
                {
                    MountHOPBCs(nq, kinvis,Q[i],Advection[i]);
                }

                Pvals = (m_pressureHBCs[m_intSteps-1]) + cnt;
                Uvals = (m_acceleration[m_intSteps]) + cnt;

                // Getting values on the edge and filling the pressure boundary
                // expansion and the acceleration term. Multiplication by the
                // normal is required
                for(int i = 0; i < m_bnd_dim; i++)
                {
                    m_fields[0]->ExtractElmtToBndPhys(n, Q[i],BndValues[i]);
                }
                m_PBndExp[n]->NormVectorIProductWRTBase(BndValues, Pvals);

                for(int i = 0; i < m_bnd_dim; i++)
                {
                    m_fields[0]->ExtractElmtToBndPhys(n, Velocity[i],BndValues[i]);
                }
                m_PBndExp[n]->NormVectorIProductWRTBase(BndValues, Uvals);

                // Get offset for next terms
                cnt += m_PBndExp[n]->GetNcoeffs();
            }
        }
    }

    void Extrapolate::CalcOutflowBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        const Array<OneD, const Array<OneD, NekDouble> >  &N,
        NekDouble kinvis)
    {
        if(m_outHBCnumber == 0)
        {
           return;
        }

        NekDouble U0,delta;
        m_session->LoadParameter("U0_HighOrderBC",U0,1.0);
        m_session->LoadParameter("Delta_HighOrderBC",delta,1/20.0);

        Array<OneD, Array<OneD, NekDouble> > Velocity(m_curl_dim);

        // Get velocity boundary conditions
        Array<OneD, Array<OneD, const
                    SpatialDomains::BoundaryConditionShPtr> >
                                            UBndConds(m_curl_dim);
        Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >
                                            UBndExp(m_curl_dim);

        for (int i = 0; i < m_curl_dim; ++i)
        {
            UBndConds[i] = m_fields[m_velocity[i]]->GetBndConditions();
            UBndExp[i]   = m_fields[m_velocity[i]]->GetBndCondExpansions();
        }

        MultiRegions::ExpListSharedPtr BndElmtExp;
        int cnt    = 0;
        for(int n  = 0; n < m_PBndConds.num_elements(); ++n)
        {
            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"HOutflow"))
            {
                // Get expansion with element on this boundary
                m_fields[0]->GetBndElmtExpansion(n, BndElmtExp, false);
                int nqb = m_PBndExp[n]->GetTotPoints();
                int nq  = BndElmtExp->GetTotPoints();

                // Get velocity and extrapolate
                for(int i = 0; i < m_curl_dim; i++)
                {
                    m_fields[0]->ExtractPhysToBndElmt(n, fields[i],
                                        m_outflowVel[cnt][i][m_intSteps-1]);
                    ExtrapolateArray(m_outflowVel[cnt][i]);
                    Velocity[i] = m_outflowVel[cnt][i][m_intSteps-1];
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

                // Calculate u.n and u^2
                Array<OneD, NekDouble> un (nqb, 0.0);
                Array<OneD, NekDouble> u2 (nqb, 0.0);
                for( int i = 0; i < m_curl_dim; i++)
                {
                    m_fields[0]->ExtractElmtToBndPhys(n, Velocity[i],bndVal);
                    Vmath::Vvtvp(nqb, normals[i], 1, bndVal, 1,
                                                 un, 1, un, 1);
                    Vmath::Vvtvp(nqb, bndVal, 1, bndVal, 1,
                                                 u2, 1, u2, 1);
                }
                // Calculate -1/2*u^2*S_0(u.n)
                for( int i = 0; i < nqb; i++)
                {
                    NekDouble fac = 0.5*(1.0-tanh(un[i]/(U0*delta)));
                    u2[i] = -0.5*u2[i]*fac;
                }

                // Calculate pressure boundary condition
                Array<OneD, NekDouble> pbc (nqb, 0.0);
                Vmath::Svtvp( nqb, kinvis, nGradUn, 1, u2, 1, pbc, 1);
                Vmath::Vsub( nqb, pbc, 1, m_PBndExp[n]->GetPhys(), 1,
                                             pbc, 1);
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
                // Calculate velocity boundary conditions
                Vmath::Vsub( nqb, pbc, 1, u2, 1, bndVal, 1);
                Vmath::Smul(nqb, kinvis, divU, 1, divU, 1);
                Vmath::Vsub( nqb, bndVal, 1, divU, 1, bndVal, 1);
                for( int i = 0; i < m_curl_dim; i++)
                {
                    if(boost::iequals(UBndConds[i][n]->GetUserDefined(),"HOutflow"))
                    {
                        // Reuse divU
                        Vmath::Vmul( nqb, normals[i], 1, bndVal, 1, divU, 1);
                        Vmath::Vadd( nqb, divU, 1, UBndExp[i][n]->GetPhys(), 1,
                                             divU, 1);
                        Vmath::Smul(nqb, 1.0/kinvis, divU, 1, divU, 1);
                        if ( m_fields[0]->GetWaveSpace())
                        {
                            UBndExp[i][n]->HomogeneousFwdTrans(divU, divU);
                        }
                        UBndExp[i][n]->IProductWRTBase(divU,
                                                UBndExp[i][n]->UpdateCoeffs());
                    }
                }

                // Get offset for next terms
                cnt++;
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
        int  nlevels = input.num_elements();
        
        Array<OneD, NekDouble> tmp;
    
        tmp = input[nlevels-1];
    
        for(int n = nlevels-1; n > 0; --n)
        {
            input[n] = input[n-1];
        }
    
        input[0] = tmp;
    }
    
    
    /**
     * Map to directly locate HOPBCs position and offsets in all scenarios
     */
    void Extrapolate::GenerateHOPBCMap()
    {
        m_PBndConds   = m_pressure->GetBndConditions();
        m_PBndExp     = m_pressure->GetBndCondExpansions();
    
        // Set up mapping from pressure boundary condition to pressure element
        // details.
        m_pressure->GetBoundaryToElmtMap(m_pressureBCtoElmtID,
                                         m_pressureBCtoTraceID);

        // find the maximum values of points  for pressure BC evaluation
        m_pressureBCsMaxPts = 0; 
        m_pressureBCsElmtMaxPts = 0; 
        int cnt, n;
        for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
        {
            for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i)
            {
                m_pressureBCsMaxPts = max(m_pressureBCsMaxPts,
                                m_PBndExp[n]->GetExp(i)->GetTotPoints());
                m_pressureBCsElmtMaxPts = max(m_pressureBCsElmtMaxPts,
                                m_pressure->GetExp(m_pressureBCtoElmtID[cnt++])
                                                            ->GetTotPoints());
            }
        }
    
        // Storage array for high order pressure BCs
        m_pressureHBCs = Array<OneD, Array<OneD, NekDouble> > (m_intSteps);
        m_acceleration = Array<OneD, Array<OneD, NekDouble> > (m_intSteps + 1);

        // Get number of expansions and points (or coeffs) for HOBCs
        m_HBCnumber = 0;
        m_numHBCDof = 0;
        m_outHBCnumber = 0;
        m_numOutHBCPts = 0;
        for( n = 0; n < m_PBndConds.num_elements(); ++n)
        {
            // High order boundary condition;
            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
            {
                m_numHBCDof += m_PBndExp[n]->GetNcoeffs();
                m_HBCnumber += m_PBndExp[n]->GetExpSize();
            }
            // High order outflow boundary condition;
            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"HOutflow"))
            {
                MultiRegions::ExpListSharedPtr BndElmtExp;
                m_pressure->GetBndElmtExpansion(n, BndElmtExp, false);
                m_numOutElmtPts += BndElmtExp->GetTotPoints();

                m_numOutHBCPts += m_PBndExp[n]->GetTotPoints();
                m_outHBCnumber += m_PBndExp[n]->GetTotPoints();;
            }
        }

        //int checkHBC = m_HBCnumber;
        //m_comm->AllReduce(checkHBC,LibUtilities::ReduceSum);
        //ASSERTL0(checkHBC > 0 ,"At least one high-order pressure boundary "
        //                       "condition is required for scheme "
        //                       "consistency");

        m_acceleration[0] = Array<OneD, NekDouble>(m_numHBCDof, 0.0);
        for(n = 0; n < m_intSteps; ++n)
        {
            m_pressureHBCs[n]   = Array<OneD, NekDouble>(m_numHBCDof, 0.0);
            m_acceleration[n+1] = Array<OneD, NekDouble>(m_numHBCDof, 0.0);
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
        if(m_numOutHBCPts > 0)
        {
            MultiRegions::ExpListSharedPtr BndElmtExp;
            m_outflowVel = Array<OneD,
                               Array<OneD,
                               Array<OneD,
                               Array<OneD, NekDouble> > > > (m_outHBCnumber);
            for(n = 0, cnt = 0; n < m_PBndConds.num_elements(); ++n)
            {
                if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"HOutflow"))
                {
                    m_outflowVel[cnt] = Array<OneD,
                                          Array<OneD,
                                          Array<OneD, NekDouble> > > (m_curl_dim);

                    m_fields[0]->GetBndElmtExpansion(n, BndElmtExp);
                    int nq  = BndElmtExp->GetTotPoints();
                    for(int j = 0; j < m_curl_dim; ++j)
                    {
                        m_outflowVel[cnt][j] = Array<OneD,
                                             Array<OneD, NekDouble> > (m_intSteps);
                        for(int k = 0; k < m_intSteps; ++k)
                        {
                            m_outflowVel[cnt][j][k] =
                                Array<OneD, NekDouble>(nq,0.0);
                        }
                    }
                    cnt++;
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
        int nvel            = inarray.num_elements();
        int cnt; 
        
        NekDouble pntVelocity;
        
        // Getting the standard velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(nvel);
        Array<OneD, NekDouble> maxV(n_element, 0.0);
        LibUtilities::PointsKeyVector ptsKeys;
        
        for (int i = 0; i < nvel; ++i)
        {
            stdVelocity[i] = Array<OneD, NekDouble>(n_points_0);
        }
        
        if (nvel == 2)
        {
            cnt = 0.0;
            for (int el = 0; el < n_element; ++el)
            { 
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();
                
                // reset local space if necessary
                if(n_points != n_points_0)
                {
                    for (int j = 0; j < nvel; ++j)
                    {
                        stdVelocity[j] = Array<OneD, NekDouble>(n_points);
                    }
                    n_points_0 = n_points;
                }        
                
                Array<TwoD, const NekDouble> gmat = 
                    m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetDerivFactors(ptsKeys);
                
                if (m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetGtype()
                    == SpatialDomains::eDeformed)
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i+cnt] 
                            + gmat[2][i]*inarray[1][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i+cnt] 
                            + gmat[3][i]*inarray[1][i+cnt];
                    }
                }
                else
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i+cnt] 
                            + gmat[2][0]*inarray[1][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i+cnt] 
                            + gmat[3][0]*inarray[1][i+cnt];
                    }
                }
                
                cnt += n_points;
                
                
                for (int i = 0; i < n_points; i++)
                {
                    pntVelocity = stdVelocity[0][i]*stdVelocity[0][i] 
                        + stdVelocity[1][i]*stdVelocity[1][i];
                    
                    if (pntVelocity>maxV[el])
                    {
                        maxV[el] = pntVelocity;
                    }
                }
                maxV[el] = sqrt(maxV[el]);
            }
        }
        else
        {
            cnt = 0;
            for (int el = 0; el < n_element; ++el)
            { 
                
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                ptsKeys = m_fields[0]->GetExp(el)->GetPointsKeys();
                
                // reset local space if necessary
                if(n_points != n_points_0)
                {
                    for (int j = 0; j < nvel; ++j)
                    {
                        stdVelocity[j] = Array<OneD, NekDouble>(n_points);
                    }
                    n_points_0 = n_points;
                }        
                
                Array<TwoD, const NekDouble> gmat =
                    m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetDerivFactors(ptsKeys);
                
                if (m_fields[0]->GetExp(el)->GetGeom()->GetMetricInfo()->GetGtype()
                    == SpatialDomains::eDeformed)
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i+cnt] 
                            + gmat[3][i]*inarray[1][i+cnt] 
                            + gmat[6][i]*inarray[2][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i+cnt] 
                            + gmat[4][i]*inarray[1][i+cnt] 
                            + gmat[7][i]*inarray[2][i+cnt];
                        
                        stdVelocity[2][i] = gmat[2][i]*inarray[0][i+cnt] 
                            + gmat[5][i]*inarray[1][i+cnt] 
                            + gmat[8][i]*inarray[2][i+cnt];
                    }
                }
                else
                {
                    for (int i = 0; i < n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i+cnt] 
                            + gmat[3][0]*inarray[1][i+cnt] 
                            + gmat[6][0]*inarray[2][i+cnt];
                        
                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i+cnt] 
                            + gmat[4][0]*inarray[1][i+cnt] 
                            + gmat[7][0]*inarray[2][i+cnt];
                        
                        stdVelocity[2][i] = gmat[2][0]*inarray[0][i+cnt] 
                            + gmat[5][0]*inarray[1][i+cnt] 
                            + gmat[8][0]*inarray[2][i+cnt];
                    }
                }
                
                cnt += n_points;
                
                for (int i = 0; i < n_points; i++)
                {
                    pntVelocity = stdVelocity[0][i]*stdVelocity[0][i] 
                        + stdVelocity[1][i]*stdVelocity[1][i] 
                        + stdVelocity[2][i]*stdVelocity[2][i];
                    
                    if (pntVelocity > maxV[el])
                    {
                        maxV[el] = pntVelocity;
                    }
                }

                maxV[el] = sqrt(maxV[el]);
                //cout << maxV[el]*maxV[el] << endl;
            }
        }
        
        return maxV;
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
        int nlevels  = array.num_elements();
        int nPts     = array[0].num_elements();

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
     *    At the end, the acceleration from BDF is stored in array[nlevels-1]
     *        and the storage has been updated to included the new value
     */
    void Extrapolate::AccelerationBDF(
            Array<OneD, Array<OneD, NekDouble> > &array)
    {
        int nlevels  = array.num_elements();
        int nPts     = array[0].num_elements();

        // Update array
        RollOver(array);

        // Calculate acceleration using Backward Differentiation Formula
        Array<OneD, NekDouble> accelerationTerm (nPts, 0.0);
        if (m_pressureCalls > 2)
        {
            int acc_order = min(m_pressureCalls-2,m_intSteps);
            Vmath::Smul(nPts,
                             StifflyStable_Gamma0_Coeffs[acc_order-1],
                             m_acceleration[0], 1,
                             accelerationTerm,  1);

            for(int i = 0; i < acc_order; i++)
            {
                Vmath::Svtvp(nPts,
                            -1*StifflyStable_Alpha_Coeffs[acc_order-1][i],
                             m_acceleration[i+1], 1,
                             accelerationTerm,    1,
                             accelerationTerm,    1);
            }
        }
        m_acceleration[nlevels-1] = accelerationTerm;
    }

}

