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
        if(m_HBCdata.num_elements()>0)
        {
            Array<OneD, NekDouble> tmp;
            Array<OneD, NekDouble> accelerationTerm;

            m_pressureCalls++;
            int  n,cnt;
            int  nint    = min(m_pressureCalls,m_intSteps);
            int  nlevels = m_pressureHBCs.num_elements();

            int acc_order = 0;

            accelerationTerm =
                Array<OneD, NekDouble>(m_acceleration[0].num_elements(), 0.0);

            // Rotate HOPBCs storage
            RollOver(m_pressureHBCs);

            // Rotate acceleration term
            RollOver(m_acceleration);

            // Calculate Neumann BCs at current level
            CalcNeumannPressureBCs(fields, N, kinvis);

            // Copy High order values into storage array 
            for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
            {
                // High order boundary condition;
                if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
                {
                    int nq = m_PBndExp[n]->GetNcoeffs();
                    Vmath::Vcopy(nq, &(m_PBndExp[n]->GetCoeffs()[0]), 1,
                                     &(m_pressureHBCs[0])[cnt], 1);
                    cnt += nq;
                }
            }

            //Calculate acceleration term at level n based on previous steps
            if (m_pressureCalls > 2)
            {
                acc_order = min(m_pressureCalls-2,m_intSteps);
                Vmath::Smul(cnt, StifflyStable_Gamma0_Coeffs[acc_order-1],
                                 m_acceleration[0], 1,
                                 accelerationTerm,  1);

                for(int i = 0; i < acc_order; i++)
                {
                    Vmath::Svtvp(cnt, 
                                -1*StifflyStable_Alpha_Coeffs[acc_order-1][i],
                                 m_acceleration[i+1], 1,
                                 accelerationTerm,    1,
                                 accelerationTerm,    1);
                }
            }

            // Adding acceleration term to HOPBCs
            Vmath::Svtvp(cnt, -1.0/m_timestep,
                         accelerationTerm,  1,
                         m_pressureHBCs[0], 1,
                         m_pressureHBCs[0], 1);


            // Extrapolate to n+1
            Vmath::Smul(cnt, StifflyStable_Betaq_Coeffs[nint-1][nint-1],
                             m_pressureHBCs[nint-1],    1,
                             m_pressureHBCs[nlevels-1], 1);

            for(n = 0; n < nint-1; ++n)
            {
                Vmath::Svtvp(cnt, StifflyStable_Betaq_Coeffs[nint-1][n],
                             m_pressureHBCs[n],1,m_pressureHBCs[nlevels-1],1,
                             m_pressureHBCs[nlevels-1],1);
            }

            // Copy values of [dP/dn]^{n+1} in the pressure bcs storage.
            // m_pressureHBCS[nlevels-1] will be cancelled at next time step
            for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
            {
                if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
                {
                    int nq = m_PBndExp[n]->GetNcoeffs();
                    Vmath::Vcopy(nq, &(m_pressureHBCs[nlevels-1])[cnt],  1,
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
    void Extrapolate::CalcNeumannPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        const Array<OneD, const Array<OneD, NekDouble> >  &N,
        NekDouble kinvis)
    {    
        Array<OneD, NekDouble> Pvals;
        Array<OneD, NekDouble> Uvals;
        StdRegions::StdExpansionSharedPtr Pbc;
        StdRegions::StdExpansionSharedPtr elmt;

        Array<OneD, Array<OneD, const NekDouble> > Velocity(m_curl_dim);
        Array<OneD, Array<OneD, const NekDouble> > Advection(m_bnd_dim);
        
        Array<OneD, Array<OneD, NekDouble> > BndValues(m_bnd_dim);
        Array<OneD, Array<OneD, NekDouble> > Q(m_bnd_dim);
        
        for(int i = 0; i < m_bnd_dim; i++)
        {
            BndValues[i] = Array<OneD, NekDouble> (m_pressureBCsMaxPts,0.0);
            Q[i]         = Array<OneD, NekDouble> (m_pressureBCsElmtMaxPts,0.0);
        }
        
        for(int j = 0 ; j < m_HBCdata.num_elements() ; j++)
        {
            /// Casting the boundary expansion to the specific case
            Pbc = boost::dynamic_pointer_cast<StdRegions::StdExpansion> 
                        (m_PBndExp[m_HBCdata[j].m_bndryElmtID]
                            ->GetExp(m_HBCdata[j].m_bndElmtOffset));

            /// Picking up the element where the HOPBc is located
            elmt = m_pressure->GetExp(m_HBCdata[j].m_globalElmtID);

            /// Assigning
            for(int i = 0; i < m_bnd_dim; i++)
            {
                Velocity[i]  = fields[i] + m_HBCdata[j].m_physOffset;
                Advection[i] = N[i]      + m_HBCdata[j].m_physOffset;
            }

            // for the 3DH1D case we need to grab the conjugate mode
            if(m_pressure->GetExpType() == MultiRegions::e3DH1D)
            {
                Velocity[2]  = fields[2] + m_HBCdata[j].m_assPhysOffset;
            }

            /// Calculating the curl-curl and storing it in Q
            CurlCurl(Velocity,Q,j);

            // Mounting advection component into the high-order condition
            for(int i = 0; i < m_bnd_dim; i++)
            {
                MountHOPBCs(m_HBCdata[j].m_ptsInElmt,kinvis,Q[i],Advection[i]);
            }

            Pvals = m_PBndExp[m_HBCdata[j].m_bndryElmtID]->UpdateCoeffs()
                        + m_PBndExp[m_HBCdata[j].m_bndryElmtID]
                            ->GetCoeff_Offset(m_HBCdata[j].m_bndElmtOffset);
            Uvals = (m_acceleration[0]) + m_HBCdata[j].m_coeffOffset;

            // Getting values on the edge and filling the pressure boundary
            // expansion and the acceleration term. Multiplication by the
            // normal is required
            switch(m_pressure->GetExpType())
            {
                case MultiRegions::e2D:
                case MultiRegions::e3DH1D:
                {
                    elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                          Q[0], BndValues[0]);
                    elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                          Q[1], BndValues[1]);
                    Pbc->NormVectorIProductWRTBase(BndValues[0], BndValues[1],
                                                   Pvals);

                    elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                          Velocity[0], BndValues[0]);
                    elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                          Velocity[1], BndValues[1]);
                    Pbc->NormVectorIProductWRTBase(BndValues[0], BndValues[1],
                                                   Uvals);
                }
                break;
                case MultiRegions::e3DH2D:
                {
                    if(m_HBCdata[j].m_elmtTraceID == 0)
                    {
                        (m_PBndExp[m_HBCdata[j].m_bndryElmtID]->UpdateCoeffs()
                            + m_PBndExp[m_HBCdata[j].m_bndryElmtID]
                                ->GetCoeff_Offset(
                                    m_HBCdata[j].m_bndElmtOffset))[0]
                                                                = -1.0*Q[0][0];
                    }
                    else if (m_HBCdata[j].m_elmtTraceID == 1)
                    {
                        (m_PBndExp[m_HBCdata[j].m_bndryElmtID]->UpdateCoeffs()
                            + m_PBndExp[m_HBCdata[j].m_bndryElmtID]
                                ->GetCoeff_Offset(
                                    m_HBCdata[j].m_bndElmtOffset))[0] 
                                            = Q[0][m_HBCdata[j].m_ptsInElmt-1];
                    }
                    else
                    {
                        ASSERTL0(false,
                                 "In the 3D homogeneous 2D approach BCs edge "
                                 "ID can be just 0 or 1 ");
                    }
                }
                break;
                case MultiRegions::e3D:
                {
                    elmt->GetFacePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc, 
                                          Q[0], BndValues[0]);
                    elmt->GetFacePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                          Q[1], BndValues[1]);
                    elmt->GetFacePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                          Q[2], BndValues[2]);
                    Pbc->NormVectorIProductWRTBase(BndValues[0], BndValues[1],
                                          BndValues[2], Pvals);

                    elmt->GetFacePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                          Velocity[0], BndValues[0]);
                    elmt->GetFacePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                          Velocity[1], BndValues[1]);
                    elmt->GetFacePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                          Velocity[2], BndValues[2]);
                    Pbc->NormVectorIProductWRTBase(BndValues[0], BndValues[1],
                                          BndValues[2], Uvals);
                }
                break;
            default:
                ASSERTL0(0,"Dimension not supported");
                break;
            }
        }
    }



    void Extrapolate::CalcOutflowBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        const Array<OneD, const Array<OneD, NekDouble> >  &N,
        NekDouble kinvis)
    {
        static bool init = true;
        static bool noHOBC = false;

        if(noHOBC == true)
        {
           return;
        }
        
        if(init) // set up storage for boundary velocity at outflow
        {
            init = false;
            int totbndpts = 0;
            for(int n = 0; n < m_PBndConds.num_elements(); ++n)
            {
                if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"HOutflow"))
                {
                    totbndpts += m_PBndExp[n]->GetTotPoints();
                }
            }
          
            if(totbndpts == 0)
            {
                noHOBC = true;
                return;
            }
  
            m_outflowVel = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (m_curl_dim);
            for(int i = 0; i < m_curl_dim; ++i)
            {
                m_outflowVel[i] = Array<OneD, Array<OneD, NekDouble> >(m_curl_dim);
                for(int j = 0; j < m_curl_dim; ++j)
                {             
                    // currently just set up for 2nd order extrapolation 
                    m_outflowVel[i][j] = Array<OneD, NekDouble>(totbndpts,0.0);
                }
            }

            if (m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
            {
                m_PhyoutfVel = Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (m_curl_dim);

                for(int i = 0; i < m_curl_dim; ++i)
                {
                    m_PhyoutfVel[i] = Array<OneD, Array<OneD, NekDouble> > (m_curl_dim);
                    for(int j = 0; j < m_curl_dim; ++j)
                    {
                        // currently just set up for 2nd order extrapolation
                        m_PhyoutfVel[i][j] = Array<OneD, NekDouble> (totbndpts,0.0);
                    }
                }

                m_nonlinearterm_phys   = Array<OneD, NekDouble> (totbndpts,0.0);
                m_nonlinearterm_coeffs = Array<OneD, NekDouble> (totbndpts,0.0);

                m_PBndCoeffs = Array<OneD, NekDouble> (totbndpts,0.0);
                m_UBndCoeffs = Array<OneD, Array<OneD, NekDouble> > (m_curl_dim);
                for(int i = 0; i < m_curl_dim; ++i)
                {
                    m_UBndCoeffs[i] = Array<OneD, NekDouble> (totbndpts);   
                }
                Array<OneD, unsigned int> planes;
                planes = m_pressure->GetZIDs();
                int num_planes = planes.num_elements();
                m_expsize_per_plane = Array<OneD, unsigned int> (m_PBndConds.num_elements());
                for(int n = 0; n < m_PBndConds.num_elements(); ++n)
                {
                    int exp_size = m_PBndExp[n]->GetExpSize();
                    m_expsize_per_plane[n] = exp_size/num_planes;
                }
                m_totexps_per_plane = 0;
                for(int n = 0; n < m_PBndConds.num_elements(); ++n)
                {
                    m_totexps_per_plane += m_PBndExp[n]->GetExpSize()/num_planes;
                }
            }
        }
        
        StdRegions::StdExpansionSharedPtr Bc,Pbc; 
        Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr> >
                                                        UBndConds(m_curl_dim);
        Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >
                                                        UBndExp(m_curl_dim);

        for (int i = 0; i < m_curl_dim; ++i)
        {
            UBndConds[i] = m_fields[m_velocity[i]]->GetBndConditions();
            UBndExp[i]   = m_fields[m_velocity[i]]->GetBndCondExpansions();
        }

        Array<OneD, Array<OneD, NekDouble> > BndValues(m_curl_dim);
        Array<OneD, Array<OneD, NekDouble> > BndElmt  (m_curl_dim);
        Array<OneD, Array<OneD, NekDouble> > nGradu   (m_curl_dim);
        Array<OneD, NekDouble > gradtmp (m_pressureBCsElmtMaxPts),
                                fgradtmp(m_pressureBCsElmtMaxPts);
        
        nGradu[0] = Array<OneD, NekDouble>(m_curl_dim*m_pressureBCsMaxPts);
        for(int i = 0; i < m_curl_dim; ++i)
        {
            BndElmt[i]   = Array<OneD, NekDouble> (m_pressureBCsElmtMaxPts,
                                                   0.0);
            BndValues[i] = Array<OneD, NekDouble> (m_pressureBCsMaxPts,0.0);
            nGradu[i]  = nGradu[0] + i*m_pressureBCsMaxPts;
            RollOver(m_outflowVel[i]);

            if (m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
            {
                RollOver(m_PhyoutfVel[i]);
            }
        }
            
        int nbc,cnt,cnt_start;
        int veloffset = 0;
        int  nint    = min(m_pressureCalls,m_intSteps);

        StdRegions::StdExpansionSharedPtr elmt;
        Array<OneD, NekDouble> PBCvals, UBCvals;
        Array<OneD, Array<OneD, NekDouble> > ubc(m_curl_dim);
        Array<OneD, Array<OneD, NekDouble> > normals;

        cnt = 0;
        for(int n = 0; n < m_PBndConds.num_elements(); ++n)
        {
            // Do outflow boundary conditions if they exist
            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"HOutflow"))
            {
                for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i,cnt++)
                {
                    cnt = max(cnt,m_PBndExp[n]->GetTotPoints());
                }
            }
        }

        for(int i =0; i < m_curl_dim; ++i)
        {
            ubc[i] = Array<OneD, NekDouble>(cnt);
        }

        NekDouble U0,delta;
        m_session->LoadParameter("U0_HighOrderBC",U0,1.0);
        m_session->LoadParameter("Delta_HighOrderBC",delta,1/20.0);

        cnt = 0;
        int count = 0;
        for(int n = 0; n < m_PBndConds.num_elements(); ++n)
        {
            cnt_start = cnt;

            // Do outflow boundary conditions if they exist
            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"HOutflow"))
            {

                if (m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
                {
                    int cnt_exp   = 0; int cnt_plane = 0;
                    int veloffset = 0;
                    for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i, cnt_exp++)
                    {
                        // count the expansion list in each plane for e3DH1D case
                        if(cnt_exp == m_expsize_per_plane[n])
                        {
                            cnt_exp = 0; cnt_plane++;
                        }
                        int cnt = cnt_plane * m_totexps_per_plane + cnt_exp + count;

                        // find element and edge of this expansion.
                        Bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion>
                            (m_PBndExp[n]->GetExp(i));

                        int elmtid = m_pressureBCtoElmtID[cnt];
                        elmt       = m_fields[0]->GetExp(elmtid);
                        int offset = m_fields[0]->GetPhys_Offset(elmtid);

                        int boundary = m_pressureBCtoTraceID[cnt];

                        // Determine extrapolated U,V values
                        int nq = elmt->GetTotPoints();
                        int nbc = m_PBndExp[n]->GetExp(i)->GetTotPoints();
                        // currently just using first order approximation here.
                        // previously have obtained value from m_integrationSoln
                        Array<OneD, NekDouble> veltmp;

                        for(int j = 0; j < m_curl_dim; ++j)
                        {
                            Vmath::Vcopy(nq, &fields[m_velocity[j]][offset], 1,
                                         &BndElmt[j][0],                 1);
                            elmt->GetTracePhysVals(boundary,Bc,BndElmt[j],
                                         veltmp = m_outflowVel[j][0] + veloffset);
                        }
                        veloffset += nbc;
                    }

                    // for velocity on the outflow boundary in e3DH1D,
                    // we need to make a backward fourier transformation
                    // to get the physical coeffs at the outflow BCs.
                    for(int j = 0; j < m_curl_dim; ++j)
                    {
                        m_PBndExp[n]->HomogeneousBwdTrans(
                                m_outflowVel[j][0],
                                m_PhyoutfVel[j][0]);
                    }

                    cnt_plane = 0; cnt_exp = 0;
                    veloffset = 0;
                    for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i, cnt_exp++)
                    {
                        // count the expansion list for each plane for e3DH1D
                        if(cnt_exp == m_expsize_per_plane[n])
                        {
                            cnt_exp = 0; cnt_plane++;
                        }
                        cnt = cnt_plane * m_totexps_per_plane + cnt_exp + count;

                        int elmtid = m_pressureBCtoElmtID[cnt];
                        elmt       = m_fields[0]->GetExp(elmtid);
                        int nbc = m_PBndExp[n]->GetExp(i)->GetTotPoints();

                        Array<OneD, NekDouble>  veltmp(nbc,0.0),
                                                normDotu(nbc,0.0), utot(nbc,0.0);
                        int boundary = m_pressureBCtoTraceID[cnt];
                        normals=elmt->GetSurfaceNormal(boundary);

                        // extrapolate velocity
                        if(nint <= 1)
                        {
                            for(int j = 0; j < m_curl_dim; ++j)
                            {
                                Vmath::Vcopy(nbc,
                                        veltmp = m_PhyoutfVel[j][0] +veloffset, 1,
                                        BndValues[j],                           1);
                            }
                        }
                        else // only set up for 2nd order extrapolation
                        {
                            for(int j = 0; j < m_curl_dim; ++j)
                            {
                                Vmath::Smul(nbc, 2.0,
                                        veltmp = m_PhyoutfVel[j][0] + veloffset, 1,
                                        BndValues[j],                            1);
                                Vmath::Svtvp(nbc, -1.0,
                                        veltmp = m_PhyoutfVel[j][1] + veloffset, 1,
                                        BndValues[j],                            1,
                                        BndValues[j],                            1);
                            }    
                        }

                        // Set up |u|^2, n.u in physical space
                        for(int j = 0; j < m_curl_dim; ++j)
                        {
                            Vmath::Vvtvp(nbc, BndValues[j], 1, BndValues[j], 1,
                                              utot,         1, utot,         1);
                        }
                        for(int j = 0; j < m_bnd_dim; ++j)
                        {
                            Vmath::Vvtvp(nbc, normals[j],   1, BndValues[j], 1,
                                              normDotu,     1, normDotu,     1);
                        }
                        
                        int Offset = m_PBndExp[n]->GetPhys_Offset(i);

                        for(int k = 0; k < nbc; ++k)
                        {
                            // calculate the nonlinear term (kinetic energy
                            //  multiplies step function) in physical space
                            NekDouble fac = 0.5*(1.0-tanh(normDotu[k]/(U0*delta)));
                            m_nonlinearterm_phys[k + Offset] =  0.5 * utot[k] * fac;
                        }

                        veloffset += nbc;
                    }

                    // for e3DH1D, we need to make a forward fourier transformation
                    // for the kinetic energy term (nonlinear)
                    UBndExp[0][n]->HomogeneousFwdTrans(
                            m_nonlinearterm_phys,
                            m_nonlinearterm_coeffs);

                    // for e3DH1D, we need to make a forward fourier transformation
                    // for Dirichlet pressure boundary condition that is from input file
                    m_PBndExp[n]->HomogeneousFwdTrans(
                            m_PBndExp[n]->UpdatePhys(),
                            m_PBndCoeffs);
                    // for e3DH1D, we need to make a forward fourier transformation
                    // for Neumann velocity boundary condition that is from input file
                    for (int j = 0; j < m_curl_dim; ++j)
                    {
                        UBndExp[j][n]->HomogeneousFwdTrans(
                            UBndExp[j][n]->UpdatePhys(),
                            m_UBndCoeffs[j]);
                    }
                }

                veloffset = 0;
                int cnt_exp = 0; int cnt_plane = 0; //only useful in e3DH1D case
                for(int i = 0; i < m_PBndExp[n]->GetExpSize(); ++i,cnt++)
                {
                    if (m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
                    {
                        // count the expansion list for e3DH1D
                        if(cnt_exp == m_expsize_per_plane[n])
                        {
                            cnt_exp = 0; cnt_plane++;
                        }
                        cnt = cnt_plane * m_totexps_per_plane + cnt_exp + count;
                        cnt_exp++;
                    }

                    // find element and edge of this expansion. 
                    Bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion> 
                        (m_PBndExp[n]->GetExp(i));

                    int elmtid = m_pressureBCtoElmtID[cnt];
                    elmt       = m_fields[0]->GetExp(elmtid);
                    int offset = m_fields[0]->GetPhys_Offset(elmtid);

                    // Determine extrapolated U,V values
                    int nq = elmt->GetTotPoints();

                    // currently just using first order approximation here. 
                    // previously have obtained value from m_integrationSoln
                    for(int j = 0; j < m_bnd_dim; ++j)
                    {
                        Vmath::Vcopy(nq, &fields[m_velocity[j]][offset], 1,
                                         &BndElmt[j][0],                 1);
                    }

                    int nbc      = m_PBndExp[n]->GetExp(i)->GetTotPoints();
                    int boundary = m_pressureBCtoTraceID[cnt];

                    Array<OneD, NekDouble>  ptmp(nbc,0.0),
                                            divU(nbc,0.0);

                    normals=elmt->GetSurfaceNormal(boundary);
                    Vmath::Zero(m_bnd_dim*m_pressureBCsMaxPts,nGradu[0],1);

                    for (int j = 0; j < m_bnd_dim; j++)
                    {
                        // Calculate Grad u =  du/dx, du/dy, du/dz, etc. 
                        for (int k = 0; k< m_bnd_dim; k++)
                        {
                            elmt->PhysDeriv(MultiRegions::DirCartesianMap[k],
                                            BndElmt[j], gradtmp);
                            elmt->GetTracePhysVals(boundary, Bc, gradtmp,
                                                   fgradtmp);
                            Vmath::Vvtvp(nbc,normals[k], 1, fgradtmp, 1,
                                             nGradu[j],  1, nGradu[j],1);
                            if(j == k)
                            {
                                Vmath::Vadd(nbc,fgradtmp, 1, divU, 1, divU, 1);
                            }
                        }
                    }

                    if (m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
                    {
                        // Set up |u|^2, n.u, div(u), and (n.grad(u) . n) for
                        // pressure condition
                        for(int j = 0; j < m_bnd_dim; ++j)
                        {
                            Vmath::Vvtvp(nbc, normals[j],   1, nGradu[j],    1,
                                              ptmp,         1, ptmp,         1);
                        }
                        int p_offset = m_PBndExp[n]->GetPhys_Offset(i);

                        for(int k = 0; k < nbc; ++k)
                        {
                            // Set up Dirichlet pressure condition and
                            // store in ptmp (m_UBndCoeffs contains Fourier Coeffs of the
                            // function from the input file )

                            ptmp[k] =  kinvis * ptmp[k] 
                                        - m_nonlinearterm_coeffs[k + p_offset]
                                                  - m_PBndCoeffs[k + p_offset];
                        }

                        int u_offset = UBndExp[0][n]->GetPhys_Offset(i);

                        for(int j = 0; j < m_bnd_dim; ++j)
                        {
                            for(int k = 0; k < nbc; ++k)
                            {
                                ubc[j][k + u_offset] = (1.0 / kinvis)
                                                * (m_UBndCoeffs[j][k + u_offset]
                                          + m_nonlinearterm_coeffs[k + u_offset]
                                                                * normals[j][k]);
                            }
                        }

                        // boundary condition for velocity in homogenous direction
                        for(int k = 0; k < nbc; ++k)
                        {
                            ubc[m_bnd_dim][k + u_offset] = (1.0 / kinvis)
                                                * m_UBndCoeffs[m_bnd_dim][k + u_offset];
                        }

                        u_offset = UBndExp[m_bnd_dim][n]->GetPhys_Offset(i);
                        UBCvals  = UBndExp[m_bnd_dim][n]->UpdateCoeffs()
                                    + UBndExp[m_bnd_dim][n]->GetCoeff_Offset(i);
                        Bc->IProductWRTBase(ubc[m_bnd_dim] + u_offset, UBCvals);
                    }
                    else
                    {

                        Array<OneD, NekDouble>  veltmp, utot(nbc,0.0),
                                                normDotu(nbc,0.0);
                        // extract velocity and store
                        for(int j = 0; j < m_bnd_dim; ++j)
                        {
                            elmt->GetTracePhysVals(boundary,Bc,BndElmt[j],
                                           veltmp = m_outflowVel[j][0] + veloffset);
                        }

                        // extrapolate velocity
                        if(nint <= 1)
                        {
                            for(int j = 0; j < m_bnd_dim; ++j)
                            {
                                Vmath::Vcopy(nbc,
                                        veltmp = m_outflowVel[j][0] +veloffset, 1,
                                        BndValues[j],                           1);
                            }
                        }
                        else // only set up for 2nd order extrapolation
                        {
                            for(int j = 0; j < m_bnd_dim; ++j)
                            {
                                Vmath::Smul(nbc, 2.0,
                                        veltmp = m_outflowVel[j][0] + veloffset, 1,
                                        BndValues[j],                            1);
                                Vmath::Svtvp(nbc, -1.0,
                                        veltmp = m_outflowVel[j][1] + veloffset, 1,
                                        BndValues[j],                            1,
                                        BndValues[j],                            1);
                            }
                        }

                        // Set up |u|^2, n.u, div(u), and (n.grad(u) . n) for
                        // pressure condition
                        for(int j = 0; j < m_bnd_dim; ++j)
                        {
                            Vmath::Vvtvp(nbc, BndValues[j], 1, BndValues[j], 1,
                                              utot,         1, utot,         1);
                            Vmath::Vvtvp(nbc, normals[j],   1, BndValues[j], 1,
                                              normDotu,     1, normDotu,     1);
                            Vmath::Vvtvp(nbc, normals[j],   1, nGradu[j],    1,
                                              ptmp,         1, ptmp,         1);
                        }

                        PBCvals = m_PBndExp[n]->GetPhys() +
                                        m_PBndExp[n]->GetPhys_Offset(i);

                        for(int k = 0; k < nbc; ++k)
                        {
                            NekDouble fac = 0.5*(1.0-tanh(normDotu[k]/(U0*delta)));

                            // Set up Dirichlet pressure condition and
                            // store in ptmp (PBCvals contains a
                            // function from the input file )
                            ptmp[k] =  kinvis * ptmp[k] - 0.5 * utot[k] * fac
                                                        - PBCvals[k];
                        }

                        int u_offset = UBndExp[0][n]->GetPhys_Offset(i);

                        for(int j = 0; j < m_bnd_dim; ++j)
                        {
                            UBCvals = UBndExp[j][n]->GetPhys()
                                        + UBndExp[j][n]->GetPhys_Offset(i);

                            for(int k = 0; k < nbc; ++k)
                            {
                                NekDouble fac        = 0.5 * (1.0 - tanh(normDotu[k]
                                                                / (U0 * delta)));
                                ubc[j][k + u_offset] = (1.0 / kinvis)
                                                * (UBCvals[k] + 0.5 * utot[k] * fac
                                                                * normals[j][k]);
                            }
                        }
                    }

                    // set up pressure boundary condition
                    PBCvals = m_PBndExp[n]->UpdateCoeffs()
                            + m_PBndExp[n]->GetCoeff_Offset(i);
                    m_PBndExp[n]->GetExp(i)->FwdTrans(ptmp,PBCvals);

                    veloffset += nbc;
                }

                // Now set up Velocity conditions.
                for(int j = 0; j < m_bnd_dim; j++)
                {
                    if(boost::iequals(UBndConds[j][n]->GetUserDefined(),"HOutflow"))
                    {
                        cnt = cnt_start;

                        int cnt_exp = 0; int cnt_plane = 0; //only useful in e3DH1D case    
                        for(int i = 0; i < UBndExp[0][n]->GetExpSize();
                                       ++i, cnt++)
                        {
                            if(m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
                            {
                            
                                // count the expansion order for e3DH1D
                                if(cnt_exp == m_expsize_per_plane[n])
                                {
                                    cnt_exp = 0; cnt_plane++;
                                }
                                cnt = cnt_plane * m_totexps_per_plane + cnt_exp + count;
                                cnt_exp++;
                            }

                            Pbc =  StdRegions::StdExpansionSharedPtr
                                            (m_PBndExp[n]->GetExp(i));
                            Bc  =  StdRegions::StdExpansionSharedPtr
                                            (UBndExp[0][n]->GetExp(i));

                            nbc = UBndExp[0][n]->GetExp(i)->GetTotPoints();
                            int boundary = m_pressureBCtoTraceID[cnt];

                            Array<OneD, NekDouble> pb(nbc), ub(nbc);
                            int elmtid = m_pressureBCtoElmtID[cnt];

                            elmt   = m_fields[0]->GetExp(elmtid);

                            normals = elmt->GetSurfaceNormal(boundary);

                            // Get p from projected boundary condition
                            PBCvals = m_PBndExp[n]->UpdateCoeffs()
                                    + m_PBndExp[n]->GetCoeff_Offset(i);
                            Pbc->BwdTrans(PBCvals,pb);

                            int u_offset = UBndExp[j][n]->GetPhys_Offset(i);

                            for(int k = 0; k < nbc; ++k)
                            {
                                ub[k] = ubc[j][k + u_offset]
                                      + pb[k] * normals[j][k] / kinvis;
                            }

                            UBCvals = UBndExp[j][n]->UpdateCoeffs()
                                    + UBndExp[j][n]->GetCoeff_Offset(i);
                            Bc->IProductWRTBase(ub,UBCvals);
                        }
                    }
                }
            }
            else
            {
                cnt += m_PBndExp[n]->GetExpSize();
                if(m_fields[0]->GetExpType() == MultiRegions::e3DH1D)
                {
                    count  += m_expsize_per_plane[n];
                }
            }
        }
    }
    
    
    /**
     * Curl Curl routine - dimension dependent
     */
    void Extrapolate::CurlCurl(
        Array<OneD, Array<OneD, const NekDouble> > &Vel,
        Array<OneD, Array<OneD, NekDouble> > &Q,
        const int j)
    {
        StdRegions::StdExpansionSharedPtr elmt
                        = m_fields[0]->GetExp(m_HBCdata[j].m_globalElmtID);

        Array<OneD,NekDouble> Vx(m_pressureBCsElmtMaxPts);
        Array<OneD,NekDouble> Uy(m_pressureBCsElmtMaxPts);

        switch(m_fields[0]->GetExpType())
        {
            case MultiRegions::e2D:
            {
                Array<OneD,NekDouble> Dummy(m_pressureBCsElmtMaxPts);

                elmt->PhysDeriv(MultiRegions::DirCartesianMap[0], Vel[1], Vx);
                elmt->PhysDeriv(MultiRegions::DirCartesianMap[1], Vel[0], Uy);

                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Vx, 1, Uy, 1, Dummy, 1);

                elmt->PhysDeriv(Dummy,Q[1],Q[0]);

                Vmath::Smul(m_HBCdata[j].m_ptsInElmt, -1.0, Q[1], 1, Q[1], 1);
            }
            break;

            case MultiRegions::e3DH1D:
            {
                Array<OneD,NekDouble> Wz(m_pressureBCsElmtMaxPts);

                Array<OneD,NekDouble> Dummy1(m_pressureBCsElmtMaxPts);
                Array<OneD,NekDouble> Dummy2(m_pressureBCsElmtMaxPts);

                elmt->PhysDeriv(MultiRegions::DirCartesianMap[0], Vel[1], Vx);
                elmt->PhysDeriv(MultiRegions::DirCartesianMap[1], Vel[0], Uy);
                Vmath::Smul(m_HBCdata[j].m_ptsInElmt, m_wavenumber[j],
                                                      Vel[2], 1, Wz,     1);

                elmt->PhysDeriv(MultiRegions::DirCartesianMap[1], Vx, Dummy1);
                elmt->PhysDeriv(MultiRegions::DirCartesianMap[1], Uy, Dummy2);
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Dummy1, 1, Dummy2, 1,
                                                      Q[0],   1);
                Vmath::Smul(m_HBCdata[j].m_ptsInElmt, m_negWavenumberSq[j],
                                                      Vel[0], 1, Dummy1, 1);
                elmt->PhysDeriv(MultiRegions::DirCartesianMap[0], Wz, Dummy2);
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Q[0],   1, Dummy1, 1,
                                                      Q[0],   1);
                Vmath::Vadd(m_HBCdata[j].m_ptsInElmt, Q[0],   1, Dummy2, 1,
                                                      Q[0],   1);

                elmt->PhysDeriv(MultiRegions::DirCartesianMap[1], Wz, Dummy1);
                Vmath::Smul(m_HBCdata[j].m_ptsInElmt, m_negWavenumberSq[j],
                                                      Vel[1], 1, Dummy2, 1);
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Dummy1, 1, Dummy2, 1,
                                                      Q[1],   1);
                elmt->PhysDeriv(MultiRegions::DirCartesianMap[0], Vx, Dummy1);
                elmt->PhysDeriv(MultiRegions::DirCartesianMap[0], Uy, Dummy2);
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Q[1],   1, Dummy1, 1,
                                                      Q[1],   1);
                Vmath::Vadd(m_HBCdata[j].m_ptsInElmt, Q[1],   1, Dummy2, 1,
                                                      Q[1],   1);
            }
            break;

            case MultiRegions::e3DH2D:
            {
                Array<OneD,NekDouble> Wx(m_pressureBCsElmtMaxPts);
                Array<OneD,NekDouble> Wz(m_pressureBCsElmtMaxPts);
                Array<OneD,NekDouble> Uz(m_pressureBCsElmtMaxPts);
                Array<OneD,NekDouble> qz(m_pressureBCsElmtMaxPts);
                Array<OneD,NekDouble> qy(m_pressureBCsElmtMaxPts);

                elmt->PhysDeriv(MultiRegions::DirCartesianMap[0], Vel[2], Wx);
                elmt->PhysDeriv(MultiRegions::DirCartesianMap[0], Vel[1], Vx);

                Vmath::Smul(m_HBCdata[j].m_ptsInElmt, m_negWavenumberSq[j],
                                                      Vel[0], 1, Uy, 1);

                Vmath::Smul(m_HBCdata[j].m_ptsInElmt, m_wavenumber[j],
                                                      Vel[0], 1, Uz, 1);

                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Wz,     1, Wx, 1,
                                                      qy,     1);
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Vx,     1, Uy, 1,
                                                      qz,     1);

                Vmath::Smul(m_HBCdata[j].m_ptsInElmt, m_negWavenumberSq[j],
                                                      qz,     1, Uy, 1);

                Vmath::Smul(m_HBCdata[j].m_ptsInElmt, m_wavenumber[j],
                                                      qy,     1, Uz, 1);

                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Uy,     1, Uz, 1,
                                                      Q[0],   1);
            }
            break;

            case MultiRegions::e3D:
            {
                Array<OneD,NekDouble> Dummy(m_pressureBCsElmtMaxPts);
                Array<OneD,NekDouble> Vz(m_pressureBCsElmtMaxPts);
                Array<OneD,NekDouble> Uz(m_pressureBCsElmtMaxPts);
                Array<OneD,NekDouble> Wx(m_pressureBCsElmtMaxPts);
                Array<OneD,NekDouble> Wy(m_pressureBCsElmtMaxPts);

                elmt->PhysDeriv(Vel[0], Dummy, Uy, Uz);
                elmt->PhysDeriv(Vel[1], Vx, Dummy, Vz);
                elmt->PhysDeriv(Vel[2], Wx, Wy, Dummy);

                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Wy, 1, Vz, 1, Q[0], 1);
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Uz, 1, Wx, 1, Q[1], 1);
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Vx, 1, Uy, 1, Q[2], 1);

                elmt->PhysDeriv(Q[0], Dummy, Wy, Vx);
                elmt->PhysDeriv(Q[1], Wx, Dummy, Uz);
                elmt->PhysDeriv(Q[2], Vz, Uy, Dummy);

                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Uy, 1, Uz, 1, Q[0], 1);
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Vx, 1, Vz, 1, Q[1], 1);
                Vmath::Vsub(m_HBCdata[j].m_ptsInElmt, Wx, 1, Wy, 1, Q[2], 1);
            }
            break;
            default:
                ASSERTL0(0,"Dimension not supported");
                break;
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
    
        int HBCnumber = 0;
        for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
        {
            // High order boundary condition;
            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
            {
                cnt += m_PBndExp[n]->GetNcoeffs();
                HBCnumber += m_PBndExp[n]->GetExpSize();
            }
        }
    
        int checkHBC = HBCnumber;
        m_comm->AllReduce(checkHBC,LibUtilities::ReduceSum);
        //ASSERTL0(checkHBC > 0 ,"At least one high-order pressure boundary "
        //                       "condition is required for scheme "
        //                       "consistency");

        m_acceleration[0] = Array<OneD, NekDouble>(cnt, 0.0);
        for(n = 0; n < m_intSteps; ++n)
        {
            m_pressureHBCs[n]   = Array<OneD, NekDouble>(cnt, 0.0);
            m_acceleration[n+1] = Array<OneD, NekDouble>(cnt, 0.0);
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
    
        
        m_HBCdata = Array<OneD, HBCInfo>(HBCnumber);
        StdRegions::StdExpansionSharedPtr elmt;
    
        switch(m_pressure->GetExpType())
        {
            case MultiRegions::e2D:
            case MultiRegions::e3D:
            {
                int coeff_count = 0;
                int exp_size;
                int j=0;
                int cnt = 0;
                for(int n = 0 ; n < m_PBndConds.num_elements(); ++n)
                {
                    exp_size = m_PBndExp[n]->GetExpSize();

                    if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
                    {
                        for(int i = 0; i < exp_size; ++i,cnt++)
                        {
                            m_HBCdata[j].m_globalElmtID =
                                                m_pressureBCtoElmtID[cnt];
                            elmt = m_pressure->GetExp(
                                                m_HBCdata[j].m_globalElmtID);
                            m_HBCdata[j].m_ptsInElmt    =
                                                elmt->GetTotPoints();
                            m_HBCdata[j].m_physOffset   =
                                                m_pressure->GetPhys_Offset(
                                                m_HBCdata[j].m_globalElmtID);
                            m_HBCdata[j].m_bndElmtOffset = i;
                            m_HBCdata[j].m_elmtTraceID  =
                                                m_pressureBCtoTraceID[cnt];
                            m_HBCdata[j].m_bndryElmtID  = n;
                            m_HBCdata[j].m_coeffOffset  = coeff_count;
                            coeff_count += elmt->GetEdgeNcoeffs(
                                                m_HBCdata[j].m_elmtTraceID);
                            j = j+1;
                        }
                    }
                    else // setting if just standard BC no High order
                    {
                        cnt += exp_size;
                    }
                }
            }
            break;
            
            case MultiRegions::e3DH1D:
            {
                Array<OneD, unsigned int> planes;
                planes = m_pressure->GetZIDs();
                int num_planes = planes.num_elements();            
                int num_elm_per_plane = (m_pressure->GetExpSize())/num_planes;
                
                m_wavenumber      = Array<OneD, NekDouble>(HBCnumber);
                m_negWavenumberSq = Array<OneD, NekDouble>(HBCnumber);
        
                int exp_size, exp_size_per_plane;
                int j=0;
                int K;
                NekDouble sign = -1.0;
                int cnt = 0;
                
                m_session->MatchSolverInfo("ModeType", "SingleMode", 
                                           m_SingleMode, false);
                m_session->MatchSolverInfo("ModeType", "HalfMode", 
                                           m_HalfMode, false);
                m_session->MatchSolverInfo("ModeType", "MultipleModes", 
                                           m_MultipleModes, false);
                m_session->LoadParameter("LZ", m_LhomZ);

                // Stability Analysis flags
                if(m_session->DefinesSolverInfo("ModeType"))
                {
                    if(m_SingleMode)
                    {
                        m_npointsZ = 2;
                    }
                    else if(m_HalfMode)
                    {
                        m_npointsZ = 1;
                    }
                    else if(m_MultipleModes)
                    {
                        m_npointsZ = m_session->GetParameter("HomModesZ");
                    }
                    else
                    {
                        ASSERTL0(false, "SolverInfo ModeType not valid");
                    }
                }
                else 
                {
                    m_npointsZ = m_session->GetParameter("HomModesZ");
                }

                   Array<OneD, int> coeff_count(m_PBndConds.num_elements(),0);
                   Array<OneD, int> coeffPlaneOffset(m_PBndConds.num_elements(),0);

                   cnt = 0;
                   for(int n = 0 ; n < m_PBndConds.num_elements(); ++n)
                   {
                       coeffPlaneOffset[n] = cnt;
                       if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
                       {
                           cnt += m_PBndExp[n]->GetNcoeffs();
                       }
                   }

                cnt = 0;
                for(int k = 0; k < num_planes; k++)
                {
                    K = planes[k]/2;
                    for(int n = 0 ; n < m_PBndConds.num_elements(); ++n)
                    {
                        exp_size = m_PBndExp[n]->GetExpSize();
                        exp_size_per_plane = exp_size/num_planes;
            
                        if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
                        {
                            for(int i = 0; i < exp_size_per_plane; ++i,cnt++)
                            {
                                m_HBCdata[j].m_globalElmtID = m_pressureBCtoElmtID[cnt];   
                                elmt        = m_pressure->GetExp(m_HBCdata[j].m_globalElmtID);
                                m_HBCdata[j].m_ptsInElmt = elmt->GetTotPoints();         
                                m_HBCdata[j].m_physOffset = m_pressure->GetPhys_Offset(m_HBCdata[j].m_globalElmtID);
                                m_HBCdata[j].m_bndElmtOffset = i+k*exp_size_per_plane;       
                                m_HBCdata[j].m_elmtTraceID = m_pressureBCtoTraceID[cnt];      
                                m_HBCdata[j].m_bndryElmtID = n;
                                m_HBCdata[j].m_coeffOffset = coeffPlaneOffset[n] + coeff_count[n];
                                   coeff_count[n] += elmt->GetEdgeNcoeffs(m_HBCdata[j].m_elmtTraceID);
    
                                if(m_SingleMode)
                                {
                                    m_wavenumber[j]      = -2*M_PI/m_LhomZ;       
                                    m_negWavenumberSq[j] = -1.0*m_wavenumber[j]*m_wavenumber[j];
                                }
                                else if(m_HalfMode || m_MultipleModes)
                                {
                                    m_wavenumber[j]      = 2*M_PI/m_LhomZ;       
                                    m_negWavenumberSq[j] = -1.0*m_wavenumber[j]*m_wavenumber[j];
                                }
                                else
                                {
                                    m_wavenumber[j]     = 2*M_PI*sign*(NekDouble(K))/m_LhomZ; 
                                    m_negWavenumberSq[j] = -1.0*m_wavenumber[j]*m_wavenumber[j];
                                }
                                
                                int assElmtID;
                
                                if(k%2==0)
                                {
                                    if(m_HalfMode)
                                    {
                                        assElmtID = m_HBCdata[j].m_globalElmtID;
                                        
                                    }
                                    else
                                    {
                                        assElmtID = m_HBCdata[j].m_globalElmtID + num_elm_per_plane;
                                    }
                                }
                                else 
                                {
                                    assElmtID = m_HBCdata[j].m_globalElmtID - num_elm_per_plane;
                                }
                
                                m_HBCdata[j].m_assPhysOffset = m_pressure->GetPhys_Offset(assElmtID);
                
                                j = j+1;
                            }
                        }
                        else // setting if just standard BC no High order
                        {
                            cnt += exp_size_per_plane;
                        }
                    }
                    sign = -1.0*sign;
                }
            }
            break;
            
            case MultiRegions::e3DH2D:
            {
                int cnt = 0;
                int exp_size, exp_size_per_line;
                int j=0;
        
                for(int k1 = 0; k1 < m_npointsZ; k1++)
                {
                    for(int k2 = 0; k2 < m_npointsY; k2++)
                    {
                        for(int n = 0 ; n < m_PBndConds.num_elements(); ++n)
                        {
                            exp_size = m_PBndExp[n]->GetExpSize();
                            
                            exp_size_per_line = exp_size/(m_npointsZ*m_npointsY);
                            
                            if(boost::iequals(m_PBndConds[n]->GetUserDefined(),"H"))
                            {
                                for(int i = 0; i < exp_size_per_line; ++i,cnt++)
                                {
                                    // find element and edge of this expansion. 
                                    m_HBCdata[j].m_globalElmtID = m_pressureBCtoElmtID[cnt];
                                    elmt      = m_pressure->GetExp(m_HBCdata[j].m_globalElmtID);
                                    m_HBCdata[j].m_ptsInElmt = elmt->GetTotPoints();
                                    m_HBCdata[j].m_physOffset = m_pressure->GetPhys_Offset(m_HBCdata[j].m_globalElmtID);
                                    m_HBCdata[j].m_bndElmtOffset = i+(k1*m_npointsY+k2)*exp_size_per_line;
                                    m_HBCdata[j].m_elmtTraceID = m_pressureBCtoTraceID[cnt];                
                                    m_HBCdata[j].m_bndryElmtID = n;
                                }
                            }
                            else
                            {
                                cnt += exp_size_per_line;
                            }
                        }
                    }
                }
            }
            break;
            default:
                ASSERTL0(0,"Dimension not supported");
                break;
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

}

