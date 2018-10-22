///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionIP.cpp
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
// Description: IP diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/DiffusionIP.h>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string DiffusionIP::type = GetDiffusionFactory().
            RegisterCreatorFunction("IP", DiffusionIP::create);

        DiffusionIP::DiffusionIP()
        {
        }

        void DiffusionIP::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            // m_session = pSession;

            // m_session->LoadSolverInfo("ShockCaptureType",
            //                       m_shockCaptureType,    "Off");		

            // Setting up the normals
            int i;
            int nDim = pFields[0]->GetCoordim(0);
            int nTracePts = pFields[0]->GetTrace()->GetTotPoints();
            
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            pFields[0]->GetTrace()->GetNormals(m_traceNormals);
        }
        
        void DiffusionIP::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            
            int nCoeffs   = fields[0]->GetNcoeffs();
            Array<OneD, Array<OneD, NekDouble> > tmp(nConvectiveFields);
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);
            }
            DiffusionIP::v_Diffuse_coeff(nConvectiveFields,fields,inarray,tmp,pFwd,pBwd);
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                fields[i]->BwdTrans             (tmp[i], outarray[i]);
            }
        }

        void DiffusionIP::v_Diffuse_coeff(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            int nBndEdgePts, i, j, k, e;
            int nDim      = fields[0]->GetCoordim(0);
            int nPts      = fields[0]->GetTotPoints();
            int nCoeffs   = fields[0]->GetNcoeffs();
            int nTracePts = fields[0]->GetTrace()->GetTotPoints();
            Array<OneD, NekDouble>  tmp(nCoeffs);

            // NOTE: the 1st&2nd indexes are exchanged compared with DiffusionLDG;
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > flux  (nDim);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > qfield(nDim);
            for (j = 0; j < nDim; ++j)
            {
                qfield[j] = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                flux[j]   = Array<OneD, Array<OneD, NekDouble> >(nConvectiveFields);
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    qfield[j][i] = Array<OneD, NekDouble>(nPts, 0.0);
                    flux[j][i]   = Array<OneD, NekDouble>(nTracePts, 0.0);
                }
            }

            // Compute q_{\eta} and q_{\xi}
            // Obtain numerical fluxes
            Array< Array<OneD, NekDouble> > qtmp(3);
            for(int nd=0; nd<3, nd++)
            {
                qtmp[nd]    =   NullNekDouble1DArray;
            }
            for(int i = 0; i < nConvectiveFields; ++i)
            {
                for(int nd=0; nd<nDim, nd++)
                {
                    qtmp[nd]    =   qfield[nd][i];
                }
                m_fields[i]->PhysDeriv(inarray[i], qtmp[0], qtmp[1], qtmp[2]);
            }

            Array<OneD, NekDouble> muvar        =   NullNekDouble1DArray;
            Array<OneD, NekDouble> FwdMuVar     =   NullNekDouble1DArray;
            Array<OneD, NekDouble> BwdMuVar     =   NullNekDouble1DArray;
            Array<OneD, NekDouble> MuVarTrace   =   NullNekDouble1DArray;
            if (m_ArtificialDiffusionVector)
            {
                muvar   =   Array<OneD, NekDouble>(nPts, 0.0);
                m_ArtificialDiffusionVector(inarray, muvar);

                FwdMuVar =  Array<OneD, NekDouble>(nTracePts, 0.0);
                BwdMuVar =  Array<OneD, NekDouble>(nTracePts, 0.0);

                // TODO: CHECK !!
                fields[0]->GetFwdBwdTracePhysInterior(muvar,FwdMuVar,BwdMuVar);
                fields[0]->FillBwdWITHBoundZero(FwdMuVar,BwdMuVar);

                for(k = 0; k < nTracePts; ++k)
                {
                    FwdMuVar[k] = 0.5 * (FwdMuVar[k] + BwdMuVar[k]) ;
                }
                MuVarTrace  =   FwdMuVar;
                FwdMuVar    =   NullNekDouble1DArray;
                BwdMuVar    =   NullNekDouble1DArray;
            }
            
            Array<OneD, NekDouble>    Fwd   =   NullNekDouble1DArray;
            Array<OneD, NekDouble>    Bwd   =   NullNekDouble1DArray;
            Bwd = Array<OneD, NekDouble>(nTracePts,0.0);
            
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > >    numDeriv(nDim);
            for (int nd = 0; nd < nDim; ++nd)
            {
                numDeriv[nd]     =   Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    numDeriv[nd][i]    = Array<OneD, NekDouble>(nTracePts,0.0);
                    Fwd =  numDeriv[nd][i]; 
                    Vmath::Fill(nTracePts, 1.0, Bwd,1);
                    fields[i]->GetFwdBwdTracePhysDeriv(qfield[nd][i], Fwd, Bwd);
                    for (int nt = 0; nt < nTracePts; ++nt)
                    {
                        numDeriv[nd][i][nt] =   0.5*( Fwd[nt] + Bwd[nt] );  
                    }
                }
            }
            // release storage
            Bwd   =   NullNekDouble1DArray;

            // volume intergration 
            // release qfield and muvar;
            
            Array<OneD, Array<OneD, NekDouble> >    solution_jump(nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> >    solution_Aver(nConvectiveFields);
            for (i = 0; i < nConvectiveFields; ++i)
            {
                solution_jump[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
                solution_Aver[i]    =   Array<OneD, NekDouble>(nTracePts,0.0);
            }

            // Careful Fwd and solution_jump, solution_Aver and Bwd share same storage
            if (pFwd == NullNekDoubleArrayofArray ||
                pBwd == NullNekDoubleArrayofArray)
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd       =    solution_jump[i];  
                    Bwd       =    solution_Aver[i]
                    fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd, Bwd);
                    for (int nt = 0; nt < nTracePts; ++nt)
                    {
                        solution_jump[i][nt]   =   Fwd[nt] - Bwd[nt];  
                        solution_Aver[i][nt]   =   0.5*solution_jump[nt] + Bwd[nt];  
                    }
                }
            }
            else
            {
                Fwd   =   NullNekDouble1DArray;
                Bwd   =   NullNekDouble1DArray;
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd  = pFwd[i];
                    Bwd  = pBwd[i];
                    for (int nt = 0; nt < nTracePts; ++nt)
                    {
                        solution_jump[i][nt]   =   Fwd[nt] - Bwd[nt];  
                        solution_Aver[i][nt]   =   0.5*solution_jump[nt] + Bwd[nt];  
                    }
                }
            }

            Array<OneD, NekDouble>  PenaltyFactor(nTracePts,0.0);
            getPenaltyFactor(PenaltyFactor);
            for (i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Vmul(nTracePts,solution_jump[i],1, PenaltyFactor,1,solution_jump[i],1);
            }

            Array<OneD, NekDouble> ElmtLength;
            ElmtLength  =   PenaltyFactor; 
            PenaltyFactor = NullNekDouble1DArray;
            Vmath::Fill(nTracePts, 1.0, ElmtLength,1);
            // getElmtLength(ElmtLength);
            for (i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Vdiv(nTracePts,solution_jump[i],1, ElmtLength,1,solution_jump[i],1);
            }

            for (int nd = 0; nd < nDim; ++nd)
            {
                for (i = 0; i < nConvectiveFields; ++i)
                {
                    Vmath::Vvtvp(nTracePts, m_traceNormals[nd],1,solution_jump[i],1, numDeriv[nd][i],1, numDeriv[nd][i],1);
                }
            }

            // Calculate normal viscous flux
            // viscousFlux(fields,solution_Aver,numDeriv,m_traceNormals,traceflux);

            for(i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Neg                      (nCoeffs, outarray[i], 1);
                fields[i]->AddTraceIntegral     (traceflux[i], outarray[i]);
                fields[i]->MultiplyByElmtInvMass(outarray[i], outarray[i]);
            }
        }
        
        
        void DiffusionIP::v_WeakPenaltyforScalar(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const Array<OneD, const NekDouble>                &ufield,
            const Array<OneD, const NekDouble>                &uplus,
                  Array<OneD,       NekDouble>                &penaltyflux)
        {
            int i, e, id1, id2;
            
            // Number of boundary regions
            int nBndEdgePts, nBndEdges;
            int cnt         = 0;
            int nBndRegions = fields[var]->GetBndCondExpansions().num_elements();

            for (i = 0; i < nBndRegions; ++i)
            {
                // Number of boundary expansion related to that region
                nBndEdges = fields[var]->
                GetBndCondExpansions()[i]->GetExpSize();

                // Weakly impose boundary conditions by modifying flux values
                for (e = 0; e < nBndEdges ; ++e)
                {
                    nBndEdgePts = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetTotPoints();
                    
                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));

                    // For Dirichlet boundary condition: uflux = g_D
                    if (fields[var]->GetBndConditions()[i]->
                        GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vcopy(nBndEdgePts, 
                                     &(fields[var]->
                                       GetBndCondExpansions()[i]->
                                       GetPhys())[id1], 1, 
                                     &penaltyflux[id2], 1);
                    }
                    // For Neumann boundary condition: uflux = u+
                    else if ((fields[var]->GetBndConditions()[i])->
                        GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vcopy(nBndEdgePts, 
                                     &uplus[id2], 1, 
                                     &penaltyflux[id2], 1);
                    }
                }
            }
        }
        
        
        
        void DiffusionIP::v_NumFluxforVector(
            const Array<OneD, MultiRegions::ExpListSharedPtr>        &fields,
            const Array<OneD, Array<OneD, NekDouble> >               &ufield,
                  Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &qfield,
                  Array<OneD, Array<OneD, NekDouble> >               &qflux)
        {
            int i, j;
            int nTracePts  = fields[0]->GetTrace()->GetTotPoints();
            int nvariables = fields.num_elements();
            int nDim       = qfield.num_elements();
            
            NekDouble C11 = 0.0;
            Array<OneD, NekDouble > Fwd(nTracePts);
            Array<OneD, NekDouble > Bwd(nTracePts);
            Array<OneD, NekDouble > Vn (nTracePts, 0.0);
            Array<OneD, NekDouble > qFwd     (nTracePts);
            Array<OneD, NekDouble > qBwd     (nTracePts);
            Array<OneD, NekDouble > qfluxtemp(nTracePts, 0.0);
            Array<OneD, NekDouble > uterm(nTracePts);
            /*
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            */
            
            // Get the normal velocity Vn
            for(i = 0; i < nDim; ++i)
            {
                Vmath::Svtvp(nTracePts, 1.0, m_traceNormals[i], 1, 
                             Vn, 1, Vn, 1);
            }
            // Evaulate upwind flux:
            // qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)
            for (i = 0; i < nvariables; ++i)
            {
                // Generate Stability term = - C11 ( u- - u+ )
                fields[i]->GetFwdBwdTracePhys(ufield[i], Fwd, Bwd);

                Vmath::Vsub(nTracePts,
                            Fwd, 1, Bwd, 1,
                            uterm, 1);

                Vmath::Smul(nTracePts,
                            -1.0 * C11, uterm, 1,
                            uterm, 1);

                qflux[i] = Array<OneD, NekDouble> (nTracePts, 0.0);
                for (j = 0; j < nDim; ++j)
                {
                    //  Compute Fwd and Bwd value of ufield of jth direction
                    fields[i]->GetFwdBwdTracePhys(qfield[j][i], qFwd, qBwd);
                    
                    // if Vn >= 0, flux = uFwd, i.e.,
                    // edge::eForward, if V*n>=0 <=> V*n_F>=0, pick 
                    // qflux = qBwd = q+
                    // edge::eBackward, if V*n>=0 <=> V*n_B<0, pick 
                    // qflux = qBwd = q-
                    
                    // else if Vn < 0, flux = uBwd, i.e.,
                    // edge::eForward, if V*n<0 <=> V*n_F<0, pick 
                    // qflux = qFwd = q-
                    // edge::eBackward, if V*n<0 <=> V*n_B>=0, pick 
                    // qflux = qFwd = q+
                    
                    fields[i]->GetTrace()->Upwind(/*m_traceNormals[j]*/Vn, 
                                                  qBwd, qFwd, 
                                                  qfluxtemp);
                    
                    Vmath::Vmul(nTracePts, 
                                m_traceNormals[j], 1, 
                                qfluxtemp, 1, 
                                qfluxtemp, 1);

                    // Flux = {Fwd, Bwd} * (nx, ny, nz) + uterm * (nx, ny)
                    Vmath::Vadd(nTracePts, 
                                uterm, 1, 
                                qfluxtemp, 1, 
                                qfluxtemp, 1);
                    
                    Array<OneD, NekDouble > qtemp(nTracePts);
                    fields[i]->ExtractTracePhys(qfield[j][i], qtemp);

                    // Imposing weak boundary condition with flux
                    if (fields[0]->GetBndCondExpansions().num_elements())
                    {
                        v_WeakPenaltyforVector(fields, i, j, 
                                               qfield[j][i], qtemp,
                                               qfluxtemp, C11);
                    }
                    
                    // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
                    // n_xi = n_x * tan_xi_x + n_y * tan_xi_y + n_z * tan_xi_z
                    // n_xi = n_x * tan_eta_x + n_y * tan_eta_y + n_z*tan_eta_z
                    Vmath::Vadd(nTracePts, 
                                qfluxtemp, 1, 
                                qflux[i], 1, 
                                qflux[i], 1);
                }
            }
        }
         
        
        
        /**
         * Diffusion: Imposing weak boundary condition for q with flux
         *  uflux = g_D  on Dirichlet boundary condition
         *  uflux = u_Fwd  on Neumann boundary condition
         */
        void DiffusionIP::v_WeakPenaltyforVector(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const int                                          var,
            const int                                          dir,
            const Array<OneD, const NekDouble>                &qfield,
            const Array<OneD, const NekDouble>                &qtemp,
                  Array<OneD,       NekDouble>                &penaltyflux,
            NekDouble                                          C11)
        {
            int i, e, id1, id2;
            int nBndEdges, nBndEdgePts;
            int nBndRegions = fields[var]->GetBndCondExpansions().num_elements();
            int cnt = 0;

            /*
            // Setting up the normals
            m_traceNormals = Array<OneD, Array<OneD, NekDouble> >(nDim);
            for(i = 0; i < nDim; ++i)
            {
                m_traceNormals[i] = Array<OneD, NekDouble> (nTracePts);
            }
            fields[0]->GetTrace()->GetNormals(m_traceNormals);
            */
            
            for (i = 0; i < nBndRegions; ++i)
            {
                nBndEdges = fields[var]->
                    GetBndCondExpansions()[i]->GetExpSize();
                                
                // Weakly impose boundary conditions by modifying flux values
                for (e = 0; e < nBndEdges ; ++e)
                {
                    nBndEdgePts = fields[var]->
                    GetBndCondExpansions()[i]->GetExp(e)->GetTotPoints();
                    
                    id1 = fields[var]->
                    GetBndCondExpansions()[i]->GetPhys_Offset(e);
                    
                    id2 = fields[0]->GetTrace()->
                    GetPhys_Offset(fields[0]->GetTraceMap()->
                                   GetBndCondTraceToGlobalTraceMap(cnt++));
                    
                    // For Dirichlet boundary condition: 
                    //qflux = q+ - C_11 (u+ -    g_D) (nx, ny)
                    if(fields[var]->GetBndConditions()[i]->
                    GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        Vmath::Vmul(nBndEdgePts, 
                                    &m_traceNormals[dir][id2], 1, 
                                    &qtemp[id2], 1, 
                                    &penaltyflux[id2], 1);
                    }
                    // For Neumann boundary condition: qflux = g_N
                    else if((fields[var]->GetBndConditions()[i])->
                    GetBoundaryConditionType() == SpatialDomains::eNeumann)
                    {
                        Vmath::Vmul(nBndEdgePts,
                                    &m_traceNormals[dir][id2], 1, 
                                    &(fields[var]->
                                      GetBndCondExpansions()[i]->
                                      GetPhys())[id1], 1, 
                                    &penaltyflux[id2], 1);
                    }
                }
            }
        }

    }
}
