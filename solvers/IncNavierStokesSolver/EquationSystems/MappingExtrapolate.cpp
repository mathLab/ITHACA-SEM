///////////////////////////////////////////////////////////////////////////////
//
// File: StandardExtrapolate.cpp
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
// Description: Abstract base class for StandardExtrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/MappingExtrapolate.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
    /**
     * Registers the class with the Factory.
     */
    std::string MappingExtrapolate::className = GetExtrapolateFactory().RegisterCreatorFunction(
        "Mapping",
        MappingExtrapolate::create,
        "Mapping");

    MappingExtrapolate::MappingExtrapolate(
        const LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
        MultiRegions::ExpListSharedPtr pPressure,
        const Array<OneD, int> pVel,
        const SolverUtils::AdvectionSharedPtr advObject)
        : StandardExtrapolate(pSession,pFields,pPressure,pVel,advObject)
    {
        m_mapping = GlobalMapping::Mapping::Load(m_session, m_fields);
        
        // Load solve parameters related to the mapping
        // Flags determining if pressure/viscous terms should be treated implicitly
        m_session->MatchSolverInfo("MappingImplicitPressure","True",m_implicitPressure,false);
        m_session->MatchSolverInfo("MappingImplicitViscous","True",m_implicitViscous,false);
        
        // Relaxation parameter for pressure system
        m_session->LoadParameter("MappingPressureRelaxation",m_pressureRelaxation,1.0);
    }

    MappingExtrapolate::~MappingExtrapolate()
    {
    }
    
    /** 
     * 
     */
    void MappingExtrapolate::v_CorrectPressureBCs(
        const Array<OneD, NekDouble>  &pressure)
    {
        if(m_HBCdata.num_elements()>0)
        {
            int cnt, n;
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_fields.num_elements()-1;
            
            Array<OneD, NekDouble> Vals;
            // Remove previous correction
            for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
            {
                if(m_PBndConds[n]->GetUserDefined() == "H")
                {
                    int nq = m_PBndExp[n]->GetNcoeffs();
                    Vmath::Vsub(nq, &(m_PBndExp[n]->GetCoeffs()[0]),  1,
                                    &(m_bcCorrection[cnt]), 1, 
                                    &(m_PBndExp[n]->UpdateCoeffs()[0]), 1);
                    cnt += nq;
                }
            }
            
            // Calculate new correction
            Array<OneD, NekDouble> Jac(physTot, 0.0);
            m_mapping->GetJacobian(Jac);
            
            Array<OneD, Array<OneD, NekDouble> > correction(nvel);
            Array<OneD, Array<OneD, NekDouble> > gradP(nvel);
            Array<OneD, Array<OneD, NekDouble> > wk(nvel);
            Array<OneD, Array<OneD, NekDouble> > wk2(nvel);
            for (int i=0; i<nvel; i++)
            {
                wk[i] = Array<OneD, NekDouble> (physTot, 0.0);
                gradP[i] = Array<OneD, NekDouble> (physTot, 0.0);
                correction[i] = Array<OneD, NekDouble> (physTot, 0.0);
            }

            // Calculate G(p)
            for(int i = 0; i < nvel; ++i)
            {
                m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i], pressure, gradP[i]);
                if(m_fields[0]->GetWaveSpace())
                {
                    m_fields[0]->HomogeneousBwdTrans(gradP[i], wk[i]);
                }
                else
                {
                    Vmath::Vcopy(physTot, gradP[i], 1, wk[i], 1);
                }
            }
            m_mapping->RaiseIndex(wk, correction);   // G(p)

            // alpha*J*(G(p))
            if (!m_mapping->HasConstantJacobian())
            {
                for(int i = 0; i < nvel; ++i)
                {
                    Vmath::Vmul(physTot, correction[i], 1, Jac, 1, correction[i], 1);
                }                
            }   
            for(int i = 0; i < nvel; ++i)
            {
                Vmath::Smul(physTot, m_pressureRelaxation, correction[i], 1, correction[i], 1); 
            }
            
            if(m_pressure->GetWaveSpace())
            {
                for(int i = 0; i < nvel; ++i)
                {
                    m_pressure->HomogeneousFwdTrans(correction[i], correction[i]);
                }
            }            
            // p_i - alpha*J*div(G(p))    
            for (int i = 0; i < nvel; ++i)
            {
                Vmath::Vsub(physTot, gradP[i], 1, correction[i], 1, correction[i], 1);
            }
            
            // Get value at boundary and calculate Inner product
            StdRegions::StdExpansionSharedPtr Pbc;
            StdRegions::StdExpansionSharedPtr elmt;
            Array<OneD, Array<OneD, const NekDouble> > correctionElmt(m_bnd_dim);
            Array<OneD, Array<OneD, NekDouble> > BndValues(m_bnd_dim);
            for(int i = 0; i < m_bnd_dim; i++)
            {
                BndValues[i] = Array<OneD, NekDouble> (m_pressureBCsMaxPts,0.0);
            }
            for(int j = 0 ; j < m_HBCdata.num_elements() ; j++)
            {
                /// Casting the boundary expansion to the specific case
                Pbc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion> 
                            (m_PBndExp[m_HBCdata[j].m_bndryID]
                                ->GetExp(m_HBCdata[j].m_bndElmtID));

                /// Picking up the element where the HOPBc is located
                elmt = m_pressure->GetExp(m_HBCdata[j].m_globalElmtID);

                /// Assigning
                for(int i = 0; i < m_bnd_dim; i++)
                {
                    correctionElmt[i]  = correction[i] + m_HBCdata[j].m_physOffset;
                }
                Vals = m_bcCorrection + m_HBCdata[j].m_coeffOffset;
                // Getting values on the edge and filling the correction
                switch(m_pressure->GetExpType())
                {
                    case MultiRegions::e2D:
                    case MultiRegions::e3DH1D:
                    {                                                         
                        elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                              correctionElmt[0], BndValues[0]);                    
                        elmt->GetEdgePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                              correctionElmt[1], BndValues[1]);

                        // InnerProduct 
                        Pbc->NormVectorIProductWRTBase(BndValues[0], BndValues[1],
                                                       Vals);
                    }
                    break;

                    case MultiRegions::e3D:
                    {
                        elmt->GetFacePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc, 
                                              correctionElmt[0], BndValues[0]);
                        elmt->GetFacePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                              correctionElmt[1], BndValues[1]);
                        elmt->GetFacePhysVals(m_HBCdata[j].m_elmtTraceID, Pbc,
                                              correctionElmt[2], BndValues[2]);
                        Pbc->NormVectorIProductWRTBase(BndValues[0], BndValues[1],
                                              BndValues[2], Vals);
                    }
                    break;
                default:
                    ASSERTL0(0,"Dimension not supported");
                    break;
                }
            }      
            
            // Apply new correction
            for(cnt = n = 0; n < m_PBndConds.num_elements(); ++n)
            {
                if(m_PBndConds[n]->GetUserDefined() == "H")
                {
                    int nq = m_PBndExp[n]->GetNcoeffs();
                    Vmath::Vadd(nq, &(m_PBndExp[n]->GetCoeffs()[0]),  1,
                                    &(m_bcCorrection[cnt]), 1, 
                                    &(m_PBndExp[n]->UpdateCoeffs()[0]), 1);
                    cnt += nq;
                }
            }                
        }        
    }

    void MappingExtrapolate::v_CalcNeumannPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        const Array<OneD, const Array<OneD, NekDouble> >  &N,
        NekDouble kinvis)
    {
        if (m_mapping->HasConstantJacobian() && !m_implicitViscous)
        {
            Extrapolate::v_CalcNeumannPressureBCs( fields, N, kinvis);
        }
        else
        {
            int physTot = m_fields[0]->GetTotPoints();
            int nvel = m_fields.num_elements()-1;
            
            Array<OneD, NekDouble> Pvals;
            StdRegions::StdExpansionSharedPtr Pbc;
            StdRegions::StdExpansionSharedPtr elmt;

            Array<OneD, Array<OneD, const NekDouble> > Velocity(m_bnd_dim);
            Array<OneD, Array<OneD, const NekDouble> > Advection(m_bnd_dim);
            // Get transformation Jacobian
            Array<OneD, NekDouble> Jac(physTot,0.0);
            m_mapping->GetJacobian(Jac);
            // Declare variables
            Array<OneD, Array<OneD, NekDouble> > BndValues(m_bnd_dim);
            Array<OneD, Array<OneD, NekDouble> > Q(m_bnd_dim);
            Array<OneD, Array<OneD, NekDouble> > Q_field(nvel);
            Array<OneD, Array<OneD, NekDouble> > fields_new(nvel);
            Array<OneD, Array<OneD, NekDouble> > N_new(m_bnd_dim);
            // Temporary variables
            Array<OneD, NekDouble> tmp(physTot,0.0);
            Array<OneD, NekDouble> tmp2(physTot,0.0);
            for(int i = 0; i < m_bnd_dim; i++)
            {
                BndValues[i] = Array<OneD, NekDouble> (m_pressureBCsMaxPts,0.0);
                Q[i]         = Array<OneD, NekDouble> (m_pressureBCsElmtMaxPts,0.0);
                N_new[i]   = Array<OneD, NekDouble> (physTot,0.0);
            }
            for(int i = 0; i < nvel; i++)
            {
                Q_field[i]   = Array<OneD, NekDouble> (physTot,0.0);
                fields_new[i]   = Array<OneD, NekDouble> (physTot,0.0);
            }
            
            // Multiply convective terms by Jacobian
            for(int i = 0; i < m_bnd_dim; i++)
            {
                if (m_fields[0]->GetWaveSpace())
                {
                    m_fields[0]->HomogeneousBwdTrans(N[i],N_new[i]);
                }
                else
                {
                    Vmath::Vcopy(physTot, N[i], 1, N_new[i], 1);
                }
                Vmath::Vmul(physTot, Jac, 1, N_new[i], 1, N_new[i], 1);
                if (m_fields[0]->GetWaveSpace())
                {
                    m_fields[0]->HomogeneousFwdTrans(N_new[i],N_new[i]);
                }                          
            }
            
            // Get velocity in physical space
            for(int i = 0; i < nvel; i++)
            {
                if (m_fields[0]->GetWaveSpace())
                {
                    m_fields[0]->HomogeneousBwdTrans(fields[i],fields_new[i]);
                }        
                else
                {
                    Vmath::Vcopy(physTot, fields[i], 1, fields_new[i], 1);
                }
            }
            
            // Calculate appropriate form of the CurlCurl operator
            m_mapping->CurlCurlField(fields_new, Q_field, m_implicitViscous);
            
            // If viscous terms are treated explicitly,
            //     add grad(U/J \dot grad J) to CurlCurl
            if ( !m_implicitViscous)
            {
                m_mapping->DotGradJacobian(fields_new, tmp);
                Vmath::Vdiv(physTot, tmp, 1, Jac, 1, tmp, 1);
                
                bool wavespace = m_fields[0]->GetWaveSpace();
                m_fields[0]->SetWaveSpace(false);
                for(int i = 0; i < m_bnd_dim; i++)
                {
                    m_fields[0]->PhysDeriv(MultiRegions::DirCartesianMap[i], 
                                            tmp, tmp2);
                    Vmath::Vadd(physTot, Q_field[i], 1, tmp2, 1, Q_field[i], 1);
                }     
                m_fields[0]->SetWaveSpace(wavespace);
            }        
            
            // Multiply by Jacobian and convert to wavespace (if necessary)
            for(int i = 0; i < m_bnd_dim; i++)
            {
                Vmath::Vmul(physTot, Jac, 1, fields_new[i], 1, fields_new[i], 1);
                Vmath::Vmul(physTot, Jac, 1, Q_field[i]   , 1, Q_field[i]   , 1);
                if (m_fields[0]->GetWaveSpace())
                {
                    m_fields[0]->HomogeneousFwdTrans(fields_new[i],fields_new[i]);
                    m_fields[0]->HomogeneousFwdTrans(Q_field[i],Q_field[i]);
                }                          
            }            

            for(int j = 0 ; j < m_HBCdata.num_elements() ; j++)
            {
                /// Casting the boundary expansion to the specific case
                Pbc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion> 
                            (m_PBndExp[m_HBCdata[j].m_bndryID]
                                ->GetExp(m_HBCdata[j].m_bndElmtID));

                /// Picking up the element where the HOPBc is located
                elmt = m_pressure->GetExp(m_HBCdata[j].m_globalElmtID);

                /// Assigning
                for(int i = 0; i < m_bnd_dim; i++)
                {
                    Velocity[i]  = fields_new[i] + m_HBCdata[j].m_physOffset;
                    Advection[i] = N_new[i]      + m_HBCdata[j].m_physOffset;
                    Q[i]         = Q_field[i]   + m_HBCdata[j].m_physOffset;
                }

                // Mounting advection component into the high-order condition
                for(int i = 0; i < m_bnd_dim; i++)
                {
                    MountHOPBCs(m_HBCdata[j].m_ptsInElmt,kinvis,Q[i],Advection[i]);
                }

                Pvals = m_pressureHBCs[0] + m_HBCdata[j].m_coeffOffset;

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

                        // InnerProduct 
                        Pbc->NormVectorIProductWRTBase(BndValues[0], BndValues[1],
                                                       Pvals);
                    }
                    break;
                    case MultiRegions::e3DH2D:
                    {
                        if(m_HBCdata[j].m_elmtTraceID == 0)
                        {
                            (m_PBndExp[m_HBCdata[j].m_bndryID]->UpdateCoeffs()
                                + m_PBndExp[m_HBCdata[j].m_bndryID]
                                    ->GetCoeff_Offset(
                                        m_HBCdata[j].m_bndElmtID))[0]
                                                                    = -1.0*Q[0][0];
                        }
                        else if (m_HBCdata[j].m_elmtTraceID == 1)
                        {
                            (m_PBndExp[m_HBCdata[j].m_bndryID]->UpdateCoeffs()
                                + m_PBndExp[m_HBCdata[j].m_bndryID]
                                    ->GetCoeff_Offset(
                                        m_HBCdata[j].m_bndElmtID))[0] 
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
                    }
                    break;
                default:
                    ASSERTL0(0,"Dimension not supported");
                    break;
                }
            }            
        }
        // If pressure terms are treated implicitly, we need to multiply
        //     by the relaxation parameter, and zero the correction term
        if (m_implicitPressure)
        {
            Vmath::Smul(m_pressureHBCs[0].num_elements(), m_pressureRelaxation,
                            m_pressureHBCs[0],  1,
                            m_pressureHBCs[0], 1);
        } 
        m_bcCorrection  = Array<OneD, NekDouble> (m_pressureHBCs[0].num_elements(), 0.0);
    }

}

