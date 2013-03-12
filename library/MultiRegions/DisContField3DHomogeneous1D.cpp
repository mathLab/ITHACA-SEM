//////////////////////////////////////////////////////////////////////////////
//
// File DisContField3DHomogeneous1D.cpp
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
// Description: Field definition for 3D domain with boundary
// conditions using LDG flux and a 1D homogeneous direction
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <MultiRegions/DisContField2D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        DisContField3DHomogeneous1D::DisContField3DHomogeneous1D(void):
            ExpList3DHomogeneous1D(),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }

        DisContField3DHomogeneous1D::DisContField3DHomogeneous1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                                                 const LibUtilities::BasisKey &HomoBasis,
                                                                 const NekDouble lhom,
																 const bool useFFT,
																 const bool dealiasing):
            ExpList3DHomogeneous1D(pSession,HomoBasis,lhom,useFFT,dealiasing),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }

        DisContField3DHomogeneous1D::DisContField3DHomogeneous1D(const DisContField3DHomogeneous1D &In, const bool DeclarePlanesSetCoeffPhys):
            ExpList3DHomogeneous1D (In,false),
            m_bndCondExpansions    (In.m_bndCondExpansions),
            m_bndConditions        (In.m_bndConditions)
        {
            if(DeclarePlanesSetCoeffPhys)
            {
                bool False = false;
                DisContField2DSharedPtr zero_plane = boost::dynamic_pointer_cast<DisContField2D> (In.m_planes[0]);
                
                for(int n = 0; n < m_planes.num_elements(); ++n)
                {
                    m_planes[n] = MemoryManager<DisContField2D>::AllocateSharedPtr(*zero_plane,False);
                }
                
                SetCoeffPhys();
            }
        }

        DisContField3DHomogeneous1D::DisContField3DHomogeneous1D(
                                       const LibUtilities::SessionReaderSharedPtr &pSession,
                                       const LibUtilities::BasisKey &HomoBasis,
                                       const NekDouble lhom,
									   const bool useFFT,
									   const bool dealiasing,
                                       const SpatialDomains::MeshGraphSharedPtr &graph2D,
                                       const std::string &variable):
            ExpList3DHomogeneous1D(pSession,HomoBasis,lhom,useFFT,dealiasing),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            int i,n,nel;
            bool True  = true; 
            bool False = false; 
            DisContField2DSharedPtr plane_zero;
            SpatialDomains::BoundaryConditions bcs(m_session, graph2D);

            // note that nzplanes can be larger than nzmodes 
            m_planes[0] = plane_zero = MemoryManager<DisContField2D>::AllocateSharedPtr(pSession,graph2D,variable,True,False);

            m_exp = MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr();
            nel = m_planes[0]->GetExpSize();

            for(i = 0; i < nel; ++i)
            {
                (*m_exp).push_back(m_planes[0]->GetExp(i));
            }

            for(n = 1; n < m_planes.num_elements(); ++n)
            {
                m_planes[n] = MemoryManager<DisContField2D>::AllocateSharedPtr(*plane_zero,graph2D,variable,True,False);
                for(i = 0; i < nel; ++i)
                {
                    (*m_exp).push_back((*m_exp)[i]);
                }
            }            

            // Setup Default optimisation information. 
            nel = GetExpSize();

            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                ::AllocateSharedPtr(nel);
            
            SetCoeffPhys();

            SetupBoundaryConditions(HomoBasis,lhom,bcs,variable);
        }

        DisContField3DHomogeneous1D::~DisContField3DHomogeneous1D()
        {
        }


        void DisContField3DHomogeneous1D::SetupBoundaryConditions(const LibUtilities::BasisKey &HomoBasis, const NekDouble lhom, SpatialDomains::BoundaryConditions &bcs, const std::string variable)
        {

            int i,j,n;
            // Setup an ExpList2DHomogeneous1D expansion for boundary
            // conditions and link to class declared in m_planes.
            const SpatialDomains::BoundaryRegionCollection  &bregions = bcs.GetBoundaryRegions();
            const SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            int nbnd = bregions.size();
            
            // count the number of non-periodic boundary regions
            int cnt = 0;
            for(i = 0; i < nbnd; ++i)
            {
                SpatialDomains::BoundaryConditionShPtr boundaryCondition = GetBoundaryCondition(bconditions, i, variable);
                if( boundaryCondition->GetBoundaryConditionType() != SpatialDomains::ePeriodic )
                {
                    cnt++;
                }              
            }

            m_bndCondExpansions  = Array<OneD,MultiRegions::ExpListSharedPtr>(cnt);
            m_bndConditions = m_planes[0]->UpdateBndConditions();

            int nplanes = m_planes.num_elements();
            Array<OneD, MultiRegions::ExpListSharedPtr> PlanesBndCondExp(nplanes);
            
            cnt = 0;
            for(i = 0; i < nbnd; ++i)
            {
                SpatialDomains::BoundaryConditionShPtr boundaryCondition = GetBoundaryCondition(bconditions, i, variable);
                if(boundaryCondition->GetBoundaryConditionType() != SpatialDomains::ePeriodic)
                {
                    
                    boost::shared_ptr<StdRegions::StdExpansionVector> exp = MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr();
                    
                    for(n = 0; n < nplanes; ++n)
                    {
                        PlanesBndCondExp[n] = m_planes[n]->UpdateBndCondExpansion(i);

                        for(j = 0; j < PlanesBndCondExp[n]->GetExpSize(); ++j)
                        {
                            (*exp).push_back(PlanesBndCondExp[n]->GetExp(j));
                        }
                    }
                
                    m_bndCondExpansions[cnt++] = MemoryManager<ExpList2DHomogeneous1D>::AllocateSharedPtr(m_session,HomoBasis,lhom,m_useFFT,false,exp,PlanesBndCondExp);
                }
            }            
            EvaluateBoundaryConditions();
        }

        void DisContField3DHomogeneous1D::EvaluateBoundaryConditions(const NekDouble time)
        {
            int n;
            const Array<OneD, const NekDouble> z = m_homogeneousBasis->GetZ();
            Array<OneD, NekDouble> local_z(m_planes.num_elements());
            
            for(n = 0; n < m_planes.num_elements(); n++)
            {
                local_z[n] = z[m_transposition->GetPlaneID(n)];
            }
            
            for(n = 0; n < m_planes.num_elements(); ++n)
            {
                m_planes[n]->EvaluateBoundaryConditions(time,0.5*m_lhom*(1.0+local_z[n]));
            }
            
            // Fourier transform coefficient space boundary values
            // This will only be undertaken for time dependent
            // boundary conditions unless time == 0.0 which is the
            // case when the method is called from the constructor.
            for(n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                if(time == 0.0 || m_bndConditions[n]->GetUserDefined() == 
                   SpatialDomains::eTimeDependent)
                {
                    m_bndCondExpansions[n]->HomogeneousFwdTrans(m_bndCondExpansions[n]->GetCoeffs(),m_bndCondExpansions[n]->UpdateCoeffs());
                }
            }
        }
        
        void DisContField3DHomogeneous1D::v_HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const FlagList &flags,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff,
                const Array<OneD, const NekDouble> &dirForcing)
        {
            int n;
            int cnt = 0;
            int cnt1 = 0;
            NekDouble beta;
            StdRegions::ConstFactorMap new_factors;

            Array<OneD, NekDouble> e_out;
            Array<OneD, NekDouble> fce(inarray.num_elements());

            // Transform forcing function in half-physical space if required
            if(m_WaveSpace)
            {
                fce = inarray;
            }
            else 
            {
                HomogeneousFwdTrans(inarray,fce);
            }
            
            for(n = 0; n < m_planes.num_elements(); ++n)
            {
                if(n != 1 || m_transposition->GetK(n) != 0)
                {
                    beta = 2*M_PI*(m_transposition->GetK(n))/m_lhom;
                    new_factors = factors;
                    // add in Homogeneous Fourier direction and SVV if turned on
                    new_factors[StdRegions::eFactorLambda] += beta*beta*(1+GetSpecVanVisc(n));
                    m_planes[n]->HelmSolve(fce + cnt,
                                           e_out = outarray + cnt1,
                                           flags, new_factors, varcoeff, dirForcing);
                }
                
                cnt  += m_planes[n]->GetTotPoints();
                cnt1 += m_planes[n]->GetNcoeffs();
            }
        }
		
        void DisContField3DHomogeneous1D::v_EvaluateBoundaryConditions(const NekDouble time,const NekDouble x2_in,const NekDouble x3_in)
        {
            EvaluateBoundaryConditions(time);
        }
	
        boost::shared_ptr<ExpList> &DisContField3DHomogeneous1D::v_UpdateBndCondExpansion(int i)
        {
            return UpdateBndCondExpansion(i);
        }
	
        Array<OneD, SpatialDomains::BoundaryConditionShPtr> &DisContField3DHomogeneous1D::v_UpdateBndConditions()
        {
            return UpdateBndConditions();
        }
	
        void DisContField3DHomogeneous1D::GetBoundaryToElmtMap(Array<OneD, int> &ElmtID, Array<OneD,int> &EdgeID)
        {
            
            if(m_BCtoElmMap.num_elements() == 0)
            {
                Array<OneD, int> ElmtID_tmp;
                Array<OneD, int> EdgeID_tmp;
		
                m_planes[0]->GetBoundaryToElmtMap(ElmtID_tmp,EdgeID_tmp);
                int nel_per_plane = m_planes[0]->GetExpSize();
                int nplanes = m_planes.num_elements();
		
                int MapSize = ElmtID_tmp.num_elements();
		
                ElmtID = Array<OneD, int>(nplanes*MapSize);
                EdgeID = Array<OneD, int>(nplanes*MapSize);

                // If this mesh (or partition) has no BCs, skip this step.
                if (MapSize > 0)
                {
                    for(int i = 0; i < nplanes; i++)
                    {
                        for(int j = 0; j < MapSize; j++)
                        {
                            ElmtID[j+i*MapSize] = ElmtID_tmp[j] + i*nel_per_plane;
                            EdgeID[j+i*MapSize] = EdgeID_tmp[j];
                        }
                    }
                    m_BCtoElmMap = Array<OneD, int>(nplanes*MapSize);
                    m_BCtoEdgMap = Array<OneD, int>(nplanes*MapSize);

                    Vmath::Vcopy(nplanes*MapSize,ElmtID,1,m_BCtoElmMap,1);
                    Vmath::Vcopy(nplanes*MapSize,EdgeID,1,m_BCtoEdgMap,1);
                }
            }
            else 
            {
                int MapSize = m_BCtoElmMap.num_elements();
		
                ElmtID = Array<OneD, int>(MapSize);
                EdgeID = Array<OneD, int>(MapSize);
		
                Vmath::Vcopy(MapSize,m_BCtoElmMap,1,ElmtID,1);
                Vmath::Vcopy(MapSize,m_BCtoEdgMap,1,EdgeID,1);
            }			
        }
	
        void DisContField3DHomogeneous1D::GetBCValues(Array<OneD, NekDouble> &BndVals, 
                                                      const Array<OneD, NekDouble> &TotField, 
                                                      int BndID)
        {
            StdRegions::StdExpansionSharedPtr elmt;
            StdRegions::StdExpansion1DSharedPtr temp_BC_exp;
            
			Array<OneD, const NekDouble> tmp_Tot;
			Array<OneD, NekDouble> tmp_BC;
						
			int cnt = 0;
			int pos = 0;
			int exp_size, exp_size_per_plane, elmtID, boundaryID, offset, exp_dim;
			
			for(int k = 0; k < m_planes.num_elements(); k++)
			{
				for(int n = 0; n < m_bndConditions.num_elements(); ++n)
				{
					exp_size = m_bndCondExpansions[n]->GetExpSize();
					exp_size_per_plane = exp_size/m_planes.num_elements();
								
					for(int i = 0; i < exp_size_per_plane; i++)
					{
						if(n == BndID)
						{
							elmtID = m_BCtoElmMap[cnt];
							boundaryID = m_BCtoEdgMap[cnt];
							exp_dim = m_bndCondExpansions[n]->GetExp(i+k*exp_size_per_plane)->GetTotPoints();
							offset = GetPhys_Offset(elmtID);
							elmt = GetExp(elmtID);
							temp_BC_exp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_bndCondExpansions[n]->GetExp(i+k*exp_size_per_plane));
							elmt->GetEdgePhysVals(boundaryID,temp_BC_exp,tmp_Tot = TotField + offset,tmp_BC = BndVals + pos);
							pos += exp_dim;
						}
						cnt++;
					}
				}
			}
        }
	
        void DisContField3DHomogeneous1D::NormVectorIProductWRTBase(Array<OneD, const NekDouble> &V1,
                                                                    Array<OneD, const NekDouble> &V2,
                                                                    Array<OneD, NekDouble> &outarray,
                                                                    int BndID)
        {
            StdRegions::StdExpansionSharedPtr elmt;
            StdRegions::StdExpansion1DSharedPtr temp_BC_exp;
            
            Array<OneD, NekDouble> tmp_V1;
            Array<OneD, NekDouble> tmp_V2;
            Array<OneD, NekDouble> tmp_outarray;
            
            int cnt = 0;
            int exp_size, exp_size_per_plane, elmtID, Phys_offset, Coef_offset;
            
            for(int k = 0; k < m_planes.num_elements(); k++)
            {
                for(int n = 0; n < m_bndConditions.num_elements(); ++n)
                {
                    exp_size = m_bndCondExpansions[n]->GetExpSize();
                    exp_size_per_plane = exp_size/m_planes.num_elements();
                    
                    for(int i = 0; i < exp_size_per_plane; i++)
                    {
                        if(n == BndID)
                        {
                            elmtID = m_BCtoElmMap[cnt];
                            
                            Phys_offset = m_bndCondExpansions[n]->GetPhys_Offset(i+k*exp_size_per_plane);
                            Coef_offset = m_bndCondExpansions[n]->GetCoeff_Offset(i+k*exp_size_per_plane);
                            
                            elmt = GetExp(elmtID);
                            temp_BC_exp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_bndCondExpansions[n]->GetExp(i+k*exp_size_per_plane));
                            
                            temp_BC_exp->NormVectorIProductWRTBase(tmp_V1 = V1 + Phys_offset,tmp_V2 = V2 + Phys_offset,tmp_outarray = outarray + Coef_offset);
                        }
                        cnt++;
                    }
                }
            }
        }
    } // end of namespace
} //end of namespace
