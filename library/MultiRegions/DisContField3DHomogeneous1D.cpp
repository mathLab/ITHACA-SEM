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

        DisContField3DHomogeneous1D::DisContField3DHomogeneous1D(LibUtilities::CommSharedPtr &pComm,
                                                                 const LibUtilities::BasisKey &HomoBasis,
                                                                 const NekDouble lhom,
																 bool useFFT):
            ExpList3DHomogeneous1D(pComm,HomoBasis,lhom,useFFT),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }

        DisContField3DHomogeneous1D::DisContField3DHomogeneous1D(const DisContField3DHomogeneous1D &In, bool DeclarePlanesSetCoeffPhys):
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
                                       LibUtilities::CommSharedPtr &pComm,
                                       const LibUtilities::BasisKey &HomoBasis,
                                       const NekDouble lhom,
									   bool useFFT,
                                       SpatialDomains::MeshGraph2D &graph2D,
                                       SpatialDomains::BoundaryConditions &bcs,
                                       const int bc_loc,
                                       const GlobalSysSolnType solnType):  
            ExpList3DHomogeneous1D(pComm,HomoBasis,lhom,useFFT),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            int i,j,n,nel;
            bool True  = true; 
            bool False = false; 
            DisContField2DSharedPtr plane_zero;

            // note that nzplanes can be larger than nzmodes 
            m_planes[0] = plane_zero = MemoryManager<DisContField2D>::AllocateSharedPtr(pComm,graph2D,bcs,bc_loc,solnType,True,False);

            m_exp = MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr();
            nel = m_planes[0]->GetExpSize();

            for(j = 0; j < nel; ++j)
            {
                (*m_exp).push_back(m_planes[0]->GetExp(j));
            }

            for(n = 1; n < m_homogeneousBasis->GetNumPoints(); ++n)
            {
                m_planes[n] = MemoryManager<DisContField2D>::AllocateSharedPtr(*plane_zero,graph2D,bcs,bc_loc,True,False);
                for(j = 0; j < nel; ++j)
                {
                    (*m_exp).push_back((*m_exp)[j]);
                }
            }            

            // Setup Default optimisation information. 
            nel = GetExpSize();

            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                ::AllocateSharedPtr(nel);
            
            SetCoeffPhys();

            SetupBoundaryConditions(HomoBasis,lhom,bcs);
        }

        DisContField3DHomogeneous1D::~DisContField3DHomogeneous1D()
        {
        }


        void DisContField3DHomogeneous1D::SetupBoundaryConditions(const LibUtilities::BasisKey &HomoBasis, const NekDouble lhom, SpatialDomains::BoundaryConditions &bcs)
        {

            int i,j,n;
            // Setup an ExpList2DHomogeneous1D expansion for boundary
            // conditions and link to class declared in m_planes.
            SpatialDomains::BoundaryRegionCollection  &bregions = bcs.GetBoundaryRegions();
            int nbnd = bregions.size();
            
            m_bndCondExpansions  = Array<OneD,MultiRegions::ExpListSharedPtr>(nbnd);
            
            m_bndConditions = m_planes[0]->UpdateBndConditions();
			
			MultiRegions::ExpList2DHomogeneous1DSharedPtr locExpList;
            
            int nplanes = m_planes.num_elements();
            Array<OneD, MultiRegions::ExpListSharedPtr> PlanesBndCondExp(nplanes);
            
            for(i = 0; i < nbnd; ++i)
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
                
				m_bndCondExpansions[i] = MemoryManager<ExpList2DHomogeneous1D>::AllocateSharedPtr(m_comm,HomoBasis,lhom,m_useFFT,exp,PlanesBndCondExp);
               
				 //= locExpList;
            }
            
            EvaluateBoundaryConditions();
        }

        void DisContField3DHomogeneous1D::EvaluateBoundaryConditions(const NekDouble time)
        {
            int n;
            const Array<OneD, const NekDouble> z = m_homogeneousBasis->GetZ();

            for(n = 0; n < m_planes.num_elements(); ++n)
            {
                m_planes[n]->EvaluateBoundaryConditions(time,0.5*m_lhom*(1.0+z[n]));
            }
            
            // Fourier transform coefficient space boundary values
            for(n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                m_bndCondExpansions[n]->HomogeneousFwdTrans(m_bndCondExpansions[n]->GetCoeffs(),m_bndCondExpansions[n]->UpdateCoeffs());
            }    
        }
        
        void DisContField3DHomogeneous1D::v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,       NekDouble> &outarray,
                    NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff)
        {
            v_HelmSolveDG(inarray, outarray, lambda, varLambda, varCoeff, 1);
        }

        void DisContField3DHomogeneous1D::v_HelmSolveDG(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,       NekDouble> &outarray,
                    NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff,
                    NekDouble tau)
        {
            int n;
            int cnt = 0;
            int cnt1 = 0;
            int nhom_modes = m_homogeneousBasis->GetNumModes();
            NekDouble beta; 
            Array<OneD, NekDouble> e_out;
            Array<OneD, NekDouble> fce(inarray.num_elements());

            // Fourier transform forcing function
			if(m_FourierSpace != eCoef)
			{
				HomogeneousFwdTrans(inarray,fce);
			}
			
            for(n = 0; n < nhom_modes; ++n)
            {
                beta = 2*M_PI*(n/2)/m_lhom;
                m_planes[n]->HelmSolve(fce + cnt,
                                         e_out = outarray + cnt1,
                                         lambda + beta*beta,varLambda,
                                         varCoeff);
                
                cnt  += m_planes[n]->GetTotPoints();
                cnt1 += m_planes[n]->GetNcoeffs();
            }
        }
		
		void DisContField3DHomogeneous1D::v_EvaluateBoundaryConditions(const NekDouble time,const NekDouble x2_in,const NekDouble x3_in)
		{
			EvaluateBoundaryConditions(time);
		}
		
		const Array<OneD,const boost::shared_ptr<ExpList> > &DisContField3DHomogeneous1D::v_GetBndCondExpansions(void)
		{
			return GetBndCondExpansions();
		}
		
		const Array<OneD,const SpatialDomains::BoundaryConditionShPtr> &DisContField3DHomogeneous1D::v_GetBndConditions()
		{
			return GetBndConditions();
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
			Array<OneD, int> ElmtID_tmp;
			Array<OneD, int> EdgeID_tmp;
			
			m_planes[0]->GetBoundaryToElmtMap(ElmtID_tmp,EdgeID_tmp);
			
			int nel_per_plane = m_planes[0]->GetExpSize();
			
			int nplanes = m_planes.num_elements();
			int MapSize = ElmtID_tmp.num_elements();
			
			ElmtID = Array<OneD, int>(nplanes*MapSize);
			EdgeID = Array<OneD, int>(nplanes*MapSize);
			
			for(int i = 0; i < nplanes; i++)
			{
				for(int j = 0; j < MapSize; j++)
				{
					ElmtID[j+i*MapSize] = ElmtID_tmp[j] + i*nel_per_plane;
					EdgeID[j+i*MapSize] = EdgeID_tmp[j];
				}
			}
		}
		
		void DisContField3DHomogeneous1D::GetBCValues(Array<OneD, NekDouble> &BndVals, 
													  const Array<OneD, NekDouble> &TotField, 
													  int BndID)
		{
			Array<OneD, int> ElmtID;
			Array<OneD, int> EdgeID;
			StdRegions::StdExpansionSharedPtr elmt;
			StdRegions::StdExpansion1DSharedPtr temp_BC_exp;
			
			Array<OneD, const NekDouble> tmp_Tot;
			Array<OneD, NekDouble> tmp_BC;
			
			GetBoundaryToElmtMap(ElmtID,EdgeID);
			
			int cnt = 0;
			int exp_size, elmtID, boundaryID, offset, exp_dim, pos;
			
			for(int n = 0; n < m_bndConditions.num_elements(); ++n)
			{
				if(n == BndID)
				{
					exp_size = m_bndCondExpansions[n]->GetExpSize();
					pos = 0;
					for(int i = 0; i < exp_size; i++)
					{
						elmtID = ElmtID[cnt];
						boundaryID = EdgeID[cnt];
						
						exp_dim = m_bndCondExpansions[n]->GetExp(i)->GetTotPoints();
						offset = GetPhys_Offset(elmtID);
						
						elmt = GetExp(elmtID);
						
						temp_BC_exp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_bndCondExpansions[n]->GetExp(i));
		
						elmt->GetEdgePhysVals(boundaryID,temp_BC_exp,tmp_Tot = TotField + offset,tmp_BC = BndVals + pos);
						
						pos += exp_dim;
						cnt++;
					}
				}
				else
				{
					cnt += m_bndCondExpansions[n]->GetExpSize();
				}
			}
		}
		
		void DisContField3DHomogeneous1D::NormVectorIProductWRTBase(Array<OneD, const NekDouble> &V1,
																	Array<OneD, const NekDouble> &V2,
																	Array<OneD, NekDouble> &outarray,
																	int BndID)
		{
			Array<OneD, int> ElmtID;
			Array<OneD, int> EdgeID;
			StdRegions::StdExpansionSharedPtr elmt;
			StdRegions::StdExpansion1DSharedPtr temp_BC_exp;
			
			Array<OneD, const NekDouble> tmp_V1;
			Array<OneD, const NekDouble> tmp_V2;
			Array<OneD, NekDouble> tmp_outarray;
			
			GetBoundaryToElmtMap(ElmtID,EdgeID);
			
			bool NegateNormals;
			
			int cnt = 0;
			int exp_size, elmtID, boundaryID, offset;
			
			for(int n = 0; n < m_bndConditions.num_elements(); ++n)
			{
				if(n == BndID)
				{
					exp_size = m_bndCondExpansions[n]->GetExpSize();
					
					for(int i = 0; i < exp_size; i++)
					{
						elmtID = ElmtID[cnt];
						boundaryID = EdgeID[cnt];
						
						offset = m_bndCondExpansions[n]->GetPhys_Offset(i);
						
						elmt = GetExp(elmtID);
						
						temp_BC_exp =  boost::dynamic_pointer_cast<StdRegions::StdExpansion1D> (m_bndCondExpansions[n]->GetExp(i));
						
						// Decide if normals facing outwards
						NegateNormals = (elmt->GetEorient(boundaryID) == StdRegions::eForwards)? false:true;
						
						temp_BC_exp->NormVectorIProductWRTBase(tmp_V1 = V1 + offset,tmp_V2 = V2 + offset,tmp_outarray = outarray + offset,NegateNormals);
			
						cnt++;
					}
				}
				else
				{
					cnt += m_bndCondExpansions[n]->GetExpSize();
				}
			}
		}


    } // end of namespace
} //end of namespace
