//////////////////////////////////////////////////////////////////////////////
//
// File ContField3DHomogeneous1D.cpp
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
// conditions and a 1D homogeneous direction
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContField3DHomogeneous1D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        ContField3DHomogeneous1D::ContField3DHomogeneous1D(void):
            DisContField3DHomogeneous1D()
        {
        }

        ContField3DHomogeneous1D::ContField3DHomogeneous1D(const ContField3DHomogeneous1D &In):
            DisContField3DHomogeneous1D (In,false)
        {
            
            bool False = false;
            ContField2DSharedPtr zero_plane = boost::dynamic_pointer_cast<ContField2D> (In.m_planes[0]);
            
            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                m_planes[n] = MemoryManager<ContField2D>::AllocateSharedPtr(*zero_plane,False);
            }
            
            SetCoeffPhys();
        }
        
        ContField3DHomogeneous1D::~ContField3DHomogeneous1D()
        {
        }

        ContField3DHomogeneous1D::ContField3DHomogeneous1D(
                                       const LibUtilities::SessionReaderSharedPtr &pSession,
                                       const LibUtilities::BasisKey &HomoBasis,
                                       const NekDouble lhom,
									   const bool useFFT,
									   const bool dealiasing,
                                       const SpatialDomains::MeshGraphSharedPtr &graph2D,
                                       const std::string &variable,
									   const bool CheckIfSingularSystem):
            DisContField3DHomogeneous1D(pSession,HomoBasis,lhom,useFFT,dealiasing)
        {
            int i,n,nel;
            bool False = false;
            ContField2DSharedPtr plane_zero;
            SpatialDomains::BoundaryConditions bcs(m_session, graph2D);

            // note that nzplanes can be larger than nzmodes 
            m_planes[0] = plane_zero = MemoryManager<ContField2D>::AllocateSharedPtr(pSession,graph2D,variable,False,CheckIfSingularSystem);

            m_exp = MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr();
            nel = m_planes[0]->GetExpSize();

            for(i = 0; i < nel; ++i)
            {
                (*m_exp).push_back(m_planes[0]->GetExp(i));
            }

            for(n = 1; n < m_planes.num_elements(); ++n)
            {
                m_planes[n] = MemoryManager<ContField2D>::AllocateSharedPtr(*plane_zero,graph2D,variable,False,CheckIfSingularSystem);
                
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


        void ContField3DHomogeneous1D::SetCoeffPhys(void)
        {
            int n,cnt;
            int contncoeffs_per_plane = m_planes[0]->GetContNcoeffs();
            int nzplanes = m_planes.num_elements();

            ExpList3DHomogeneous1D::SetCoeffPhys();

            // Set total coefficients and points
            m_contNcoeffs = contncoeffs_per_plane*nzplanes;
            
            m_contCoeffs = Array<OneD, NekDouble> (m_contNcoeffs);

            Array<OneD, NekDouble> tmparray;

            for(cnt  = n = 0; n < nzplanes; ++n)
            {
                m_planes[n]->SetContCoeffsArray(tmparray= m_contCoeffs + contncoeffs_per_plane*n);
            }
        }

        const Array<OneD, const NekDouble> &ContField3DHomogeneous1D::v_GetContCoeffs(void) const
        {
            return m_contCoeffs;
        }
       

        Array<OneD, NekDouble> &ContField3DHomogeneous1D::v_UpdateContCoeffs(void)
        {
            return m_contCoeffs;
        }
       
 
        /**
         * 
         */
        void  ContField3DHomogeneous1D::v_LocalToGlobal(void) 
        {
            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                m_planes[n]->LocalToGlobal();
            }
        };


        /**
         * 
         */
        void  ContField3DHomogeneous1D::v_GlobalToLocal(void) 
        {
            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                m_planes[n]->GlobalToLocal();
            }
        };


        void ContField3DHomogeneous1D::v_HelmSolve(
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

            // Fourier transform forcing function
			if(m_WaveSpace)
			{
				fce = inarray;
			}
			else 
			{
				HomogeneousFwdTrans(inarray,fce,flags.isSet(eUseContCoeff));
			}

            for(n = 0; n < m_planes.num_elements(); ++n)
            {
				beta = 2*M_PI*(m_transposition->GetK(n))/m_lhom;
				new_factors = factors;
				new_factors[StdRegions::eFactorLambda] += beta*beta;

                m_planes[n]->HelmSolve(fce + cnt,
                                       e_out = outarray + cnt1,
                                       flags, new_factors, varcoeff, dirForcing);
                
                cnt  += m_planes[n]->GetTotPoints();
                if(flags.isSet(eUseContCoeff))
                {
                    cnt1 += m_planes[n]->GetContNcoeffs();
                }
                else
                {
                    cnt1 += m_planes[n]->GetNcoeffs();
                }
            }
        }
    } // end of namespace
} //end of namespace
