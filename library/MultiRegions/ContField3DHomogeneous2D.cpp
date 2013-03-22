//////////////////////////////////////////////////////////////////////////////
//
// File ContField3DHomogeneous2D.cpp
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
// conditions and a 2 homogeneous directions
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <MultiRegions/ContField1D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        ContField3DHomogeneous2D::ContField3DHomogeneous2D(void):
            DisContField3DHomogeneous2D()
        {
        }

        ContField3DHomogeneous2D::ContField3DHomogeneous2D(const ContField3DHomogeneous2D &In):
            DisContField3DHomogeneous2D (In,false)
        {
            
            ContField1DSharedPtr zero_line = boost::dynamic_pointer_cast<ContField1D> (In.m_lines[0]);
            
            for(int n = 0; n < m_lines.num_elements(); ++n)
            {
                m_lines[n] = MemoryManager<ContField1D>::AllocateSharedPtr(*zero_line);
            }
            
            SetCoeffPhys();
        }
        
        ContField3DHomogeneous2D::~ContField3DHomogeneous2D()
        {
        }

        ContField3DHomogeneous2D::ContField3DHomogeneous2D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                                           const LibUtilities::BasisKey &HomoBasis_y,
														   const LibUtilities::BasisKey &HomoBasis_z,
														   const NekDouble lhom_y,
														   const NekDouble lhom_z,
														   const bool useFFT,
														   const bool dealiasing,
														   const SpatialDomains::MeshGraphSharedPtr &graph1D,
														   const std::string &variable):
            DisContField3DHomogeneous2D(pSession,HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,useFFT,dealiasing)
        {
            int i,n,nel;
            ContField1DSharedPtr line_zero;
            SpatialDomains::BoundaryConditions bcs(pSession, graph1D);

            m_lines[0] = line_zero = MemoryManager<ContField1D>::AllocateSharedPtr(pSession,graph1D,variable);

            m_exp = MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr();
            nel = m_lines[0]->GetExpSize();

            for(i = 0; i < nel; ++i)
            {
                (*m_exp).push_back(m_lines[0]->GetExp(i));
            }

			int nylines = m_homogeneousBasis_y->GetNumPoints();
			int nzlines = m_homogeneousBasis_z->GetNumPoints();

            for(n = 1; n < nylines*nzlines; ++n)
            {
                m_lines[n] = MemoryManager<ContField1D>::AllocateSharedPtr(pSession,graph1D,variable);
                
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

            SetupBoundaryConditions(HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,bcs);
        }


        void ContField3DHomogeneous2D::v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray)
        {
            Array<OneD, NekDouble> tmp; 
            int ncoeffs = m_lines[0]->GetNcoeffs();

            for(int n = 0; n < m_lines.num_elements(); ++n)
            {
                m_lines[n]->ImposeDirichletConditions(tmp = outarray + 
                                                       n*ncoeffs);
            }
        }


        /**
         * 
         */
        void  ContField3DHomogeneous2D::v_LocalToGlobal(void) 
        {
            for(int n = 0; n < m_lines.num_elements(); ++n)
            {
                m_lines[n]->LocalToGlobal();
            }
        };


        /**
         * 
         */
        void  ContField3DHomogeneous2D::v_GlobalToLocal(void) 
        {
            for(int n = 0; n < m_lines.num_elements(); ++n)
            {
                m_lines[n]->GlobalToLocal();
            }
        };


        void ContField3DHomogeneous2D::v_HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const FlagList &flags,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff,
                const Array<OneD, const NekDouble> &dirForcing)
        {
            int n,m;
            int cnt = 0;
            int cnt1 = 0;
            int nhom_modes_y = m_homogeneousBasis_y->GetNumModes();
			int nhom_modes_z = m_homogeneousBasis_z->GetNumModes();
            NekDouble beta_y;
			NekDouble beta_z;
			NekDouble beta;
			StdRegions::ConstFactorMap new_factors;

            Array<OneD, NekDouble> e_out;
            Array<OneD, NekDouble> fce(inarray.num_elements());
			
			if(m_WaveSpace)
			{
				fce = inarray;
			}
			else 
			{
				// Fourier transform forcing function
                            HomogeneousFwdTrans(inarray,fce,(flags.isSet(eUseGlobal))?eGlobal:eLocal);
			}

            for(n = 0; n < nhom_modes_z; ++n)
            {
                for(m = 0; m < nhom_modes_y; ++m)
                {
                    beta_z = 2*M_PI*(n/2)/m_lhom_z;
                    beta_y = 2*M_PI*(m/2)/m_lhom_y;
                    beta = beta_y*beta_y + beta_z*beta_z;
                    new_factors = factors;
                    new_factors[StdRegions::eFactorLambda] += beta;
                    
                    m_lines[n]->HelmSolve(fce + cnt,
                                          e_out = outarray + cnt1,
                                          flags, new_factors, varcoeff, dirForcing);
                    
                    cnt  += m_lines[n]->GetTotPoints();
                    
                    cnt1 += m_lines[n]->GetNcoeffs();
                    
                }
            }
        }
        
    } // end of namespace
} //end of namespace
