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
            
            bool False = false;
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

        ContField3DHomogeneous2D::ContField3DHomogeneous2D(LibUtilities::CommSharedPtr &pComm,
                                                           const LibUtilities::BasisKey &HomoBasis_y,
														   const LibUtilities::BasisKey &HomoBasis_z,
														   const NekDouble lhom_y,
														   const NekDouble lhom_z,
														   bool useFFT,
														   SpatialDomains::MeshGraph1D &graph1D,
														   SpatialDomains::BoundaryConditions &bcs,
														   const int bc_loc,
														   const GlobalSysSolnType solnType):
            DisContField3DHomogeneous2D(pComm,HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,useFFT)
        {
            int i,j,n,nel;
            bool False = false;
            ContField1DSharedPtr line_zero;

            m_lines[0] = line_zero = MemoryManager<ContField1D>::AllocateSharedPtr(pComm,graph1D,bcs,bc_loc,solnType); 

            m_exp = MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr();
            nel = m_lines[0]->GetExpSize();

            for(j = 0; j < nel; ++j)
            {
                (*m_exp).push_back(m_lines[0]->GetExp(j));
            }

			int nylines = m_homogeneousBasis_y->GetNumPoints();
			int nzlines = m_homogeneousBasis_z->GetNumPoints();

            for(n = 1; n < nylines*nzlines; ++n)
            {
                m_lines[n] = MemoryManager<ContField1D>::AllocateSharedPtr(pComm,graph1D,bcs,bc_loc,solnType); 
                
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

            SetupBoundaryConditions(HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,bcs);
        }


        void ContField3DHomogeneous2D::SetCoeffPhys(void)
        {
            int i,n,cnt;
            int contncoeffs_per_line = m_lines[0]->GetContNcoeffs();
            int nyzlines = m_lines.num_elements();

            ExpList3DHomogeneous2D::SetCoeffPhys();

            // Set total coefficients and points
            m_contNcoeffs = contncoeffs_per_line*nyzlines;
            
            m_contCoeffs = Array<OneD, NekDouble> (m_contNcoeffs);

            Array<OneD, NekDouble> tmparray;

            for(cnt  = n = 0; n < nyzlines; ++n)
            {
                m_lines[n]->SetContCoeffsArray(tmparray= m_contCoeffs + contncoeffs_per_line*n);
            }
        }

        const Array<OneD, const NekDouble> &ContField3DHomogeneous2D::v_GetContCoeffs(void) const
        {
            return m_contCoeffs;
        }
       

        Array<OneD, NekDouble> &ContField3DHomogeneous2D::v_UpdateContCoeffs(void)
        {
            return m_contCoeffs;
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
                          NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff)
        {
            v_HelmSolveCG(inarray, outarray, lambda, varLambda, varCoeff,
                              false, NullNekDouble1DArray);
        }


        void ContField3DHomogeneous2D::v_HelmSolveCG(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,       NekDouble> &outarray,
                    NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff,
                    bool UseContCoeffs,
                    const Array<OneD, const NekDouble> &dirForcing)
        {
            int n,m;
            int cnt = 0;
            int cnt1 = 0;
			
			NekDouble beta_y;
			NekDouble beta_z;
			NekDouble beta; 
            Array<OneD, NekDouble> e_out;
            Array<OneD, NekDouble> fce(inarray.num_elements());
			
			int npoints_per_line = m_lines[0]->GetTotPoints();
			int ncoeffs_per_line;
			
			if(UseContCoeffs)
			{
				ncoeffs_per_line = m_lines[0]->GetContNcoeffs();
			}
			else 
			{
				ncoeffs_per_line = m_lines[0]->GetNcoeffs();
			}

            // Fourier transform forcing function
			Homogeneous2DFwdTrans(inarray,fce,UseContCoeffs);

            for(n = 0; n < m_nz; ++n)
            {
				for(m = 0; m < m_ny; ++m)
				{
					beta_z = 2*M_PI*(n/2)/m_lhom_z;
					beta_y = 2*M_PI*(m/2)/m_lhom_y;
					
					beta = beta_y*beta_y + beta_z*beta_z;
					
					m_lines[m+(n*m_ny)]->HelmSolve(fce + cnt,
												   e_out = outarray + cnt1,
												   lambda + beta,UseContCoeffs,
												   dirForcing, varLambda,varCoeff);
                
					cnt  += npoints_per_line;
					cnt1 += ncoeffs_per_line;
				}
			}
        }

        
    } // end of namespace
} //end of namespace
