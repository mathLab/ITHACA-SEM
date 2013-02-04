///////////////////////////////////////////////////////////////////////////////
//
// File ExpList3DHomogeneous1D.cpp
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
// Description: An ExpList which is homogeneous in 1 direction and so
// uses much of the functionality from a ExpList2D and its daughters
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declaration for typedefs
        ExpList3DHomogeneous1D::ExpList3DHomogeneous1D():
            ExpListHomogeneous1D()
        {
        }

        ExpList3DHomogeneous1D::ExpList3DHomogeneous1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                                       const LibUtilities::BasisKey &HomoBasis,
                                                       const NekDouble lhom,
                                                       const bool useFFT,
													   const bool dealiasing):
            ExpListHomogeneous1D(pSession,HomoBasis,lhom,useFFT,dealiasing)
        {
        }

        // Constructor for ExpList3DHomogeneous1D to act as a Explist2D field
        ExpList3DHomogeneous1D::ExpList3DHomogeneous1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                                       const LibUtilities::BasisKey &HomoBasis,
                                                       const NekDouble lhom,
                                                       const bool useFFT,
													   const bool dealiasing,
                                                       const SpatialDomains::MeshGraphSharedPtr &graph2D,
                                                       const std::string &var):
            ExpListHomogeneous1D(pSession,HomoBasis,lhom,useFFT,dealiasing)
        {
            GenExpList3DHomogeneous1D(graph2D->GetExpansions(var));
        }
        
        // Constructor for ExpList3DHomogeneous1D to act as a Explist2D field
        ExpList3DHomogeneous1D::ExpList3DHomogeneous1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                                       const LibUtilities::BasisKey &HomoBasis,
                                                       const NekDouble lhom,
                                                       const bool useFFT,
													   const bool dealiasing,
                                                       const SpatialDomains::ExpansionMap  &expansions):                                                      
            ExpListHomogeneous1D(pSession,HomoBasis,lhom,useFFT,dealiasing)
        {
            GenExpList3DHomogeneous1D(expansions);
        }

        void ExpList3DHomogeneous1D::GenExpList3DHomogeneous1D(const SpatialDomains::ExpansionMap &expansions)
        {
            int  n,j,nel;
            bool False = false;
            ExpList2DSharedPtr plane_zero;

            // note that nzplanes can be larger than nzmodes
            m_planes[0] = plane_zero = MemoryManager<ExpList2D>::AllocateSharedPtr(m_session, expansions, False);

            m_exp = MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr();
            nel = m_planes[0]->GetExpSize();
            
            for(j = 0; j < nel; ++j)
            {
                (*m_exp).push_back(m_planes[0]->GetExp(j));
            }
            
            for(n = 1; n < m_homogeneousBasis->GetNumPoints(); ++n)
            {
                m_planes[n] = MemoryManager<ExpList2D>::AllocateSharedPtr(*plane_zero,False);
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
        }


        /**
         * @param   In          ExpList3DHomogeneous1D object to copy.
         */
        ExpList3DHomogeneous1D::ExpList3DHomogeneous1D(const ExpList3DHomogeneous1D &In, bool DeclarePlanesSetCoeffPhys):
            ExpListHomogeneous1D(In)
        {
            if(DeclarePlanesSetCoeffPhys)
            {
                bool False = false;
                ExpList2DSharedPtr zero_plane = boost::dynamic_pointer_cast<ExpList2D> (In.m_planes[0]);

                for(int n = 0; n < m_planes.num_elements(); ++n)
                {
                    m_planes[n] = MemoryManager<ExpList2D>::AllocateSharedPtr(*zero_plane,False);
                }

                SetCoeffPhys();
            }
        }

        /**
         * Destructor
         */
        ExpList3DHomogeneous1D::~ExpList3DHomogeneous1D()
        {
        }

        void ExpList3DHomogeneous1D::SetCoeffPhys(void)
        {
            int i,n,cnt;
            int ncoeffs_per_plane = m_planes[0]->GetNcoeffs();
            int npoints_per_plane = m_planes[0]->GetTotPoints();

            int nzplanes = m_planes.num_elements();

            // Set total coefficients and points
            m_ncoeffs = ncoeffs_per_plane*nzplanes;
            m_npoints = npoints_per_plane*nzplanes;

            m_coeffs = Array<OneD, NekDouble> (m_ncoeffs,0.0);
            m_phys   = Array<OneD, NekDouble> (m_npoints,0.0);

            int nel = m_planes[0]->GetExpSize();
            m_coeff_offset   = Array<OneD,int>(nel*nzplanes);
            m_phys_offset    = Array<OneD,int>(nel*nzplanes);
            m_offset_elmt_id = Array<OneD,int>(nel*nzplanes);
            Array<OneD, NekDouble> tmparray;

            for(cnt  = n = 0; n < nzplanes; ++n)
            {
                m_planes[n]->SetCoeffsArray(tmparray= m_coeffs + ncoeffs_per_plane*n);
                m_planes[n]->SetPhysArray(tmparray = m_phys + npoints_per_plane*n);

                for(i = 0; i < nel; ++i)
                {
                    m_coeff_offset[cnt] = m_planes[n]->GetCoeff_Offset(i) + n*ncoeffs_per_plane;
                    m_phys_offset[cnt] =  m_planes[n]->GetPhys_Offset(i) + n*npoints_per_plane;
                    m_offset_elmt_id[cnt++] = m_planes[n]->GetOffset_Elmt_Id(i) + n*nel;
                }
            }
        }

        void ExpList3DHomogeneous1D::GetCoords(const int eid,
                                               Array<OneD, NekDouble> &xc0,
                                               Array<OneD, NekDouble> &xc1,
                                               Array<OneD, NekDouble> &xc2)
        {
            int n;
            Array<OneD, NekDouble> tmp_xc;
            int nzplanes = m_planes.num_elements();
            int npoints  = GetTotPoints(eid);

            (*m_exp)[eid]->GetCoords(xc0,xc1);

            // Fill z-direction
            Array<OneD, const NekDouble> pts =  m_homogeneousBasis->GetZ();
            Array<OneD, NekDouble> local_pts(m_planes.num_elements());
            
            for(n = 0; n < m_planes.num_elements(); n++)
            {
                local_pts[n] = pts[m_transposition->GetPlaneID(n)];
            }
            
            Array<OneD, NekDouble> z(nzplanes);

            Vmath::Smul(nzplanes,m_lhom/2.0,local_pts,1,z,1);
            Vmath::Sadd(nzplanes,m_lhom/2.0,z,1,z,1);

            for(n = 0; n < nzplanes; ++n)
            {
                Vmath::Fill(npoints,z[n],tmp_xc = xc2 + npoints*n,1);
                if(n)
                {
                    Vmath::Vcopy(npoints,xc0,1,tmp_xc = xc0+npoints*n,1);
                    Vmath::Vcopy(npoints,xc1,1,tmp_xc = xc1+npoints*n,1);
                }
            }
        }

        /**
         * The operation calls the 2D plane coordinates through the
         * function ExpList#GetCoords and then evaluated the third
         * coordinate using the member \a m_lhom
         *
         * @param coord_0 After calculation, the \f$x_1\f$ coordinate
         *                          will be stored in this array.
         *
         * @param coord_1 After calculation, the \f$x_2\f$ coordinate
         *                          will be stored in this array.
         *
         * @param coord_2 After calculation, the \f$x_3\f$ coordinate
         *                          will be stored in this array. This
         *                          coordinate is evaluated using the
         *                          predefined value \a m_lhom
         */
        void ExpList3DHomogeneous1D::v_GetCoords(Array<OneD, NekDouble> &xc0,
                                                 Array<OneD, NekDouble> &xc1,
                                                 Array<OneD, NekDouble> &xc2)
        {
            int n;
            Array<OneD, NekDouble> tmp_xc;
            int nzplanes = m_planes.num_elements();
            int npoints = m_planes[0]->GetTotPoints();

            m_planes[0]->GetCoords(xc0,xc1);

            // Fill z-direction
            Array<OneD, const NekDouble> pts =  m_homogeneousBasis->GetZ();
			
			Array<OneD, NekDouble> local_pts(m_planes.num_elements());
			
			for(n = 0; n < m_planes.num_elements(); n++)
			{
				local_pts[n] = pts[m_transposition->GetPlaneID(n)];
			}
			
            Array<OneD, NekDouble> z(nzplanes);

            Vmath::Smul(nzplanes,m_lhom/2.0,local_pts,1,z,1);
            Vmath::Sadd(nzplanes,m_lhom/2.0,z,1,z,1);

            for(n = 0; n < nzplanes; ++n)
            {
                Vmath::Fill(npoints,z[n],tmp_xc = xc2 + npoints*n,1);
                if(n)
                {
                    Vmath::Vcopy(npoints,xc0,1,tmp_xc = xc0+npoints*n,1);
                    Vmath::Vcopy(npoints,xc1,1,tmp_xc = xc1+npoints*n,1);
                }
            }
        }


        /**
         * Write Tecplot Files Zone
         * @param   outfile    Output file name.
         * @param   expansion  Expansion that is considered
         */
        void ExpList3DHomogeneous1D::v_WriteTecplotZone(std::ofstream &outfile, int expansion)
        {
            int i,j;

            int nquad0 = (*m_exp)[expansion]->GetNumPoints(0);
            int nquad1 = (*m_exp)[expansion]->GetNumPoints(1);
            int nquad2 = m_planes.num_elements();

            Array<OneD,NekDouble> coords[3];

            coords[0] = Array<OneD,NekDouble>(3*nquad0*nquad1*nquad2);
            coords[1] = coords[0] + nquad0*nquad1*nquad2;
            coords[2] = coords[1] + nquad0*nquad1*nquad2;

            GetCoords(expansion,coords[0],coords[1],coords[2]);

            outfile << "Zone, I=" << nquad0 << ", J=" << nquad1 <<",K="
                    << nquad2 << ", F=Block" << std::endl;

            for(j = 0; j < 3; ++j)
            {
                for(i = 0; i < nquad0*nquad1*nquad2; ++i)
                {
                    outfile << coords[j][i] << " ";
                }
                outfile << std::endl;
            }
        }


        void ExpList3DHomogeneous1D::v_WriteVtkPieceHeader(std::ofstream &outfile, int expansion)
        {
            int i,j,k;
            int nquad0 = (*m_exp)[expansion]->GetNumPoints(0);
            int nquad1 = (*m_exp)[expansion]->GetNumPoints(1);
            int nquad2 = m_planes.num_elements();
            int ntot = nquad0*nquad1*nquad2;
            int ntotminus = (nquad0-1)*(nquad1-1)*(nquad2-1);

            Array<OneD,NekDouble> coords[3];
            coords[0] = Array<OneD,NekDouble>(ntot);
            coords[1] = Array<OneD,NekDouble>(ntot);
            coords[2] = Array<OneD,NekDouble>(ntot);
            GetCoords(expansion,coords[0],coords[1],coords[2]);

            outfile << "    <Piece NumberOfPoints=\""
                    << ntot << "\" NumberOfCells=\""
                    << ntotminus << "\">" << endl;
            outfile << "      <Points>" << endl;
            outfile << "        <DataArray type=\"Float32\" "
                    << "NumberOfComponents=\"3\" format=\"ascii\">" << endl;
            outfile << "          ";
            for (i = 0; i < ntot; ++i)
            {
                for (j = 0; j < 3; ++j)
                {
                    outfile << coords[j][i] << " ";
                }
                outfile << endl;
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "      </Points>" << endl;
            outfile << "      <Cells>" << endl;
            outfile << "        <DataArray type=\"Int32\" "
                    << "Name=\"connectivity\" format=\"ascii\">" << endl;
            for (i = 0; i < nquad0-1; ++i)
            {
                for (j = 0; j < nquad1-1; ++j)
                {
                    for (k = 0; k < nquad2-1; ++k)
                    {
                        outfile << k*nquad0*nquad1 + j*nquad0 + i << " "
                                << k*nquad0*nquad1 + j*nquad0 + i + 1 << " "
                                << k*nquad0*nquad1 + (j+1)*nquad0 + i + 1 << " "
                                << k*nquad0*nquad1 + (j+1)*nquad0 + i << " "
                                << (k+1)*nquad0*nquad1 + j*nquad0 + i << " "
                                << (k+1)*nquad0*nquad1 + j*nquad0 + i + 1 << " "
                                << (k+1)*nquad0*nquad1 + (j+1)*nquad0 + i + 1 << " "
                                << (k+1)*nquad0*nquad1 + (j+1)*nquad0 + i << " " << endl;
                    }
                }
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "        <DataArray type=\"Int32\" "
                    << "Name=\"offsets\" format=\"ascii\">" << endl;
            for (i = 0; i < ntotminus; ++i)
            {
                outfile << i*8+8 << " ";
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "        <DataArray type=\"UInt8\" "
                    << "Name=\"types\" format=\"ascii\">" << endl;
            for (i = 0; i < ntotminus; ++i)
            {
                outfile << "12 ";
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "      </Cells>" << endl;
            outfile << "      <PointData>" << endl;
        }


        NekDouble ExpList3DHomogeneous1D::v_L2(const Array<OneD, const NekDouble> &soln)
        {
            int cnt = 0;
            NekDouble errL2,err = 0.0;
            Array<OneD, const NekDouble> w = m_homogeneousBasis->GetW();
            Array<OneD, NekDouble> local_w(m_planes.num_elements());
            
            for(int n = 0; n < m_planes.num_elements(); n++)
            {
                local_w[n] = w[m_transposition->GetPlaneID(n)];
            }
            
            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                errL2 = m_planes[n]->L2(soln + cnt);
                cnt  += m_planes[n]->GetTotPoints();
                err  += errL2*errL2*local_w[n]*m_lhom*0.5;
            }
            
            m_comm->GetColumnComm()->AllReduce(err, LibUtilities::ReduceSum);
            
            return sqrt(err);
        }
	
        NekDouble ExpList3DHomogeneous1D::v_L2(void)
        {
            NekDouble errL2,err = 0;
            Array<OneD, const NekDouble> w = m_homogeneousBasis->GetW();
            Array<OneD, NekDouble> local_w(m_planes.num_elements());
            
            for(int n = 0; n < m_planes.num_elements(); n++)
            {
                local_w[n] = w[m_transposition->GetPlaneID(n)];
            }

            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                errL2 = m_planes[n]->L2();
                err += errL2*errL2*local_w[n]*m_lhom*0.5;
            }
            
            m_comm->GetColumnComm()->AllReduce(err, LibUtilities::ReduceSum);
                        
            return sqrt(err);
        }
		
        Array<OneD, const NekDouble> ExpList3DHomogeneous1D::v_HomogeneousEnergy(void)
        {
            Array<OneD, NekDouble> energy(m_planes.num_elements()/2);
            NekDouble area = 0.0;

            // Calculate total area of elements.
            for (int n = 0; n < m_planes[0]->GetExpSize(); ++n)
            {
                Array<OneD, NekDouble> inarray(m_planes[0]->GetExp(n)->GetTotPoints(), 1.0);
                area += m_planes[0]->GetExp(n)->Integral(inarray);
            }
            
            m_comm->GetRowComm()->AllReduce(area, LibUtilities::ReduceSum);
            
            // Calculate L2 norm of real/imaginary planes.
            for (int n = 0; n < m_planes.num_elements(); n += 2)
            {
                NekDouble err;
                
                energy[n/2] = 0;
                
                for(int i = 0; i < m_planes[n]->GetExpSize(); ++i)
                {
                    StdRegions::StdExpansionSharedPtr exp = m_planes[n]->GetExp(i);
                    exp->BwdTrans(m_planes[n]->GetCoeffs()+m_planes[n]->GetCoeff_Offset(i),
                                  exp->UpdatePhys());
                    err = exp->L2();
                    energy[n/2] += err*err;
                    
                    exp = m_planes[n+1]->GetExp(i);
                    exp->BwdTrans(m_planes[n+1]->GetCoeffs()+m_planes[n+1]->GetCoeff_Offset(i),
                                  exp->UpdatePhys());
                    err = exp->L2();
                    energy[n/2] += err*err;
                }
                
                m_comm->GetRowComm()->AllReduce(energy[n/2], LibUtilities::ReduceSum);
                energy[n/2] /= 2.0*area;
            }
            
            return energy;
        }
    } //end of namespace
} //end of namespace


