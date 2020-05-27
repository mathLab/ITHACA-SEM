///////////////////////////////////////////////////////////////////////////////
//
// File ExpList3DHomogeneous2D.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: An ExpList which is homogeneous in 2 directions and so
// uses much of the functionality from a ExpList1D and its daughters
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <MultiRegions/ExpList1D.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declaration for typedefs
        ExpList3DHomogeneous2D::ExpList3DHomogeneous2D():
            ExpListHomogeneous2D()
        {
            SetExpType(e3DH2D);
        }

        ExpList3DHomogeneous2D::ExpList3DHomogeneous2D(
                       const LibUtilities::SessionReaderSharedPtr &pSession,
                       const LibUtilities::BasisKey &HomoBasis_y,
                       const LibUtilities::BasisKey &HomoBasis_z,
                       const NekDouble lhom_y,
                       const NekDouble lhom_z,
                       const bool useFFT,
                       const bool dealiasing,
                       const Collections::ImplementationType ImpType):
            ExpListHomogeneous2D(pSession,HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,useFFT,dealiasing)
        {
            boost::ignore_unused(ImpType);
            SetExpType(e3DH2D);
        }

        // Constructor for ExpList3DHomogeneous2D to act as a Explist1D field
        ExpList3DHomogeneous2D::ExpList3DHomogeneous2D(
                    const LibUtilities::SessionReaderSharedPtr &pSession,
                    const LibUtilities::BasisKey &HomoBasis_y,
                    const LibUtilities::BasisKey &HomoBasis_z,
                    const NekDouble lhom_y,
                    const NekDouble lhom_z,
                    const bool useFFT,
                    const bool dealiasing,
                    const SpatialDomains::MeshGraphSharedPtr &graph1D,
                    const Collections::ImplementationType ImpType):
            ExpListHomogeneous2D(pSession,HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,useFFT,dealiasing)
        {
            SetExpType(e3DH2D);

            int n,j,nel;
            bool False = false;
            ExpList1DSharedPtr line_zero;

            //
            m_lines[0] = line_zero = MemoryManager<ExpList1D>::
                AllocateSharedPtr(m_session,graph1D, False,ImpType);

            m_exp = MemoryManager<LocalRegions::ExpansionVector>::
                AllocateSharedPtr();
            nel = m_lines[0]->GetExpSize();

            for(j = 0; j < nel; ++j)
            {
                (*m_exp).push_back(m_lines[0]->GetExp(j));
            }

            int ny = m_homogeneousBasis_y->GetNumPoints();
            int nz = m_homogeneousBasis_z->GetNumPoints();

            for(n = 1; n < (ny*nz); ++n)
            {
                m_lines[n] = MemoryManager<ExpList1D>::AllocateSharedPtr(*line_zero,False);
                for(j = 0; j < nel; ++j)
                {
                    (*m_exp).push_back((*m_exp)[j]);
                }
            }

            SetCoeffPhys();
        }


        /**
         * @param   In          ExpList3DHomogeneous2D object to copy.
         */
        ExpList3DHomogeneous2D::ExpList3DHomogeneous2D(
                            const ExpList3DHomogeneous2D &In,
                            const bool DeclareLinesSetCoeffPhys):
            ExpListHomogeneous2D(In)
        {
            SetExpType(e3DH2D);

            if(DeclareLinesSetCoeffPhys)
            {
                bool False = false;
                ExpList1DSharedPtr zero_line = std::dynamic_pointer_cast<ExpList1D> (In.m_lines[0]);

                for(int n = 0; n < m_lines.size(); ++n)
                {
                    m_lines[n] = MemoryManager<ExpList1D>::AllocateSharedPtr(*zero_line,False);
                }

                SetCoeffPhys();
            }
        }

        /**
         *
         */
        ExpList3DHomogeneous2D::ExpList3DHomogeneous2D(
                       const ExpList3DHomogeneous2D &In,
                       const std::vector<unsigned int> &eIDs,
                       const bool DeclareLinesSetCoeffPhys,
                       const Collections::ImplementationType ImpType):
            ExpListHomogeneous2D(In, eIDs)
        {
            SetExpType(e3DH2D);

            if(DeclareLinesSetCoeffPhys)
            {
                bool False = false;
                std::vector<unsigned int> eIDsLine;
                int nel = eIDs.size()/m_lines.size();

                for (int i = 0; i < nel; ++i)
                {
                    eIDsLine.push_back(eIDs[i]);
                }

                ExpList1DSharedPtr zero_line_old =
                        std::dynamic_pointer_cast<ExpList1D> (In.m_lines[0]);

                ExpList1DSharedPtr zero_line =
                    MemoryManager<ExpList1D>::AllocateSharedPtr(*(zero_line_old), eIDsLine, ImpType);

                for(int n = 0; n < m_lines.size(); ++n)
                {
                    m_lines[n] = MemoryManager<ExpList1D>::AllocateSharedPtr(*zero_line,False);
                }

                SetCoeffPhys();
            }
        }

        /**
         * Destructor
         */
        ExpList3DHomogeneous2D::~ExpList3DHomogeneous2D()
        {
        }

        void ExpList3DHomogeneous2D::SetCoeffPhys(void)
        {
            int i,n,cnt;
            int ncoeffs_per_line = m_lines[0]->GetNcoeffs();
            int npoints_per_line = m_lines[0]->GetTotPoints();

            int nyzlines = m_lines.size();

            // Set total coefficients and points
            m_ncoeffs = ncoeffs_per_line*nyzlines;
            m_npoints = npoints_per_line*nyzlines;

            m_coeffs = Array<OneD, NekDouble> (m_ncoeffs,0.0);
            m_phys   = Array<OneD, NekDouble> (m_npoints,0.0);

            int nel = m_lines[0]->GetExpSize();
            m_coeff_offset   = Array<OneD,int>(nel*nyzlines);
            m_phys_offset    = Array<OneD,int>(nel*nyzlines);
            Array<OneD, NekDouble> tmparray;

            for(cnt  = n = 0; n < nyzlines; ++n)
            {
                m_lines[n]->SetCoeffsArray(tmparray= m_coeffs + ncoeffs_per_line*n);
                m_lines[n]->SetPhysArray(tmparray = m_phys + npoints_per_line*n);

                for(i = 0; i < nel; ++i)
                {
                    m_coeff_offset[cnt] = m_lines[n]->GetCoeff_Offset(i) + n*ncoeffs_per_line;
                    m_phys_offset[cnt++] =  m_lines[n]->GetPhys_Offset(i) + n*npoints_per_line;
                }
            }
        }

        void ExpList3DHomogeneous2D::GetCoords(const int eid,
                                               Array<OneD, NekDouble> &xc0,
                                               Array<OneD, NekDouble> &xc1,
                                               Array<OneD, NekDouble> &xc2)
        {
            int n,m,j;
            Array<OneD, NekDouble> tmp_xc;
			int nylines = m_homogeneousBasis_y->GetNumPoints();
			int nzlines = m_homogeneousBasis_z->GetNumPoints();

			int npoints  = GetTotPoints(eid);

			// Fill x-y-z-direction
            Array<OneD, const NekDouble> pts_y =  m_homogeneousBasis_y->GetZ();
			Array<OneD, const NekDouble> pts_z =  m_homogeneousBasis_z->GetZ();

			Array<OneD, NekDouble> x(npoints);
			Array<OneD, NekDouble> y(nylines);
			Array<OneD, NekDouble> z(nzlines);

            Vmath::Smul(nylines,m_lhom_y/2.0,pts_y,1,y,1);
            Vmath::Sadd(nylines,m_lhom_y/2.0,y,1,y,1);

			Vmath::Smul(nzlines,m_lhom_z/2.0,pts_z,1,z,1);
            Vmath::Sadd(nzlines,m_lhom_z/2.0,z,1,z,1);

			(*m_exp)[eid]->GetCoords(x);


            for(m = 0; m < nzlines; ++m)
            {
				for(j = 0; j < nylines; ++j)
				{
					for(n = 0; n < npoints; ++n)
					{
						Vmath::Fill(1,x[n],tmp_xc = xc0 + n +(j*npoints) + (m*npoints*nylines), 1);
						Vmath::Fill(1,y[j],tmp_xc = xc1 + n +(j*npoints) + (m*npoints*nylines), 1);
						Vmath::Fill(1,z[m],tmp_xc = xc2 + n +(j*npoints) + (m*npoints*nylines), 1);
					}
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
        void ExpList3DHomogeneous2D::v_GetCoords(Array<OneD, NekDouble> &xc0,
                                                 Array<OneD, NekDouble> &xc1,
                                                 Array<OneD, NekDouble> &xc2)
        {
            int n,m,j;
            Array<OneD, NekDouble> tmp_xc;
            int npoints = m_lines[0]->GetTotPoints();

            int nylines = m_homogeneousBasis_y->GetNumPoints();
			int nzlines = m_homogeneousBasis_z->GetNumPoints();

            // Fill z-direction
            Array<OneD, const NekDouble> pts_y =  m_homogeneousBasis_y->GetZ();
			Array<OneD, const NekDouble> pts_z =  m_homogeneousBasis_z->GetZ();

            Array<OneD, NekDouble> x(npoints);
			Array<OneD, NekDouble> y(nylines);
			Array<OneD, NekDouble> z(nzlines);

			m_lines[0]->GetCoords(x);

            Vmath::Smul(nylines,m_lhom_y/2.0,pts_y,1,y,1);
            Vmath::Sadd(nylines,m_lhom_y/2.0,y,1,y,1);

			Vmath::Smul(nzlines,m_lhom_z/2.0,pts_z,1,z,1);
            Vmath::Sadd(nzlines,m_lhom_z/2.0,z,1,z,1);

            for(m = 0; m < nzlines; ++m)
            {
				for(j = 0; j < nylines; ++j)
				{
					for(n = 0; n < npoints; ++n)
					{
						Vmath::Fill(1,x[n],tmp_xc = xc0 + n +(j*npoints) + (m*npoints*nylines), 1);
						Vmath::Fill(1,y[j],tmp_xc = xc1 + n +(j*npoints) + (m*npoints*nylines), 1);
						Vmath::Fill(1,z[m],tmp_xc = xc2 + n +(j*npoints) + (m*npoints*nylines), 1);
					}
				}
            }
        }


        /**
         * Write Tecplot Files Zone
         * @param   outfile    Output file name.
         * @param   expansion  Expansion that is considered
         */
        void ExpList3DHomogeneous2D::v_WriteTecplotZone(std::ostream &outfile, int expansion)
        {
            int i,j;

            int nquad0 = (*m_exp)[expansion]->GetNumPoints(0);
            int nquad1 = m_homogeneousBasis_y->GetNumPoints();
            int nquad2 = m_homogeneousBasis_z->GetNumPoints();

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


        void ExpList3DHomogeneous2D::v_WriteVtkPieceHeader(std::ostream &outfile, int expansion, int)
        {
            int i,j,k;
            int nquad0 = (*m_exp)[expansion]->GetNumPoints(0);
            int nquad1 = m_homogeneousBasis_y->GetNumPoints();
            int nquad2 = m_homogeneousBasis_z->GetNumPoints();
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
            outfile << "        <DataArray type=\"Float64\" "
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


        NekDouble ExpList3DHomogeneous2D::v_L2(
            const Array<OneD, const NekDouble> &inarray,
            const Array<OneD, const NekDouble> &soln)
        {
            int cnt = 0;
            NekDouble errL2,err = 0.0;
            Array<OneD, const NekDouble> w_y = m_homogeneousBasis_y->GetW();
            Array<OneD, const NekDouble> w_z = m_homogeneousBasis_z->GetW();

            int nylines = m_homogeneousBasis_y->GetNumPoints();
            int nzlines = m_homogeneousBasis_z->GetNumPoints();

            for(int m = 0; m < nzlines; ++m)
            {
                for(int n = 0; n < nylines; ++n)
                {
                    errL2 = m_lines[n+(m*nylines)]->L2(inarray + cnt, soln + cnt);
                    cnt  += m_lines[n+(m*nylines)]->GetTotPoints();
                    err  += errL2*errL2*w_y[n]*m_lhom_y*0.5*w_z[m]*m_lhom_z*0.5;
                }
            }

            return sqrt(err);
        }
    } //end of namespace
} //end of namespace
