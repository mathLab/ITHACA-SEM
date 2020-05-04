///////////////////////////////////////////////////////////////////////////////
//
// File ExpList1DHomogeneous2D.cpp
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
// Description: An ExpList1D which is homogeneous in 2 directions and so
// uses much of the functionality from a ExpList2D and its daughters
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/ExpList1DHomogeneous2D.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declaration for typedefs
        ExpList1DHomogeneous2D::ExpList1DHomogeneous2D():
            ExpListHomogeneous2D()
        {
        }

        // Constructor for ExpList1DHomogeneous2D to act as a Explist1D field
        ExpList1DHomogeneous2D::ExpList1DHomogeneous2D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                                       const LibUtilities::BasisKey &HomoBasis_y,
                                                       const LibUtilities::BasisKey &HomoBasis_z,
                                                       const NekDouble lhom_y,
                                                       const NekDouble lhom_z,
                                                       const bool useFFT,
                                                       const bool dealiasing,
                                                       const Array<OneD, ExpListSharedPtr> &points):
            ExpListHomogeneous2D(pSession,HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,useFFT,dealiasing)
        {
            int n;

            ASSERTL1(m_ny*m_nz == points.size(),
                    "Size of basis number of points and number of lines are "
                    "not the same");

            for(n = 0; n < points.size(); ++n)
            {
                m_lines[n] = points[n];
                (*m_exp).push_back(points[n]->GetExp(0));
            }

            SetCoeffPhys();
        }

        /**
         * @param   In          ExpList1DHomogeneous2D object to copy.
         */
        ExpList1DHomogeneous2D::ExpList1DHomogeneous2D(const ExpList1DHomogeneous2D &In):
            ExpListHomogeneous2D(In)
        {
            for(int n = 0; n < m_lines.size(); ++n)
            {
                m_lines[n] = In.m_lines[n];
            }

            SetCoeffPhys();
        }

        /**
         * Destructor
         */
        ExpList1DHomogeneous2D::~ExpList1DHomogeneous2D()
        {
        }

        void ExpList1DHomogeneous2D::SetCoeffPhys(void)
        {
            int i,n,cnt;
            int ncoeffs_per_line = m_lines[0]->GetNcoeffs();
            int npoints_per_line = m_lines[0]->GetTotPoints();

            int nyzlines = m_lines.size();

            // Set total coefficients and points
            m_ncoeffs = ncoeffs_per_line*nyzlines;
            m_npoints = npoints_per_line*nyzlines;

            m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
            m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};

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

        void ExpList1DHomogeneous2D::GetCoords(const int eid,
                                               Array<OneD, NekDouble> &xc0,
                                               Array<OneD, NekDouble> &xc1,
                                               Array<OneD, NekDouble> &xc2)
        {
            boost::ignore_unused(eid);

            int n,m,j;
            Array<OneD, NekDouble> tmp_xc;
            int nylines = m_homogeneousBasis_y->GetNumPoints();
            int nzlines = m_homogeneousBasis_z->GetNumPoints();
            int npoints = 1;

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

            m_lines[0]->GetCoords(x);


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
         *                          will be stored in this array.  This
         *                          coordinate might be  evaluated using the
         *                          predefined value \a m_lhom
         *
         * @param coord_2 After calculation, the \f$x_3\f$ coordinate
         *                          will be stored in this array. This
         *                          coordinate is evaluated using the
         *                          predefined value \a m_lhom
         */
        void ExpList1DHomogeneous2D::v_GetCoords(Array<OneD, NekDouble> &xc0,
                                                 Array<OneD, NekDouble> &xc1,
                                                 Array<OneD, NekDouble> &xc2)
        {
            int n,m,j;
            Array<OneD, NekDouble> tmp_xc;
            int npoints = 1;

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
         * Perform the 2D Forward transform of a set of points representing a plane of
         * boundary conditions which are merely the collection of the boundary conditions
         * coming from each 1D expansion.
         * @param   inarray    The value of the BC on each point of the y-z homogeneous plane.
         * @param   outarray   The value of the the coefficient of the 2D Fourier expansion
         */
        //void HomoFwdTrans2D(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        //{


        //}

        /**
         * Write Tecplot Files Zone
         * @param   outfile    Output file name.
         * @param   expansion  Expansion that is considered
         */
        void ExpList1DHomogeneous2D::v_WriteTecplotZone(std::ostream &outfile, int expansion)
        {
            int i,j;

            int nquad0 = 1;
            int nquad1 = m_homogeneousBasis_y->GetNumPoints();
            int nquad2 = m_homogeneousBasis_z->GetNumPoints();

            Array<OneD,NekDouble> coords[3];

            coords[0] = Array<OneD,NekDouble>(3*nquad0*nquad1*nquad2);
            coords[1] = coords[0] + nquad0*nquad1*nquad2;
            coords[2] = coords[1] + nquad0*nquad1*nquad2;

            GetCoords(expansion,coords[0],coords[1],coords[2]);

            outfile << "Zone, I=" << nquad1 << ", J=" << nquad0*nquad2
                    << ", F=Block" << std::endl;

            for(j = 0; j < nquad1; ++j)
            {
                for(i = 0; i < nquad2*GetCoordim(0)+1; ++i)
                {
                    outfile << coords[j][i] << " ";
                }
                outfile << std::endl;
            }
        }


        void ExpList1DHomogeneous2D::v_WriteVtkPieceHeader(std::ostream &outfile, int expansion, int istrip)
        {
            boost::ignore_unused(istrip);

            int i,j;

            int nquad0 = 1;
            int nquad1 = m_homogeneousBasis_y->GetNumPoints();
            int nquad2 = m_homogeneousBasis_z->GetNumPoints();

            int ntot = nquad0*nquad1*nquad2;
            int ntotminus = (nquad0)*(nquad1-1)*(nquad2-1);

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
            for (i = 0; i < nquad0; ++i)
            {
                for (j = 0; j < nquad1-1; ++j)
                {
                    outfile << j*nquad0 + i << " "
                            << j*nquad0 + i + 1 << " "
                            << (j+1)*nquad0 + i + 1 << " "
                            << (j+1)*nquad0 + i << endl;
                }
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "        <DataArray type=\"Int32\" "
                    << "Name=\"offsets\" format=\"ascii\">" << endl;
            for (i = 0; i < ntotminus; ++i)
            {
                outfile << i*4+4 << " ";
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "        <DataArray type=\"UInt8\" "
                    << "Name=\"types\" format=\"ascii\">" << endl;
            for (i = 0; i < ntotminus; ++i)
            {
                outfile << "9 ";
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "      </Cells>" << endl;
            outfile << "      <PointData>" << endl;
        }


    } //end of namespace
} //end of namespace
