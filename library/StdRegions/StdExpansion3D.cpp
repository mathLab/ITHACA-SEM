///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansioneD.cpp
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
// Description: Daughter of StdExpansion. This class contains routine
// which are common to 3D expansion. Typically this inolves physiocal
// space operations.
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdExpansion3D.h>

#ifdef max
#undef max
#endif

namespace Nektar

{
    namespace StdRegions
    {
    
    
    StdExpansion3D::StdExpansion3D()
    {
    }
    
    StdExpansion3D::StdExpansion3D(int numcoeffs, const LibUtilities::BasisKey &Ba, 
                       const LibUtilities::BasisKey &Bb, const LibUtilities::BasisKey &Bc):
        StdExpansion(numcoeffs,3,Ba,Bb,Bc)
    {
    }
    
    StdExpansion3D::StdExpansion3D(const StdExpansion3D &T):
        StdExpansion(T)
    {
    }
    
    StdExpansion3D::~StdExpansion3D()
    { 
    }
    
    namespace {    
        // Physical tenser terivative based on Spen's book page 151.
        void EasyDerivatives(const Array<OneD, const NekDouble> &inarray, 
                         Array<OneD, NekDouble> &outarray_dx, 
                         Array<OneD, NekDouble> &outarray_dy, 
                         Array<OneD, NekDouble> &outarray_dz, 
                         int Qx, int Qy, int Qz, DNekMatSharedPtr derivativeMatrix[3])
        {
            Array<OneD, NekDouble> u = Array<OneD, NekDouble>(Qx*Qy*Qz);
    
            // copy inarray to wsp in case inarray is used as outarray
            Vmath::Vcopy(Qx*Qy*Qz, &inarray[0], 1, &u[0], 1);
            
            double *Dx = &(derivativeMatrix[0]->GetPtr())[0];            
            double *Dy = &(derivativeMatrix[1]->GetPtr())[0];            
            double *Dz = &(derivativeMatrix[2]->GetPtr())[0];            
            
            int ijk = 0;
            for( int k = 0; k < Qz; ++k ) {         // Loop over xi_3
                for( int j = 0; j < Qy; ++j ) {     // Loop over xi_2
                    for( int i = 0; i < Qx; ++i ) { // Loop over xi_1
                        // calculate du/dxi_0
                        if(outarray_dx.num_elements() > 0 ) {
                            outarray_dx[ijk] = 0.0;
                            for( int p = 0; p < Qx; ++p ) {                            
                                int pjk = p + Qx*(j + Qy*k);
                                int ip =  i + Qx*p;
                                outarray_dx[ijk] += Dx[ip] * u[pjk];
                            }
                        }
                        
                        // calculate du/dxi_1
                        if(outarray_dy.num_elements() > 0 ) {
			                outarray_dy[ijk] = 0.0;
                            for( int q = 0; q < Qy; ++q ) {
                                int iqk = i + Qx*(q + Qy*k);
                                int jq  = j + Qy*q;
                                outarray_dy[ijk] += Dy[jq] * u[iqk];
                            }
                        }
                        
                        // calculate du/dxi_2
                        if(outarray_dz.num_elements() > 0 ) {
                            outarray_dz[ijk] = 0.0;
                            for( int r = 0; r < Qz; ++r ) {
                                int ijr = i + Qx*(j + Qy*r);
                                int kr  = k + Qz*r;
                                outarray_dz[ijk] += Dz[kr] * u[ijr];
                            }
                        }
                        ++ijk;
                    }
                }
            }
            
        }
        
                   
        // This is not working
        void AlternativeMethodForComputingTheDerivatives(const Array<OneD, const NekDouble> &inarray, 
                         Array<OneD, NekDouble> &outarray_dx, 
                         Array<OneD, NekDouble> &outarray_dy, 
                         Array<OneD, NekDouble> &outarray_dz, 
                         int Qx, int Qy, int Qz, DNekMatSharedPtr derivativeMatrix[3])
        {
            
            Array<OneD, NekDouble> wsp = Array<OneD, NekDouble>(Qx*Qy*Qz);
    
            // copy inarray to wsp in case inarray is used as outarray
            Vmath::Vcopy(Qx*Qy*Qz, &inarray[0], 1, &wsp[0], 1);
    
            DNekMatSharedPtr D0,D1,D2;
            D0 = derivativeMatrix[0];
            D1 = derivativeMatrix[1];
            D2 = derivativeMatrix[2];
            double *Dx = &(derivativeMatrix[0]->GetPtr())[0];            
            double *Dy = &(derivativeMatrix[1]->GetPtr())[0];            
            double *Dz = &(derivativeMatrix[2]->GetPtr())[0];            
   
    
            if(outarray_dx.num_elements() > 0)
            {
                Dx = &(D0->GetPtr())[0];
                for(int i=0; i < Qz; ++i)
                {
                    Blas::Dgemm('T', 'N', Qx, Qy, Qx, 1.0, Dx, Qx, &wsp[i*Qx*Qy], Qx, 0.0, &outarray_dx[i*Qx*Qy], Qx);
                }
            }
    
   
            if(outarray_dy.num_elements() > 0 ) {
                    Dy = &(D1->GetPtr())[0];                
                    for(int j=0; j<Qz; ++j){
                    
                    Blas:: Dgemm('T','N', Qx, Qy, Qy, 1.0, &wsp[j*Qx*Qy], Qx, Dy, Qy, 0.0, &outarray_dy[j*Qx*Qy], Qx);                   
                    //  Blas:: Dgemm('N','N', Qx, Qy, Qy, 1.0, &wsp[j*Qx*Qy], Qy, Dy, Qy, 0.0, &outarray_dy[j*Qx*Qy], Qy);    
                    }            
            }
     
            // calculate du/dx_2
            if(outarray_dz.num_elements() > 0) {   
                Dz = &(D2->GetPtr())[0];           
                for(int k=0; k < Qx*Qy; ++k)
                {
                        Blas:: Dgemv('N', Qz, Qz, 1.0, Dz, Qz, &wsp[0]+k, Qx*Qy, 0.0,  &outarray_dz[0]+k, Qx*Qy);
                }
            }   
        }        
    } // End of anonymous namespace

   
   
    void StdExpansion3D::PhysTensorDeriv(const Array<OneD, const NekDouble> &inarray, 
                         Array<OneD, NekDouble> &outarray_dx, 
                         Array<OneD, NekDouble> &outarray_dy, 
                         Array<OneD, NekDouble> &outarray_dz)
    {
        int    Qx = m_base[0]->GetNumPoints();
        int    Qy = m_base[1]->GetNumPoints();
        int    Qz = m_base[2]->GetNumPoints();
        DNekMatSharedPtr D[3];
        D[0] = m_base[0]->GetD();
        D[1] = m_base[1]->GetD();
        D[2] = m_base[2]->GetD();


    
        EasyDerivatives(inarray, outarray_dx, outarray_dy, outarray_dz, Qx, Qy, Qz, D); 
        
        // This is not working
        // AlternativeMethodForComputingTheDerivatives(inarray, outarray_dx, outarray_dy, outarray_dz, Qx, Qy, Qz, D); ;
    }
    
    
    NekDouble StdExpansion3D::v_PhysEvaluate(const Array<OneD, const NekDouble> &coords)
    {
        NekDouble  value;
        ASSERTL2(coords[0] < -1,"coord[0] < -1");
        ASSERTL2(coords[0] >  1,"coord[0] >  1");
        ASSERTL2(coords[1] < -1,"coord[1] < -1");
        ASSERTL2(coords[1] >  1,"coord[1] >  1");
        ASSERTL2(coords[2] < -1,"coord[2] < -1");
        ASSERTL2(coords[2] >  1,"coord[2] >  1");

        int Qx = m_base[0]->GetNumPoints();
        int Qy = m_base[1]->GetNumPoints();
        int Qz = m_base[2]->GetNumPoints();

        Array<OneD, NekDouble> sumFactorization_qr = Array<OneD, NekDouble>(Qy*Qz);
        Array<OneD, NekDouble> sumFactorization_r  = Array<OneD, NekDouble>(Qz);
        
        // Lagrangian interpolation matrix
        DNekMatSharedPtr I;
        double *interpolatingNodes = 0;

         // Interpolate first coordinate direction
        I = m_base[0]->GetI(coords);
        interpolatingNodes = &I->GetPtr()[0];
        for(int i = 0; i < Qy*Qz;++i)
        {
            sumFactorization_qr[i] =  Blas::Ddot(Qx, interpolatingNodes, 1, &m_phys[ i*Qx ], 1);
        }

        // Interpolate in second coordinate direction 
        I = m_base[1]->GetI(coords+1);
        interpolatingNodes = &I->GetPtr()[0];
        for(int j =0; j < Qz; ++j)
        {
            sumFactorization_r[j] = Blas::Ddot(Qy, interpolatingNodes, 1, &sumFactorization_qr[ j*Qy ], 1);
//            cout << "Index: " << j << endl;
//        for (int i = 0; i < Qz; ++i)
//        {
//            cout << interpolatingNodes[i] << ", " << sumFactorization_qr[i] << endl;
//        }
//        cout << endl;
        }

        // Interpolate in third coordinate direction 
        I = m_base[2]->GetI(coords+2);
        interpolatingNodes = &I->GetPtr()[0];
//        for (int i = 0; i < Qz; ++i)
//        {
//            cout << interpolatingNodes[i] << ", " << sumFactorization_r[i] << endl;
//        }
        value = Blas::Ddot(Qz, interpolatingNodes, 1, &sumFactorization_r[0], 1);            

        return value;
    }

        /**
         * Writes out the header for a <PIECE> VTK XML segment describing the
         * geometric information which comprises this element. This includes
         * vertex coordinates for each quadrature point, vertex connectivity
         * information, cell types and cell offset data.
         * 
         * @param   outfile     Output stream to write data to.
         */
        void StdExpansion3D::v_WriteVtkPieceHeader(std::ofstream &outfile)
        {
            int i,j,k;
            int coordim  = GetCoordim();
            int nquad0 = GetNumPoints(0);
            int nquad1 = GetNumPoints(1);
            int nquad2 = GetNumPoints(2);
            int ntot = nquad0*nquad1*nquad2;
            int ntotminus = (nquad0-1)*(nquad1-1)*(nquad2-1);
            
            Array<OneD,NekDouble> coords[3];
            coords[0] = Array<OneD,NekDouble>(ntot);
            coords[1] = Array<OneD,NekDouble>(ntot);
            coords[2] = Array<OneD,NekDouble>(ntot);
            GetCoords(coords[0],coords[1],coords[2]);

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
    }//end namespace
}//end namespace

/** 
 * $Log: StdExpansion3D.cpp,v $
 * Revision 1.19  2008/12/09 17:15:49  rcantao
 * Missing outarray_dy[ijk] = 0.0 prior to calculation.
 *
 * Revision 1.18  2008/11/23 00:32:57  sherwin
 * Removed dump value from PhysEvaluate
 *
 * Revision 1.17  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.16  2008/05/07 16:04:57  pvos
 * Mapping + Manager updates
 *
 * Revision 1.15  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.14  2007/11/08 14:27:07  ehan
 * Fixed PhysTensorDerivative3D matrix and improved L1 error up to 1e-15.
 *
 * Revision 1.13  2007/10/29 20:30:44  ehan
 * Fixed floating point approximation up to 1-e15 for PhysEvaluate.
 *
 * Revision 1.12  2007/10/15 20:37:59  ehan
 * Make changes of column major matrix
 *
 * Revision 1.11  2007/07/20 02:16:53  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.10  2007/05/22 02:01:49  bnelson
 * Changed Array::size to Array::num_elements.
 *
 * Fixed some compiler errors in assertions.
 *
 * Revision 1.9  2007/05/15 05:18:23  bnelson
 * Updated to use the new Array object.
 *
 * Revision 1.8  2007/04/10 14:00:45  sherwin
 * Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
 *
 * Revision 1.7  2007/04/04 20:48:17  sherwin
 * Update to handle SharedArrays
 *
 * Revision 1.6  2007/03/29 19:35:09  bnelson
 * Replaced boost::shared_array with SharedArray
 *
 * Revision 1.5  2007/03/20 16:58:43  sherwin
 * Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
 *
 * Revision 1.4  2007/01/18 18:44:45  bnelson
 * Updates to compile on Visual Studio 2005.
 *
 * Revision 1.3  2007/01/17 16:05:40  pvos
 * updated doxygen documentation
 *
 * Revision 1.2  2006/06/01 14:46:16  kirby
 * *** empty log message ***
 *
 * Revision 1.1  2006/05/04 18:58:31  kirby
 * *** empty log message ***
 *
 * Revision 1.9  2006/04/01 21:59:27  sherwin
 * Sorted new definition of ASSERT
 *
 * Revision 1.8  2006/03/21 09:21:32  sherwin
 * Introduced NekMemoryManager
 *
 * Revision 1.7  2006/02/27 23:47:23  sherwin
 *
 * Standard coding update upto compilation of StdHexExp.cpp
 *
 **/ 



