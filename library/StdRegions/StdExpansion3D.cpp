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


        NekDouble StdExpansion3D::v_PhysEvaluate(const Array<OneD, const NekDouble> &coords )
        {
            return v_PhysEvaluate(coords,m_phys);
        }

    NekDouble StdExpansion3D::v_PhysEvaluate(const Array<OneD, const NekDouble> &coords, const Array<OneD, const NekDouble> & physvals)
    {
        NekDouble  value;
        ASSERTL2(coords[0] >= -1,"coord[0] < -1");
        ASSERTL2(coords[0] <=  1,"coord[0] >  1");
        ASSERTL2(coords[1] >= -1,"coord[1] < -1");
        ASSERTL2(coords[1] <=  1,"coord[1] >  1");
        ASSERTL2(coords[2] >= -1,"coord[2] < -1");
        ASSERTL2(coords[2] <=  1,"coord[2] >  1");

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
            sumFactorization_qr[i] =  Blas::Ddot(Qx, interpolatingNodes, 1, &physvals[ i*Qx ], 1);
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

        const NormalVector & StdExpansion3D::v_GetSurfaceNormal() const
        {
            return m_surfaceNormal;
        }
        const NormalVector & StdExpansion3D::v_GetFaceNormal(const int face) const
        {
            std::map<int, NormalVector>::const_iterator x;
            x = m_faceNormals.find(face);
            ASSERTL0 (x != m_faceNormals.end(),
                        "face normal not computed.");
            return x->second;
        }

    }//end namespace
}//end namespace
