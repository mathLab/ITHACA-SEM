///////////////////////////////////////////////////////////////////////////////
//
// File StdExpansion3D.cpp
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
        
        StdExpansion3D::StdExpansion3D(int                           numcoeffs, 
                                       const LibUtilities::BasisKey &Ba,
                                       const LibUtilities::BasisKey &Bb, 
                                       const LibUtilities::BasisKey &Bc) :
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
        void StdExpansion3D::PhysTensorDeriv(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &out_dx,
                  Array<OneD,       NekDouble> &out_dy,
                  Array<OneD,       NekDouble> &out_dz)
        {
            const int nquad0 = m_base[0]->GetNumPoints();
            const int nquad1 = m_base[1]->GetNumPoints();
            const int nquad2 = m_base[2]->GetNumPoints();

            Array<OneD, NekDouble> wsp(nquad0*nquad1*nquad2);
            
            // copy inarray to wsp in case inarray is used as outarray
            Vmath::Vcopy(nquad0*nquad1*nquad2, &inarray[0], 1, &wsp[0], 1);
            
            if (out_dx.num_elements() > 0)
            {
                NekDouble  *D0 = &((m_base[0]->GetD())->GetPtr())[0];

                Blas::Dgemm('N','N', nquad0,nquad1*nquad2,nquad0,1.0,
                            D0,nquad0,&wsp[0],nquad0,0.0,&out_dx[0],nquad0);
            }

            if (out_dy.num_elements() > 0) 
            {
                NekDouble   *D1 = &((m_base[1]->GetD())->GetPtr())[0];
                for (int j = 0; j < nquad2; ++j)
                {
                    Blas::Dgemm('N', 'T', nquad0, nquad1,      nquad1,
                                1.0, &wsp[j*nquad0*nquad1],    nquad0,
                                D1,                            nquad1,
                                0.0, &out_dy[j*nquad0*nquad1], nquad0);
                }
            }

            if (out_dz.num_elements() > 0) 
            {
                NekDouble     *D2 = &((m_base[2]->GetD())->GetPtr())[0];

                Blas::Dgemm('N','T',nquad0*nquad1,nquad2,nquad2,1.0,
                            &wsp[0],nquad0*nquad1,D2,nquad2,0.0,&out_dz[0],
                            nquad0*nquad1);
            }
        }

        NekDouble StdExpansion3D::v_PhysEvaluate(
            const Array<OneD, const NekDouble> &coords)
        {
            return PhysEvaluate(coords,m_phys);
        }
        
        NekDouble StdExpansion3D::v_PhysEvaluate(
            const Array<OneD, const NekDouble> &coords, 
            const Array<OneD, const NekDouble> &physvals)
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
            NekDouble *interpolatingNodes = 0;
            
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
            }
            
            // Interpolate in third coordinate direction
            I = m_base[2]->GetI(coords+2);
            interpolatingNodes = &I->GetPtr()[0];
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

        void StdExpansion3D::v_NegateFaceNormal(const int face)
        {
            m_negatedNormals[face] = true;
            for (int i = 0; i < GetCoordim(); ++i)
            {
                Vmath::Neg(m_faceNormals[face][i].num_elements(), 
                           m_faceNormals[face][i], 1);
            }
        }
    }//end namespace
}//end namespace
