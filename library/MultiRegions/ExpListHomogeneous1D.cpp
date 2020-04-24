///////////////////////////////////////////////////////////////////////////////
//
// File ExpListHomogeneous1D.cpp
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
// Description: An ExpList which is homogeneous in 1-direction
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpListHomogeneous1D.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdPointExp.h>
#include <LocalRegions/Expansion.h>
#include <LocalRegions/Expansion2D.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declaration for typedefs
        ExpListHomogeneous1D::ExpListHomogeneous1D():
            ExpList(),
            m_homogeneousBasis(LibUtilities::NullBasisSharedPtr),
            m_lhom(1),
            m_homogeneous1DBlockMat(MemoryManager<Homo1DBlockMatrixMap>::AllocateSharedPtr())
        {
        }

        ExpListHomogeneous1D::ExpListHomogeneous1D(
                   const LibUtilities::SessionReaderSharedPtr
                   &pSession,const LibUtilities::BasisKey &HomoBasis,
                   const NekDouble lhom,
                   const bool useFFT,
                   const bool dealiasing):
            ExpList(pSession),
            m_useFFT(useFFT),
            m_lhom(lhom),
            m_homogeneous1DBlockMat(MemoryManager<Homo1DBlockMatrixMap>::AllocateSharedPtr()),
            m_dealiasing(dealiasing)
        {
            ASSERTL2(HomoBasis != LibUtilities::NullBasisKey,
                     "Homogeneous Basis is a null basis");

            m_homogeneousBasis = LibUtilities::BasisManager()[HomoBasis];

            m_StripZcomm = m_session->DefinesSolverInfo("HomoStrip") ?
                           m_comm->GetColumnComm()->GetColumnComm()  :
                           m_comm->GetColumnComm();

            ASSERTL0( m_homogeneousBasis->GetNumPoints() %
                      m_StripZcomm->GetSize() == 0,
                      "HomModesZ should be a multiple of npz.");

            if (  (m_homogeneousBasis->GetBasisType() !=
                    LibUtilities::eFourierHalfModeRe)
               && (m_homogeneousBasis->GetBasisType() !=
                    LibUtilities::eFourierHalfModeIm) )
            {
                ASSERTL0(
                    (m_homogeneousBasis->GetNumPoints() /
                     m_StripZcomm->GetSize()) % 2 == 0,
                    "HomModesZ/npz should be an even integer.");
            }

            m_transposition = MemoryManager<LibUtilities::Transposition>
                                ::AllocateSharedPtr(HomoBasis, m_comm,
                                                    m_StripZcomm);

            m_planes = Array<OneD,ExpListSharedPtr>(
                                m_homogeneousBasis->GetNumPoints() /
                                m_StripZcomm->GetSize());

            if(m_useFFT)
            {
                m_FFT = LibUtilities::GetNektarFFTFactory().CreateInstance(
                                "NekFFTW", m_homogeneousBasis->GetNumPoints());
            }

            if(m_dealiasing)
            {
                if(m_useFFT)
                {
                    int size = m_homogeneousBasis->GetNumPoints() +
                               m_homogeneousBasis->GetNumPoints() / 2;
                    m_padsize = size + (size % 2);
                    m_FFT_deal = LibUtilities::GetNektarFFTFactory()
                                    .CreateInstance("NekFFTW", m_padsize);
                }
                else
                {
                    ASSERTL0(false, "Dealiasing available just in combination "
                                    "with FFTW");
                }
            }
        }


        /**
         * @param   In          ExpListHomogeneous1D object to copy.
         */
        ExpListHomogeneous1D::ExpListHomogeneous1D(const ExpListHomogeneous1D &In):
            ExpList(In,false),
            m_transposition(In.m_transposition),
            m_StripZcomm(In.m_StripZcomm),
            m_useFFT(In.m_useFFT),
            m_FFT(In.m_FFT),
            m_tmpIN(In.m_tmpIN),
            m_tmpOUT(In.m_tmpOUT),
            m_homogeneousBasis(In.m_homogeneousBasis),
            m_lhom(In.m_lhom),
            m_homogeneous1DBlockMat(In.m_homogeneous1DBlockMat),
            m_dealiasing(In.m_dealiasing),
            m_padsize(In.m_padsize)
        {
            m_planes = Array<OneD, ExpListSharedPtr>(In.m_planes.size());
        }

        ExpListHomogeneous1D::ExpListHomogeneous1D(const ExpListHomogeneous1D &In,
                                            const std::vector<unsigned int> &eIDs):
            ExpList(In,eIDs,false),
            m_transposition(In.m_transposition),
            m_useFFT(In.m_useFFT),
            m_FFT(In.m_FFT),
            m_tmpIN(In.m_tmpIN),
            m_tmpOUT(In.m_tmpOUT),
            m_homogeneousBasis(In.m_homogeneousBasis),
            m_lhom(In.m_lhom),
            m_homogeneous1DBlockMat(MemoryManager<Homo1DBlockMatrixMap>::AllocateSharedPtr()),
            m_dealiasing(In.m_dealiasing),
            m_padsize(In.m_padsize)
        {
            m_planes = Array<OneD, ExpListSharedPtr>(In.m_planes.size());
        }

        /**
         * Destructor
         */
        ExpListHomogeneous1D::~ExpListHomogeneous1D()
        {
        }
    
        void ExpListHomogeneous1D::v_HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                                         Array<OneD, NekDouble> &outarray, 
                                                         bool Shuff,
                                                         bool UnShuff)
        {
            // Forwards trans
            Homogeneous1DTrans(inarray,outarray,true,Shuff,UnShuff);
        }
    
        void ExpListHomogeneous1D::v_HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                                         Array<OneD, NekDouble> &outarray, 
                                                         bool Shuff,
                                                         bool UnShuff)
        {
            // Backwards trans
            Homogeneous1DTrans(inarray,outarray,false,Shuff,UnShuff);
        }

        /**
         * Dealiasing routine
         *
         * @param inarray1  First term of the product
         * @param inarray2  Second term of the product
         * @param outarray  Dealiased product
         */
        void ExpListHomogeneous1D::v_DealiasedProd(const Array<OneD, NekDouble> &inarray1,
                                                   const Array<OneD, NekDouble> &inarray2,
                                                   Array<OneD, NekDouble> &outarray)
        {
            int num_dofs = inarray1.size();
            int N = m_homogeneousBasis->GetNumPoints();

            Array<OneD, NekDouble> V1(num_dofs);
            Array<OneD, NekDouble> V2(num_dofs);
            Array<OneD, NekDouble> V1V2(num_dofs);

            if(m_WaveSpace)
            {
                V1 = inarray1;
                V2 = inarray2;
            }
            else
            {
                HomogeneousFwdTrans(inarray1,V1);
                HomogeneousFwdTrans(inarray2,V2);
            }

            int num_points_per_plane = num_dofs/m_planes.size();
            int num_proc;
            if(!m_session->DefinesSolverInfo("HomoStrip"))
            {
                num_proc             = m_comm->GetColumnComm()->GetSize();
            }
            else
            {
                num_proc             = m_StripZcomm->GetSize();
            }
            int num_dfts_per_proc    = num_points_per_plane / num_proc
                                        + (num_points_per_plane % num_proc > 0);

            Array<OneD, NekDouble> ShufV1(num_dfts_per_proc*N,0.0);
            Array<OneD, NekDouble> ShufV2(num_dfts_per_proc*N,0.0);
            Array<OneD, NekDouble> ShufV1V2(num_dfts_per_proc*N,0.0);

            Array<OneD, NekDouble> ShufV1_PAD_coef(m_padsize,0.0);
            Array<OneD, NekDouble> ShufV2_PAD_coef(m_padsize,0.0);
            Array<OneD, NekDouble> ShufV1_PAD_phys(m_padsize,0.0);
            Array<OneD, NekDouble> ShufV2_PAD_phys(m_padsize,0.0);

            Array<OneD, NekDouble> ShufV1V2_PAD_coef(m_padsize,0.0);
            Array<OneD, NekDouble> ShufV1V2_PAD_phys(m_padsize,0.0);

            m_transposition->Transpose(V1, ShufV1, false, LibUtilities::eXYtoZ);
            m_transposition->Transpose(V2, ShufV2, false, LibUtilities::eXYtoZ);

            // Looping on the pencils
            for(int i = 0 ; i < num_dfts_per_proc ; i++)
            {
                // Copying the i-th pencil pf lenght N into a bigger
                // pencil of lenght 2N We are in Fourier space
                Vmath::Vcopy(N, &(ShufV1[i*N]), 1, &(ShufV1_PAD_coef[0]), 1);
                Vmath::Vcopy(N, &(ShufV2[i*N]), 1, &(ShufV2_PAD_coef[0]), 1);

                // Moving to physical space using the padded system
                m_FFT_deal->FFTBwdTrans(ShufV1_PAD_coef, ShufV1_PAD_phys);
                m_FFT_deal->FFTBwdTrans(ShufV2_PAD_coef, ShufV2_PAD_phys);

                // Perfroming the vectors multiplication in physical space on
                // the padded system
                Vmath::Vmul(m_padsize, ShufV1_PAD_phys,   1,
                                       ShufV2_PAD_phys,   1,
                                       ShufV1V2_PAD_phys, 1);

                // Moving back the result (V1*V2)_phys in Fourier space, padded
                // system
                m_FFT_deal->FFTFwdTrans(ShufV1V2_PAD_phys, ShufV1V2_PAD_coef);

                // Copying the first half of the padded pencil in the full
                // vector (Fourier space)
                Vmath::Vcopy(N, &(ShufV1V2_PAD_coef[0]), 1,
                                &(ShufV1V2[i*N]),        1);
            }

            // Moving the results to the output
            if (m_WaveSpace)
            {
                m_transposition->Transpose(ShufV1V2, outarray, false,
                                       LibUtilities::eZtoXY);
            }
            else
            {
                m_transposition->Transpose(ShufV1V2, V1V2, false,
                                       LibUtilities::eZtoXY);
                HomogeneousBwdTrans(V1V2, outarray);
            }
        }

        /**
         * Dealiasing routine for dot product
         *
         * @param inarray1  First term of the product with dimension ndim
         *                  (e.g. U)
         * @param inarray2  Second term of the product with dimension ndim*nvec
         *                  (e.g. grad U)
         * @param outarray  Dealiased product with dimension nvec
         */
        void ExpListHomogeneous1D::v_DealiasedDotProd(
                        const Array<OneD, Array<OneD, NekDouble> > &inarray1,
                        const Array<OneD, Array<OneD, NekDouble> > &inarray2,
                        Array<OneD, Array<OneD, NekDouble> > &outarray)
        {
            int ndim = inarray1.size();
            ASSERTL1( inarray2.size() % ndim == 0,
                     "Wrong dimensions for DealiasedDotProd.");
            int nvec = inarray2.size() / ndim;

            int num_dofs = inarray1[0].size();
            int N = m_homogeneousBasis->GetNumPoints();

            int num_points_per_plane = num_dofs/m_planes.size();
            int num_proc;
            if(!m_session->DefinesSolverInfo("HomoStrip"))
            {
                num_proc             = m_comm->GetColumnComm()->GetSize();
            }
            else
            {
                num_proc             = m_StripZcomm->GetSize();
            }
            int num_dfts_per_proc    = num_points_per_plane / num_proc
                                        + (num_points_per_plane % num_proc > 0);

            // Get inputs in Fourier space
            Array<OneD, Array<OneD, NekDouble> > V1(ndim);
            Array<OneD, Array<OneD, NekDouble> > V2(ndim*nvec);
            if(m_WaveSpace)
            {
                for (int i = 0; i < ndim; i++)
                {
                    V1[i] = inarray1[i];
                }
                for (int i = 0; i < ndim*nvec; i++)
                {
                    V2[i] = inarray2[i];
                }
            }
            else
            {
                for (int i = 0; i < ndim; i++)
                {
                    V1[i] = Array<OneD, NekDouble> (num_dofs);
                    HomogeneousFwdTrans(inarray1[i],V1[i]);
                }
                for (int i = 0; i < ndim*nvec; i++)
                {
                    V2[i] = Array<OneD, NekDouble> (num_dofs);
                    HomogeneousFwdTrans(inarray2[i],V2[i]);
                }
            }

            // Allocate variables for ffts
            Array<OneD, Array<OneD, NekDouble> > ShufV1(ndim);
            Array<OneD, NekDouble>               ShufV1_PAD_coef(m_padsize,0.0);
            Array<OneD, Array<OneD, NekDouble> > ShufV1_PAD_phys(ndim);
            for (int i = 0; i < ndim; i++)
            {
                ShufV1[i]          = Array<OneD, NekDouble>
                                     (num_dfts_per_proc*N,0.0);
                ShufV1_PAD_phys[i] = Array<OneD, NekDouble>
                                     (m_padsize,0.0);
            }

            Array<OneD, Array<OneD, NekDouble> > ShufV2(ndim*nvec);
            Array<OneD, NekDouble>               ShufV2_PAD_coef(m_padsize,0.0);
            Array<OneD, Array<OneD, NekDouble> > ShufV2_PAD_phys(ndim*nvec);
            for (int i = 0; i < ndim*nvec; i++)
            {
                ShufV2[i]          = Array<OneD, NekDouble>
                                     (num_dfts_per_proc*N,0.0);
                ShufV2_PAD_phys[i] = Array<OneD, NekDouble>
                                     (m_padsize,0.0);
            }

            Array<OneD, Array<OneD, NekDouble> > ShufV1V2(nvec);
            Array<OneD, NekDouble>               ShufV1V2_PAD_coef(m_padsize,0.0);
            Array<OneD, NekDouble>               ShufV1V2_PAD_phys(m_padsize,0.0);
            for (int i = 0; i < nvec; i++)
            {
                ShufV1V2[i]          = Array<OneD, NekDouble>
                                     (num_dfts_per_proc*N,0.0);
            }

            // Transpose inputs
            for (int i = 0; i < ndim; i++)
            {
                m_transposition->Transpose(V1[i], ShufV1[i], false,
                                           LibUtilities::eXYtoZ);
            }
            for (int i = 0; i < ndim*nvec; i++)
            {
                m_transposition->Transpose(V2[i], ShufV2[i], false,
                                           LibUtilities::eXYtoZ);
            }

            // Looping on the pencils
            for(int i = 0 ; i < num_dfts_per_proc ; i++)
            {
                for (int j = 0; j < ndim; j++)
                {
                    // Copying the i-th pencil pf lenght N into a bigger
                    // pencil of lenght 1.5N We are in Fourier space
                    Vmath::Vcopy(N, &(ShufV1[j][i*N]), 1,
                                    &(ShufV1_PAD_coef[0]), 1);
                    // Moving to physical space using the padded system
                    m_FFT_deal->FFTBwdTrans(ShufV1_PAD_coef, ShufV1_PAD_phys[j]);
                }
                for (int j = 0; j < ndim*nvec; j++)
                {
                    Vmath::Vcopy(N, &(ShufV2[j][i*N]), 1,
                                    &(ShufV2_PAD_coef[0]), 1);
                    m_FFT_deal->FFTBwdTrans(ShufV2_PAD_coef, ShufV2_PAD_phys[j]);
                }

                // Performing the vectors multiplication in physical space on
                // the padded system
                for (int j = 0; j < nvec; j++)
                {
                    Vmath::Zero(m_padsize, ShufV1V2_PAD_phys, 1);
                    for (int k = 0; k < ndim; k++)
                    {
                        Vmath::Vvtvp(m_padsize, ShufV1_PAD_phys[k], 1,
                                                ShufV2_PAD_phys[j*ndim+k], 1,
                                                ShufV1V2_PAD_phys, 1,
                                                ShufV1V2_PAD_phys, 1);
                    }
                    // Moving back the result (V1*V2)_phys in Fourier space,
                    // padded system
                    m_FFT_deal->FFTFwdTrans(ShufV1V2_PAD_phys, ShufV1V2_PAD_coef);
                    // Copying the first half of the padded pencil in the full
                    // vector (Fourier space)
                    Vmath::Vcopy(N, &(ShufV1V2_PAD_coef[0]), 1,
                                    &(ShufV1V2[j][i*N]),     1);
                }
            }

            // Moving the results to the output
            if (m_WaveSpace)
            {
                for (int j = 0; j < nvec; j++)
                {
                    m_transposition->Transpose(ShufV1V2[j], outarray[j],
                                           false,
                                           LibUtilities::eZtoXY);
                }
            }
            else
            {
                Array<OneD, NekDouble> V1V2(num_dofs);
                for (int j = 0; j < nvec; j++)
                {
                    m_transposition->Transpose(ShufV1V2[j], V1V2, false,
                                       LibUtilities::eZtoXY);
                    HomogeneousBwdTrans(V1V2, outarray[j]);
                }
            }
        }

        /**
         * Forward transform
         */
        void ExpListHomogeneous1D::v_FwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;

            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->FwdTrans(inarray+cnt, tmparray = outarray + cnt1);
                cnt   += m_planes[n]->GetTotPoints();

                cnt1  += m_planes[n]->GetNcoeffs(); // need to skip ncoeffs
            }
            if(!m_WaveSpace)
            {
                HomogeneousFwdTrans(outarray,outarray);
            }
        }

        /**
         * Forward transform element by element
         */
        void ExpListHomogeneous1D::v_FwdTrans_IterPerExp(const Array<OneD,
                                                         const NekDouble> &inarray,
                                                         Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;

            //spectral element FwdTrans plane by plane
            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->FwdTrans_IterPerExp(inarray+cnt, tmparray = outarray + cnt1);
                cnt   += m_planes[n]->GetTotPoints();
                cnt1  += m_planes[n]->GetNcoeffs();
            }
            if(!m_WaveSpace)
            {
                HomogeneousFwdTrans(outarray,outarray);
            }
        }

        /**
         * Forward transform element by element with boundaries constrained
         */
        void ExpListHomogeneous1D::v_FwdTrans_BndConstrained(const Array<OneD,
                                                             const NekDouble> &inarray,
                                                             Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;

            //spectral element FwdTrans plane by plane
            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->FwdTrans_BndConstrained(inarray+cnt, tmparray = outarray + cnt1);
                cnt   += m_planes[n]->GetTotPoints();
                cnt1  += m_planes[n]->GetNcoeffs();
            }
            if(!m_WaveSpace)
            {
                HomogeneousFwdTrans(outarray,outarray);
            }
        }

        /**
         * Backward transform
         */
        void ExpListHomogeneous1D::v_BwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;

            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->BwdTrans(inarray+cnt, tmparray = outarray + cnt1);
                cnt  += m_planes[n]->GetNcoeffs();
                cnt1 += m_planes[n]->GetTotPoints();
            }
            if(!m_WaveSpace)
            {
                HomogeneousBwdTrans(outarray,outarray);
            }
        }

        /**
         * Backward transform element by element
         */
        void ExpListHomogeneous1D::v_BwdTrans_IterPerExp(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;

            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->BwdTrans_IterPerExp(inarray+cnt, tmparray = outarray + cnt1);

                cnt    += m_planes[n]->GetNcoeffs();
                cnt1   += m_planes[n]->GetTotPoints();
            }
            if(!m_WaveSpace)
            {
                HomogeneousBwdTrans(outarray,outarray);
            }
        }

        /**
         * Inner product
         */
        void ExpListHomogeneous1D::v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray, tmpIn;

            if(m_WaveSpace)
            {
                tmpIn = inarray;
            }
            else
            {
                tmpIn = Array<OneD, NekDouble> (inarray.size(), 0.0);
                HomogeneousFwdTrans(inarray,tmpIn);
            }

            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->IProductWRTBase(tmpIn+cnt, tmparray = outarray + cnt1);

                cnt1    += m_planes[n]->GetNcoeffs();
                cnt   += m_planes[n]->GetTotPoints();
            }
        }

        /**
         * Inner product element by element
         */
        void ExpListHomogeneous1D::v_IProductWRTBase_IterPerExp(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray, tmpIn;

            if(m_WaveSpace)
            {
                tmpIn = inarray;
            }
            else
            {
                tmpIn = Array<OneD, NekDouble> (inarray.size(), 0.0);
                HomogeneousFwdTrans(inarray,tmpIn);
            }

            for(int n = 0; n < m_planes.size(); ++n)
            {
                m_planes[n]->IProductWRTBase_IterPerExp(tmpIn+cnt, tmparray = outarray + cnt1);

                cnt1  += m_planes[n]->GetNcoeffs();
                cnt   += m_planes[n]->GetTotPoints();
            }
        }

        /**
         * Homogeneous transform Bwd/Fwd (MVM and FFT)
         */
        void ExpListHomogeneous1D::Homogeneous1DTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, 
                                                      bool IsForwards, 
                                                      bool Shuff,
                                                      bool UnShuff)
        {
            int num_dofs;

            if(IsForwards)
            {
                num_dofs = inarray.size();
            }
            else
            {
                num_dofs = outarray.size();
            }

            if(m_useFFT)
            {
                int num_points_per_plane = num_dofs/m_planes.size();
                int num_dfts_per_proc;
                if(!m_session->DefinesSolverInfo("HomoStrip"))
                {
                    int nP = m_comm->GetColumnComm()->GetSize();
                    num_dfts_per_proc = num_points_per_plane / nP
                                      + (num_points_per_plane % nP > 0);
                }
                else
                {
                    int nP = m_StripZcomm->GetSize();
                    num_dfts_per_proc = num_points_per_plane / nP
                                      + (num_points_per_plane % nP > 0);
                }

                Array<OneD, NekDouble> fft_in (num_dfts_per_proc*m_homogeneousBasis->GetNumPoints(),0.0);
                Array<OneD, NekDouble> fft_out(num_dfts_per_proc*m_homogeneousBasis->GetNumPoints(),0.0);

                if(Shuff)
                {
                    m_transposition->Transpose(inarray,fft_in,false,LibUtilities::eXYtoZ);
                }
                else
                {
                    Vmath::Vcopy(num_dfts_per_proc*m_homogeneousBasis->GetNumPoints(),
                                 inarray,1,fft_in,1);
                }

                if(IsForwards)
                {
                    for(int i = 0 ; i < num_dfts_per_proc ; i++)
                    {
                        m_FFT->FFTFwdTrans(m_tmpIN = fft_in + i*m_homogeneousBasis->GetNumPoints(), m_tmpOUT = fft_out + i*m_homogeneousBasis->GetNumPoints());
                    }
                }
                else
                {
                    for(int i = 0 ; i < num_dfts_per_proc ; i++)
                    {
                        m_FFT->FFTBwdTrans(m_tmpIN = fft_in + i*m_homogeneousBasis->GetNumPoints(), m_tmpOUT = fft_out + i*m_homogeneousBasis->GetNumPoints());
                    }
                }

                if(UnShuff)
                {
                    m_transposition->Transpose(fft_out,outarray,false,LibUtilities::eZtoXY);
                }
                else
                {
                    Vmath::Vcopy(num_dfts_per_proc*m_homogeneousBasis->GetNumPoints(),
                                 fft_out,1,outarray,1);
                }
            }
            else
            {
                DNekBlkMatSharedPtr blkmat;

                if(num_dofs == m_npoints) //transform phys space
                {
                    if(IsForwards)
                    {
                        blkmat = GetHomogeneous1DBlockMatrix(eForwardsPhysSpace1D);
                    }
                    else
                    {
                        blkmat = GetHomogeneous1DBlockMatrix(eBackwardsPhysSpace1D);
                    }
                }
                else
                {
                    if(IsForwards)
                    {
                        blkmat = GetHomogeneous1DBlockMatrix(eForwardsCoeffSpace1D);
                    }
                    else
                    {
                        blkmat = GetHomogeneous1DBlockMatrix(eBackwardsCoeffSpace1D);
                    }
                }

                int nrows = blkmat->GetRows();
                int ncols = blkmat->GetColumns();

                Array<OneD, NekDouble> sortedinarray(ncols,0.0);
                Array<OneD, NekDouble> sortedoutarray(nrows,0.0);

                if(Shuff)
                {
                    m_transposition->Transpose(inarray,sortedinarray,!IsForwards,LibUtilities::eXYtoZ);
                }
                else
                {
                    Vmath::Vcopy(ncols,inarray,1,sortedinarray,1);
                }

                // Create NekVectors from the given data arrays
                NekVector<NekDouble> in (ncols,sortedinarray,eWrapper);
                NekVector<NekDouble> out(nrows,sortedoutarray,eWrapper);

                // Perform matrix-vector multiply.
                out = (*blkmat)*in;

                if(UnShuff)
                {
                    m_transposition->Transpose(sortedoutarray,outarray,IsForwards,LibUtilities::eZtoXY);
                }
                else
                {
                    Vmath::Vcopy(nrows,sortedoutarray,1,outarray,1);
                }

            }
        }

        DNekBlkMatSharedPtr ExpListHomogeneous1D::GetHomogeneous1DBlockMatrix(Homogeneous1DMatType mattype) const
        {
            auto matrixIter = m_homogeneous1DBlockMat->find(mattype);

            if(matrixIter == m_homogeneous1DBlockMat->end())
            {
                return ((*m_homogeneous1DBlockMat)[mattype] =
                        GenHomogeneous1DBlockMatrix(mattype));
            }
            else
            {
                return matrixIter->second;
            }
        }


        DNekBlkMatSharedPtr ExpListHomogeneous1D::GenHomogeneous1DBlockMatrix(Homogeneous1DMatType mattype) const
        {
            DNekMatSharedPtr    loc_mat;
            DNekBlkMatSharedPtr BlkMatrix;
            int n_exp = 0;
            int num_trans_per_proc = 0;

            if((mattype == eForwardsCoeffSpace1D)
               ||(mattype == eBackwardsCoeffSpace1D)) // will operate on m_coeffs
            {
                n_exp = m_planes[0]->GetNcoeffs();
            }
            else
            {
                n_exp = m_planes[0]->GetTotPoints(); // will operatore on m_phys
            }

            num_trans_per_proc = n_exp/m_comm->GetColumnComm()->GetSize() + (n_exp%m_comm->GetColumnComm()->GetSize() > 0);

            Array<OneD,unsigned int> nrows(num_trans_per_proc);
            Array<OneD,unsigned int> ncols(num_trans_per_proc);

            if((mattype == eForwardsCoeffSpace1D)||(mattype == eForwardsPhysSpace1D))
            {
                nrows = Array<OneD, unsigned int>(num_trans_per_proc,m_homogeneousBasis->GetNumModes());
                ncols = Array<OneD, unsigned int>(num_trans_per_proc,m_homogeneousBasis->GetNumPoints());
            }
            else
            {
                nrows = Array<OneD, unsigned int>(num_trans_per_proc,m_homogeneousBasis->GetNumPoints());
                ncols = Array<OneD, unsigned int>(num_trans_per_proc,m_homogeneousBasis->GetNumModes());
            }

            MatrixStorage blkmatStorage = eDIAGONAL;
            BlkMatrix = MemoryManager<DNekBlkMat>
                ::AllocateSharedPtr(nrows,ncols,blkmatStorage);

            //Half Mode
            if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeRe || m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeIm)
            {
                StdRegions::StdPointExp StdPoint(m_homogeneousBasis->GetBasisKey());

                if((mattype == eForwardsCoeffSpace1D)||(mattype == eForwardsPhysSpace1D))
                {
                    StdRegions::StdMatrixKey matkey(StdRegions::eFwdTrans,
                                                    StdPoint.DetShapeType(),
                                                    StdPoint);

                    loc_mat = StdPoint.GetStdMatrix(matkey);
                }
                else
                {
                    StdRegions::StdMatrixKey matkey(StdRegions::eBwdTrans,
                                                    StdPoint.DetShapeType(),
                                                    StdPoint);

                    loc_mat = StdPoint.GetStdMatrix(matkey);
                }
            }
            //other cases
            else
            {
                StdRegions::StdSegExp StdSeg(m_homogeneousBasis->GetBasisKey());

                if((mattype == eForwardsCoeffSpace1D)||(mattype == eForwardsPhysSpace1D))
                {
                    StdRegions::StdMatrixKey matkey(StdRegions::eFwdTrans,
                                                    StdSeg.DetShapeType(),
                                                    StdSeg);

                    loc_mat = StdSeg.GetStdMatrix(matkey);
                }
                else
                {
                    StdRegions::StdMatrixKey matkey(StdRegions::eBwdTrans,
                                                    StdSeg.DetShapeType(),
                                                    StdSeg);

                    loc_mat = StdSeg.GetStdMatrix(matkey);
                }

            }

            // set up array of block matrices.
            for(int i = 0; i < num_trans_per_proc; ++i)
            {
                BlkMatrix->SetBlock(i,i,loc_mat);
            }

            return BlkMatrix;
        }

        std::vector<LibUtilities::FieldDefinitionsSharedPtr> ExpListHomogeneous1D::v_GetFieldDefinitions()
        {
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> returnval;

            // Set up Homogeneous length details.
            Array<OneD,LibUtilities::BasisSharedPtr> HomoBasis(1,m_homogeneousBasis);

            std::vector<NekDouble> HomoLen;
            HomoLen.push_back(m_lhom);

            std::vector<unsigned int> StripsIDs;

            bool strips;
            m_session->MatchSolverInfo("HomoStrip","True",strips,false);
            if (strips)
            {
                StripsIDs.push_back(m_transposition->GetStripID());
            }

            std::vector<unsigned int> PlanesIDs;
            int IDoffset = 0;

            // introduce a 2 plane offset for single mode case so can
            // be post-processed or used in MultiMode expansion.
            if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierSingleMode)
            {
                IDoffset  = 2;
            }

            for(int i = 0; i < m_planes.size(); i++)
            {
                PlanesIDs.push_back(m_transposition->GetPlaneID(i)+IDoffset);
            }

            m_planes[0]->GeneralGetFieldDefinitions(returnval, 1, HomoBasis,
                    HomoLen, strips, StripsIDs, PlanesIDs);
            return returnval;
        }

        void  ExpListHomogeneous1D::v_GetFieldDefinitions(std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef)
        {
            // Set up Homogeneous length details.
            Array<OneD,LibUtilities::BasisSharedPtr> HomoBasis(1,m_homogeneousBasis);

            std::vector<NekDouble> HomoLen;
            HomoLen.push_back(m_lhom);

            std::vector<unsigned int> StripsIDs;

            bool strips;
            m_session->MatchSolverInfo("HomoStrip","True",strips,false);
            if (strips)
            {
                StripsIDs.push_back(m_transposition->GetStripID());
            }

            std::vector<unsigned int> PlanesIDs;
            int IDoffset = 0;

            if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierSingleMode)
            {
                IDoffset = 2;
            }

            for(int i = 0; i < m_planes.size(); i++)
            {
                PlanesIDs.push_back(m_transposition->GetPlaneID(i)+IDoffset);
            }

            // enforce NumHomoDir == 1 by direct call
            m_planes[0]->GeneralGetFieldDefinitions(fielddef, 1, HomoBasis,
                    HomoLen, strips, StripsIDs, PlanesIDs);
        }


        /** This routine appends the data from the expansion list into
            the output format where each element is given by looping
            over its Fourier modes where as data in the expandion is
            stored with all consecutive elements and then the Fourier
            modes
         */
        void ExpListHomogeneous1D::v_AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, Array<OneD, NekDouble> &coeffs)
        {
            int i,n;
            int ncoeffs_per_plane = m_planes[0]->GetNcoeffs();

            // Determine mapping from element ids to location in
            // expansion list
            if (m_elmtToExpId.size() == 0)
            {
                for(i = 0; i < m_planes[0]->GetExpSize(); ++i)
                {
                    m_elmtToExpId[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
                }
            }

            for(i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                int eid     = m_elmtToExpId[fielddef->m_elementIDs[i]];
                int datalen = (*m_exp)[eid]->GetNcoeffs();

                for(n = 0; n < m_planes.size(); ++n)
                {
                    fielddata.insert(fielddata.end(),&coeffs[m_coeff_offset[eid]+n*ncoeffs_per_plane],&coeffs[m_coeff_offset[eid]+n*ncoeffs_per_plane]+datalen);
                }
            }
        }

        void ExpListHomogeneous1D::v_AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata)
        {
           v_AppendFieldData(fielddef,fielddata,m_coeffs);
        }

        //Extract the data in fielddata into the m_coeff list
        void ExpListHomogeneous1D::v_ExtractDataToCoeffs(
            LibUtilities::FieldDefinitionsSharedPtr    &fielddef,
            std::vector<NekDouble>       &fielddata,
            std::string                  &field,
            Array<OneD, NekDouble>       &coeffs)
        {
            int i,n;
            int offset  = 0;
            int nzmodes = 1;
            int datalen = fielddata.size()/fielddef->m_fields.size();
            std::vector<unsigned int> fieldDefHomoZids;


            // Find data location according to field definition
            for(i = 0; i < fielddef->m_fields.size(); ++i)
            {
                if(fielddef->m_fields[i] == field)
                {
                    break;
                }
                offset += datalen;
            }

            if(i == fielddef->m_fields.size())
            {
                cout << "Field "<< field<< "not found in data file. "  << endl;
            }
            else
            {

                int modes_offset = 0;
                int planes_offset = 0;
                Array<OneD, NekDouble> coeff_tmp;

                // Build map of plane IDs lying on this processor and determine
                // mapping from element ids to location in expansion list.
                if (m_zIdToPlane.size() == 0)
                {
                    for (i = 0; i < m_planes.size(); ++i)
                    {
                        m_zIdToPlane[m_transposition->GetPlaneID(i)] = i;
                    }

                    for (i = 0; i < m_planes[0]->GetExpSize(); ++i)
                    {
                        m_elmtToExpId[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
                    }
                }

                if(fielddef->m_numHomogeneousDir)
                {
                    nzmodes = fielddef->m_homogeneousZIDs.size();
                    fieldDefHomoZids = fielddef->m_homogeneousZIDs;
                }
                else // input file is 2D and so set nzmodes to 1
                {
                    nzmodes = 1;
                    fieldDefHomoZids.push_back(0);
                }

                // calculate number of modes in the current partition
                int ncoeffs_per_plane = m_planes[0]->GetNcoeffs();

                for(i = 0; i < fielddef->m_elementIDs.size(); ++i)
                {
                    if(fielddef->m_uniOrder == true) // reset modes_offset to zero
                    {
                        modes_offset = 0;
                    }

                    int datalen = LibUtilities::GetNumberOfCoefficients(fielddef->m_shapeType,
                                                                        fielddef->m_numModes,
                                                                        modes_offset);

                    auto it = m_elmtToExpId.find(fielddef->m_elementIDs[i]);

                    // ensure element is on this partition for parallel case.
                    if(it == m_elmtToExpId.end())
                    {
                        // increase offset for correct FieldData access
                        offset += datalen*nzmodes;
                        modes_offset += (*m_exp)[0]->GetNumBases() +
                                        fielddef->m_numHomogeneousDir;
                        continue;
                    }

                    int eid = it->second;
                    bool sameBasis = true;
                    for (int j = 0; j < fielddef->m_basis.size()-1; ++j)
                    {
                        if (fielddef->m_basis[j] != (*m_exp)[eid]->GetBasisType(j))
                        {
                            sameBasis = false;
                            break;
                        }
                    }

                    for(n = 0; n < nzmodes; ++n, offset += datalen)
                    {

                        it = m_zIdToPlane.find(fieldDefHomoZids[n]);

                        // Check to make sure this mode number lies in this field.
                        if (it == m_zIdToPlane.end())
                        {
                            continue;
                        }

                        planes_offset = it->second;
                        if(datalen == (*m_exp)[eid]->GetNcoeffs() && sameBasis)
                        {
                            Vmath::Vcopy(datalen,&fielddata[offset],1,&coeffs[m_coeff_offset[eid]+planes_offset*ncoeffs_per_plane],1);
                        }
                        else // unpack data to new order
                        {
                            (*m_exp)[eid]->ExtractDataToCoeffs(&fielddata[offset], fielddef->m_numModes,modes_offset,&coeffs[m_coeff_offset[eid] + planes_offset*ncoeffs_per_plane], fielddef->m_basis);
                        }
                    }
                    modes_offset += (*m_exp)[0]->GetNumBases() + fielddef->m_numHomogeneousDir;
                }
            }
        }

        //Extract the data in fielddata into the m_coeff list
        void ExpListHomogeneous1D::v_ExtractCoeffsToCoeffs(
                                                           const std::shared_ptr<ExpList> &fromExpList,const  Array<OneD, const NekDouble> &fromCoeffs, Array<OneD, NekDouble> &toCoeffs)
        {
            int i;
            int fromNcoeffs_per_plane = fromExpList->GetPlane(0)->GetNcoeffs();
            int toNcoeffs_per_plane = m_planes[0]->GetNcoeffs();
            Array<OneD, NekDouble> tocoeffs_tmp, fromcoeffs_tmp;

            for(i = 0; i < m_planes.size(); ++i)
            {
                m_planes[i]->ExtractCoeffsToCoeffs(fromExpList->GetPlane(i),fromcoeffs_tmp =  fromCoeffs + fromNcoeffs_per_plane*i, tocoeffs_tmp = toCoeffs + toNcoeffs_per_plane*i);
            }
        }

        void ExpListHomogeneous1D::v_WriteVtkPieceData(std::ostream &outfile, int expansion,
                                        std::string var)
        {
            // If there is only one plane (e.g. HalfMode), we write a 2D plane.
            if (m_planes.size() == 1)
            {
                m_planes[0]->WriteVtkPieceData(outfile, expansion, var);
                return;
            }

            int i;
            int nq = (*m_exp)[expansion]->GetTotPoints();
            int npoints_per_plane = m_planes[0]->GetTotPoints();

            // If we are using Fourier points, output extra plane to fill domain
            int outputExtraPlane = 0;
            Array<OneD, NekDouble> extraPlane;
            if ( m_homogeneousBasis->GetBasisType()   == LibUtilities::eFourier
               && m_homogeneousBasis->GetPointsType() ==
                    LibUtilities::eFourierEvenlySpaced)
            {
                outputExtraPlane = 1;
                // Get extra plane data
                if (m_StripZcomm->GetSize() == 1)
                {
                    extraPlane = m_phys + m_phys_offset[expansion];
                }
                else
                {
                    // Determine to and from rank for communication
                    int size     = m_StripZcomm->GetSize();
                    int rank     = m_StripZcomm->GetRank();
                    int fromRank = (rank+1) % size;
                    int toRank   = (rank == 0) ? size-1 : rank-1;
                    // Communicate using SendRecv
                    extraPlane = Array<OneD, NekDouble>(nq);
                    Array<OneD, NekDouble> send (nq,
                            m_phys + m_phys_offset[expansion]);
                    m_StripZcomm->SendRecv(toRank, send,
                                           fromRank, extraPlane);
                }
            }

            // printing the fields of that zone
            outfile << "        <DataArray type=\"Float64\" Name=\""
                    << var << "\">" << endl;
            outfile << "          ";
            for (int n = 0; n < m_planes.size(); ++n)
            {
                const Array<OneD, NekDouble> phys = m_phys + m_phys_offset[expansion] + n*npoints_per_plane;
                for(i = 0; i < nq; ++i)
                {
                    outfile << (fabs(phys[i]) < NekConstants::kNekZeroTol ? 0 : phys[i]) << " ";
                }
            }
            if (outputExtraPlane)
            {
                for(i = 0; i < nq; ++i)
                {
                    outfile << (fabs(extraPlane[i]) < NekConstants::kNekZeroTol ?
                                0 : extraPlane[i]) << " ";
                }
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
        }

        void ExpListHomogeneous1D::v_PhysInterp1DScaled(const NekDouble scale, const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt,cnt1;
            Array<OneD, NekDouble> tmparray;
            cnt  = m_planes[0]->GetTotPoints();
            cnt1 = m_planes[0]->Get1DScaledTotPoints(scale);

            ASSERTL1(m_planes.size()*cnt1 <= outarray.size(),"size of outarray does not match internal estimage");


            for(int i = 0; i < m_planes.size(); i++)
            {

                m_planes[i]->PhysInterp1DScaled(scale,inarray+i*cnt,
                                                 tmparray = outarray+i*cnt1);
            }
        }


        void ExpListHomogeneous1D::v_PhysGalerkinProjection1DScaled(const NekDouble scale, const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt,cnt1;
            Array<OneD, NekDouble> tmparray;
            cnt  = m_planes[0]->Get1DScaledTotPoints(scale);
            cnt1 = m_planes[0]->GetTotPoints();

            ASSERTL1(m_planes.size()*cnt <= inarray.size(),"size of outarray does not match internal estimage");


            for(int i = 0; i < m_planes.size(); i++)
            {
                m_planes[i]->PhysGalerkinProjection1DScaled(scale,inarray+i*cnt,
                                                 tmparray = outarray+i*cnt1);
            }

        }
        void ExpListHomogeneous1D::v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                               Array<OneD, NekDouble> &out_d0,
                                               Array<OneD, NekDouble> &out_d1,
                                               Array<OneD, NekDouble> &out_d2)
        {
            int nT_pts = inarray.size();          //number of total points = n. of Fourier points * n. of points per plane (nT_pts)
            int nP_pts = nT_pts/m_planes.size();    //number of points per plane = n of Fourier transform required (nP_pts)

            Array<OneD, NekDouble> temparray(nT_pts);
            Array<OneD, NekDouble> outarray(nT_pts);
            Array<OneD, NekDouble> tmp1;
            Array<OneD, NekDouble> tmp2;
            Array<OneD, NekDouble> tmp3;

            for(int i = 0; i < m_planes.size(); i++)
            {
                m_planes[i]->PhysDeriv(inarray + i*nP_pts ,tmp2 = out_d0 + i*nP_pts , tmp3 = out_d1 + i*nP_pts );
            }

            if(out_d2 != NullNekDouble1DArray)
            {
                if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourier || m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierSingleMode ||
                   m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeRe || m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeIm)
                {
                    if(m_WaveSpace)
                    {
                        temparray = inarray;
                    }
                    else
                    {
                        HomogeneousFwdTrans(inarray,temparray);
                    }

                    NekDouble sign = -1.0;
                    NekDouble beta;

                    //Half Mode
                    if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeRe)
                    {
                        beta = sign*2*M_PI*(m_transposition->GetK(0))/m_lhom;

                        Vmath::Smul(nP_pts,beta,temparray,1,outarray,1);
                    }
                    else if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeIm)
                    {
                        beta = -sign*2*M_PI*(m_transposition->GetK(0))/m_lhom;

                        Vmath::Smul(nP_pts,beta,temparray,1,outarray,1);
                    }

                    //Fully complex
                    else
                    {
                        for(int i = 0; i < m_planes.size(); i++)
                        {
                            beta = -sign*2*M_PI*(m_transposition->GetK(i))/m_lhom;

                            Vmath::Smul(nP_pts,beta,tmp1 = temparray + i*nP_pts,1,tmp2 = outarray + (i-int(sign))*nP_pts,1);

                            sign = -1.0*sign;
                        }
                    }

                    if(m_WaveSpace)
                    {
                        out_d2 = outarray;
                    }
                    else
                    {
                        HomogeneousBwdTrans(outarray,out_d2);
                    }
                }
                else
                {
                    if(!m_session->DefinesSolverInfo("HomoStrip"))
                    {
                        ASSERTL0(m_comm->GetColumnComm()->GetSize() == 1,
                                 "Parallelisation in the homogeneous direction "
                                 "implemented just for Fourier basis");
                    }
                    else
                    {
                        ASSERTL0(m_StripZcomm->GetSize()            == 1,
                                 "Parallelisation in the homogeneous direction "
                                 "implemented just for Fourier basis");
                    }

                    if(m_WaveSpace)
                    {
                        ASSERTL0(false, "Semi-phyisical time-stepping not "
                                        "implemented yet for non-Fourier "
                                        "basis");
                    }
                    else
                    {
                        StdRegions::StdSegExp StdSeg(m_homogeneousBasis->GetBasisKey());

                        m_transposition->Transpose(inarray,temparray,false,LibUtilities::eXYtoZ);

                        for(int i = 0; i < nP_pts; i++)
                        {
                            StdSeg.PhysDeriv(temparray + i*m_planes.size(), tmp2 = outarray + i*m_planes.size());
                        }

                        m_transposition->Transpose(outarray,out_d2,false,LibUtilities::eZtoXY);

                        Vmath::Smul(nT_pts,2.0/m_lhom,out_d2,1,out_d2,1);
                    }
                }
            }
        }

        void ExpListHomogeneous1D::v_PhysDeriv(Direction edir,
                                               const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &out_d)

        {
            int nT_pts = inarray.size();        //number of total points = n. of Fourier points * n. of points per plane (nT_pts)
            int nP_pts = nT_pts/m_planes.size();  //number of points per plane = n of Fourier transform required (nP_pts)

            int dir= (int)edir;

            Array<OneD, NekDouble> temparray(nT_pts);
            Array<OneD, NekDouble> outarray(nT_pts);
            Array<OneD, NekDouble> tmp1;
            Array<OneD, NekDouble> tmp2;

            if (dir < 2)
            {
                for(int i=0; i < m_planes.size(); i++)
                {
                    m_planes[i]->PhysDeriv(edir, inarray + i*nP_pts ,tmp2 = out_d + i*nP_pts);
                }
            }
            else
            {
                if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourier || m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierSingleMode ||
                   m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeRe || m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeIm)
                {
                    if(m_WaveSpace)
                    {
                        temparray = inarray;
                    }
                    else
                    {
                        HomogeneousFwdTrans(inarray,temparray);
                    }

                    NekDouble sign = -1.0;
                    NekDouble beta;

                    //HalfMode
                    if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeRe)
                    {
                        beta = 2*sign*M_PI*(m_transposition->GetK(0))/m_lhom;

                        Vmath::Smul(nP_pts,beta,temparray,1,outarray,1);
                    }
                    else if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierHalfModeIm)
                    {
                        beta = -2*sign*M_PI*(m_transposition->GetK(0))/m_lhom;

                        Vmath::Smul(nP_pts,beta,temparray,1,outarray,1);
                    }
                    //Fully complex
                    else
                    {
                        for(int i = 0; i < m_planes.size(); i++)
                        {
                            beta = -sign*2*M_PI*(m_transposition->GetK(i))/m_lhom;

                            Vmath::Smul(nP_pts,beta,tmp1 = temparray + i*nP_pts,1,tmp2 = outarray + (i-int(sign))*nP_pts,1);

                            sign = -1.0*sign;
                        }
                    }
                    if(m_WaveSpace)
                    {
                        out_d = outarray;
                    }
                    else
                    {
                        HomogeneousBwdTrans(outarray,out_d);
                    }
                }
                else
                {
                    if(!m_session->DefinesSolverInfo("HomoStrip"))
                    {
                        ASSERTL0(m_comm->GetColumnComm()->GetSize() == 1,
                                 "Parallelisation in the homogeneous direction "
                                 "implemented just for Fourier basis");
                    }
                    else
                    {
                        ASSERTL0(m_StripZcomm->GetSize()            == 1,
                                 "Parallelisation in the homogeneous direction "
                                 "implemented just for Fourier basis");
                    }

                    if(m_WaveSpace)
                    {
                        ASSERTL0(false,"Semi-phyisical time-stepping not implemented yet for non-Fourier basis");
                    }
                    else
                    {
                        StdRegions::StdSegExp StdSeg(m_homogeneousBasis->GetBasisKey());

                        m_transposition->Transpose(inarray,temparray,false,LibUtilities::eXYtoZ);

                        for(int i = 0; i < nP_pts; i++)
                        {
                            StdSeg.PhysDeriv(temparray + i*m_planes.size(), tmp2 = outarray + i*m_planes.size());
                        }

                        m_transposition->Transpose(outarray,out_d,false,LibUtilities::eZtoXY);

                        Vmath::Smul(nT_pts,2.0/m_lhom,out_d,1,out_d,1);
                    }
                }
            }
        }

        void ExpListHomogeneous1D::PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD, NekDouble> &out_d0,
                                             Array<OneD, NekDouble> &out_d1,
                                             Array<OneD, NekDouble> &out_d2)

        {
            v_PhysDeriv(inarray,out_d0,out_d1,out_d2);
        }

        void ExpListHomogeneous1D::PhysDeriv(Direction edir,
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD, NekDouble> &out_d)
        {
            v_PhysDeriv(edir,inarray,out_d);
        }

        LibUtilities::TranspositionSharedPtr ExpListHomogeneous1D::v_GetTransposition(void)
        {
            return m_transposition;
        }

        NekDouble ExpListHomogeneous1D::v_GetHomoLen(void)
        {
            return m_lhom;
        }

        void ExpListHomogeneous1D::v_SetHomoLen(const NekDouble lhom)
        {
            m_lhom = lhom;
        }

        Array<OneD, const unsigned int> ExpListHomogeneous1D::v_GetZIDs(void)
        {
            return m_transposition->GetPlanesIDs();
        }
    } //end of namespace
} //end of namespace
