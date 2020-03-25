///////////////////////////////////////////////////////////////////////////////
//
// File ExpListHomogeneous2D.cpp
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
// Description: An ExpList which is homogeneous in 2-directions
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/ExpListHomogeneous2D.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdQuadExp.h>
#include <LocalRegions/Expansion.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declaration for typedefs
        ExpListHomogeneous2D::ExpListHomogeneous2D():
            ExpList(),
            m_homogeneousBasis_y(LibUtilities::NullBasisSharedPtr),
            m_homogeneousBasis_z(LibUtilities::NullBasisSharedPtr),
            m_lhom_y(1),
            m_lhom_z(1),
            m_homogeneous2DBlockMat(MemoryManager<Homo2DBlockMatrixMap>::AllocateSharedPtr())
        {
        }

        ExpListHomogeneous2D::ExpListHomogeneous2D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                                   const LibUtilities::BasisKey &HomoBasis_y,
                                                   const LibUtilities::BasisKey &HomoBasis_z,
                                                   const NekDouble lhom_y,
                                                   const NekDouble lhom_z,
                                                   const bool useFFT,
                                                   const bool dealiasing):
            ExpList(pSession),
            m_useFFT(useFFT),
            m_lhom_y(lhom_y),
            m_lhom_z(lhom_z),
            m_homogeneous2DBlockMat(MemoryManager<Homo2DBlockMatrixMap>::AllocateSharedPtr()),
            m_dealiasing(dealiasing)
        {
            ASSERTL2(HomoBasis_y != LibUtilities::NullBasisKey,
                     "Homogeneous Basis in y direction is a null basis");
            ASSERTL2(HomoBasis_z != LibUtilities::NullBasisKey,
                     "Homogeneous Basis in z direction is a null basis");

            m_homogeneousBasis_y = LibUtilities::BasisManager()[HomoBasis_y];
            m_homogeneousBasis_z = LibUtilities::BasisManager()[HomoBasis_z];

            m_transposition = MemoryManager<LibUtilities::Transposition>::AllocateSharedPtr(HomoBasis_y,HomoBasis_z,m_comm->GetColumnComm());

            m_Ycomm = m_comm->GetColumnComm()->GetRowComm();
            m_Zcomm = m_comm->GetColumnComm()->GetRowComm();

            m_ny = m_homogeneousBasis_y->GetNumPoints()/m_Ycomm->GetSize();
            m_nz = m_homogeneousBasis_z->GetNumPoints()/m_Zcomm->GetSize();

            m_lines = Array<OneD,ExpListSharedPtr>(m_ny*m_nz);

            if(m_useFFT)
            {
                m_FFT_y = LibUtilities::GetNektarFFTFactory().CreateInstance("NekFFTW", m_ny);
                m_FFT_z = LibUtilities::GetNektarFFTFactory().CreateInstance("NekFFTW", m_nz);
            }

            if(m_dealiasing)
            {
                ASSERTL0(m_comm->GetColumnComm()->GetSize() == 1,"Remove dealiasing if you want to run in parallel");
                SetPaddingBase();
            }
        }


        /**
         * @param   In          ExpListHomogeneous2D object to copy.
         */
        ExpListHomogeneous2D::ExpListHomogeneous2D(const ExpListHomogeneous2D &In):
            ExpList(In,false),
            m_useFFT(In.m_useFFT),
            m_FFT_y(In.m_FFT_y),
            m_FFT_z(In.m_FFT_z),
            m_transposition(In.m_transposition),
            m_Ycomm(In.m_Ycomm),
            m_Zcomm(In.m_Ycomm),
            m_homogeneousBasis_y(In.m_homogeneousBasis_y),
            m_homogeneousBasis_z(In.m_homogeneousBasis_z),
            m_lhom_y(In.m_lhom_y),
            m_lhom_z(In.m_lhom_z),
            m_homogeneous2DBlockMat(In.m_homogeneous2DBlockMat),
            m_ny(In.m_ny),
            m_nz(In.m_nz),
            m_dealiasing(In.m_dealiasing),
            m_padsize_y(In.m_padsize_y),
            m_padsize_z(In.m_padsize_z),
            MatFwdPAD(In.MatFwdPAD),
            MatBwdPAD(In.MatBwdPAD)
        {
            m_lines = Array<OneD, ExpListSharedPtr>(In.m_lines.size());
        }

        ExpListHomogeneous2D::ExpListHomogeneous2D(const ExpListHomogeneous2D &In,
                                            const std::vector<unsigned int> &eIDs):
            ExpList(In,eIDs,false),
            m_useFFT(In.m_useFFT),
            m_FFT_y(In.m_FFT_y),
            m_FFT_z(In.m_FFT_z),
            m_transposition(In.m_transposition),
            m_Ycomm(In.m_Ycomm),
            m_Zcomm(In.m_Ycomm),
            m_homogeneousBasis_y(In.m_homogeneousBasis_y),
            m_homogeneousBasis_z(In.m_homogeneousBasis_z),
            m_lhom_y(In.m_lhom_y),
            m_lhom_z(In.m_lhom_z),
            m_homogeneous2DBlockMat(MemoryManager<Homo2DBlockMatrixMap>::AllocateSharedPtr()),
            m_ny(In.m_ny),
            m_nz(In.m_nz),
            m_dealiasing(In.m_dealiasing),
            m_padsize_y(In.m_padsize_y),
            m_padsize_z(In.m_padsize_z),
            MatFwdPAD(In.MatFwdPAD),
            MatBwdPAD(In.MatBwdPAD)
        {
            m_lines = Array<OneD, ExpListSharedPtr>(In.m_lines.size());
        }

        /**
         * Destructor
         */
        ExpListHomogeneous2D::~ExpListHomogeneous2D()
        {
        }
        
        void ExpListHomogeneous2D::v_HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                                         Array<OneD, NekDouble> &outarray, 
                                                         bool Shuff,
                                                         bool UnShuff)
        {
            // Forwards trans
            Homogeneous2DTrans(inarray,outarray,true,Shuff,UnShuff);
        }
    
        void ExpListHomogeneous2D::v_HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                                         Array<OneD, NekDouble> &outarray, 
                                                         bool Shuff,
                                                         bool UnShuff)
        {
            // Backwards trans
            Homogeneous2DTrans(inarray,outarray,false,Shuff,UnShuff);
        }

        void ExpListHomogeneous2D::v_DealiasedProd(const Array<OneD, NekDouble> &inarray1,
                                                   const Array<OneD, NekDouble> &inarray2,
                                                   Array<OneD, NekDouble> &outarray)
        {
            int npoints = outarray.size(); // number of total physical points
            int nlines  = m_lines.size();  // number of lines == number of Fourier modes = number of Fourier coeff = number of points per slab
            int nslabs  = npoints/nlines;          // number of slabs = numebr of physical points per line

            Array<OneD, NekDouble> V1(npoints);
            Array<OneD, NekDouble> V2(npoints);
            Array<OneD, NekDouble> V1V2(npoints);
            Array<OneD, NekDouble> ShufV1(npoints);
            Array<OneD, NekDouble> ShufV2(npoints);
            Array<OneD, NekDouble> ShufV1V2(npoints);

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

            m_transposition->Transpose(V1,ShufV1,false,LibUtilities::eXtoYZ);
            m_transposition->Transpose(V2,ShufV2,false,LibUtilities::eXtoYZ);

            Array<OneD, NekDouble> PadV1_slab_coeff(m_padsize_y*m_padsize_z,0.0);
            Array<OneD, NekDouble> PadV2_slab_coeff(m_padsize_y*m_padsize_z,0.0);
            Array<OneD, NekDouble> PadRe_slab_coeff(m_padsize_y*m_padsize_z,0.0);

            Array<OneD, NekDouble> PadV1_slab_phys(m_padsize_y*m_padsize_z,0.0);
            Array<OneD, NekDouble> PadV2_slab_phys(m_padsize_y*m_padsize_z,0.0);
            Array<OneD, NekDouble> PadRe_slab_phys(m_padsize_y*m_padsize_z,0.0);

            NekVector<NekDouble> PadIN_V1(m_padsize_y*m_padsize_z,PadV1_slab_coeff,eWrapper);
            NekVector<NekDouble> PadOUT_V1(m_padsize_y*m_padsize_z,PadV1_slab_phys,eWrapper);

            NekVector<NekDouble> PadIN_V2(m_padsize_y*m_padsize_z,PadV2_slab_coeff,eWrapper);
            NekVector<NekDouble> PadOUT_V2(m_padsize_y*m_padsize_z,PadV2_slab_phys,eWrapper);

            NekVector<NekDouble> PadIN_Re(m_padsize_y*m_padsize_z,PadRe_slab_phys,eWrapper);
            NekVector<NekDouble> PadOUT_Re(m_padsize_y*m_padsize_z,PadRe_slab_coeff,eWrapper);

            //Looping on the slabs
            for(int j = 0 ; j< nslabs ; j++)
            {
                //Copying the j-th slab of size N*M into a bigger slab of lenght 2*N*M
                //We are in Fourier space
                for(int i = 0 ; i< m_nz ; i++)
                {
                    Vmath::Vcopy(m_ny,&(ShufV1[i*m_ny + j*nlines]),1,&(PadV1_slab_coeff[i*2*m_ny]),1);
                    Vmath::Vcopy(m_ny,&(ShufV2[i*m_ny + j*nlines]),1,&(PadV2_slab_coeff[i*2*m_ny]),1);
                }

                //Moving to physical space using the padded system
                PadOUT_V1 = (*MatBwdPAD)*PadIN_V1;
                PadOUT_V2 = (*MatBwdPAD)*PadIN_V2;

                //Perfroming the vectors multiplication in physical
                //space on the padded system
                Vmath::Vmul(m_padsize_y*m_padsize_z,PadV1_slab_phys,1,PadV2_slab_phys,1,PadRe_slab_phys,1);

                //Moving back the result (V1*V2)_phys in Fourier
                //space, padded system
                PadOUT_Re = (*MatFwdPAD)*PadIN_Re;

                //Copying the first half of the padded pencil in the
                //full vector (Fourier space)
                for (int i = 0; i < m_nz; i++)
                {
                    Vmath::Vcopy(m_ny,&(PadRe_slab_coeff[i*2*m_ny]),1,&(ShufV1V2[i*m_ny + j*nlines]),1);
                }
            }

            if(m_WaveSpace)
            {
                m_transposition->Transpose(ShufV1V2,outarray,false,LibUtilities::eYZtoX);
            }
            else
            {
                m_transposition->Transpose(ShufV1V2,V1V2,false,LibUtilities::eYZtoX);

                //Moving the results in physical space for the output
                HomogeneousBwdTrans(V1V2,outarray);
            }
        }

        void ExpListHomogeneous2D::v_DealiasedDotProd(
                        const Array<OneD, Array<OneD, NekDouble> > &inarray1,
                        const Array<OneD, Array<OneD, NekDouble> > &inarray2,
                        Array<OneD, Array<OneD, NekDouble> > &outarray)
        {
            // TODO Proper implementation of this
            int ndim = inarray1.size();
            ASSERTL1( inarray2.size() % ndim == 0,
                     "Wrong dimensions for DealiasedDotProd.");
            int nvec = inarray2.size() % ndim;
            int npts = inarray1[0].size();

            Array<OneD, NekDouble> out(npts);
            for (int i = 0; i < nvec; i++)
            {
                Vmath::Zero(npts, outarray[i], 1);
                for (int j = 0; j < ndim; j++)
                {
                    DealiasedProd(inarray1[j], inarray2[i*ndim+j], out);
                    Vmath::Vadd(npts, outarray[i], 1, out, 1, outarray[i], 1);
                }
            }
        }

        void ExpListHomogeneous2D::v_FwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            int nlines = m_lines.size();

            for(int n = 0; n < nlines; ++n)
            {
                m_lines[n]->FwdTrans(inarray+cnt, tmparray = outarray + cnt1);
                cnt   += m_lines[n]->GetTotPoints();
                cnt1  += m_lines[n]->GetNcoeffs();
            }
            if(!m_WaveSpace)
            {
                HomogeneousFwdTrans(outarray,outarray);
            }
        }

        void ExpListHomogeneous2D::v_FwdTrans_IterPerExp(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            int nlines = m_lines.size();

            for(int n = 0; n < nlines; ++n)
            {
                m_lines[n]->FwdTrans_IterPerExp(inarray+cnt, tmparray = outarray + cnt1);

                cnt   += m_lines[n]->GetTotPoints();
                cnt1  += m_lines[n]->GetNcoeffs();
            }
            if(!m_WaveSpace)
            {
                HomogeneousFwdTrans(outarray,outarray);
            }
        }
        
        void ExpListHomogeneous2D::v_BwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            int nlines = m_lines.size();

            for(int n = 0; n < nlines; ++n)
            {
                m_lines[n]->BwdTrans(inarray+cnt, tmparray = outarray + cnt1);
                cnt    += m_lines[n]->GetNcoeffs();
                cnt1   += m_lines[n]->GetTotPoints();
            }
            if(!m_WaveSpace)
            {
                HomogeneousBwdTrans(outarray,outarray);
            }
        }

        void ExpListHomogeneous2D::v_BwdTrans_IterPerExp(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            int nlines = m_lines.size();

            for(int n = 0; n < nlines; ++n)
            {
                m_lines[n]->BwdTrans_IterPerExp(inarray+cnt, tmparray = outarray + cnt1);

                cnt    += m_lines[n]->GetNcoeffs();
                cnt1   += m_lines[n]->GetTotPoints();
            }
            if(!m_WaveSpace)
            {
                HomogeneousBwdTrans(outarray,outarray);
            }
        }


        void ExpListHomogeneous2D::v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            int nlines = m_lines.size();

            for(int n = 0; n < nlines; ++n)
            {
                m_lines[n]->IProductWRTBase(inarray+cnt, tmparray = outarray + cnt1);

                cnt    += m_lines[n]->GetNcoeffs();
                cnt1   += m_lines[n]->GetTotPoints();
            }
        }

        void ExpListHomogeneous2D::v_IProductWRTBase_IterPerExp(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            int nlines = m_lines.size();

            for(int n = 0; n < nlines; ++n)
            {
                m_lines[n]->IProductWRTBase_IterPerExp(inarray+cnt, tmparray = outarray + cnt1);

                cnt    += m_lines[n]->GetNcoeffs();
                cnt1   += m_lines[n]->GetTotPoints();
            }
        }
        
        void ExpListHomogeneous2D::Homogeneous2DTrans(const Array<OneD, const NekDouble> &inarray, 
                                                      Array<OneD, NekDouble> &outarray, 
                                                      bool IsForwards, 
                                                      bool Shuff,
                                                      bool UnShuff)
        {
            boost::ignore_unused(Shuff, UnShuff);

            if(m_useFFT)
            {

                int n  = m_lines.size();   //number of Fourier points in the Fourier directions (x-z grid)
                int s  = inarray.size();   //number of total points = n. of Fourier points * n. of points per line
                int p  = s/n;                      //number of points per line = n of Fourier transform required

                Array<OneD, NekDouble> fft_in(s);
                Array<OneD, NekDouble> fft_out(s);

                m_transposition->Transpose(inarray,fft_in,false,LibUtilities::eXtoYZ);

                if(IsForwards)
                {
                    for(int i=0;i<(p*m_nz);i++)
                    {
                        m_FFT_y->FFTFwdTrans(m_tmpIN = fft_in + i*m_ny, m_tmpOUT = fft_out + i*m_ny);
                    }

                }
                else
                {
                    for(int i=0;i<(p*m_nz);i++)
                    {
                        m_FFT_y->FFTBwdTrans(m_tmpIN = fft_in + i*m_ny, m_tmpOUT = fft_out + i*m_ny);
                    }
                }

                m_transposition->Transpose(fft_out,fft_in,false,LibUtilities::eYZtoZY);

                if(IsForwards)
                {
                    for(int i=0;i<(p*m_ny);i++)
                    {
                        m_FFT_z->FFTFwdTrans(m_tmpIN = fft_in + i*m_nz, m_tmpOUT = fft_out + i*m_nz);
                    }

                }
                else
                {
                    for(int i=0;i<(p*m_ny);i++)
                    {
                        m_FFT_z->FFTBwdTrans(m_tmpIN = fft_in + i*m_nz, m_tmpOUT = fft_out + i*m_nz);
                    }
                }

                //TODO: required ZYtoX routine
                m_transposition->Transpose(fft_out,fft_in,false,LibUtilities::eZYtoYZ);

                m_transposition->Transpose(fft_in,outarray,false,LibUtilities::eYZtoX);

            }
            else
            {
                DNekBlkMatSharedPtr blkmatY;
                DNekBlkMatSharedPtr blkmatZ;

                if(inarray.size() == m_npoints) //transform phys space
                {
                    if(IsForwards)
                    {
                        blkmatY = GetHomogeneous2DBlockMatrix(eForwardsPhysSpaceY1D);
                        blkmatZ = GetHomogeneous2DBlockMatrix(eForwardsPhysSpaceZ1D);
                    }
                    else
                    {
                        blkmatY = GetHomogeneous2DBlockMatrix(eBackwardsPhysSpaceY1D);
                        blkmatZ = GetHomogeneous2DBlockMatrix(eBackwardsPhysSpaceZ1D);
                    }
                }
                else
                {
                    if(IsForwards)
                    {
                        blkmatY = GetHomogeneous2DBlockMatrix(eForwardsCoeffSpaceY1D);
                        blkmatZ = GetHomogeneous2DBlockMatrix(eForwardsCoeffSpaceZ1D);
                    }
                    else
                    {
                        blkmatY = GetHomogeneous2DBlockMatrix(eBackwardsCoeffSpaceY1D);
                        blkmatZ = GetHomogeneous2DBlockMatrix(eBackwardsCoeffSpaceZ1D);
                    }
                }

                int nrowsY = blkmatY->GetRows();
                int ncolsY = blkmatY->GetColumns();

                Array<OneD, NekDouble> sortedinarrayY(ncolsY);
                Array<OneD, NekDouble> sortedoutarrayY(nrowsY);

                int nrowsZ = blkmatZ->GetRows();
                int ncolsZ = blkmatZ->GetColumns();

                Array<OneD, NekDouble> sortedinarrayZ(ncolsZ);
                Array<OneD, NekDouble> sortedoutarrayZ(nrowsZ);

                NekVector<NekDouble> inY (ncolsY,sortedinarrayY,eWrapper);
                NekVector<NekDouble> outY(nrowsY,sortedoutarrayY,eWrapper);

                NekVector<NekDouble> inZ (ncolsZ,sortedinarrayZ,eWrapper);
                NekVector<NekDouble> outZ(nrowsZ,sortedoutarrayZ,eWrapper);

                m_transposition->Transpose(inarray,sortedinarrayY,!IsForwards,LibUtilities::eXtoYZ);

                outY = (*blkmatY)*inY;

                m_transposition->Transpose(sortedoutarrayY,sortedinarrayZ,false,LibUtilities::eYZtoZY);

                outZ = (*blkmatZ)*inZ;

                m_transposition->Transpose(sortedoutarrayZ,sortedoutarrayY,false,LibUtilities::eZYtoYZ);

                m_transposition->Transpose(sortedoutarrayY,outarray,false,LibUtilities::eYZtoX);
            }
        }
        
        DNekBlkMatSharedPtr ExpListHomogeneous2D::GetHomogeneous2DBlockMatrix(Homogeneous2DMatType mattype) const
        {
            auto matrixIter = m_homogeneous2DBlockMat->find(mattype);

            if(matrixIter == m_homogeneous2DBlockMat->end())
            {
                return ((*m_homogeneous2DBlockMat)[mattype] =
                        GenHomogeneous2DBlockMatrix(mattype));
            }
            else
            {
                return matrixIter->second;
            }
        }

        DNekBlkMatSharedPtr ExpListHomogeneous2D::GenHomogeneous2DBlockMatrix(Homogeneous2DMatType mattype) const
        {
            int i;
            int n_exp = 0;

            DNekMatSharedPtr    loc_mat;
            DNekBlkMatSharedPtr BlkMatrix;

            LibUtilities::BasisSharedPtr Basis;

            int NumPoints = 0;
            int NumModes = 0;
            int NumPencils = 0;

            if((mattype == eForwardsCoeffSpaceY1D) || (mattype == eBackwardsCoeffSpaceY1D)
               ||(mattype == eForwardsPhysSpaceY1D) || (mattype == eBackwardsPhysSpaceY1D))
            {
                Basis = m_homogeneousBasis_y;
                NumPoints  = m_homogeneousBasis_y->GetNumModes();
                NumModes   = m_homogeneousBasis_y->GetNumPoints();
                NumPencils = m_homogeneousBasis_z->GetNumPoints();
            }
            else
            {
                Basis = m_homogeneousBasis_z;
                NumPoints  = m_homogeneousBasis_z->GetNumModes();
                NumModes   = m_homogeneousBasis_z->GetNumPoints();
                NumPencils = m_homogeneousBasis_y->GetNumPoints();
            }

            if((mattype == eForwardsCoeffSpaceY1D) || (mattype == eForwardsCoeffSpaceZ1D)
               ||(mattype == eBackwardsCoeffSpaceY1D)||(mattype == eBackwardsCoeffSpaceZ1D))
            {
                n_exp = m_lines[0]->GetNcoeffs();
            }
            else
            {
                n_exp = m_lines[0]->GetTotPoints(); // will operatore on m_phys
            }

            Array<OneD,unsigned int> nrows(n_exp);
            Array<OneD,unsigned int> ncols(n_exp);

            if((mattype == eForwardsCoeffSpaceY1D)||(mattype == eForwardsPhysSpaceY1D) ||
               (mattype == eForwardsCoeffSpaceZ1D)||(mattype == eForwardsPhysSpaceZ1D))
            {
                nrows = Array<OneD, unsigned int>(n_exp*NumPencils,NumModes);
                ncols = Array<OneD, unsigned int>(n_exp*NumPencils,NumPoints);
            }
            else
            {
                nrows = Array<OneD, unsigned int>(n_exp*NumPencils,NumPoints);
                ncols = Array<OneD, unsigned int>(n_exp*NumPencils,NumModes);
            }

            MatrixStorage blkmatStorage = eDIAGONAL;
            BlkMatrix = MemoryManager<DNekBlkMat>::AllocateSharedPtr(nrows,ncols,blkmatStorage);

            StdRegions::StdSegExp StdSeg(Basis->GetBasisKey());

            if((mattype == eForwardsCoeffSpaceY1D)||(mattype == eForwardsPhysSpaceY1D) ||
               (mattype == eForwardsCoeffSpaceZ1D)||(mattype == eForwardsPhysSpaceZ1D))
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

            // set up array of block matrices.
            for(i = 0; i < (n_exp*NumPencils); ++i)
            {
                BlkMatrix->SetBlock(i,i,loc_mat);
            }

            return BlkMatrix;
        }

        std::vector<LibUtilities::FieldDefinitionsSharedPtr> ExpListHomogeneous2D::v_GetFieldDefinitions()
        {
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> returnval;
            // Set up Homogeneous length details.
            Array<OneD,LibUtilities::BasisSharedPtr> HomoBasis(2);
            HomoBasis[0] = m_homogeneousBasis_y;
            HomoBasis[1] = m_homogeneousBasis_z;

            std::vector<NekDouble> HomoLen(2);
            HomoLen[0] = m_lhom_y;
            HomoLen[1] = m_lhom_z;

            int nhom_modes_y = m_homogeneousBasis_y->GetNumModes();
            int nhom_modes_z = m_homogeneousBasis_z->GetNumModes();

            std::vector<unsigned int> sIDs
                = LibUtilities::NullUnsignedIntVector;

            std::vector<unsigned int> yIDs;
            std::vector<unsigned int> zIDs;

            for(int n = 0; n < nhom_modes_z; ++n)
            {
                for(int m = 0; m < nhom_modes_y; ++m)
                {
                    zIDs.push_back(n);
                    yIDs.push_back(m);
                }
            }

            m_lines[0]->GeneralGetFieldDefinitions(returnval, 2, HomoBasis,
                                                     HomoLen, false,
                                                     sIDs, zIDs, yIDs);
            return returnval;
        }

        void  ExpListHomogeneous2D::v_GetFieldDefinitions(std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef)
        {
            // Set up Homogeneous length details.
            Array<OneD,LibUtilities::BasisSharedPtr> HomoBasis(2);
            HomoBasis[0] = m_homogeneousBasis_y;
            HomoBasis[1] = m_homogeneousBasis_z;
            std::vector<NekDouble> HomoLen(2);
            HomoLen[0] = m_lhom_y;
            HomoLen[1] = m_lhom_z;

            int nhom_modes_y = m_homogeneousBasis_y->GetNumModes();
            int nhom_modes_z = m_homogeneousBasis_z->GetNumModes();

            std::vector<unsigned int> sIDs
                =LibUtilities::NullUnsignedIntVector;

            std::vector<unsigned int> yIDs;
            std::vector<unsigned int> zIDs;

            for(int n = 0; n < nhom_modes_z; ++n)
            {
                for(int m = 0; m < nhom_modes_y; ++m)
                {
                    zIDs.push_back(n);
                    yIDs.push_back(m);
                }
            }

            // enforce NumHomoDir == 1 by direct call
             m_lines[0]->GeneralGetFieldDefinitions(fielddef, 2, HomoBasis,
                                                    HomoLen, false,
                                                    sIDs, zIDs, yIDs);
        }

        void ExpListHomogeneous2D::v_AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, Array<OneD, NekDouble> &coeffs)
        {
            int i,k;

            int NumMod_y = m_homogeneousBasis_y->GetNumModes();
            int NumMod_z = m_homogeneousBasis_z->GetNumModes();

            int ncoeffs_per_line = m_lines[0]->GetNcoeffs();

            // Determine mapping from element ids to location in
            // expansion list
            map<int, int> ElmtID_to_ExpID;
            for(i = 0; i < m_lines[0]->GetExpSize(); ++i)
            {
                ElmtID_to_ExpID[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
            }

            for(i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                int eid     = ElmtID_to_ExpID[fielddef->m_elementIDs[i]];
                int datalen = (*m_exp)[eid]->GetNcoeffs();

                for(k = 0; k < (NumMod_y*NumMod_z); ++k)
                {
                    fielddata.insert(fielddata.end(),&coeffs[m_coeff_offset[eid]+k*ncoeffs_per_line],&coeffs[m_coeff_offset[eid]+k*ncoeffs_per_line]+datalen);
                }
            }
        }

        void ExpListHomogeneous2D::v_AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata)
        {
            v_AppendFieldData(fielddef,fielddata,m_coeffs);
        }

        //Extract the data in fielddata into the m_coeff list
        void ExpListHomogeneous2D::v_ExtractDataToCoeffs(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, std::string &field, Array<OneD, NekDouble> &coeffs)
        {
            int i,k;
            int offset = 0;
            int datalen = fielddata.size()/fielddef->m_fields.size();
            int ncoeffs_per_line = m_lines[0]->GetNcoeffs();
            int NumMod_y = m_homogeneousBasis_y->GetNumModes();
            int NumMod_z = m_homogeneousBasis_z->GetNumModes();

            // Find data location according to field definition
            for(i = 0; i < fielddef->m_fields.size(); ++i)
            {
                if(fielddef->m_fields[i] == field)
                {
                    break;
                }
                offset += datalen;
            }

            ASSERTL0(i!= fielddef->m_fields.size(),"Field not found in data file");

            // Determine mapping from element ids to location in
            // expansion list
            map<int, int> ElmtID_to_ExpID;
            for(i = 0; i < m_lines[0]->GetExpSize(); ++i)
            {
                ElmtID_to_ExpID[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
            }

            for(i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                int eid = ElmtID_to_ExpID[fielddef->m_elementIDs[i]];
                int datalen = (*m_exp)[eid]->GetNcoeffs();

                for(k = 0; k < (NumMod_y*NumMod_z); ++k)
                {
                    Vmath::Vcopy(datalen,&fielddata[offset],1,&coeffs[m_coeff_offset[eid] + k*ncoeffs_per_line],1);
                    offset += datalen;
                }
            }
        }

        void ExpListHomogeneous2D::v_WriteVtkPieceData(std::ostream &outfile, int expansion,
                                        std::string var)
        {
            int i;
            int nq = (*m_exp)[expansion]->GetTotPoints();
            int npoints_per_line = m_lines[0]->GetTotPoints();

            // printing the fields of that zone
            outfile << "        <DataArray type=\"Float64\" Name=\""
                    << var << "\">" << endl;
            outfile << "          ";
            for (int n = 0; n < m_lines.size(); ++n)
            {
                const Array<OneD, NekDouble> phys = m_phys + m_phys_offset[expansion] + n*npoints_per_line;
                for(i = 0; i < nq; ++i)
                {
                    outfile << (fabs(phys[i]) < NekConstants::kNekZeroTol ? 0 : phys[i]) << " ";
                }
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
        }

        void ExpListHomogeneous2D::v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                               Array<OneD, NekDouble> &out_d0,
                                               Array<OneD, NekDouble> &out_d1,
                                               Array<OneD, NekDouble> &out_d2)

        {
            int nyzlines      = m_lines.size();   //number of Fourier points in the Fourier directions (nF_pts)
            int npoints       = inarray.size();   //number of total points = n. of Fourier points * n. of points per line (nT_pts)
            int n_points_line = npoints/nyzlines;         //number of points per line

            Array<OneD, NekDouble> temparray(npoints);
            Array<OneD, NekDouble> temparray1(npoints);
            Array<OneD, NekDouble> temparray2(npoints);
            Array<OneD, NekDouble> tmp1;
            Array<OneD, NekDouble> tmp2;
            Array<OneD, NekDouble> tmp3;

            for( int i=0 ; i<nyzlines ; i++ )
            {
                m_lines[i]->PhysDeriv( tmp1 = inarray + i*n_points_line ,tmp2 = out_d0 + i*n_points_line);
            }

            if(m_homogeneousBasis_y->GetBasisType() == LibUtilities::eFourier && m_homogeneousBasis_z->GetBasisType() == LibUtilities::eFourier)
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

                //along y
                for(int i = 0; i < m_ny; i++)
                {
                    beta = -sign*2*M_PI*(i/2)/m_lhom_y;

                    for(int j = 0; j < m_nz; j++)
                    {
                        Vmath::Smul(n_points_line,beta,tmp1 = temparray + n_points_line*(i+j*m_ny),1, tmp2 = temparray1 + n_points_line*((i-int(sign))+j*m_ny),1);
                    }

                    sign = -1.0*sign;
                }

                //along z
                sign = -1.0;
                for(int i = 0; i < m_nz; i++)
                {
                    beta = -sign*2*M_PI*(i/2)/m_lhom_z;
                    Vmath::Smul(m_ny*n_points_line,beta,tmp1 = temparray + i*m_ny*n_points_line,1,tmp2 = temparray2 + (i-int(sign))*m_ny*n_points_line,1);
                    sign = -1.0*sign;
                }
                if(m_WaveSpace)
                {
                    out_d1 = temparray1;
                    out_d2 = temparray2;
                }
                else
                {
                    HomogeneousBwdTrans(temparray1,out_d1);
                    HomogeneousBwdTrans(temparray2,out_d2);
                }
            }
            else
            {
                if(m_WaveSpace)
                {
                    ASSERTL0(false,"Semi-phyisical time-stepping not implemented yet for non-Fourier basis")
                        }
                else
                {
                    StdRegions::StdQuadExp StdQuad(m_homogeneousBasis_y->GetBasisKey(),m_homogeneousBasis_z->GetBasisKey());

                    m_transposition->Transpose(inarray,temparray,false,LibUtilities::eXtoYZ);

                    for(int i = 0; i < n_points_line; i++)
                    {
                        StdQuad.PhysDeriv(tmp1 = temparray + i*nyzlines, tmp2 = temparray1 + i*nyzlines, tmp3 = temparray2 + i*nyzlines);
                    }

                    m_transposition->Transpose(temparray1,out_d1,false,LibUtilities::eYZtoX);
                    m_transposition->Transpose(temparray2,out_d2,false,LibUtilities::eYZtoX);
                    Vmath::Smul(npoints,2.0/m_lhom_y,out_d1,1,out_d1,1);
                    Vmath::Smul(npoints,2.0/m_lhom_z,out_d2,1,out_d2,1);
                }
            }
        }

        void ExpListHomogeneous2D::v_PhysDeriv(Direction edir,
                                               const Array<OneD, const NekDouble> &inarray,
                                               Array<OneD, NekDouble> &out_d)

        {
            int nyzlines      = m_lines.size();   //number of Fourier points in the Fourier directions (nF_pts)
            int npoints       = inarray.size();   //number of total points = n. of Fourier points * n. of points per line (nT_pts)
            int n_points_line = npoints/nyzlines;         //number of points per line
            //convert enum into int
            int dir = (int)edir;

            Array<OneD, NekDouble> temparray(npoints);
            Array<OneD, NekDouble> temparray1(npoints);
            Array<OneD, NekDouble> temparray2(npoints);
            Array<OneD, NekDouble> tmp1;
            Array<OneD, NekDouble> tmp2;
            Array<OneD, NekDouble> tmp3;

            if (dir < 1)
            {
                for( int i=0 ; i<nyzlines ; i++)
                {
                    m_lines[i]->PhysDeriv( tmp1 = inarray + i*n_points_line ,tmp2 = out_d + i*n_points_line);
                }
            }
            else
            {
                if(m_homogeneousBasis_y->GetBasisType() == LibUtilities::eFourier && m_homogeneousBasis_z->GetBasisType() == LibUtilities::eFourier)
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

                    if (dir == 1)
                    {
                        //along y
                        for(int i = 0; i < m_ny; i++)
                        {
                            beta = -sign*2*M_PI*(i/2)/m_lhom_y;

                            for(int j = 0; j < m_nz; j++)
                            {
                                Vmath::Smul(n_points_line,beta,tmp1 = temparray + n_points_line*(i+j*m_ny),1, tmp2 = temparray1 + n_points_line*((i-int(sign))+j*m_ny),1);
                            }
                            sign = -1.0*sign;
                        }
                        if(m_WaveSpace)
                        {
                            out_d = temparray1;
                        }
                        else
                        {
                            HomogeneousBwdTrans(temparray1,out_d);
                        }
                    }
                    else
                    {
                        //along z
                        for(int i = 0; i < m_nz; i++)
                        {
                            beta = -sign*2*M_PI*(i/2)/m_lhom_z;
                            Vmath::Smul(m_ny*n_points_line,beta,tmp1 = temparray + i*m_ny*n_points_line,1,tmp2 = temparray2 + (i-int(sign))*m_ny*n_points_line,1);
                            sign = -1.0*sign;
                        }
                        if(m_WaveSpace)
                        {
                                                    out_d = temparray2;
                        }
                        else
                        {
                            HomogeneousBwdTrans(temparray2,out_d);
                        }
                    }
                }
                else
                {
                    if(m_WaveSpace)
                    {
                        ASSERTL0(false,"Semi-phyisical time-stepping not implemented yet for non-Fourier basis")
                            }
                    else
                    {
                        StdRegions::StdQuadExp StdQuad(m_homogeneousBasis_y->GetBasisKey(),m_homogeneousBasis_z->GetBasisKey());

                        m_transposition->Transpose(inarray,temparray,false,LibUtilities::eXtoYZ);

                        for(int i = 0; i < n_points_line; i++)
                        {
                            StdQuad.PhysDeriv(tmp1 = temparray + i*nyzlines, tmp2 = temparray1 + i*nyzlines, tmp3 = temparray2 + i*nyzlines);
                        }

                        if (dir == 1)
                        {
                            m_transposition->Transpose(temparray1,out_d,false,LibUtilities::eYZtoX);
                            Vmath::Smul(npoints,2.0/m_lhom_y,out_d,1,out_d,1);
                        }
                        else
                        {
                            m_transposition->Transpose(temparray2,out_d,false,LibUtilities::eYZtoX);
                            Vmath::Smul(npoints,2.0/m_lhom_z,out_d,1,out_d,1);
                        }
                    }
                }
            }
        }

        void ExpListHomogeneous2D::PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD, NekDouble> &out_d0,
                                             Array<OneD, NekDouble> &out_d1,
                                             Array<OneD, NekDouble> &out_d2)

        {
            v_PhysDeriv(inarray,out_d0,out_d1,out_d2);
        }

        void ExpListHomogeneous2D::PhysDeriv(Direction edir,
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD, NekDouble> &out_d)
        {
            //convert int into enum
            v_PhysDeriv(edir,inarray,out_d);
        }

        void ExpListHomogeneous2D::SetPaddingBase(void)
        {
            NekDouble size_y = 1.5*m_ny;
            NekDouble size_z = 1.5*m_nz;
            m_padsize_y = int(size_y);
            m_padsize_z = int(size_z);

            const LibUtilities::PointsKey Ppad_y(m_padsize_y,LibUtilities::eFourierEvenlySpaced);
            const LibUtilities::BasisKey  Bpad_y(LibUtilities::eFourier,m_padsize_y,Ppad_y);

            const LibUtilities::PointsKey Ppad_z(m_padsize_z,LibUtilities::eFourierEvenlySpaced);
            const LibUtilities::BasisKey  Bpad_z(LibUtilities::eFourier,m_padsize_z,Ppad_z);

            m_paddingBasis_y = LibUtilities::BasisManager()[Bpad_y];
            m_paddingBasis_z = LibUtilities::BasisManager()[Bpad_z];

            StdRegions::StdQuadExp StdQuad(m_paddingBasis_y->GetBasisKey(),m_paddingBasis_z->GetBasisKey());

            StdRegions::StdMatrixKey matkey1(StdRegions::eFwdTrans,StdQuad.DetShapeType(),StdQuad);
            StdRegions::StdMatrixKey matkey2(StdRegions::eBwdTrans,StdQuad.DetShapeType(),StdQuad);

            MatFwdPAD = StdQuad.GetStdMatrix(matkey1);
            MatBwdPAD = StdQuad.GetStdMatrix(matkey2);
        }
    } //end of namespace
} //end of namespace
