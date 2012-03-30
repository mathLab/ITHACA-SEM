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
// Description: An ExpList which is homogeneous in 1-direction
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpListHomogeneous1D.h>

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

        ExpListHomogeneous1D::ExpListHomogeneous1D(const LibUtilities::SessionReaderSharedPtr
                &pSession,const LibUtilities::BasisKey &HomoBasis, const NekDouble lhom, const bool useFFT, const bool dealiasing):
            ExpList(pSession),
            m_lhom(lhom),
            m_useFFT(useFFT),
		    m_dealiasing(dealiasing),
            m_homogeneous1DBlockMat(MemoryManager<Homo1DBlockMatrixMap>::AllocateSharedPtr())
        {
            ASSERTL2(HomoBasis != LibUtilities::NullBasisKey,"Homogeneous Basis is a null basis");
            m_homogeneousBasis = LibUtilities::BasisManager()[HomoBasis];
						
			SetParalleInfo();
			
			m_planes = Array<OneD,ExpListSharedPtr>(m_num_planes_per_proc);
            
            if(m_useFFT)
            {
                m_FFT = LibUtilities::GetNektarFFTFactory().CreateInstance("NekFFTW", m_num_fourier_points);
            }
			
			if(m_dealiasing)
			{
				ASSERTL0(m_num_processes == 1,"Remove dealiasing if you want to run in parallel");
				SetPaddingBase();
			}
		}


        /**
         * @param   In          ExpListHomogeneous1D object to copy.
         */
        ExpListHomogeneous1D::ExpListHomogeneous1D(const ExpListHomogeneous1D &In):
            ExpList(In,false),
            m_homogeneousBasis(In.m_homogeneousBasis),
            m_homogeneous1DBlockMat(In.m_homogeneous1DBlockMat),
            m_lhom(In.m_lhom),
            m_useFFT(In.m_useFFT),
		    m_FFT(In.m_FFT),
		    m_dealiasing(In.m_dealiasing),
		    m_padsize(In.m_padsize),
            MatBwdPAD(In.MatBwdPAD),
		    MatFwdPAD(In.MatFwdPAD),
            m_tmpIN(In.m_tmpIN),
            m_tmpOUT(In.m_tmpOUT)
        {
            m_planes = Array<OneD, ExpListSharedPtr>(In.m_planes.num_elements());
			
			SetParalleInfo();
        }

        /**
         * Destructor
         */
        ExpListHomogeneous1D::~ExpListHomogeneous1D()
        {
        }
		
        void ExpListHomogeneous1D::v_HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs)
        {
            // Forwards trans
            Homogeneous1DTrans(inarray,outarray,true, UseContCoeffs);
        }
		
        void ExpListHomogeneous1D::v_HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs)
        {
            // Backwards trans
            Homogeneous1DTrans(inarray,outarray,false, UseContCoeffs);
        }
		
		/**
		 * Dealiasing routine
		 */
		void ExpListHomogeneous1D::v_DealiasedProd(const Array<OneD, NekDouble> &inarray1,
												   const Array<OneD, NekDouble> &inarray2,
												   Array<OneD, NekDouble> &outarray, 
												   bool UseContCoeffs)
		{
			// inarray1 = first term of the product
			// inarray2 = second term of the product
			// dealiased product stored in outarray
						
			int npoints  = outarray.num_elements(); // number of total physical points
			int nplanes  = m_planes.num_elements(); // number of planes == number of Fourier modes = number of Fourier coeff
			int npencils = npoints/nplanes;         // number of pencils = numebr of physical points per plane

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
				HomogeneousFwdTrans(inarray1,V1,UseContCoeffs);
				HomogeneousFwdTrans(inarray2,V2,UseContCoeffs);
			}
			
			ShuffleIntoHomogeneous1DClosePacked(V1,ShufV1,false);
			ShuffleIntoHomogeneous1DClosePacked(V2,ShufV2,false);
			

			/////////////////////////////////////////////////////////////////////////////
			// Creating padded vectors for each pencil
			Array<OneD, NekDouble> PadV1_pencil_coeff(m_padsize,0.0);
			Array<OneD, NekDouble> PadV2_pencil_coeff(m_padsize,0.0);
			Array<OneD, NekDouble> PadRe_pencil_coeff(m_padsize,0.0);

			Array<OneD, NekDouble> PadV1_pencil_phys(m_padsize,0.0);
			Array<OneD, NekDouble> PadV2_pencil_phys(m_padsize,0.0);
			Array<OneD, NekDouble> PadRe_pencil_phys(m_padsize,0.0);

			NekVector<NekDouble> PadIN_V1(m_padsize,PadV1_pencil_coeff,eWrapper);
			NekVector<NekDouble> PadOUT_V1(m_padsize,PadV1_pencil_phys,eWrapper);

			NekVector<NekDouble> PadIN_V2(m_padsize,PadV2_pencil_coeff,eWrapper);
			NekVector<NekDouble> PadOUT_V2(m_padsize,PadV2_pencil_phys,eWrapper);

			NekVector<NekDouble> PadIN_Re(m_padsize,PadRe_pencil_phys,eWrapper);
			NekVector<NekDouble> PadOUT_Re(m_padsize,PadRe_pencil_coeff,eWrapper);

			//Looping on the pencils
			for(int i = 0 ; i< npencils ; i++)
			{
				//Copying the i-th pencil pf lenght N into a bigger pencil of lenght 2N
				//We are in Fourier space
				Vmath::Vcopy(nplanes,&(ShufV1[i*nplanes]),1,&(PadV1_pencil_coeff[0]),1);
				Vmath::Vcopy(nplanes,&(ShufV2[i*nplanes]),1,&(PadV2_pencil_coeff[0]),1);
				//Moving to physical space using the padded system
				PadOUT_V1 = (*MatBwdPAD)*PadIN_V1;
				PadOUT_V2 = (*MatBwdPAD)*PadIN_V2;

				//Perfroming the vectors multiplication in physical space on the padded system
				Vmath::Vmul(m_padsize,PadV1_pencil_phys,1,PadV2_pencil_phys,1,PadRe_pencil_phys,1);

				//Moving back the result (V1*V2)_phys in Fourier space, padded system
				PadOUT_Re = (*MatFwdPAD)*PadIN_Re;

				//Copying the first half of the padded pencil in the full vector (Fourier space)
				Vmath::Vcopy(nplanes,&(PadRe_pencil_coeff[0]),1,&(ShufV1V2[i*nplanes]),1);
				
			}
			
			if(m_WaveSpace)
			{

				//Unshuffle the dealiased result vector (still in Fourier space) in the original ordering
				UnshuffleFromHomogeneous1DClosePacked(ShufV1V2,outarray,false);
				
			}
			else 
			{
				//Unshuffle the dealiased result vector (still in Fourier space) in the original ordering
				UnshuffleFromHomogeneous1DClosePacked(ShufV1V2,V1V2,false);
				//Moving the results in physical space for the output
				HomogeneousBwdTrans(V1V2,outarray,UseContCoeffs);
			}

		}
		

		/**
		 * Forward transform
		 */
        void ExpListHomogeneous1D::v_FwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            
            for(int n = 0; n < m_num_planes_per_proc; ++n)
            {
                m_planes[n]->FwdTrans(inarray+cnt, tmparray = outarray + cnt1,
                                      UseContCoeffs);
                cnt   += m_planes[n]->GetTotPoints();

                if(UseContCoeffs)
                {
                    cnt1  += m_planes[n]->GetContNcoeffs();
                }
                else
                {
                    cnt1  += m_planes[n]->GetNcoeffs();
                }
            }
            if(!m_WaveSpace)
			{
				HomogeneousFwdTrans(outarray,outarray,UseContCoeffs);
			}
        }

		/**
		 * Forward transform element by element
		 */
        void ExpListHomogeneous1D::v_FwdTrans_IterPerExp(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            
			//spectral element FwdTrans plane by plane
            for(int n = 0; n < m_num_planes_per_proc; ++n)
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
		 * Backward transform
		 */
        void ExpListHomogeneous1D::v_BwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            
            for(int n = 0; n < m_num_planes_per_proc; ++n)
            {
                m_planes[n]->BwdTrans(inarray+cnt, tmparray = outarray + cnt1,
                                      UseContCoeffs);
                if(UseContCoeffs)
                {
                    cnt    += m_planes[n]->GetContNcoeffs();
                }
                else
                {
                    cnt    += m_planes[n]->GetNcoeffs();
                }
                cnt1   += m_planes[n]->GetTotPoints();
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
			
            for(int n = 0; n < m_num_planes_per_proc; ++n)
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
        void ExpListHomogeneous1D::v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;

            for(int n = 0; n < m_num_planes_per_proc; ++n)
            {
                m_planes[n]->IProductWRTBase(inarray+cnt, tmparray = outarray + cnt1,UseContCoeffs);

                if(UseContCoeffs)
                {
                    cnt1    += m_planes[n]->GetContNcoeffs();
                }
                else
                {
                    cnt1    += m_planes[n]->GetNcoeffs();
                }
                cnt   += m_planes[n]->GetTotPoints();
            }
        }
		
		/**
		 * Inner product element by element
		 */
		void ExpListHomogeneous1D::v_IProductWRTBase_IterPerExp(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        { 
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
			
            for(int n = 0; n < m_num_planes_per_proc; ++n)
            {
                m_planes[n]->IProductWRTBase_IterPerExp(inarray+cnt, tmparray = outarray + cnt1);
				
				cnt1  += m_planes[n]->GetNcoeffs();
				cnt   += m_planes[n]->GetTotPoints();
            } 
        }
		
		/**
		 * Homogeneous transform Bwd/Fwd (MVM and FFT)
		 */
        void ExpListHomogeneous1D::Homogeneous1DTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool IsForwards, bool UseContCoeffs)
        {
            if(m_useFFT)
            {				
				int num_dofs             = inarray.num_elements();
				int num_points_per_plane = num_dofs/m_num_planes_per_proc;
				int num_dfts_per_proc    = num_points_per_plane/m_num_processes + (num_points_per_plane%m_num_processes > 0);
		
                Array<OneD, NekDouble> fft_in(num_dfts_per_proc*m_num_fourier_points);
                Array<OneD, NekDouble> fft_out(num_dfts_per_proc*m_num_fourier_points);
		
                ShuffleIntoHomogeneous1DClosePacked(inarray,fft_in,false);
		
                if(IsForwards)
                {
                    for(int i = 0 ; i < num_dfts_per_proc ; i++)
                    {
                        m_FFT->FFTFwdTrans(m_tmpIN = fft_in + i*m_num_fourier_points, m_tmpOUT = fft_out + i*m_num_fourier_points);
                    }
                }
                else 
                {
                    for(int i = 0 ; i < num_dfts_per_proc ; i++)
                    {
                        m_FFT->FFTBwdTrans(m_tmpIN = fft_in + i*m_num_fourier_points, m_tmpOUT = fft_out + i*m_num_fourier_points);
                    }
                }
		
                UnshuffleFromHomogeneous1DClosePacked(fft_out,outarray,false);
            }
            else 
            {
                DNekBlkMatSharedPtr blkmat;
		
                if(inarray.num_elements() == m_npoints) //transform phys space
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
                        blkmat = GetHomogeneous1DBlockMatrix(eForwardsCoeffSpace1D,UseContCoeffs);
                    }
                    else
                    {
                        blkmat = GetHomogeneous1DBlockMatrix(eBackwardsCoeffSpace1D,UseContCoeffs);
                    }
                }
		
                int nrows = blkmat->GetRows();
                int ncols = blkmat->GetColumns();
		
                Array<OneD, NekDouble> sortedinarray(ncols);
                Array<OneD, NekDouble> sortedoutarray(nrows);
		
                ShuffleIntoHomogeneous1DClosePacked(inarray,sortedinarray,!IsForwards);
                
                // Create NekVectors from the given data arrays
                NekVector<NekDouble> in (ncols,sortedinarray,eWrapper);
                NekVector<NekDouble> out(nrows,sortedoutarray,eWrapper);
		
                // Perform matrix-vector multiply.
                out = (*blkmat)*in;
				
                UnshuffleFromHomogeneous1DClosePacked(sortedoutarray,outarray,IsForwards);
            }
        }

		/*
		 * Shuffle routine into Fourier pencils
		 */
        void ExpListHomogeneous1D::ShuffleIntoHomogeneous1DClosePacked(
                              const Array<OneD, const NekDouble> &inarray,
                              Array<OneD, NekDouble> &outarray,
                              bool UseNumModes)
        {
			int i,j;
            int num_dofs             = inarray.num_elements();
			int num_points_per_plane = num_dofs/m_num_planes_per_proc;
			int num_dfts_per_proc    = num_points_per_plane/m_num_processes + (num_points_per_plane%m_num_processes > 0);
			int copy_len             = num_dfts_per_proc;
			
			Array<OneD, int> SizeMap(m_num_processes,0);
			Array<OneD, int> OffsetMap(m_num_processes,0);
			
			Array< OneD, NekDouble> tmp_outarray(num_dfts_per_proc*m_num_fourier_points,0.0);
			
			for(i = 0; i < m_num_processes ; i++)
			{
				if(i == m_num_processes-1)
				{
					copy_len = num_dfts_per_proc - ((num_dfts_per_proc*m_num_processes) - num_points_per_plane);
				}
				
				for(j = 0; j < m_num_planes_per_proc;j++)
				{
					Vmath::Vcopy(copy_len,
								 &(inarray[i*num_dfts_per_proc+j*num_points_per_plane]),1,
								 &(outarray[i*num_dfts_per_proc*m_num_planes_per_proc+j*num_dfts_per_proc]),1);
				}
				
				SizeMap[i]   = num_dfts_per_proc*m_num_planes_per_proc;
				OffsetMap[i] = i*num_dfts_per_proc*m_num_planes_per_proc;				
			}
			
			m_comm->GetColumnComm()->AlltoAllv(outarray,SizeMap,OffsetMap,tmp_outarray,SizeMap,OffsetMap);
			
			// reshuffle
			int packed_len;
			
            if(UseNumModes)
            {
                packed_len = m_num_fourier_coeffs;
            }
            else
            {
                packed_len = m_num_fourier_points;
            }

            ASSERTL1(&inarray[0] != &outarray[0],"Inarray and outarray cannot be the same");

            for(i = 0; i < packed_len; ++i)
            {
                Vmath::Vcopy(num_dfts_per_proc,&(tmp_outarray[i*num_dfts_per_proc]),1,
                             &(outarray[i]),packed_len);
            }
        }

		/*
		 * UnShuffle routine into spectral elements planes
		 */
        void ExpListHomogeneous1D::UnshuffleFromHomogeneous1DClosePacked(
                              const Array<OneD, const NekDouble> &inarray,
                              Array<OneD, NekDouble> &outarray,
                              bool UseNumModes)
        {
			int i,j;
            int num_dofs             = inarray.num_elements();
			int num_points_per_plane = num_dofs/m_num_planes_per_proc;
			int num_dfts_per_proc    = num_points_per_plane/m_num_processes + (num_points_per_plane%m_num_processes > 0);
			int copy_len             = num_dfts_per_proc;
			
			Array<OneD, int> SizeMap(m_num_processes,0);
			Array<OneD, int> OffsetMap(m_num_processes,0);
			
			Array< OneD, NekDouble> tmp_inarray(num_dfts_per_proc*m_num_fourier_points,0.0);
			Array< OneD, NekDouble> tmp_outarray(num_dfts_per_proc*m_num_fourier_points,0.0);
			
			int packed_len;
			
            if(UseNumModes)
            {
                packed_len = m_num_fourier_coeffs;
            }
            else
            {
                packed_len = m_num_fourier_points;
            }
			
            ASSERTL1(&inarray[0] != &outarray[0],"Inarray and outarray cannot be the same");
			
            for(i = 0; i < packed_len; ++i)
            {
                Vmath::Vcopy(num_dfts_per_proc,&(inarray[i]),packed_len,
                             &(tmp_inarray[i*num_dfts_per_proc]),1);
            }
			
			for(i = 0; i < m_num_processes ; i++)
			{				
				SizeMap[i]   = num_dfts_per_proc*m_num_planes_per_proc;
				OffsetMap[i] = i*num_dfts_per_proc*m_num_planes_per_proc;				
			}
			
			m_comm->GetColumnComm()->AlltoAllv(tmp_inarray,SizeMap,OffsetMap,tmp_outarray,SizeMap,OffsetMap);
			
			for(i = 0; i < m_num_processes;i++)
			{	
				if(i == m_num_processes-1)
				{
					copy_len = num_dfts_per_proc - ((num_dfts_per_proc*m_num_processes)-num_points_per_plane);
				}
				
				for(j = 0; j < m_num_planes_per_proc;j++)
				{
					Vmath::Vcopy(copy_len,
								 &(tmp_outarray[i*num_dfts_per_proc*m_num_planes_per_proc+j*num_dfts_per_proc]),1,
								 &(outarray[i*num_dfts_per_proc+j*num_points_per_plane]),1);
				}
			}
        }

        DNekBlkMatSharedPtr ExpListHomogeneous1D::GetHomogeneous1DBlockMatrix(Homogeneous1DMatType mattype, bool UseContCoeffs) const
        {
            Homo1DBlockMatrixMap::iterator matrixIter = m_homogeneous1DBlockMat->find(mattype);

            if(matrixIter == m_homogeneous1DBlockMat->end())
            {
                return ((*m_homogeneous1DBlockMat)[mattype] =
                        GenHomogeneous1DBlockMatrix(mattype,UseContCoeffs));
            }
            else
            {
                return matrixIter->second;
            }
        }


        DNekBlkMatSharedPtr ExpListHomogeneous1D::GenHomogeneous1DBlockMatrix(Homogeneous1DMatType mattype, bool UseContCoeffs) const
        {
            DNekMatSharedPtr    loc_mat;
            DNekBlkMatSharedPtr BlkMatrix;
			int n_exp = 0;
			int num_dfts_per_proc = 0;
			

            if((mattype == eForwardsCoeffSpace1D)
               ||(mattype == eBackwardsCoeffSpace1D)) // will operate on m_coeffs
            {
                if(UseContCoeffs)
                {
                    n_exp = m_planes[0]->GetContNcoeffs();
					
                }
                else
                {
                    n_exp = m_planes[0]->GetNcoeffs();
                }
            }
            else
            {
                n_exp = m_planes[0]->GetTotPoints(); // will operatore on m_phys
            }
			
			num_dfts_per_proc = n_exp/m_num_processes + (n_exp%m_num_processes > 0);

            Array<OneD,unsigned int> nrows(num_dfts_per_proc);
            Array<OneD,unsigned int> ncols(num_dfts_per_proc);

            if((mattype == eForwardsCoeffSpace1D)||(mattype == eForwardsPhysSpace1D))
            {
                nrows = Array<OneD, unsigned int>(num_dfts_per_proc,m_num_fourier_coeffs);
                ncols = Array<OneD, unsigned int>(num_dfts_per_proc,m_num_fourier_points);
            }
            else
            {
                nrows = Array<OneD, unsigned int>(num_dfts_per_proc,m_num_fourier_points);
                ncols = Array<OneD, unsigned int>(num_dfts_per_proc,m_num_fourier_coeffs);
            }

            MatrixStorage blkmatStorage = eDIAGONAL;
            BlkMatrix = MemoryManager<DNekBlkMat>
                ::AllocateSharedPtr(nrows,ncols,blkmatStorage);

            StdRegions::StdSegExp StdSeg(m_homogeneousBasis->GetBasisKey());

            if((mattype == eForwardsCoeffSpace1D)||(mattype == eForwardsPhysSpace1D))
            {
                StdRegions::StdMatrixKey matkey(StdRegions::eFwdTrans,
                                                StdSeg.DetExpansionType(),
                                                StdSeg);

                loc_mat = StdSeg.GetStdMatrix(matkey);
				
            }
            else
            {
                StdRegions::StdMatrixKey matkey(StdRegions::eBwdTrans,
                                                StdSeg.DetExpansionType(),
                                                StdSeg);

                loc_mat = StdSeg.GetStdMatrix(matkey);
				
            }

            // set up array of block matrices.
            for(int i = 0; i < num_dfts_per_proc; ++i)
            {
                BlkMatrix->SetBlock(i,i,loc_mat);
            }

            return BlkMatrix;
        }

        std::vector<SpatialDomains::FieldDefinitionsSharedPtr> ExpListHomogeneous1D::v_GetFieldDefinitions()
        {
            std::vector<SpatialDomains::FieldDefinitionsSharedPtr> returnval;
            
			// Set up Homogeneous length details.
            Array<OneD,LibUtilities::BasisSharedPtr> HomoBasis(1,m_homogeneousBasis);
            
			std::vector<NekDouble> HomoLen;
            HomoLen.push_back(m_lhom);
			
			std::vector<unsigned int> PlanesIDs;
			for(int i = 0; i < m_num_planes_per_proc; i++)
			{
				PlanesIDs.push_back(m_planes_IDs[i]);
			}

            m_planes[0]->GeneralGetFieldDefinitions(returnval, 1, HomoBasis, HomoLen, PlanesIDs);
            
			return returnval;
        }

        void  ExpListHomogeneous1D::v_GetFieldDefinitions(std::vector<SpatialDomains::FieldDefinitionsSharedPtr> &fielddef)
        {
            // Set up Homogeneous length details.
            Array<OneD,LibUtilities::BasisSharedPtr> HomoBasis(1,m_homogeneousBasis);
            
			std::vector<NekDouble> HomoLen;
            HomoLen.push_back(m_lhom);
			
			std::vector<unsigned int> PlanesIDs;
			for(int i = 0; i < m_num_planes_per_proc; i++)
			{
				PlanesIDs.push_back(m_planes_IDs[i]);
			}

             // enforce NumHomoDir == 1 by direct call
            m_planes[0]->GeneralGetFieldDefinitions(fielddef,1, HomoBasis,HomoLen,PlanesIDs);
        }


        void ExpListHomogeneous1D::v_AppendFieldData(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, Array<OneD, NekDouble> &coeffs)
        {
            int i,n;
            int ncoeffs_per_plane = m_planes[0]->GetNcoeffs();

            // Determine mapping from element ids to location in
            // expansion list
            map<int, int> ElmtID_to_ExpID;
            for(i = 0; i < m_planes[0]->GetExpSize(); ++i)
            {
                ElmtID_to_ExpID[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
            }

            for(i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                int eid     = ElmtID_to_ExpID[fielddef->m_elementIDs[i]];
                int datalen = (*m_exp)[eid]->GetNcoeffs();

                for(n = 0; n < m_num_planes_per_proc; ++n)
                {
                    fielddata.insert(fielddata.end(),&coeffs[m_coeff_offset[eid]+n*ncoeffs_per_plane],&coeffs[m_coeff_offset[eid]+n*ncoeffs_per_plane]+datalen);
                }
            }
        }
		
		void ExpListHomogeneous1D::v_AppendFieldData(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata)
        {
           v_AppendFieldData(fielddef,fielddata,m_coeffs);
        }

        //Extract the data in fielddata into the m_coeff list
        void ExpListHomogeneous1D::v_ExtractDataToCoeffs(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, std::string &field)
        {
            int i,n;
            int offset = 0;
            int nzmodes; 
            int datalen = fielddata.size()/fielddef->m_fields.size();
            int ncoeffs_per_plane = m_planes[0]->GetNcoeffs();
            

            for(i = 0; i < fielddef->m_basis.size(); ++i)
            {
                if(fielddef->m_basis[i] == m_homogeneousBasis->GetBasisType())
                {
                    nzmodes = fielddef->m_homogeneousZIDs.size();
                    break;
                }
            }
            ASSERTL1(i != fielddef->m_basis.size()," Failed to determine number of modes");
            
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
            for(i = 0; i < m_planes[0]->GetExpSize(); ++i)
            {
                ElmtID_to_ExpID[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
            }

            int modes_offset = 0;
            Array<OneD, NekDouble> coeff_tmp;             
            for(i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                int eid = ElmtID_to_ExpID[fielddef->m_elementIDs[i]];
                int datalen = (*m_exp)[eid]->CalcNumberOfCoefficients(fielddef->m_numModes,modes_offset);
                if(fielddef->m_uniOrder == true) // reset modes_offset to zero
                {
                    modes_offset = 0;
                }

                for(n = 0; n < nzmodes; ++n)
                {
                    if(datalen == (*m_exp)[eid]->GetNcoeffs())
                    {
                        Vmath::Vcopy(datalen,&fielddata[offset],1,&m_coeffs[m_coeff_offset[eid]+ n*ncoeffs_per_plane],1);
                    }
                    else // unpack data to new order
                    {
                        (*m_exp)[eid]->ExtractDataToCoeffs(fielddata, offset, fielddef->m_numModes,modes_offset,coeff_tmp = m_coeffs + m_coeff_offset[eid] + n*ncoeffs_per_plane);
                    }
                    
                    offset += datalen;
                }
            }
        }
		
		
		
        //Extract the data in fielddata into the m_coeff list (for 2D files into 3D cases)
        void ExpListHomogeneous1D::v_ExtractDataToCoeffs(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, std::string &field, bool BaseFlow3D)
        {
            int i,n;
            int offset = 0;
            int nzmodes = m_homogeneousBasis->GetNumModes();
            int datalen = fielddata.size()/fielddef->m_fields.size();
            int ncoeffs_per_plane = m_planes[0]->GetNcoeffs();
		
			
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
            for(i = 0; i < m_planes[0]->GetExpSize(); ++i)
            {
                ElmtID_to_ExpID[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
            }
			

            Vmath::Vcopy(datalen,&fielddata[offset],1,&m_coeffs[0],1);
                          
        }
		


        /**
         * Write Tecplot Files Header
         * @param   outfile Output file name.
         * @param   var                 variables names
         */
        void ExpListHomogeneous1D::v_WriteTecplotHeader(std::ofstream &outfile, std::string var)
        {
            if(GetExp(0)->GetCoordim() == 1)
            {
                outfile << "Variables = x, y";
            }
            else
            {
                outfile << "Variables = x, y, z";
            }
            outfile << ", "<< var << std::endl << std::endl;
        }

        /**
         * Write Tecplot Files Field
         * @param   outfile    Output file name.
         * @param   expansion  Expansion that is considered
         */
        void ExpListHomogeneous1D::v_WriteTecplotField(std::ofstream &outfile, int expansion)
        {
            int npoints_per_plane = m_planes[0]->GetTotPoints();

            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                (*m_exp)[expansion]->SetPhys(m_phys+m_phys_offset[expansion]+
                                             n*npoints_per_plane);
                (*m_exp)[expansion]->WriteTecplotField(outfile);
            }
        }

        void ExpListHomogeneous1D::v_WriteVtkPieceData(std::ofstream &outfile, int expansion,
                                        std::string var)
        {
            int i;
            int nq = (*m_exp)[expansion]->GetTotPoints();
            int npoints_per_plane = m_planes[0]->GetTotPoints();

            // printing the fields of that zone
            outfile << "        <DataArray type=\"Float32\" Name=\""
                    << var << "\">" << endl;
            outfile << "          ";
            for (int n = 0; n < m_planes.num_elements(); ++n)
            {
                const Array<OneD, NekDouble> phys = m_phys + m_phys_offset[expansion] + n*npoints_per_plane;
                for(i = 0; i < nq; ++i)
                {
                    outfile << (fabs(phys[i]) < NekConstants::kNekZeroTol ? 0 : phys[i]) << " ";
                }
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
        }
		
		void ExpListHomogeneous1D::v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
											   Array<OneD, NekDouble> &out_d0,
											   Array<OneD, NekDouble> &out_d1, 
											   Array<OneD, NekDouble> &out_d2, bool UseContCoeffs)
		
		{
			int nT_pts = inarray.num_elements();          //number of total points = n. of Fourier points * n. of points per plane (nT_pts)
			int nP_pts = nT_pts/m_num_planes_per_proc;    //number of points per plane = n of Fourier transform required (nP_pts)
			
			Array<OneD, NekDouble> temparray(nT_pts);
			Array<OneD, NekDouble> outarray(nT_pts);
			Array<OneD, NekDouble> tmp1;
			Array<OneD, NekDouble> tmp2;
			Array<OneD, NekDouble> tmp3;            
			
			for(int i = 0; i < m_num_planes_per_proc; i++)
			{
				m_planes[i]->PhysDeriv(tmp1 = inarray + i*nP_pts ,tmp2 = out_d0 + i*nP_pts , tmp3 = out_d1 + i*nP_pts );
			}
			
			if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourier || m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierSingleMode)
			{
				if(m_WaveSpace)
				{
					temparray = inarray;
				}
				else 
				{ 
					HomogeneousFwdTrans(inarray,temparray,UseContCoeffs);
				}

				NekDouble sign = -1.0;
				NekDouble beta;
				
				for(int i = 0; i < m_num_planes_per_proc; i++)
				{
					beta = sign*2*M_PI*m_K[i]/m_lhom;
					
					Vmath::Smul(nP_pts,beta,tmp1 = temparray + i*nP_pts,1,tmp2 = outarray + (i-int(sign))*nP_pts,1);
					
					sign = -1.0*sign;
				}
				
				if(m_WaveSpace)
				{
					out_d2 = outarray;
				}
				else 
				{
					HomogeneousBwdTrans(outarray,out_d2,UseContCoeffs);
				}
			}
			else 
			{
				ASSERTL0(m_num_processes == 1,"Parallelisation in the homogeneous direction implemented just for Fourier basis");
				
				if(m_WaveSpace)
				{

					ASSERTL0(false,"Semi-phyisical time-stepping not implemented yet for non-Fourier basis");
				}
				else 
				{
					StdRegions::StdSegExp StdSeg(m_homogeneousBasis->GetBasisKey());
					
					ShuffleIntoHomogeneous1DClosePacked(inarray,temparray,false);
					
					for(int i = 0; i < nP_pts; i++)
					{
						StdSeg.PhysDeriv(tmp1 = temparray + i*m_num_planes_per_proc, tmp2 = outarray + i*m_num_planes_per_proc);
					}
					
					UnshuffleFromHomogeneous1DClosePacked(outarray,out_d2,false);
					
					Vmath::Smul(nT_pts,2.0/m_lhom,out_d2,1,out_d2,1);
					
				}
			}
		}
		
		void ExpListHomogeneous1D::v_PhysDeriv(Direction edir,
											   const Array<OneD, const NekDouble> &inarray,
											   Array<OneD, NekDouble> &out_d, bool UseContCoeffs)
		
		{
			int nT_pts = inarray.num_elements();        //number of total points = n. of Fourier points * n. of points per plane (nT_pts)
			int nP_pts = nT_pts/m_num_planes_per_proc;  //number of points per plane = n of Fourier transform required (nP_pts)
			
			int dir= (int)edir;
			
			Array<OneD, NekDouble> temparray(nT_pts);
			Array<OneD, NekDouble> outarray(nT_pts);
			Array<OneD, NekDouble> tmp1;
			Array<OneD, NekDouble> tmp2;
            			
			if (dir < 2)
			{
				for(int i=0; i<m_num_planes_per_proc; i++)
				{
					m_planes[i]->PhysDeriv(edir, tmp1 = inarray + i*nP_pts ,tmp2 = out_d + i*nP_pts);
				}
			}
			else
			{
				if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourier || m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierSingleMode)
				{
					if(m_WaveSpace)
					{
						temparray = inarray;
					}
					else 
					{ 
						HomogeneousFwdTrans(inarray,temparray,UseContCoeffs);
					}
					
					NekDouble sign = -1.0;
					NekDouble beta;
					
					for(int i = 0; i < m_num_planes_per_proc; i++)
					{
						beta = sign*2*M_PI*m_K[i]/m_lhom;
						
						Vmath::Smul(nP_pts,beta,tmp1 = temparray + i*nP_pts,1,tmp2 = outarray + (i-int(sign))*nP_pts,1);
						
						sign = -1.0*sign;
					}
					if(m_WaveSpace)
					{
						out_d = outarray;
					}
					else 
					{
						HomogeneousBwdTrans(outarray,out_d,UseContCoeffs);
					}
				}
				else 
				{
					ASSERTL0(m_num_processes == 1,"Parallelisation in the homogeneous direction implemented just for Fourier basis");
					
					if(m_WaveSpace)
					{
						ASSERTL0(false,"Semi-phyisical time-stepping not implemented yet for non-Fourier basis");
					}
					else 
					{
						StdRegions::StdSegExp StdSeg(m_homogeneousBasis->GetBasisKey());
						
						ShuffleIntoHomogeneous1DClosePacked(inarray,temparray,false);
						
						for(int i = 0; i < nP_pts; i++)
						{
							StdSeg.PhysDeriv(tmp1 = temparray + i*m_num_planes_per_proc, tmp2 = outarray + i*m_num_planes_per_proc);
						}
						
						UnshuffleFromHomogeneous1DClosePacked(outarray,out_d,false);
						
						Vmath::Smul(nT_pts,2.0/m_lhom,out_d,1,out_d,1);
					}
				}
			}
		}
        
        void ExpListHomogeneous1D::PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD, NekDouble> &out_d0,
                                             Array<OneD, NekDouble> &out_d1, 
                                             Array<OneD, NekDouble> &out_d2, bool UseContCoeffs)
            
        {
            v_PhysDeriv(inarray,out_d0,out_d1,out_d2,UseContCoeffs);
        }
	
        void ExpListHomogeneous1D::PhysDeriv(Direction edir,
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD, NekDouble> &out_d, bool UseContCoeffs)
        {
            v_PhysDeriv(edir,inarray,out_d,UseContCoeffs);
        }
		
		/*
		 * Setting the Padding base for dealisaing
		 */
		void ExpListHomogeneous1D::SetPaddingBase(void)
		{
			NekDouble size = 1.5*m_homogeneousBasis->GetNumPoints();
			m_padsize = int(size);
			
			const LibUtilities::PointsKey Ppad(m_padsize,LibUtilities::eFourierEvenlySpaced);
            const LibUtilities::BasisKey  Bpad(LibUtilities::eFourier,m_padsize,Ppad);
            
            m_paddingBasis = LibUtilities::BasisManager()[Bpad];
			
			StdRegions::StdSegExp StdSeg(m_paddingBasis->GetBasisKey());
			
			StdRegions::StdMatrixKey matkey1(StdRegions::eFwdTrans,StdSeg.DetExpansionType(),StdSeg);
			StdRegions::StdMatrixKey matkey2(StdRegions::eBwdTrans,StdSeg.DetExpansionType(),StdSeg);
			
			MatFwdPAD = StdSeg.GetStdMatrix(matkey1);
			MatBwdPAD = StdSeg.GetStdMatrix(matkey2);
		}
		
		/*
		 * Setting parallelisation dimension and offsets
		 */
		void ExpListHomogeneous1D::SetParalleInfo(void)
		{
			m_num_fourier_points  = m_homogeneousBasis->GetNumPoints();
			m_num_fourier_coeffs  = m_homogeneousBasis->GetNumModes();
			m_num_processes       = m_comm->GetColumnComm()->GetSize();
			m_num_planes_per_proc = m_num_fourier_points/m_num_processes;
			m_rank_id             = m_comm->GetColumnComm()->GetRank();
			
			m_planes_IDs = Array<OneD,unsigned int>(m_num_planes_per_proc);
			m_K          = Array<OneD,unsigned int>(m_num_planes_per_proc);
			
			for(int i = 0 ; i < m_num_planes_per_proc ; i++)
			{
				m_planes_IDs[i] = m_rank_id*m_num_planes_per_proc + i;
				m_K[i]          = m_planes_IDs[i]/2;
			}
			
			if(m_homogeneousBasis->GetBasisType() == LibUtilities::eFourierSingleMode)
			{
				m_K[0] = 1;
				m_K[1] = 1;
			}
		}
		
		Array<OneD, unsigned int> ExpListHomogeneous1D::v_GetZIDs(void)
		{
			return m_planes_IDs;
		}
    } //end of namespace
} //end of namespace


/**
* $Log: v $
*
**/

