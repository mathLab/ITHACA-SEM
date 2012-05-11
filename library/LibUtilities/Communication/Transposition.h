///////////////////////////////////////////////////////////////////////////////
//
// File Transposition.h
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
// Description: Manager class for homogeneous data transposition 
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_UTILITIES_HOMOGENEOUS1D_H
#define NEKTAR_LIB_UTILITIES_HOMOGENEOUS1D_H

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Foundations/Basis.h>

namespace Nektar
{
    namespace LibUtilities
    {
		enum TranspositionDir
		{
			eXYtoZ,
			eZtoXY,
			eXtoYZ,
			eYZtoX,
			eYZtoZY,
			eZYtoYZ,
			eXtoY,
			eYtoZ,
			eZtoX,
			eNoTrans
		};
		
        class Transposition
        {
		public:
			
			LIB_UTILITIES_EXPORT Transposition(const LibUtilities::BasisKey &HomoBasis0, 
											   LibUtilities::CommSharedPtr hcomm);
			
			LIB_UTILITIES_EXPORT Transposition(const LibUtilities::BasisKey &HomoBasis0,
											   const LibUtilities::BasisKey &HomoBasis1,
											   LibUtilities::CommSharedPtr hcomm);
			
			LIB_UTILITIES_EXPORT Transposition(const LibUtilities::BasisKey &HomoBasis0,
											   const LibUtilities::BasisKey &HomoBasis1,
											   const LibUtilities::BasisKey &HomoBasis2,
											   LibUtilities::CommSharedPtr hcomm);

			LIB_UTILITIES_EXPORT ~Transposition();
			
			LIB_UTILITIES_EXPORT unsigned int GetK(int i);
			
			LIB_UTILITIES_EXPORT Array<OneD, unsigned int> GetKs(void);
			
			LIB_UTILITIES_EXPORT unsigned int GetPlaneID(int i);
			
			LIB_UTILITIES_EXPORT Array<OneD, unsigned int> GetPlanesIDs(void);
			
			LIB_UTILITIES_EXPORT void Transpose(const Array<OneD,const NekDouble> &inarray,
												Array<OneD, NekDouble> &outarray, 
												bool UseNumMode = false, 
												TranspositionDir dir = eNoTrans);
			
		protected:
			
			CommSharedPtr         m_hcomm;
			
		private:
			
			LIB_UTILITIES_EXPORT void TransposeXYtoZ(const Array<OneD,const NekDouble> &inarray,
													 Array<OneD, NekDouble> &outarray, 
													 bool UseNumMode = false);
			
			LIB_UTILITIES_EXPORT void TransposeZtoXY(const Array<OneD,const NekDouble> &inarray,
													 Array<OneD, NekDouble> &outarray, 
													 bool UseNumMode = false);
			
			LIB_UTILITIES_EXPORT void TransposeXtoYZ(const Array<OneD,const NekDouble> &inarray,
													 Array<OneD, NekDouble> &outarray, 
													 bool UseNumMode = false);
			
			LIB_UTILITIES_EXPORT void TransposeYZtoX(const Array<OneD,const NekDouble> &inarray,
													 Array<OneD, NekDouble> &outarray, 
													 bool UseNumMode = false);
			
			LIB_UTILITIES_EXPORT void TransposeYZtoZY(const Array<OneD,const NekDouble> &inarray,
													 Array<OneD, NekDouble> &outarray, 
													 bool UseNumMode = false);
			
			LIB_UTILITIES_EXPORT void TransposeZYtoYZ(const Array<OneD,const NekDouble> &inarray,
													  Array<OneD, NekDouble> &outarray, 
													  bool UseNumMode = false);
			
			int m_num_homogeneous_directions;
			
			Array<OneD,int> m_num_points_per_proc;     // number of homogenous points on each process for each of the three homogeneous direction
			
			Array<OneD,int> m_num_homogeneous_points;  // total number homogeneous point per homogeneous direction
						
			Array<OneD,int> m_num_homogeneous_coeffs;  // total number homogeneous coefficients 
			
			Array<OneD,int> m_num_processes;
			
			int m_rank_id;                            // rank of the process
			
			Array<OneD, unsigned int> m_planes_IDs;    // IDs of the planes on the process
			
			Array<OneD, unsigned int> m_K;             // Fourier wave numbers associated with the planes on the process
			
			Array<OneD,int> m_SizeMap;                // MPI_Alltoallv map containing the size of send/rec buffer
			  
			Array<OneD,int> m_OffsetMap;              // MPI_Alltoallv map containing the offset of send/rec buffer in the global vector
			
        };
		
		typedef boost::shared_ptr<Transposition>      TranspositionSharedPtr;
	}
}
#endif
