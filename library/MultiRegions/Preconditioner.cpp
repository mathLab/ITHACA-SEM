///////////////////////////////////////////////////////////////////////////////
//
// File Preconditioner.cpp
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
// Description: Preconditioner definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/Preconditioner.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <LocalRegions/MatrixKey.h>
#include <math.h>

namespace Nektar
{
    namespace MultiRegions
    {
        std::string Preconditioner::lookupIds[7] = {
                LibUtilities::SessionReader::RegisterEnumValue("Preconditioner","Null",MultiRegions::eNull),
                LibUtilities::SessionReader::RegisterEnumValue("Preconditioner","Diagonal",MultiRegions::eDiagonal),
                LibUtilities::SessionReader::RegisterEnumValue("Preconditioner","InverseLinear",MultiRegions::eInverseLinear),
                LibUtilities::SessionReader::RegisterEnumValue("Preconditioner","LowEnergy",MultiRegions::eLowEnergy),
                LibUtilities::SessionReader::RegisterEnumValue("Preconditioner","LinearLowEnergy",MultiRegions::eLinearLowEnergy),
                LibUtilities::SessionReader::RegisterEnumValue("Preconditioner","Block",MultiRegions::eBlock),
                LibUtilities::SessionReader::RegisterEnumValue("Preconditioner","LocalLowEnergy",MultiRegions::eLocalLowEnergy),
        };
        std::string Preconditioner::def = LibUtilities::SessionReader::RegisterDefaultSolverInfo("Preconditioner","Diagonal");

        /**
         * @class Preconditioner
         *
         * This class implements preconditioning for the conjugate 
	 * gradient matrix solver.
	 */

         Preconditioner::Preconditioner(
                         const boost::shared_ptr<GlobalLinSys> &plinsys,
	                 const AssemblyMapSharedPtr &pLocToGloMap):
	   m_linsys(plinsys),
           m_preconType(pLocToGloMap->GetPreconType()),
           m_locToGloMap(pLocToGloMap)
         {
         }

        /**
         *
         */
        PreconFactory& GetPreconFactory()
        {
            typedef Loki::SingletonHolder<PreconFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy > Type;
            return Type::Instance();
        }

        void Preconditioner::v_InitObject()
        {
            NEKERROR(ErrorUtil::efatal,"Method does not exist" );
	}



        void Preconditioner::v_DoPreconditioner(
                      const Array<OneD, NekDouble>& pInput,
		      Array<OneD, NekDouble>& pOutput)
        {
            NEKERROR(ErrorUtil::efatal,"Method does not exist" );
	}

        /**
         * \brief Get the transformation matrix \f$\mathbf{R}\f$
         */
        const DNekMatSharedPtr& Preconditioner::v_GetTransformationMatrix() const
	{
            NEKERROR(ErrorUtil::efatal,"Method does not exist" );
            return NullDNekMatSharedPtr;
	}

        /**
         * \brief Get the transposed transformation matrix \f$\mathbf{R}^{T}\f$
         */
        const DNekMatSharedPtr& Preconditioner::v_GetTransposedTransformationMatrix() const
	{
            NEKERROR(ErrorUtil::efatal,"Method does not exist" );
            return NullDNekMatSharedPtr;
	}

        /**
         *  Performs global assembly of diagonal entries
         *  to global schur complement matrix.
         */
        Array<OneD, NekDouble> Preconditioner::AssembleStaticCondGlobalDiagonals()
        {
            int nGlobalBnd = m_locToGloMap->GetNumGlobalBndCoeffs();
            int nDirBnd = m_locToGloMap->GetNumGlobalDirBndCoeffs();
            int rows = nGlobalBnd - nDirBnd;

            DNekScalBlkMatSharedPtr loc_mat;
            DNekScalMatSharedPtr    bnd_mat;
            int sign1, sign2, gid1, gid2, i, j, n, cnt;
            Array<OneD, NekDouble> diagonals(rows,0.0);

            // Extract diagonal contributions of globally assembled
            // schur complement matrix
            for(cnt = n = 0; n < m_linsys.lock()->GetNumBlocks(); ++n)
            {
                //Get statically condensed local matrix
                loc_mat = (m_linsys.lock())->GetStaticCondBlock(n);

                //Extract boundary block
                bnd_mat=loc_mat->GetBlock(0,0);

                //offset by number of rows
                int bnd_lda = bnd_mat->GetRows();

                for(i = 0; i < bnd_lda; ++i)
                {
                    gid1  = m_locToGloMap->GetLocalToGlobalBndMap (cnt + i) - nDirBnd;
                    sign1 = m_locToGloMap->GetLocalToGlobalBndSign(cnt + i);

                    if(gid1 >= 0)
                    {
                        for(j = 0; j < bnd_lda; ++j)
                        {
                            gid2  = m_locToGloMap->GetLocalToGlobalBndMap(cnt+j) - nDirBnd;
                            sign2 = m_locToGloMap->GetLocalToGlobalBndSign(cnt+j);

                            if(gid2 == gid1)
                            {
                                diagonals[gid1] += sign1*sign2*(*bnd_mat)(i,j);
                            }
                        }
                   }
                }
                cnt += bnd_lda;
            }
            return diagonals;
        }
    }
}






