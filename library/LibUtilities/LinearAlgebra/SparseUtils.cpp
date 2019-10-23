///////////////////////////////////////////////////////////////////////////////
//
// File: SparseUtils.hpp
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
// Description: common utility functions for sparse matrices
//
///////////////////////////////////////////////////////////////////////////////

#include <utility>
#include <map>

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/LinearAlgebra/SparseMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/SparseUtils.hpp>
#include <LibUtilities/LinearAlgebra/SparseMatrix.hpp>
#include <LibUtilities/LinearAlgebra/SparseDiagBlkMatrix.hpp>
#include <LibUtilities/LinearAlgebra/StorageSmvBsr.hpp>


namespace Nektar{

    void convertCooToBco(
                    const unsigned int  blkDim,
                    const COOMatType&   cooMat,
                          BCOMatType&   bcoMat)
    {
        COOMatTypeConstIt  entry;
        BCOMatTypeConstIt  blk;
        int rowcoord, localRow, blkRowCoord;
        int colcoord, localCol, blkColCoord;

        for(entry = cooMat.begin(); entry != cooMat.end(); entry++)
        {
            rowcoord = (entry->first).first;
            colcoord = (entry->first).second;

            blkRowCoord = rowcoord / blkDim;
            blkColCoord = colcoord / blkDim;

            CoordType blkCoords = std::make_pair(blkRowCoord,blkColCoord);
            blk = bcoMat.find(blkCoords);
            if (blk == bcoMat.end())
            {
                BCOEntryType b(blkDim*blkDim, 0.0);
                bcoMat[blkCoords] = b;
            }

            localRow    = rowcoord % blkDim;
            localCol    = colcoord % blkDim;

            // transpose it: NIST SpBLAS expects Fortran ordering
            // of dense subblocks
            const unsigned int localoffset = localRow + localCol*blkDim;
            (bcoMat[blkCoords])[localoffset] = entry->second;
        }
    }

    template<typename SparseStorageType>
    std::ostream& operator<<(std::ostream& os, const NekSparseMatrix<SparseStorageType>& rhs)
    {
        int oswidth = 9;
        int osprecision = 6;

        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            os << "[";
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                os.width(oswidth);
                os.precision(osprecision);
                os << rhs(i,j);
                if( j != rhs.GetColumns() - 1 )
                {
                    os << ", ";
                }
            }
            os << "]";
            if( i != rhs.GetRows()-1 )
            {
                os << std::endl;
            }
        }
        return os;
    }

    template<typename SparseStorageType>
    std::ostream& operator<<(std::ostream& os, const NekSparseDiagBlkMatrix<SparseStorageType>& rhs)
    {
        int oswidth = 9;
        int osprecision = 6;

        for(unsigned int i = 0; i < rhs.GetRows(); ++i)
        {
            os << "[";
            for(unsigned int j = 0; j < rhs.GetColumns(); ++j)
            {
                os.width(oswidth);
                os.precision(osprecision);
                os << rhs(i,j);
                if( j != rhs.GetColumns() - 1 )
                {
                    os << ", ";
                }
            }
            os << "]";
            if( i != rhs.GetRows()-1 )
            {
                os << std::endl;
            }
        }
        return os;
    }

    template std::ostream& operator<<(std::ostream& os, const NekSparseMatrix<StorageSmvBsr<NekDouble> >& rhs);
    template std::ostream& operator<<(std::ostream& os, const NekSparseDiagBlkMatrix<StorageSmvBsr<NekDouble> >& rhs);

} // namespace
