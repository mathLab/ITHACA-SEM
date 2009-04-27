///////////////////////////////////////////////////////////////////////////////
//
// File: SparseStandardMatrix.hpp
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
// Description: sparse matrix class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSE_STANDARD_MATRIX_HPP
#define NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSE_STANDARD_MATRIX_HPP

#include <map>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <boost/lexical_cast.hpp>

namespace Nektar
{    
    template<typename DataType>
    class NekSparseMatrix
    {

    public:
        typedef std::pair<  int,  int> CoordType;
        typedef std::map< CoordType, DataType >        COOMatType;
        typedef typename COOMatType::const_iterator    COOMatTypeConstIt;

    public:
        
        // construct a CSR sparse matrix based on a matrix in 
        // coordinate storage (COO) sparse matrix.
        // This COO sparse matrix is given as a map< pair<int,int>, NekDouble>
        // where the pair refers to the coordinate of the non-zero entry
        // and the NekDouble contains its value. 
        // The constructor now convers from COO storage to CSR storage.
        NekSparseMatrix(unsigned int rows, 
                        unsigned int columns,
                        const COOMatType &cooMat):
            m_size(),
            m_nnz(cooMat.size()),
            m_val(m_nnz),
            m_indx(m_nnz),
            m_pntr(rows+1)
        {
            unsigned int i;
            COOMatTypeConstIt entry;
            int rowcoord;
            int colcoord;
            DataType     value;

            Array<OneD,  int> tmp(rows+1,0);

            m_size[0] = rows;
            m_size[1] = columns;

            // calculate the number of entries on each row
            // and store the result in tmp
            for(entry = cooMat.begin(); entry != cooMat.end(); entry++)
            {
                rowcoord = (entry->first).first;
                tmp[rowcoord]++;
            }
            // Based upon this information, fill the array m_pntr
            // which basically contains the offset of each row's 
            // first entry in the other arrays m_val and m_indx
            m_pntr[0] = 0;
            for(i = 0; i < rows; i++)
            {
                m_pntr[i+1] = m_pntr[i] + tmp[i]; 
            }
            // Copy the values of m_pntr into tmp as this will be needed 
            // in the following step
            std::copy(m_pntr.get(),m_pntr.get()+rows+1,tmp.get());
            // Now, fill in index and value entries.
            for(entry = cooMat.begin(); entry != cooMat.end(); entry++)
            {
                rowcoord = (entry->first).first;
                colcoord = (entry->first).second;
                value    =  entry->second;

                m_val [ tmp[rowcoord] ] = value;
                m_indx[ tmp[rowcoord] ] = colcoord;
                tmp[rowcoord]++;
            }
        }
                  
        unsigned int GetRows() const
        {
            return m_size[0]; 
        }
        
        unsigned int GetColumns() const
        {
            return m_size[1];
        }

        unsigned int GetNumNonZeroEntries() const
        {
            return m_nnz;
        }
        
        const unsigned int* GetSize() const { return m_size; }

        const Array<OneD, const NekDouble>& GetVal() const
        {
            return m_val;
        }

        const Array<OneD, const int>& GetIndx() const
        {
            return m_indx;
        }

        const Array<OneD, const int>& GetPntr() const
        {
            return m_pntr;
        }

        typename boost::call_traits<DataType>::const_reference operator()(unsigned int row, unsigned int column) const
        {                
                return this->GetValue(row, column);
        }
        
        typename boost::call_traits<DataType>::const_reference GetValue(unsigned int row, unsigned int column) const
        {
            ASSERTL1(row < GetRows(), std::string("Row ") + boost::lexical_cast<std::string>(row) + 
                     std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(GetRows()) +
                     std::string(" rows"));
            ASSERTL1(column < GetColumns(), std::string("Column ") + boost::lexical_cast<std::string>(column) + 
                     std::string(" requested in a matrix with a maximum of ") + boost::lexical_cast<std::string>(GetColumns()) +
                     std::string(" columns"));

            unsigned int i;
            static DataType defaultReturnValue;
            for( i = m_pntr[row]; i < m_pntr[row+1]; i++)
            {
                if(column == m_indx[i])
                {
                    return m_val[i];
                }
            }

            return defaultReturnValue;
        }
        
        
    protected:
        
        unsigned int m_size[2];
        unsigned int m_nnz; // number of entries in the sparse matrix
        
        Array<OneD, DataType>     m_val;  // values of non-zero entries
        Array<OneD, int> m_indx; // column indices of non-zero entries
        Array<OneD, int> m_pntr; // m_pntr(i) contains index in m_val of first non-zero element in row i
        
    private:
        
    };

    template<typename DataType>
    std::ostream& operator<<(std::ostream& os, const NekSparseMatrix<DataType>& rhs)
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
        
}

#endif //NEKTAR_LIB_UTILITIES_LINEAR_ALGEBRA_SPARSE_STANDARD_MATRIX_HPP
