///////////////////////////////////////////////////////////////////////////////
//
// File GlobalMatrixKey.cpp
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
// Description: Definition of GlobalMatrixKey 
//
///////////////////////////////////////////////////////////////////////////////


#include <MultiRegions/GlobalMatrixKey.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Register Mass Matrix creator. 
        GlobalMatrixKey::GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap):
            m_matrixType(matrixType),
            m_locToGloMap(locToGloMap),
            m_nconstants(0),
            m_constant(m_nconstants),
            m_nvariablecoefficients(0),
            m_variablecoefficient(m_nvariablecoefficients)
        {
        }
        
        GlobalMatrixKey::GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                                         const NekDouble factor,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap):
            m_matrixType(matrixType),
            m_locToGloMap(locToGloMap),
            m_nconstants(1),
            m_constant(m_nconstants),
            m_nvariablecoefficients(0),
            m_variablecoefficient(m_nvariablecoefficients)
        {
            m_constant[0] = factor;
        }
        
        GlobalMatrixKey::GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                                         const NekDouble factor1,
                                         const NekDouble factor2,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap):
            m_matrixType(matrixType),
            m_locToGloMap(locToGloMap),
            m_nconstants(2),
            m_constant(m_nconstants),
            m_nvariablecoefficients(0),
            m_variablecoefficient(m_nvariablecoefficients)
        {
            m_constant[0] = factor1;
            m_constant[1] = factor2;
        }

        GlobalMatrixKey::GlobalMatrixKey(const StdRegions::MatrixType matrixType, 
                                         const Array<OneD,NekDouble>& varcoeffs,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap):
            m_matrixType(matrixType),
            m_locToGloMap(locToGloMap),
            m_nconstants(0),
            m_constant(m_nconstants),
            m_nvariablecoefficients(1)
        {
            m_variablecoefficient[0] = varcoeffs;
        }  

        GlobalMatrixKey::GlobalMatrixKey(const StdRegions::MatrixType matrixType, 
                                         const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap):
            m_matrixType(matrixType),
            m_locToGloMap(locToGloMap),
            m_nconstants(0),
            m_constant(m_nconstants),
            m_nvariablecoefficients(varcoeffs.num_elements()),
            m_variablecoefficient(varcoeffs)
        {
        }          
            
        GlobalMatrixKey::GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                                         const NekDouble factor,
                                         const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap):
            m_matrixType(matrixType),
            m_locToGloMap(locToGloMap),
            m_nconstants(1),
            m_constant(m_nconstants),
            m_nvariablecoefficients(varcoeffs.num_elements()),
            m_variablecoefficient(varcoeffs)
        {
            m_constant[0] = factor;
        }     

        GlobalMatrixKey::GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                                         const NekDouble factor1,
                                         const NekDouble factor2,
                                         const Array<OneD, Array<OneD,NekDouble> >& varcoeffs,
                                         const LocalToGlobalBaseMapSharedPtr &locToGloMap):
            m_matrixType(matrixType),
            m_locToGloMap(locToGloMap),
            m_nconstants(2),
            m_constant(m_nconstants),
            m_nvariablecoefficients(varcoeffs.num_elements()),
            m_variablecoefficient(varcoeffs)
        {
            m_constant[0] = factor1;
            m_constant[1] = factor2;
        }

        GlobalMatrixKey::GlobalMatrixKey(const GlobalMatrixKey &key):
            m_matrixType(key.m_matrixType),
            m_locToGloMap(key.m_locToGloMap),
            m_nconstants(key.m_nconstants),
            m_constant(key.m_constant),
            m_nvariablecoefficients(key.m_nvariablecoefficients),
            m_variablecoefficient(key.m_variablecoefficient)
            
        {
        }

        bool operator<(const GlobalMatrixKey &lhs, const GlobalMatrixKey &rhs)
        {
            if(lhs.m_matrixType < rhs.m_matrixType)
            {
                return true;
            }

            if(lhs.m_matrixType > rhs.m_matrixType)
            {
                return false;
            }

            if(lhs.m_nconstants < rhs.m_nconstants)
            {
                return true;
            }
            else if(lhs.m_nconstants > rhs.m_nconstants)
            {
                return false;
            }
            else 
            {
                for(unsigned int i = 0; i < lhs.m_nconstants; ++i)
                {
                    if(lhs.m_constant[i] < rhs.m_constant[i])
                    {
                        return true;
                    }
                    
                    if(lhs.m_constant[i] > rhs.m_constant[i])
                    {
                        return false;
                    }
                }
            }

            if(lhs.m_nvariablecoefficients < rhs.m_nvariablecoefficients)
            {
                return true;
            }
            else if(lhs.m_nvariablecoefficients > rhs.m_nvariablecoefficients)
            {
                return false;
            }
            else 
            {
                for(unsigned int i = 0; i < lhs.m_nvariablecoefficients; ++i)
                {
                    if((lhs.m_variablecoefficient[i]).get() < (rhs.m_variablecoefficient[i]).get())
                    {
                        return true;
                    }
                    
                    if((lhs.m_variablecoefficient[i]).get() > (rhs.m_variablecoefficient[i]).get())
                    {
                        return false;
                    }
                }
            }


            if(lhs.m_locToGloMap.get() < rhs.m_locToGloMap.get())
            {
                return true;
            }

            return false;
        }

        std::ostream& operator<<(std::ostream& os, const GlobalMatrixKey& rhs)
        {
            int i;
            os << "MatrixType: " << rhs.GetMatrixType() << endl;
            os << "Number of constants: " << rhs.GetNconstants() << endl;
            for(i = 0; i < rhs.GetNconstants();i++) 
            {
                os << "  Constant " << i << ": " << rhs.GetConstant(i) << endl;
            }
            os << "Number of variable coefficients: " << rhs.GetNvariableCoefficients() << endl;
            
            return os;
        }
    }
}

/**
* $Log: GlobalMatrixKey.cpp,v $
* Revision 1.1  2009/03/23 10:46:54  pvos
* Added GlobalMatrixKey
*
**/

