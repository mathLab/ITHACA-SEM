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

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {
        GlobalMatrixKey::GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varCoeffs) :
            m_matrixType(matrixType),
            m_shapeType(LibUtilities::eNoShapeType),
            m_constFactors(factors),
            m_varCoeffs(varCoeffs),
            m_locToGloMap(locToGloMap)
        {
        }

        GlobalMatrixKey::GlobalMatrixKey(const GlobalMatrixKey &key,
                                         const LibUtilities::ShapeType shapeType):
            m_matrixType(key.m_matrixType),
            m_shapeType(shapeType),
            m_constFactors(key.m_constFactors),
            m_varCoeffs(key.m_varCoeffs),
            m_locToGloMap(key.m_locToGloMap)
        {
        }

        GlobalMatrixKey::GlobalMatrixKey(const GlobalMatrixKey &key):
            m_matrixType(key.m_matrixType),
            m_shapeType(key.m_shapeType),
            m_constFactors(key.m_constFactors),
            m_varCoeffs(key.m_varCoeffs),
            m_locToGloMap(key.m_locToGloMap)
        {
        }

        GlobalMatrixKey::~GlobalMatrixKey()
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


            if(lhs.m_shapeType < rhs.m_shapeType)
            {
                return true;
            }
            

            if(lhs.m_shapeType > rhs.m_shapeType)
            {
                return false;
            }
            
            if(lhs.m_constFactors.size() < rhs.m_constFactors.size())
            {
                return true;
            }
            else if(lhs.m_constFactors.size() > rhs.m_constFactors.size())
            {
                return false;
            }
            else
            {
                StdRegions::ConstFactorMap::const_iterator x, y;
                for(x = lhs.m_constFactors.begin(), y = rhs.m_constFactors.begin();
                    x != lhs.m_constFactors.end(); ++x, ++y)
                {
                    if (x->second < y->second)
                    {
                        return true;
                    }
                    if (x->second > y->second)
                    {
                        return false;
                    }
                }
            }

            if(lhs.m_varCoeffs.size() < rhs.m_varCoeffs.size())
            {
                return true;
            }
            else if(lhs.m_varCoeffs.size() > rhs.m_varCoeffs.size())
            {
                return false;
            }
//            else
//            {
//                StdRegions::VarCoeffMap::const_iterator x, y;
//                for (x = lhs.m_varCoeffs.begin(), y = rhs.m_varCoeffs.begin();
//                     x != lhs.m_varCoeffs.end(); ++x, ++y)
//                {
//                    if (x->second.get() < y->second.get())
//                    {
//                        return true;
//                    }
//                    if (x->second.get() > y->second.get())
//                    {
//                        return false;
//                    }
//                }
//            }

            if(!rhs.m_locToGloMap.lock().get())
            {
                return false;
            }
            else if(!lhs.m_locToGloMap.lock().get() && rhs.m_locToGloMap.lock().get() )
            {
                return true;
            }
            if(lhs.m_locToGloMap.lock()->GetHash() < rhs.m_locToGloMap.lock()->GetHash())
            {
                return true;
            }

            return false;
        }

        std::ostream& operator<<(std::ostream& os, const GlobalMatrixKey& rhs)
        {
            os << "MatrixType: " << rhs.GetMatrixType() << endl;
            os << "Number of constants: " << rhs.GetNConstFactors() << endl;
            for(auto &x : rhs.GetConstFactors())
            {
                os << "  Constant " << StdRegions::ConstFactorTypeMap[x.first]
                   << ": " << x.second << endl;
            }
            os << "Number of variable coefficients: " 
               << rhs.GetNVarCoeffs() << endl;

            return os;
        }
    }
}
