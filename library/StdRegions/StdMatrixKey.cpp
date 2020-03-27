///////////////////////////////////////////////////////////////////////////////
//
// File StdMatrixKey.cpp
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
// Description: Definition of StdMatrixKey
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdMatrixKey.h>
#include <StdRegions/StdExpansion.h>
#include <LibUtilities/BasicUtils/HashUtils.hpp>

using namespace std;

namespace Nektar
{
    namespace StdRegions
    {
        StdMatrixKey::StdMatrixKey(const MatrixType matrixType,
                                   const LibUtilities::ShapeType shapeType,
                                   const StdExpansion &stdExpansion,
                                   const ConstFactorMap &factorMap,
                                   const VarCoeffMap &varCoeffMap,
                                   LibUtilities::PointsType nodalType) :
            m_shapeType(shapeType),
            m_base(stdExpansion.GetBase()),
            m_ncoeffs(stdExpansion.GetNcoeffs()),
            m_matrixType(matrixType),
            m_nodalPointsType(nodalType),
            m_factors(factorMap),
            m_varcoeffs(varCoeffMap),
            m_varcoeff_hashes(varCoeffMap.size())
        {
            // Create hash
            int i = 0;
            for (auto &x : varCoeffMap)
            {
                m_varcoeff_hashes[i] = hash_range(x.second.begin(), x.second.begin() + stdExpansion.GetTotPoints());
                hash_combine(m_varcoeff_hashes[i], (int)x.first);
                i++;
            }
        }

        StdMatrixKey::StdMatrixKey(const StdMatrixKey& rhs,
                      const StdRegions::MatrixType matrixType) :
            m_shapeType(rhs.m_shapeType),
            m_base(rhs.m_base),
            m_ncoeffs(rhs.m_ncoeffs),
            m_matrixType(matrixType),
            m_nodalPointsType(rhs.m_nodalPointsType),
            m_factors(rhs.m_factors),
            m_varcoeffs(rhs.m_varcoeffs),
            m_varcoeff_hashes(rhs.m_varcoeff_hashes)
        {
        }

        StdMatrixKey::StdMatrixKey(const StdMatrixKey& rhs) :
            m_shapeType(rhs.m_shapeType),
            m_base(rhs.m_base),
            m_ncoeffs(rhs.m_ncoeffs),
            m_matrixType(rhs.m_matrixType),
            m_nodalPointsType(rhs.m_nodalPointsType),
            m_factors(rhs.m_factors),
            m_varcoeffs(rhs.m_varcoeffs),
            m_varcoeff_hashes(rhs.m_varcoeff_hashes)
        {
        }
        

        bool StdMatrixKey::opLess::operator()(const StdMatrixKey &lhs, const StdMatrixKey &rhs) const
        {        
            return (lhs.m_matrixType < rhs.m_matrixType);
        }

        bool operator<(const StdMatrixKey &lhs, const StdMatrixKey &rhs)
        {
            if(lhs.m_matrixType < rhs.m_matrixType)
            {
                return true;
            }
            
            if(lhs.m_matrixType > rhs.m_matrixType)
            {
                return false;
            }
            
            if(lhs.m_ncoeffs < rhs.m_ncoeffs) // probably do not need this check since checking the m_base below? 
            {
                return true;
            }
            
            if(lhs.m_ncoeffs > rhs.m_ncoeffs)
            {
                return false;
            }

            if( LibUtilities::ShapeTypeDimMap[lhs.m_shapeType] <
                LibUtilities::ShapeTypeDimMap[rhs.m_shapeType])
            {
                return true;
            }

            if( LibUtilities::ShapeTypeDimMap[lhs.m_shapeType] >
                LibUtilities::ShapeTypeDimMap[rhs.m_shapeType])
            {
                return false;
            }

            for(unsigned int i = 0; i < LibUtilities::ShapeTypeDimMap[lhs.m_shapeType]; ++i)
            {
                if(lhs.m_base[i].get() < rhs.m_base[i].get())
                {
                    return true;
                }
                
                if(lhs.m_base[i].get() > rhs.m_base[i].get())
                {
                    return false;
                }
            }

            if(lhs.m_factors.size() < rhs.m_factors.size())
            {
                return true;
            }
            else if(lhs.m_factors.size() > rhs.m_factors.size())
            {
                return false;
            }
            else
            {
                for(auto x = lhs.m_factors.begin(), y = rhs.m_factors.begin();
                    x != lhs.m_factors.end(); ++x, ++y)
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

            if(lhs.m_varcoeffs.size() < rhs.m_varcoeffs.size())
            {
                return true;
            }

            if(lhs.m_varcoeffs.size() > rhs.m_varcoeffs.size())
            {
                return false;
            }

            for (unsigned int i = 0; i < lhs.m_varcoeff_hashes.size(); ++i)
            {
                if(lhs.m_varcoeff_hashes[i] < rhs.m_varcoeff_hashes[i])
                {
                    return true;
                }
                if(lhs.m_varcoeff_hashes[i] > rhs.m_varcoeff_hashes[i])
                {
                    return false;
                }
            }
            
            if(lhs.m_nodalPointsType < rhs.m_nodalPointsType)
            {
                return true;
            }
            
            if(lhs.m_nodalPointsType > rhs.m_nodalPointsType)
            {
                return false;
            }
            
            return false;
        }

        bool operator==(const StdMatrixKey &lhs, const StdMatrixKey &rhs)
        {
            if(lhs.m_matrixType != rhs.m_matrixType)
            {
                return false;
            }

            if(lhs.m_ncoeffs != rhs.m_ncoeffs)
            {
                return false;
            }

            for(unsigned int i = 0; i < LibUtilities::ShapeTypeDimMap[lhs.m_shapeType]; ++i)
            {
                if(lhs.m_base[i].get() != rhs.m_base[i].get())
                {
                    return false;
                }
            }

            if(lhs.m_factors.size() != rhs.m_factors.size())
            {
                return false;
            }
            else
            {
                return lhs.m_factors.size() == rhs.m_factors.size() &&
                    std::equal(lhs.m_factors.begin(), lhs.m_factors.end(),
                               rhs.m_factors.begin());
            }

            if(lhs.m_nodalPointsType != rhs.m_nodalPointsType)
            {
                return false;
            }

            if(lhs.m_varcoeffs.size() != rhs.m_varcoeffs.size())
            {
                return false;
            }

            for (unsigned int i = 0; i < lhs.m_varcoeff_hashes.size(); ++i)
            {
                if(lhs.m_varcoeff_hashes[i] != rhs.m_varcoeff_hashes[i])
                {
                    return false;
                }
            }

            for (auto &x : lhs.m_varcoeffs)
            {
                auto y = rhs.m_varcoeffs.find(x.first);
                // Check var coeff is found
                if (y == rhs.m_varcoeffs.end())
                {
                    return false;
                }

                if (x.second != y->second)
                {
                    return false;
                }
            }
            for (unsigned int i = 0; i < lhs.m_varcoeffs.size(); ++i)
            {
                if(lhs.m_varcoeff_hashes[i] != rhs.m_varcoeff_hashes[i])
                {
                    return false;
                }
            }

            return true;
        }

        std::ostream& operator<<(std::ostream& os, const StdMatrixKey& rhs)
        {
            os << "MatrixType: " << MatrixTypeMap[rhs.GetMatrixType()] << ", ShapeType: " 
                << LibUtilities::ShapeTypeMap[rhs.GetShapeType()] << ", Ncoeffs: " << rhs.GetNcoeffs() 
                << std::endl;

            if(rhs.GetConstFactors().size())
            {
                os << "Constants: " << endl;
                for(auto &x : rhs.GetConstFactors())
                {
                    os << "\t value " << ConstFactorTypeMap[x.first] <<" : " << x.second << endl;
                }
            }
            if(rhs.GetVarCoeffs().size())
            {
                os << "Variable coefficients: " << endl;
                unsigned int i = 0;
                for (auto &x : rhs.GetVarCoeffs())
                {
                    os << "\t Coeff defined: " << VarCoeffTypeMap[x.first] << endl;
                    os << "\t Hash:          " << rhs.GetVarCoeffHashes()[i++] << endl;
                }
            }
            
            for(unsigned int i = 0; i < LibUtilities::ShapeTypeDimMap[rhs.GetShapeType()]; ++i)
            {
                os << rhs.GetBase()[i]->GetBasisKey();
            }

            return os;
        }
    }
}
