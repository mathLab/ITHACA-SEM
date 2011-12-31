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
// Description: Definition of StdMatrixKey
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/functional/hash.hpp>
#include "StdRegions/StdExpansion.h"
#include "StdRegions/StdMatrixKey.h"

namespace Nektar
{
    namespace StdRegions
    {
    
        static int s_matrixcnt = 99;

        StdMatrixKey::StdMatrixKey(const MatrixType matrixType,
                                   const ExpansionType expansionType,
                                   const StdExpansion &stdExpansion,
                                   const ConstFactorMap &factorMap,
                                   const VarCoeffMap &varCoeffMap,
                                   LibUtilities::PointsType nodalType) :
            m_expansionType(expansionType),
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
            for (VarCoeffMap::const_iterator x = varCoeffMap.begin(); x != varCoeffMap.end(); ++x)
            {
                m_varcoeff_hashes[i] = boost::hash_range(x->second.begin(), x->second.end());
                boost::hash_combine(m_varcoeff_hashes[i], x->first);
				i++;
            }
        }


        StdMatrixKey::StdMatrixKey(const StdMatrixKey& rhs) :
            m_expansionType(rhs.m_expansionType),
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
            
            if(lhs.m_ncoeffs < rhs.m_ncoeffs)
            {
                return true;
            }
            
            if(lhs.m_ncoeffs > rhs.m_ncoeffs)
            {
                return false;
            }
            
            for(unsigned int i = 0; i < ExpansionTypeDimMap[lhs.m_expansionType]; ++i)
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
                ConstFactorMap::const_iterator x, y;
                for(x = lhs.m_factors.begin(), y = rhs.m_factors.begin();
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
        
        std::ostream& operator<<(std::ostream& os, const StdMatrixKey& rhs)
        {
            os << "MatrixType: " << MatrixTypeMap[rhs.GetMatrixType()] << ", ShapeType: " 
                << ExpansionTypeMap[rhs.GetExpansionType()] << ", Ncoeffs: " << rhs.GetNcoeffs() 
                << std::endl;

            if(rhs.GetConstFactors().size())
            {
                os << "Constants: " << endl;
                ConstFactorMap::const_iterator x;
                for(x = rhs.GetConstFactors().begin(); x != rhs.GetConstFactors().end(); ++x)
                {
                    os << "\t value " << ConstFactorTypeMap[x->first] <<" : " << x->second << endl;
                }
            }
            if(rhs.GetVarCoeffs().size())
            {
                os << "Variable coefficients: " << endl;
                VarCoeffMap::const_iterator x;
                unsigned int i = 0;
                for (x = rhs.GetVarCoeffs().begin(); x != rhs.GetVarCoeffs().end(); ++x)
                {
                    os << "\t Coeff defined: " << VarCoeffTypeMap[x->first] << endl;
                    os << "\t Hash:          " << rhs.GetVarCoeffHashes()[i++] << endl;
                }
            }
            
            for(unsigned int i = 0; i < ExpansionTypeDimMap[rhs.GetExpansionType()]; ++i)
            {
                os << rhs.GetBase()[i]->GetBasisKey();
            }

            return os;
        }
    }
}

/**
* $Log: StdMatrixKey.cpp,v $
* Revision 1.17  2009/11/07 21:09:11  sehunchun
* Add more functions with various parameters
*
* Revision 1.16  2009/09/06 22:23:18  sherwin
* Somehow deleted opless operator in previous submission
*
* Revision 1.15  2009/09/06 21:55:26  sherwin
* Updates related to Navier Stokes Solver
*
* Revision 1.14  2008/11/19 16:02:47  pvos
* Added functionality for variable Laplacian coeffcients
*
* Revision 1.13  2008/05/30 00:33:49  delisi
* Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
*
* Revision 1.12  2008/04/06 06:04:15  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.11  2007/12/17 13:03:51  sherwin
* Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
*
* Revision 1.10  2007/08/11 23:42:26  sherwin
* A few changes
*
* Revision 1.9  2007/07/26 02:39:21  bnelson
* Fixed Visual C++ compiler errors when compiling in release mode.
*
* Revision 1.8  2007/07/20 02:16:54  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.7  2007/07/09 15:19:15  sherwin
* Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
*
* Revision 1.6  2007/05/15 05:18:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.5  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.4  2007/03/29 19:35:09  bnelson
* Replaced boost::shared_array with SharedArray
*
* Revision 1.3  2007/03/20 16:58:43  sherwin
* Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
*
* Revision 1.2  2007/03/05 08:07:13  sherwin
* Modified so that StdMatrixKey has const calling arguments in its constructor.
*
* Revision 1.1  2007/02/28 19:05:11  sherwin
* Moved key definitions to their own files to make things more transparent
*
* Revision 1.6  2007/02/28 09:53:17  sherwin
* Update including adding GetBasis call to StdExpansion
*
* Revision 1.5  2007/02/24 09:07:25  sherwin
* Working version of stdMatrixManager and stdLinSysMatrix
*
* Revision 1.4  2007/02/23 19:26:04  jfrazier
* General bug fix and formatting.
*
* Revision 1.3  2007/02/22 22:02:27  sherwin
* Update with executing StdMatManager
*
* Revision 1.2  2007/02/22 18:11:31  sherwin
* Version with some create functions introduced for StdMatManagers
*
* Revision 1.1  2007/02/21 22:55:16  sherwin
* First integration of StdMatrixManagers
*
***/

