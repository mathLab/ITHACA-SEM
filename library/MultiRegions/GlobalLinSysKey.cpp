///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysKey.cpp
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
// Description: Definition of GlobalLinSysKey 
//
///////////////////////////////////////////////////////////////////////////////


#include <MultiRegions/GlobalLinSysKey.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysKey
         * 
         * This class represents a global linear system and is in essence a
         * wrapper around a global matrix key, augmenting it with a specific
         * solution type from GlobalSysSolnType. Each constructor accepts a
         * MatrixType, describing the matrix to be constructed, a
         * AssemblyMap, defining the mapping from the local elemental
         * expansions to a global system, and a GlobalSysSolnType, defining the
         * type of solution (e.g. full matrix, static condenstation). Some
         * constructors include additional parameters for customising the 
         * global operator matrix.
         */

        GlobalLinSysKey::GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const StdRegions::ConstFactorMap &factors,
                    const StdRegions::VarCoeffMap &varCoeffs) :
            GlobalMatrixKey(matrixType, locToGloMap, factors, varCoeffs),
            m_solnType(locToGloMap->GetGlobalSysSolnType())
        {

        }


        /**
         * @param   key         Existing key to duplicate.
         */
        GlobalLinSysKey::GlobalLinSysKey(const GlobalLinSysKey &key):
            GlobalMatrixKey(key),
            m_solnType(key.m_solnType)
        {
        }


        /**
         *
         */
        GlobalLinSysKey::~GlobalLinSysKey()
        {
        }
        
        
        /**
         * Compares two GlobalLinSysKeys by comparing their solution types and
         * matrix keys.
         * @param   lhs         First operand.
         * @param   rhs         Second operand.
         * @returns true if the first operand is considered less than the
         *          second operand.
         */
        bool operator<(const GlobalLinSysKey &lhs, const GlobalLinSysKey &rhs)
        {
            if(lhs.m_solnType < rhs.m_solnType)
            {
                return true;
            }

            if(lhs.m_solnType > rhs.m_solnType)
            {
                return false;
            }

            return (*dynamic_cast<const GlobalMatrixKey*>(&lhs)
                    < *dynamic_cast<const GlobalMatrixKey*>(&rhs));
        }


        /**
         * Writes the vital statistics of a global linear system to a stream.
         * @param   os          Output stream.
         * @param   rhs         GlobalLinSys object to use.
         * @returns Reference to the output stream \a os.
         */
        std::ostream& operator<<(std::ostream& os, const GlobalLinSysKey& rhs)
        {
            os << "MatrixType: " << StdRegions::MatrixTypeMap[rhs.GetMatrixType()] << ", ShapeType: "
                            << StdRegions::ExpansionTypeMap[rhs.GetExpansionType()]
                            << std::endl;
            os << "Solution Type: " 
               << GlobalSysSolnTypeMap[rhs.GetGlobalSysSolnType()] << endl;
            os << "Number of constants: " << rhs.GetNConstFactors() << endl;
            StdRegions::ConstFactorMap::const_iterator x;
            for (x = rhs.GetConstFactors().begin(); x != rhs.GetConstFactors().end(); ++x)
            {
                os << "  Constant " << StdRegions::ConstFactorTypeMap[x->first]
                   << ": " << x->second << endl;
            }
            os << "Number of variable coefficients: " 
               << rhs.GetNVarCoeffs() << endl;
            
            return os;
        }
    }
}
