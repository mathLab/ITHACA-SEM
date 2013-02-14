///////////////////////////////////////////////////////////////////////////////
//
// File GlobalMatrixKey.h
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
// Description: Headers for GlobalMatrixKey
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_GLOBALMATRIXKEY_H
#define NEKTAR_LIBS_MULTIREGIONS_GLOBALMATRIXKEY_H

#include <MultiRegions/AssemblyMap/AssemblyMap.h>
#include <StdRegions/StdRegions.hpp>

namespace Nektar
{
    namespace MultiRegions
    {

        /// Describes a matrix with ordering defined by a local to global map.
        class GlobalMatrixKey
        {
        public:
            MULTI_REGIONS_EXPORT GlobalMatrixKey(const StdRegions::MatrixType matrixType,
                            const AssemblyMapSharedPtr &locToGloMap
                                    = NullAssemblyMapSharedPtr,
                            const StdRegions::ConstFactorMap &factors = StdRegions::NullConstFactorMap,
                            const StdRegions::VarCoeffMap &varCoeffs = StdRegions::NullVarCoeffMap);

            /// Copy constructor with change in expansion type
            GlobalMatrixKey(const GlobalMatrixKey &key,
                            const StdRegions::ExpansionType expType);
            /// Copy constructor.
            MULTI_REGIONS_EXPORT GlobalMatrixKey(const GlobalMatrixKey &key);

            /// Destructor
            MULTI_REGIONS_EXPORT virtual ~GlobalMatrixKey();

            /// Provides ordering of GlobalMatrixKey objects.
            MULTI_REGIONS_EXPORT friend bool operator<(const GlobalMatrixKey &lhs,
                                  const GlobalMatrixKey &rhs);

            /// Return the matrix type.
            MULTI_REGIONS_EXPORT StdRegions::MatrixType GetMatrixType() const;
            /// Return the expansion type associated with key
            MULTI_REGIONS_EXPORT StdRegions::ExpansionType GetExpansionType()  const;
            /// Returns true if a local to global map is defined.
            MULTI_REGIONS_EXPORT bool LocToGloMapIsDefined() const;
            /// Returns the number of constants defined for this matrix.
            MULTI_REGIONS_EXPORT int GetNConstFactors() const;
            /// Returns the requested constant.
            MULTI_REGIONS_EXPORT NekDouble GetConstFactor(const StdRegions::ConstFactorType & factor) const;
            /// Returns all the constants.
            MULTI_REGIONS_EXPORT const StdRegions::ConstFactorMap& GetConstFactors() const;

            MULTI_REGIONS_EXPORT int GetNVarCoeffs() const;

            MULTI_REGIONS_EXPORT const Array<OneD, const NekDouble> & GetVarCoeff(const StdRegions::VarCoeffType& coeff) const;
            MULTI_REGIONS_EXPORT const StdRegions::VarCoeffMap & GetVarCoeffs() const;

        protected:
            /// Default constructor.
            GlobalMatrixKey();

            /// Stores the matrix type based on the enum StdRegions::MatrixType.
            StdRegions::MatrixType m_matrixType;
            
            /// Stores the expansion/shape type that the matrix is to
            /// be based on
            StdRegions::ExpansionType m_expansionType;

            StdRegions::ConstFactorMap  m_constFactors;
            StdRegions::VarCoeffMap     m_varCoeffs;

            /// Pointer to the local to global mapping.
            AssemblyMapSharedPtr m_locToGloMap;

        private:

        };

        /// Writes statistics about the matrix key to an output stream.
        MULTI_REGIONS_EXPORT std::ostream& operator<<(std::ostream& os, const GlobalMatrixKey& rhs);

        /// A pointer to a GlobalMatrixKey.
        typedef  boost::shared_ptr<GlobalMatrixKey> GlobalMatrixKeySharedPtr;

        inline StdRegions::MatrixType
                        GlobalMatrixKey::GetMatrixType() const
        {
            return m_matrixType;
        }

        inline StdRegions::ExpansionType
                        GlobalMatrixKey::GetExpansionType() const
        {
            return m_expansionType;
        }

        inline bool GlobalMatrixKey::LocToGloMapIsDefined(void) const
        {
            if( m_locToGloMap.get() == 0) //NullAssemblyMapSharedPtr)
            {
                return false;
            }

            return true;
        }

        inline int GlobalMatrixKey::GetNConstFactors() const
        {
            return m_constFactors.size();
        }

        /// @Todo error checking
        inline NekDouble GlobalMatrixKey::GetConstFactor(const StdRegions::ConstFactorType &factor) const
        {
            StdRegions::ConstFactorMap::const_iterator found = m_constFactors.find(factor);
            return (*found).second;
        }

        inline const StdRegions::ConstFactorMap&
                        GlobalMatrixKey::GetConstFactors() const
        {
            return m_constFactors;
        }

        inline int GlobalMatrixKey::GetNVarCoeffs() const
        {
            return m_varCoeffs.size();
        }

        inline const Array<OneD, const NekDouble> & GlobalMatrixKey::GetVarCoeff(const StdRegions::VarCoeffType &coeff) const
        {
            StdRegions::VarCoeffMap::const_iterator found = m_varCoeffs.find(coeff);
            return (*found).second;
        }

        inline const StdRegions::VarCoeffMap & GlobalMatrixKey::GetVarCoeffs() const
        {
            return m_varCoeffs;
        }
    }
}

#endif
