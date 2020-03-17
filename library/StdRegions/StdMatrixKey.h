///////////////////////////////////////////////////////////////////////////////
//
// File StdMatrixKeys.h
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
// Description: Headers for StdMatrixKey
//
///////////////////////////////////////////////////////////////////////////////

#ifndef STDMATRIXKEY_H
#define STDMATRIXKEY_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdRegionsDeclspec.h>
#include <LibUtilities/Foundations/Foundations.hpp>  // for PointsType, etc
#include <LibUtilities/Foundations/Basis.h>


namespace Nektar
{
    namespace StdRegions
    {

        class StdExpansion;

        class StdMatrixKey
        {
        public:
            STD_REGIONS_EXPORT StdMatrixKey( const StdRegions::MatrixType matrixType,
                          const LibUtilities::ShapeType shapeType,
                          const StdRegions::StdExpansion &stdExpansion,
                          const ConstFactorMap &factorMap = NullConstFactorMap,
                          const VarCoeffMap &varCoeffMap = NullVarCoeffMap,
                          LibUtilities::PointsType nodalType = LibUtilities::eNoPointsType);

            STD_REGIONS_EXPORT StdMatrixKey(const StdMatrixKey& rhs,
                          const StdRegions::MatrixType matrixType);

            STD_REGIONS_EXPORT StdMatrixKey(const StdMatrixKey& rhs);

            virtual ~StdMatrixKey()
            {
            }

            /// Used to lookup the create function in NekManager.
            struct opLess
            {
                STD_REGIONS_EXPORT bool operator()(const StdMatrixKey &lhs, const StdMatrixKey &rhs) const;
            };

            /// Used for finding value given the key in NekManager.
            STD_REGIONS_EXPORT friend bool operator<(const StdMatrixKey &lhs, const StdMatrixKey &rhs);
            STD_REGIONS_EXPORT friend bool operator==(const StdMatrixKey &lhs, const StdMatrixKey &rhs);
            STD_REGIONS_EXPORT friend bool opLess::operator()(const StdMatrixKey &lhs, const StdMatrixKey &rhs) const;

            MatrixType GetMatrixType() const
            {
                return m_matrixType;
            }

            LibUtilities::ShapeType GetShapeType() const
            {
                return m_shapeType;
            }

            LibUtilities::PointsType GetNodalPointsType() const
            {
                return m_nodalPointsType;
            }
           
            int GetNcoeffs() const
            {
                return m_ncoeffs;
            }

            inline const Array<OneD, const LibUtilities::BasisSharedPtr>& GetBase() const
            {
                return m_base;
            }

            std::vector<std::size_t> GetVarCoeffHashes() const
            {
                return m_varcoeff_hashes;
            }

            inline const LibUtilities::BasisSharedPtr GetBasis(int dir) const
            {
                return(m_base[dir]);
            }

            inline int GetNConstFactors() const
            {
                return m_factors.size();
            }

            inline NekDouble GetConstFactor(const ConstFactorType& factor) const
            {
                auto x = m_factors.find(factor);
                ASSERTL1(x != m_factors.end(),
                        "Constant factor not defined: "
                        + std::string(StdRegions::ConstFactorTypeMap[factor]));
                return x->second;
            }

            inline  bool ConstFactorExists(const ConstFactorType& factor) const
            {
                return m_factors.find(factor) != m_factors.end();
            }

            inline const ConstFactorMap& GetConstFactors() const
            {
                return m_factors;
            }

            inline int GetNVarCoeff() const
            {
                return m_varcoeffs.size();
            }

            inline const Array<OneD, const NekDouble> &GetVarCoeff(const StdRegions::VarCoeffType & coeff) const
            {
                auto x = m_varcoeffs.find(coeff);
                ASSERTL1(x != m_varcoeffs.end(),
                        "Variable coefficient not defined: "
                        + std::string(StdRegions::VarCoeffTypeMap[coeff]));
                return x->second;
            }

            inline const VarCoeffMap GetVarCoeffAsMap(const VarCoeffType & coeff) const
            {
                VarCoeffMap m;
                m[coeff] = GetVarCoeff(coeff);
                return m;
            }

            inline const VarCoeffMap& GetVarCoeffs() const
            {
                return m_varcoeffs;
            }

            inline bool HasVarCoeff(const StdRegions::VarCoeffType & coeff) const
            {
                return (m_varcoeffs.find(coeff) != m_varcoeffs.end());
            }

        protected:
            LibUtilities::ShapeType m_shapeType;
            Array<OneD, const LibUtilities::BasisSharedPtr> m_base;

            unsigned int m_ncoeffs;
            MatrixType   m_matrixType;
            LibUtilities::PointsType m_nodalPointsType;
            
            ConstFactorMap m_factors;
            VarCoeffMap m_varcoeffs;

            std::vector<std::size_t>     m_varcoeff_hashes;
        private:
            StdMatrixKey();
        };

        STD_REGIONS_EXPORT std::ostream& operator<<(std::ostream& os, const StdMatrixKey& rhs);

        typedef  std::shared_ptr<StdMatrixKey> StdMatrixKeySharedPtr;

    } // end of namespace
} // end of namespace

#endif
