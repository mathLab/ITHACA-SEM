///////////////////////////////////////////////////////////////////////////////
//
// File Basis.h
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
// Description: Header file of Basis definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILIITIES_FOUNDATIONS_BASIS_H
#define NEKTAR_LIB_UTILIITIES_FOUNDATIONS_BASIS_H

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        /// Describes the specification for a Basis.
        class BasisKey
        {
        public:
            // Use for looking up the creator. The creator for number of points
            // can generate for any number, so we want the same creator called
            // for all number.
            struct opLess
            {
                LIB_UTILITIES_EXPORT  bool operator()(const BasisKey &lhs, const BasisKey &rhs) const;
            };

            /// Constructor
            BasisKey(const BasisType btype, const int nummodes,
                     const PointsKey pkey):
                m_nummodes(nummodes),
                m_basistype(btype),
                m_pointsKey(pkey)
            {
            }

            /// Copy constructor
            BasisKey(const BasisKey &B):
                m_nummodes(B.m_nummodes),
                m_basistype(B.m_basistype),
                m_pointsKey(B.m_pointsKey)
            {
            }

            /// Destructor
            ~BasisKey()
            {
            }

            /// Assignment operator
            BasisKey& operator=(const BasisKey &) = default;

            /// Returns the order of the basis.
            inline int GetNumModes() const
            {
                return m_nummodes;
            }

            /// \todo Generalise to arbitrary polynomials
            inline int GetTotNumModes() const
            {
                int value = 0;

                switch(m_basistype)
                {
                case eOrtho_B:
                case eModified_B:
                    value = m_nummodes*(m_nummodes+1)/2;
                    break;

                case eModified_C:
                case eOrtho_C:
                    value = m_nummodes*(m_nummodes+1)*(m_nummodes+2)/6;
                    break;
                case eModifiedPyr_C:
                case eOrthoPyr_C:
                    value = m_nummodes*(m_nummodes+1)*(2*m_nummodes+1)/6;
                    break;
                case eOrtho_A:
                case eModified_A:
                case eFourier:
                case eGLL_Lagrange:
                case eGauss_Lagrange:
                case eLegendre:
                case eChebyshev:
                case eMonomial:
                case eFourierSingleMode:
                case eFourierHalfModeRe:
                case eFourierHalfModeIm:
                    value = m_nummodes;
                    break;
                    
                default:
                    NEKERROR(ErrorUtil::efatal,"Unknown basis being used");
                }
                return value;
            }


            /// Return points order at which  basis is defined
            inline int GetNumPoints() const
            {
                return m_pointsKey.GetNumPoints();
            }

            inline int GetTotNumPoints() const
            {
                return m_pointsKey.GetTotNumPoints();
            }

            /// Return type of expansion basis.
            inline BasisType GetBasisType() const
            {
                return m_basistype;
            }

            /// Return distribution of points.
            inline PointsKey GetPointsKey() const
            {
                return m_pointsKey;
            }

            /// Return type of quadrature.
            inline PointsType GetPointsType() const
            {
                return m_pointsKey.GetPointsType();
            }

            /// Determine if quadrature of expansion \a x matches this.
            inline bool SamePoints(const BasisKey &x) const
            {
                return (x.m_pointsKey == m_pointsKey);
            }

            /// Determine if basis expansion \a x matches this.
            inline bool SameExp(const BasisKey &x) const
            {
                return ((x.m_nummodes == m_nummodes)
                                &&(x.m_basistype == m_basistype));
            }

            /// Determine if basis has exact integration for inner product.
            LIB_UTILITIES_EXPORT bool ExactIprodInt() const;

            /// Determine if basis has collocation properties.
            LIB_UTILITIES_EXPORT  bool  Collocation() const;

            /// Overloaded Operators
            LIB_UTILITIES_EXPORT friend bool operator  == (const BasisKey& x, const BasisKey& y);
            LIB_UTILITIES_EXPORT friend bool operator  == (const BasisKey* x, const BasisKey& y);
            LIB_UTILITIES_EXPORT friend bool operator  == (const BasisKey& x, const BasisKey *y);
            LIB_UTILITIES_EXPORT friend bool operator  != (const BasisKey& x, const BasisKey& y);
            LIB_UTILITIES_EXPORT friend bool operator  != (const BasisKey* x, const BasisKey& y);
            LIB_UTILITIES_EXPORT friend bool operator  != (const BasisKey& x, const BasisKey *y);

            LIB_UTILITIES_EXPORT friend bool operator<(const BasisKey &lhs, const BasisKey &rhs);
            LIB_UTILITIES_EXPORT friend bool opLess::operator()( const BasisKey &lhs,
                                            const BasisKey &rhs) const;

        protected:
            int        m_nummodes;   ///< Expansion order.
            BasisType  m_basistype;  ///< Expansion type.
            PointsKey  m_pointsKey;  ///< Points specification.

        private:
            BasisKey():m_pointsKey(NullPointsKey)
            {
                NEKERROR(ErrorUtil::efatal,
                         "Default Constructor BasisKey should never be called");
            }

        };

        /// Defines a null basis with no type or points.
        static const BasisKey NullBasisKey(eNoBasisType, 0, NullPointsKey);


        /// Represents a basis of a given type.
        class Basis
        {
        public:
            /// Returns a new instance of a Basis with given BasisKey.
            static std::shared_ptr<Basis> Create(const BasisKey &bkey);

            /// Destructor.
            virtual ~Basis()
            {
            };

            /// Return order of basis from the basis specification.
            inline int GetNumModes() const
            {
                return m_basisKey.GetNumModes();
            }

            /// Return total number of modes from the basis specification.
            inline int GetTotNumModes() const
            {
                return m_basisKey.GetTotNumModes();
            }

            /// Return the number of points from the basis specification.
            inline int GetNumPoints() const
            {
                return m_basisKey.GetNumPoints();
            }

            /// Return total number of points from the basis specification.
            inline int GetTotNumPoints() const
            {
                return m_basisKey.GetTotNumPoints();
            }

            /// Return the type of expansion basis.
            inline BasisType GetBasisType() const
            {
                return m_basisKey.GetBasisType();
            }

            /// Return the points specification for the basis.
            inline PointsKey GetPointsKey() const
            {
                return m_basisKey.GetPointsKey();
            }

            /// Return the type of quadrature.
            inline PointsType GetPointsType() const
            {
                return m_basisKey.GetPointsType();
            }

            inline const Array<OneD, const NekDouble>& GetZ() const
            {
                return m_points->GetZ();
            }

            inline const Array<OneD, const NekDouble>& GetW() const
            {
                return m_points->GetW();
            }

            inline void GetZW(Array<OneD, const NekDouble> &z,
                              Array<OneD, const NekDouble> &w) const
            {
                m_points->GetZW(z,w);
            }

            inline const std::shared_ptr<NekMatrix<NekDouble> > & GetD(
                              Direction dir = xDir) const
            {
                return m_points->GetD(dir);
            }

            const std::shared_ptr<NekMatrix<NekDouble> > GetI(
                                const Array<OneD, const NekDouble>& x)
            {
                return m_points->GetI(x);
            }

            const std::shared_ptr<NekMatrix<NekDouble> > GetI(
                                const BasisKey &bkey)
            {
                ASSERTL0(bkey.GetPointsKey().GetPointsDim()==1,
                         "Interpolation only to other 1d basis");
                return m_InterpManager[bkey];
            }

            /// Determine if basis has exact integration for inner product.
            inline bool ExactIprodInt() const
            {
                return m_basisKey.ExactIprodInt();
            }

            /// Determine if basis has collocation properties.
            inline bool  Collocation() const
            {
                return m_basisKey.Collocation();
            }

            /// Return basis definition array m_bdata.
            inline const Array<OneD, const NekDouble>& GetBdata() const
            {
                return m_bdata;
            }

            /// Return basis definition array m_dbdata.
            inline const Array<OneD, const NekDouble>& GetDbdata() const
            {
                return m_dbdata;
            }

            /// Returns the specification for the Basis.
            inline const BasisKey GetBasisKey() const
            {
                return m_basisKey;
            }

            LIB_UTILITIES_EXPORT virtual void Initialize();

        protected:
            BasisKey                m_basisKey; ///< Basis specification.
            PointsSharedPtr         m_points;   ///< Set of points.
            Array<OneD, NekDouble>  m_bdata;    ///< Basis definition.
            Array<OneD, NekDouble>  m_dbdata;   ///< Derivative Basis definition.
            NekManager<BasisKey, NekMatrix<NekDouble>, BasisKey::opLess>
                                    m_InterpManager;

        private:
            static bool initBasisManager[];

            /// Private constructor with BasisKey.
            Basis(const BasisKey &bkey);

            /// Private default constructor.
            Basis():m_basisKey(NullBasisKey)
            {
                NEKERROR(ErrorUtil::efatal,
                         "Default Constructor for Basis should not be called");
            }

            std::shared_ptr< NekMatrix<NekDouble> > CalculateInterpMatrix(
                                const BasisKey &tbasis0);

            /// Generate appropriate basis and their derivatives.
            void GenBasis();
        };

        LIB_UTILITIES_EXPORT bool operator<(const BasisKey &lhs, const BasisKey &rhs);
        LIB_UTILITIES_EXPORT bool operator>(const BasisKey &lhs, const BasisKey &rhs);

        LIB_UTILITIES_EXPORT std::ostream& operator<<(std::ostream& os, const BasisKey& rhs);

        static BasisSharedPtr NullBasisSharedPtr;
        static Array<OneD, BasisSharedPtr> NullBasisSharedPtr1DArray;
        
    } // end of namespace LibUtilities
} // end of namespace Nektar

#endif //NEKTAR_LIB_UTILIITIES_FOUNDATIONS_BASIS_H


