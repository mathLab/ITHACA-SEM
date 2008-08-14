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
// Description: Header file of Basis definition 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILIITIES_FOUNDATIONS_BASIS_H
#define NEKTAR_LIB_UTILIITIES_FOUNDATIONS_BASIS_H

#include <math.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Points.h>

namespace Nektar
{
    namespace LibUtilities 
    {
        class BasisKey
        {
        public:
            // Use for looking up the creator. The creator for number of points
            // can generate for any number, so we want the same creator called
            // for all number.
            struct opLess
            {
                bool operator()(const BasisKey &lhs, const BasisKey &rhs) const;
            };


            BasisKey(const BasisType btype, const int nummodes, const PointsKey pkey):
                m_nummodes(nummodes),
                m_basistype(btype), 
                m_pointsKey(pkey)
            {
            }

            BasisKey(const BasisKey &B): 
                m_nummodes(B.m_nummodes),
                m_basistype(B.m_basistype),
                m_pointsKey(B.m_pointsKey)
            {
            }

            ~BasisKey()
            {
            }

            /** \brief return order of basis */
            inline int GetNumModes() const
            {
                return m_nummodes;
            }

            inline int GetTotNumModes() const
            {
                int value = 0;

                switch(m_basistype)
                {
                case eOrtho_B:
                case eModified_B:
                case eModified_C:
                    value = m_nummodes*(m_nummodes+1)/2;
                    break;

                case eOrtho_C:
                    value = m_nummodes*(m_nummodes+1)*(m_nummodes+2)/6;
                    break;

                case eOrtho_A:      
                case eModified_A:   
                case eFourier:      
                case eGLL_Lagrange: 
                case eLegendre:    
                case eChebyshev:
                case eMonomial:
                    value = m_nummodes;
                    break;

                default:
                    NEKERROR(ErrorUtil::efatal,"Unknown basis being used");
                }
                return value;
            }


            /** \brief return points order at which  basis is defined */
            inline int GetNumPoints() const
            {
                return m_pointsKey.GetNumPoints();
            }

            inline int GetTotNumPoints() const
            {
                return m_pointsKey.GetTotNumPoints();
            }

            /** \brief return type of expansion basis */
            inline BasisType GetBasisType() const
            {
                return m_basistype;
            }

            inline PointsKey GetPointsKey() const
            {
                return m_pointsKey;
            }

            /** \brief return type of quadrature */
            inline PointsType GetPointsType() const
            {
                return m_pointsKey.GetPointsType();
            }    

            /** \brief Check to see if the quadrature of expansions x is the same as the calling basis */

            inline bool SamePoints(const BasisKey &x) const
            {
                return (x.m_pointsKey == m_pointsKey);
            }

            /** \brief Check to see if basis expansions x is the same as the calling basis */
            inline bool SameExp(const BasisKey &x) const
            {
                return ((x.m_nummodes == m_nummodes)&&(x.m_basistype == m_basistype));
            }


            /** \brief determine if basis definition has exact integration for
            *  inner product
            */
            bool ExactIprodInt() const;

            /** \brief Determine if basis has collocation properties,
            *  i.e. GLL_Lagrange with appropriate quadrature
            */
            bool  Collocation() const;

            //Overloaded Operators
            friend bool operator  == (const BasisKey& x, const BasisKey& y);
            friend bool operator  == (const BasisKey* x, const BasisKey& y);
            friend bool operator  == (const BasisKey& x, const BasisKey *y);
            friend bool operator  != (const BasisKey& x, const BasisKey& y);
            friend bool operator  != (const BasisKey* x, const BasisKey& y);
            friend bool operator  != (const BasisKey& x, const BasisKey *y);

            friend bool operator<(const BasisKey &lhs, const BasisKey &rhs);
            friend bool opLess::operator()(const BasisKey &lhs, const BasisKey &rhs) const;

        protected:
            int        m_nummodes;   /**< Expansion Order */
            BasisType  m_basistype;  /**< Expansion Type */
            PointsKey  m_pointsKey;

        private:
        BasisKey():m_pointsKey(NullPointsKey)
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor BasisKey should never be called");
            }

        };

        static const BasisKey NullBasisKey(eNoBasisType, 0, NullPointsKey);

        /////////////////////////////////////////////////////////////////////////
        class Basis
        {
        public:

            static boost::shared_ptr<Basis> Create(const BasisKey &bkey);

            // default destructor()
            virtual ~Basis()
            {
            };

            /** \brief return order of basis */
            inline int GetNumModes() const
            {
                return m_basisKey.GetNumModes();
            }

            inline int GetTotNumModes() const
            {
                return m_basisKey.GetTotNumModes();
            }

            /** \brief return points order at which  basis is defined */
            inline int GetNumPoints() const
            {
                return m_basisKey.GetNumPoints();
            }

            /** \brief return points order at which  basis is defined */
            inline int GetTotNumPoints() const
            {
                return m_basisKey.GetTotNumPoints();
            }

            /** \brief return type of expansion basis */
            inline BasisType GetBasisType() const
            {
                return m_basisKey.GetBasisType();
            }

            inline PointsKey GetPointsKey() const
            {
                return m_basisKey.GetPointsKey();
            }

            /** \brief return type of quadrature */
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

            inline const  boost::shared_ptr<NekMatrix<NekDouble> > GetD(Direction dir = xDir) const
            {
                return m_points->GetD(dir);
            }

            const boost::shared_ptr<NekMatrix<NekDouble> > GetI(const Array<OneD, const NekDouble>& x)
            {
                return m_points->GetI(x);
            }
            
            /** \brief determine if basis definition has exact integration for
            *  inner product
            */
            inline bool ExactIprodInt() const
            {
                return m_basisKey.ExactIprodInt();
            }

            /** \brief Determine if basis has collocation properties,
            *  i.e. GLL_Lagrange with appropriate quadrature
            */
            inline bool  Collocation() const
            {
                return m_basisKey.Collocation();
            }

            /** \brief return basis definition array m_bdata */
            inline const Array<OneD, const NekDouble>& GetBdata() const 
            {
                return m_bdata;
            }

            /** \brief return basis definition array m_dbdata */

            inline const Array<OneD, const NekDouble>& GetDbdata() const
            {
                return m_dbdata;
            }

            inline const BasisKey GetBasisKey() const
            {
                return m_basisKey;
            }

            virtual void Initialize();

        protected:
            BasisKey        m_basisKey;
            PointsSharedPtr m_points;
            Array<OneD, NekDouble> m_bdata; /**< Basis definition */
            Array<OneD, NekDouble> m_dbdata; /**< Derivative Basis definition */

        private:

            Basis(const BasisKey &bkey);

            Basis():m_basisKey(NullBasisKey)
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for Basis should not be called");
            }

            void GenBasis();

        };

        bool operator<(const BasisKey &lhs, const BasisKey &rhs);
        bool operator>(const BasisKey &lhs, const BasisKey &rhs);

        std::ostream& operator<<(std::ostream& os, const BasisKey& rhs);

        typedef boost::shared_ptr<Basis> BasisSharedPtr;
        typedef std::vector< BasisSharedPtr > BasisVector; 
        typedef std::vector< BasisSharedPtr >::iterator BasisVectorIter; 
        
        static BasisSharedPtr NullBasisSharedPtr;

    } // end of namespace
} // end of namespace

#endif //NEKTAR_LIB_UTILIITIES_FOUNDATIONS_BASIS_H


