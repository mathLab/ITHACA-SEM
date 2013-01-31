///////////////////////////////////////////////////////////////////////////////
//
// File: Points.hpp
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
// Description: Header file of Points definition 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_FOUNDATIONS_POINTS_H
#define NEKTAR_LIB_UTILITIES_FOUNDATIONS_POINTS_H

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>

//#include <LibUtilities/BasicUtils/BasicUtilsFwd.hpp>  // for NekManager
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>


namespace Nektar
{

    namespace LibUtilities
    {
        // Need to add method to compute total number of points given dimension
        // and number of points.

        /// Defines a specification for a set of points.
        class PointsKey
        {
        public:
            // Used for looking up the creator.  The creator for number of points
            // can generate for any number, so we want the same creator called
            // for all number.
            struct opLess
            {
                LIB_UTILITIES_EXPORT bool operator()(const PointsKey &lhs, const PointsKey &rhs) const;
            };
            
            /// Default constructor.
            PointsKey(void):
                m_numpoints(0), 
                m_pointstype(eNoPointsType)
            {
            }

            /// Constructor defining the number and distribution of points.
            PointsKey(const int &numpoints, const PointsType &pointstype): 
                m_numpoints(numpoints), 
                m_pointstype(pointstype)
            {
            }

            /// Destructor.
            virtual ~PointsKey()
            {
            }

            /// Copy constructor.
            PointsKey(const PointsKey &key)
            {
                *this = key; // defer to assignment operator
            }

            PointsKey& operator=(const PointsKey &key)
            {
                m_numpoints = key.m_numpoints;
                m_pointstype  = key.m_pointstype;

                return *this;
            }

            inline unsigned int GetNumPoints() const
            {
                return m_numpoints;
            }

            inline PointsType GetPointsType() const
            {
                return m_pointstype;
            }

            inline bool operator==(const PointsKey &key)
            {
                return (m_numpoints == key.m_numpoints &&
                    m_pointstype == key.m_pointstype);
            }

            inline bool operator== (const PointsKey *y)
            {
                return (*this == *y);
            }

            inline bool operator != (const PointsKey& y)
            {
                return (!(*this == y));
            }

            inline bool operator != (const PointsKey *y)
            {
                return (!(*this == *y));
            }

            // If new points are added, this function must be modified
            inline unsigned int GetPointsDim() const
            {
                int dimpoints = 1;

                switch(m_pointstype)
                {
                case eNodalTriElec:
                case eNodalTriFekete:
                case eNodalTriEvenlySpaced:
                    dimpoints = 2;
                    break;

                case eNodalTetElec:
                case eNodalTetEvenlySpaced:
                case eNodalPrismEvenlySpaced:
                    dimpoints = 3;
                    break;
                    
                default:
                    break;
                }

                return dimpoints;
            }

            // If new points are added, this function must be modified
            inline unsigned int GetTotNumPoints() const
            {
                int totpoints = m_numpoints;

                switch(m_pointstype)
                {
                case eNodalTriElec:
                case eNodalTriFekete:
                case eNodalTriEvenlySpaced:
                    totpoints = m_numpoints*(m_numpoints+1)/2;
                    break;

                case eNodalTetElec:
                case eNodalTetEvenlySpaced:
                    totpoints = m_numpoints*(m_numpoints+1)*(m_numpoints+2)/6;
                    break;
                    
                case eNodalPrismEvenlySpaced:
                    totpoints = m_numpoints*m_numpoints*(m_numpoints+1)/2;
                    break;
                    
                default:
                    break;
                }

                return totpoints;
            }

            LIB_UTILITIES_EXPORT friend bool operator==(const PointsKey &lhs, const PointsKey &rhs);
            LIB_UTILITIES_EXPORT friend bool operator<(const PointsKey &lhs, const PointsKey &rhs);
            LIB_UTILITIES_EXPORT friend bool opLess::operator()(const PointsKey &lhs, const PointsKey &rhs) const;

        protected:
            unsigned int m_numpoints;     //!< number of the points (as appropriately defined for PointsType)
            PointsType m_pointstype;      //!< Type of Points

        private:
        };

        static const PointsKey NullPointsKey(0, eNoPointsType);

        LIB_UTILITIES_EXPORT bool operator==(const PointsKey &lhs, const PointsKey &rhs);
        LIB_UTILITIES_EXPORT bool operator<(const PointsKey &lhs, const PointsKey &rhs);
        LIB_UTILITIES_EXPORT std::ostream& operator<<(std::ostream& os, const PointsKey& rhs);

        /// Stores a set of points of datatype DataT, defined by a PointKey.
        template<typename DataT>
        class Points
        {
        public:
            typedef DataT DataType;
            typedef boost::shared_ptr<NekMatrix<DataType> > MatrixSharedPtrType;

            virtual ~Points()
            {
            }

            virtual void Initialize(void)
            {
                CalculatePoints();
                CalculateWeights();
                CalculateDerivMatrix();
            }

            inline unsigned int GetPointsDim() const
            {
                return m_pointsKey.GetPointsDim();
            }

            inline unsigned int GetNumPoints() const
            {
                return m_pointsKey.GetNumPoints();
            }

            inline unsigned int GetTotNumPoints() const
            {
                return m_pointsKey.GetTotNumPoints();
            }

            inline PointsType GetPointsType() const
            {
                return m_pointsKey.GetPointsType();
            }

            inline const Array<OneD, const DataType>& GetZ() const
            {
                return m_points[0];
            }

            inline const Array<OneD, const DataType>& GetW() const 
            {
                return m_weights; 
            } 

            inline void GetZW(Array<OneD, const DataType> &z,
                Array<OneD, const DataType> &w) const 
            {
                z = m_points[0];
                w = m_weights;
            }

            inline void GetPoints(Array<OneD, const DataType> &x) const
            {
                x = m_points[0];
            }

            inline void GetPoints(Array<OneD, const DataType> &x,
                                  Array<OneD, const DataType> &y) const
            {
                x = m_points[0];
                y = m_points[1];
            }

            inline void GetPoints(Array<OneD, const DataType> &x,
                                  Array<OneD, const DataType> &y,
                                  Array<OneD, const DataType> &z) const
            {
                x = m_points[0];
                y = m_points[1];
                z = m_points[2];
            }

            inline const MatrixSharedPtrType& GetD(Direction dir = xDir) const
            {
                return m_derivmatrix[(int)dir];
            }

            virtual const MatrixSharedPtrType GetI(const PointsKey &pkey)
            {
                ASSERTL0(false, "Method not implemented ");
                boost::shared_ptr<NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr());
                return returnval;
            }
			
            virtual const MatrixSharedPtrType GetI(const Array<OneD, const DataType>& x)
            {
                ASSERTL0(false, "Method not implemented");
                boost::shared_ptr<NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr());
                return returnval;
            }
			
            virtual const MatrixSharedPtrType GetI(unsigned int numpoints, const Array<OneD, const DataType>& x)
            {
                ASSERTL0(false, "Method not implemented");
                boost::shared_ptr<NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr());
                return returnval;
            }
			
            virtual const MatrixSharedPtrType GetI(const Array<OneD, const DataType>& x, const Array<OneD, const DataType>& y)
            {
                ASSERTL0(false, "Method not implemented");
                boost::shared_ptr<NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr());
                return returnval;
            }

            virtual const MatrixSharedPtrType GetI(const Array<OneD, const DataType>& x, const Array<OneD, const DataType>& y,
                                                   const Array<OneD, const DataType>& z)
            {
                ASSERTL0(false, "Method not implemented");
                boost::shared_ptr<NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr());
                return returnval;
            }

            virtual const MatrixSharedPtrType GetGalerkinProjection(const PointsKey &pkey)
            {
                ASSERTL0(false, "Method not implemented ");
                boost::shared_ptr<NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr());
                return returnval;
            }

        protected:
            PointsKey             m_pointsKey;
            Array<OneD, DataType> m_points[3];
            Array<OneD, DataType> m_weights;
            MatrixSharedPtrType   m_derivmatrix[3];
            NekManager<PointsKey, NekMatrix<DataType>, PointsKey::opLess> m_InterpManager;
            NekManager<PointsKey, NekMatrix<DataType>, PointsKey::opLess> m_GalerkinProjectionManager;
            
            virtual void CalculatePoints()
            {
                unsigned int pointsDim = GetPointsDim();
                unsigned int totNumPoints = GetTotNumPoints();

                for (unsigned int i=0; i<pointsDim; ++i)
                {
                    m_points[i] = Array<OneD, DataType>(totNumPoints);
                }
            }

            virtual void CalculateWeights()
            {
                m_weights = Array<OneD, DataType>(GetTotNumPoints());
            }

            virtual void CalculateDerivMatrix()
            {
                int totNumPoints = GetTotNumPoints();
                for(unsigned int i = 0; i < m_pointsKey.GetPointsDim(); ++i)
                {
                    m_derivmatrix[i] = MemoryManager<NekMatrix<DataType> >::AllocateSharedPtr(totNumPoints,totNumPoints);
                }
            }
            Points(const PointsKey &key):m_pointsKey(key)
            {
            }

        private:
            // These should never be called
            Points(const Points &pts)
            {
                NEKERROR(ErrorUtil::efatal,"Copy Constructor for Points should not be called");
            }
            Points()
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for Points should not be called");
            }
        };

    }; // end of namespace
} // end of namespace 

#endif //NEKTAR_LIB_UTILITIES_FOUNDATIONS_POINTS_H
