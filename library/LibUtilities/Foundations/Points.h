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

#include <math.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        // Need to add method to compute total number of points given dimension
        // and number of points.

        class PointsKey
        {
        public:
            // Used for looking up the creator.  The creator for number of points
            // can generate for any number, so we want the same creator called
            // for all number.
            struct opLess
            {
                bool operator()(const PointsKey &lhs, const PointsKey &rhs) const;
            };

            PointsKey(const int &numpoints, const PointsType &pointstype): 
                m_numpoints(numpoints), 
                m_pointstype(pointstype)
            {
            }

            virtual ~PointsKey()
            {
            }

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
                    dimpoints = 2;
                    break;

                case eNodalTetElec:
                    dimpoints = 3;
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
                    totpoints = m_numpoints*(m_numpoints+1)/2;
                    break;

                case eNodalTetElec:
                    totpoints = m_numpoints*(m_numpoints+1)*(m_numpoints+2)/6;
                    break;
                }

                return totpoints;
            }

            friend bool operator==(const PointsKey &lhs, const PointsKey &rhs);
            friend bool operator<(const PointsKey &lhs, const PointsKey &rhs);
            friend bool opLess::operator()(const PointsKey &lhs, const PointsKey &rhs) const;

        protected:
            unsigned int m_numpoints;     //!< number of the points (as appropriately defined for PointsType)
            PointsType m_pointstype;      //!< Type of Points

        private:

            // This should never be called
            PointsKey()
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for PointsKey should not be called");
            }

        };

        static const PointsKey NullPointsKey(0, eNoPointsType);

        bool operator==(const PointsKey &lhs, const PointsKey &rhs);
        bool operator<(const PointsKey &lhs, const PointsKey &rhs);
        std::ostream& operator<<(std::ostream& os, const PointsKey& rhs);

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

            inline const ConstArray<OneD, DataType>& GetZ() const
            {
                return m_points[0];
            }

            inline const ConstArray<OneD, DataType>& GetW() const 
            {
                return m_weights; 
            } 

            inline void GetZW(ConstArray<OneD, DataType> &z,
                ConstArray<OneD, DataType> &w) const 
            {
                z = m_points[0];
                w = m_weights;
            }

            inline void GetPoints(ConstArray<OneD, DataType> &x) const
            {
                x = m_points[0];
            }

            inline void GetPoints(ConstArray<OneD, DataType> &x,
                                  ConstArray<OneD, DataType> &y) const
            {
                x = m_points[0];
                y = m_points[1];
            }

            inline void GetPoints(ConstArray<OneD, DataType> &x,
                                  ConstArray<OneD, DataType> &y,
                                  ConstArray<OneD, DataType> &z) const
            {
                x = m_points[0];
                y = m_points[1];
                z = m_points[2];
            }

            inline const MatrixSharedPtrType GetD() const
            {
                return m_derivmatrix;
            }

            virtual const MatrixSharedPtrType GetI(const PointsKey &pkey)=0;
            virtual const MatrixSharedPtrType GetI(ConstArray<OneD, DataType>& x) = 0;
            virtual const MatrixSharedPtrType GetI(unsigned int numpoints, ConstArray<OneD, DataType>& x) = 0;

        protected:
            PointsKey m_pointsKey;
            Array<OneD, DataType> m_points[3];
            Array<OneD, DataType> m_weights;
            MatrixSharedPtrType m_derivmatrix;
            NekManager<PointsKey, NekMatrix<DataType>, PointsKey::opLess> m_InterpManager;

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
                m_derivmatrix.reset(MemoryManager<NekMatrix<DataType> >::Allocate(totNumPoints,totNumPoints));
            }

            Points(const PointsKey &key):m_pointsKey(key)
            {
            }

        private:
            // These should never be called
            Points(const Points &pts);
            Points();
        };

    }; // end of namespace
} // end of namespace 

#endif //NEKTAR_LIB_UTILITIES_FOUNDATIONS_POINTS_H
