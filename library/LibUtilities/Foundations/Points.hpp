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

#ifndef POINTS_HPP
#define POINTS_HPP

#include <math.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
//#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        // Need to add method to compute total number of points given dimension
        // and number of points.

        class PointsKey
        {
        public:
             // Use for looking up the creator.  The creator for number of points
            // can generate for any number, so we want the same creator called
            // for all number.
            struct opLess
            {
                bool operator()(const PointsKey &lhs, const PointsKey &rhs);
            };

           PointsKey(const int &numpoints, const PointsType &pointstype): 
                    m_numpoints(numpoints), 
                    m_pointstype(pointstype);
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

            inline unsigned int GetTotNumPoints() const
            {
                exit(0);//return m_numpoints;
            }

            inline PointsType GetPointsType() const
            {
                return m_pointstype;
            }


            bool operator==(const PointsKey &key)
            {
                return (m_numpoints == key.m_numpoints &&
                        m_pointstype == key.m_pointstype);
            }


            bool operator == (const PointsKey *y)
            {
                return (*this == *y);
            }

            bool operator != (const PointsKey& y)
            {
                return (!(*this == y));
            }

            bool operator != (const PointsKey *y)
            {
                return (!(*this == *y));
            }

            friend bool operator<(const PointsKey &lhs, const PointsKey &rhs);
            friend bool opLess::operator()(const PointsKey &lhs, const PointsKey &rhs);

        protected:
            unsigned int m_numpoints;     //!< number of the points (as appropriately defined for PointsType)
            PointsType m_pointstype;      //!< Type of Points

        private:
            PointsKey()
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for PointsKey should not be called");
            }
        };

        bool operator<(const PointsKey &lhs, const PointsKey &rhs);

        std::ostream& operator<<(std::ostream& os, const PointsKey& rhs);

        template<typename DataT>
        class Points
        {
        public:
            typedef DataT DataType;

            Points(const PointsKey &key): m_pointsKey(key)
            {
            }

            virtual ~Points()
            {
                //if(m_pointsKey.GetNumPoints())
                //{
                //    unsigned int dim = m_pkey.GetPointsDim();

                //    for(unsigned int i = 0; i < dim; ++i)
                //    {
                //        delete[] m_points[i];
                //    }

                //    delete[] m_points;
                //    delete[] m_weights;
                //}
            }

            void Initialize(void)
            {
                CalculatePoints();
                CalculateWeights();
                CalculateDerivMatrix();
            }

            inline unsigned int GetPointsDim() const
            {
                return 1;
            }

            inline unsigned int GetNumPoints() const
            {
                return m_pointsKey.GetNumPoints();
            }

            inline PointsType GetPointsType() const
            {
                return m_pointsKey.GetPointsType();
            }


            inline double *GetZ() const
            {
                return m_points[0];
            }

            inline double *GetW() const 
            {
                return m_weights; 
            } 


            inline void GetZW(const double *&z, const double *&w) const 
            {
                z = m_points[0];
                w = m_weights;
            }

            inline void GetPoints(const double *&x) const
            {
                x = m_points[0];
            }

            inline void GetPoints(const double *&x, const double *&y) const
            {
                x = m_points[0];
                y = m_points[1];
            }

            inline void GetPoints(const double *&x, const double *&y, const double *&z) const
            {
                x = m_points[0];
                y = m_points[1];
                z = m_points[2];
            }


            //inline const boost::shared_ptr<NekMatrix<DataType> > GetI(PointsKey)
            //{
            //    return m_derivmatrix;
            //}


//            inline const boost::shared_ptr<NekMatrix<DataType> > GetI(PointsKey)
//            {
//                return m_derivmatrix;
//            }

            inline const boost::shared_ptr<NekMatrix<DataType> > GetD() const
            {
                return m_derivmatrix;
            }

            inline const boost::shared_ptr<NekMatrix<DataType> > GetI(double x) const
            {
                return m_derivmatrix;
            }

            inline const boost::shared_ptr<NekMatrix<DataType> > GetI(unsigned int numpoints, double *x) const
            {
                return m_derivmatrix;
            }

            inline const boost::shared_ptr<NekMatrix<DataType> > GetI(PoinstsKey pkey) const
            {
                return m_derivmatrix;
            }

        protected:
            PointsKey m_pointsKey;
            DataType **m_points;
            DataType *m_weights;
//            boost::shared_ptr<NekMatrix<DataType> > m_derivmatrix;

            virtual void CalculatePoints()
            {
 /*               unsigned int pointsDim = GetPointsDim();

                m_points = new Points::DataType*[pointsDim];

                for (unsigned int i=0; i<pointsDim; ++i)
                {
                    m_points[i] = new Points::DataType[GetNumPoints()];
                }*/
            }

            virtual void CalculateWeights()
            {
                //m_weights = new Points::DataType[GetNumPoints()];
            }

            virtual void CalculateDerivMatrix()
            {
            }

        private:
            // This should never be called.
            Points()
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for Points should not be called");
            }
        };

    }; // end of namespace
} // end of namespace 

#endif //POINTS_HPP
