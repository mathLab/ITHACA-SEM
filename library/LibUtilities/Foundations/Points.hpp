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
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{
    namespace LibUtilities
    {

        typedef boost::shared_ptr<NekMatrix> NekMatrixSharedPtr;

        class PointsKey
        {
        public:

            PointsKey()
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for PointsKey should not be called");
            }

            PointsKey(int numpoints, PointsType pointstype, PointsIdentifier pointsid)
            {
                m_numpoints = numpoints;
                m_pointstype = pointstype;
                m_pointsid = pointsid;
            }

            virtual ~PointsKey()
            {
            }

            PointsKey(const PointsKey &key)
            {
                *this = key;
            }

            PointsKey& operator=(const PointsKey &key)
            {
                m_numpoints = key.m_numpoints;
                m_pointstype = key.m_pointstype;
                m_pointsid = key.m_pointsid;
                return *this;
            }

            inline int GetNumPoints() const
            {
                return m_numpoints;
            }

            inline PointsType GetPointsType() const
            {
                return m_pointstype;
            }

            inline PointsIdentifier GetPointsId() const
            {
                return m_pointsid;
            }

        protected:
            int m_numpoints;              //!< Number of points
            PointsType m_pointstype       //!< Type of Points
            PointsIdentifier m_pointsid;  //!< Unique indentifier (when needed) 
        };

        template<typename DataType, unsigned int dim>
        class Points
        {
        public:
            Points()
            {
                NEKERROR(ErrorUtil::efatal,"Default Constructor for Points should not be called");
            }

            Points(const PointsKey &key)
            {
                m_pkey = key;

                     for(unsigned int i = 0; i < dim; ++i)
                    {
                        delete[] m_points[i];
                    }
                    delete[] m_weights;
                    key.m_numpoints = 0;
                }

            }


           friend bool operator  == (const BasisKey& x, const BasisKey& y);

            virtual ~Points()
            {
                if(m_numpoints)
                {
                    for(unsigned int i = 0; i < dim; ++i)
                    {
                        delete[] m_points[i];
                    }
                    delete[] m_weights;
                    key.m_numpoints = 0;
                }
            }
  
            inline int GetNumPoints() const
            {
                return key.GetNumPoints();
            }

            inline PointsType GetPointsType() const
            {
                return key.GetPointsType();
            }

            inline PointsIdentifier GetPointsId() const
            {
                return key.GetPointsId();
            }

            inline double *GetZ() const
            {
                BOOST_STATIC_ASSERT(dim == 1);
                return m_points[0];
            }

            inline double *GetW() const 
            {
                return m_weights; 
            } 
      
            inline void GetZW(const double *&z, const double *&w) const 
            {
                BOOST_STATIC_ASSERT(dim == 1);
	            z = m_points[0];
            	w = m_weights;
            }

            inline void GetPoints(const double *&x) const
            {
                BOOST_STATIC_ASSERT(dim == 1);
                x = m_points[0];
            }

            inline void GetPoints(const double *&x, const double *&y) const
            {
                BOOST_STATIC_ASSERT(dim == 2);
                x = m_points[0];
                y = m_points[1];
            }

            inline void GetPoints(const double *&x, const double *&y, const double *&z) const
            {
                BOOST_STATIC_ASSERT(dim == 3);
                x = m_points[0];
                y = m_points[1];
                z = m_points[2];
            }

            inline const NekMatrix<DataType,FULL> *GetD() const
            {
                return m_derivmatrix;
            }
        
        protected:
            PointsKey m_pointskey;
            DataType * m_points[dim];
            DataType * m_weights;
            NekMatrix<DataType,eFull> m_derivmatrix;

        private:
            virtual void CalculatePoints() = 0;
            virtual void CalculateWeights() = 0;
            virtual void CalculateDerivMatrix() = 0;
        };


        bool operator == (const PointsKey& x, const PointsKey& y)
        {
            if( (x.m_numpoints == y.m_numpoints) &&
                (x.m_pointstype == y.m_pointstype) &&
                (x.m_pointsid == y.m_pointsid) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }


        bool operator == (const PointsKey* x, const PointsKey& y)
        {
            if( ((*x).m_numpoints == y.m_numpoints) &&
                ((*x).m_pointstype == y.m_pointstype) &&
                ((*x).m_pointsid == y.m_pointsid) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }


        bool operator == (const PointsKey& x, const PointsKey *y)
        {
            if( (x.m_numpoints == (*y).m_numpoints) &&
                (x.m_pointstype == (*y).m_pointstype) &&
                (x.m_pointsid == (*y).m_pointsid) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }


        bool operator != (const PointsKey& x, const PointsKey& y)
        {
            if( (x.m_numpoints == y.m_numpoints) &&
                (x.m_pointstype == y.m_pointstype) &&
                (x.m_pointsid == y.m_pointsid) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }


        bool operator != (const PointsKey* x, const PointsKey& y)
        {
            if( ((*x).m_numpoints == y.m_numpoints) &&
                ((*x).m_pointstype == y.m_pointstype) &&
                ((*x).m_pointsid == y.m_pointsid) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }


        bool operator != (const PoinstKey& x, const PointsKey *y)
        {
            if( (x.m_numpoints == (*y).m_numpoints) &&
                (x.m_pointstype == (*y).m_pointstype) &&
                (x.m_pointsid == (*y).m_pointsid) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }

    } // end of namespace
} // end of namespace 

#endif //POINTS_HPP
