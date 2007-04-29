///////////////////////////////////////////////////////////////////////////////
//
// File FourierPoints.cpp
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
// Description: 1D Evenly-Spaced Points 
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Foundations.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Foundations/FourierPoints.h>

namespace Nektar
{
    namespace LibUtilities 
    {

        template<typename T>
        void FourierPoints<T>::CalculatePoints()
        {
            // Allocate the storage for points
            PointsBaseType::CalculatePoints();

            unsigned int npts = m_pointsKey.GetNumPoints();
            ASSERTL0(!(npts%2), "Fourier points need to be of even order");

            if(npts==1)
            {
                (*m_points[0])[0] = 0.0;
            }
            else
            {
                double dx = 2.0/(double)(npts);
                for(unsigned int i=0;i<npts;++i)
                {
                    (*m_points[0])[i] = -1.0 + i*dx;
                }
            }
        }

        template<typename T>
        void FourierPoints<T>::CalculateWeights()
        {
            // Allocate the storage for points
            PointsBaseType::CalculateWeights();

            unsigned int npts = m_pointsKey.GetNumPoints();
            if(npts==1)
            {
                (*m_weights)[0] = 2.0; //midpoint rule
            }
            else
            {
                for(unsigned int i=0; i<npts; ++i)
                {
                    (*m_weights)[i] =  2.0/(double)npts;
                }
            }
        }

        template<typename T>
        void FourierPoints<T>::CalculateDerivMatrix()
        {
            //// Allocate the derivative matrix.
            Points<double>::CalculateDerivMatrix();

            unsigned int npts = m_pointsKey.GetNumPoints();

            for(unsigned int i=1;i<npts;++i)
            {
                (*m_derivmatrix)(i,i) = 0.0;
            }

            for(unsigned int i=1;i<npts;++i)
            {
                (*m_derivmatrix)(0,i) = -0.5*M_PI*pow(-1.0,double(i))*cos(M_PI*i/npts)/sin(M_PI*i/npts);
            }

            for(unsigned int i=1;i<npts;++i)
            {
                for(unsigned int j=0;j<npts;++j)
                {
                    (*m_derivmatrix)(i,j) = (*m_derivmatrix)(0,(j-i+npts)%npts);
                }
            }           
        }

        template<typename T>
        boost::shared_ptr<typename FourierPoints<T>::PointsBaseType> FourierPoints<T>::Create(const PointsKey &key)
        {
            boost::shared_ptr<PointsBaseType> returnval(MemoryManager::AllocateSharedPtr<FourierPoints<T> >(key));

            returnval->Initialize();

            return returnval;
        }


        template<typename T>
        boost::shared_ptr< NekMatrix<T> > FourierPoints<T>::CreateMatrix(const PointsKey &pkey)
        {
            int numpoints = pkey.GetNumPoints();
            typename Nek1DConstSharedArray<T>::type xpoints;

            PointsManager()[pkey]->GetPoints(xpoints);

            /// Delegate to function below.
            return GetI(numpoints, xpoints);
        }

        template<typename T>
        const boost::shared_ptr<NekMatrix<T> > FourierPoints<T>::GetI(const PointsKey& pkey)
        {
            ASSERTL0(pkey.GetPointsDim()==1, "Fourier Points can only interp to other 1d point distributions");

            return m_InterpManager[pkey];
        }

        template<typename T>
        const boost::shared_ptr<NekMatrix<T> > FourierPoints<T>::GetI(typename Nek1DConstSharedArray<T>::type &x)
        {
            int numpoints = 1;

            /// Delegate to function below.
            return GetI(numpoints, x);
        }

        template<typename T>
        const boost::shared_ptr<NekMatrix<T> > FourierPoints<T>::GetI(unsigned int numpoints, typename Nek1DConstSharedArray<T>::type &x)
        {
            Nek1DSharedArray<T> interp = MemoryManager::AllocateSharedPtr<Nek1DArray<T> >(boost::extents[GetNumPoints()*numpoints]);

            CalculateInterpMatrix(numpoints, x, interp);

            boost::shared_ptr< NekMatrix<T> > returnval(MemoryManager::AllocateSharedPtr<NekMatrix<T> >(numpoints,GetNumPoints(),interp.get()));

            return returnval;
        }

        template<typename T>
        void FourierPoints<T>::CalculateInterpMatrix(unsigned int npts, typename Nek1DConstSharedArray<T>::type xpoints, typename Nek1DSharedArray<T>::type interp)
        {
            const T h = 2.0/m_pointsKey.GetNumPoints();
            for(unsigned int i=0;i<npts;++i)
            {
                for(unsigned int j=0;j<m_pointsKey.GetNumPoints();++j)
                {
                    (*interp)[i*m_pointsKey.GetNumPoints()+j] = PeriodicSincFunction(M_PI*(xpoints[i]-(*m_points[0])[j]),h);
                }
            }
        }

        template<typename T>
        double FourierPoints<T>::PeriodicSincFunction(const T x, const T h)
        {
            // This formula is based upon a mapped version of 
            // the periodic sinc presented in Trefethen's "Spectral Methods
            // in Matlab"

            double y = 1.0;
            
            if(fabs(x)>1.0e-12)
            {
                y = sin(M_PI*x/(M_PI*h))/((2.0/h)*tan(0.5*x));
            }

            return y;
        }
            
    } // end of namespace LibUtilities
} // end of namespace Nektar

