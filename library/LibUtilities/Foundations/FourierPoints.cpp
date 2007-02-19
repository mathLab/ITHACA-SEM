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

        void FourierPoints::CalculatePoints()
        {
            // Allocate the storage for points
            PointsBaseType::CalculatePoints();

            unsigned int npts = m_pointsKey.GetNumPoints();
            ASSERTL0(!(npts%2), "Fourier points need to be of even order");

            if(npts==1)
            {
                m_points[0][0] = 0.0;
            }
            else
            {
                double dx = 2.0/(double)(npts);
                for(int i=0;i<npts;++i)
                {
                    m_points[0][i] = -1.0 + i*dx;
                }
            }
        }

        void FourierPoints::CalculateWeights()
        {
            // Allocate the storage for points
            PointsBaseType::CalculateWeights();

            unsigned int npts = m_pointsKey.GetNumPoints();
            if(npts==1)
            {
                m_weights[0] = 2.0; //midpoint rule
            }
            else
            {
                for(unsigned int i=0; i<npts; ++i)
                {
                    m_weights[i] =  2.0/(double)npts;
                }
            }
        }

        void FourierPoints::CalculateDerivMatrix()
        {
            //// Allocate the derivative matrix.
            Points<double>::CalculateDerivMatrix();

            int npts = m_pointsKey.GetNumPoints();

            for(unsigned int i=0;i<npts;++i)
            {
                for(unsigned int j=0;j<npts;++j)
                {
                    if(i==j)
                    {
                        (*m_derivmatrix)(i,j) = 0.0;
                    }
                    else
                    {
                        int k = (-2*((i+j)%2)) + 1;
                        (*m_derivmatrix)(i,j) = M_PI*(double(k)/2.0)*(cos(M_PI*double(i-j)/double(npts)))/sin(M_PI*double(i-j)/double(npts));
                    }
                }
            }
        }



        boost::shared_ptr< PointsBaseType > FourierPoints::Create(const PointsKey &key)
        {
            boost::shared_ptr< PointsBaseType > returnval(new FourierPoints(key));

            returnval->Initialize();

            return returnval;
        }

        const boost::shared_ptr<NekMatrix<double> > FourierPoints::GetI(const PointsKey &pkey)
        {

            boost::shared_ptr< NekMatrix<DataType> > returnval(new NekMatrix<DataType>(pkey.GetNumPoints(),GetNumPoints()));
            return returnval;
        }

        const boost::shared_ptr<NekMatrix<double> > FourierPoints::GetI(const double * x)
        {
            int numpoints = 1;

            /// Delegate to function below.
            return GetI(numpoints, x);
        }

        const boost::shared_ptr<NekMatrix<double> > FourierPoints::GetI(unsigned int numpoints, const double *x)
        {

            boost::shared_ptr< NekMatrix<DataType> > returnval(new NekMatrix<DataType>(numpoints,GetNumPoints()));

            return returnval;
        }       


        void FourierPoints::CalculateInterpMatrix(unsigned int npts, const double * xpoints, double * interp)
        {
            const double h = 2.0/m_pointsKey.GetNumPoints();
            for(unsigned int i=0;i<npts;++i)
            {
                for(unsigned int j=0;j<m_pointsKey.GetNumPoints();++j)
                {
                    interp[i*m_pointsKey.GetNumPoints()+j] = PeriodicSincFunction(xpoints[i]-m_points[0][j],h);
                }
            }
        }

        double FourierPoints::PeriodicSincFunction(const double h, const double x)
        {
            // This formula is based upon a mapped version of 
            // the periodic sinc presented in Trefethen's "Spectral Methods
            // in Matlab"

            double xi = -M_PI*(x+1.0)+2*M_PI;
            return sin(M_PI*xi/h)/((2.0*M_PI/h)*tan(x/2.0));
        }
    } // end of namespace LibUtilities
} // end of namespace Nektar



