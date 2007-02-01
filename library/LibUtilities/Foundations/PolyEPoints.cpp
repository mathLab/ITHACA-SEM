///////////////////////////////////////////////////////////////////////////////
//
// File PolyEPoints.cpp
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
#include <LibUtilities/Foundations/PolyEPoints.h>

namespace Nektar
{
    namespace LibUtilities
    {

        void PolyEPoints::CalculatePoints()
        {
            // Allocate the storage for points
            PointsBaseType::CalculatePoints();

            unsigned int npts = m_pointsKey.GetNumPoints();
            if(npts==1)
            {
                m_points[0][0] = 0.0;
            }
            else
            {
                double dx = 2.0/(double)(npts-1);
                for(unsigned int i=0;i<npts;++i)
                {
                    m_points[0][i] = -1.0 + i*dx;
                }
            }
        }

        void PolyEPoints::CalculateWeights()
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
                PointsKey gaussKey(npts, eGaussLobattoLegendre);
                boost::shared_ptr<Points<double> > ptr = PointsManager()[gaussKey];
                const double *z,*w;
                ptr->GetZW(z,w);
                for(unsigned int i=0; i<npts;++i)
                {
                    m_weights[i] = 0.0;
                    for(unsigned j=0;j<npts;++j)
                    {
                        m_weights[i] += w[j]*LagrangePoly(z[j],i,npts,m_points[0]);
                    }
                }
            }
        }

        void PolyEPoints::CalculateDerivMatrix()
        {
            // Allocate the derivative matrix.
            Points<double>::CalculateDerivMatrix();

            for(unsigned int i=0;i<m_pointsKey.GetNumPoints();++i)
            {
                for(unsigned int j=0;j<m_pointsKey.GetNumPoints();++j)
                {
                    (*m_derivmatrix)(i,j) = LagrangePolyDeriv(m_points[0][i],j,m_pointsKey.GetNumPoints(),m_points[0]);
                }
            }

        }

        void PolyEPoints::CalculateInterpMatrix(unsigned int npts, const double * xpoints, double * interp)
        {
            for(unsigned int i=0;i<npts;++i)
            {
                for(unsigned int j=0;j<m_pointsKey.GetNumPoints();++j)
                {
                    interp[i*m_pointsKey.GetNumPoints()+j] = LagrangePoly(xpoints[i],j,m_pointsKey.GetNumPoints(),m_points[0]);
                }
            }

        }

        boost::shared_ptr< PointsBaseType > PolyEPoints::Create(const PointsKey &key)
        {
            boost::shared_ptr< PointsBaseType > returnval(new PolyEPoints(key));

            returnval->Initialize();

            return returnval;
        }

        const boost::shared_ptr<NekMatrix<double> > PolyEPoints::GetI(const PointsKey &pkey)
        {

            ASSERTL0(pkey.GetPointsDim()==1, "Gauss Points can only interp to other 1d point distributions");

            int numpoints = pkey.GetNumPoints();
            const double * xpoints;

            PointsManager()[pkey]->GetPoints(xpoints);

            return GetI(numpoints, xpoints);
        }

        const boost::shared_ptr<NekMatrix<double> > PolyEPoints::GetI(const double * x)
        {
            int numpoints = 1;

            /// Delegate to function below.
            return GetI(numpoints, x);
        }

        const boost::shared_ptr<NekMatrix<double> > PolyEPoints::GetI(unsigned int numpoints, const double *x)
        {
            double * interp = new double[GetNumPoints()*numpoints];

            CalculateInterpMatrix(numpoints, x, interp);

            boost::shared_ptr< NekMatrix<DataType> > returnval(new NekMatrix<DataType>(numpoints,GetNumPoints(),interp));

            delete[] interp;

            return returnval;
        }


        double PolyEPoints::LagrangeInterpolant(double x, int npts, double *xpts, double * funcvals)
        {
            double sum = 0.0;

            for(int i=0;i<npts;++i)
            {
                sum = sum + funcvals[i]*LagrangePoly(x,i,npts,xpts);
            }
            return sum;
        }


        double PolyEPoints::LagrangePoly(double x, int pt, int npts, double * xpts)
        {
            double h=1.0;

            for(int i=0;i<pt; ++i)
            {
                h = h * (x - xpts[i])/(xpts[pt]-xpts[i]);
            }

            for(int i=pt+1;i<npts;++i)
            {
                h = h * (x - xpts[i])/(xpts[pt]-xpts[i]);
            }

            return h;
        }

        double PolyEPoints::LagrangePolyDeriv(double x, int pt, int npts, double * xpts)
        {
            double h;
            double y=0.0;

            for(int j=0;j<npts;++j)
            {
                if(j!=pt)
                {
                    h=1.0;
                    for(int i=0;i<npts;++i)
                    {
                        if(i!=pt)
                        {
                            if(i!=j)
                            {
                                h = h*(x-xpts[i]);
                            }
                            h = h/(xpts[pt]-xpts[i]);
                        }
                        y = y + h;
                    }
                }
            }
            return y;
        }
    } // end of namespace LibUtilities
} // end of namespace Nektar



