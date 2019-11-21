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

#include <LibUtilities/Foundations/PolyEPoints.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        bool PolyEPoints::initPointsManager[] = {
            PointsManager().RegisterCreator(PointsKey(0, ePolyEvenlySpaced), PolyEPoints::Create)
        };

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
                NekDouble dx = 2.0/(NekDouble)(npts-1);
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
                std::shared_ptr<PointsBaseType> ptr = PointsManager()[gaussKey];
                Array<OneD, const NekDouble> z;
                Array<OneD, const NekDouble> w;

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
            PointsBaseType::CalculateDerivMatrix();

            for(unsigned int i=0;i<m_pointsKey.GetNumPoints();++i)
            {
                for(unsigned int j=0;j<m_pointsKey.GetNumPoints();++j)
                {
                    (*m_derivmatrix[0])(i,j) = LagrangePolyDeriv(m_points[0][i],j,m_pointsKey.GetNumPoints(),m_points[0]);
                }
            }
        }

        void PolyEPoints::CalculateInterpMatrix(unsigned int npts, const Array<OneD, const NekDouble>& xpoints, Array<OneD, NekDouble>& interp)
        {
            for(unsigned int i=0;i<npts;++i)
            {
                for(unsigned int j=0;j<m_pointsKey.GetNumPoints();++j)
                {
                    interp[i + j*npts] = LagrangePoly(xpoints[i],j,m_pointsKey.GetNumPoints(),m_points[0]);
                }
            }
        }

        std::shared_ptr<PolyEPoints::PointsBaseType> PolyEPoints::Create(const PointsKey &key)
        {
            std::shared_ptr<PointsBaseType> returnval(MemoryManager<PolyEPoints>::AllocateSharedPtr(key));

            returnval->Initialize();

            return returnval;
        }

        const std::shared_ptr<NekMatrix<NekDouble> > PolyEPoints::GetI(const PointsKey &pkey)
        {
            ASSERTL0(pkey.GetPointsDim()==1, "Gauss Points can only interp to other 1d point distributions");

            int numpoints = pkey.GetNumPoints();
            Array<OneD, const NekDouble> xpoints;

            PointsManager()[pkey]->GetPoints(xpoints);

            return GetI(numpoints, xpoints);
        }

        const std::shared_ptr<NekMatrix<NekDouble> > PolyEPoints::GetI(const Array<OneD, const NekDouble>& x)
        {
            int numpoints = 1;

            /// Delegate to function below.
            return GetI(numpoints, x);
        }

        const std::shared_ptr<NekMatrix<NekDouble> > PolyEPoints::GetI(unsigned int numpoints, const Array<OneD, const NekDouble>& x)
        {
            Array<OneD, NekDouble> interp(GetNumPoints()*numpoints);

            CalculateInterpMatrix(numpoints, x, interp);

            unsigned int np = GetTotNumPoints();
            NekDouble* d = interp.data();
            std::shared_ptr< NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(numpoints,np, d));

            return returnval;
        }


        NekDouble PolyEPoints::LagrangeInterpolant(NekDouble x, int npts, const Array<OneD, const NekDouble>& xpts,
            const Array<OneD, const NekDouble>& funcvals)
        {
            NekDouble sum = 0.0;

            for(int i=0;i<npts;++i)
            {
                sum += funcvals[i]*LagrangePoly(x,i,npts,xpts);
            }
            return sum;
        }


        NekDouble PolyEPoints::LagrangePoly(NekDouble x, int pt, int npts, const Array<OneD, const NekDouble>& xpts)
        {
            NekDouble h=1.0;

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

        NekDouble PolyEPoints::LagrangePolyDeriv(NekDouble x, int pt, int npts, const Array<OneD, const NekDouble>& xpts)
        {
            NekDouble h;
            NekDouble y=0.0;

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
                    }
                    y = y + h;
                }
            }
            return y;
        }
    } // end of namespace LibUtilities
} // end of namespace Nektar



