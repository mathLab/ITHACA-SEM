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

#include <LibUtilities/Foundations/FourierPoints.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>


#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc

namespace Nektar
{
    namespace LibUtilities 
    {
        bool FourierPoints::initPointsManager[] = {
            PointsManager().RegisterCreator(PointsKey(0, eFourierEvenlySpaced),     FourierPoints::Create)
        };

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
                NekDouble dx = 2.0/(NekDouble)(npts);
                for(unsigned int i=0;i<npts;++i)
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
                m_weights[0] = 1.0; //midpoint rule
            }
            else
            {
                for(unsigned int i=0; i<npts; ++i)
                {
                    m_weights[i] =  1.0/(NekDouble)npts;
                }
            }
        }

        void FourierPoints::CalculateDerivMatrix()
        {
            //// Allocate the derivative matrix.
            Points<NekDouble>::CalculateDerivMatrix();

            unsigned int npts = m_pointsKey.GetNumPoints();

            for(unsigned int i=1;i<npts;++i)
            {
                m_derivmatrix[0]->SetValue(i,i, 0.0);
            }

            for(unsigned int i=1;i<npts;++i)
            {
                m_derivmatrix[0]->SetValue(0,i, -0.5*M_PI*pow(-1.0,NekDouble(i))*cos(M_PI*i/npts)/sin(M_PI*i/npts));
            }

            for(unsigned int i=1;i<npts;++i)
            {
                for(unsigned int j=0;j<npts;++j)
                {
                    m_derivmatrix[0]->SetValue(i,j, (*m_derivmatrix[0])(0,(j-i+npts)%npts));
                }
            }           
        }

        std::shared_ptr<Points<NekDouble> > FourierPoints::Create(const PointsKey &key)
        {
            std::shared_ptr<Points<NekDouble> > returnval(MemoryManager<FourierPoints>::AllocateSharedPtr(key));

            returnval->Initialize();

            return returnval;
        }


        std::shared_ptr< NekMatrix<NekDouble> > FourierPoints::CreateMatrix(const PointsKey &pkey)
        {
            int numpoints = pkey.GetNumPoints();
            Array<OneD, const NekDouble> xpoints;

            PointsManager()[pkey]->GetPoints(xpoints);

            /// Delegate to function below.
            return GetI(numpoints, xpoints);
        }

        const std::shared_ptr<NekMatrix<NekDouble> > FourierPoints::GetI(const PointsKey& pkey)
        {
            ASSERTL0(pkey.GetPointsDim()==1, "Fourier Points can only interp to other 1d point distributions");

            return m_InterpManager[pkey];
        }

        const std::shared_ptr<NekMatrix<NekDouble> > FourierPoints::GetI(const Array<OneD, const NekDouble>& x)
        {
            int numpoints = 1;

            /// Delegate to function below.
            return GetI(numpoints, x);
        }

        const std::shared_ptr<NekMatrix<NekDouble> > FourierPoints::GetI(unsigned int numpoints, const Array<OneD, const NekDouble>& x)
        {
            Array<OneD, NekDouble> interp(GetNumPoints()*numpoints);

            CalculateInterpMatrix(numpoints, x, interp);

            NekDouble* t = interp.data();
            unsigned int np = GetNumPoints();
            std::shared_ptr< NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(numpoints,np,t));

            return returnval;
        }

        void FourierPoints::CalculateInterpMatrix(unsigned int npts, const Array<OneD, const NekDouble>& xpoints, Array<OneD, NekDouble>& interp)
        {
            const NekDouble h = 2.0/m_pointsKey.GetNumPoints();
            for(unsigned int i=0;i<npts;++i)
            {
                for(unsigned int j=0;j<m_pointsKey.GetNumPoints();++j)
                {
                    interp[i*m_pointsKey.GetNumPoints()+j] = PeriodicSincFunction(M_PI*(xpoints[i]-m_points[0][j]),h);
                }
            }
        }

        NekDouble FourierPoints::PeriodicSincFunction(const NekDouble x, const NekDouble h)
        {
            // This formula is based upon a mapped version of 
            // the periodic sinc presented in Trefethen's "Spectral Methods
            // in Matlab"

            NekDouble y = 1.0;
            
            if(fabs(x)>1.0e-12)
            {
                y = sin(M_PI*x/(M_PI*h))/((2.0/h)*tan(0.5*x));
            }

            return y;
        }
            
    } // end of namespace LibUtilities
} // end of namespace Nektar

