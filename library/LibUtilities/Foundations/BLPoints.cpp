///////////////////////////////////////////////////////////////////////////////
//
// File BLPoints.cpp
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
// Description: 1D boundary layer points
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Foundations/BLPoints.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

namespace Nektar
{
    namespace LibUtilities 
    {
        bool BLPoints::initPointsManager[] = {
            PointsManager().RegisterCreator(PointsKey(0, eBoundaryLayerPoints),    BLPoints::Create),
            PointsManager().RegisterCreator(PointsKey(0, eBoundaryLayerPointsRev), BLPoints::Create)
        };

        void BLPoints::CalculatePoints()
        {
            // Allocate the storage for points.
            PointsBaseType::CalculatePoints();
            unsigned int npts = m_pointsKey.GetNumPoints(); 

	    // Derived power coefficient.
            NekDouble r = m_pointsKey.GetFactor();
            ASSERTL0(r != NekConstants::kNekUnsetDouble,
                     "Must set factor in BLPoints key");

            if (fabs(r-1.0) < 1e-6)
            {
                NekDouble tmp = 2.0/(npts-1.0);
                for (unsigned int i = 0; i < npts; ++i)
                {
                    m_points[0][i] = -1.0 + i * tmp;
                }
            }
            else
            {
                NekDouble a = 2.0 * (1.0-r) / (1.0 - pow(r,(double)(npts-1)));
                m_points[0][0] = -1.0;
                
                for (unsigned int i = 1; i < npts; ++i)
                {
                    m_points[0][i] = m_points[0][i-1] + a*pow(r,(double)(i-1));
                }

                m_points[0][npts-1] = 1.0;
            }
            
            if (m_pointsKey.GetPointsType() == eBoundaryLayerPointsRev)
            {
                std::vector<NekDouble> tmp(npts);
                for (unsigned int i = 0; i < npts; ++i)
                {
                    tmp[i] = - m_points[0][npts-1-i];
                }
                
                for (unsigned int i = 0; i < npts; ++i)
                {
                    m_points[0][i] = tmp[i];
                }
            }
        }

        void BLPoints::CalculateWeights()
        {
            
        }

        void BLPoints::CalculateDerivMatrix()
        {
            //// Allocate the derivative matrix.
            Points<NekDouble>::CalculateDerivMatrix();
        }

        std::shared_ptr<Points<NekDouble> > BLPoints::Create(const PointsKey &key)
        {
            std::shared_ptr<Points<NekDouble> > returnval(MemoryManager<BLPoints>::AllocateSharedPtr(key));

            returnval->Initialize();

            return returnval;
        }


        std::shared_ptr< NekMatrix<NekDouble> > BLPoints::CreateMatrix(const PointsKey &pkey)
        {
            int numpoints = pkey.GetNumPoints();
            Array<OneD, const NekDouble> xpoints;

            PointsManager()[pkey]->GetPoints(xpoints);

            /// Delegate to function below.
            return GetI(numpoints, xpoints);
        }

        const std::shared_ptr<NekMatrix<NekDouble> > BLPoints::GetI(const PointsKey& pkey)
        {
            ASSERTL0(pkey.GetPointsDim()==1, "Fourier Points can only interp to other 1d point distributions");

            return m_InterpManager[pkey];
        }

        const std::shared_ptr<NekMatrix<NekDouble> > BLPoints::GetI(const Array<OneD, const NekDouble>& x)
        {
            int numpoints = 1;

            /// Delegate to function below.
            return GetI(numpoints, x);
        }

        const std::shared_ptr<NekMatrix<NekDouble> > BLPoints::GetI(unsigned int numpoints, const Array<OneD, const NekDouble>& x)
        {
            Array<OneD, NekDouble> interp(GetNumPoints()*numpoints);

            CalculateInterpMatrix(numpoints, x, interp);

            NekDouble* t = interp.data();
            unsigned int np = GetNumPoints();
            std::shared_ptr< NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(numpoints,np,t));

            return returnval;
        }

        void BLPoints::CalculateInterpMatrix(unsigned int npts, const Array<OneD, const NekDouble>& xpoints, Array<OneD, NekDouble>& interp)
        {
            boost::ignore_unused(npts, xpoints, interp);
        }
    } // end of namespace LibUtilities
} // end of namespace Nektar

