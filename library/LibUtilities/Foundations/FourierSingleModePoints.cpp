///////////////////////////////////////////////////////////////////////////////
//
// File FourierSingleModePoints.cpp
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
// Description: 1D Non Evenly-Spaced Fourier Points for stability analysis
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/Foundations/FourierSingleModePoints.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{
    namespace LibUtilities 
    {
        bool FourierSingleModePoints::initPointsManager[] = {
            PointsManager().RegisterCreator(PointsKey(0, eFourierSingleModeSpaced), FourierSingleModePoints::Create)
        };

        void FourierSingleModePoints::CalculatePoints()
        {
            // Allocate the storage for points
            PointsBaseType::CalculatePoints();
            unsigned int npts = m_pointsKey.GetNumPoints();

			if(npts==1)
            {
                m_points[0][0] = 0.25;
            }
			else
			
			{
            ASSERTL0(npts==2, "Fourier points for single mode analysis should be 2");
			
			m_points[0][0] = 0.0;
			m_points[0][1] = 0.5;
			}
		}

        void FourierSingleModePoints::CalculateWeights()
        {
            // Allocate the storage for points
            PointsBaseType::CalculateWeights();

            unsigned int npts = m_pointsKey.GetNumPoints();
				//Here I need to insert the weight for the new point distribution
                for(unsigned int i=0; i<npts; ++i)
                {
                    m_weights[i] =  1.0;
                }
        }

        void FourierSingleModePoints::CalculateDerivMatrix()
        {

			Points<NekDouble>::CalculateDerivMatrix();
		
		}

        std::shared_ptr<Points<NekDouble> > FourierSingleModePoints::Create(const PointsKey &key)
        {
            std::shared_ptr<Points<NekDouble> > returnval(MemoryManager<FourierSingleModePoints>::AllocateSharedPtr(key));

            returnval->Initialize();

            return returnval;
        }


        std::shared_ptr< NekMatrix<NekDouble> > FourierSingleModePoints::CreateMatrix(const PointsKey &pkey)
        {
            int numpoints = pkey.GetNumPoints();
            Array<OneD, const NekDouble> xpoints;

            PointsManager()[pkey]->GetPoints(xpoints);

            /// Delegate to function below.
            return GetI(numpoints, xpoints);
        }

        const std::shared_ptr<NekMatrix<NekDouble> > FourierSingleModePoints::GetI(const PointsKey& pkey)
        {
            ASSERTL0(pkey.GetPointsDim()==1, "Fourier Points can only interp to other 1d point distributions");

            return m_InterpManager[pkey];
			
        }

        const std::shared_ptr<NekMatrix<NekDouble> > FourierSingleModePoints::GetI(const Array<OneD, const NekDouble>& x)
        {
            int numpoints = 1;

            /// Delegate to function below.
            return GetI(numpoints, x);
        }

        const std::shared_ptr<NekMatrix<NekDouble> > FourierSingleModePoints::GetI(unsigned int numpoints, const Array<OneD, const NekDouble>& x)
        {
            Array<OneD, NekDouble> interp(GetNumPoints()*numpoints);

            CalculateInterpMatrix(numpoints, x, interp);

            NekDouble* t = interp.data();
            unsigned int np = GetNumPoints();
            std::shared_ptr< NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(numpoints,np,t));

            return returnval;
        }

        void FourierSingleModePoints::CalculateInterpMatrix(unsigned int npts, const Array<OneD, const NekDouble>& xpoints, Array<OneD, NekDouble>& interp)
        {
            boost::ignore_unused(npts, xpoints, interp);
        }


            
    } // end of namespace LibUtilities
} // end of namespace Nektar

