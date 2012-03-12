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
// Description: 1D boundary layer points
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LibUtilities.h>
#include <iostream>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Foundations.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Foundations/BLPoints.h>

namespace Nektar
{
    namespace LibUtilities 
    {
        void BLPoints::CalculatePoints()
        {
            // Allocate the storage for points
            PointsBaseType::CalculatePoints();
            unsigned int npts = m_pointsKey.GetNumPoints();
	    double rr = pow(10,log10(6000)/(npts-2));
	    //cerr << "power coefficiant :  " << rr << endl;
	    for(unsigned int i=0;i<npts;++i)
            {
		if ( i == 0 )
                {
		    m_points[0][i] = -1.0;
		    //cerr << m_points[0][i] +1.0 << endl;
                }		 
		else
                {
		    m_points[0][i] = -1 + (1.0/3000.0)*pow(rr,(double(i-1)));
		    //cerr << m_points[0][i] +1.0 << endl;
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

        boost::shared_ptr<Points<NekDouble> > BLPoints::Create(const PointsKey &key)
        {
            boost::shared_ptr<Points<NekDouble> > returnval(MemoryManager<BLPoints>::AllocateSharedPtr(key));

            returnval->Initialize();

            return returnval;
        }


        boost::shared_ptr< NekMatrix<NekDouble> > BLPoints::CreateMatrix(const PointsKey &pkey)
        {
            int numpoints = pkey.GetNumPoints();
            Array<OneD, const NekDouble> xpoints;

            PointsManager()[pkey]->GetPoints(xpoints);

            /// Delegate to function below.
            return GetI(numpoints, xpoints);
        }

        const boost::shared_ptr<NekMatrix<NekDouble> > BLPoints::GetI(const PointsKey& pkey)
        {
            ASSERTL0(pkey.GetPointsDim()==1, "Fourier Points can only interp to other 1d point distributions");

            return m_InterpManager[pkey];
        }

        const boost::shared_ptr<NekMatrix<NekDouble> > BLPoints::GetI(const Array<OneD, const NekDouble>& x)
        {
            int numpoints = 1;

            /// Delegate to function below.
            return GetI(numpoints, x);
        }

        const boost::shared_ptr<NekMatrix<NekDouble> > BLPoints::GetI(unsigned int numpoints, const Array<OneD, const NekDouble>& x)
        {
            Array<OneD, NekDouble> interp(GetNumPoints()*numpoints);

            CalculateInterpMatrix(numpoints, x, interp);

            NekDouble* t = interp.data();
            unsigned int np = GetNumPoints();
            boost::shared_ptr< NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr(numpoints,np,t));

            return returnval;
        }

        void BLPoints::CalculateInterpMatrix(unsigned int npts, const Array<OneD, const NekDouble>& xpoints, Array<OneD, NekDouble>& interp)
        {

        }

        double BLPoints::PeriodicSincFunction(const NekDouble x, const NekDouble h)
        {

        }
            
    } // end of namespace LibUtilities
} // end of namespace Nektar

