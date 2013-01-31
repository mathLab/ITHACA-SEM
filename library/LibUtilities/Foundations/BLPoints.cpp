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


#include <LibUtilities/Foundations/BLPoints.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
//#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{
    namespace LibUtilities 
    {
        // Default value.
        NekDouble BLPoints::delta_star = 1.0;
        
        void BLPoints::CalculatePoints()
        {
            // Allocate the storage for points.
            PointsBaseType::CalculatePoints();
            unsigned int npts = m_pointsKey.GetNumPoints(); 

	    // Derived power coefficient.
	    NekDouble rn = powf(2.0/BLPoints::delta_star,1.0/(npts-2.0));
            
            m_points[0][0] = -1.0;
            
	    for (unsigned int i = 1; i < npts; ++i)
            {
                m_points[0][i] = -1.0 + delta_star*pow(rn,(NekDouble)i-1.0);
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
    } // end of namespace LibUtilities
} // end of namespace Nektar

