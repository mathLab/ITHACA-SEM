///////////////////////////////////////////////////////////////////////////////
//
// File FourierPoints.h
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
// Description: Header file of 1D Evenly-Spaced Point Definitions 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FOURIERPOINTS_H
#define FOURIERPOINTS_H

#include <math.h>
#include <boost/shared_ptr.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

namespace Nektar
{
    namespace LibUtilities 
    {
        typedef Points<double> PointsBaseType;

        class FourierPoints: public PointsBaseType
        {
        public:
            virtual ~FourierPoints()
            {
            }            

            static boost::shared_ptr< PointsBaseType > Create(const PointsKey &key);

            const boost::shared_ptr<NekMatrix<double> > GetI(const PointsKey &pkey);
            const boost::shared_ptr<NekMatrix<double> > GetI(const double * x);
            const boost::shared_ptr<NekMatrix<double> > GetI(unsigned int numpoints, const double *x);

        protected:
            FourierPoints(const PointsKey &key):PointsBaseType(key)
            {
            }

        private:
            /// Default constructor should not be called except by Create method.
            FourierPoints():PointsBaseType(NullPointsKey)
            {
            }

            /// Copy constructor should not be called.
            FourierPoints(const FourierPoints &points):PointsBaseType(points.m_pointsKey)
            {
            }
            
            void CalculatePoints();
            void CalculateWeights();
            void CalculateDerivMatrix();

        }; // class FourierPoints

        namespace
        {
            const bool FourierPointsInited1 = PointsManager().RegisterCreator(PointsKey(0, ePolyEvenlySpaced), FourierPoints::Create);
        }
    } // end of namespace
} // end of namespace 


#endif //FOURIERPOINTS_H
