///////////////////////////////////////////////////////////////////////////////
//
// File GaussPoints.h
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
// Description: Header file of GaussPoints Distributions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef GAUSSPOINTS_H
#define GAUSSPOINTS_H

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

        class GaussPoints: public PointsBaseType
        {
        public:
            virtual ~GaussPoints()
            {
            }            

            static boost::shared_ptr< PointsBaseType > Create(const PointsKey &key);
            static boost::shared_ptr< NekMatrix<DataType> > CreateMatrix(const PointsKey &key);

            const boost::shared_ptr<NekMatrix<double> > GetI(const PointsKey &pkey);
            const boost::shared_ptr<NekMatrix<double> > GetI(const double *x);
            const boost::shared_ptr<NekMatrix<double> > GetI(unsigned int numpoints, const double *x);

        protected:
            GaussPoints(const PointsKey &key):PointsBaseType(key)
            {
            }

        private:
            /// Default constructor should not be called except by Create method.
            GaussPoints():PointsBaseType(NullPointsKey)
            {
            }

            /// Copy constructor should not be called.
            GaussPoints(const GaussPoints &points):PointsBaseType(points.m_pointsKey)
            {
            }

            void CalculatePoints();
            void CalculateWeights();
            void CalculateDerivMatrix();
            void CalculateInterpMatrix(unsigned int npts, const double * xpoints, double * interp);
        }; // class GaussPoints

        namespace
        {
            const bool gaussInited1 = PointsManager().RegisterCreator(PointsKey(0, eGaussGaussLegendre), GaussPoints::Create);
            const bool gaussInited2 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMLegendre), GaussPoints::Create);
            const bool gaussInited3 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauPLegendre), GaussPoints::Create);
            const bool gaussInited4 = PointsManager().RegisterCreator(PointsKey(0, eGaussLobattoLegendre), GaussPoints::Create);
            const bool gaussInited5 = PointsManager().RegisterCreator(PointsKey(0, eGaussGaussChebyshev), GaussPoints::Create);
            const bool gaussInited6 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMChebyshev), GaussPoints::Create);
            const bool gaussInited7 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauPChebyshev), GaussPoints::Create);
            const bool gaussInited8 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta1), GaussPoints::Create);
            const bool gaussInited9 = PointsManager().RegisterCreator(PointsKey(0, eGaussGaussLegendre), GaussPoints::Create);
            const bool gaussInited10 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta2), GaussPoints::Create);
        }
    } // end of namespace
} // end of namespace 

#endif //GAUSSPOINTS_H
