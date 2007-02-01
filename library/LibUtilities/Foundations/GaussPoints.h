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
#include <boost/bind.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/Points.h>

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

            static boost::shared_ptr< PointsBaseType > Create(const PointsKey &pkey);
            static boost::shared_ptr< NekMatrix<double> > CreateMatrix(const PointsKey &ownpkey, const PointsKey &pkey);

            const boost::shared_ptr<NekMatrix<double> > GetI(const PointsKey &pkey);
            const boost::shared_ptr<NekMatrix<double> > GetI(const double *x);
            const boost::shared_ptr<NekMatrix<double> > GetI(unsigned int numpoints, const double *x);

        protected:
            GaussPoints(const PointsKey &pkey):PointsBaseType(pkey)
            {
            //eGaussGaussLegendre,     //!< 1D Gauss-Gauss-Legendre quadrature points
            //eGaussRadauMLegendre,    //!< 1D Gauss-Radau-Legendre quadrature points, pinned at x=-1
            //eGaussRadauPLegendre,    //!< 1D Gauss-Radau-Legendre quadrature points, pinned at x=1
            //eGaussLobattoLegendre,   //!< 1D Gauss-Lobatto-Legendre quadrature points
            //eGaussGaussChebyshev,    //!< 1D Gauss-Gauss-Chebyshev quadrature points
            //eGaussRadauMChebyshev,   //!< 1D Gauss-Radau-Chebyshev quadrature points, pinned at x=-1
            //eGaussRadauPChebyshev,   //!< 1D Gauss-Radau-Chebyshev quadrature points, pinned at x=1
            //eGaussLobattoChebyshev,  //!< 1D Gauss-Lobatto-Legendre quadrature points
            //eGaussRadauMAlpha0Beta1, //!< Gauss Radau pinned at x=-1, \f$ \alpha =    0, \beta =    1 \f$
            //eGaussRadauMAlpha0Beta2, //!< Gauss Radau pinned at x=-1, \f$ \alpha =    0, \beta =    2 \f$
            //ePolyEvenlySpaced,       //!< 1D Evenly-spaced points using Lagrange polynomial

                //m_InterpManager.RegisterCreator(PointsKey(0, eGaussGaussLegendre),
                //    boost::bind(&GaussPoints::CreateMatrix, this, _1, _2)(pkey));
            }

        private:
            /// These should not be called.  All creation is done using the constructor
            /// requiring the key, declared above.
            GaussPoints();
            GaussPoints(const GaussPoints &points);

            void CalculatePoints();
            void CalculateWeights();
            void CalculateDerivMatrix();
            void CalculateInterpMatrix(unsigned int npts, const double * xpoints, double * interp);
        }; // class GaussPoints
    } // end of namespace
} // end of namespace 

#endif //GAUSSPOINTS_H
