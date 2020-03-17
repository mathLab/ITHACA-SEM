///////////////////////////////////////////////////////////////////////////////
//
// File PolyEPoints.h
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
// Description: Header file of 1D Evenly-Spaced Point Definitions 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef POLYEPOINTS_H
#define POLYEPOINTS_H

#include <memory>

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    namespace LibUtilities 
    {



        class PolyEPoints: public Points<NekDouble>
        {
        public:
            typedef Points<NekDouble> PointsBaseType;

            virtual ~PolyEPoints()
            {
            }            

            LIB_UTILITIES_EXPORT static std::shared_ptr< PointsBaseType > Create(const PointsKey &key);

            LIB_UTILITIES_EXPORT const std::shared_ptr<NekMatrix<NekDouble> > GetI(const PointsKey &pkey);
            LIB_UTILITIES_EXPORT const std::shared_ptr<NekMatrix<NekDouble> > GetI(const Array<OneD, const NekDouble>& x);
            LIB_UTILITIES_EXPORT const std::shared_ptr<NekMatrix<NekDouble> > GetI(unsigned int numpoints, const Array<OneD, const NekDouble>& x);

            PolyEPoints(const PointsKey &key):PointsBaseType(key)
            {
            }

        private:
            static bool initPointsManager[];

            /// Default constructor should not be called except by Create method.
            PolyEPoints():PointsBaseType(NullPointsKey)
            {
            }

            /// Copy constructor should not be called.
            PolyEPoints(const PolyEPoints &points):PointsBaseType(points.m_pointsKey)
            {
            }
            
            void CalculatePoints();
            void CalculateWeights();
            void CalculateDerivMatrix();
            void CalculateInterpMatrix(unsigned int npts, const Array<OneD, const NekDouble>& xpoints, Array<OneD, NekDouble>&  interp);

            NekDouble LagrangeInterpolant(NekDouble x, int npts, const Array<OneD, const NekDouble>& xpts, const Array<OneD, const NekDouble>& funcvals);
            NekDouble LagrangePoly(NekDouble x, int pt, int npts, const Array<OneD, const NekDouble>& xpts);     
            NekDouble LagrangePolyDeriv(NekDouble x, int pt, int npts, const Array<OneD, const NekDouble>& xpts);

        }; // class PolyEPoints
    } // end of namespace
} // end of namespace 


#endif //POLYEPOINTS_H
