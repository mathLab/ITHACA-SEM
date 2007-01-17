///////////////////////////////////////////////////////////////////////////////
//
// File Points1D.h
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
// Description: Header file of 1D Points definition 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef POINTS1D_H
#define POINTS1D_H

#include <math.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/Points.hpp>

namespace Nektar
{
    namespace LibUtilities 
    {
        class GaussPolyPoints: public Points<double>
        {
        public:
        GaussPolyPoints::GaussPolyPoints(const PointsKey &key, 
            const double alpha, const double beta);

            virtual ~GaussPolyPoints()
            {
            }

            // Not sure what to return here for the create function.
            // Since this function will be registered in a manager with 
            // other points types, it needs to return a base type.
            // I am not sure if Points<> will work because of its type
            // obtained from the template parameters.  That would make
            // it so the given manager can only hold, say, Points<double>.
            static boost::shared_ptr< Points<double> > Create(const PointsKey &key);

        protected:
            double m_alpha;
            double m_beta;

        private:
            void CalculatePoints();
            void CalculateWeights();
            void CalculateDerivMatrix();
        };  
    } // end of namespace
} // end of namespace 

#endif //POINTS1D_H
