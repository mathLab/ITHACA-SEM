///////////////////////////////////////////////////////////////////////////////
//
// File NodalTriFekete.h
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
// Description: Header for 2D Fekete Points on a Tri (in Barocentric Coordinates)
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALTRIFEKETE_H
#define NODALTRIFEKETE_H

#include <math.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/Points.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

const int NodalTriFeketeAvailable = 16;
static const int NodalTriFeketeNPTS[] = {1,2,3,4,5,7,8,10,12,14,16,19,21,24,27,30};

namespace Nektar
{
    namespace LibUtilities 
    {
        class NodalTriFekete: public Points<double,2>
        {
        public:
            NodalTriFekete()
            {
                NEKERROR(efatal, "This constructor should not be called");
            }
 
            NodalTriFekete(const int &order): 
                Points<double,2>(order,eNodalTriFekete,eWildcard)
            {
            }

            virtual ~NodalTriFekete()
            {
            }


        protected:


        private:
            virtual void CalculateNumPoints();
            virtual void CalculatePoints();
            virtual void CalculateWeights();
            virtual void CalculateDerivMatrix();
        };  
        

    } // end of namespace
} // end of namespace 

#endif