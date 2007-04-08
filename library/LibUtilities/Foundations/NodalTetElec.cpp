///////////////////////////////////////////////////////////////////////////////
//
// File NodalTetElec.cpp
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
// Description: 3D Nodal Tet Electrostatic Point Definitions
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <algorithm>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Foundations.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Foundations/NodalTetElec.h>
#include <LibUtilities/Foundations/NodalTetElecData.h>

namespace Nektar
{
    namespace LibUtilities 
    {
        void NodalTetElec::CalculatePoints()
        {
            // Allocate the storage for points
            Points<double>::CalculatePoints();
            ASSERTL0(false, "3D Point Expansion Not Implemented Yet");
        }

        void NodalTetElec::CalculateWeights()
        {
            // No weights computed
        }

        void NodalTetElec::CalculateDerivMatrix()
        {
            // No derivative matrix computed
        }

        boost::shared_ptr< Points<double> > NodalTetElec::Create(const PointsKey &key)
        {
            boost::shared_ptr< Points<double> > returnval(new NodalTetElec(key));
            returnval->Initialize();
            return returnval;
        }

        void NodalTetElec::NodalPointReorder3d()
        {
        }     

    } // end of namespace stdregion
} // end of namespace stdregion


/**
* %Log%
*/
