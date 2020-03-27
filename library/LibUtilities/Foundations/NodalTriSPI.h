///////////////////////////////////////////////////////////////////////////////
//
// File NodalTriSPI.h
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
// Description: 2D Nodal Triangle SPI points
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALTRISPI_H
#define NODALTRISPI_H

#include <memory>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

namespace Nektar
{
namespace LibUtilities
{

class NodalTriSPI : public Points<NekDouble>
{
public:
    virtual ~NodalTriSPI()
    {
    }

    LIB_UTILITIES_EXPORT static std::shared_ptr<PointsBaseType> Create(
        const PointsKey &key);

    NodalTriSPI(const PointsKey &key) : PointsBaseType(key)
    {
    }

private:
    static bool initPointsManager[];

    NodalTriSPI() : PointsBaseType(NullPointsKey)
    {
    }

    void CalculatePoints();
    void CalculateWeights();
    void CalculateDerivMatrix();
};

}
}

#endif // NODALTRIELEC_H
