///////////////////////////////////////////////////////////////////////////////
//
// File FoundationsFwd.hpp
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
// Description: Forward declarations of Foundation classes and typedefs
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FOUNDATIONS_FWD_H
#define FOUNDATIONS_FWD_H

#include <vector>
#include <memory>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/NekManager.hpp>

namespace Nektar
{
    namespace LibUtilities
    {

        class BLPoints;
        class Basis;
        class BasisKey;
        class FourierPoints;
        class FourierSingleModePoints;
        class GaussPoints;
        class GraphVertexObject;
        class GraphEdgeObject;
        class Graph;
        class NodalTetEvenlySpaced;
        class NodalTetElec;
        class NodalPrismEvenlySpaced;
        class NodalPrismElec;
        class NodalTriElec;
        class NodalTriEvenlySpaced;
        class NodalTriFekete;
        class PointsKey;
        class PolyEPoints;

        template<typename DataT>
        class Points;



        /// Name for a vector of BasisKeys.
        typedef std::vector< BasisKey > BasisKeyVector;
        typedef std::shared_ptr<Basis> BasisSharedPtr;
        typedef std::vector< BasisSharedPtr > BasisVector;

        typedef Points<NekDouble> PointsBaseType;
        typedef std::shared_ptr<Points<NekDouble> > PointsSharedPtr;
        typedef int GraphVertexID;

    } // end of namespace
} // end of namespace

#endif // FOUNDATIONS_FWD_H
