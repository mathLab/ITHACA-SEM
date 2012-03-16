///////////////////////////////////////////////////////////////////////////////
//
// File: LibUtilities.h
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_LIB_UTILITIES_H
#define NEKTAR_LIB_UTILITIES_LIB_UTILITIES_H

#ifdef NEKTAR_USE_PRECOMPILED_HEADERS

#include <LibUtilities/BasicConst/NektarUnivConsts.hpp>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/Polylib/Polylib.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/Memory/ThreadSpecificPool.hpp>

#include <LibUtilities/LinearAlgebra/Arpack.hpp>
#include <LibUtilities/LinearAlgebra/BandedMatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <LibUtilities/LinearAlgebra/SparseBlas.hpp>
#include <LibUtilities/LinearAlgebra/BlockMatrix.hpp>
#include <LibUtilities/LinearAlgebra/DiagonalMatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/FullMatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/LinearAlgebra/MatrixBase.hpp>
#include <LibUtilities/LinearAlgebra/MatrixOperations.hpp>
#include <LibUtilities/LinearAlgebra/MatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/MatrixStorageType.h>
#include <LibUtilities/LinearAlgebra/NekLinAlgAlgorithms.hpp>
#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrixMetadata.hpp>
#include <LibUtilities/LinearAlgebra/NekPoint.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/LinearAlgebra/NekVector.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorCommon.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorConstantSized.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorMetadata.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorObject.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorTypeTraits.hpp>
#include <LibUtilities/LinearAlgebra/NekVectorVariableSized.hpp>
#include <LibUtilities/LinearAlgebra/PointerWrapper.h>
#include <LibUtilities/LinearAlgebra/ScaledMatrix.hpp>
#include <LibUtilities/LinearAlgebra/StandardMatrix.hpp>
#include <LibUtilities/LinearAlgebra/SymmetricMatrixStoragePolicy.hpp>
#include <LibUtilities/LinearAlgebra/TransF77.hpp>
#include <LibUtilities/LinearAlgebra/TriangularMatrixStoragePolicy.hpp>
//#include <LibUtilities/LinearAlgebra/SparseStandardMatrix.hpp>

#include <LibUtilities/Interpreter/AnalyticExpressionEvaluator.hpp>

#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/FourierPoints.h>
#include <LibUtilities/Foundations/FourierSingleModePoints.h>
#include <LibUtilities/Foundations/GaussPoints.h>
#include <LibUtilities/Foundations/Graph.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/NodalTetElec.h>
#include <LibUtilities/Foundations/NodalTriElec.h>
#include <LibUtilities/Foundations/NodalTriFekete.h>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/PolyEPoints.h>


#include <LibUtilities/BasicUtils/ArrayPolicies.hpp>
#include <LibUtilities/BasicUtils/BoostUtil.hpp>
#include <LibUtilities/BasicUtils/Concepts.hpp>
#include <LibUtilities/BasicUtils/ConsistentObjectAccess.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/BasicUtils/OperatorGenerators.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/Kernel/kernel.h>


#endif //NEKTAR_USE_PRECOMPILED_HEADERS

#endif //NEKTAR_LIB_UTILITIES_LIB_UTILITIES_H


