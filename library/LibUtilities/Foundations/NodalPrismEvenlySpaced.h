///////////////////////////////////////////////////////////////////////////////
//
// File NodalPrismEvenlySpaced.h
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
// Description: Header file of 3D Nodal Prism Evenly Spaced Points
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALPRISMEVENLYSPACED_H
#define NODALPRISMEVENLYSPACED_H

#include <memory>

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Nektar
{
    namespace LibUtilities
    {
        class NodalPrismEvenlySpaced: public Points<NekDouble>
        {
        public:
            virtual ~NodalPrismEvenlySpaced()
            {

            }

            NodalPrismEvenlySpaced(const PointsKey &key):PointsBaseType(key)
            {

            }

            LIB_UTILITIES_EXPORT static std::shared_ptr<PointsBaseType>
                Create(const PointsKey &key);

            const MatrixSharedPtrType GetI(const PointsKey &pkey)
            {
                ASSERTL0(pkey.GetPointsDim() == 3,
                         "NodalPrismEvenlySpaced Points can only interp to "
                         "other 3d point distributions");
                Array<OneD, const NekDouble> x, y, z;
                PointsManager()[pkey]->GetPoints(x, y, z);
                return GetI(x, y, z);
            }

            const MatrixSharedPtrType GetI(
                const Array<OneD, const NekDouble> &x,
                const Array<OneD, const NekDouble> &y,
                const Array<OneD, const NekDouble> &z)
            {
                size_t       numpoints = x.size();
                unsigned int np        = GetTotNumPoints();

                Array<OneD, NekDouble> interp(GetTotNumPoints()*numpoints);
                CalculateInterpMatrix(x, y, z, interp);

                NekDouble* d = interp.data();
                return MemoryManager<NekMatrix<NekDouble> >
                    ::AllocateSharedPtr(numpoints, np, d);
            }

        private:
            static bool initPointsManager[];

            std::shared_ptr<NodalUtilPrism> m_util;

            /// Default constructor should not be called except by Create matrix
            NodalPrismEvenlySpaced():PointsBaseType(NullPointsKey)
            {
            }

            void CalculatePoints();
            void CalculateWeights();
            void CalculateDerivMatrix();
            void NodalPointReorder3d();
            void CalculateInterpMatrix(
                const Array<OneD, const NekDouble> &xi,
                const Array<OneD, const NekDouble> &yi,
                const Array<OneD, const NekDouble> &zi,
                      Array<OneD,       NekDouble> &interp);
        }; // end of NodalPrismEvenlySpaced
   } // end of namespace
} // end of namespace

#endif //NODALPRISMEVENLYSPACED_H
