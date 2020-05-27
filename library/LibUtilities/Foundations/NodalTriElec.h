///////////////////////////////////////////////////////////////////////////////
//
// File NodalTriElec.h
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
// Description: Header file of 2D Nodal Triangle Fekete Points
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALTRIELEC_H
#define NODALTRIELEC_H

#include <memory>

#include <LibUtilities/Foundations/FoundationsFwd.hpp>
#include <LibUtilities/Foundations/NodalUtil.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar
{
    namespace LibUtilities
    {

        class NodalTriElec: public Points<NekDouble>
        {
        public:

            virtual ~NodalTriElec()
            {
            }

            LIB_UTILITIES_EXPORT static std::shared_ptr<PointsBaseType>
                Create(const PointsKey &key);

            const MatrixSharedPtrType GetI(const PointsKey &pkey)
            {
                ASSERTL0(pkey.GetPointsDim()==2,
                         "NodalTriElec Points can only interpolate to other 2d "
                         "point distributions");
                Array<OneD, const NekDouble> x, y;
                PointsManager()[pkey]->GetPoints(x, y);
                return GetI(x, y);
            }

            const MatrixSharedPtrType GetI(
                const Array<OneD, const NekDouble>& x,
                const Array<OneD, const NekDouble>& y)
            {
                size_t numpoints = x.size();
                Array<OneD, NekDouble> interp(GetTotNumPoints()*numpoints);
                CalculateInterpMatrix(x, y, interp);

                unsigned int np = GetTotNumPoints();
                NekDouble* d = interp.data();
                return MemoryManager<NekMatrix<NekDouble> >
                    ::AllocateSharedPtr(numpoints, np, d);
            }

            NodalTriElec(const PointsKey &key):PointsBaseType(key)
            {
            }

        private:
            static bool initPointsManager[];

            std::shared_ptr<NodalUtilTriangle> m_util;

            NodalTriElec():PointsBaseType(NullPointsKey)
            {
            }

            void CalculatePoints();
            void CalculateWeights();
            void CalculateDerivMatrix();
            void NodalPointReorder2d();

            void CalculateInterpMatrix(
                const Array<OneD, const NekDouble> &xia,
                const Array<OneD, const NekDouble> &yia,
                      Array<OneD, NekDouble> &interp);
        };
   } // end of namespace
} // end of namespace

#endif //NODALTRIELEC_H
