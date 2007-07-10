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
// Description: Header file of 2D Nodal Triangle Fekete Points 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NODALTRIFEKETE_H
#define NODALTRIFEKETE_H

#include <math.h>
#include <boost/shared_ptr.hpp>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/ManagerAccess.h>

namespace Nektar
{
    namespace LibUtilities 
    {
 
        class NodalTriFekete: public Points<NekDouble>
        {
        public:
            typedef Points<NekDouble> PointsBaseType;

            virtual ~NodalTriFekete()
            {
            }

            static boost::shared_ptr<PointsBaseType> Create(const PointsKey &key);

            const boost::shared_ptr<NekMatrix<NekDouble> > GetI(const PointsKey &pkey)
            {
                ASSERTL0(false, "NodalTriFekete Method not implemented");
                boost::shared_ptr<NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr());

                return returnval;
            }

            const boost::shared_ptr<NekMatrix<NekDouble> > GetI(const ConstArray<OneD, NekDouble>& x)
            {
                ASSERTL0(false, "NodalTriFekete Method not implemented");

                boost::shared_ptr<NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr());

                return returnval;
            }

            const boost::shared_ptr<NekMatrix<NekDouble> > GetI(unsigned int numpoints, const ConstArray<OneD, NekDouble>& x)
            {
                ASSERTL0(false, "NodalTriFekete Method not implemented");

                boost::shared_ptr<NekMatrix<NekDouble> > returnval(MemoryManager<NekMatrix<NekDouble> >::AllocateSharedPtr());

                return returnval;
            }

            NodalTriFekete(const PointsKey &key):PointsBaseType(key)
            {
            }

        private:
            NodalTriFekete():PointsBaseType(NullPointsKey)
            {
            }

            void CalculatePoints();
            void CalculateWeights();
            void CalculateDerivMatrix();
            void NodalPointReorder2d();

           NekDouble LagrangePoly(NekDouble x, int pt, int npts, const ConstArray<OneD, NekDouble>& xpts);

        }; // end of NodalTriFekete
   } // end of namespace
} // end of namespace 

#endif //NODALTRIFEKETE_H
