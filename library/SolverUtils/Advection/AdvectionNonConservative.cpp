///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionNonConservative.cpp
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
// Description: Non-conservative advection class.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Advection/AdvectionNonConservative.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string AdvectionNonConservative::type = GetAdvectionFactory().
            RegisterCreatorFunction("NonConservative",
                                    AdvectionNonConservative::create);

        AdvectionNonConservative::AdvectionNonConservative()
        {

        }

        /**
         * @brief Initialise AdvectionNonConservative objects and store them
         * before starting the time-stepping.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void AdvectionNonConservative::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            Advection::v_InitObject(pSession, pFields);
        }

        void AdvectionNonConservative::v_Advect(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const NekDouble                                   &time,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            boost::ignore_unused(time, pFwd, pBwd);

            int nDim       = advVel.size();
            int nPointsTot = fields[0]->GetNpoints();
            Array<OneD, NekDouble> grad0,grad1,grad2;

            grad0 = Array<OneD, NekDouble> (nPointsTot);

            if (nDim > 1)
            {
                grad1 = Array<OneD,NekDouble>(nPointsTot);
            }

            if (nDim > 2)
            {
                grad2 = Array<OneD,NekDouble>(nPointsTot);
            }


            for (int i = 0; i < nConvectiveFields; ++i)
            {
                // Evaluate V \cdot Grad(u)
                switch(nDim)
                {
                    case 1:
                        fields[0]->PhysDeriv(inarray[i], grad0);

                        Vmath::Vmul(nPointsTot,
                                    grad0,          1,
                                    advVel[0],      1,
                                    outarray[i],    1);
                        break;
                    case 2:
                        fields[0]->PhysDeriv(inarray[i], grad0, grad1);


                        // Calculate advection terms
                        Vmath::Vmul (nPointsTot,
                                     grad0, 1,
                                     advVel[0], 1,
                                     outarray[i], 1);

                        Vmath::Vvtvp(nPointsTot,
                                     grad1, 1,
                                     advVel[1], 1,
                                     outarray[i], 1,
                                     outarray[i], 1);

                        break;
                      case 3:
                        fields[0]->PhysDeriv(inarray[i], grad0, grad1, grad2);

                        // Calculate advection terms
                        Vmath::Vmul (nPointsTot,
                                     grad0, 1,
                                     advVel[0], 1,
                                     outarray[i], 1);

                        Vmath::Vvtvp(nPointsTot,
                                     grad1, 1,
                                     advVel[1], 1,
                                     outarray[i], 1,
                                     outarray[i], 1);

                        Vmath::Vvtvp(nPointsTot,
                                     grad2, 1,
                                     advVel[2], 1,
                                     outarray[i], 1,
                                     outarray[i], 1);
                        break;
                    default:
                        ASSERTL0(false,"dimension unknown");
                }
            }
        }
    }
}
