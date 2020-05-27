///////////////////////////////////////////////////////////////////////////////
//
// File InterpCoeff.cpp
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
// Description: Definition of Interpolation methods for Coefficients 
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/InterpCoeff.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/Foundations/Basis.h>

namespace Nektar
{    
    namespace LibUtilities
    {     
        void InterpCoeff1D(const BasisKey                      &fbasis0,
                           const Array<OneD, const NekDouble>&  from,
                           const BasisKey                      &tbasis0,
                                 Array<OneD, NekDouble>        &to)
        {
            ASSERTL0(fbasis0.GetNumModes() == tbasis0.GetNumModes(),
                     "Number of modes must be the same for "
                     "interpolating coefficients");

            // Check to see if the same basis
            if (fbasis0.GetBasisType() == tbasis0.GetBasisType())
            {
                Vmath::Vcopy(fbasis0.GetNumModes(), from, 1, to, 1);
            }
            else
            {
                // interpolate
                DNekMatSharedPtr ftB = BasisManager()[fbasis0]->GetI(tbasis0);
	    
                NekVector<NekDouble> in (fbasis0.GetNumModes(), from, eWrapper);
                NekVector<NekDouble> out(tbasis0.GetNumModes(), to,   eWrapper);
                
                out = (*ftB)*in;
            }
        }

        void InterpCoeff2D(const BasisKey                      &fbasis0,
                           const BasisKey                      &fbasis1,
                           const Array<OneD, const NekDouble>&  from,
                           const BasisKey                      &tbasis0,
                           const BasisKey                      &tbasis1,
                                 Array<OneD, NekDouble>        &to)
        {
            InterpCoeff2D(fbasis0, fbasis1, from.data(),
                          tbasis0, tbasis1, to.  data());
        }

        void InterpCoeff2D(const BasisKey  &fbasis0,
                           const BasisKey  &fbasis1,
                           const NekDouble *from,
                           const BasisKey  &tbasis0,
                           const BasisKey  &tbasis1,
                                 NekDouble *to)
        {
            const int fnm0 = fbasis0.GetNumModes();
            const int fnm1 = fbasis1.GetNumModes();
            const int tnm0 = tbasis0.GetNumModes();
            const int tnm1 = tbasis1.GetNumModes();

            Array<OneD, NekDouble> wsp(tnm1 * fnm0);

            if (fbasis1.GetBasisType() == tbasis1.GetBasisType())
            {
                Vmath::Vcopy(fnm0*tnm1, from, 1, wsp.get(), 1);
            }
            else
            {
                // interpolate
                DNekMatSharedPtr ft1 = BasisManager()[fbasis1]->GetI(tbasis1);

                Blas::Dgemm('N', 'T', fnm0, tnm1, fnm1, 1.0, from, fnm0,
                            ft1->GetPtr().get(), tnm1, 0.0,  wsp.get(), fnm0);
            }

            if (fbasis0.GetBasisType() == tbasis0.GetBasisType())
            {
                Vmath::Vcopy(tnm0*tnm1, wsp.get(), 1, to, 1);
            }
            else
            {
                // interpolate
                DNekMatSharedPtr ft0 = BasisManager()[fbasis0]->GetI(tbasis0);

                Blas::Dgemm('N', 'N', tnm0, tnm1, fnm0, 1.0, ft0->GetPtr().get(),
                            tnm0, wsp.get(), fnm0, 0.0, to, tnm0);
            }
        }

        void InterpCoeff3D(const BasisKey                      &fbasis0,
                           const BasisKey                      &fbasis1,
                           const BasisKey                      &fbasis2,
                           const Array<OneD, const NekDouble>&  from,
                           const BasisKey                      &tbasis0,
                           const BasisKey                      &tbasis1,
                           const BasisKey                      &tbasis2,
                                 Array<OneD, NekDouble>        &to)
        {
            InterpCoeff3D(fbasis0, fbasis1, fbasis2, from.data(),
                          tbasis0, tbasis1, tbasis2, to.  data());
        }

        void InterpCoeff3D(const BasisKey  &fbasis0,
                           const BasisKey  &fbasis1,
                           const BasisKey  &fbasis2,
                           const NekDouble *from,
                           const BasisKey  &tbasis0,
                           const BasisKey  &tbasis1,
                           const BasisKey  &tbasis2,
                           NekDouble       *to)
        {
            const int fnm0 = fbasis0.GetNumModes();
            const int fnm1 = fbasis1.GetNumModes();
            const int fnm2 = fbasis2.GetNumModes();
            const int tnm0 = tbasis0.GetNumModes();
            const int tnm1 = tbasis1.GetNumModes();
            const int tnm2 = tbasis2.GetNumModes();

            Array<OneD, NekDouble> wsp1(tnm0 * tnm1 * fnm2);
            Array<OneD, NekDouble> wsp2(tnm0 * fnm1 * fnm2);

            DNekMatSharedPtr ft0 = BasisManager()[fbasis0]->GetI(tbasis0);
            DNekMatSharedPtr ft1 = BasisManager()[fbasis1]->GetI(tbasis1);
            DNekMatSharedPtr ft2 = BasisManager()[fbasis2]->GetI(tbasis2);

            Blas::Dgemm('N', 'N', tnm0, fnm1*fnm2, fnm0, 1.0,
                        ft0->GetPtr().get(), tnm0, from, fnm0, 0.0,
                        wsp2.get(), tnm0);

            for (int i = 0; i < fnm2; i++)
            {
                Blas::Dgemm('N', 'T', tnm0,  tnm1,  fnm1, 1.0,
                            wsp2.get()+i*tnm0*fnm1, tnm0, ft1->GetPtr().get(),
                            tnm1, 0.0, wsp1.get()+i*tnm0*tnm1, tnm0);
            }

            Blas::Dgemm('N', 'T', tnm0*tnm1, tnm2, fnm2, 1.0, wsp1.get(),
                        tnm0*tnm1, ft2->GetPtr().get(), tnm2,
                        0.0, to, tnm0*tnm1);
        }
    }
}
