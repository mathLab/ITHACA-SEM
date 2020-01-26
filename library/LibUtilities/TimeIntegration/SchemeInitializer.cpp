///////////////////////////////////////////////////////////////////////////////
//
// File SchemeInitializor.cpp
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
// Description: This file isused to add each of the Time Integration Schemes
//              to the NekFactory.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/TimeIntegration/AdamsBashforthTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/AdamsMoultonTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/BDFImplicitTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/EulerTimeIntegrationSchemes.h>

#include <LibUtilities/TimeIntegration/CNABTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/DIRKTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/IMEXGearTimeIntegrationScheme.h>
#include <LibUtilities/TimeIntegration/IMEXTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/IMEXdirkTimeIntegrationSchemes.h>
#include <LibUtilities/TimeIntegration/MCNABTimeIntegrationScheme.h>

#include <LibUtilities/TimeIntegration/RungeKuttaTimeIntegrationSchemes.h>

#include <LibUtilities/TimeIntegration/EulerExponentialTimeIntegrationSchemes.h>

namespace Nektar
{
namespace LibUtilities
{

// Register all the schemes with the Time Integration Scheme Facatory...
//
#define REGISTER(scheme)                                                       \
    std::string scheme##TimeIntegrationScheme::className =                     \
        GetTimeIntegrationSchemeFactory().RegisterCreatorFunction(             \
            #scheme, scheme##TimeIntegrationScheme::create)

// AdamsBashforthTimeIntegrationSchemes.h
REGISTER(AdamsBashforth);
REGISTER(AdamsBashforthOrder1);
REGISTER(AdamsBashforthOrder2);
REGISTER(AdamsBashforthOrder3);
REGISTER(AdamsBashforthOrder4);

// AdamsMoultonTimeIntegrationSchemes.h
REGISTER(AdamsMoulton);
REGISTER(AdamsMoultonOrder1);
REGISTER(AdamsMoultonOrder2);
REGISTER(AdamsMoultonOrder3);
REGISTER(AdamsMoultonOrder4);

// BDFImplicitTimeIntegrationSchemes.h
REGISTER(BDFImplicit);
REGISTER(BDFImplicitOrder1);
REGISTER(BDFImplicitOrder2);

// EulerTimeIntegrationSchemes.h
REGISTER(BackwardEuler);
REGISTER(ForwardEuler);

// EulerExponentialTimeIntegrationSchemes.h
REGISTER(EulerExponential);

// CNABTimeIntegrationScheme.h
REGISTER(CNAB);

// DIRKTimeIntegrationSchemes.h
REGISTER(DIRK);
REGISTER(DIRKOrder2);
REGISTER(DIRKOrder3);

// IMEXdirkTimeIntegrationSchemes.h
REGISTER(IMEXdirk);
REGISTER(IMEXdirk_1_1_1);
REGISTER(IMEXdirk_1_2_1);
REGISTER(IMEXdirk_1_2_2);
REGISTER(IMEXdirk_2_2_2);
REGISTER(IMEXdirk_2_3_2);
REGISTER(IMEXdirk_2_3_3);
REGISTER(IMEXdirk_3_4_3);
REGISTER(IMEXdirk_4_4_3);

// IMEXGearTimeIntegrationScheme.h
REGISTER(IMEXGear);

// IMEXTimeIntegrationSchemes.h
REGISTER(IMEX);
REGISTER(IMEXOrder1);
REGISTER(IMEXOrder2);
REGISTER(IMEXOrder3);
REGISTER(IMEXOrder4);

// MCNABTimeIntegrationScheme.h
REGISTER(MCNAB);

// RungeKuttaTimeIntegrationSchemes.h
REGISTER(RungeKutta);
REGISTER(RungeKutta2);
REGISTER(RungeKutta2_ImprovedEuler);
REGISTER(RungeKutta2_SSP);
REGISTER(RungeKutta3_SSP);
REGISTER(ClassicalRungeKutta4);
REGISTER(RungeKutta5);

} // end namespace LibUtilities
} // end namespace NekTar
