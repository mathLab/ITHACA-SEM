///////////////////////////////////////////////////////////////////////////////
//
// File: TimeIntegratorBase.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2018 Division of Applied Mathematics, Brown University (USA),
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
// Description: Time integration scheme wrappers; one class per time integration
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/TimeIntegration/TimeIntegratorBase.h>

#include <LibUtilities/TimeIntegration/AdamsBashforthOrder2TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/AdamsBashforthOrder3TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/AdamsMoultonOrder2TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/BackwardEulerTimeIntegrator.h>
#include <LibUtilities/TimeIntegration/BDFImplicitOrder1TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/BDFImplicitOrder2TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/ClassicalRungeKutta4TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/CNABTimeIntegrator.h>
#include <LibUtilities/TimeIntegration/DirkOrder2TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/DirkOrder3TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/ForwardEulerTimeIntegrator.h>
#include <LibUtilities/TimeIntegration/IMEXDirk_1_2_2TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/IMEXDirk_2_2_2TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/IMEXDirk_2_3_2TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/IMEXDirk_3_4_3TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/IMEXGearTimeIntegrator.h>
#include <LibUtilities/TimeIntegration/IMEXOrder1TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/IMEXOrder2TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/IMEXOrder3TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/MCNABTimeIntegrator.h>
#include <LibUtilities/TimeIntegration/RungeKutta2TimeIntegrator.h>
#include <LibUtilities/TimeIntegration/RungeKutta2_ImprovedEulerTimeIntegrator.h>
#include <LibUtilities/TimeIntegration/RungeKutta2_SSPTimeIntegrator.h>
#include <LibUtilities/TimeIntegration/RungeKutta3_SSPTimeIntegrator.h>

///////////////////////////////////////////////////////////////////////////////

using std::string;

///////////////////////////////////////////////////////////////////////////////

namespace Nektar {
  namespace LibUtilities {

    TimeIntegratorFactory & GetTimeIntegratorFactory()
    {
        static TimeIntegratorFactory instance;
        return instance;
    }

#   define GTIF GetTimeIntegratorFactory
    string AdamsBashforthOrder2TimeIntegrator::className = GTIF().RegisterCreatorFunction( "AdamsBashforthOrder2", AdamsBashforthOrder2TimeIntegrator::create );
    string AdamsBashforthOrder3TimeIntegrator::className = GTIF().RegisterCreatorFunction( "AdamsBashforthOrder3", AdamsBashforthOrder3TimeIntegrator::create );
    string AdamsMoultonOrder2TimeIntegrator::className   = GTIF().RegisterCreatorFunction( "AdamsMoultonOrder2",   AdamsMoultonOrder2TimeIntegrator::create );
    string BackwardEulerTimeIntegrator::className        = GTIF().RegisterCreatorFunction( "BackwardEuler",        BackwardEulerTimeIntegrator::create );
    string BDFImplicitOrder1TimeIntegrator::className    = GTIF().RegisterCreatorFunction( "BDFImplicitOrder1",    BDFImplicitOrder1TimeIntegrator::create );
    string BDFImplicitOrder2TimeIntegrator::className    = GTIF().RegisterCreatorFunction( "BDFImplicitOrder2",    BDFImplicitOrder2TimeIntegrator::create );
    string ClassicalRungeKutta4TimeIntegrator::className = GTIF().RegisterCreatorFunction( "ClassicalRungeKutta4", ClassicalRungeKutta4TimeIntegrator::create );
    string CNABTimeIntegrator::className                 = GTIF().RegisterCreatorFunction( "CNAB",                 CNABTimeIntegrator::create );
    string DirkOrder2TimeIntegrator::className           = GTIF().RegisterCreatorFunction( "DIRKOrder2",           DirkOrder2TimeIntegrator::create );
    string DirkOrder3TimeIntegrator::className           = GTIF().RegisterCreatorFunction( "DIRKOrder3",           DirkOrder3TimeIntegrator::create );
    string ForwardEulerTimeIntegrator::className         = GTIF().RegisterCreatorFunction( "ForwardEuler",         ForwardEulerTimeIntegrator::create );
    string IMEXDirk_1_2_2TimeIntegrator::className       = GTIF().RegisterCreatorFunction( "IMEXdirk_1_2_2",       IMEXDirk_1_2_2TimeIntegrator::create );
    string IMEXDirk_2_2_2TimeIntegrator::className       = GTIF().RegisterCreatorFunction( "IMEXdirk_2_2_2",       IMEXDirk_2_2_2TimeIntegrator::create );
    string IMEXDirk_2_3_2TimeIntegrator::className       = GTIF().RegisterCreatorFunction( "IMEXdirk_2_3_2",       IMEXDirk_2_3_2TimeIntegrator::create );
    string IMEXDirk_3_4_3TimeIntegrator::className       = GTIF().RegisterCreatorFunction( "IMEXdirk_3_4_3",       IMEXDirk_3_4_3TimeIntegrator::create );
    string IMEXGearTimeIntegrator::className             = GTIF().RegisterCreatorFunction( "IMEXGear",             IMEXGearTimeIntegrator::create );
    string IMEXOrder1TimeIntegrator::className           = GTIF().RegisterCreatorFunction( "IMEXOrder1",           IMEXOrder1TimeIntegrator::create );
    string IMEXOrder2TimeIntegrator::className           = GTIF().RegisterCreatorFunction( "IMEXOrder2",           IMEXOrder2TimeIntegrator::create );
    string IMEXOrder3TimeIntegrator::className           = GTIF().RegisterCreatorFunction( "IMEXOrder3",           IMEXOrder3TimeIntegrator::create );
    string MCNABTimeIntegrator::className                = GTIF().RegisterCreatorFunction( "MCNAB",                MCNABTimeIntegrator::create );
    string RungeKutta2TimeIntegrator::className          = GTIF().RegisterCreatorFunction( "RungeKutta2",          RungeKutta2TimeIntegrator::create );
    string RungeKutta2_ImprovedEulerTimeIntegrator::className = GTIF().RegisterCreatorFunction( "RungeKutta2_ImprovedEuler", RungeKutta2_ImprovedEulerTimeIntegrator::create );
    string RungeKutta2_SSPTimeIntegrator::className      = GTIF().RegisterCreatorFunction( "RungeKutta2_SSP", RungeKutta2_SSPTimeIntegrator::create );
    string RungeKutta3_SSPTimeIntegrator::className      = GTIF().RegisterCreatorFunction( "RungeKutta3_SSP", RungeKutta3_SSPTimeIntegrator::create );

  } // end namespace LibUtilities
} // end namespace Nektar
