///////////////////////////////////////////////////////////////////////////////
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
/// The above copyright notice and this permission notice shall be included
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

#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>


#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

namespace Nektar
{
namespace RiemannTests
{


    LibUtilities::SessionReaderSharedPtr dummySession;

    SolverUtils::RiemannSolverSharedPtr riemannSolver;
    std::string riemannName = "Roe";
    riemannSolver = SolverUtils::GetRiemannSolverFactory()
                                .CreateInstance(riemannName);

    // // Setting up parameters for advection operator Riemann solver
    // NekDouble gamma = 1.4;
    // riemannSolver->SetParam ("gamma",   gamma,   this);
    // riemannSolver->SetAuxVec("vecLocs", &CompressibleFlowSystem::GetVecLocs, this);
    // riemannSolver->SetVector("N",       &CompressibleFlowSystem::GetNormals, this);

    BOOST_AUTO_TEST_CASE(Roe)
    {
        // size_t spaceDim = 1;
        // size_t npts = 10;
        // Array<OneD, Array<OneD, NekDouble> > fwd(spaceDim), bwd(spaceDim),
        //     flx(spaceDim);
        // for (size_t i = 0, i < spaceDim, ++i)
        // {
        //     fwd[i] = Array<OneD, NekDouble>(npts);
        //     bwd[i] = Array<OneD, NekDouble>(npts);
        //     flx[i] = Array<OneD, NekDouble>(npts);

        // }
        // m_riemann->Solve(spaceDim, fwd, bwd, flx);

        BOOST_CHECK_CLOSE(1., 1., 1e-10);
    }




}
}