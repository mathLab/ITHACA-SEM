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

// #include <SolverUtils/RiemannSolvers/RiemannSolver.h>
#include "../RoeSolver.h"


namespace Nektar
{
namespace RiemannTests
{
    BOOST_AUTO_TEST_CASE(RoeAlongXconstSolution)
    {
        // LibUtilities::SessionReaderSharedPtr dummySession;
        // SolverUtils::RiemannSolverSharedPtr riemannSolver;
        // std::string riemannName = "Roe";
        // riemannSolver = SolverUtils::GetRiemannSolverFactory()
                                    // .CreateInstance(riemName, m_session);
        auto riemannSolver = RoeSolver();
        // Setting up parameters for Riemann solver
        NekDouble gamma = 1.4;
        riemannSolver.SetParam("gamma", [&gamma]()
            -> NekDouble& {return gamma;});

        size_t spaceDim = 3;
        size_t npts = 1;

        // Set up locations of velocity vector.
        Array<OneD, Array<OneD, NekDouble>> vecLocs(1);
        vecLocs[0] = Array<OneD, NekDouble>(spaceDim);
        for (size_t i = 0; i < spaceDim; ++i)
        {
            vecLocs[0][i] = 1 + i;
        }
        riemannSolver.SetAuxVec("vecLocs", [&vecLocs]()
            -> const Array<OneD, const Array<OneD, NekDouble>>&
            {return vecLocs;});

        // setup normals
        Array<OneD, Array<OneD, NekDouble>> normals(spaceDim);
        for (size_t i = 0; i < spaceDim; ++i)
        {
           normals[i] = Array<OneD, NekDouble>(npts);
        }
        riemannSolver.SetVector("N", [&normals]()
            -> const Array<OneD, const Array<OneD, NekDouble>>&
            {return normals;});

        size_t nFields = spaceDim + 2;
        Array<OneD, Array<OneD, NekDouble>> fwd(nFields), bwd(nFields),
            flx(nFields), flxRef(nFields);
        for (size_t i = 0; i < nFields; ++i)
        {
            fwd[i] = Array<OneD, NekDouble>(npts);
            bwd[i] = Array<OneD, NekDouble>(npts);
            flx[i] = Array<OneD, NekDouble>(npts);
            flxRef[i] = Array<OneD, NekDouble>(npts);
        }

        // up to now it is "boiler plate" code to set up the test
        // below are the conditions for the test

        // density
        NekDouble rho = 0.9;
        fwd[0][0] = rho;
        bwd[0][0] = rho;
        // x-momentum
        NekDouble rhou = rho*1.0;
        fwd[1][0] = rhou;
        bwd[1][0] = rhou;
        // y-momentum
        NekDouble rhov = rho*2.0;
        fwd[2][0] = rhov;
        bwd[2][0] = rhov;
        // z-momentum
        NekDouble rhow = rho*3.0;
        fwd[3][0] = rhow;
        bwd[3][0] = rhow;
        // energy
        NekDouble p = 1.0;
        NekDouble rhoe = p / (gamma - 1.0);
        NekDouble E = rhoe + 0.5*(rhou*rhou + rhov*rhov + rhow*rhow) / rho;
        fwd[nFields-1][0] = E;
        bwd[nFields-1][0] = E;
        // set face normal along x
        normals[0][0] = 1.0;
        // Ref solution
        flxRef[0][0] = rhou;
        flxRef[1][0] = rhou*rhou/rho + p;
        flxRef[2][0] = rhou*rhov/rho;
        flxRef[3][0] = rhou*rhow/rho;
        flxRef[nFields-1][0] = (E+p)*rhou/rho;

        riemannSolver.Solve(spaceDim, fwd, bwd, flx);

        // check fluxes
        BOOST_CHECK_CLOSE(flxRef[0][0], flx[0][0], 1e-10);
        BOOST_CHECK_CLOSE(flxRef[1][0], flx[1][0], 1e-10);
        BOOST_CHECK_CLOSE(flxRef[2][0], flx[2][0], 1e-10);
        BOOST_CHECK_CLOSE(flxRef[3][0], flx[3][0], 1e-10);
        BOOST_CHECK_CLOSE(flxRef[4][0], flx[4][0], 1e-10);

    }


    BOOST_AUTO_TEST_CASE(RoeAlongYconstSolution)
    {
        // LibUtilities::SessionReaderSharedPtr dummySession;
        // SolverUtils::RiemannSolverSharedPtr riemannSolver;
        // std::string riemannName = "Roe";
        // riemannSolver = SolverUtils::GetRiemannSolverFactory()
                                    // .CreateInstance(riemName, m_session);
        auto riemannSolver = RoeSolver();
        // Setting up parameters for Riemann solver
        NekDouble gamma = 1.4;
        riemannSolver.SetParam("gamma", [&gamma]()
            -> NekDouble& {return gamma;});

        size_t spaceDim = 3;
        size_t npts = 1;

        // Set up locations of velocity vector.
        Array<OneD, Array<OneD, NekDouble>> vecLocs(1);
        vecLocs[0] = Array<OneD, NekDouble>(spaceDim);
        for (size_t i = 0; i < spaceDim; ++i)
        {
            vecLocs[0][i] = 1 + i;
        }
        riemannSolver.SetAuxVec("vecLocs", [&vecLocs]()
            -> const Array<OneD, const Array<OneD, NekDouble>>&
            {return vecLocs;});

        // setup normals
        Array<OneD, Array<OneD, NekDouble>> normals(spaceDim);
        for (size_t i = 0; i < spaceDim; ++i)
        {
           normals[i] = Array<OneD, NekDouble>(npts);
        }
        riemannSolver.SetVector("N", [&normals]()
            -> const Array<OneD, const Array<OneD, NekDouble>>&
            {return normals;});

        size_t nFields = spaceDim + 2;
        Array<OneD, Array<OneD, NekDouble>> fwd(nFields), bwd(nFields),
            flx(nFields), flxRef(nFields);
        for (size_t i = 0; i < nFields; ++i)
        {
            fwd[i] = Array<OneD, NekDouble>(npts);
            bwd[i] = Array<OneD, NekDouble>(npts);
            flx[i] = Array<OneD, NekDouble>(npts);
            flxRef[i] = Array<OneD, NekDouble>(npts);
        }

        // up to now it is "boiler plate" code to set up the test
        // below are the conditions for the test

        // density
        NekDouble rho = 0.9;
        fwd[0][0] = rho;
        bwd[0][0] = rho;
        // x-momentum
        NekDouble rhou = rho*1.0;
        fwd[1][0] = rhou;
        bwd[1][0] = rhou;
        // y-momentum
        NekDouble rhov = rho*2.0;
        fwd[2][0] = rhov;
        bwd[2][0] = rhov;
        // z-momentum
        NekDouble rhow = rho*3.0;
        fwd[3][0] = rhow;
        bwd[3][0] = rhow;
        // energy
        NekDouble p = 1.0;
        NekDouble rhoe = p / (gamma - 1.0);
        NekDouble E = rhoe + 0.5*(rhou*rhou + rhov*rhov + rhow*rhow) / rho;
        fwd[nFields-1][0] = E;
        bwd[nFields-1][0] = E;
        // set face normal along y
        normals[1][0] = 1.0;
        // Ref solution
        flxRef[0][0] = rhov;
        flxRef[1][0] = rhov*rhou/rho;
        flxRef[2][0] = rhov*rhov/rho + p;
        flxRef[3][0] = rhov*rhow/rho;
        flxRef[nFields-1][0] = (E+p)*rhov/rho;

        riemannSolver.Solve(spaceDim, fwd, bwd, flx);

        // check fluxes
        BOOST_CHECK_CLOSE(flxRef[0][0], flx[0][0], 1e-10);
        BOOST_CHECK_CLOSE(flxRef[1][0], flx[1][0], 1e-10);
        BOOST_CHECK_CLOSE(flxRef[2][0], flx[2][0], 1e-10);
        BOOST_CHECK_CLOSE(flxRef[3][0], flx[3][0], 1e-10);
        BOOST_CHECK_CLOSE(flxRef[4][0], flx[4][0], 1e-10);

    }

    BOOST_AUTO_TEST_CASE(RoeAlongZconstSolution)
    {
        // LibUtilities::SessionReaderSharedPtr dummySession;
        // SolverUtils::RiemannSolverSharedPtr riemannSolver;
        // std::string riemannName = "Roe";
        // riemannSolver = SolverUtils::GetRiemannSolverFactory()
                                    // .CreateInstance(riemName, m_session);
        auto riemannSolver = RoeSolver();
        // Setting up parameters for Riemann solver
        NekDouble gamma = 1.4;
        riemannSolver.SetParam("gamma", [&gamma]()
            -> NekDouble& {return gamma;});

        size_t spaceDim = 3;
        size_t npts = 1;

        // Set up locations of velocity vector.
        Array<OneD, Array<OneD, NekDouble>> vecLocs(1);
        vecLocs[0] = Array<OneD, NekDouble>(spaceDim);
        for (size_t i = 0; i < spaceDim; ++i)
        {
            vecLocs[0][i] = 1 + i;
        }
        riemannSolver.SetAuxVec("vecLocs", [&vecLocs]()
            -> const Array<OneD, const Array<OneD, NekDouble>>&
            {return vecLocs;});

        // setup normals
        Array<OneD, Array<OneD, NekDouble>> normals(spaceDim);
        for (size_t i = 0; i < spaceDim; ++i)
        {
           normals[i] = Array<OneD, NekDouble>(npts);
        }
        riemannSolver.SetVector("N", [&normals]()
            -> const Array<OneD, const Array<OneD, NekDouble>>&
            {return normals;});

        size_t nFields = spaceDim + 2;
        Array<OneD, Array<OneD, NekDouble>> fwd(nFields), bwd(nFields),
            flx(nFields), flxRef(nFields);
        for (size_t i = 0; i < nFields; ++i)
        {
            fwd[i] = Array<OneD, NekDouble>(npts);
            bwd[i] = Array<OneD, NekDouble>(npts);
            flx[i] = Array<OneD, NekDouble>(npts);
            flxRef[i] = Array<OneD, NekDouble>(npts);
        }

        // up to now it is "boiler plate" code to set up the test
        // below are the conditions for the test

        // density
        NekDouble rho = 0.9;
        fwd[0][0] = rho;
        bwd[0][0] = rho;
        // x-momentum
        NekDouble rhou = rho*1.0;
        fwd[1][0] = rhou;
        bwd[1][0] = rhou;
        // y-momentum
        NekDouble rhov = rho*2.0;
        fwd[2][0] = rhov;
        bwd[2][0] = rhov;
        // z-momentum
        NekDouble rhow = rho*3.0;
        fwd[3][0] = rhow;
        bwd[3][0] = rhow;
        // energy
        NekDouble p = 1.0;
        NekDouble rhoe = p / (gamma - 1.0);
        NekDouble E = rhoe + 0.5*(rhou*rhou + rhov*rhov + rhow*rhow) / rho;
        fwd[nFields-1][0] = E;
        bwd[nFields-1][0] = E;
        // set face normal along y
        normals[2][0] = 1.0;
        // Ref solution
        flxRef[0][0] = rhow;
        flxRef[1][0] = rhow*rhou/rho;
        flxRef[2][0] = rhow*rhov/rho;
        flxRef[3][0] = rhow*rhow/rho + p;
        flxRef[nFields-1][0] = (E+p)*rhow/rho;

        riemannSolver.Solve(spaceDim, fwd, bwd, flx);

        // check fluxes
        BOOST_CHECK_CLOSE(flxRef[0][0], flx[0][0], 1e-10);
        BOOST_CHECK_CLOSE(flxRef[1][0], flx[1][0], 1e-10);
        BOOST_CHECK_CLOSE(flxRef[2][0], flx[2][0], 1e-10);
        BOOST_CHECK_CLOSE(flxRef[3][0], flx[3][0], 1e-10);
        BOOST_CHECK_CLOSE(flxRef[4][0], flx[4][0], 1e-10);

    }



}
}