///////////////////////////////////////////////////////////////////////////////
//
// File: APEUpwindSolver.cpp
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
// Description: Upwind Riemann solver for the APE equations.
//
///////////////////////////////////////////////////////////////////////////////

#include <APESolver/RiemannSolvers/APEUpwindSolver.h>

using namespace std;

namespace Nektar
{

std::string APEUpwindSolver::solverName = SolverUtils::GetRiemannSolverFactory().
        RegisterCreatorFunction("APEUpwind", APEUpwindSolver::create, "APE Upwind solver");

/**
*
*/

/**
*
*/
APEUpwindSolver::APEUpwindSolver() :
    RiemannSolver()
{
    // we need out own rotation logic
    m_requiresRotation = false;
}

/**
*
*/
void APEUpwindSolver::v_Solve(
        const int                                         nDim,
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux)
{

    ASSERTL1(CheckVectors("N"),         "N not defined.");
    ASSERTL1(CheckAuxVec ("vecLocs"),   "vecLoc not defined.");
    ASSERTL1(CheckVectors("basefield"), "basefield not defined.");
    const Array<OneD, const Array<OneD, NekDouble> > normals =
            m_vectors["N"]();
    const Array<OneD, const Array<OneD, NekDouble> > vecLocs =
            m_auxVec["vecLocs"]();
    const Array<OneD, const Array<OneD, NekDouble> > basefield =
            m_vectors["basefield"]();

    int nFields = Fwd   .num_elements();
    int nPts    = Fwd[0].num_elements();

    // rotate and store basefield
    m_rotBasefield = Array<OneD, Array<OneD, NekDouble> > (nDim+1);
    for (int i = 0; i < nDim + 1; i++)
    {
        m_rotBasefield[i] = Array<OneD, NekDouble>(nPts);
    }
    Array<OneD, Array<OneD, NekDouble> > baseVecLocs(1);
    baseVecLocs[0] = Array<OneD, NekDouble>(nDim);
    for (int i = 0; i < nDim; ++i)
    {
        baseVecLocs[0][i] = 1+i;
    }

    rotateToNormal(basefield, normals, baseVecLocs, m_rotBasefield);

    if (m_rotStorage[0].num_elements()    != nFields ||
            m_rotStorage[0][0].num_elements() != nPts)
    {
        for (int i = 0; i < 3; ++i)
        {
            m_rotStorage[i] =
                    Array<OneD, Array<OneD, NekDouble> >(nFields);
            for (int j = 0; j < nFields; ++j)
            {
                m_rotStorage[i][j] = Array<OneD, NekDouble>(nPts);
            }
        }
    }

    rotateToNormal(Fwd, normals, vecLocs, m_rotStorage[0]);
    rotateToNormal(Bwd, normals, vecLocs, m_rotStorage[1]);
    Solve1D(m_rotStorage[0], m_rotStorage[1], m_rotStorage[2]);
    rotateFromNormal(m_rotStorage[2], normals, vecLocs, flux);
}



/**
*
*/
void APEUpwindSolver::Solve1D(
        const Array<OneD, const Array<OneD, NekDouble> > &Fwd,
        const Array<OneD, const Array<OneD, NekDouble> > &Bwd,
              Array<OneD,       Array<OneD, NekDouble> > &flux)
{
    // fetch params
    ASSERTL1(CheckParams("Gamma"), "Gamma not defined.");
    ASSERTL1(CheckParams("Rho"), "Rho not defined.");
    ASSERTL1(CheckVectors("N"), "N not defined.");
    const NekDouble &gamma= m_params["Gamma"]();
    const NekDouble &rho = m_params["Rho"]();
    const Array<OneD, const Array<OneD, NekDouble> > normals = m_vectors["N"]();
    const Array<OneD, const Array<OneD, NekDouble> > basefield = m_rotBasefield;

    int dim  = normals.num_elements();
    int nvar = dim +1;
    ASSERTL1(Fwd.num_elements() == nvar, "Fwd malformed.");
    ASSERTL1(Bwd.num_elements() == nvar, "Bwd malformed.");

    int nTracePts = Fwd[0].num_elements();

    Array<OneD, Array<OneD, NekDouble> > upphysfield(2);
    for (int i = 0; i < 2; i++)
    {
        upphysfield[i] = Array<OneD, NekDouble>(nTracePts);
    }

    for (int i = 0; i < nTracePts; i++)
    {
        // assign variables
        NekDouble pL = Fwd[0][i];
        NekDouble uL = Fwd[1][i];

        NekDouble pR = Bwd[0][i];
        NekDouble uR = Bwd[1][i];

        NekDouble p0 = basefield[0][i];
        NekDouble u0 = basefield[1][i];

        Array<OneD, NekDouble> characteristic(4);
        Array<OneD, NekDouble> W(2);
        Array<OneD, NekDouble> lambda(2);

        // compute the wave speeds
        lambda[0]=u0 + sqrt(p0*gamma*rho)/rho;
        lambda[1]=u0 - sqrt(p0*gamma*rho)/rho;

        // calculate the caracteristic variables
        //left characteristics
        characteristic[0] = pL/2 + uL*sqrt(p0*gamma*rho)/2;
        characteristic[1] = pL/2 - uL*sqrt(p0*gamma*rho)/2;
        //right characteristics
        characteristic[2] = pR/2 + uR*sqrt(p0*gamma*rho)/2;
        characteristic[3] = pR/2 - uR*sqrt(p0*gamma*rho)/2;

        //take left or right value of characteristic variable
        for (int j = 0; j < 2; j++)
        {
            if (lambda[j]>=0)
            {
                W[j]=characteristic[j];
            }
            if(lambda[j]<0)
            {
                W[j]=characteristic[j+2];
            }
        }

        //calculate conservative variables from characteristics
        upphysfield[0][i] = W[0]+W[1];
        // do not divide by zero
        if (p0*gamma*rho != 0)
        {
            upphysfield[1][i] = (W[0]-W[1])/sqrt(p0*gamma*rho);
        }
        else
        {
            upphysfield[1][i] = 0.0;
        }
    }

    // compute the fluxes

    // flux[0][i] = gamma*p0 * upphysfield[1] + u0*upphysfield[0]
    Vmath::Smul(nTracePts, gamma, basefield[0], 1, flux[0], 1);
    Vmath::Vmul(nTracePts, upphysfield[1], 1, flux[0], 1, flux[0], 1);
    Vmath::Vvtvp(nTracePts, basefield[1], 1, upphysfield[0], 1, flux[0], 1, flux[0], 1);

    // flux[1][i] = upphysfield[0]/rho + u0*upphysfield[1];
    Vmath::Smul(nTracePts, 1/rho, upphysfield[0], 1, flux[1], 1);
    Vmath::Vvtvp(nTracePts, basefield[1], 1, upphysfield[1], 1, flux[1], 1, flux[1], 1);
    for (int j = 2; j < nvar; j++)
    {
        // flux[1][i] = upphysfield[0]/rho + u0*upphysfield[1] + v0 * vL + w0 + wL;
        Vmath::Vvtvp(nTracePts, basefield[j], 1, Fwd[j], 1, flux[1], 1, flux[1], 1);

        // flux[2/3][i] = 0.0;
        Vmath::Zero(nTracePts, flux[j], 1);
    }
}
}
