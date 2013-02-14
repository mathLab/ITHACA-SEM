///////////////////////////////////////////////////////////////////////////////
//
// File FitzhughNagumo.cpp
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
// Description: FitzhughNagumo phenomological cell model.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <CardiacEPSolver/CellModels/FitzhughNagumo.h>

namespace Nektar
{
    /**
     * Registers the class with the Factory.
     */
    std::string CellModelFitzHughNagumo::className
            = GetCellModelFactory().RegisterCreatorFunction(
                "FitzHughNagumo",
                CellModelFitzHughNagumo::create,
                "Phenomological model of squid nerve cell.");

    CellModelFitzHughNagumo::CellModelFitzHughNagumo(
                    const LibUtilities::SessionReaderSharedPtr& pSession,
                    const MultiRegions::ExpListSharedPtr& pField)
            : CellModel(pSession, pField)
    {
        pSession->LoadParameter("beta",          m_beta,         0.0);
        pSession->LoadParameter("epsilon",       m_epsilon,      1.0);

        m_uuu  = Array<OneD, NekDouble>(m_nq, 0.0);

        m_nvar = 2;
        m_concentrations.push_back(1);
    }


    void CellModelFitzHughNagumo::v_Update(
                    const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                    const NekDouble time)
    {
        NekDouble m_gamma = 0.5;

        // compute u^2: m_u = u*u
        Vmath::Vmul(m_nq, &inarray[0][0], 1, &inarray[0][0], 1, &m_uuu[0], 1);

        // compute u^3: m_u = u*u*u
        Vmath::Vmul(m_nq, &inarray[0][0], 1, &m_uuu[0], 1, &m_uuu[0], 1);

        // For u: (1/m_epsilon)*( u*-u*u*u/3 - v )
        // physfield = u - (1.0/3.0)*u*u*u
        Vmath::Svtvp(m_nq, (-1.0/3.0), &m_uuu[0], 1, &inarray[0][0], 1, &outarray[0][0], 1);

        Vmath::Vsub(m_nq, &inarray[1][0], 1, &outarray[0][0], 1, &outarray[0][0], 1);
        Vmath::Smul(m_nq, -1.0/m_epsilon, &outarray[0][0], 1, &outarray[0][0], 1);

        // For v: m_epsilon*( u + m_beta - m_gamma*v )
        Vmath::Svtvp(m_nq, -1.0*m_gamma, &inarray[1][0], 1, &inarray[0][0], 1, &outarray[1][0], 1);
        Vmath::Sadd(m_nq, m_beta, &outarray[1][0], 1, &outarray[1][0], 1);
        Vmath::Smul(m_nq, m_epsilon, &outarray[1][0], 1, &outarray[1][0], 1);
    }

    /**
     *
     */
    void CellModelFitzHughNagumo::v_PrintSummary(std::ostream &out)
    {
        out << "\tCell model      : FitzHugh-Nagumo" << std::endl;
        out << "\tBeta            : " << m_beta << std::endl;
    }


    /**
     *
     */
    void CellModelFitzHughNagumo::v_SetInitialConditions()
    {
        Vmath::Fill(m_nq, 0.0,        m_cellSol[0],  1);
        Vmath::Fill(m_nq, 0.0,        m_cellSol[1],  1);
    }

}
