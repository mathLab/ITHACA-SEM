///////////////////////////////////////////////////////////////////////////////
//
// File AlievPanfilov.cpp
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
// Description: Aliev-Panfilov phenomological cell model.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>

#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <CardiacEPSolver/CellModels/AlievPanfilov.h>

namespace Nektar
{
    /**
     * Registers the class with the Factory.
     */
    std::string CellModelAlievPanfilov::className
            = GetCellModelFactory().RegisterCreatorFunction(
                "AlievPanfilov",
                CellModelAlievPanfilov::create,
                "Phenomological model of canine cardiac electrophysiology.");

    CellModelAlievPanfilov::CellModelAlievPanfilov(
                    const LibUtilities::SessionReaderSharedPtr& pSession,
                    const MultiRegions::ExpListSharedPtr& pField)
            : CellModel(pSession, pField)
    {
        pSession->LoadParameter("k",          m_k,         0.0);
        pSession->LoadParameter("a",          m_a,         0.0);
        pSession->LoadParameter("mu1",        m_mu1,       0.0);
        pSession->LoadParameter("mu2",        m_mu2,       0.0);
        pSession->LoadParameter("eps",        m_eps,       0.0);

        m_uu   = Array<OneD, NekDouble>(m_nq, 0.0);
        m_uuu  = Array<OneD, NekDouble>(m_nq, 0.0);
        m_tmp1 = Array<OneD, NekDouble>(m_nq, 0.0);
        m_tmp2 = Array<OneD, NekDouble>(m_nq, 0.0);

        m_nvar = 2;
        m_concentrations.push_back(1);
    }


    void CellModelAlievPanfilov::v_Update(
                    const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                          Array<OneD,        Array<OneD, NekDouble> >&outarray,
                    const NekDouble time)
    {
        // inarray[0] holds initial physical u values throughout
        // inarray[1] holds initial physical v values throughout

        // compute u^2: m_u = u*u
        Vmath::Vmul(m_nq, &inarray[0][0], 1, &inarray[0][0], 1, &m_uu[0], 1);

        // compute u^3: m_u = u*u*u
        Vmath::Vmul(m_nq, &inarray[0][0], 1, &m_uu[0], 1, &m_uuu[0], 1);

        // --------------------------------------
        // Compute reaction term f(u,v)
        // --------------------------------------
//        if (m_spatialParameters->Exists("a"))
//        {
//          Vmath::Vmul(m_nq,  &m_spatialParameters->GetData("a")->GetPhys()[0], 1,
//                           &inarray[0][0], 1, &m_tmp1[0], 1);
//
//          Vmath::Vvtvm(m_nq, &m_spatialParameters->GetData("a")->GetPhys()[0], 1,
//                           &m_uu[0], 1, &m_tmp1[0], 1, &m_tmp1[0], 1);
//
//          Vmath::Svtvm(m_nq, -1.0, &m_uu[0], 1, &m_tmp1[0], 1, &m_tmp1[0], 1);
//        }
//        else
//        {
          // Ru = au
          Vmath::Smul(m_nq, m_a, &inarray[0][0], 1, &m_tmp1[0], 1);
          // Ru = (-1-a)u*u + au
          Vmath::Svtvp(m_nq, (-1.0-m_a), &m_uu[0], 1, &m_tmp1[0], 1,
                                       &m_tmp1[0], 1);
//        }
        // Ru = u*u*u - (1+a)u*u + au
        Vmath::Vadd(m_nq, &m_uuu[0], 1, &m_tmp1[0], 1, &m_tmp1[0], 1);
        // Ru = k(u*u*u - (1+a)u*u + au)
//        if (m_spatialParameters->Exists("k"))
//        {
//          Vmath::Vmul(m_nq, &m_spatialParameters->GetData("k")->GetPhys()[0], 1,
//                          &m_tmp1[0], 1, &m_tmp1[0], 1);
//        }
//        else
//        {
          Vmath::Smul(m_nq, m_k, &m_tmp1[0], 1, &m_tmp1[0], 1);
//        }

        // Ru = k(u*u*u - (1+a)u*u + au) + I_stim
        Vmath::Vadd(m_nq, &outarray[0][0], 1, &m_tmp1[0], 1, &outarray[0][0], 1);

        // Ru = k(u*u*u - (1+a)u*u + au) + uv + I_stim
        Vmath::Vvtvp(m_nq, &inarray[0][0], 1, &inarray[1][0], 1, &m_tmp1[0], 1,
                         &outarray[0][0], 1);
        // Ru = -k(u*u*u - (1+a)u*u + au) - uv - I_stim
        Vmath::Neg(m_nq, &outarray[0][0], 1);


        // --------------------------------------
        // Compute reaction term g(u,v)
        // --------------------------------------
        // tmp2 = mu2 + u
        Vmath::Sadd(m_nq, m_mu2, &inarray[0][0], 1, &m_tmp2[0], 1);

        // tmp2 = v/(mu2 + u)
        Vmath::Vdiv(m_nq, &inarray[1][0], 1, &m_tmp2[0], 1, &m_tmp2[0], 1);

        // tmp2 = mu1*v/(mu2 + u)
        Vmath::Smul(m_nq, m_mu1, &m_tmp2[0], 1, &m_tmp2[0], 1);

        // tmp1 = Eps + mu1*v/(mu2+u)
        Vmath::Sadd(m_nq, m_eps, &m_tmp2[0], 1, &m_tmp2[0], 1);

        // tmp1 = (-a-1) + u
//        if (m_spatialParameters->Exists("a"))
//        {
//          Vmath::Vsub(m_nq, &inarray[0][0], 1,
//                          &m_spatialParameters->GetData("a")->GetPhys()[0], 1,
//                          &m_tmp1[0], 1);
//
//          Vmath::Sadd(m_nq, -1.0, &inarray[0][0], 1, &m_tmp1[0], 1);
//        }
//        else
//        {
          Vmath::Sadd(m_nq, (-m_a-1), &inarray[0][0], 1, &m_tmp1[0], 1);
//        }

        // tmp1 = k(u-a-1)
//        if (m_spatialParameters->Exists("k"))
//        {
//          Vmath::Vmul(m_nq, &m_spatialParameters->GetData("k")->GetPhys()[0], 1,
//                          &m_tmp1[0], 1, &m_tmp1[0], 1);
//        }
//        else
//        {
          Vmath::Smul(m_nq, m_k, &m_tmp1[0], 1, &m_tmp1[0], 1);
//        }

        // tmp1 = ku(u-a-1) + v
        Vmath::Vvtvp(m_nq, &inarray[0][0], 1, &m_tmp1[0], 1, &inarray[1][0], 1,
                         &m_tmp1[0], 1);

        // tmp1 = -ku(u-a-1)-v
        Vmath::Neg(m_nq, &m_tmp1[0], 1);

        // outarray = [Eps + mu1*v/(mu2+u)] * [-ku(u-a-1)-v]
        Vmath::Vmul(m_nq, &m_tmp1[0], 1, &m_tmp2[0], 1, &outarray[1][0], 1);
    }

    /**
     *
     */
    void CellModelAlievPanfilov::v_GenerateSummary(SummaryList& s)
    {
        SolverUtils::AddSummaryItem(s, "Cell model","Aliev-Panfilov");
        SolverUtils::AddSummaryItem(s, "k", m_k);
        SolverUtils::AddSummaryItem(s, "a", m_a);
        SolverUtils::AddSummaryItem(s, "eps", m_eps);
        SolverUtils::AddSummaryItem(s, "mu1", m_mu1);
        SolverUtils::AddSummaryItem(s, "mu2", m_mu2);
    }


    /**
     *
     */
    void CellModelAlievPanfilov::v_SetInitialConditions()
    {
        Vmath::Fill(m_nq, 0.0,  m_cellSol[0],  1);
        Vmath::Fill(m_nq, 0.0,  m_cellSol[1],  1);
    }
}
