///////////////////////////////////////////////////////////////////////////////
//
// File CellModel.cpp
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
// Description: Cell model base class.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VmathArray.hpp>

#include <CardiacEPSolver/CellModels/CellModel.h>

namespace Nektar
{
    CellModelFactory& GetCellModelFactory()
    {
        typedef Loki::SingletonHolder<CellModelFactory,
            Loki::CreateUsingNew,
            Loki::NoDestroy > Type;
        return Type::Instance();
    }

    /**
     * @class CellModel
     *
     * The CellModel class and derived classes implement a range of cell model
     * ODE systems. A cell model comprises a system of ion concentration
     * variables and zero or more gating variables. Gating variables are
     * time-integrated using the Rush-Larsen method and for each variable y,
     * the corresponding y_inf and tau_y value is computed by Update(). The tau
     * values are stored in separate storage to inarray/outarray, #m_gates_tau.
     */

    /**
     * Cell model base class constructor.
     */
    CellModel::CellModel(const LibUtilities::SessionReaderSharedPtr& pSession,
                         const MultiRegions::ExpListSharedPtr& pField)
    {
        m_session = pSession;
        m_field = pField;
        m_timestep = pSession->GetParameter("TimeStep");
        m_nvar = 0;
//        m_spatialParameters = MemoryManager<SpatialDomains::SpatialParameters>
//                                          ::AllocateSharedPtr(nq);


        //m_spatialParameters->Read(m_filename);

        //Array<OneD, NekDouble> x(nq), y(nq), z(nq);
        //m_fields[0]->GetCoords(x,y,z);
        //m_spatialParameters->EvaluateParameters(x,y,z);
    }


    /**
     * Initialise the cell model. Allocate workspace and variable storage.
     */
    void CellModel::Initialise()
    {
        ASSERTL1(m_nvar > 0, "Cell model must have at least 1 variable.");

        int nq = m_field->GetNpoints();
        Array<OneD, NekDouble> xc0(nq), xc1(nq), xc2(nq);
        m_field->GetCoords(xc0, xc1, xc2);

        m_cellSol = Array<OneD, Array<OneD, NekDouble> >(m_nvar);
        m_wsp = Array<OneD, Array<OneD, NekDouble> >(m_nvar);
        for (unsigned int i = 0; i < m_nvar; ++i)
        {
//            LibUtilities::EquationSharedPtr f = m_session->GetFunction("InitialConditions", i);
            m_cellSol[i] = Array<OneD, NekDouble>(nq);
//            for (unsigned int j = 0; j < nq; ++j)
//            {
//                m_cellSol[i][j] = f->Evaluate(xc0[j], xc1[j], xc2[j]);
//            }
            m_wsp[i] = Array<OneD, NekDouble>(nq);
        }
        m_gates_tau = Array<OneD, Array<OneD, NekDouble> >(m_gates.size());
        for (unsigned int i = 0; i < m_gates.size(); ++i)
        {
            m_gates_tau[i] = Array<OneD, NekDouble>(nq);
        }

        v_SetInitialConditions();
    }

    /**
     * Integrates the cell model for one PDE time-step. Cell model is
     * sub-stepped.
     *
     * Ion concentrations and membrane potential are integrated using forward
     * Euler, while gating variables are integrated using the Rush-Larsen
     * scheme.
     */
    void CellModel::TimeIntegrate(
            const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                  Array<OneD,       Array<OneD, NekDouble> > &outarray,
            const NekDouble time)
    {
        int nq = m_nq;
        int nvar = inarray.num_elements();
        int substeps = 10;
        NekDouble delta_t = m_timestep/substeps;

        // Copy new transmembrane potential into cell model
        Vmath::Vcopy(nq, inarray[0], 1, m_cellSol[0], 1);

        // Perform substepping
        for (unsigned int i = 0; i < substeps - 1; ++i)
        {
            Update(m_cellSol, m_wsp, time);
            // Voltage
            Vmath::Svtvp(nq, delta_t, m_wsp[0], 1, m_cellSol[0], 1, m_cellSol[0], 1);
            // Ion concentrations
            for (unsigned int j = 0; j < m_concentrations.size(); ++j)
            {
                Vmath::Svtvp(nq, delta_t, m_wsp[m_concentrations[j]], 1, m_cellSol[m_concentrations[j]], 1, m_cellSol[m_concentrations[j]], 1);
            }
            // Gating variables: Rush-Larsen scheme
            for (unsigned int j = 0; j < m_gates.size(); ++j)
            {
                Vmath::Sdiv(nq, -delta_t, m_gates_tau[j], 1, m_gates_tau[j], 1);
                Vmath::Vexp(nq, m_gates_tau[j], 1, m_gates_tau[j], 1);
                Vmath::Vsub(nq, m_cellSol[m_gates[j]], 1, m_wsp[m_gates[j]], 1, m_cellSol[m_gates[j]], 1);
                Vmath::Vvtvp(nq, m_cellSol[m_gates[j]], 1, m_gates_tau[j], 1, m_wsp[m_gates[j]], 1, m_cellSol[m_gates[j]], 1);
            }
        }

        // Perform final cell model step
        Update(m_cellSol, m_wsp, time);

        // Output dV/dt from last step but integrate remaining cell model vars
        Vmath::Vcopy(nq, m_wsp[0], 1, outarray[0], 1);

        // Ion concentrations
        for (unsigned int j = 0; j < m_concentrations.size(); ++j)
        {
            Vmath::Svtvp(nq, delta_t, m_wsp[m_concentrations[j]], 1, m_cellSol[m_concentrations[j]], 1, m_cellSol[m_concentrations[j]], 1);
        }

        // Gating variables: Rush-Larsen scheme
        for (unsigned int j = 0; j < m_gates.size(); ++j)
        {
            Vmath::Sdiv(nq, -delta_t, m_gates_tau[j], 1, m_gates_tau[j], 1);
            Vmath::Vexp(nq, m_gates_tau[j], 1, m_gates_tau[j], 1);
            Vmath::Vsub(nq, m_cellSol[m_gates[j]], 1, m_wsp[m_gates[j]], 1, m_cellSol[m_gates[j]], 1);
            Vmath::Vvtvp(nq, m_cellSol[m_gates[j]], 1, m_gates_tau[j], 1, m_wsp[m_gates[j]], 1, m_cellSol[m_gates[j]], 1);
        }
    }
}
