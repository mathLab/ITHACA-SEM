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

#include <StdRegions/StdNodalTriExp.h>
//#include <LibUtilities/LinearAlgebra/Blas.hpp>

using namespace std;

namespace Nektar
{
    CellModelFactory& GetCellModelFactory()
    {
        static CellModelFactory instance;
        return instance;
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
        m_lastTime = 0.0;
        m_substeps = pSession->GetParameter("Substeps");
        m_nvar = 0;
        m_useNodal = false;

        // Number of points in nodal space is the number of coefficients
        // in modified basis
        std::set<enum LibUtilities::ShapeType> s;
        for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
        {
            s.insert(m_field->GetExp(i)->DetShapeType());
        }

        // Use nodal projection if only triangles
        if (s.size() == 1 && (s.count(LibUtilities::eTriangle) == 1 ||
                              s.count(LibUtilities::eTetrahedron) == 1))
        {
            // This is disabled for now as it causes problems at high order.
            // m_useNodal = true;
        }

        // ---------------------------
        // Move to nodal points
        if (m_useNodal)
        {
            m_nq = pField->GetNcoeffs();
            int order = m_field->GetExp(0)->GetBasis(0)->GetNumModes();

            // Set up a nodal tri
            LibUtilities::BasisKey B0(
                LibUtilities::eModified_A, order,
                LibUtilities::PointsKey(order, LibUtilities::eGaussLobattoLegendre));
            LibUtilities::BasisKey B1(
                LibUtilities::eModified_B, order,
                LibUtilities::PointsKey(order, LibUtilities::eGaussRadauMAlpha1Beta0));
            LibUtilities::BasisKey B2(
                LibUtilities::eModified_C, order,
                LibUtilities::PointsKey(order, LibUtilities::eGaussRadauMAlpha2Beta0));

            m_nodalTri = MemoryManager<StdRegions::StdNodalTriExp>::AllocateSharedPtr(
                        B0, B1, LibUtilities::eNodalTriEvenlySpaced);
            m_nodalTet = MemoryManager<StdRegions::StdNodalTetExp>::AllocateSharedPtr(
                        B0, B1, B2, LibUtilities::eNodalTetEvenlySpaced);
        }
        else
        {
            m_nq = pField->GetTotPoints();
        }
    }


    /**
     * Initialise the cell model. Allocate workspace and variable storage.
     */
    void CellModel::Initialise()
    {
        ASSERTL1(m_nvar > 0, "Cell model must have at least 1 variable.");

        m_cellSol = Array<OneD, Array<OneD, NekDouble> >(m_nvar);
        m_wsp = Array<OneD, Array<OneD, NekDouble> >(m_nvar);
        for (unsigned int i = 0; i < m_nvar; ++i)
        {
            m_cellSol[i] = Array<OneD, NekDouble>(m_nq);
            m_wsp[i] = Array<OneD, NekDouble>(m_nq);
        }
        m_gates_tau = Array<OneD, Array<OneD, NekDouble> >(m_gates.size());
        for (unsigned int i = 0; i < m_gates.size(); ++i)
        {
            m_gates_tau[i] = Array<OneD, NekDouble>(m_nq);
        }

        if (m_session->DefinesFunction("CellModelInitialConditions"))
        {
            LoadCellModel();
        }
        else
        {
            v_SetInitialConditions();
        }
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
        int phys_offset = 0;
        int coef_offset = 0;
        int nvar = inarray.size();
        Array<OneD, NekDouble> tmp;

        // ---------------------------
        // Check nodal temp array set up
        if (m_useNodal)
        {
            if (!m_nodalTmp.size())
            {
                m_nodalTmp = Array<OneD, Array<OneD, NekDouble> >(nvar);
                for (unsigned int k = 0; k < nvar; ++k)
                {
                    m_nodalTmp[k] = Array<OneD, NekDouble>(m_nq);
                }
            }

            // Move to nodal points
            Array<OneD, NekDouble> tmpCoeffs(max(m_nodalTri->GetNcoeffs(), m_nodalTet->GetNcoeffs()));

            for (unsigned int k = 0; k < nvar; ++k)
            {
                for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
                {
                    phys_offset = m_field->GetPhys_Offset(i);
                    coef_offset = m_field->GetCoeff_Offset(i);
                    if (m_field->GetExp(0)->DetShapeType() == LibUtilities::eTriangle)
                    {
                        m_field->GetExp(0)->FwdTrans(inarray[k] + phys_offset, tmpCoeffs);
                        m_nodalTri->ModalToNodal(tmpCoeffs, tmp=m_nodalTmp[k]+coef_offset);
                    }
                    else
                    {
                        m_field->GetExp(0)->FwdTrans(inarray[k] + phys_offset, tmpCoeffs);
                        m_nodalTet->ModalToNodal(tmpCoeffs, tmp=m_nodalTmp[k]+coef_offset);
                    }
                }
            }
            // Copy new transmembrane potential into cell model
            Vmath::Vcopy(m_nq, m_nodalTmp[0], 1, m_cellSol[0], 1);
        }
        else
        {
            // Copy new transmembrane potential into cell model
            Vmath::Vcopy(m_nq, inarray[0], 1, m_cellSol[0], 1);
        }
        // -------------------------

        NekDouble delta_t = (time - m_lastTime)/m_substeps;


        // Perform substepping
        for (unsigned int i = 0; i < m_substeps - 1; ++i)
        {
            Update(m_cellSol, m_wsp, time);
            // Voltage
            Vmath::Svtvp(m_nq, delta_t, m_wsp[0], 1, m_cellSol[0], 1, m_cellSol[0], 1);
            // Ion concentrations
            for (unsigned int j = 0; j < m_concentrations.size(); ++j)
            {
                Vmath::Svtvp(m_nq, delta_t, m_wsp[m_concentrations[j]], 1, m_cellSol[m_concentrations[j]], 1, m_cellSol[m_concentrations[j]], 1);
            }
            // Gating variables: Rush-Larsen scheme
            for (unsigned int j = 0; j < m_gates.size(); ++j)
            {
                Vmath::Sdiv(m_nq, -delta_t, m_gates_tau[j], 1, m_gates_tau[j], 1);
                Vmath::Vexp(m_nq, m_gates_tau[j], 1, m_gates_tau[j], 1);
                Vmath::Vsub(m_nq, m_cellSol[m_gates[j]], 1, m_wsp[m_gates[j]], 1, m_cellSol[m_gates[j]], 1);
                Vmath::Vvtvp(m_nq, m_cellSol[m_gates[j]], 1, m_gates_tau[j], 1, m_wsp[m_gates[j]], 1, m_cellSol[m_gates[j]], 1);
            }
        }

        // Perform final cell model step
        Update(m_cellSol, m_wsp, time);

        // Output dV/dt from last step but integrate remaining cell model vars
        // Transform cell model I_total from nodal to modal space
        if (m_useNodal)
        {
            Array<OneD, NekDouble> tmpCoeffs(max(m_nodalTri->GetNcoeffs(), m_nodalTet->GetNcoeffs()));

            for (unsigned int k = 0; k < nvar; ++k)
            {
                for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
                {
                    int phys_offset = m_field->GetPhys_Offset(i);
                    int coef_offset = m_field->GetCoeff_Offset(i);
                    if (m_field->GetExp(0)->DetShapeType() == LibUtilities::eTriangle)
                    {
                        m_nodalTri->NodalToModal(m_wsp[k]+coef_offset, tmpCoeffs);
                        m_field->GetExp(0)->BwdTrans(tmpCoeffs, tmp=outarray[k] + phys_offset);
                    }
                    else
                    {
                        m_nodalTet->NodalToModal(m_wsp[k]+coef_offset, tmpCoeffs);
                        m_field->GetExp(0)->BwdTrans(tmpCoeffs, tmp=outarray[k] + phys_offset);
                    }
                }
            }
        }
        else
        {
            Vmath::Vcopy(m_nq, m_wsp[0], 1, outarray[0], 1);
        }

        // Ion concentrations
        for (unsigned int j = 0; j < m_concentrations.size(); ++j)
        {
            Vmath::Svtvp(m_nq, delta_t, m_wsp[m_concentrations[j]], 1, m_cellSol[m_concentrations[j]], 1, m_cellSol[m_concentrations[j]], 1);
        }

        // Gating variables: Rush-Larsen scheme
        for (unsigned int j = 0; j < m_gates.size(); ++j)
        {
            Vmath::Sdiv(m_nq, -delta_t, m_gates_tau[j], 1, m_gates_tau[j], 1);
            Vmath::Vexp(m_nq, m_gates_tau[j], 1, m_gates_tau[j], 1);
            Vmath::Vsub(m_nq, m_cellSol[m_gates[j]], 1, m_wsp[m_gates[j]], 1, m_cellSol[m_gates[j]], 1);
            Vmath::Vvtvp(m_nq, m_cellSol[m_gates[j]], 1, m_gates_tau[j], 1, m_wsp[m_gates[j]], 1, m_cellSol[m_gates[j]], 1);
        }

        m_lastTime = time;
    }

    Array<OneD, NekDouble> CellModel::GetCellSolutionCoeffs(unsigned int idx)
    {
        ASSERTL0(idx < m_nvar, "Index out of range for cell model.");

        Array<OneD, NekDouble> outarray(m_field->GetNcoeffs());
        Array<OneD, NekDouble> tmp;

        if (m_useNodal)
        {
            for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
            {
                int coef_offset = m_field->GetCoeff_Offset(i);
                if (m_field->GetExp(0)->DetShapeType() == LibUtilities::eTriangle)
                {
                    m_nodalTri->NodalToModal(m_cellSol[idx]+coef_offset, tmp=outarray+coef_offset);
                }
                else
                {
                    m_nodalTet->NodalToModal(m_cellSol[idx]+coef_offset, tmp=outarray+coef_offset);
                }
            }
        }
        else
        {
            m_field->FwdTrans_IterPerExp(m_cellSol[idx], outarray);
        }

        return outarray;
    }

    Array<OneD, NekDouble> CellModel::GetCellSolution(unsigned int idx)
    {
        return m_cellSol[idx];
    }

    void CellModel::LoadCellModel()
    {
        const bool root = (m_session->GetComm()->GetRank() == 0);
        const std::string fncName = "CellModelInitialConditions";
        const int nvar = m_cellSol[0].size();
        std::string varName;
        Array<OneD, NekDouble> coeffs(m_field->GetNcoeffs());
        Array<OneD, NekDouble> tmp;
        int j = 0;

        SpatialDomains::MeshGraphSharedPtr vGraph = m_field->GetGraph();

        if (root)
        {
            cout << "Cell model initial conditions: " << endl;
        }

        // First determine all the files we need to load
        std::set<std::string> filelist;
        for (j = 1; j < nvar; ++j)
        {
            // Get the name of the jth variable
            varName = GetCellVarName(j);

            if (m_session->GetFunctionType(fncName, varName) ==
                    LibUtilities::eFunctionTypeFile)
            {
                filelist.insert(m_session->GetFunctionFilename(fncName,
                                                               varName));
            }
        }

        // Read files
        typedef std::vector<LibUtilities::FieldDefinitionsSharedPtr> FDef;
        typedef std::vector<std::vector<NekDouble> > FData;
        std::map<std::string, FDef>  FieldDef;
        std::map<std::string, FData> FieldData;
        LibUtilities::FieldMetaDataMap fieldMetaDataMap;

        for (auto &setIt : filelist)
        {
            if (root)
            {
                cout << "  - Reading file: " << setIt << endl;
            }
            FieldDef[setIt] = FDef(0);
            FieldData[setIt] = FData(0);
            LibUtilities::FieldIOSharedPtr fld =
                LibUtilities::FieldIO::CreateForFile(m_session, setIt);
            fld->Import(setIt, FieldDef[setIt], FieldData[setIt],
                        fieldMetaDataMap);
        }

        // Get time of checkpoint from file if available
        auto iter = fieldMetaDataMap.find("Time");
        if(iter != fieldMetaDataMap.end())
        {
             m_lastTime = boost::lexical_cast<NekDouble>(iter->second);
        }

        // Load each cell model variable
        // j=0 and j=1 are for transmembrane or intra/extra-cellular volt.
        Vmath::Zero(m_nq, m_cellSol[0], 1);
        for(j = 1; j < m_cellSol.size(); ++j)
        {
            // Get the name of the jth variable
            varName = GetCellVarName(j);

            // Check if this variable is defined in a file or analytically
            if (m_session->GetFunctionType(fncName, varName) ==
                    LibUtilities::eFunctionTypeFile)
            {
                const std::string file =
                        m_session->GetFunctionFilename(fncName, varName);

                if (root)
                {
                    cout << "  - Field " << varName << ": from file "
                         << file << endl;
                }

                // Extract the data into the modal coefficients
                for(int i = 0; i < FieldDef[file].size(); ++i)
                {
                    m_field->ExtractDataToCoeffs(FieldDef[file][i],
                                                 FieldData[file][i],
                                                 varName,
                                                 coeffs);
                }

                // If using nodal cell model then we do a modal->nodal transform
                // otherwise we do a backward transform onto physical points.
                if (m_useNodal)
                {
                    for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
                    {
                        int coef_offset = m_field->GetCoeff_Offset(i);
                        if (m_field->GetExp(0)->DetShapeType() ==
                                LibUtilities::eTriangle)
                        {
                            m_nodalTri->ModalToNodal(coeffs+coef_offset,
                                                tmp=m_cellSol[j]+coef_offset);
                        }
                        else
                        {
                            m_nodalTet->ModalToNodal(coeffs+coef_offset,
                                                tmp=m_cellSol[j]+coef_offset);
                        }
                    }
                }
                else
                {
                    m_field->BwdTrans(coeffs, m_cellSol[j]);
                }
            }
            else if (m_session->GetFunctionType(fncName, varName) ==
                    LibUtilities::eFunctionTypeExpression)
            {
                LibUtilities::EquationSharedPtr equ =
                        m_session->GetFunction(fncName, varName);

                if (root)
                {
                    cout << "  - Field " << varName << ": "
                         << equ->GetExpression() << endl;
                }

                const unsigned int nphys = m_field->GetNpoints();
                Array<OneD, NekDouble> x0(nphys);
                Array<OneD, NekDouble> x1(nphys);
                Array<OneD, NekDouble> x2(nphys);
                m_field->GetCoords(x0,x1,x2);

                if (m_useNodal)
                {
                    Array<OneD, NekDouble> phys(nphys);
                    Array<OneD, NekDouble> tmpCoeffs(max(m_nodalTri->GetNcoeffs(), m_nodalTet->GetNcoeffs()));

                    equ->Evaluate(x0, x1, x2, phys);
                    for (unsigned int i = 0; i < m_field->GetNumElmts(); ++i)
                    {
                        int phys_offset = m_field->GetPhys_Offset(i);
                        int coef_offset = m_field->GetCoeff_Offset(i);
                        if (m_field->GetExp(0)->DetShapeType() ==
                                LibUtilities::eTriangle)
                        {
                            m_field->GetExp(0)->FwdTrans(
                                            phys + phys_offset, tmpCoeffs);
                            m_nodalTri->ModalToNodal(
                                            tmpCoeffs,
                                            tmp = m_cellSol[j] + coef_offset);
                        }
                        else
                        {
                            m_field->GetExp(0)->FwdTrans(
                                            phys + phys_offset, tmpCoeffs);
                            m_nodalTet->ModalToNodal(
                                            tmpCoeffs,
                                            tmp = m_cellSol[j] + coef_offset);
                        }
                    }
                }
                else
                {
                    equ->Evaluate(x0, x1, x2, m_cellSol[j]);
                }
            }
        }
    }
}
