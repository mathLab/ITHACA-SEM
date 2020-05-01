///////////////////////////////////////////////////////////////////////////////
//
// File IterativeElasticSystem.cpp
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
// Description: IterativeElasticSystem solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/format.hpp>

#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTriExp.h>
#include <LocalRegions/MatrixKey.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>
#include <GlobalMapping/Deform.h>

#include <LinearElasticSolver/EquationSystems/IterativeElasticSystem.h>

using namespace std;

namespace Nektar
{

string IterativeElasticSystem::className = GetEquationSystemFactory().
    RegisterCreatorFunction("IterativeElasticSystem",
                            IterativeElasticSystem::create);

IterativeElasticSystem::IterativeElasticSystem(
    const LibUtilities::SessionReaderSharedPtr& pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : LinearElasticSystem(pSession, pGraph)
{
}

void IterativeElasticSystem::v_InitObject()
{
    LinearElasticSystem::v_InitObject();

    const int nVel = m_fields[0]->GetCoordim(0);

    // Read in number of steps to take.
    m_session->LoadParameter("NumSteps", m_numSteps, 0);
    ASSERTL0(m_numSteps > 0, "You must specify at least one step");

    // Read in whether to repeatedly apply boundary conditions (for e.g.
    // rotation purposes).
    string bcType;
    m_session->LoadSolverInfo("BCType", bcType, "Normal");
    m_repeatBCs = bcType != "Normal";

    if (!m_repeatBCs)
    {
        // Loop over BCs, identify which ones we need to deform.
        const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
            &bndCond = m_fields[0]->GetBndConditions();

        for (int i = 0; i < bndCond.size(); ++i)
        {
            if (boost::iequals(bndCond[i]->GetUserDefined(), "Wall"))
            {
                m_toDeform.push_back(i);
            }
        }

        int numDeform = m_toDeform.size();
        m_comm->AllReduce(numDeform, LibUtilities::ReduceMax);
        ASSERTL0(numDeform > 0, "You must specify at least one WALL tag on"
                 "a boundary condition");

        m_savedBCs  = Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(
            m_toDeform.size());

        for (int i = 0; i < m_toDeform.size(); ++i)
        {
            m_savedBCs[i] = Array<OneD, Array<OneD, NekDouble> >(nVel);
            for (int j = 0; j < nVel; ++j)
            {
                const int id = m_toDeform[i];
                MultiRegions::ExpListSharedPtr bndCondExp =
                    m_fields[j]->GetBndCondExpansions()[id];
                int nCoeffs = bndCondExp->GetNcoeffs();

                m_savedBCs[i][j] = Array<OneD, NekDouble>(nCoeffs);
                Vmath::Smul(nCoeffs, 1.0/m_numSteps,
                            bndCondExp->GetCoeffs(),    1,
                            bndCondExp->UpdateCoeffs(), 1);
                Vmath::Vcopy(nCoeffs, bndCondExp->GetCoeffs(), 1,
                             m_savedBCs[i][j], 1);
            }
        }
    }
}

void IterativeElasticSystem::v_GenerateSummary(SolverUtils::SummaryList& s)
{
    LinearElasticSystem::v_GenerateSummary(s);
}

void IterativeElasticSystem::v_DoSolve()
{
    int i, j, k;

    // Write initial geometry for consistency/script purposes
    WriteGeometry(0);

    // Now loop over desired number of steps
    for (i = 1; i <= m_numSteps; ++i)
    {
        int invalidElmtId = -1;

        // Perform solve for this iteration and update geometry accordingly.
        LinearElasticSystem::v_DoSolve();
        GlobalMapping::UpdateGeometry(m_graph, m_fields);
        WriteGeometry(i);

        // Check for invalid elements.
        for (j = 0; j < m_fields[0]->GetExpSize(); ++j)
        {
            SpatialDomains::GeomFactorsSharedPtr geomFac =
                m_fields[0]->GetExp(j)->GetGeom()->GetGeomFactors();

            if (!geomFac->IsValid())
            {
                invalidElmtId =
                    m_fields[0]->GetExp(j)->GetGeom()->GetGlobalID();
                break;
            }
        }

        m_session->GetComm()->AllReduce(invalidElmtId, LibUtilities::ReduceMax);

        // If we found an invalid element, exit loop without writing output.
        if (invalidElmtId >= 0)
        {
            if (m_session->GetComm()->GetRank() == 0)
            {
                cout << "- Detected negative Jacobian in element "
                     << invalidElmtId << "; terminating at"
                    " step: "<<i<< endl;
            }

            break;
        }

        if (m_session->GetComm()->GetRank() == 0)
        {
            cout << "Step: " << i << endl;
        }

        // Update boundary conditions
        if (m_repeatBCs)
        {
            for (j = 0; j < m_fields.size(); ++j)
            {
                string varName = m_session->GetVariable(j);
                m_fields[j]->EvaluateBoundaryConditions(m_time, varName);
            }
        }
        else
        {
            for (j = 0; j < m_fields.size(); ++j)
            {
                const Array<OneD, const MultiRegions::ExpListSharedPtr>
                    &bndCondExp = m_fields[j]->GetBndCondExpansions();

                for (k = 0; k < m_toDeform.size(); ++k)
                {
                    const int id = m_toDeform[k];
                    const int nCoeffs = bndCondExp[id]->GetNcoeffs();
                    Vmath::Vcopy(nCoeffs,
                                 m_savedBCs[k][j],               1,
                                 bndCondExp[id]->UpdateCoeffs(), 1);
                }
            }
        }
    }
}

/**
 * @brief Write out a file in serial or directory in parallel containing new
 * mesh geometry.
 */
void IterativeElasticSystem::WriteGeometry(const int i)
{
    fs::path filename;
    stringstream s;
    s << m_session->GetSessionName() << "-" << i;

    if (m_session->GetComm()->GetSize() > 1)
    {
        s << "_xml";

        string ss = s.str();
        if(!fs::is_directory(ss))
        {
            fs::create_directory(ss);
        }

        boost::format pad("P%1$07d.xml");
        pad % m_session->GetComm()->GetRank();
        filename = fs::path(ss) / fs::path(pad.str());
    }
    else
    {
        s << ".xml";
        filename = fs::path(s.str());
    }

    string fname = LibUtilities::PortablePath(filename);
    m_fields[0]->GetGraph()->WriteGeometry(fname);
}

}
