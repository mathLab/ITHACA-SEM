////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessC0Projection.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Computes C0 projection.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessC0Projection.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessC0Projection::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "C0Projection"),
        ProcessC0Projection::create,
        "Computes C0 projection.");

ProcessC0Projection::ProcessC0Projection(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["fields"] = ConfigOption(false, "All", "Start field to project");
    m_config["localtoglobalmap"] = ConfigOption(
        true, "0", "Just perform a local to global mapping and back");
    m_config["usexmlbcs"] = ConfigOption(
        true, "0", "Use boundary conditions given in xml file. Requires all "
                   "projected fields to be defined in xml file");
    m_config["helmsmoothing"] = ConfigOption(
        false, "Not Set", "Use a Helmholtz smoother to remove high frequency "
                          "components above specified L");

    f->m_declareExpansionAsContField = true;
}

ProcessC0Projection::~ProcessC0Projection()
{
}

void ProcessC0Projection::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    // ensure not using diagonal preconditioner since tends not to converge fo
    // mass matrix
    if (m_f->m_graph->GetMeshDimension() == 3)
    {
        if (boost::iequals(m_f->m_session->GetSolverInfo("GLOBALSYSSOLN"),
                           "IterativeStaticCond"))
        {
            if (boost::iequals(m_f->m_session->GetSolverInfo("PRECONDITIONER"),
                               "Diagonal"))
            {
                m_f->m_session->SetSolverInfo("PRECONDITIONER",
                                              "LowEnergyBlock");
            }
            if (boost::iequals(m_f->m_session->GetSolverInfo("PRECONDITIONER"),
                               "FullLinearSpaceWithDiagonal"))
            {
                m_f->m_session->SetSolverInfo(
                    "PRECONDITIONER", "FullLinearSpaceWithLowEnergyBlock");
            }

            if (m_f->m_verbose)
            {
                if (m_f->m_comm->GetRank() == 0)
                {
                    cout << "Resetting diagonal precondition to low energy "
                            "block "
                         << endl;
                }
            }
        }
    }
    bool JustPerformLocToGloMap = m_config["localtoglobalmap"].as<bool>();
    bool HelmSmoother =
        (boost::iequals(m_config["helmsmoothing"].as<string>(), "Not Set"))
            ? false
            : true;
    int nfields = m_f->m_exp.size();
    Array<OneD, MultiRegions::ExpListSharedPtr> C0ProjectExp(nfields);
    if (m_config["usexmlbcs"].as<bool>())
    {
        for (int i = 0; i < nfields; ++i)
        {
            C0ProjectExp[i] = m_f->m_exp[i];
        }
    }
    else
    {
        // generate a C0 expansion field with no boundary conditions.
        bool savedef                       = m_f->m_declareExpansionAsContField;
        bool savedef2                      = m_f->m_requireBoundaryExpansion;
        m_f->m_declareExpansionAsContField = true;
        m_f->m_requireBoundaryExpansion    = false;
        C0ProjectExp[0]                    = m_f->AppendExpList(
            m_f->m_numHomogeneousDir, "DefaultVar", true);
        m_f->m_declareExpansionAsContField = savedef;
        m_f->m_requireBoundaryExpansion    = savedef2;
        for (int i = 1; i < nfields; ++i)
        {
            C0ProjectExp[i] = C0ProjectExp[0];
        }
    }

    string fields = m_config["fields"].as<string>();
    vector<unsigned int> processFields;

    if (fields.compare("All") == 0)
    {
        for (int i = 0; i < nfields; ++i)
        {
            processFields.push_back(i);
        }
    }
    else
    {
        ASSERTL0(
            ParseUtils::GenerateVector(fields, processFields),
            "Failed to interpret field string in C0Projection");
    }

    for (int i = 0; i < processFields.size(); ++i)
    {
        ASSERTL0(processFields[i] < nfields,
                 "Attempt to process field that is larger than then number of "
                 "fields available");

        if (m_f->m_verbose)
        {
            if (m_f->m_comm->GetRank() == 0)
            {
                cout << "\t Processing field: " << processFields[i] << endl;
            }
        }

        if (JustPerformLocToGloMap)
        {
            int ncoeffs = m_f->m_exp[0]->GetNcoeffs();
            Vmath::Vcopy(ncoeffs, m_f->m_exp[processFields[i]]->GetCoeffs(), 1,
                         C0ProjectExp[processFields[i]]->UpdateCoeffs(), 1);
            C0ProjectExp[processFields[i]]->LocalToGlobal();
            C0ProjectExp[processFields[i]]->GlobalToLocal();
            Vmath::Vcopy(ncoeffs, C0ProjectExp[processFields[i]]->GetCoeffs(),
                         1, m_f->m_exp[processFields[i]]->UpdateCoeffs(), 1);
        }
        else
        {
            if (HelmSmoother)
            {
                int dim          = m_f->m_graph->GetSpaceDimension();
                int npoints      = m_f->m_exp[0]->GetNpoints();
                NekDouble lambda = m_config["helmsmoothing"].as<NekDouble>();
                lambda           = 2 * M_PI / lambda;
                lambda           = lambda * lambda;

                if (m_f->m_verbose)
                {
                    cout << "Setting up Helmholtz smoother with lambda = "
                         << lambda << endl;
                }

                StdRegions::ConstFactorMap factors;
                Array<OneD, NekDouble> forcing(npoints);
                factors[StdRegions::eFactorLambda] = -lambda;

                Array<OneD, Array<OneD, NekDouble> > Velocity(dim);
                for (int j = 0; j < dim; ++j)
                {
                    Velocity[j] = Array<OneD, NekDouble>(npoints, 0.0);
                }

                Vmath::Smul(npoints, -lambda,
                            m_f->m_exp[processFields[i]]->GetPhys(), 1, forcing,
                            1);

                // Note we are using the
                // LinearAdvectionDiffusionReaction solver here
                // instead of HelmSolve since lambda is negative and
                // so matrices are not positive definite. Ideally
                // should allow for negative lambda coefficient in
                // HelmSolve
                C0ProjectExp[processFields[i]]
                    ->LinearAdvectionDiffusionReactionSolve(
                        Velocity, forcing,
                        m_f->m_exp[processFields[i]]->UpdateCoeffs(), -lambda);
            }
            else
            {
                C0ProjectExp[processFields[i]]->FwdTrans(
                    m_f->m_exp[processFields[i]]->GetPhys(),
                    m_f->m_exp[processFields[i]]->UpdateCoeffs());
            }
        }
        C0ProjectExp[processFields[i]]->BwdTrans(
                    m_f->m_exp[processFields[i]]->GetCoeffs(),
                    m_f->m_exp[processFields[i]]->UpdatePhys());
    }

}
}
}
