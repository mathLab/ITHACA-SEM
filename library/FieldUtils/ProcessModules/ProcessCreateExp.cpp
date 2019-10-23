////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessCreateExp.cpp
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
//  Description: Dummy module to create m_exp.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include "ProcessCreateExp.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/Timer.h>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessCreateExp::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "createExp"),
        ProcessCreateExp::create,
        "dummy module used to create m_exp.");

ProcessCreateExp::ProcessCreateExp(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessCreateExp::~ProcessCreateExp()
{
}

void ProcessCreateExp::Process(po::variables_map &vm)
{
    if(m_f->m_graph)
    {
        int i, j;
        LibUtilities::Timer timerpart;
        if (m_f->m_verbose)
        {
            if (m_f->m_comm->TreatAsRankZero())
            {
                timerpart.Start();
            }
        }

        // check to see if fld file defined so can use in
        // expansion defintion if required
        bool fldfilegiven = (m_f->m_fielddef.size() != 0);
        bool expFromFld   = fldfilegiven  && !vm.count("useSessionExpansion");

        // load fielddef header if fld file is defined. This gives
        // precedence to Homogeneous definition in fld file
        m_f->m_numHomogeneousDir = 0;
        if (expFromFld)
        {
            m_f->m_numHomogeneousDir = m_f->m_fielddef[0]->m_numHomogeneousDir;

            // Set up Expansion information to use mode order from field
            m_f->m_graph->SetExpansions(m_f->m_fielddef);
        }
        else
        {
            if (m_f->m_session->DefinesSolverInfo("HOMOGENEOUS"))
            {
                std::string HomoStr = 
                        m_f->m_session->GetSolverInfo("HOMOGENEOUS");

                if ((HomoStr == "HOMOGENEOUS1D") ||
                    (HomoStr == "Homogeneous1D") ||
                    (HomoStr == "1D") || (HomoStr == "Homo1D"))
                {
                    m_f->m_numHomogeneousDir = 1;
                }
                if ((HomoStr == "HOMOGENEOUS2D") ||
                    (HomoStr == "Homogeneous2D") ||
                    (HomoStr == "2D") || (HomoStr == "Homo2D"))
                {
                    m_f->m_numHomogeneousDir = 2;
                }
            }
        }

        m_f->m_exp.resize(1);

        // Check  if there are any elements to process
        vector<int> IDs;
        auto domain = m_f->m_graph->GetDomain();
        for(int d = 0; d < domain.size(); ++d)
        {
            for (auto &compIter : domain[d])
            {
                for (auto &x : compIter.second->m_geomVec)
                {
                    IDs.push_back(x->GetGlobalID());
                }
            }
        }

        // if Range has been specified it is possible to have a
        // partition which is empty so check this and return with empty
        // expansion if no elements present.
        if (!IDs.size())
        {
            m_f->m_exp[0] = MemoryManager<MultiRegions::ExpList>::
                            AllocateSharedPtr();
            return;
        }

        // Adjust number of quadrature points
        if (vm.count("output-points"))
        {
            int nPointsNew = vm["output-points"].as<int>();
            m_f->m_graph->SetExpansionsToPointOrder(nPointsNew);
        }

        if (m_f->m_verbose)
        {
            if (m_f->m_comm->TreatAsRankZero())
            {
                timerpart.Stop();
                NekDouble cpuTime = timerpart.TimePerTest(1);

                stringstream ss;
                ss << cpuTime << "s";
                cout << "\t ProcessCreateExp setexpansion CPU Time: "
                     << setw(8) << left
                     << ss.str() << endl;
                timerpart.Start();
            }
        }

        // Override number of planes with value from cmd line
        if (m_f->m_numHomogeneousDir == 1 && vm.count("output-points-hom-z"))
        {
            int expdim = m_f->m_graph->GetMeshDimension();
            m_f->m_fielddef[0]->m_numModes[expdim] =
                vm["output-points-hom-z"].as<int>();
        }

        m_f->m_exp[0] = m_f->SetUpFirstExpList(m_f->m_numHomogeneousDir,
                                                expFromFld);

        if (m_f->m_verbose)
        {
            if (m_f->m_comm->TreatAsRankZero())
            {
                timerpart.Stop();
                NekDouble cpuTime = timerpart.TimePerTest(1);

                stringstream ss1;

                ss1 << cpuTime << "s";
                cout << "\t ProcessCreateExp set first exp CPU Time: "
                     << setw(8)   << left
                     << ss1.str() << endl;
            }
        }

        if (fldfilegiven)
        {
            int nfields, nstrips;

            m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

            vector<string> vars = m_f->m_session->GetVariables();
            if (vm.count("useSessionVariables"))
            {
                m_f->m_variables = vars;
            }
            nfields = m_f->m_variables.size();

            m_f->m_exp.resize(nfields * nstrips);

            // declare other fields;
            for (int s = 0; s < nstrips; ++s) // homogeneous strip varient
            {
                for (i = 0; i < nfields; ++i)
                {
                    if (i < vars.size())
                    {
                        // check to see if field already defined
                        if (!m_f->m_exp[s * nfields + i])
                        {
                            m_f->m_exp[s * nfields + i] = m_f->AppendExpList(
                              m_f->m_numHomogeneousDir, vars[i]);
                        }
                    }
                    else
                    {
                        if (vars.size())
                        {
                            m_f->m_exp[s * nfields + i] = m_f->AppendExpList(
                              m_f->m_numHomogeneousDir, vars[0]);
                        }
                        else
                        {
                            m_f->m_exp[s * nfields + i] = m_f->AppendExpList(
                                m_f->m_numHomogeneousDir);
                        }
                    }
                }
            }

            // Extract data to coeffs and bwd transform
            for (int s = 0; s < nstrips; ++s) // homogeneous strip varient
            {
                for (j = 0; j < nfields; ++j)
                {
                    for (i = 0; i < m_f->m_data.size() / nstrips; ++i)
                    {
                        int n = i * nstrips + s;
                        // In case of multiple flds, we might not have a
                        //   variable in this m_data[n] -> skip in this case
                        auto it = find (m_f->m_fielddef[n]->m_fields.begin(),
                                        m_f->m_fielddef[n]->m_fields.end(),
                                        m_f->m_variables[j]);
                        if(it !=m_f->m_fielddef[n]->m_fields.end())
                        {
                            m_f->m_exp[s * nfields + j]->ExtractDataToCoeffs(
                                m_f->m_fielddef[n],
                                m_f->m_data[n],
                                m_f->m_variables[j],
                                m_f->m_exp[s * nfields + j]->UpdateCoeffs());
                        }
                    }
                    m_f->m_exp[s * nfields + j]->BwdTrans(
                        m_f->m_exp[s * nfields + j]->GetCoeffs(),
                        m_f->m_exp[s * nfields + j]->UpdatePhys());
                }
            }
            // Clear fielddef and data
            //    (they should not be used after running this module)
            m_f->m_fielddef = vector<LibUtilities::FieldDefinitionsSharedPtr>();
            m_f->m_data     = vector<std::vector<NekDouble> >();
        }
    }

}
}
}
