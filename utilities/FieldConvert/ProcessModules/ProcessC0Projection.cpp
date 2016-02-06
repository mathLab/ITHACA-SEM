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
//  License for the specific language governing rights and limitations under
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

#include <string>
#include <iostream>
using namespace std;

#include "ProcessC0Projection.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessC0Projection::className =
GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "C0Projection"),
    ProcessC0Projection::create, "Computes C0 projection.");

ProcessC0Projection::ProcessC0Projection(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["fields"] = ConfigOption(false,"All","Start field to project");
}

ProcessC0Projection::~ProcessC0Projection()
{
}

void ProcessC0Projection::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessC0Projection: Projecting field into C0 space..."
             << endl;
    }

    // ensure not using diagonal preconditioner since tends not to converge fo
    // mass matrix
    if(m_f->m_graph->GetMeshDimension() == 3)
    {
        if(boost::iequals(m_f->m_session->GetSolverInfo("GLOBALSYSSOLN"),
                          "IterativeStaticCond"))
        {
            if(boost::iequals(m_f->m_session->GetSolverInfo("PRECONDITIONER"),
                              "Diagonal"))
            {
                m_f->m_session->SetSolverInfo("PRECONDITIONER","LowEnergyBlock");
            }
            if(boost::iequals(m_f->m_session->GetSolverInfo("PRECONDITIONER"),
                              "FullLinearSpaceWithDiagonal"))
            {
                m_f->m_session->SetSolverInfo("PRECONDITIONER","FullLinearSpaceWithLowEnergyBlock");
            }
            
            if(m_f->m_verbose)
            {
                cout << "Resetting diagonal precondition to low energy block " << endl;;
            }
        }
    }
    

    // generate an C0 expansion field with no boundary conditions.
    bool savedef = m_f->m_declareExpansionAsContField;
    m_f->m_declareExpansionAsContField = true;
    MultiRegions::ExpListSharedPtr C0ProjectExp = 
        m_f->AppendExpList(m_f->m_fielddef[0]->m_numHomogeneousDir,
                           "DefaultVar",true);
    m_f->m_declareExpansionAsContField = savedef;

    int nfields = m_f->m_exp.size();

    string fields = m_config["fields"].as<string>();
    vector<unsigned int> processFields;

    if(fields.compare("All") == 0)
    {
        for(int i = 0; i < nfields; ++i)
        {
            processFields.push_back(i);
        }
    }
    else
    {
        ASSERTL0(ParseUtils::GenerateOrderedVector(fields.c_str(),
                                                   processFields),
                 "Failed to interpret field string in C0Projection");
    }

    for (int i = 0; i < processFields.size(); ++i)
    {
        ASSERTL0(processFields[i] < nfields,
                 "Attempt to process field that is larger than then number of "
                 "fields available");

        if (m_f->m_verbose)
        {
            cout << "\t Processing field: " << processFields[i] << endl;
        }
        C0ProjectExp->FwdTrans(m_f->m_exp[processFields[i]]->GetPhys(),
                                 m_f->m_exp[processFields[i]]->UpdateCoeffs());
        C0ProjectExp->BwdTrans(m_f->m_exp[processFields[i]]->GetCoeffs(),
                                 m_f->m_exp[processFields[i]]->UpdatePhys());
    }

    // reset FieldDef in case of serial input and parallel output
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = m_f->m_exp[0]->GetFieldDefinitions();
    // reset up FieldData with new values before projecting
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for(int i = 0; i < nfields; ++i)
    {
        for (int j = 0; j < FieldDef.size(); ++j)
        {
            FieldDef[j]->m_fields.push_back(m_f->m_fielddef[0]->m_fields[i]);
            m_f->m_exp[i]->AppendFieldData(FieldDef[j], FieldData[j]);
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}

}
}
