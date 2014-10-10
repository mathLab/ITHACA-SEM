////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessInterpPointDataToFld.cpp
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
//  Description: Interpolate, using finite different approximation,
//  given data to a fld file
//
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <iostream>
using namespace std;

#include "ProcessInterpPointDataToFld.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessInterpPointDataToFld::className =
    GetModuleFactory().RegisterCreatorFunction(
       ModuleKey(eProcessModule, "interppointdatatofld"),
       ProcessInterpPointDataToFld::create,
       "Interpolates given discrete data using a finite difference "
       "approximation to a fld file given a xml file");


ProcessInterpPointDataToFld::ProcessInterpPointDataToFld(FieldSharedPtr f)
        : ProcessModule(f)
{

    m_config["interpcoord"] = ConfigOption(false, "0",
                                    "coordinate id ot use for interpolation");

}

ProcessInterpPointDataToFld::~ProcessInterpPointDataToFld()
{
}

void ProcessInterpPointDataToFld::Process(po::variables_map &vm)
{
    int i,j;
    if(m_f->m_verbose)
    {
        cout << "Processing point data interpolation (linear)" << endl;
    }

    // Check for command line point specification if no .pts file specified
    ASSERTL0(m_f->m_fieldPts != NullFieldPts,
             "No input points found");

    ASSERTL0(m_f->m_fieldPts->m_nFields > 0,
             "No field values provided in input");

    // assume one field is already defined from input file.
    m_f->m_exp.resize(m_f->m_fieldPts->m_nFields+1);
    for(i = 1; i < m_f->m_fieldPts->m_nFields; ++i)
    {
        m_f->m_exp[i] = m_f->AppendExpList(0);
    }

    if(m_f->m_session->GetComm()->GetRank() == 0)
    {
        cout << "Interpolating [" << flush;
    }

    int totpoints = m_f->m_exp[0]->GetTotPoints();
    Array<OneD,NekDouble> coords[3];

    coords[0] = Array<OneD,NekDouble>(totpoints);
    coords[1] = Array<OneD,NekDouble>(totpoints);
    coords[2] = Array<OneD,NekDouble>(totpoints);

    m_f->m_exp[0]->GetCoords(coords[0],coords[1],coords[2]);

    int coord_id = m_config["interpcoord"].as<int>();

    // interpolate points and transform
    Array<OneD, NekDouble> intfields(m_f->m_fieldPts->m_nFields);
    for(i = 0; i < totpoints; ++i)
    {
        m_f->m_fieldPts->Interp1DPts(coords[coord_id][i],intfields);
        for(j = 0; j < m_f->m_fieldPts->m_nFields; ++j)
        {
            m_f->m_exp[j]->SetPhys(i,intfields[j]);
        }
    }

    if(m_f->m_session->GetComm()->GetRank() == 0)
    {
        cout << "]" << endl;
    }

    // forward transform fields
    for(i = 0; i < m_f->m_fieldPts->m_nFields; ++i)
    {
        m_f->m_exp[i]->FwdTrans_IterPerExp(m_f->m_exp[i]->GetPhys(),
                                           m_f->m_exp[i]->UpdateCoeffs());
    }


    // set up output fld file.
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for (j = 0; j < m_f->m_fieldPts->m_nFields; ++j)
    {
        for (i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back(m_f->m_fieldPts->m_fields[j]);

            m_f->m_exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;

}

}
}


