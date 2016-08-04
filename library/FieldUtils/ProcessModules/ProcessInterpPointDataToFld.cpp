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
#include <iostream>
#include <string>

using namespace std;

#include "ProcessInterpPointDataToFld.h"

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
namespace Nektar
{
namespace FieldUtils
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

    m_config["interpcoord"] =
        ConfigOption(false, "-1", "coordinate id to use for interpolation");
}

ProcessInterpPointDataToFld::~ProcessInterpPointDataToFld()
{
}

void ProcessInterpPointDataToFld::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            cout
                << "ProcessInterpPointDataToFld: interpolating data to field..."
                << endl;
        }
    }

    int i, j;

    // Check for command line point specification if no .pts file specified
    ASSERTL0(m_f->m_fieldPts != LibUtilities::NullPtsField,
             "No input points found");

    int nFields = m_f->m_fieldPts->GetNFields();
    ASSERTL0(nFields > 0, "No field values provided in input");

    // assume one field is already defined from input file.
    m_f->m_exp.resize(nFields + 1);
    for (i = 1; i < nFields; ++i)
    {
        m_f->m_exp[i] = m_f->AppendExpList(0);
    }

    int totpoints = m_f->m_exp[0]->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble> > intFields(3 + nFields);
    for (int i = 0; i < 3 + nFields; ++i)
    {
        intFields[i] = Array<OneD, NekDouble>(totpoints);
    }
    m_f->m_exp[0]->GetCoords(intFields[0], intFields[1], intFields[2]);
    LibUtilities::PtsFieldSharedPtr outPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(3, intFields);

    int coord_id = m_config["interpcoord"].as<int>();
    ASSERTL0(coord_id <= m_f->m_fieldPts->GetDim() - 1,
             "interpcoord is bigger than the Pts files dimension");

    Interpolator interp(eNoMethod, coord_id);

    if (m_f->m_comm->GetRank() == 0)
    {
        interp.SetProgressCallback(
            &ProcessInterpPointDataToFld::PrintProgressbar, this);
    }
    interp.Interpolate(m_f->m_fieldPts, outPts);
    if (m_f->m_comm->GetRank() == 0)
    {
        cout << endl;
    }

    for (i = 0; i < totpoints; ++i)
    {
        for (j = 0; j < nFields; ++j)
        {
            m_f->m_exp[j]->SetPhys(i, outPts->GetPointVal(j, i));
        }
    }

    // forward transform fields
    for (i = 0; i < nFields; ++i)
    {
        m_f->m_exp[i]->FwdTrans_IterPerExp(m_f->m_exp[i]->GetPhys(),
                                           m_f->m_exp[i]->UpdateCoeffs());
    }

    // set up output fld file.
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef =
        m_f->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for (j = 0; j < nFields; ++j)
    {
        for (i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back(m_f->m_fieldPts->GetFieldName(j));

            m_f->m_exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}
}
}
