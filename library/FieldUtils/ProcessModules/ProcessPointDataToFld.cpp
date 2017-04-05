////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessPointDataToFld.cpp
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
//  Description: Load a file of interpolated vales at physical
//  quadrature points and project to a fld file
//
////////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <string>

using namespace std;

#include "ProcessPointDataToFld.h"

#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessPointDataToFld::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "pointdatatofld"),
        ProcessPointDataToFld::create,
        "Given discrete data at quadrature points project them onto an "
        "expansion"
        "basis and output fld file. Requires .pts .xml and .fld files.");

ProcessPointDataToFld::ProcessPointDataToFld(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_requireEquiSpaced = true;

    m_config["setnantovalue"] = ConfigOption(
        false, "NotSet", "reset any nan value to prescribed value");

    if ((f->m_inputfiles.count("pts") == 0))
    {
        cout << endl
             << "A pts input file must be specified for the boundary "
                "extraction module"
             << endl;

        cout
            << "Usage: Fieldconvert -m pointdatatofld file.pts file.xml out.fld"
            << endl;
        exit(3);
    }

    if ((f->m_inputfiles.count("xml") == 0) &&
        (f->m_inputfiles.count("xml.gz") == 0))
    {
        cout << endl
             << "An xml or xml.gz input file must be specified for the "
                "boundary extraction module"
             << endl;
        cout
            << "Usage: Fieldconvert -m pointdatatofld file.pts file.xml out.fld"
            << endl;
        exit(3);
    }
}

ProcessPointDataToFld::~ProcessPointDataToFld()
{
}

void ProcessPointDataToFld::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            cout << "ProcessPointDataToFld: projecting data to expansion..."
                 << endl;
        }
    }

    int i, j;
    bool setnantovalue = false;
    NekDouble defvalue=0.0;

    if (!boost::iequals(m_config["setnantovalue"].as<string>(), "NotSet"))
    {
        setnantovalue = true;
        defvalue      = m_config["setnantovalue"].as<NekDouble>();
    }

    // Check for command line point specification if no .pts file specified
    ASSERTL0(m_f->m_fieldPts != LibUtilities::NullPtsField,
             "No input points found");

    int nFields = m_f->m_fieldPts->GetNFields();
    ASSERTL0(nFields > 0, "No field values provided in input");

    int dim = m_f->m_fieldPts->GetDim();

    // assume one field is already defined from input file.
    m_f->m_exp.resize(nFields);
    for (i = 1; i < nFields; ++i)
    {
        m_f->m_exp[i] = m_f->AppendExpList(0);
    }

    Array<OneD, Array<OneD, NekDouble> > pts;
    m_f->m_fieldPts->GetPts(pts);

    // set any nan values to default value if requested
    if (setnantovalue)
    {
        for (int i = 0; i < pts[0].num_elements(); ++i)
        {
            for (int j = 0; j < nFields; ++j)
            {
                if ((boost::math::isnan)(pts[j + dim][i]))
                {
                    pts[j + dim][i] = defvalue;
                }
            }
        }
    }

    if (m_f->m_fieldPts->m_ptsInfo.count(LibUtilities::eIsEquiSpacedData) != 0)
    {
        int totcoeffs = m_f->m_exp[0]->GetNcoeffs();

        ASSERTL0(pts[0].num_elements() != totcoeffs,
                 "length of points in .pts file is different "
                 "to the number of coefficients in expansion ");

        for (int i = 0; i < nFields; ++i)
        {
            Array<OneD, NekDouble> coeffs = m_f->m_exp[i]->UpdateCoeffs(), tmp;
            int cnt = 0;
            for (int e = 0; e < m_f->m_exp[0]->GetNumElmts(); ++e)
            {
                int ncoeffs = m_f->m_exp[i]->GetExp(e)->GetNcoeffs();
                int offset  = m_f->m_exp[i]->GetCoeff_Offset(e);
                Vmath::Vcopy(ncoeffs, &pts[i + dim][cnt], 1, &coeffs[offset],
                             1);

                m_f->m_exp[i]->GetExp(e)->EquiSpacedToCoeffs(
                    coeffs + offset, tmp = coeffs + offset);
                cnt += ncoeffs;
            }
        }
    }
    else
    {
        int totpoints = m_f->m_exp[0]->GetTotPoints();
        Array<OneD, Array<OneD, NekDouble> > coords(3);

        coords[0] = Array<OneD, NekDouble>(totpoints);
        coords[1] = Array<OneD, NekDouble>(totpoints);
        coords[2] = Array<OneD, NekDouble>(totpoints);

        m_f->m_exp[0]->GetCoords(coords[0], coords[1], coords[2]);

        if (pts[0].num_elements() != totpoints)
        {
            WARNINGL0(false, "length of points in .pts file is different to "
                             "the number of quadrature points in xml file");
            totpoints = min(totpoints, (int)pts[0].num_elements());
        }

        int mismatch = 0;
        for (i = 0; i < totpoints; ++i)
        {
            for (j = 0; j < dim; ++j)
            {
                if (fabs(coords[j][i] - pts[j][i]) > 1e-4)
                {
                    string outstring =
                        "Coordinates do not match within 1e-4: " +
                        boost::lexical_cast<string>(coords[j][i]) + " versus " +
                        boost::lexical_cast<string>(pts[j][i]) + " diff " +
                        boost::lexical_cast<string>(
                            fabs(coords[j][i] - pts[j][i]));
                    ;
                    WARNINGL0(false, outstring);
                    mismatch += 1;
                }
            }

            for (j = 0; j < nFields; ++j)
            {
                m_f->m_exp[j]->SetPhys(i, pts[j + dim][i]);
            }
        }

        if (m_f->m_session->GetComm()->GetRank() == 0)
        {
            cout << endl;
        }

        // forward transform fields
        for (i = 0; i < nFields; ++i)
        {
            m_f->m_exp[i]->FwdTrans_IterPerExp(m_f->m_exp[i]->GetPhys(),
                                               m_f->m_exp[i]->UpdateCoeffs());
        }
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
