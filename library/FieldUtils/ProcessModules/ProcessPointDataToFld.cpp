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

#include <boost/core/ignore_unused.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessPointDataToFld.h"

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
        "basis and output fld file. Requires frompts and .xml and .fld files.");

ProcessPointDataToFld::ProcessPointDataToFld(FieldSharedPtr f)
    : ProcessModule(f)
{
    m_config["setnantovalue"] = ConfigOption(
        false, "NotSet", "reset any nan value to prescribed value");

    m_config["frompts"] = ConfigOption(
        false, "NotSet", "Pts file from which to interpolate field");
}

ProcessPointDataToFld::~ProcessPointDataToFld()
{
}

void ProcessPointDataToFld::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    int i, j;
    bool setnantovalue = false;
    NekDouble defvalue=0.0;

    if (!boost::iequals(m_config["setnantovalue"].as<string>(), "NotSet"))
    {
        setnantovalue = true;
        defvalue      = m_config["setnantovalue"].as<NekDouble>();
    }

    // Check for command line point specification if no .pts file specified
    // Load pts file
    LibUtilities::PtsFieldSharedPtr fieldPts;
    ASSERTL0( m_config["frompts"].as<string>().compare("NotSet") != 0,
            "ProcessInterpPointDataToFld requires frompts parameter");
    string inFile = m_config["frompts"].as<string>().c_str();
    LibUtilities::CommSharedPtr  c     =
            LibUtilities::GetCommFactory().CreateInstance("Serial", 0, 0);
    LibUtilities::PtsIOSharedPtr ptsIO =
            MemoryManager<LibUtilities::PtsIO>::AllocateSharedPtr(c);
    ptsIO->Import(inFile, fieldPts);

    int nFields = fieldPts->GetNFields();
    ASSERTL0(nFields > 0, "No field values provided in input");

    int dim = fieldPts->GetDim();

    // assume one field is already defined from input file.
    ASSERTL0(m_f->m_numHomogeneousDir == 0,
        "ProcessInterpPointDataToFld does not support homogeneous expansion");

    m_f->m_exp.resize(nFields);
    for (i = 1; i < nFields; ++i)
    {
        m_f->m_exp[i] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
    }
    Array<OneD, Array<OneD, NekDouble> > pts;
    fieldPts->GetPts(pts);

    // set any nan values to default value if requested
    if (setnantovalue)
    {
        for (int i = 0; i < pts[0].size(); ++i)
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

    if (fieldPts->m_ptsInfo.count(LibUtilities::eIsEquiSpacedData) != 0)
    {
        int totcoeffs = m_f->m_exp[0]->GetNcoeffs();

        ASSERTL0(pts[0].size() != totcoeffs,
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
            m_f->m_exp[i]->BwdTrans(m_f->m_exp[i]->GetCoeffs(),
                                    m_f->m_exp[i]->UpdatePhys());
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

        if (pts[0].size() != totpoints)
        {
            WARNINGL0(false, "length of points in .pts file is different to "
                             "the number of quadrature points in xml file");
            totpoints = min(totpoints, (int)pts[0].size());
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

        // forward transform fields
        for (i = 0; i < nFields; ++i)
        {
            m_f->m_exp[i]->FwdTrans_IterPerExp(m_f->m_exp[i]->GetPhys(),
                                               m_f->m_exp[i]->UpdateCoeffs());
        }
    }

    // save field names
    for (int j = 0; j < fieldPts->GetNFields(); ++j)
    {
        m_f->m_variables.push_back(fieldPts->GetFieldName(j));
    }
}
}
}
