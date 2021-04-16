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

#include <boost/core/ignore_unused.hpp>
#include <boost/geometry.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include <FieldUtils/Interpolator.h>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/CsvIO.h>

#include "ProcessInterpPointDataToFld.h"

namespace bg  = boost::geometry;
namespace bgi = boost::geometry::index;

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
    m_config["frompts"] = ConfigOption(
        false, "NotSet", "Pts file from which to interpolate field");

    m_config["interpcoord"] =
        ConfigOption(false, "-1", "coordinate id to use for interpolation");
}

ProcessInterpPointDataToFld::~ProcessInterpPointDataToFld()
{
}

void ProcessInterpPointDataToFld::Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);

    int i, j;
    LibUtilities::PtsFieldSharedPtr fieldPts;
    // Load pts file
    ASSERTL0( m_config["frompts"].as<string>().compare("NotSet") != 0,
            "ProcessInterpPointDataToFld requires frompts parameter");
    string inFile = m_config["frompts"].as<string>().c_str();

    int totpoints = m_f->m_exp[0]->GetTotPoints();


    Array<OneD, Array<OneD, NekDouble> > intFields(3);
    for (int i = 0; i < 3; ++i)
    {
        intFields[i] = Array<OneD, NekDouble>(totpoints, 0.);
    }
    m_f->m_exp[0]->GetCoords(intFields[0], intFields[1], intFields[2]);

    if (boost::filesystem::path(inFile).extension() == ".pts")
    {
        LibUtilities::PtsIOSharedPtr ptsIO =
            MemoryManager<LibUtilities::PtsIO>::AllocateSharedPtr(m_f->m_comm);

        ptsIO->Import(inFile, fieldPts);
    }
    else if (boost::filesystem::path(inFile).extension() == ".csv")
    {
        LibUtilities::CsvIOSharedPtr csvIO =
            MemoryManager<LibUtilities::CsvIO>::AllocateSharedPtr(m_f->m_comm);

        LibUtilities::DomainRangeShPtr Range = MemoryManager<LibUtilities::DomainRange>::
            AllocateSharedPtr();
        
        NekDouble vmax, vmin, margin = 0.1;
        vmax = intFields[0][Vmath::Imax(totpoints, intFields[0],1)];
        vmin = intFields[0][Vmath::Imin(totpoints, intFields[0],1)];
        Range->m_xmax = (vmax - vmin)*margin + vmax;
        Range->m_xmin =-(vmax - vmin)*margin + vmin;

        vmax = intFields[1][Vmath::Imax(totpoints, intFields[1],1)];
        vmin = intFields[1][Vmath::Imin(totpoints, intFields[1],1)];
        Range->m_ymax = (vmax - vmin)*margin + vmax;
        Range->m_ymin =-(vmax - vmin)*margin + vmin;

        vmax = intFields[2][Vmath::Imax(totpoints, intFields[2],1)];
        vmin = intFields[2][Vmath::Imin(totpoints, intFields[2],1)];
        Range->m_zmax = (vmax - vmin)*margin + vmax;
        Range->m_zmin =-(vmax - vmin)*margin + vmin;
        
        csvIO->Import(inFile, fieldPts, LibUtilities::NullFieldMetaDataMap, Range);
    }
    else
    {
        ASSERTL0(false, "unknown frompts file type");
    }

    int nFields = fieldPts->GetNFields();
    ASSERTL0(nFields > 0, "No field values provided in input");

    // Define new expansions.
    ASSERTL0(m_f->m_numHomogeneousDir == 0,
        "ProcessInterpPointDataToFld does not support homogeneous expansion");

    m_f->m_exp.resize(nFields);
    for (i = 1; i < nFields; ++i)
    {
        m_f->m_exp[i] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
    }

    Array<OneD, Array<OneD, NekDouble> > intFields1(3+nFields);

    for (int i = 0; i < 3; ++i)
    {
        intFields1[i] = intFields[i];
    }

    // Declare space for interpolated fields
    for (int i = 3; i < 3 + nFields; ++i)
    {
        intFields1[i] = Array<OneD, NekDouble>(totpoints);
    }

    LibUtilities::PtsFieldSharedPtr outPts =
        MemoryManager<LibUtilities::PtsField>::AllocateSharedPtr(3, intFields1);

    int coord_id = m_config["interpcoord"].as<int>();
    ASSERTL0(coord_id <= static_cast<int>(outPts->GetDim()) - 1,
             "interpcoord is bigger than the Pts files dimension");

    Interpolator interp(LibUtilities::eNoMethod, coord_id);

    if (m_f->m_verbose && m_f->m_comm->TreatAsRankZero())
    {
        interp.SetProgressCallback(
            &ProcessInterpPointDataToFld::PrintProgressbar, this);
    }
    interp.Interpolate(fieldPts, outPts);
    if (m_f->m_verbose && m_f->m_comm->TreatAsRankZero())
    {
        cout << endl;
    }

    for (i = 0; i < totpoints; ++i)
    {
        for (j = 0; j < nFields; ++j)
        {
            m_f->m_exp[j]->SetPhys(i, outPts->GetPointVal(3 + j, i));
        }
    }

    // forward transform fields
    for (i = 0; i < nFields; ++i)
    {
        m_f->m_exp[i]->FwdTrans_IterPerExp(m_f->m_exp[i]->GetPhys(),
                                           m_f->m_exp[i]->UpdateCoeffs());
    }

    // save field names
    for (int j = 0; j < fieldPts->GetNFields(); ++j)
    {
        m_f->m_variables.push_back(fieldPts->GetFieldName(j));
    }
}
}
}
