////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessCombineAvg.cpp
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
//  Description: Combines two fld files containing average fields.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessCombineAvg.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessCombineAvg::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "combineAvg"),
        ProcessCombineAvg::create,
        "combine two fields containing averages (and possibly Reynolds "
        "stresses). Must specify fromfld.");

ProcessCombineAvg::ProcessCombineAvg(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["fromfld"] =
        ConfigOption(false, "NotSet", "Fld file form which to add field");
}

ProcessCombineAvg::~ProcessCombineAvg()
{
}

void ProcessCombineAvg::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    ASSERTL0(m_config["fromfld"].as<string>().compare("NotSet") != 0,
             "Need to specify fromfld=file.fld ");

    int nfields  = m_f->m_variables.size();
    int nq       = m_f->m_exp[0]->GetTotPoints();
    int expdim   = m_f->m_graph->GetMeshDimension();
    int spacedim = expdim;
    if ((m_f->m_numHomogeneousDir) == 1 || (m_f->m_numHomogeneousDir) == 2)
    {
        spacedim += m_f->m_numHomogeneousDir;
    }

    // Allocate storage for new field and correction (for Reynolds stress)
    Array<OneD, Array<OneD, NekDouble> > fromPhys(nfields);
    Array<OneD, Array<OneD, NekDouble> > correction(nfields);
    for (int j = 0; j < nfields; ++j)
    {
        fromPhys[j]   = Array<OneD, NekDouble>(nq, 0.0);
        correction[j] = Array<OneD, NekDouble>(nq, 0.0);
    }

    string fromfld           = m_config["fromfld"].as<string>();
    FieldSharedPtr fromField = std::shared_ptr<Field>(new Field());
    LibUtilities::FieldMetaDataMap fromFieldMetaDataMap;

    // Set up ElementGIDs in case of parallel processing
    Array<OneD, int> ElementGIDs(m_f->m_exp[0]->GetExpSize());
    for (int i = 0; i < m_f->m_exp[0]->GetExpSize(); ++i)
    {
        ElementGIDs[i] = m_f->m_exp[0]->GetExp(i)->GetGeom()->GetGlobalID();
    }
    // Import fromfld file
    m_f->FieldIOForFile(fromfld)->Import(
        fromfld, fromField->m_fielddef, fromField->m_data, fromFieldMetaDataMap,
        ElementGIDs);
    ASSERTL0(fromField->m_fielddef[0]->m_fields.size() == nfields,
             "Mismatch in number of fields");
    // Extract data to fromPhys
    for (int j = 0; j < nfields; ++j)
    {
        ASSERTL0(fromField->m_fielddef[0]->m_fields[j] ==
                     m_f->m_variables[j],
                 "Field names do not match.");

        // load new field (overwrite m_f->m_exp coeffs for now)
        for (int i = 0; i < fromField->m_data.size(); ++i)
        {
            m_f->m_exp[j]->ExtractDataToCoeffs(
                fromField->m_fielddef[i], fromField->m_data[i],
                m_f->m_variables[j],
                m_f->m_exp[j]->UpdateCoeffs());
        }
        m_f->m_exp[j]->BwdTrans(m_f->m_exp[j]->GetCoeffs(), fromPhys[j]);
    }

    // Load number of samples in each file
    ASSERTL0(m_f->m_fieldMetaDataMap.count("NumberOfFieldDumps") != 0,
             "Missing NumberOfFieldDumps metadata.");
    ASSERTL0(fromFieldMetaDataMap.count("NumberOfFieldDumps") != 0,
             "Missing NumberOfFieldDumps metadata.");
    string s_num;
    s_num  = m_f->m_fieldMetaDataMap["NumberOfFieldDumps"];
    int na = atoi(s_num.c_str());
    s_num  = fromFieldMetaDataMap["NumberOfFieldDumps"];
    int nb = atoi(s_num.c_str());

    // Look for Reynolds stresses
    int stress = -1;
    for (int j = 0; j < nfields; ++j)
    {
        if (m_f->m_variables[j] == "uu")
        {
            stress = j;
            break;
        }
    }

    // Calculate correction for Reynolds stresses
    if (stress != -1)
    {
        Array<OneD, NekDouble> tmp(nq, 0.0);
        int n = stress;
        // Follow same numbering as FilterReynoldsStresses
        for (int i = 0; i < spacedim; ++i)
        {
            for (int j = i; j < spacedim; ++j, ++n)
            {
                // correction is zero for averages and
                //      = (\bar{x_a}-\bar{x_b})*(\bar{y_a}-\bar{y_b})*na*nb/N
                //      for Reynolds stresses
                NekDouble fac = ((NekDouble)(na * nb)) / ((NekDouble)(na + nb));
                Vmath::Vsub(nq, m_f->m_exp[i]->GetPhys(), 1, fromPhys[i], 1,
                            correction[n], 1);
                Vmath::Vsub(nq, m_f->m_exp[j]->GetPhys(), 1, fromPhys[j], 1,
                            tmp, 1);
                Vmath::Vmul(nq, correction[n], 1, tmp, 1, correction[n], 1);
                Vmath::Smul(nq, fac, correction[n], 1, correction[n], 1);
            }
        }
    }
    // Combine fields
    for (int j = 0; j < nfields; ++j)
    {
        // The new value is: (x_a*na + x_b*nb + correction)/N
        Vmath::Smul(nq, 1.0 * na, m_f->m_exp[j]->GetPhys(), 1,
                    m_f->m_exp[j]->UpdatePhys(), 1);
        Vmath::Svtvp(nq, 1.0 * nb, fromPhys[j], 1, m_f->m_exp[j]->GetPhys(), 1,
                     m_f->m_exp[j]->UpdatePhys(), 1);
        Vmath::Vadd(nq, m_f->m_exp[j]->GetPhys(), 1, correction[j], 1,
                    m_f->m_exp[j]->UpdatePhys(), 1);
        Vmath::Smul(nq, 1.0 / (na + nb), m_f->m_exp[j]->GetPhys(), 1,
                    m_f->m_exp[j]->UpdatePhys(), 1);

        m_f->m_exp[j]->FwdTrans_IterPerExp(m_f->m_exp[j]->GetPhys(),
                                           m_f->m_exp[j]->UpdateCoeffs());
    }

    // Update metadata
    m_f->m_fieldMetaDataMap["NumberOfFieldDumps"] =
        boost::lexical_cast<std::string>(na + nb);
    NekDouble t0      = -1;
    NekDouble finTime = -1;
    if (m_f->m_fieldMetaDataMap.count("InitialTime"))
    {
        string s_t  = m_f->m_fieldMetaDataMap["InitialTime"];
        NekDouble t = atof(s_t.c_str());

        t0 = t;
    }
    if (fromFieldMetaDataMap.count("InitialTime"))
    {
        string s_t  = fromFieldMetaDataMap["InitialTime"];
        NekDouble t = atof(s_t.c_str());

        if (t0 == -1)
        {
            t0 = t;
        }
        else
        {
            t0 = std::min(t0, t);
        }
    }
    if (m_f->m_fieldMetaDataMap.count("FinalTime"))
    {
        string s_t  = m_f->m_fieldMetaDataMap["FinalTime"];
        NekDouble t = atof(s_t.c_str());

        finTime = std::max(t0, t);
    }
    if (fromFieldMetaDataMap.count("FinalTime"))
    {
        string s_t  = fromFieldMetaDataMap["FinalTime"];
        NekDouble t = atof(s_t.c_str());

        finTime = std::max(t0, t);
    }
    if (t0 != -1)
    {
        m_f->m_fieldMetaDataMap["InitialTime"] =
            boost::lexical_cast<std::string>(t0);
    }
    if (finTime != -1)
    {
        m_f->m_fieldMetaDataMap["FinalTime"] =
            boost::lexical_cast<std::string>(finTime);
    }

}
}
}
