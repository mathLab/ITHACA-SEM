////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessBoundaryExtract.cpp
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
//  Description: Set up boundary to be extracted when writing fld file.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessBoundaryExtract.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey ProcessBoundaryExtract::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "extract"),
        ProcessBoundaryExtract::create,
        "Extract Boundary field");

ProcessBoundaryExtract::ProcessBoundaryExtract(FieldSharedPtr f)
    : ProcessModule(f)
{
    // set up dafault values.
    m_config["bnd"] = ConfigOption(false, "All", "Boundary to be processed");
    m_config["addnormals"] =
        ConfigOption(true, "0", "Add normals to output");

    f->m_writeBndFld                 = true;
    f->m_declareExpansionAsContField = true;
    f->m_requireBoundaryExpansion    = true;
}

ProcessBoundaryExtract::~ProcessBoundaryExtract()
{
}

void ProcessBoundaryExtract::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    m_f->m_addNormals = m_config["addnormals"].as<bool>();

    // Set up Field options to output boundary fld
    string bvalues = m_config["bnd"].as<string>();

    vector<unsigned int> bndRegions;
    if (boost::iequals(bvalues, "All"))
    {
        int numBndExp = 0;

        SpatialDomains::BoundaryConditions bcs(m_f->m_session,
                                               m_f->m_exp[0]->GetGraph());
        const SpatialDomains::BoundaryRegionCollection bregions =
            bcs.GetBoundaryRegions();

        for (auto &breg_it : bregions)
        {
            numBndExp = max(numBndExp, breg_it.first);
        }
        // assuming all boundary regions are consecutive number if
        // regions is one more than maximum id
        numBndExp++;

        // not all partitions in parallel touch all boundaries so
        // find maximum number of boundaries
        m_f->m_comm->AllReduce(numBndExp, LibUtilities::ReduceMax);

        // THis presumes boundary regions are numbered consecutively
        for (int i = 0; i < numBndExp; ++i)
        {
            bndRegions.push_back(i);
        }
    }
    else
    {
        ASSERTL0(ParseUtils::GenerateVector(bvalues, bndRegions),
                 "Failed to interpret bnd values string");
    }

    if(m_f->m_bndRegionsToWrite.size())
    {
        // This was already called. Just check if the bnd option is the same
        ASSERTL0(m_f->m_bndRegionsToWrite == bndRegions,
                "Incompatible bnd parameters.");
    }
    else
    {
        m_f->m_bndRegionsToWrite = bndRegions;

        if (m_f->m_exp[0]->GetNumElmts() != 0)
        {
            for (int i = 0; i < m_f->m_exp.size(); ++i)
            {
                m_f->m_exp[i]->FillBndCondFromField();
            }
        }
    }
}
}
}
