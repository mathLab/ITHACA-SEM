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
//  Description: Set up boundary to be extracted when writing fld file.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessBoundaryExtract.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessBoundaryExtract::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "extract"),
        ProcessBoundaryExtract::create, "Extract Boundary field");

ProcessBoundaryExtract::ProcessBoundaryExtract(FieldSharedPtr f) : ProcessModule(f)
{
    // set up dafault values.
    m_config["bnd"] = ConfigOption(false,"All","Boundary to be extracted");
    m_config["fldtoboundary"] = ConfigOption(true,"NotSet","Extract fld values to boundary");
    m_config["addnormals"] = ConfigOption(true,"NotSet","Add normals to output");


    f->m_writeBndFld = true;
    f->m_declareExpansionAsContField = true;

}

ProcessBoundaryExtract::~ProcessBoundaryExtract()
{
}

void ProcessBoundaryExtract::Process(po::variables_map &vm)
{
    Timer timer;

    if (m_f->m_verbose)
    {
        if(m_f->m_comm->GetRank() == 0)
        {
            cout << "Process Boundary Extract: Setting up boundary extraction..." << endl;
            timer.Start();
        }
    }


    // check for correct input files
    if((m_f->m_inputfiles.count("xml") == 0)&&(m_f->m_inputfiles.count("xml.gz") == 0))
    {
        cout << "An xml or xml.gz input file must be specified for the boundary extraction module" << endl;
        exit(3);
    }

    if((m_f->m_inputfiles.count("fld") == 0)&&(m_f->m_inputfiles.count("chk") == 0)&&(m_f->m_inputfiles.count("rst") == 0))
    {
        cout << "A fld or chk or rst input file must be specified for the boundary extraction module" << endl;

        exit(3);
    }

    // Set up Field options to output boundary fld
    string bvalues =  m_config["bnd"].as<string>();

    if(boost::iequals(bvalues,"All"))
    {
        int numBndExp = 0; 

        SpatialDomains::BoundaryConditions bcs(m_f->m_session,
                                               m_f->m_exp[0]->GetGraph());
        const SpatialDomains::BoundaryRegionCollection bregions  =
                                                    bcs.GetBoundaryRegions();

        SpatialDomains::BoundaryRegionCollection::const_iterator breg_it;
        for(breg_it = bregions.begin(); breg_it != bregions.end();
            ++breg_it)
        {
            numBndExp = max(numBndExp,breg_it->first);
        }
        // assuming all boundary regions are consecutive number if
        // regions is one more tham maximum id
        numBndExp++; 

        // not all partitions in parallel touch all boundaries so 
        // find maximum number of boundaries
        m_f->m_session->GetComm()->AllReduce(numBndExp,LibUtilities::ReduceMax);

        // THis presumes boundary regions are numbered consecutively
        for(int i = 0; i < numBndExp; ++i)
        {
            m_f->m_bndRegionsToWrite.push_back(i);
        }
    }
    else
    {
        ASSERTL0(ParseUtils::GenerateOrderedVector(bvalues.c_str(),
                                                   m_f->m_bndRegionsToWrite),
                 "Failed to interpret bnd values string");
    }

    m_f->m_fldToBnd = m_config["fldtoboundary"].m_beenSet;
    m_f->m_addNormals = m_config["addnormals"].m_beenSet;

    if(m_f->m_verbose)
    {
        if(m_f->m_comm->GetRank() == 0)
        {
            timer.Stop();
            NekDouble cpuTime = timer.TimePerTest(1);
            
            stringstream ss;
            ss << cpuTime << "s";
            cout << "Process Boundary Extract CPU Time: " << setw(8) << left
                 << ss.str() << endl;
            cpuTime = 0.0;
        }
    }
}

}
}


