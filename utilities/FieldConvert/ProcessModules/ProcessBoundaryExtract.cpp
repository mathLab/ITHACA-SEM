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

    // check for correct input files
    if((f->m_inputfiles.count("xml") == 0)&&(f->m_inputfiles.count("xml.gz") == 0))
    {
        cout << "An xml or xml.gz input file must be specified for the boundary extraction module" << endl;
        exit(3);
    }

    if((f->m_inputfiles.count("fld") == 0)&&(f->m_inputfiles.count("chk") == 0)&&(f->m_inputfiles.count("rst") == 0))
    {
        cout << "A fld or chk or rst input file must be specified for the boundary extraction module" << endl;

        exit(3);
    }

}

ProcessBoundaryExtract::~ProcessBoundaryExtract()
{
}

void ProcessBoundaryExtract::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessBoundaryExtract: Setting up boundary extraction..." << endl;
    }

    // Set up Field options to output boundary fld
    string bvalues =  m_config["bnd"].as<string>();

    if(bvalues.compare("All") == 0)
    {
        Array<OneD, const MultiRegions::ExpListSharedPtr>
            BndExp = m_f->m_exp[0]->GetBndCondExpansions();

        for(int i = 0; i < BndExp.num_elements(); ++i)
        {
            m_f->m_bndRegionsToWrite.push_back(i);
        }
    }
    else
    {
        ASSERTL0(ParseUtils::GenerateOrderedVector(bvalues.c_str(),
                                                   m_f->m_bndRegionsToWrite),"Failed to interpret range string");
    }

    m_f->m_fldToBnd = m_config["fldtoboundary"].m_beenSet;
    m_f->m_addNormals = m_config["addnormals"].m_beenSet;

}

}
}


