////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputInfo.cpp
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
//  Description: Generates a Nektar++ info XML file.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>
#include <boost/format.hpp>

#include <LibUtilities/BasicUtils/FieldIOXml.h>
#include <SpatialDomains/MeshPartition.h>

#include "OutputInfo.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey OutputInfo::m_className =
    GetModuleFactory().RegisterCreatorFunction(ModuleKey(eOutputModule, "info"),
                                               OutputInfo::create,
                                               "Writes an Info file.");

OutputInfo::OutputInfo(FieldSharedPtr f) : OutputModule(f)
{
    m_config["nparts"] = ConfigOption(
        false, "NotSet",
        "Number of partitions over which to create the info file");
}

OutputInfo::~OutputInfo()
{
}

void OutputInfo::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();

    // partition mesh
    ASSERTL0(m_config["nparts"].as<string>().compare("NotSet") != 0,
             "Need to specify nparts for info output");
    int nparts = m_config["nparts"].as<int>();

    // Input/output file
    LibUtilities::CommSharedPtr c = m_f->m_comm;
    std::shared_ptr<LibUtilities::FieldIOXml> fldXml =
        std::static_pointer_cast<LibUtilities::FieldIOXml>(
            LibUtilities::GetFieldIOFactory().CreateInstance("Xml", c, true));

    // open file and setup meta data.
    fs::path pinfilename(filename);
    std::vector<std::string> filenames;
    std::vector<std::vector<unsigned int> > ElementIDs;

    for (int p = 0; p < nparts; ++p)
    {
        boost::format pad("P%1$07d.%2$s");
        pad % p % "fld";
        std::string s = pad.str();

        fs::path fullpath              = pinfilename / s;
        string fname                   = LibUtilities::PortablePath(fullpath);
        LibUtilities::DataSourceSharedPtr dataSource =
            LibUtilities::XmlDataSource::create(fname);

        std::vector<LibUtilities::FieldDefinitionsSharedPtr> fielddefs;
        std::vector<unsigned int>              PartElmtIDs;

        // read in header of partition if it exists
        fldXml->ImportFieldDefs(dataSource,fielddefs,false);

        // create ElmenetIDs list then use
        for(int i = 0; i < fielddefs.size(); ++i)
        {
            for(int j = 0; j < fielddefs[i]->m_elementIDs.size(); ++j)
            {
                PartElmtIDs.push_back(fielddefs[i]->m_elementIDs[j]);
            }
        }

        ElementIDs.push_back(PartElmtIDs);
        filenames.push_back(s);
    }

    // Write the Info.xml file
    string infofile =
        LibUtilities::PortablePath(pinfilename / fs::path("Info.xml"));

    fldXml->WriteMultiFldFileIDs(infofile,filenames, ElementIDs);
}

}
}
