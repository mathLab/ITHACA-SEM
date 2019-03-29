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
//  Description: Generates a Nektar++ info XML file.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "OutputInfo.h"

#include <LibUtilities/BasicUtils/FieldIOXml.h>
#include <SpatialDomains/MeshPartition.h>
#include <boost/format.hpp>

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
    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();
    int i = 0;

    // partition mesh
    ASSERTL0(m_config["nparts"].as<string>().compare("NotSet") != 0,
             "Need to specify nparts for info output");
    const int nparts = m_config["nparts"].as<int>();

    std::vector<std::string> files;
    // load .xml ending
    for (auto &x : m_f->m_inputfiles["xml"]) {
        files.push_back(x);
    }

    // load any .xml.gz endings
    for (auto &x: m_f->m_inputfiles["xml.gz"])
    {
        files.push_back(x);
    }

    ASSERTL0(m_f->m_comm->GetSize() == 1,
             "OutputInfo module should be run in serial.");

    // Default partitioner to use is Scotch. Override default with
    // command-line flags if they are set.
    string vPartitionerName = "Scotch";
    if (m_f->m_session->DefinesCmdLineArgument("use-metis"))
    {
        vPartitionerName = "Metis";
    }
    if (m_f->m_session->DefinesCmdLineArgument("use-scotch"))
    {
        vPartitionerName = "Scotch";
    }

    // Construct mesh partitioning.
    SpatialDomains::MeshPartitionSharedPtr meshPartition =
        SpatialDomains::GetMeshPartitionFactory().CreateInstance(
            vPartitionerName, m_f->m_session, m_f->m_graph->GetMeshDimension(),
            m_f->m_graph->CreateMeshEntities(),
            m_f->m_graph->CreateCompositeDescriptor());

    meshPartition->PartitionMesh(nparts, true);

    // get hold of local partition ids
    std::vector<std::vector<unsigned int> > ElementIDs(nparts);

    // Populate the list of element ID lists from all processes
    for (i = 0; i < nparts; ++i)
    {
        std::vector<unsigned int> tmp;
        meshPartition->GetElementIDs(i, tmp);
        ElementIDs[i] = tmp;
    }

    // Set up output names
    std::vector<std::string> filenames;
    for (int i = 0; i < nparts; ++i)
    {
        boost::format pad("P%1$07d.fld");
        pad % i;
        filenames.push_back(pad.str());
    }

    // Write the output file
    LibUtilities::CommSharedPtr c = m_f->m_comm;
    std::shared_ptr<LibUtilities::FieldIOXml> fldXml =
        std::static_pointer_cast<LibUtilities::FieldIOXml>(
            LibUtilities::GetFieldIOFactory().CreateInstance("Xml", c, true));
    fldXml->WriteMultiFldFileIDs(filename, filenames, ElementIDs);
}
}
}
