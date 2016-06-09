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
#include <LibUtilities/BasicUtils/MeshPartition.h>
#include <boost/format.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey OutputInfo::m_className = GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eOutputModule, "info"), OutputInfo::create,
        "Writes an Info file.");

OutputInfo::OutputInfo(FieldSharedPtr f) : OutputModule(f)
{
}

OutputInfo::~OutputInfo()
{
}

void OutputInfo::Process(po::variables_map &vm)
{
    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();
    int i;

    if (m_f->m_verbose)
    {
        cout << "OutputInfo: Writing Info file..." << endl;
    }

    // partition mesh
    ASSERTL0(vm.count("nprocs") > 0,
             "--nprocs nust be specified with info output");

    int nprocs = vm["nprocs"].as<int>();

    LibUtilities::CommSharedPtr vComm = boost::shared_ptr<FieldConvertComm>(
            new FieldConvertComm(0, NULL, nprocs,0));
    vComm->SplitComm(1,nprocs);

    // define new session with psuedo parallel communicator
    string xml_ending    = "xml";
    string xml_gz_ending = "xml.gz";

    std::vector<std::string> files;
    // load .xml ending
    for (int i = 0; i < m_f->m_inputfiles[xml_ending].size(); ++i)
    {
        files.push_back(m_f->m_inputfiles[xml_ending][i]);
    }

    // load any .xml.gz endings
    for (int j =0; j < m_f->m_inputfiles[xml_gz_ending].size(); ++j)
    {
        files.push_back(m_f->m_inputfiles[xml_gz_ending][j]);
    }

    LibUtilities::SessionReaderSharedPtr vSession =
                    boost::shared_ptr<LibUtilities::SessionReader>(
                            new LibUtilities::SessionReader(0,0,files,vComm));
    vSession->SetUpXmlDoc();

    // Default partitioner to use is Metis. Use Scotch as default
    // if it is installed. Override default with command-line flags
    // if they are set.
    string vPartitionerName = "Metis";
    if (LibUtilities::GetMeshPartitionFactory().ModuleExists("Scotch"))
    {
        vPartitionerName = "Scotch";
    }
    if (vSession->DefinesCmdLineArgument("use-metis"))
    {
        vPartitionerName = "Metis";
    }
    if (vSession->DefinesCmdLineArgument("use-scotch"))
    {
        vPartitionerName = "Scotch";
    }

    LibUtilities::MeshPartitionSharedPtr vMeshPartition =
        LibUtilities::GetMeshPartitionFactory().CreateInstance(
                                                 vPartitionerName, vSession);

    vMeshPartition->PartitionMesh(nprocs, true);

    // get hold of local partition ids
    std::vector<std::vector<unsigned int> > ElementIDs(nprocs);

    // Populate the list of element ID lists from all processes
    for (i = 0; i < nprocs; ++i)
    {
        std::vector<unsigned int> tmp;
        vMeshPartition->GetElementIDs(i,tmp);
        ElementIDs[i] = tmp;
    }

    // Set up output names
    std::vector<std::string> filenames;
    for(int i = 0; i < nprocs; ++i)
    {
        boost::format pad("P%1$07d.fld");
        pad % i;
        filenames.push_back(pad.str());
    }

    // Write the output file
    m_f->m_fld->WriteMultiFldFileIDs(filename,filenames, ElementIDs);
}

}
}
