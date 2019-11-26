////////////////////////////////////////////////////////////////////////////////
//
//  File: InputXml.cpp
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
//  Description: Read xml file and set up expansions
//
////////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>
#include <string>
using namespace std;

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/Timer.h>

#include "InputXml.h"
using namespace Nektar;

namespace Nektar
{
namespace FieldUtils
{

ModuleKey InputXml::m_className[5] = {
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "xml"), InputXml::create, "Reads Xml file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "xml.gz"), InputXml::create, "Reads Xml file."),
};

/**
 * @brief Set up InputXml object.
 *
 */
InputXml::InputXml(FieldSharedPtr f) : InputModule(f)
{
    m_allowedFiles.insert("xml");
    m_allowedFiles.insert("xml.gz");
}

/**
 *
 */
InputXml::~InputXml()
{
}

/**
 *
 */
void InputXml::Process(po::variables_map &vm)
{
    LibUtilities::Timer timerpart;
    if (m_f->m_verbose)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            timerpart.Start();
        }
    }

    // check for multiple calls to inputXml due to split xml
    // files. If so just return
    if (m_f->m_graph)
    {
        return;
    }

    string xml_ending    = "xml";
    string xml_gz_ending = "xml.gz";

    std::vector<std::string> files;
    // load .xml ending
    for (int i = 0; i < m_f->m_inputfiles[xml_ending].size(); ++i)
    {
        files.push_back(m_f->m_inputfiles[xml_ending][i]);
    }

    // load any .xml.gz endings
    for (int j = 0; j < m_f->m_inputfiles[xml_gz_ending].size(); ++j)
    {
        files.push_back(m_f->m_inputfiles[xml_gz_ending][j]);
    }

    SpatialDomains::DomainRangeShPtr rng = SpatialDomains::NullDomainRangeShPtr;

    // define range to process output
    if (vm.count("range"))
    {
        vector<NekDouble> values;
        ASSERTL0(ParseUtils::GenerateVector(vm["range"].as<string>(), values),
                 "Failed to interpret range string");

        ASSERTL0(values.size() > 1, "Do not have minimum values of xmin,xmax");
        ASSERTL0(values.size() % 2 == 0,
                 "Do not have an even number of range values");

        int nvalues = values.size() / 2;
        rng = MemoryManager<SpatialDomains::DomainRange>::AllocateSharedPtr();

        rng->m_doZrange   = false;
        rng->m_doYrange   = false;
        rng->m_checkShape = false;

        switch (nvalues)
        {
            case 3:
                rng->m_doZrange = true;
                rng->m_zmin     = values[4];
                rng->m_zmax     = values[5];
                /* Falls through. */
            case 2:
                rng->m_doYrange = true;
                rng->m_ymin     = values[2];
                rng->m_ymax     = values[3];
                /* Falls through. */
            case 1:
                rng->m_doXrange = true;
                rng->m_xmin     = values[0];
                rng->m_xmax     = values[1];
                break;
            default:
                NEKERROR(ErrorUtil::efatal, "too many values specfied in range");
        }
    }

    // define range to only take a single shape.
    if (vm.count("onlyshape"))
    {
        if (rng == SpatialDomains::NullDomainRangeShPtr)
        {
            rng =
                MemoryManager<SpatialDomains::DomainRange>::AllocateSharedPtr();
            rng->m_doXrange = false;
            rng->m_doYrange = false;
            rng->m_doZrange = false;
        }

        rng->m_checkShape = true;

        string shapematch = boost::to_upper_copy(vm["onlyshape"].as<string>());
        int i;
        for (i = 0; i < LibUtilities::SIZE_ShapeType; ++i)
        {
            string shapeval = LibUtilities::ShapeTypeMap[i];
            boost::to_upper(shapeval);
            if (shapematch.compare(shapeval) == 0)
            {
                rng->m_shapeType = (LibUtilities::ShapeType)i;
                break;
            }
        }
        ASSERTL0(i != LibUtilities::SIZE_ShapeType,
                 "Failed to find shape type in -onlyshape command line "
                 "argument");
    }

    // Set up command lines options that will be passed through to SessionReader
    vector<string> cmdArgs;
    cmdArgs.push_back("FieldConvert");

    if (m_f->m_verbose)
    {
        cmdArgs.push_back("--verbose");
    }

    if (vm.count("part-only"))
    {
        cmdArgs.push_back("--part-only");
        cmdArgs.push_back(
            boost::lexical_cast<string>(vm["part-only"].as<int>()));
    }

    if (vm.count("part-only-overlapping"))
    {
        cmdArgs.push_back("--part-only-overlapping");
        cmdArgs.push_back(
            boost::lexical_cast<string>(vm["part-only-overlapping"].as<int>()));
    }

    if (vm.count("npz"))
    {
        cmdArgs.push_back("--npz");
        cmdArgs.push_back(boost::lexical_cast<string>(vm["npz"].as<int>()));
    }

    int argc          = cmdArgs.size();
    const char **argv = new const char *[argc];
    for (int i = 0; i < argc; ++i)
    {
        argv[i] = cmdArgs[i].c_str();
    }

    m_f->m_session = LibUtilities::SessionReader::CreateInstance(
        argc, (char **)argv, files, m_f->m_comm);

    if (vm.count("nparts"))
    {
        // make sure have pre-partitioned mesh for nparts option
        ASSERTL0(boost::icontains(files[0],"_xml"),
                 "Expect the mesh to have been pre-partitioned when "
                 " using the\"--nparts\" option. Please use \"--part-only\" "
                 "option to prepartition xml file.");
    }

    // Free up memory.
    delete[] argv;

    if (m_f->m_verbose)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            timerpart.Stop();
            NekDouble cpuTime = timerpart.TimePerTest(1);

            stringstream ss;
            ss << cpuTime << "s";
            cout << "\t InputXml session reader CPU Time: " << setw(8) << left
                 << ss.str() << endl;
            timerpart.Start();
        }
    }

    m_f->m_graph = SpatialDomains::MeshGraph::Read(m_f->m_session, rng);

    if (m_f->m_verbose)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            timerpart.Stop();
            NekDouble cpuTime = timerpart.TimePerTest(1);

            stringstream ss;
            ss << cpuTime << "s";
            cout << "\t InputXml mesh graph setup  CPU Time: " << setw(8)
                 << left << ss.str() << endl;
            timerpart.Start();
        }
    }
}
}
}
