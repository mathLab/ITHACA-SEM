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
//  Description: Read xml file and set up expansions
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <iomanip>
using namespace std;

#include <LibUtilities/BasicUtils/Timer.h>

#include "InputXml.h"
using namespace Nektar;

static std::string npts = LibUtilities::SessionReader::RegisterCmdLineArgument(
                "NumberOfPoints","n","Define number of points to dump output");

namespace Nektar
{
namespace Utilities
{

ModuleKey InputXml::m_className[5] = {
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "xml"), InputXml::create,
        "Reads Xml file."),
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "xml.gz"), InputXml::create,
        "Reads Xml file."),
};


/**
 * @brief Set up InputXml object.
 *
 */
InputXml::InputXml(FieldSharedPtr f) : InputModule(f)
{
    m_allowedFiles.insert("xml");
    m_allowedFiles.insert("xml.gz");
    m_allowedFiles.insert("fld"); // these files could be allowed with xml files
    m_allowedFiles.insert("chk");
    m_allowedFiles.insert("rst");
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
    Timer timer, timerpart;

    //check for multiple calls to inputXml due to split xml
    //files. If so just return
    int expsize = m_f->m_exp.size();
    m_f->m_comm->AllReduce(expsize,LibUtilities::ReduceMax);
    
    if(expsize != 0)
    {
        return; 
    }

    if(m_f->m_verbose)
    {
        if(m_f->m_comm->GetRank() == 0)
        {
            cout << "Processing input xml file" << endl;
            timer.Start();
            timerpart.Start();
        }
    }

    // check to see if fld file defined so can use in
    // expansion defintion if required
    string fldending;
    bool fldfilegiven = true;

    //Determine appropriate field input
    if(m_f->m_inputfiles.count("fld") != 0)
    {
        fldending = "fld";
    }
    else if(m_f->m_inputfiles.count("chk") != 0)
    {
        fldending = "chk";
    }
    else if (m_f->m_inputfiles.count("rst") != 0)
    {
        fldending = "rst";
    }
    else
    {
        fldfilegiven = false;
    }

    string xml_ending = "xml";
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

    SpatialDomains::DomainRangeShPtr rng =
                                SpatialDomains::NullDomainRangeShPtr;


    // define range to process output
    if(vm.count("range"))
    {
        vector<NekDouble> values;
        ASSERTL0(ParseUtils::GenerateUnOrderedVector(
                    vm["range"].as<string>().c_str(), values),
                 "Failed to interpret range string");

        ASSERTL0(values.size() > 1,
                 "Do not have minimum values of xmin,xmax");
        ASSERTL0(values.size() % 2 == 0,
                 "Do not have an even number of range values");

        int nvalues = values.size()/2;
        rng = MemoryManager<SpatialDomains::DomainRange>::
                                                AllocateSharedPtr();

        rng->m_doZrange     = false;
        rng->m_doYrange     = false;
        rng->m_checkShape   = false;

        switch(nvalues)
        {
        case 3:
            rng->m_doZrange = true;
            rng->m_zmin     = values[4];
            rng->m_zmax     = values[5];
        case 2:
            rng->m_doYrange = true;
            rng->m_ymin     = values[2];
            rng->m_ymax     = values[3];
        case 1:
            rng->m_doXrange = true;
            rng->m_xmin     = values[0];
            rng->m_xmax     = values[1];
            break;
        default:
            ASSERTL0(false,"too many values specfied in range");
        }
    }

    // define range to only take a single shape.
    if(vm.count("onlyshape"))
    {
        if(rng == SpatialDomains::NullDomainRangeShPtr)
        {
            rng = MemoryManager<SpatialDomains::DomainRange>::
                                                    AllocateSharedPtr();
            rng->m_doXrange = false;
            rng->m_doYrange = false;
            rng->m_doZrange = false;
        }

        rng->m_checkShape = true;

        string shapematch =
                    boost::to_upper_copy(vm["onlyshape"].as<string>());
        int i;
        for(i = 0; i < LibUtilities::SIZE_ShapeType; ++i)
        {
            string shapeval = LibUtilities::ShapeTypeMap[i];
            boost::to_upper(shapeval);
            if(shapematch.compare(shapeval) == 0)
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

    if(m_f->m_verbose)
    {
        cmdArgs.push_back("--verbose");
    }

    if(vm.count("shared-filesystem"))
    {
        cmdArgs.push_back("--shared-filesystem");
    }

    if(vm.count("part-only"))
    {
        cmdArgs.push_back("--part-only");
        cmdArgs.push_back(
            boost::lexical_cast<string>(vm["part-only"].as<int>()));
    }

    if(vm.count("part-only-overlapping"))
    {
        cmdArgs.push_back("--part-only-overlapping");
        cmdArgs.push_back(
            boost::lexical_cast<string>(vm["part-only-overlapping"].as<int>()));
    }

    int argc = cmdArgs.size();
    const char **argv = new const char*[argc];
    for (int i = 0; i < argc; ++i)
    {
        argv[i] = cmdArgs[i].c_str();
    }

    m_f->m_session = LibUtilities::SessionReader::
        CreateInstance(argc, (char **) argv, files, m_f->m_comm);

    // Free up memory.
    delete [] argv;

    if(m_f->m_verbose)
    {
        if(m_f->m_comm->GetRank() == 0)
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
    
    m_f->m_graph = SpatialDomains::MeshGraph::Read(m_f->m_session,rng);
    m_f->m_fld = MemoryManager<LibUtilities::FieldIO>
                    ::AllocateSharedPtr(m_f->m_session->GetComm());

    if(m_f->m_verbose)
    {
        if(m_f->m_comm->GetRank() == 0)
        {
            timerpart.Stop();
            NekDouble cpuTime = timerpart.TimePerTest(1);
            
            stringstream ss;
            ss << cpuTime << "s";
            cout << "\t InputXml mesh graph setup  CPU Time: " << setw(8) << left
                 << ss.str() << endl;
            timerpart.Start();
        }
    }

    // currently load all field (possibly could read data from
    // expansion list but it is re-arranged in expansion)
    const SpatialDomains::ExpansionMap &expansions = m_f->m_graph->GetExpansions();

    // if Range has been speficied it is possible to have a
    // partition which is empty so ccheck this and return if
    // no elements present.
    if(!expansions.size())
    {
        return;
    }

    m_f->m_exp.resize(1);

    // load fielddef header if fld file is defined. This gives
    // precedence to Homogeneous definition in fld file
    int NumHomogeneousDir = 0;
    if(fldfilegiven)
    {
        // use original expansion to identify which elements are in
        // this partition/subrange

        Array<OneD,int> ElementGIDs(expansions.size());
        SpatialDomains::ExpansionMap::const_iterator expIt;

        int i = 0; 
        for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
        {
            ElementGIDs[i++] = expIt->second->m_geomShPtr->GetGlobalID();
        }
        
        m_f->m_fld->Import(m_f->m_inputfiles[fldending][0],m_f->m_fielddef,
                           LibUtilities::NullVectorNekDoubleVector,
                           LibUtilities::NullFieldMetaDataMap,
                           ElementGIDs);
        NumHomogeneousDir = m_f->m_fielddef[0]->m_numHomogeneousDir;

        //----------------------------------------------
        // Set up Expansion information to use mode order from field
        m_f->m_graph->SetExpansions(m_f->m_fielddef);
    }
    else
    {
        if(m_f->m_session->DefinesSolverInfo("HOMOGENEOUS"))
        {
            std::string HomoStr = m_f->m_session->GetSolverInfo("HOMOGENEOUS");

            if((HomoStr == "HOMOGENEOUS1D") || (HomoStr == "Homogeneous1D")
               || (HomoStr == "1D") || (HomoStr == "Homo1D"))
            {
                NumHomogeneousDir = 1;
            }
            if((HomoStr == "HOMOGENEOUS2D") || (HomoStr == "Homogeneous2D")
               || (HomoStr == "2D") || (HomoStr == "Homo2D"))
            {
                NumHomogeneousDir = 2;
            }
        }
    }

    // reset expansion defintion to use equispaced points if required.
    if(m_requireEquiSpaced || vm.count("output-points"))
    {
        int nPointsNew = 0;

        if(vm.count("output-points"))
        {
            nPointsNew = vm["output-points"].as<int>();
        }


        m_f->m_graph->SetExpansionsToEvenlySpacedPoints(nPointsNew);
    }
    else
    {
        if(vm.count("output-points"))
        {
            int nPointsNew = vm["output-points"].as<int>();
            m_f->m_graph->SetExpansionsToPointOrder(nPointsNew);
        }
    }

    if(m_f->m_verbose)
    {
        if(m_f->m_comm->GetRank() == 0)
        {
            timerpart.Stop();
            NekDouble cpuTime = timerpart.TimePerTest(1);
            
            stringstream ss;
            ss << cpuTime << "s";
            cout << "\t InputXml setexpansion CPU Time: " << setw(8) << left
                 << ss.str() << endl;
            timerpart.Start();
        }
    }

    if(m_f->m_verbose)
    {
        if(m_f->m_comm->GetRank() == 0)
        {
            timerpart.Stop();
            NekDouble cpuTime = timerpart.TimePerTest(1);
            
            stringstream ss;
            ss << cpuTime << "s";
            cout << "\t InputXml setexpansion CPU Time: " << setw(8) << left
                 << ss.str() << endl;
            timerpart.Start();
        }
    }

    // Override number of planes with value from cmd line
    if(NumHomogeneousDir == 1 && vm.count("output-points-hom-z"))
    {
        int expdim = m_f->m_graph->GetMeshDimension();
        m_f->m_fielddef[0]->m_numModes[expdim] = vm["output-points-hom-z"].as<int>();
    }

    m_f->m_exp[0] = m_f->SetUpFirstExpList(NumHomogeneousDir,fldfilegiven);
    
    if(m_f->m_verbose)
    {
        if(m_f->m_comm->GetRank() == 0)
        {
            timerpart.Stop();
            NekDouble cpuTime = timerpart.TimePerTest(1);
            
            stringstream ss1;

            ss1 << cpuTime << "s";
            cout << "\t InputXml set first exp CPU Time: " << setw(8) << left
                 << ss1.str() << endl;

            
            timer.Stop();
            cpuTime = timer.TimePerTest(1);
            
            stringstream ss;
            ss << cpuTime << "s";
            cout << "InputXml  CPU Time: " << setw(8) << left
                 << ss.str() << endl;

        }
    }
}
}
}
