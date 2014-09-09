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
using namespace std;

#include "InputXml.h"

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

        InputXml::~InputXml()
        {
        }

        /**
         *
         */
        void InputXml::Process(po::variables_map &vm)
        {

            if(m_f->m_verbose)
            {
                cout << "Processing input xml file" << endl;
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
                rng = MemoryManager<SpatialDomains::DomainRange>::AllocateSharedPtr();

                rng->doZrange = false;
                rng->doYrange = false;

                switch(nvalues)
                {
                case 3:
                    rng->doZrange = true;
                    rng->zmin = values[4];
                    rng->zmax = values[5];
                case 2:
                    rng->doYrange = true;
                    rng->ymin = values[2];
                    rng->ymax = values[3];
                case 1:
                    rng->doXrange = true;
                    rng->xmin = values[0];
                    rng->xmax = values[1];
                    break;
                default:
                    ASSERTL0(false,"too many values specfied in range");
                }    
            }

            if(m_f->m_verbose)
            {
                string firstarg = "FieldConvert";
                string verbose = "-v";
                char **argv; 
                argv = (char**)malloc(2*sizeof(char*));
                argv[0] = (char *)malloc(firstarg.size()*sizeof(char));
                argv[1] = (char *)malloc(verbose.size()*sizeof(char));

                sprintf(argv[0],"%s",firstarg.c_str());
                sprintf(argv[1],"%s",verbose.c_str());

                m_f->m_session = LibUtilities::SessionReader::
                    CreateInstance(2, (char **)argv, files, m_f->m_comm);
            }
            else
            {
                m_f->m_session = LibUtilities::SessionReader::
                    CreateInstance(0, 0, files, m_f->m_comm);
            }




            m_f->m_graph = SpatialDomains::MeshGraph::Read(m_f->m_session,rng);
            m_f->m_fld = MemoryManager<LibUtilities::FieldIO>
                            ::AllocateSharedPtr(m_f->m_session->GetComm());

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

            // load fielddef if fld file is defined This gives
            // precedence to Homogeneous definition in fld file
            int NumHomogeneousDir = 0;
            if(fldfilegiven)
            {
                m_f->m_fld->Import(m_f->m_inputfiles[fldending][0],m_f->m_fielddef);
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
            if(m_requireEquiSpaced) // set up points to be equispaced 
            {
                int nPointsNew = 0;
                
                if(vm.count("output-points"))
                {
                    nPointsNew = vm["output-points"].as<int>();
                }
                

                m_f->m_graph->SetExpansionsToEvenlySpacedPoints(nPointsNew);
            }

            m_f->m_exp[0] = m_f->SetUpFirstExpList(NumHomogeneousDir,fldfilegiven);
        }
    }
}
