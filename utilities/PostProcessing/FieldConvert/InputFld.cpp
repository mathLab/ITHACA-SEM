////////////////////////////////////////////////////////////////////////////////
//
//  File: InputFld.cpp
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
//  Description: GMSH converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "InputFld.h"


static std::string npts = LibUtilities::SessionReader::RegisterCmdLineArgument(
                "NumberOfPoints","n","Define number of points to dump output");

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey InputFld::m_className[5] = {
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "xml"), InputFld::create,
                "Reads Xml file."),
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "xml.gz"), InputFld::create,
                "Reads Xml file."),
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "fld"), InputFld::create,
                "Reads Fld file."),
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "chk"), InputFld::create,
                "Reads Fld file."),
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "rst"), InputFld::create,
                "Reads Fld file."),
        };

        /**
         * @brief Set up InputFld object.
         *
         */
        InputFld::InputFld(FieldSharedPtr f) : InputModule(f)
        {
            m_allowedFiles.insert("xml");
            m_allowedFiles.insert("xml.gz");
            m_allowedFiles.insert("fld");
            m_allowedFiles.insert("chk");
            m_allowedFiles.insert("rst");
        }

        InputFld::~InputFld()
        {
            m_f->m_session->GetComm()->Finalise();
        }

        /**
         *
         */
        void InputFld::Process(po::variables_map &vm)
        {
            string fldending;
            //Determine appropriate field input 
            if(m_files.count("fld") != 0)
            {
                fldending = "fld";
            }
            else if(m_files.count("chk") != 0)
            {
                fldending = "chk";
            }
            else if (m_files.count("rst") != 0)
            {
                fldending = "rst";
            }
            else
            {
                ASSERTL0(false,"no input file found");
            }
            
            if (m_files.count("xml") == 0 && m_files.count("xml.gz") == 0 )
            {
                return;
            }
            
            string xml_ending = "xml";
            string xml_gz_ending = "xml.gz";

            if(vm.count("boundary-region"))
            {
                m_f->m_declareExpansionAsContField = true;
            }

            int   argc = m_files[xml_ending].size()+m_files[xml_gz_ending].size()+1;
            char *argv[argc];
            const char *instring = "ProcessField";
            argv[0] = strdup(instring);
            // load .xml ending
            int i;
            for (i = 0; i < m_files[xml_ending].size(); ++i)
            {
                argv[i+1] = strdup(m_files[xml_ending][i].c_str());
            }

            // load any .xml.gz endings
            for (int j =0; j < m_files[xml_gz_ending].size(); ++j)
            {
                argv[i+j+1] = strdup(m_files[xml_gz_ending][j].c_str());
            }


            SpatialDomains::DomainRangeShPtr rng = SpatialDomains::NullDomainRangeShPtr; 
            // define range to process output 
            if(vm.count("range"))
            {
                vector<NekDouble> values;
                ASSERTL0(ParseUtils::GenerateUnOrderedVector(vm["range"].as<string>().c_str(),values),"Failed to interpret range string");
                
                ASSERTL0(values.size() > 1,"Do not have minimum values of xmin,xmax");
                ASSERTL0(values.size()%2 == 0,"Do not have an even number of range values");
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


            m_f->m_session = LibUtilities::SessionReader::
                CreateInstance(argc, argv);
            m_f->m_session->GetComm();
            m_f->m_graph = SpatialDomains::MeshGraph::Read(m_f->m_session,rng);
            m_f->m_fld = MemoryManager<LibUtilities::FieldIO>
                ::AllocateSharedPtr(m_f->m_session->GetComm());
            

            // Set up expansion list
            int expdim  = m_f->m_graph->GetMeshDimension();
            
            if(m_requireEquiSpaced) // set up points to be equispaced 
            {
                int nPointsNew = 0;
                
                if(vm.count("output-points"))
                {
                    LibUtilities::Equation expession(m_f->m_session, 
                                                     vm["output-points"].as<string>());
                    nPointsNew = expession.Evaluate();
                }
                

                m_f->m_graph->SetExpansionsToEvenlySpacedPoints(nPointsNew);
            }
            
            bool useFFT     = false;
            bool dealiasing = false;
            
            // currently load all field (possibly could read data from expansion list
            // but it is re-arranged in expansion) 
            
            const SpatialDomains::ExpansionMap &expansions = m_f->m_graph->GetExpansions();
  
            
            // if Range has been speficied it is possible to have a
            // partition which is empty so ccheck this and return if
            // no elements present.
            if(!expansions.size())
            {
                return;
            }

            m_f->m_exp.resize(1);

            Array<OneD,int> ElementGIDs(expansions.size());
            SpatialDomains::ExpansionMap::const_iterator expIt;

            i = 0;
            for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
            {
                ElementGIDs[i++] = expIt->second->m_geomShPtr->GetGlobalID();
            }

            m_f->m_fld->Import(m_files[fldending][0],m_f->m_fielddef,m_f->m_data,
                               LibUtilities::NullFieldMetaDataMap,
                               ElementGIDs);

            switch (expdim)
            {
            case 1:
                {
                    ASSERTL0(m_f->m_fielddef[0]->m_numHomogeneousDir <= 2,
                             "Quasi-3D approach is only set up for 1 or 2 "
                             "homogeneous directions");
                    
                    if (m_f->m_fielddef[0]->m_numHomogeneousDir == 1)
                    {
                        MultiRegions::ExpList2DHomogeneous1DSharedPtr Exp2DH1;
                        
                        // Define Homogeneous expansion
                        int nplanes;
                        m_f->m_session->LoadParameter("HomModesZ", nplanes, 
                                             m_f->m_fielddef[0]->m_numModes[1]);
                        
                        // Choose points to be at evenly spaced points at
                        // nplanes points
                        const LibUtilities::PointsKey 
                        Pkey(nplanes, LibUtilities::ePolyEvenlySpaced);
                        const LibUtilities::BasisKey  
                            Bkey(m_f->m_fielddef[0]->m_basis[1], nplanes, Pkey);
                        NekDouble ly = m_f->m_fielddef[0]->m_homogeneousLengths[0];
                        
                        
                        if(m_f->m_declareExpansionAsContField||
                            m_f->m_declareExpansionAsDisContField)
                        {
                            ASSERTL0(false,"ContField2DHomogeneous1D or DisContField2DHomogenenous1D has not been implemented");
                        }

                        Exp2DH1 = MemoryManager<MultiRegions::
                        ExpList2DHomogeneous1D>::
                            AllocateSharedPtr(m_f->m_session, Bkey, ly, useFFT, 
                                              dealiasing, m_f->m_graph);
                        m_f->m_exp[0] = Exp2DH1;
                    }
                    else if (m_f->m_fielddef[0]->m_numHomogeneousDir == 2)
                    {
                        MultiRegions::ExpList3DHomogeneous2DSharedPtr Exp3DH2;
                        
                        int nylines;
                        int nzlines;
                        m_f->m_session->LoadParameter("HomModesY", nylines, 
                                            m_f->m_fielddef[0]->m_numModes[1]);
                        m_f->m_session->LoadParameter("HomModesZ", nzlines, 
                                            m_f->m_fielddef[0]->m_numModes[2]);
                        
                        // Choose points to be at evenly spaced points at
                        // nplanes points
                        const LibUtilities::PointsKey 
                        PkeyY(nylines, LibUtilities::ePolyEvenlySpaced);
                        const LibUtilities::BasisKey  
                            BkeyY(m_f->m_fielddef[0]->m_basis[1], nylines, PkeyY);
                        
                        const LibUtilities::PointsKey 
                        PkeyZ(nzlines, LibUtilities::ePolyEvenlySpaced);
                        const LibUtilities::BasisKey  
                        BkeyZ(m_f->m_fielddef[0]->m_basis[2], nzlines, PkeyZ);
                        
                        NekDouble ly = m_f->m_fielddef[0]->m_homogeneousLengths[0];
                        NekDouble lz = m_f->m_fielddef[0]->m_homogeneousLengths[1];

                        if(m_f->m_declareExpansionAsContField)
                        {
                            Exp3DH2 = MemoryManager<MultiRegions::
                                ContField3DHomogeneous2D>::
                                AllocateSharedPtr(m_f->m_session, BkeyY, BkeyZ, 
                                                  ly, lz, useFFT, dealiasing, 
                                                  m_f->m_graph,
                                                  m_f->m_session->GetVariable(0));
                        }
                        else if(m_f->m_declareExpansionAsContField)
                        {
                            Exp3DH2 = MemoryManager<MultiRegions::
                                DisContField3DHomogeneous2D>::
                                AllocateSharedPtr(m_f->m_session, BkeyY, BkeyZ, 
                                                  ly, lz, useFFT, dealiasing, 
                                                  m_f->m_graph,
                                                  m_f->m_session->GetVariable(0));
                        }
                        else
                        {
                            Exp3DH2 = MemoryManager<MultiRegions::
                                ExpList3DHomogeneous2D>::
                                AllocateSharedPtr(m_f->m_session, BkeyY, BkeyZ, 
                                                  ly, lz, useFFT, dealiasing, 
                                                  m_f->m_graph);
                        }


                        m_f->m_exp[0] = Exp3DH2;
                    }
                    else
                    {
                        MultiRegions::ExpList1DSharedPtr Exp1D;

                        if(m_f->m_declareExpansionAsContField)
                        {
                            Exp1D = MemoryManager<MultiRegions::ContField1D>
                                ::AllocateSharedPtr(m_f->m_session, m_f->m_graph,
                                                    m_f->m_session->GetVariable(0));
                        }
                        else if(m_f->m_declareExpansionAsContField)
                        {
                            Exp1D = MemoryManager<MultiRegions::DisContField1D>
                                ::AllocateSharedPtr(m_f->m_session, m_f->m_graph,
                                                  m_f->m_session->GetVariable(0));
                        }
                        else 
                        {
                            Exp1D = MemoryManager<MultiRegions::ExpList1D>
                                ::AllocateSharedPtr(m_f->m_session, m_f->m_graph);
                        }


                        m_f->m_exp[0] = Exp1D;
                    }
                }
                break;
            case 2:
                {
                    ASSERTL0(m_f->m_fielddef[0]->m_numHomogeneousDir <= 1, 
                             "NumHomogeneousDir is only set up for 1");
                    
                    if (m_f->m_fielddef[0]->m_numHomogeneousDir == 1)
                    {
                        MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;
                        
                        int nplanes;
                        m_f->m_session->LoadParameter("HomModesZ", nplanes, 
                                           m_f->m_fielddef[0]->m_numModes[2]);
                        
                        // Choose points to be at evenly spaced points at
                        // nplanes points
                        const LibUtilities::PointsKey 
                            Pkey(nplanes, LibUtilities::ePolyEvenlySpaced);
                        const LibUtilities::BasisKey  
                            Bkey(m_f->m_fielddef[0]->m_basis[2], nplanes, Pkey);
                        NekDouble lz = m_f->m_fielddef[0]->m_homogeneousLengths[0];
                        
                        if(m_f->m_declareExpansionAsContField)
                        {
                            Exp3DH1 = MemoryManager<MultiRegions::
                                ContField3DHomogeneous1D>::
                                AllocateSharedPtr(m_f->m_session, Bkey, lz, useFFT, 
                                                  dealiasing, m_f->m_graph,
                                                  m_f->m_session->GetVariable(0));
                        }
                        else if (m_f->m_declareExpansionAsContField)
                        {
                            Exp3DH1 = MemoryManager<MultiRegions::
                                DisContField3DHomogeneous1D>::
                                AllocateSharedPtr(m_f->m_session,
                                                  Bkey, lz, useFFT, 
                                                  dealiasing, m_f->m_graph,
                                                  m_f->m_session->GetVariable(0));
                        }
                        else
                        {
                            Exp3DH1 = MemoryManager<MultiRegions::
                                ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(m_f->m_session, Bkey, lz, useFFT, 
                                                  dealiasing, m_f->m_graph);
                        }
                        m_f->m_exp[0] = Exp3DH1;
                    }
                    else
                    {
                        MultiRegions::ExpList2DSharedPtr Exp2D;

                        if(m_f->m_declareExpansionAsContField)
                        {
                            Exp2D = MemoryManager<MultiRegions::ContField2D>
                                ::AllocateSharedPtr(m_f->m_session,m_f->m_graph,
                                                    m_f->m_session->GetVariable(0));
                        }
                        else if(m_f->m_declareExpansionAsDisContField)
                        {
                            Exp2D = MemoryManager<MultiRegions::DisContField2D>
                                ::AllocateSharedPtr(m_f->m_session,m_f->m_graph,
                                                    m_f->m_session->GetVariable(0));
                        }
                        else
                        {
                            Exp2D = MemoryManager<MultiRegions::ExpList2D>
                                ::AllocateSharedPtr(m_f->m_session,m_f->m_graph);
                        }
                        m_f->m_exp[0] = Exp2D;
                        
                    }
                }
                break;
            case 3:
                {
                    MultiRegions::ExpList3DSharedPtr Exp3D;

                    if(m_f->m_declareExpansionAsContField)
                    {
                        Exp3D = MemoryManager<MultiRegions::ContField3D>
                            ::AllocateSharedPtr(m_f->m_session,m_f->m_graph,
                                                m_f->m_session->GetVariable(0));
                    }
                    else if(m_f->m_declareExpansionAsDisContField)
                    {
                        Exp3D = MemoryManager<MultiRegions::DisContField3D>
                            ::AllocateSharedPtr(m_f->m_session,m_f->m_graph,
                                                m_f->m_session->GetVariable(0));
                    }
                    else
                    {
                        Exp3D = MemoryManager<MultiRegions::ExpList3D>
                            ::AllocateSharedPtr(m_f->m_session, m_f->m_graph);
                    }
                    m_f->m_exp[0] = Exp3D;
                }
                break;
            default:
                ASSERTL0(false, "Expansion dimension not recognised");
                break;
            }
            
            
            int nfields = m_f->m_fielddef[0]->m_fields.size();
            m_f->m_exp.resize(nfields);
        
            vector<string> vars = m_f->m_session->GetVariables();

            // declare other fields; 
            for (int i = 1; i < nfields; ++i)
            {
                if(i < vars.size())
                {
                    m_f->m_exp[i] = m_f->AppendExpList(vars[i]);
                }
                else
                {
                    m_f->m_exp[i] = m_f->AppendExpList();
                }                   
            }
            
            for (int j = 0; j < nfields; ++j)
            {
                for (int i = 0; i < m_f->m_data.size(); ++i)
                {
                    m_f->m_exp[j]->ExtractDataToCoeffs(m_f->m_fielddef[i], 
                                                       m_f->m_data[i],
                                                       m_f->m_fielddef[i]->m_fields[j],
                                                       m_f->m_exp[j]->UpdateCoeffs());
                }
                m_f->m_exp[j]->BwdTrans(m_f->m_exp[j]->GetCoeffs(), 
                                        m_f->m_exp[j]->UpdatePhys());
            }


            // if range is defined reset up output field in case or
            // reducing fld definition
            if(vm.count("range"))
            {
                std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                    = m_f->m_exp[0]->GetFieldDefinitions();
                std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
                
                for (int j = 0; j < nfields; ++j)
                {
                    for (i = 0; i < FieldDef.size(); ++i)
                    {   
                        FieldDef[i]->m_fields.push_back(m_f->m_fielddef[0]->m_fields[j]);
                        m_f->m_exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
                    }
                }   
                m_f->m_fielddef = FieldDef;
                m_f->m_data     = FieldData;
            }
        }
    }
}
