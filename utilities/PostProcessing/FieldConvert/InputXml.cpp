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
            m_f->m_session->GetComm()->Finalise();
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

            // Boundary output requires field file
            if(m_f->m_writeBndFld)
            {
                m_f->m_declareExpansionAsContField = true;
            }

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
                CreateInstance(0, 0, files, m_f->m_comm);
            m_f->m_graph = SpatialDomains::MeshGraph::Read(m_f->m_session,rng);
            m_f->m_fld = MemoryManager<LibUtilities::FieldIO>
                ::AllocateSharedPtr(m_f->m_session->GetComm());
            

            // Set up expansion list
            int expdim  = m_f->m_graph->GetMeshDimension();
                        
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


            // load fielddef if fld file is defined This gives
            // precedence to Homogeneous definition in fld file
            int NumHomogeneousDir = 0;
            if(fldfilegiven)
            {
                m_f->m_fld = MemoryManager<LibUtilities::FieldIO>
                    ::AllocateSharedPtr(m_f->m_session->GetComm());
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
                    LibUtilities::Equation expession(m_f->m_session, 
                                                     vm["output-points"].as<string>());
                    nPointsNew = expession.Evaluate();
                }
                

                m_f->m_graph->SetExpansionsToEvenlySpacedPoints(nPointsNew);
            }

            switch (expdim)
            {
            case 1:
                {
                    ASSERTL0(NumHomogeneousDir <= 2,
                             "Quasi-3D approach is only set up for 1 or 2 "
                             "homogeneous directions");
                    
                    if (NumHomogeneousDir == 1)
                    {
                        MultiRegions::ExpList2DHomogeneous1DSharedPtr Exp2DH1;
                        
                        // Define Homogeneous expansion
                        int nplanes;
                        NekDouble ly;
                        LibUtilities::BasisType btype;

                        if(fldfilegiven)
                        {
                            nplanes =  m_f->m_fielddef[0]->m_numModes[1];
                            ly      = m_f->m_fielddef[0]->m_homogeneousLengths[0];
                            btype   = m_f->m_fielddef[0]->m_basis[1];
                        }
                        else
                        {
                            m_f->m_session->LoadParameter("HomModesZ", nplanes);
                            m_f->m_session->LoadParameter("LY",ly);
                            btype = LibUtilities::eFourier;
                        }
                                                          
                        // Choose points to be at evenly spaced points at
                        // nplanes points
                        const LibUtilities::PointsKey 
                        Pkey(nplanes, LibUtilities::ePolyEvenlySpaced);

                        const LibUtilities::BasisKey Bkey(btype, nplanes, Pkey);

                        
                        
                        if(m_f->m_declareExpansionAsContField||
                            m_f->m_declareExpansionAsDisContField)
                        {
                            ASSERTL0(false,"ContField2DHomogeneous1D or "
                                     "DisContField2DHomogenenous1D has "
                                     "not been implemented");
                        }

                        Exp2DH1 = MemoryManager<MultiRegions::
                        ExpList2DHomogeneous1D>::
                            AllocateSharedPtr(m_f->m_session, Bkey, ly, useFFT, 
                                              dealiasing, m_f->m_graph);
                        m_f->m_exp[0] = Exp2DH1;
                    }
                    else if (NumHomogeneousDir == 2)
                    {
                        MultiRegions::ExpList3DHomogeneous2DSharedPtr Exp3DH2;
                        
                        int nylines,nzlines;
                        NekDouble ly,lz;
                        LibUtilities::BasisType btype1,btype2;

                        if(fldfilegiven)
                        {
                            nylines=  m_f->m_fielddef[0]->m_numModes[1];
                            nzlines=  m_f->m_fielddef[0]->m_numModes[2];
                            ly      = m_f->m_fielddef[0]->m_homogeneousLengths[0];
                            lz      = m_f->m_fielddef[0]->m_homogeneousLengths[1];
                            btype1   = m_f->m_fielddef[0]->m_basis[1];
                            btype2   = m_f->m_fielddef[0]->m_basis[2];
                        }
                        else
                        {
                            m_f->m_session->LoadParameter("HomModesY", nylines);
                            m_f->m_session->LoadParameter("HomModesZ", nzlines);
                            m_f->m_session->LoadParameter("LY",ly);
                            m_f->m_session->LoadParameter("LZ",lz);
                            btype1 = LibUtilities::eFourier;
                            btype2 = LibUtilities::eFourier;
                        }

                        // Choose points to be at evenly spaced points at
                        // nplanes points
                        const LibUtilities::PointsKey 
                        PkeyY(nylines, LibUtilities::ePolyEvenlySpaced);
                        const LibUtilities::BasisKey BkeyY(btype1, nylines, PkeyY);
                        
                        const LibUtilities::PointsKey 
                        PkeyZ(nzlines, LibUtilities::ePolyEvenlySpaced);
                        const LibUtilities::BasisKey BkeyZ(btype2, nzlines, PkeyZ);
                        
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
                    ASSERTL0(NumHomogeneousDir <= 1, 
                             "NumHomogeneousDir is only set up for 1");
                    
                    if (NumHomogeneousDir == 1)
                    {
                        MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;
                        
                        // Define Homogeneous expansion
                        int nplanes;
                        NekDouble lz;
                        LibUtilities::BasisType btype;

                        if(fldfilegiven)
                        {
                            nplanes =  m_f->m_fielddef[0]->m_numModes[2];
                            lz      = m_f->m_fielddef[0]->m_homogeneousLengths[0];
                            btype   = m_f->m_fielddef[0]->m_basis[1];
                        }
                        else
                        {
                            m_f->m_session->LoadParameter("HomModesZ", nplanes);
                            m_f->m_session->LoadParameter("LZ",lz);
                            btype = LibUtilities::eFourier;
                        }
                        
                        // Choose points to be at evenly spaced points at
                        // nplanes points
                        const LibUtilities::PointsKey 
                            Pkey(nplanes, LibUtilities::ePolyEvenlySpaced);
                        
                        const LibUtilities::BasisKey  Bkey(btype, nplanes, Pkey);
                        
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
            
        }
    }
}
