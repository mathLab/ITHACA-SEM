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

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>


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

            LibUtilities::Import(m_files[fldending][0], m_f->m_fielddef, m_f->m_data);
            
            if (m_files.count("xml") == 0 && m_files.count("xml.gz") == 0 )
            {
                return;
            }
            
            string xml_ending; 
            if(m_files.count("xml") == 0)
            {
                xml_ending = "xml.gz";
            }
            else
            {
                xml_ending = "xml";
            }

            int argc = m_files[xml_ending].size()+1;
            char *argv[argc];
            argv[0] = "ProcessField";
            for (int i = 0; i < m_files[xml_ending].size(); ++i)
            {
                argv[i+1] = strdup(m_files[xml_ending][i].c_str());
            }
            m_f->m_session = LibUtilities::SessionReader::
                CreateInstance(argc, argv);
            m_f->m_session->GetComm();
            m_f->m_graph = SpatialDomains::MeshGraph::Read(m_f->m_session);
            
            // Set up expansion list
            int expdim  = m_f->m_graph->GetMeshDimension();
            int nfields = m_f->m_fielddef[0]->m_fields.size();
            

            if(m_requireEquiSpaced) // set up points to be equispaced 
            {
                int i,j;
                int nPointsNew = -1;
                
                if(vm.count("output-points"))
                {
                    LibUtilities::Equation expession(m_f->m_session, 
                                                     vm["output-points"].as<string>());
                    nPointsNew = expession.Evaluate();
                }
                
                
                


                for(i = 0; i < m_f->m_fielddef.size(); ++i)
                {
                    vector<LibUtilities::PointsType> ptype;
                    for(j = 0; j < 3; ++j)
                    {
                        ptype.push_back(LibUtilities::ePolyEvenlySpaced);
                    }
                    
                    m_f->m_fielddef[i]->m_pointsDef = true;
                    m_f->m_fielddef[i]->m_points    = ptype;
                    
                    vector<unsigned int> porder;
                    if(nPointsNew != -1)
                    {
                        for(j = 0; j < m_f->m_fielddef[i]->m_numModes.size(); ++j)
                        {
                            porder.push_back(nPointsNew);
                        }
                    }
                    else
                    {
                        for(j = 0; j < m_f->m_fielddef[i]->m_numModes.size(); ++j)
                        {
                            porder.push_back(m_f->m_fielddef[i]->m_numModes[j]);
                        }
                    }
                    
                    m_f->m_fielddef[i]->m_numPointsDef = true;
                    m_f->m_fielddef[i]->m_numPoints = porder;
                }
                
                m_f->m_graph->SetExpansions(m_f->m_fielddef);
            }

            m_f->m_exp.resize(nfields);
        
            bool useFFT     = false;
            bool dealiasing = false;

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
                        
                        Exp2DH1 = MemoryManager<MultiRegions::
                        ExpList2DHomogeneous1D>::
                            AllocateSharedPtr(m_f->m_session, Bkey, ly, useFFT, 
                                              dealiasing, m_f->m_graph);
                        m_f->m_exp[0] = Exp2DH1;
                        
                        for (int i = 1; i < nfields; ++i)
                        {
                            m_f->m_exp[i] = MemoryManager<MultiRegions::
                            ExpList2DHomogeneous1D>::
                            AllocateSharedPtr(*Exp2DH1);
                        }
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
                        
                        Exp3DH2 = MemoryManager<MultiRegions::
                            ExpList3DHomogeneous2D>::
                            AllocateSharedPtr(m_f->m_session, BkeyY, BkeyZ, 
                                              ly, lz, useFFT, dealiasing, 
                                              m_f->m_graph);
                        m_f->m_exp[0] = Exp3DH2;
                        
                        for (int i = 1; i < nfields; ++i)
                        {
                            m_f->m_exp[i] = MemoryManager<MultiRegions::
                                ExpList3DHomogeneous2D>::
                                AllocateSharedPtr(*Exp3DH2);
                        }
                    }
                    else
                    {
                        MultiRegions::ExpList1DSharedPtr Exp1D;
                        Exp1D = MemoryManager<MultiRegions::ExpList1D>
                        ::AllocateSharedPtr(m_f->m_session, m_f->m_graph);
                        m_f->m_exp[0] = Exp1D;
                        for (int i = 1; i < nfields; ++i)
                        {
                            m_f->m_exp[i] = MemoryManager<MultiRegions::ExpList1D>
                                ::AllocateSharedPtr(*Exp1D);
                        }
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
                        
                        Exp3DH1 = MemoryManager<MultiRegions::
                            ExpList3DHomogeneous1D>::
                            AllocateSharedPtr(m_f->m_session, Bkey, lz, useFFT, 
                                              dealiasing, m_f->m_graph);
                        m_f->m_exp[0] = Exp3DH1;
                        
                        for (int i = 1; i < nfields; ++i)
                        {
                            m_f->m_exp[i] = MemoryManager<MultiRegions::
                                ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*Exp3DH1);   
                        }
                    }
                    else
                    {
                        MultiRegions::ExpList2DSharedPtr Exp2D;
                        Exp2D = MemoryManager<MultiRegions::ExpList2D>
                            ::AllocateSharedPtr(m_f->m_session,m_f->m_graph);
                        m_f->m_exp[0] = Exp2D;
                        
                        for (int i = 1; i < nfields; ++i)
                        {
                            m_f->m_exp[i] = MemoryManager<MultiRegions::ExpList2D>
                                ::AllocateSharedPtr(*Exp2D);
                        }
                    }
                }
                break;
            case 3:
                {
                    MultiRegions::ExpList3DSharedPtr Exp3D;
                    Exp3D = MemoryManager<MultiRegions::ExpList3D>
                        ::AllocateSharedPtr(m_f->m_session, m_f->m_graph);
                    m_f->m_exp[0] = Exp3D;
                    
                    for (int i = 1; i < nfields; ++i)
                    {
                        m_f->m_exp[i] = MemoryManager<MultiRegions::ExpList3D>
                            ::AllocateSharedPtr(*Exp3D);
                    }
                }
                break;
            default:
                ASSERTL0(false, "Expansion dimension not recognised");
                break;
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
        }
    }
}
