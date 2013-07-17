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

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey InputFld::className[4] = {
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "xml"), InputFld::create,
                "Reads Fld file."),
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "fld"), InputFld::create,
                "Reads Fld file."),
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "chk"), InputFld::create,
                "Reads Fld file."),
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "rst"), InputFld::create,
                "Reads Fld file.")
        };

        /**
         * @brief Set up InputFld object.
         *
         */
        InputFld::InputFld(FieldSharedPtr f) : InputModule(f)
        {
            allowedFiles.insert("xml");
            allowedFiles.insert("fld");
            allowedFiles.insert("chk");
            allowedFiles.insert("rst");
        }

        InputFld::~InputFld()
        {
            
        }

        /**
         *
         */
        void InputFld::Process()
        {
            map<string, vector<string> >::iterator it;

            string fieldfile();
            
            LibUtilities::Import(files["fld"][0], f->fielddef, f->data);
            
            if (files.count("xml") == 0)
            {
                return;
            }
            int argc = files["xml"].size()+1;
            char *argv[argc];
            argv[0] = "ProcessField";
            for (int i = 0; i < files["xml"].size(); ++i)
            {
                argv[i+1] = strdup(files["xml"][i].c_str());
            }
            f->session = LibUtilities::SessionReader::
                CreateInstance(argc, argv);
            f->graph = SpatialDomains::MeshGraph::Read(f->session);
            
            // Set up expansion list
            int expdim  = f->graph->GetMeshDimension();
            int nfields = f->fielddef[0]->m_fields.size();
            
            f->exp.resize(nfields);
            
            bool useFFT     = false;
            bool dealiasing = false;

            switch (expdim)
            {
                case 1:
                {
                    ASSERTL0(f->fielddef[0]->m_numHomogeneousDir <= 2,
                             "Quasi-3D approach is only set up for 1 or 2 "
                             "homogeneous directions");
                    
                    if (f->fielddef[0]->m_numHomogeneousDir == 1)
                    {
                        MultiRegions::ExpList2DHomogeneous1DSharedPtr Exp2DH1;
                        
                        // Define Homogeneous expansion
                        int nplanes;
                        f->session->LoadParameter("HomModesZ", nplanes, 
                                                f->fielddef[0]->m_numModes[1]);
                        
                        // Choose points to be at evenly spaced points at
                        // nplanes points
                        const LibUtilities::PointsKey 
                        Pkey(nplanes, LibUtilities::ePolyEvenlySpaced);
                        const LibUtilities::BasisKey  
                        Bkey(f->fielddef[0]->m_basis[1], nplanes, Pkey);
                        NekDouble ly = f->fielddef[0]->m_homogeneousLengths[0];
                        
                        Exp2DH1 = MemoryManager<MultiRegions::
                        ExpList2DHomogeneous1D>::
                        AllocateSharedPtr(f->session, Bkey, ly, useFFT, 
                                          dealiasing, f->graph);
                        f->exp[0] = Exp2DH1;
                        
                        for (int i = 1; i < nfields; ++i)
                        {
                            f->exp[i] = MemoryManager<MultiRegions::
                            ExpList2DHomogeneous1D>::
                            AllocateSharedPtr(*Exp2DH1);
                        }
                    }
                    else if (f->fielddef[0]->m_numHomogeneousDir == 2)
                    {
                        MultiRegions::ExpList3DHomogeneous2DSharedPtr Exp3DH2;
                        
                        int nylines;
                        int nzlines;
                        f->session->LoadParameter("HomModesY", nylines, 
                                                f->fielddef[0]->m_numModes[1]);
                        f->session->LoadParameter("HomModesZ", nzlines, 
                                                f->fielddef[0]->m_numModes[2]);
                        
                        // Choose points to be at evenly spaced points at
                        // nplanes points
                        const LibUtilities::PointsKey 
                        PkeyY(nylines, LibUtilities::ePolyEvenlySpaced);
                        const LibUtilities::BasisKey  
                        BkeyY(f->fielddef[0]->m_basis[1], nylines, PkeyY);
                        
                        const LibUtilities::PointsKey 
                        PkeyZ(nzlines, LibUtilities::ePolyEvenlySpaced);
                        const LibUtilities::BasisKey  
                        BkeyZ(f->fielddef[0]->m_basis[2], nzlines, PkeyZ);
                        
                        NekDouble ly = f->fielddef[0]->m_homogeneousLengths[0];
                        NekDouble lz = f->fielddef[0]->m_homogeneousLengths[1];
                        
                        Exp3DH2 = MemoryManager<MultiRegions::
                        ExpList3DHomogeneous2D>::
                        AllocateSharedPtr(f->session, BkeyY, BkeyZ, 
                                          ly, lz, useFFT, dealiasing, 
                                          f->graph);
                        f->exp[0] = Exp3DH2;
                        
                        for (int i = 1; i < nfields; ++i)
                        {
                            f->exp[i] = MemoryManager<MultiRegions::
                            ExpList3DHomogeneous2D>::
                            AllocateSharedPtr(*Exp3DH2);
                        }
                    }
                    else
                    {
                        MultiRegions::ExpList1DSharedPtr Exp1D;
                        Exp1D = MemoryManager<MultiRegions::ExpList1D>
                        ::AllocateSharedPtr(f->session, f->graph);
                        f->exp[0] = Exp1D;
                        for (int i = 1; i < nfields; ++i)
                        {
                            f->exp[i] = MemoryManager<MultiRegions::ExpList1D>
                            ::AllocateSharedPtr(*Exp1D);
                        }
                    }
                }
                    break;
                case 2:
                {
                    ASSERTL0(f->fielddef[0]->m_numHomogeneousDir <= 1, 
                             "NumHomogeneousDir is only set up for 1");
                    
                    if (f->fielddef[0]->m_numHomogeneousDir == 1)
                    {
                        MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;
                        
                        int nplanes;
                        f->session->LoadParameter("HomModesZ", nplanes, 
                                                f->fielddef[0]->m_numModes[2]);
                        
                        // Choose points to be at evenly spaced points at
                        // nplanes points
                        const LibUtilities::PointsKey 
                        Pkey(nplanes, LibUtilities::ePolyEvenlySpaced);
                        const LibUtilities::BasisKey  
                        Bkey(f->fielddef[0]->m_basis[2], nplanes, Pkey);
                        NekDouble lz = f->fielddef[0]->m_homogeneousLengths[0];
                        
                        Exp3DH1 = MemoryManager<MultiRegions::
                        ExpList3DHomogeneous1D>::
                        AllocateSharedPtr(f->session, Bkey, lz, useFFT, 
                                          dealiasing, f->graph);
                        f->exp[0] = Exp3DH1;
                        
                        for (int i = 1; i < nfields; ++i)
                        {
                            f->exp[i] = MemoryManager<MultiRegions::
                            ExpList3DHomogeneous1D>::
                            AllocateSharedPtr(*Exp3DH1);   
                        }
                    }
                    else
                    {
                        MultiRegions::ExpList2DSharedPtr Exp2D;
                        Exp2D = MemoryManager<MultiRegions::ExpList2D>
                        ::AllocateSharedPtr(f->session,f->graph);
                        f->exp[0] = Exp2D;
                        
                        for (int i = 1; i < nfields; ++i)
                        {
                            f->exp[i] = MemoryManager<MultiRegions::ExpList2D>
                            ::AllocateSharedPtr(*Exp2D);
                        }
                    }
                }
                    break;
                case 3:
                {
                    MultiRegions::ExpList3DSharedPtr Exp3D;
                    Exp3D = MemoryManager<MultiRegions::ExpList3D>
                    ::AllocateSharedPtr(f->session, f->graph);
                    f->exp[0] = Exp3D;
                    
                    for (int i = 1; i < nfields; ++i)
                    {
                        f->exp[i] = MemoryManager<MultiRegions::ExpList3D>
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
                for (int i = 0; i < f->data.size(); ++i)
                {
                    f->exp[j]->ExtractDataToCoeffs(f->fielddef[i], f->data[i],
                                                   f->fielddef[i]->m_fields[j],
                                                   f->exp[j]->UpdateCoeffs());
                }
                f->exp[j]->BwdTrans(f->exp[j]->GetCoeffs(), 
                                    f->exp[j]->UpdatePhys());
            }
        }
    }
}