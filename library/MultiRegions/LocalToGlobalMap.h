///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalMap.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Local to Global mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_LOC2GLOMAP_H
#define NEKTAR_LIB_MULTIREGIONS_LOC2GLOMAP_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/LocalToGlobalBndryMap.h>
#include <vector>
#include <LibUtilities/Memory/NekMemoryManager.hpp>


namespace Nektar
{
    namespace MultiRegions
    {
        
    class LocalToGlobalMap:
        public LocalToGlobalBndryMap
        {
        public:
            LocalToGlobalMap();
            LocalToGlobalMap(const int totdata, const int totbnddata, 
                             Array<OneD,int> &map);
            ~LocalToGlobalMap();
            
            inline int GetMap(const int i) const
            {
                return m_locToContMap[i];
            }
            
            inline NekDouble GetSign(int i) 
            {
                return v_GetSign(i);
            }
            
            inline void LocalToCont(const ConstArray<OneD, NekDouble> &loc, 
                                    Array<OneD, NekDouble> &cont)
            {
                v_LocalToCont(loc,cont);
            }                            
            
            inline void ContToLocal(const ConstArray<OneD, NekDouble> &cont, 
                                    Array<OneD, NekDouble> &loc)
            {
                v_ContToLocal(cont,loc);
            }
            
            inline void Assemble(const ConstArray<OneD, NekDouble> &loc, 
                                 Array<OneD, NekDouble> &cont)
            {
                v_Assemble(loc,cont);
            }
            
            inline int GetTotGloDofs()
            {
                return m_totGloDofs;
            }

            inline int GetTotLocDofs()
            {
                return m_totLocDofs;
            }
            
        protected:
            int             m_totLocDofs;         //< number of local dofs
            int             m_totGloDofs;         //< number of global dofs 
            Array<OneD,int> m_locToContMap;       //< integer map
            
        private:
            virtual void v_LocalToCont(const ConstArray<OneD, NekDouble> &loc, 
                                       Array<OneD, NekDouble> &cont)
            {
                ASSERTL0(false,"LocalToCont needs defining");
            }                
            
            
            virtual void v_ContToLocal(const ConstArray<OneD, NekDouble> &cont, 
                                       Array<OneD, NekDouble> &loc)
            {
                ASSERTL0(false,"ContToLocal needs defining");
            }
            
            virtual void v_Assemble(const ConstArray<OneD, NekDouble> &loc, 
                                    Array<OneD, NekDouble> &cont)
            {
                ASSERTL0(false,"Assemble needs defining");
            }

            virtual NekDouble v_GetSign(int i) 
            {
                ASSERTL0(false,"GetSign needs defining");
                return 0.0;
            }

        };
        
        typedef boost::shared_ptr<LocalToGlobalMap>  LocalToGlobalMapSharedPtr;
        
    } // end of namespace
} // end of namespace

#endif //LOC2GLOMAP_H


/** $Log: LocalToGlobalMap.h,v $
/** Revision 1.15  2007/12/06 22:52:30  pvos
/** 2D Helmholtz solver updates
/**
/** Revision 1.14  2007/11/20 16:27:16  sherwin
/** Zero Dirichlet version of UDG Helmholtz solver
/**
/** Revision 1.13  2007/10/04 13:57:01  pvos
/** fixed some more errors
/**
/** Revision 1.12  2007/10/03 11:37:50  sherwin
/** Updates relating to static condensation implementation
/**
/** Revision 1.11  2007/09/25 14:25:29  pvos
/** Update for helmholtz1D with different expansion orders
/**
/** Revision 1.10  2007/07/22 23:04:21  bnelson
/** Backed out Nektar::ptr.
/**
/** Revision 1.9  2007/07/20 02:04:13  bnelson
/** Replaced boost::shared_ptr with Nektar::ptr
/**
/** Revision 1.8  2007/07/16 18:28:43  sherwin
/** Modification to introduce non-zero Dirichlet boundary conditions into the Helmholtz1D Demo
/**
/** Revision 1.7  2007/06/08 12:58:27  sherwin
/** Added ContField1D and remove previous structure using Fields
/**
/** Revision 1.6  2007/05/28 16:15:00  sherwin
/** Updated files in MultiRegions to make 1D demos work
/**
/** Revision 1.5  2007/05/27 16:09:43  bnelson
/** Update to new Array type.
/**
/** Revision 1.4  2007/04/26 15:00:16  sherwin
/** SJS compiling working version using SHaredArrays
/**
/** Revision 1.3  2007/03/20 16:58:42  sherwin
/** Update to use Array<OneD, NekDouble> storage and NekDouble usage, compiling and executing up to Demos/StdRegions/Project1D
/**
/** Revision 1.2  2007/03/14 21:24:08  sherwin
/** Update for working version of MultiRegions up to ExpList1D
/**
/** Revision 1.1  2006/07/02 17:16:17  sherwin
/**
/** Modifications to make MultiRegions work for a connected domain in 2D (Tris)
/**
/** Revision 1.3  2006/06/05 00:14:33  bnelson
/** Fixed a compiler error (couldn't find boost::shared_ptr<) and a couple of formatting updates for the standard.
/** */

