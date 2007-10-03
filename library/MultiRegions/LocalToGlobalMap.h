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
#include <vector>
#include <LibUtilities/Memory/NekMemoryManager.hpp>


namespace Nektar
{
    namespace MultiRegions
    {
    
    class LocalToGlobalMap
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
        
        inline void LocalToCont(const ConstArray<OneD, NekDouble> &loc, 
                                Array<OneD, NekDouble> &cont)
        {
            Vmath::Scatr(m_totLocDofs, &loc[0],&m_locToContMap[0],&cont[0]);
        }                

        
        inline void ContToLocal(const ConstArray<OneD, NekDouble> &cont, 
                                Array<OneD, NekDouble> &loc)
        {
            Vmath::Gathr(m_totLocDofs,&cont[0],&m_locToContMap[0], &loc[0]);
        }
        
        inline void Assemble(const ConstArray<OneD, NekDouble> &loc, 
                             Array<OneD, NekDouble> &cont)
        {
        Vmath::Zero(m_totGloDofs,&cont[0],1);
        
                Vmath::Assmb(m_totLocDofs,&loc[0],&m_locToContMap[0],&cont[0]);
        }


        inline void ContToLocalBnd(const DNekVec &cont, DNekVec &loc, int offset = 0)
        {
            ASSERTL1(loc.GetDimension() >= m_totLocBndDofs,"Local vector is not of correct dimension");
            ASSERTL1(cont.GetDimension() >= m_totGloBndDofs-offset,"Global vector is not of correct dimension");

            // offset input data by length "offset" for Dirichlet boundary conditions.
            Array<OneD,NekDouble> tmp(cont.GetDimension()+offset,0.0);
            Vmath::Vcopy(cont.GetDimension(),cont.GetPtr(),1,&tmp[offset],1);
            
            Vmath::Gathr(m_totLocBndDofs,&tmp[0],&m_locToContBndMap[0], loc.GetPtr());
        }


        inline void AssembleBnd(const DNekVec &loc, DNekVec &cont, int offset = 0)
        {
            ASSERTL1(loc.GetDimension() >= m_totLocBndDofs,"Local vector is not of correct dimension");
            ASSERTL1(cont.GetDimension() >= m_totGloBndDofs-offset,"Global vector is not of correct dimension");
            Array<OneD,NekDouble> tmp(cont.GetDimension()+offset,0.0);

            Vmath::Assmb(m_totLocBndDofs,loc.GetPtr(), &m_locToContBndMap[0], &tmp[0]);
            Vmath::Vcopy(cont.GetDimension(),&tmp[offset],1,cont.GetPtr(),1);
        }
        
        inline int GetTotGloDofs()
        {
            return m_totGloDofs;
        }
        
        inline int GetTotGloBndDofs()
        {
            return m_totGloBndDofs;
        }
        
        inline void SetNumDirichetBCs(const int NumDirichlet)
        {
            m_numDirichletBCs = NumDirichlet;
        }

        inline void SetNumNeumannBCs(const int NumNeumann)
        {
            m_numNeumannBCs = NumNeumann;
        }

        inline void SetNumRobinBCs(const int NumNeumann)
        {
            m_numNeumannBCs = NumNeumann;
        }
        
        inline int GetNumDirichletBCs()
        {
            return m_numDirichletBCs;
        }

        inline int GetNumNeumannBCs()
        {
            return m_numNeumannBCs;
        }

        inline int GetNumRobinBCs()
        {
            return m_numRobinBCs;
        }
        
        void ResetMapping(const int NumDirichlet,
                          SpatialDomains::BoundaryConditions &bcs,
                          const std::string variable)
        {
            v_ResetMapping(NumDirichlet,bcs,variable);
        }
        
        inline const ConstArray<OneD, int> &GetBndCondGlobalID() const 
        {
            return m_bndCondGlobalID;
        }

        protected:
            int             m_totLocDofs;      //< number of local dofs
            int             m_totLocBndDofs;   //< number of local dofs
            int             m_totGloDofs;      //< number of global dofs 
            int             m_totGloBndDofs;   //< number of global boundary dofs;
            int             m_numDirichletBCs; //< number of Dirichlet conditions 
            int             m_numNeumannBCs;   //< number of Neumann conditions 
            int             m_numRobinBCs;     //< number of Robin conditions 
            Array<OneD,int> m_locToContMap;    //< integer map
            Array<OneD,int> m_locToContBndMap; //< integer maps of boundary dofs
            Array<OneD,int> m_bndCondGlobalID; //< global id of all boundary conditions

        private:
            virtual void v_ResetMapping(const int NumDirichlet, 
                                        SpatialDomains::BoundaryConditions &bcs,
                                        const std::string variable)
            {
                ASSERTL0(false,"ResetMapping needs defining");
            }
    };
    
    typedef boost::shared_ptr<LocalToGlobalMap>  LocalToGlobalMapSharedPtr;

    } // end of namespace
} // end of namespace

#endif //LOC2GLOMAP_H


/** $Log: LocalToGlobalMap.h,v $
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

