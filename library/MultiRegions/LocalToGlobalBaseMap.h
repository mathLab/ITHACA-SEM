///////////////////////////////////////////////////////////////////////////////
//
// File LocalToGlobalBaseMap.h
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
// Description: Local to Global base mapping routines, header file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef MULTIREGIONS_LOCALTOGLOBALBASEMAP_H
#define MULTIREGIONS_LOCALTOGLOBALBASEMAP_H

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/SubStructuredGraph.h>
#include <vector>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        class LocalToGlobalBaseMap;
        typedef boost::shared_ptr<LocalToGlobalBaseMap>  LocalToGlobalBaseMapSharedPtr;
        static LocalToGlobalBaseMapSharedPtr NullLocalToGlobalBaseMapSharedPtr;

        class LocalToGlobalBaseMap
        {
        public:
            LocalToGlobalBaseMap(); 

            LocalToGlobalBaseMap(LocalToGlobalBaseMap* oldLevelMap, 
                                 const BottomUpSubStructuredGraphSharedPtr& multiLevelGraph);
            
            virtual ~LocalToGlobalBaseMap(); 
                
            int GetLocalToGlobalBndMap(const int i) const
            {
                return m_localToGlobalBndMap[i];
            }
            
            inline const Array<OneD,const int>&  GetLocalToGlobalBndMap(void)
            {
                return m_localToGlobalBndMap;
            }

            bool GetSignChange()
            {
                return m_signChange;
            }

            Array<OneD, const NekDouble> GetLocalToGlobalBndSign(void) const
            {
                return m_localToGlobalBndSign;
            }

            NekDouble GetLocalToGlobalBndSign(const int i) const
            {
                if(m_signChange)
                {
                    return m_localToGlobalBndSign[i];
                }
                else
                {
                    return 1.0;
                }
            }
            
            inline int GetBndCondCoeffsToGlobalCoeffsMap(const int i)
            {
                return m_bndCondCoeffsToGlobalCoeffsMap[i];
            }

	    /**
             * \brief It returns the global index of the boundary trace giving the 
	     * index on the boundary  expansion
             */ 
	    inline int GetBndCondTraceToGlobalTraceMap(const int i)
            {
                return m_bndCondTraceToGlobalTraceMap[i];
            }
            
            inline NekDouble GetBndCondCoeffsToGlobalCoeffsSign(const int i)
            {
                if(m_signChange)
                {
                    return m_bndCondCoeffsToGlobalCoeffsSign[i];
                }
                else
                {
                    return 1.0;
                }
            }
                        
            inline const Array<OneD,const int>& GetBndCondCoeffsToGlobalCoeffsMap(void)
            {
                return m_bndCondCoeffsToGlobalCoeffsMap;
            }

            inline int GetNumGlobalDirBndCoeffs() const
            {
                return m_numGlobalDirBndCoeffs;
            }

            inline int GetNumLocalDirBndCoeffs() const
            {
                return m_numLocalDirBndCoeffs;
            }

            inline int GetNumLocalBndCoeffs() const
            {
                return m_numLocalBndCoeffs;
            }

            inline int GetNumGlobalBndCoeffs() const
            {
                return m_numGlobalBndCoeffs;
            }

            inline int GetNumLocalCoeffs() const
            {
                return m_numLocalCoeffs;
            }

            inline int GetNumGlobalCoeffs() const
            {
                return m_numGlobalCoeffs;
            }

            inline void GlobalToLocalBnd(const NekVector<const NekDouble>& global, NekVector<NekDouble>& loc, int offset)
            {
                ASSERTL1(loc.GetDimension() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
                ASSERTL1(global.GetDimension() >= m_numGlobalBndCoeffs-offset,"Global vector is not of correct dimension");
                
                // offset input data by length "offset" for Dirichlet boundary conditions.
                Array<OneD,NekDouble> tmp(global.GetDimension()+offset,0.0);
                Vmath::Vcopy(global.GetDimension(), global.GetRawPtr(), 1, tmp.get() + offset, 1);

                if(m_signChange)
                {
                    Vmath::Gathr(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(), tmp.get(), m_localToGlobalBndMap.get(), loc.GetRawPtr());
                }
                else
                {
                    Vmath::Gathr(m_numLocalBndCoeffs, tmp.get(), m_localToGlobalBndMap.get(), loc.GetRawPtr());
                }
            }     

            inline void GlobalToLocalBnd(const NekVector<const NekDouble>& global, 
                                         NekVector<NekDouble>& loc)
            {
                ASSERTL1(loc.GetDimension() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
                ASSERTL1(global.GetDimension() >= m_numGlobalBndCoeffs,"Global vector is not of correct dimension");

                if(m_signChange)
                {
                    Vmath::Gathr(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(), global.GetRawPtr(), m_localToGlobalBndMap.get(), loc.GetRawPtr());
                }
                else
                {
                    Vmath::Gathr(m_numLocalBndCoeffs, global.GetRawPtr(), m_localToGlobalBndMap.get(), loc.GetRawPtr());
                }
            } 

            inline void GlobalToLocalBnd(const Array<OneD, const NekDouble>& global, 
                                         Array<OneD,NekDouble>& loc, int offset)
            {
                ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
                ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs-offset,"Global vector is not of correct dimension");
                
                // offset input data by length "offset" for Dirichlet boundary conditions.
                Array<OneD,NekDouble> tmp(m_numGlobalBndCoeffs,0.0);
                Vmath::Vcopy(m_numGlobalBndCoeffs-offset, global.get(), 1, tmp.get() + offset, 1);

                if(m_signChange)
                {
                    Vmath::Gathr(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(), tmp.get(), m_localToGlobalBndMap.get(), loc.get());
                }
                else
                {
                    Vmath::Gathr(m_numLocalBndCoeffs, tmp.get(), m_localToGlobalBndMap.get(), loc.get());
                }
            }     

            inline void GlobalToLocalBnd(const Array<OneD, const NekDouble>& global, 
                                         Array<OneD,NekDouble>& loc)
            {
                ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
                ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs,"Global vector is not of correct dimension");
                
                if(m_signChange)
                { 
                    Vmath::Gathr(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(), global.get(), m_localToGlobalBndMap.get(), loc.get()); 
                }
                else
                {
                    Vmath::Gathr(m_numLocalBndCoeffs, global.get(), m_localToGlobalBndMap.get(), loc.get()); 
                }
            } 
            
            inline void AssembleBnd(const NekVector<const NekDouble>& loc, 
                                    NekVector<NekDouble>& global, int offset)
            {
                ASSERTL1(loc.GetDimension() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
                ASSERTL1(global.GetDimension() >= m_numGlobalBndCoeffs-offset,"Global vector is not of correct dimension");
                Array<OneD,NekDouble> tmp(global.GetDimension()+offset,0.0);

                if(m_signChange)
                {
                    Vmath::Assmb(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(), loc.GetRawPtr(), m_localToGlobalBndMap.get(), tmp.get());
                }
                else
                {                
                    Vmath::Assmb(m_numLocalBndCoeffs,loc.GetRawPtr(), m_localToGlobalBndMap.get(), tmp.get());
                }
                Vmath::Vcopy(global.GetDimension(), tmp.get() + offset, 1, global.GetRawPtr(), 1);
            }   
            
            inline void AssembleBnd(const NekVector<const NekDouble>& loc, 
                                    NekVector<NekDouble>& global)
            {
                ASSERTL1(loc.GetDimension() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
                ASSERTL1(global.GetDimension() >= m_numGlobalBndCoeffs,"Global vector is not of correct dimension");

                Vmath::Zero(m_numGlobalBndCoeffs, global.GetRawPtr(), 1);

                if(m_signChange)
                {
                    Vmath::Assmb(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(), loc.GetRawPtr(), m_localToGlobalBndMap.get(), global.GetRawPtr());
                }
                else
                {                
                    Vmath::Assmb(m_numLocalBndCoeffs,loc.GetRawPtr(), m_localToGlobalBndMap.get(), global.GetRawPtr());
                }
            }    

            inline void AssembleBnd(const Array<OneD,const NekDouble>& loc, 
                                    Array<OneD, NekDouble>& global, int offset)
            {
                ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local array is not of correct dimension");
                ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs-offset,"Global array is not of correct dimension");
                Array<OneD,NekDouble> tmp(m_numGlobalBndCoeffs,0.0);

                if(m_signChange)
                {
                    Vmath::Assmb(m_numLocalBndCoeffs, m_localToGlobalBndSign.get(),loc.get(), m_localToGlobalBndMap.get(), tmp.get());
                }
                else
                {                
                    Vmath::Assmb(m_numLocalBndCoeffs,loc.get(), m_localToGlobalBndMap.get(), tmp.get());
                }
                Vmath::Vcopy(m_numGlobalBndCoeffs-offset, tmp.get() + offset, 1, global.get(), 1);
            }   
            
            inline void AssembleBnd(const Array<OneD, const NekDouble>& loc, 
                                    Array<OneD, NekDouble>& global)
            {
                ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
                ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs,"Global vector is not of correct dimension");

                Vmath::Zero(m_numGlobalBndCoeffs, global.get(), 1);

                if(m_signChange)
                {
                    Vmath::Assmb(m_numLocalBndCoeffs,m_localToGlobalBndSign.get(), 
                                 loc.get(), m_localToGlobalBndMap.get(), global.get());
                }
                else
                {                
                    Vmath::Assmb(m_numLocalBndCoeffs,loc.get(), m_localToGlobalBndMap.get(), global.get());
                }
            }   

            inline int GetBndSystemBandWidth() const
            {
                return m_bndSystemBandWidth;
            }

            inline int GetStaticCondLevel() const
            {
                return m_staticCondLevel;
            }

            inline int GetNumPatches() const
            {
                return m_numPatches;
            }   
      
            inline const Array<OneD,const unsigned int>& GetNumLocalBndCoeffsPerPatch(void)
            {
                return m_numLocalBndCoeffsPerPatch;
            }  
      
            inline const Array<OneD,const unsigned int>& GetNumLocalIntCoeffsPerPatch(void)
            {
                return m_numLocalIntCoeffsPerPatch;
            }

            inline const LocalToGlobalBaseMapSharedPtr GetNextLevelLocalToGlobalMap() const
            {
                return  m_nextLevelLocalToGlobalMap;
            }

            inline const PatchMapSharedPtr& GetPatchMapFromPrevLevel(const int i) const
            {
                return m_patchMapFromPrevLevel[i];
            }

            inline bool AtLastLevel() const
            {
                return !( (bool) m_nextLevelLocalToGlobalMap.get() );
            }

            inline const GlobalSysSolnType  GetGlobalSysSolnType() const
            {
                return m_solnType; 
            }

        protected:
            // ---- Data members ----
            int m_numLocalBndCoeffs;     //< number of local Bnd coefficients
            int m_numGlobalBndCoeffs;    //< Total number of global boundary coefficients
            int m_numLocalDirBndCoeffs;  //< Number of Local Dirichlet Boundary Coefficient
            int m_numGlobalDirBndCoeffs; //< Number of Global Dirichlet Boundary Coefficient

            // Both data members below correspond to the number of total coefficients
            // - for CG
            //   This correpsonds to the total of bnd + int dofs
            // - for DG
            //   This corresponds to the number of bnd dofs 
            //   This means that
            //    m_numLocalCoeffs  = m_numLocalBndCoeffs
            //    m_numGlobalCoeffs = m_numGlobalBndCoeffs
            //   This way, we can consider the trace-system solve as
            //   a satically condensed solve without interior dofs
            //   This allows us to use the same global system solver for both
            //   cases.
            int m_numLocalCoeffs;      //< number of local coefficients
            int m_numGlobalCoeffs;     //< Total number of global coefficients

            bool m_signChange;
            Array<OneD,int>       m_localToGlobalBndMap;  //< integer map of local boundary coeffs to global space 
            Array<OneD,NekDouble> m_localToGlobalBndSign; //< integer sign of local boundary coeffs to global space 
 
            Array<OneD,int>       m_bndCondCoeffsToGlobalCoeffsMap;  //< integer map of bnd cond coeffs to global coefficients
            Array<OneD,NekDouble> m_bndCondCoeffsToGlobalCoeffsSign; //< integer map of bnd cond coeffs to global coefficients
	    Array<OneD,int>       m_bndCondTraceToGlobalTraceMap;  //< integer map of bnd cond trace number to global trace number

            GlobalSysSolnType m_solnType; //< The solution type of the global system
            int m_bndSystemBandWidth;     //< the bandwith of the global bnd system

            // The data below are introduced to allow a multilevel static condensation implementation
            int m_staticCondLevel;  //< The level of recursion
            int m_numPatches;       //< the number of patches (~elements) in the current level
            Array<OneD, unsigned int> m_numLocalBndCoeffsPerPatch; //< the number of bnd dofs per patch
            Array<OneD, unsigned int> m_numLocalIntCoeffsPerPatch; //< the number of int dofs per patch

            Array<OneD, PatchMapSharedPtr> m_patchMapFromPrevLevel; //< map from the patches of the previous level
                                                                    //< to the patches of the current level

            LocalToGlobalBaseMapSharedPtr m_nextLevelLocalToGlobalMap; //< The local to global mapping of 
                                                                       //< the next level of recursion

            // ---- End Data members ----

            void CalculateBndSystemBandWidth();

            inline void GlobalToLocalBndWithoutSign(const Array<OneD, const NekDouble>& global, 
                                                    Array<OneD,NekDouble>& loc)
            {
                ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
                ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs,"Global vector is not of correct dimension");
                
                Vmath::Gathr(m_numLocalBndCoeffs, global.get(), m_localToGlobalBndMap.get(), loc.get()); 
            } 
            

        private:
        }; 

        
    } // end of namespace
} // end of namespace

#endif //LOC2GLOMAP_H


/** 
 $Log: LocalToGlobalBaseMap.h,v $
 Revision 1.16  2009/10/30 14:02:55  pvos
 Multi-level static condensation updates

 Revision 1.15  2009/09/06 22:28:45  sherwin
 Updates for Navier-Stokes solver

 Revision 1.14  2009/04/27 11:42:14  sherwin
 Updated to make DG helmsolve efficient

 Revision 1.13  2009/04/08 06:38:55  sherwin
 Put eigensolve into NekMatrix. Some bandwidths mods

 Revision 1.12  2009/04/02 13:06:42  sherwin
 Modified to take symmetric banded system for HDH solver

 Revision 1.11  2009/03/04 14:17:38  pvos
 Removed all methods that take and Expansion as argument

 Revision 1.10  2009/02/08 09:10:47  sherwin
 Added NUllLocalToGlobalBaseMapSharedPtr definition

 Revision 1.9  2008/12/19 15:12:50  pvos
 Updates for precomputed dirichlet forcing functionality

 Revision 1.8  2008/12/17 17:08:22  pvos
 Performance updates

 Revision 1.7  2008/12/16 14:08:43  pvos
 Performance updates

 Revision 1.6  2008/11/01 22:36:06  bnelson
 Removed uneeded files.


 Revision 1.5  2008/11/01 22:07:46  bnelson
 Fixed compiler warning

 Revision 1.4  2008/09/23 18:21:00  pvos
 Updates for working ProjectContField3D demo

 Revision 1.3  2008/09/17 13:46:40  pvos
 Added LocalToGlobalC0ContMap for 3D expansions

 Revision 1.2  2008/09/16 13:36:06  pvos
 Restructured the LocalToGlobalMap classes

 Revision 1.1  2008/08/18 08:16:23  sherwin
 First version of this new class container for mappings

 */

