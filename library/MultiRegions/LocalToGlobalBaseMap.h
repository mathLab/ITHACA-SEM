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
#include <MultiRegions/MultiRegionsDeclspec.h>
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

        /// Base class for constructing local to global mapping of degrees of
        /// freedom.
        class LocalToGlobalBaseMap
        {
        public:
        	/// Default constructor.
            MULTI_REGIONS_EXPORT LocalToGlobalBaseMap();
            /// Constructor for next level in multi-level static condensation.
            MULTI_REGIONS_EXPORT LocalToGlobalBaseMap(LocalToGlobalBaseMap* oldLevelMap,
                    const BottomUpSubStructuredGraphSharedPtr& multiLevelGraph);
            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~LocalToGlobalBaseMap();

            /// Retrieve the global index of a given local boundary mode.
            inline int GetLocalToGlobalBndMap(const int i) const;
            /// Retrieve the global indices of the local boundary modes.
            inline const Array<OneD,const int>&  GetLocalToGlobalBndMap();
            /// Set the global indices of the local boundary modes.
            inline  void SetLocalToGlobalBndMap(Array<OneD,int> inarray);

            /// Returns true if using a modal expansion requiring a change of
            /// sign of some modes.
            inline bool GetSignChange();            
            /// Set  true if using a modal expansion requiring a change of
            /// sign of some modes.
            inline void SetSignChange(bool signChange);

            /// Retrieve the sign change of a given local boundary mode.
            inline NekDouble GetLocalToGlobalBndSign(const int i) const;
            /// Retrieve the sign change for all local boundary modes.
            inline Array<OneD, const NekDouble> GetLocalToGlobalBndSign() const;
            /// Sets the sign change for all local boundary modes.
            inline void SetLocalToGlobalBndSign(Array<OneD, NekDouble> inarray);
            /// Retrieves the global index corresponding to a boundary expansion
            /// mode.
            inline int GetBndCondCoeffsToGlobalCoeffsMap(const int i);
            /// Sets the global index corresponding to a boundary expansion
            /// mode.
            inline void SetBndCondCoeffsToGlobalCoeffsMap(Array<OneD, int> inarray);
            /// Retrieves the global indices corresponding to the boundary
            /// expansion modes.
            inline const Array<OneD,const int>&
                    GetBndCondCoeffsToGlobalCoeffsMap();
            /// Returns the modal sign associated with a given boundary
            /// expansion mode.
            inline NekDouble GetBndCondCoeffsToGlobalCoeffsSign(const int i);

            /// Sets the modal sign associated with a given boundary
            inline void SetBndCondCoeffsToGlobalCoeffsSign(Array<OneD, NekDouble> inarray);

            /// Returns the global index of the boundary trace giving the
            /// index on the boundary  expansion
            inline int GetBndCondTraceToGlobalTraceMap(const int i);
 
            /// Returns the number of global Dirichlet boundary coefficients.
            inline int GetNumGlobalDirBndCoeffs() const;
            /// Set the number of global Dirichlet boundary coefficients.
            inline void SetNumGlobalDirBndCoeffs(const int n);
            /// Returns the number of local Dirichlet boundary coefficients.
            inline int GetNumLocalDirBndCoeffs() const;
            /// Set the number of local Dirichlet boundary coefficients.
            inline void SetNumLocalDirBndCoeffs(const int n);
            /// Returns the total number of global boundary coefficients.
            inline int GetNumGlobalBndCoeffs() const;
            /// Set the total number of global boundary coefficients.
            inline void SetNumGlobalBndCoeffs(const int n);
            /// Returns the total number of local boundary coefficients.
            inline int GetNumLocalBndCoeffs() const;
            /// Sets the total number of local boundary coefficients.
            inline void SetNumLocalBndCoeffs(const int n);
            /// Returns the total number of local coefficients.
            inline int GetNumLocalCoeffs() const;
            /// Sets the total number of local coefficients.
            inline void SetNumLocalCoeffs(const int n);
            /// Returns the total number of global coefficients.
            inline int GetNumGlobalCoeffs() const;
            /// Sets the total number of global coefficients.
            inline void SetNumGlobalCoeffs(const int n);

            ///
            inline void GlobalToLocalBnd(
                    const NekVector<const NekDouble>& global,
                    NekVector<NekDouble>& loc,
                    int offset);


            inline void GlobalToLocalBnd(
                    const NekVector<const NekDouble>& global,
                    NekVector<NekDouble>& loc);

            inline void GlobalToLocalBnd(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc,
                    int offset);

            inline void GlobalToLocalBnd(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc);

            inline void AssembleBnd(const NekVector<const NekDouble>& loc,
                    NekVector<NekDouble>& global, int offset);

            inline void AssembleBnd(const NekVector<const NekDouble>& loc,
                    NekVector<NekDouble>& global);

            inline void AssembleBnd(const Array<OneD,const NekDouble>& loc,
                    Array<OneD, NekDouble>& global, int offset);

            inline void AssembleBnd(const Array<OneD, const NekDouble>& loc,
                    Array<OneD, NekDouble>& global);

            /// Returns the bandwidth of the boundary system.
            inline int GetBndSystemBandWidth() const;
            /// Sets the bandwidth of the boundary system.
            inline void SetBndSystemBandWidth(const int n);
            /// Returns the level of static condensation for this map.
            inline int GetStaticCondLevel() const;
            /// Sets the level of static condensation for this map.
            inline void SetStaticCondLevel(const int n);
            /// Returns the number of patches in this static condensation level.
            inline int GetNumPatches() const;
            /// Sets the number of patches in this static condensation level.
            inline void SetNumPatches(const int n);
            /// Returns the number of local boundary coefficients in each patch.
            inline const Array<OneD,const unsigned int>&
                    GetNumLocalBndCoeffsPerPatch();
            /// Sets the number of local boundary coefficients in each patch.
            inline void SetNumLocalBndCoeffsPerPatch(Array<OneD,unsigned int> inarray);
            /// Returns the number of local interior coefficients in each patch.
            inline const Array<OneD,const unsigned int>&
                    GetNumLocalIntCoeffsPerPatch();
            /// Sets the number of local interior coefficients in each patch.
            inline void SetNumLocalIntCoeffsPerPatch(Array<OneD,unsigned int> inarray);
            /// Returns the local to global mapping for the next level in the
            /// multi-level static condensation.
            inline const LocalToGlobalBaseMapSharedPtr
                    GetNextLevelLocalToGlobalMap() const;

            inline void SetNextLevelLocalToGlobalMap( LocalToGlobalBaseMapSharedPtr  pNextLevelLocalToGlobalMap);

            /// Returns the patch map from the previous level 
            /// of the multi-level static condensation.
            inline const PatchMapSharedPtr&
                    GetPatchMapFromPrevLevel(const int i) const;
            /// Returns true if this is the last level in the multi-level
            /// static condensation.
            inline bool AtLastLevel() const;
            /// Returns the method of solving global systems.
            inline const GlobalSysSolnType  GetGlobalSysSolnType() const;
            inline void  SetGlobalSysSolnType(GlobalSysSolnType stype);
            
        protected:
            /// Number of local boundary coefficients
            int m_numLocalBndCoeffs;
            /// Total number of global boundary coefficients
            int m_numGlobalBndCoeffs;
            /// Number of Local Dirichlet Boundary Coefficients
            int m_numLocalDirBndCoeffs;
            /// Number of Global Dirichlet Boundary Coefficients
            int m_numGlobalDirBndCoeffs;

            /// Total number of local coefficients
            /** This corresponds to the number of total number of coefficients
             *  - For CG this correpsonds to the total of bnd + int DOFs
             *  - For DG this corresponds to the number of bnd DOFs.
             *    This means that #m_numLocalCoeffs = #m_numLocalBndCoeffs
             *    This way, we can consider the trace-system solve as a
             *    statically condensed solve without interior DOFs. This allows
             *    us to use the same global system solver for both cases.
             */
            int m_numLocalCoeffs;

            /// Total number of global coefficients
            /** This corresponds to the number of total number of coefficients
             *  - For CG this correpsonds to the total of bnd + int DOFs.
             *  - For DG this corresponds to the number of bnd DOFs.
             *    This means that #m_numGlobalCoeffs = #m_numGlobalBndCoeffs
             *    This way, we can consider the trace-system solve as a
             *    statically condensed solve without interior DOFs. This allows
             *    us to use the same global system solver for both cases.
             */
            int m_numGlobalCoeffs;

            /// Flag indicating if modes require sign reversal.
            bool m_signChange;

            /// Integer map of local boundary coeffs to global space
            Array<OneD,int>       m_localToGlobalBndMap;
            /// Integer sign of local boundary coeffs to global space
            Array<OneD,NekDouble> m_localToGlobalBndSign;
            /// Integer map of bnd cond coeffs to global coefficients
            Array<OneD,int>       m_bndCondCoeffsToGlobalCoeffsMap;
            /// Integer map of bnd cond coeffs to global coefficients
            Array<OneD,NekDouble> m_bndCondCoeffsToGlobalCoeffsSign;
            /// Integer map of bnd cond trace number to global trace number
            Array<OneD,int>       m_bndCondTraceToGlobalTraceMap;

            /// The solution type of the global system
            GlobalSysSolnType m_solnType;
            /// The bandwith of the global bnd system
            int m_bndSystemBandWidth;

            /// The level of recursion in the case of multi-level static
            /// condensation.
            int m_staticCondLevel;
            /// The number of patches (~elements) in the current level
            int m_numPatches;
            /// The number of bnd dofs per patch
            Array<OneD, unsigned int> m_numLocalBndCoeffsPerPatch;
            /// The number of int dofs per patch
            Array<OneD, unsigned int> m_numLocalIntCoeffsPerPatch;
            /// Map from the patches of the previous level to the patches of
            /// the current level
            Array<OneD, PatchMapSharedPtr> m_patchMapFromPrevLevel;
            /// The local to global mapping of the next level of recursion
            LocalToGlobalBaseMapSharedPtr m_nextLevelLocalToGlobalMap;

            /// Calculates the bandwidth of the boundary system.
            void CalculateBndSystemBandWidth();

            inline void GlobalToLocalBndWithoutSign(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc);

        private:

        };


        inline int LocalToGlobalBaseMap::GetLocalToGlobalBndMap(const int i)
                                                                        const
        {
            return m_localToGlobalBndMap[i];
        }


        inline const Array<OneD,const int>&
                    LocalToGlobalBaseMap::GetLocalToGlobalBndMap(void)
        {
            return m_localToGlobalBndMap;
        }

        inline void LocalToGlobalBaseMap::SetLocalToGlobalBndMap(const Array<OneD, int> inarray)
        {
            m_localToGlobalBndMap = inarray;
        }

        inline bool LocalToGlobalBaseMap::GetSignChange()
        {
            return m_signChange;
        }

        inline void  LocalToGlobalBaseMap::SetSignChange(bool signChange)
        {
            m_signChange = signChange;
        }


        inline Array<OneD, const NekDouble>
                    LocalToGlobalBaseMap::GetLocalToGlobalBndSign(void) const
        {
            return m_localToGlobalBndSign;
        }


        inline void LocalToGlobalBaseMap::SetLocalToGlobalBndSign(Array<OneD, NekDouble> inarray)
        {
            m_localToGlobalBndSign = inarray;
        }


        inline NekDouble LocalToGlobalBaseMap::GetLocalToGlobalBndSign(
                    const int i) const
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


        inline int LocalToGlobalBaseMap::GetBndCondCoeffsToGlobalCoeffsMap(
                    const int i)
        {
            return m_bndCondCoeffsToGlobalCoeffsMap[i];
        }


        inline int LocalToGlobalBaseMap::GetBndCondTraceToGlobalTraceMap(
                    const int i)
        {
            return m_bndCondTraceToGlobalTraceMap[i];
        }


        inline NekDouble
                    LocalToGlobalBaseMap::GetBndCondCoeffsToGlobalCoeffsSign(
                    const int i)
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

        inline void LocalToGlobalBaseMap::SetBndCondCoeffsToGlobalCoeffsSign(
                                                                             Array<OneD, NekDouble> inarray )
        {
            m_bndCondCoeffsToGlobalCoeffsSign = inarray;
        }


        inline const Array<OneD,const int>&
                    LocalToGlobalBaseMap::GetBndCondCoeffsToGlobalCoeffsMap()
        {
            return m_bndCondCoeffsToGlobalCoeffsMap;
        }

        inline void LocalToGlobalBaseMap::SetBndCondCoeffsToGlobalCoeffsMap(Array<OneD,int> inarray)
        {
            m_bndCondCoeffsToGlobalCoeffsMap = inarray;
        }


        inline int LocalToGlobalBaseMap::GetNumGlobalDirBndCoeffs() const
        {
            return m_numGlobalDirBndCoeffs;
        }

        inline void LocalToGlobalBaseMap::SetNumGlobalDirBndCoeffs(const int n)
        {
            m_numGlobalDirBndCoeffs = n;
        }


        inline int LocalToGlobalBaseMap::GetNumLocalDirBndCoeffs() const
        {
            return m_numLocalDirBndCoeffs;
        }

        inline void LocalToGlobalBaseMap::SetNumLocalDirBndCoeffs(const int n)
        {
            m_numLocalDirBndCoeffs = n;
        }


        inline int LocalToGlobalBaseMap::GetNumLocalBndCoeffs() const
        {
            return m_numLocalBndCoeffs;
        }

        inline void LocalToGlobalBaseMap::SetNumLocalBndCoeffs(const int n)
        {
            m_numLocalBndCoeffs = n;
        }


        inline int LocalToGlobalBaseMap::GetNumGlobalBndCoeffs() const
        {
            return m_numGlobalBndCoeffs;
        }

        inline void LocalToGlobalBaseMap::SetNumGlobalBndCoeffs(const int n)
        {
            m_numGlobalBndCoeffs = n;
        }


        inline int LocalToGlobalBaseMap::GetNumLocalCoeffs() const
        {
            return m_numLocalCoeffs;
        }

        inline void LocalToGlobalBaseMap::SetNumLocalCoeffs(const int n)
        {
            m_numLocalCoeffs = n;
        }


        inline int LocalToGlobalBaseMap::GetNumGlobalCoeffs() const
        {
            return m_numGlobalCoeffs;
        }

        inline void LocalToGlobalBaseMap::SetNumGlobalCoeffs(const int n)
        {
            m_numGlobalCoeffs = n;
        }


        inline void LocalToGlobalBaseMap::GlobalToLocalBnd(
                    const NekVector<const NekDouble>& global,
                    NekVector<NekDouble>& loc,
                    int offset)
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


        inline void LocalToGlobalBaseMap::GlobalToLocalBnd(
                    const NekVector<const NekDouble>& global,
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


        inline void LocalToGlobalBaseMap::GlobalToLocalBnd(
                    const Array<OneD, const NekDouble>& global,
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


        inline void LocalToGlobalBaseMap::GlobalToLocalBnd(
                    const Array<OneD, const NekDouble>& global,
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


        inline void LocalToGlobalBaseMap::AssembleBnd(
                    const NekVector<const NekDouble>& loc,
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


        inline void LocalToGlobalBaseMap::AssembleBnd(
                    const NekVector<const NekDouble>& loc,
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


        inline void LocalToGlobalBaseMap::AssembleBnd(
                    const Array<OneD,const NekDouble>& loc,
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


        inline void LocalToGlobalBaseMap::AssembleBnd(
                    const Array<OneD, const NekDouble>& loc,
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


        inline int LocalToGlobalBaseMap::GetBndSystemBandWidth() const
        {
            return m_bndSystemBandWidth;
        }

        inline void LocalToGlobalBaseMap::SetBndSystemBandWidth(const int n)
        {
            m_bndSystemBandWidth = n;
        }

        inline int LocalToGlobalBaseMap::GetStaticCondLevel() const
        {
            return m_staticCondLevel;
        }
        
        inline void LocalToGlobalBaseMap::SetStaticCondLevel(const int n)
        {
            m_staticCondLevel = n;
        }


        inline int LocalToGlobalBaseMap::GetNumPatches() const
        {
            return m_numPatches;
        }

        inline void LocalToGlobalBaseMap::SetNumPatches(const int n)
        {
            m_numPatches = n;
        }


        inline const Array<OneD,const unsigned int>&
                    LocalToGlobalBaseMap::GetNumLocalBndCoeffsPerPatch()
        {
            return m_numLocalBndCoeffsPerPatch;
        }

        inline void LocalToGlobalBaseMap::SetNumLocalBndCoeffsPerPatch(Array<OneD,unsigned int> inarray)
        {
            m_numLocalBndCoeffsPerPatch = inarray;
        }

        inline const Array<OneD,const unsigned int>&
            LocalToGlobalBaseMap::GetNumLocalIntCoeffsPerPatch()
        {
            return m_numLocalIntCoeffsPerPatch;
        }

        inline void LocalToGlobalBaseMap::SetNumLocalIntCoeffsPerPatch(Array<OneD,unsigned int> inarray)
        {
            m_numLocalIntCoeffsPerPatch = inarray;
        }


        inline const LocalToGlobalBaseMapSharedPtr
                    LocalToGlobalBaseMap::GetNextLevelLocalToGlobalMap() const
        {
            return  m_nextLevelLocalToGlobalMap;
        }

        inline void LocalToGlobalBaseMap::SetNextLevelLocalToGlobalMap(
                  LocalToGlobalBaseMapSharedPtr pNextLevelLocalToGlobalMap )
        {
            m_nextLevelLocalToGlobalMap = pNextLevelLocalToGlobalMap;
        }

        inline const PatchMapSharedPtr&
                    LocalToGlobalBaseMap::GetPatchMapFromPrevLevel(const int i)
                                                                        const
        {
            return m_patchMapFromPrevLevel[i];
        }


        inline bool LocalToGlobalBaseMap::AtLastLevel() const
        {
            return !( (bool) m_nextLevelLocalToGlobalMap.get() );
        }


        inline const GlobalSysSolnType
                    LocalToGlobalBaseMap::GetGlobalSysSolnType() const
        {
            return m_solnType;
        }


        inline void LocalToGlobalBaseMap::SetGlobalSysSolnType(GlobalSysSolnType stype) 
        {
            m_solnType = stype;
        }


        inline void LocalToGlobalBaseMap::GlobalToLocalBndWithoutSign(
                    const Array<OneD, const NekDouble>& global,
                    Array<OneD,NekDouble>& loc)
        {
            ASSERTL1(loc.num_elements() >= m_numLocalBndCoeffs,"Local vector is not of correct dimension");
            ASSERTL1(global.num_elements() >= m_numGlobalBndCoeffs,"Global vector is not of correct dimension");

            Vmath::Gathr(m_numLocalBndCoeffs, global.get(), m_localToGlobalBndMap.get(), loc.get());
        }
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

