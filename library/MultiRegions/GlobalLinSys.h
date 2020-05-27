///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSys.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: GlobalLinSys header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYS_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYS_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <MultiRegions/GlobalLinSysKey.h>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;
        class GlobalLinSys;
        class Preconditioner;

        /// Pointer to a GlobalLinSys object.
        typedef std::shared_ptr<GlobalLinSys> GlobalLinSysSharedPtr;
        /// Mapping between GlobalLinSys objects and their associated keys.
        typedef std::map<GlobalLinSysKey,GlobalLinSysSharedPtr> GlobalLinSysMap;
        /// Pointer to a GlobalLinSys/key map.
        typedef std::shared_ptr<GlobalLinSysMap> GlobalLinSysMapShPtr;

        // Forward declaration
        typedef std::shared_ptr<Preconditioner> PreconditionerSharedPtr;

        /// Datatype of the NekFactory used to instantiate classes derived from
        /// the EquationSystem class.
        typedef LibUtilities::NekFactory< std::string, GlobalLinSys, 
            const GlobalLinSysKey&,
            const std::weak_ptr<ExpList>&,
            const std::shared_ptr<AssemblyMap>& > GlobalLinSysFactory;
        GlobalLinSysFactory& GetGlobalLinSysFactory();


        /// A global linear system.
        class GlobalLinSys: public std::enable_shared_from_this<GlobalLinSys>
        {
        public:
            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSys(
                const GlobalLinSysKey                &pKey,
                const std::weak_ptr<ExpList>       &pExpList,
                const std::shared_ptr<AssemblyMap> &pLocToGloMap);

            MULTI_REGIONS_EXPORT
            virtual ~GlobalLinSys() {}

            /// Returns the key associated with the system.
            const inline GlobalLinSysKey &GetKey(void) const;

            //Returns the local matrix associated with the system
            const inline std::weak_ptr<ExpList> &GetLocMat(void) const;

            inline void InitObject();
            inline void Initialise(
                const std::shared_ptr<AssemblyMap>& pLocToGloMap);

            /// Solve the linear system for given input and output vectors
            /// using a specified local to global map.
            MULTI_REGIONS_EXPORT
            inline void Solve(
                const Array<OneD, const NekDouble> &in,
                      Array<OneD,       NekDouble> &out,
                const AssemblyMapSharedPtr         &locToGloMap,
                const Array<OneD, const NekDouble> &dirForcing
                    = NullNekDouble1DArray);

            /// Returns a shared pointer to the current object.
            std::shared_ptr<GlobalLinSys> GetSharedThisPtr()
            {
                return shared_from_this();
            }

            inline int                     GetNumBlocks      ();
            inline DNekScalMatSharedPtr    GetBlock          (unsigned int n);
            inline DNekScalBlkMatSharedPtr GetStaticCondBlock(unsigned int n);
            inline void                    DropStaticCondBlock(unsigned int n);

            /// Solve the linear system for given input and output vectors.
            inline void SolveLinearSystem(
                const int                          pNumRows,
                const Array<OneD,const NekDouble> &pInput,
                      Array<OneD,      NekDouble> &pOutput,
                const AssemblyMapSharedPtr        &locToGloMap,
                const int                          pNumDir = 0);

        protected:
            /// Key associated with this linear system.
            const GlobalLinSysKey                m_linSysKey;
            /// Local Matrix System
            const std::weak_ptr<ExpList>       m_expList;
            /// Robin boundary info
            const std::map<int, RobinBCInfoSharedPtr> m_robinBCInfo;
            // Provide verbose output
            bool                                 m_verbose;

            virtual int                     v_GetNumBlocks      ();
            virtual DNekScalMatSharedPtr    v_GetBlock          (unsigned int n);
            virtual DNekScalBlkMatSharedPtr v_GetStaticCondBlock(unsigned int n);
            virtual void                    v_DropStaticCondBlock(unsigned int n);

            PreconditionerSharedPtr CreatePrecon(AssemblyMapSharedPtr asmMap);

        private:
            LocalRegions::MatrixKey GetBlockMatrixKey(unsigned int n);
            
            /// Solve a linear system based on mapping.
            virtual void v_Solve(
                const Array<OneD, const NekDouble> &in,
                      Array<OneD,       NekDouble> &out,
                const AssemblyMapSharedPtr         &locToGloMap,
                const Array<OneD, const NekDouble> &dirForcing
                    = NullNekDouble1DArray) = 0;

            /// Solve a basic matrix system.
            virtual void v_SolveLinearSystem(
                const int                          pNumRows,
                const Array<OneD,const NekDouble> &pInput,
                      Array<OneD,      NekDouble> &pOutput,
                const AssemblyMapSharedPtr        &locToGloMap,
                const int                          pNumDir) = 0;

            virtual void v_InitObject();
            virtual void v_Initialise(
                const std::shared_ptr<AssemblyMap>& pLocToGloMap);

            static std::string lookupIds[];
            static std::string def;
        };


        /**
         *
         */
        const inline GlobalLinSysKey &GlobalLinSys::GetKey(void) const
        {
            return m_linSysKey;
        }

        /**
         *
         */
        const inline std::weak_ptr<ExpList> &GlobalLinSys::GetLocMat(void) const
        {
            return m_expList;
        }


        /**
         *
         */
        inline void GlobalLinSys::Solve(
                    const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const Array<OneD, const NekDouble> &dirForcing)
        {
            v_Solve(in,out,locToGloMap,dirForcing);
        }


        /**
         *
         */
        inline void GlobalLinSys::SolveLinearSystem(
                const int pNumRows,
                const Array<OneD,const NekDouble> &pInput,
                      Array<OneD,      NekDouble> &pOutput,
                const AssemblyMapSharedPtr &locToGloMap,
                const int pNumDir)
        {
            v_SolveLinearSystem(pNumRows, pInput, pOutput, locToGloMap, pNumDir);
        }

        inline void GlobalLinSys::InitObject()
        {
            v_InitObject();
        }

        inline void GlobalLinSys::Initialise(
            const std::shared_ptr<AssemblyMap>& pLocToGloMap)
        {
            v_Initialise(pLocToGloMap);
        }

        inline DNekScalMatSharedPtr GlobalLinSys::GetBlock(unsigned int n)
        {
            return v_GetBlock(n);
        }
        
        inline DNekScalBlkMatSharedPtr GlobalLinSys::GetStaticCondBlock(unsigned int n)
        {
            return v_GetStaticCondBlock(n);
        }

        inline void GlobalLinSys::DropStaticCondBlock(unsigned int n)
        {
            return v_DropStaticCondBlock(n);
        }

        inline int GlobalLinSys::GetNumBlocks()
        {
            return v_GetNumBlocks();
        }
    } //end of namespace
} //end of namespace

#endif
