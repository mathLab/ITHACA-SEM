///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSys.h
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
// Description: GlobalLinSysStaticCond header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSDIRECTSTATICCOND_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSDIRECTSTATICCOND_H

#include <MultiRegions/GlobalLinSysDirect.h>
#include <MultiRegions/GlobalLinSysStaticCond.h>
#include <MultiRegions/MultiRegionsDeclspec.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;
        class GlobalLinSysDirectStaticCond;

        typedef std::shared_ptr<GlobalLinSysDirectStaticCond>
                                        GlobalLinSysDirectStaticCondSharedPtr;

        /// A global linear system.
        class GlobalLinSysDirectStaticCond : virtual public GlobalLinSysDirect,
                                             virtual public GlobalLinSysStaticCond
        {
        public:
            /// Creates an instance of this class
            static GlobalLinSysSharedPtr create(
                        const GlobalLinSysKey                &pLinSysKey,
                        const std::weak_ptr<ExpList>         &pExpList,
                        const std::shared_ptr<AssemblyMap>   &pLocToGloMap)
            {
                GlobalLinSysDirectStaticCondSharedPtr ret =
                    MemoryManager<GlobalLinSysDirectStaticCond>
                    ::AllocateSharedPtr(pLinSysKey, pExpList, pLocToGloMap);
                ret->InitObject();
                return ret;
            }

            /// Name of class
            MULTI_REGIONS_EXPORT static std::string className;
            static std::string className2;

            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysDirectStaticCond(
                        const GlobalLinSysKey                &mkey,
                        const std::weak_ptr<ExpList>         &pExpList,
                        const std::shared_ptr<AssemblyMap>   &locToGloMap);

            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysDirectStaticCond(
                        const GlobalLinSysKey                &mkey,
                        const std::weak_ptr<ExpList>         &pExpList,
                        const DNekScalBlkMatSharedPtr         pSchurCompl,
                        const DNekScalBlkMatSharedPtr         pBinvD,
                        const DNekScalBlkMatSharedPtr         pC,
                        const DNekScalBlkMatSharedPtr         pInvD,
                        const std::shared_ptr<AssemblyMap>   &locToGloMap);

            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysDirectStaticCond();

        protected:
            virtual void v_AssembleSchurComplement(
                std::shared_ptr<AssemblyMap> pLocToGloMap);
            virtual GlobalLinSysStaticCondSharedPtr v_Recurse(
                const GlobalLinSysKey                &mkey,
                const std::weak_ptr<ExpList>         &pExpList,
                const DNekScalBlkMatSharedPtr         pSchurCompl,
                const DNekScalBlkMatSharedPtr         pBinvD,
                const DNekScalBlkMatSharedPtr         pC,
                const DNekScalBlkMatSharedPtr         pInvD,
                const std::shared_ptr<AssemblyMap>   &l2gMap);

        private:
            /// Matrix Storage type for known matrices
            MatrixStorage DetermineMatrixStorage(
                   const std::shared_ptr<AssemblyMap>& locToGloMap);
        };
    }
}

#endif
