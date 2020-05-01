///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysDirectXxt.h
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
// Description: GlobalLinSysDirectXxt header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSPETSCFULL_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSPETSCFULL_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/GlobalLinSysPETSc.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations

        //class AssemblyMapDG;
        class ExpList;

        /// A global linear system.
        class GlobalLinSysPETScFull : public GlobalLinSysPETSc
        {
        public:

            /// Creates an instance of this class
            static GlobalLinSysSharedPtr create(
                const GlobalLinSysKey                &pLinSysKey,
                const std::weak_ptr<ExpList>         &pExpList,
                const std::shared_ptr<AssemblyMap>   &pLocToGloMap)
            {
                return MemoryManager<GlobalLinSysPETScFull>
                    ::AllocateSharedPtr(pLinSysKey, pExpList, pLocToGloMap);
            }

            /// Name of class
            MULTI_REGIONS_EXPORT static std::string className;

            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysPETScFull(
                const GlobalLinSysKey                &pLinSysKey,
                const std::weak_ptr<ExpList>         &pExpList,
                const std::shared_ptr<AssemblyMap>   &pLocToGloMap);

            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysPETScFull();

        private:
            // Local to global map.
            std::shared_ptr<AssemblyMap>     m_locToGloMap;

            /// Solve the linear system for given input and output vectors
            /// using a specified local to global map.
            virtual void v_Solve(
                const Array<OneD, const NekDouble> &in,
                      Array<OneD,       NekDouble> &out,
                const AssemblyMapSharedPtr         &locToGloMap,
                const Array<OneD, const NekDouble> &dirForcing
                                                        = NullNekDouble1DArray);
            virtual void v_DoMatrixMultiply(
                const Array<OneD, const NekDouble> &input,
                      Array<OneD,       NekDouble> &output);
        };
    }
}

#endif
