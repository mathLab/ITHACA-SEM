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
// Description: GlobalLinSysDirect header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSDIRECT_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSDIRECT_H

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/GlobalLinSys.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;

        /// A global linear system.
        class GlobalLinSysDirect : virtual public GlobalLinSys
        {
        public:
            MULTI_REGIONS_EXPORT GlobalLinSysDirect(
                const GlobalLinSysKey &pKey,
                const std::weak_ptr<ExpList> &pExp,
                const std::shared_ptr<AssemblyMap> &pLocToGloMap);
                
                MULTI_REGIONS_EXPORT virtual ~GlobalLinSysDirect();

        protected:
            /// Basic linear system object.
            DNekLinSysSharedPtr m_linSys;
            /// Solve the linear system for given input and output vectors.
            virtual void v_SolveLinearSystem(
                const int pNumRows,
                const Array<OneD,const NekDouble> &pInput,
                      Array<OneD,      NekDouble> &pOutput,
                const AssemblyMapSharedPtr &locToGloMap,
                const int pNumDir);
        };
    }
}

#endif
