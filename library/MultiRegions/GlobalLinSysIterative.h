///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysIterative.h
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
// Description: GlobalLinSysIterative header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVE_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSITERATIVE_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/GlobalLinSysKey.h>
#include <MultiRegions/GlobalLinSys.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class LocalToGlobalC0ContMap;
        class ExpList;

        /// A global linear system.
        class GlobalLinSysIterative : public GlobalLinSys
        {
        public:
            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysIterative(const GlobalLinSysKey &pKey,
                         const boost::shared_ptr<ExpList> &pExpList,
                         const boost::shared_ptr<LocalToGlobalBaseMap>
                                                                &pLocToGloMap);

            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysIterative();

            /// Solve the linear system for given input and output vectors.
            MULTI_REGIONS_EXPORT virtual void Solve( const Array<OneD,const NekDouble> &in,
                              Array<OneD,      NekDouble> &out);

            /// Solve the linear system for given input and output vectors
            /// using a specified local to global map.
            MULTI_REGIONS_EXPORT virtual void Solve( const Array<OneD, const NekDouble> &in,
                              Array<OneD,       NekDouble> &out,
                        const LocalToGlobalBaseMapSharedPtr &pLocToGloMap,
                        const Array<OneD, const NekDouble> &pDirForcing
                                                        = NullNekDouble1DArray);

        protected:

        };
    }
}

#endif
