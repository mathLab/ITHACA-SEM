///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysKeys.h
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
// Description: Headers for GlobalLinSysKey
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_GLOBALLINSYSKEY_H
#define NEKTAR_LIBS_MULTIREGIONS_GLOBALLINSYSKEY_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/GlobalMatrixKey.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /// Describe a linear system.
        class GlobalLinSysKey : public GlobalMatrixKey
        {
        public:
            MULTI_REGIONS_EXPORT GlobalLinSysKey(const StdRegions::MatrixType matrixType,
                            const AssemblyMapSharedPtr &locToGloMap
                                        = NullAssemblyMapSharedPtr,
                            const StdRegions::ConstFactorMap &factors = StdRegions::NullConstFactorMap,
                            const StdRegions::VarCoeffMap &varCoeffs = StdRegions::NullVarCoeffMap);

            /// Copy constructor.
            MULTI_REGIONS_EXPORT GlobalLinSysKey(const GlobalLinSysKey &key);
            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysKey();

            /// Less-than operator for GlobalLinSysKey comparison.
            MULTI_REGIONS_EXPORT friend bool operator<(const GlobalLinSysKey &lhs, 
                                  const GlobalLinSysKey &rhs);
            /// Return the associated solution type.
            inline GlobalSysSolnType  GetGlobalSysSolnType() const;

        protected:
            /// Store the solution type associated with the linear system. This
            /// may be none, full matrix, static condensation or multi-level
            /// static condensation.
            GlobalSysSolnType        m_solnType;
            
        private:
        
        };

        /// Writes information about the object to a given stream.
        MULTI_REGIONS_EXPORT std::ostream& operator<<(std::ostream& os, const GlobalLinSysKey& rhs);

        inline GlobalSysSolnType GlobalLinSysKey::GetGlobalSysSolnType() 
                                                                        const
        {
            return m_solnType; 
        }

    }
}

#endif
