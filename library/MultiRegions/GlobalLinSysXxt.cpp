///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysXxt.cpp
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
// Description: GlobalLinSysXxt definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Communication/Xxt.hpp>
#include <MultiRegions/GlobalLinSysXxt.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysXxt
         *
         * Solves a linear system using direct methods.
         */

        /// Constructor for full direct matrix solve.
        GlobalLinSysXxt::GlobalLinSysXxt(
                const GlobalLinSysKey &pKey,
                const boost::weak_ptr<ExpList> &pExp,
                const boost::shared_ptr<AssemblyMap>
                                                        &pLocToGloMap)
                : GlobalLinSys(pKey, pExp, pLocToGloMap)
        {
        }

        GlobalLinSysXxt::~GlobalLinSysXxt()
        {
            Xxt::Finalise(m_crsData);
        }

        /// Solve the linear system for given input and output vectors.
        void GlobalLinSysXxt::v_SolveLinearSystem(
                const int pNumRows,
                const Array<OneD,const NekDouble> &pInput,
                      Array<OneD,      NekDouble> &pOutput,
                const AssemblyMapSharedPtr &pLocToGloMap,
                const int pNumDir)
        {
            unsigned int nCoeffs = pLocToGloMap->GetNumLocalCoeffs();
            Array<OneD, NekDouble> vLocalIn(m_rank, 0.0);
            Array<OneD, NekDouble> vLocalOut(m_rank, 0.0);
            GlobalToLocalNonDir(pInput, vLocalIn, pLocToGloMap);
            Xxt::Solve(vLocalOut, m_crsData, vLocalIn);
            LocalNonDirToGlobal(vLocalOut, pOutput, pLocToGloMap);
        }

        /// Solve the linear system for given input and output vectors
        /// using a specified local to global map.
        void GlobalLinSysXxt::v_Solve( const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out,
                    const AssemblyMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble> &pDirForcing)
        {
            ASSERTL0(false, "Not implemented for this GlobalLinSys type.");
        }

        void GlobalLinSysXxt::GlobalToLocalNonDir(const Array<OneD, const NekDouble> &global,
                                       Array<OneD, NekDouble> &local,
                               const boost::shared_ptr<AssemblyMap>
                                                      &pLocToGloMap)
        {
            int n = m_locToGloMap.num_elements();
            Vmath::Gathr(n, m_locToGloSignMult.get(), global.get(), m_locToGloMap.get(), local.get());
        }

        void GlobalLinSysXxt::LocalNonDirToGlobal(const Array<OneD, const NekDouble> &local,
                                       Array<OneD, NekDouble> &global,
                               const boost::shared_ptr<AssemblyMap>
                                                      &pLocToGloMap)
        {
            int n = m_locToGloMap.num_elements();
            Vmath::Scatr(n, m_locToGloSign.get(), local.get(), m_locToGloMap.get(), global.get());
        }

    }
}

