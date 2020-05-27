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

#include <boost/core/ignore_unused.hpp>

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
                const std::weak_ptr<ExpList> &pExp,
                const std::shared_ptr<AssemblyMap>
                                                        &pLocToGloMap)
                : GlobalLinSys(pKey, pExp, pLocToGloMap)
        {
            m_crsData = 0;
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
            boost::ignore_unused(pNumRows, pLocToGloMap, pNumDir);

            int nLocal = m_map.size();
            Array<OneD, NekDouble> vLocalIn(nLocal, 0.0);
            Array<OneD, NekDouble> vLocalOut(nLocal, 0.0);
            GlobalToLocalNoSign(pInput, vLocalIn);
            Xxt::Solve(vLocalOut, m_crsData, vLocalIn);
            LocalToGlobalNoSign(vLocalOut, pOutput);
        }

        void GlobalLinSysXxt::GlobalToLocalNoSign(const Array<OneD, const NekDouble> &global,
                                       Array<OneD, NekDouble> &local)
        {
            Vmath::Gathr(m_map.size(), m_locToGloSignMult.get(), global.get(), m_map.get(), local.get());
        }

        void GlobalLinSysXxt::LocalToGlobalNoSign(const Array<OneD, const NekDouble> &local,
                                       Array<OneD, NekDouble> &global)
        {
            Vmath::Scatr(m_map.size(), local.get(), m_map.get(), global.get());
        }

    }
}

