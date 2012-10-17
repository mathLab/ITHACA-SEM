///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSys.cpp
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
// Description: GlobalLinSys definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GlobalLinSysDirect.h>
//#include <MultiRegions/AssemblyMapCG.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class GlobalLinSysDirect
         *
         * Solves a linear system using direct methods.
         */


//        GlobalLinSysDirect::GlobalLinSysDirect(void)
//        {
//        }

        /// Constructor for full direct matrix solve.
        GlobalLinSysDirect::GlobalLinSysDirect(
                const GlobalLinSysKey &pKey,
                const boost::weak_ptr<ExpList> &pExp,
                const boost::shared_ptr<AssemblyMap>
                                                        &pLocToGloMap)
                : GlobalLinSys(pKey, pExp, pLocToGloMap)
        {
        }

        GlobalLinSysDirect::~GlobalLinSysDirect()
        {
        }

        /// Solve the linear system for given input and output vectors.
        void GlobalLinSysDirect::v_SolveLinearSystem(
                const int pNumRows,
                const Array<OneD,const NekDouble> &pInput,
                      Array<OneD,      NekDouble> &pOutput,
                const AssemblyMapSharedPtr &pLocToGloMap,
                const int pNumDir)
        {
            DNekVec Vin(pInput.num_elements(),pInput);
            //  DNekVec Vout(pOutput.num_elements(),pOutput,eWrapper);
            ASSERTL1(pInput.num_elements() <= pOutput.num_elements(),
                     "output array must be at least as long as input array");
            DNekVec Vout(pInput.num_elements(),pOutput,eWrapper);
            m_linSys->Solve(Vin,Vout);
        }

        /// Solve the linear system for given input and output vectors
        /// using a specified local to global map.
        void GlobalLinSysDirect::v_Solve( const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out,
                    const AssemblyMapSharedPtr &pLocToGloMap,
                    const Array<OneD, const NekDouble> &pDirForcing)
        {
            ASSERTL0(false, "Not implemented for this GlobalLinSys type.");
        }
    }
}
