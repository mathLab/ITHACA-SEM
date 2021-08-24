///////////////////////////////////////////////////////////////////////////////
//
// File: GlobalLinSysStaticCond.h
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
// Description: A collection of routines common to statically condensed systems.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSSTATICCOND_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSSTATICCOND_H

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/GlobalLinSysIterative.h>
#include <LibUtilities/LinearAlgebra/SparseMatrixFwd.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;
        class GlobalLinSysStaticCond;

        typedef std::shared_ptr<GlobalLinSysStaticCond>
            GlobalLinSysStaticCondSharedPtr;

        /// A global linear system.
        class GlobalLinSysStaticCond : virtual public GlobalLinSys
        {
        public:
            /// Constructor for full direct matrix solve.
            GlobalLinSysStaticCond(
                const GlobalLinSysKey                &mkey,
                const std::weak_ptr<ExpList>         &pExpList,
                const std::shared_ptr<AssemblyMap>   &locToGloMap);

            virtual ~GlobalLinSysStaticCond();

        protected:
            virtual void v_PreSolve(int                   scLevel,
                                    Array<OneD,NekDouble> &F_bnd)
            {
                boost::ignore_unused(scLevel, F_bnd);
            }

            virtual void v_BasisFwdTransform(
                Array<OneD, NekDouble>& pInOut)
            {
                boost::ignore_unused(pInOut);
            }

            virtual void v_CoeffsBwdTransform(
                Array<OneD, NekDouble>& pInOut)
            {
                boost::ignore_unused(pInOut);
            }

            virtual void v_CoeffsFwdTransform(
                const Array<OneD, NekDouble>& pInput,
                Array<OneD, NekDouble>& pOutput)
            {
                boost::ignore_unused(pInput,pOutput);
            }

            virtual void v_AssembleSchurComplement(
                std::shared_ptr<AssemblyMap> pLoctoGloMap)
            {
                boost::ignore_unused(pLoctoGloMap);
            }

            virtual int v_GetNumBlocks();

            virtual GlobalLinSysStaticCondSharedPtr v_Recurse(
                const GlobalLinSysKey                &mkey,
                const std::weak_ptr<ExpList>         &pExpList,
                const DNekScalBlkMatSharedPtr         pSchurCompl,
                const DNekScalBlkMatSharedPtr         pBinvD,
                const DNekScalBlkMatSharedPtr         pC,
                const DNekScalBlkMatSharedPtr         pInvD,
                const std::shared_ptr<AssemblyMap>   &locToGloMap) = 0;
            
            /// Schur complement for Direct Static Condensation.
            GlobalLinSysStaticCondSharedPtr          m_recursiveSchurCompl;
            /// Block Schur complement matrix.
            DNekScalBlkMatSharedPtr                  m_schurCompl;
            /// Block \f$ BD^{-1} \f$ matrix.
            DNekScalBlkMatSharedPtr                  m_BinvD;
            /// Block \f$ C \f$ matrix.
            DNekScalBlkMatSharedPtr                  m_C;
            /// Block \f$ D^{-1} \f$ matrix.
            DNekScalBlkMatSharedPtr                  m_invD;
            /// Local to global map.
            std::weak_ptr<AssemblyMap>               m_locToGloMap;
            /// Workspace array for matrix multiplication
            Array<OneD, NekDouble>                   m_wsp;

            /// Solve the linear system for given input and output vectors
            /// using a specified local to global map.
            virtual void v_Solve(
                const Array<OneD, const NekDouble> &in,
                      Array<OneD,       NekDouble> &out,
                const AssemblyMapSharedPtr         &locToGloMap,
                const Array<OneD, const NekDouble> &dirForcing
                    = NullNekDouble1DArray);

            virtual void v_InitObject();

            /// Initialise this object
            virtual void v_Initialise(
                    const std::shared_ptr<AssemblyMap>& locToGloMap);

            /// Set up the storage for the Schur complement or the top level
            /// of the multi-level Schur complement.
            void SetupTopLevel(
                    const std::shared_ptr<AssemblyMap>& locToGloMap);

            ///
            void ConstructNextLevelCondensedSystem(
                    const std::shared_ptr<AssemblyMap>& locToGloMap);
        };
    }
}

#endif
