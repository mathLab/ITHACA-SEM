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
// Description: GlobalLinSysStaticCond header
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSDIRECTSTATICCOND_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSDIRECTSTATICCOND_H

#include <MultiRegions/GlobalLinSysKey.h>
#include <MultiRegions/GlobalLinSysDirect.h>
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;
        class GlobalLinSysDirectStaticCond;

        typedef boost::shared_ptr<GlobalLinSysDirectStaticCond>
                                        GlobalLinSysDirectStaticCondSharedPtr;

        /// A global linear system.
        class GlobalLinSysDirectStaticCond : public GlobalLinSysDirect
        {
        public:
            /// Creates an instance of this class
            static GlobalLinSysSharedPtr create(
                        const GlobalLinSysKey &pLinSysKey,
                        const boost::weak_ptr<ExpList> &pExpList,
                        const boost::shared_ptr<AssemblyMap>
                                                               &pLocToGloMap)
            {
                return MemoryManager<GlobalLinSysDirectStaticCond>
                    ::AllocateSharedPtr(pLinSysKey, pExpList, pLocToGloMap);
            }

            /// Name of class
            MULTI_REGIONS_EXPORT static std::string className;
            static std::string className2;

            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysDirectStaticCond(
                        const GlobalLinSysKey &mkey,
                        const boost::weak_ptr<ExpList> &pExpList,
                        const boost::shared_ptr<AssemblyMap>
                                                                &locToGloMap);

            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysDirectStaticCond(
                        const GlobalLinSysKey &mkey,
                        const boost::weak_ptr<ExpList> &pExpList,
                        const DNekScalBlkMatSharedPtr pSchurCompl,
                        const DNekScalBlkMatSharedPtr pBinvD,
                        const DNekScalBlkMatSharedPtr pC,
                        const DNekScalBlkMatSharedPtr pInvD,
                        const boost::shared_ptr<AssemblyMap>
                                                                &locToGloMap);

//            MULTI_REGIONS_EXPORT GlobalLinSysDirectStaticCond(
//                        const DNekScalBlkMatSharedPtr pSchurCompl,
//                        const DNekScalBlkMatSharedPtr pBinvD,
//                        const DNekScalBlkMatSharedPtr pC,
//                        const DNekScalBlkMatSharedPtr pInvD,
//                        const boost::shared_ptr<AssemblyMap>
//                        &pLocToGloMap);

            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysDirectStaticCond();

        private:
            /// Schur complement for Direct Static Condensation.
            GlobalLinSysDirectStaticCondSharedPtr m_recursiveSchurCompl;

            /// Block matrices at this level
            DNekScalBlkMatSharedPtr m_schurCompl;
            DNekScalBlkMatSharedPtr m_BinvD;
            DNekScalBlkMatSharedPtr m_C;
            DNekScalBlkMatSharedPtr m_invD;

            /// Solve the linear system for given input and output vectors
            /// using a specified local to global map.
            virtual void v_Solve(
                        const Array<OneD, const NekDouble>  &in,
                              Array<OneD,       NekDouble>  &out,
                        const AssemblyMapSharedPtr &locToGloMap,
                        const Array<OneD, const NekDouble>  &dirForcing
                                                        = NullNekDouble1DArray);

            /// Initialise this object
            void Initialise(
                  const boost::shared_ptr<AssemblyMap>& locToGloMap,
                  MatrixStorage matStorage);

            /// Matrix Storage type for known matrices
            MatrixStorage DetermineMatrixStorage(
                   const boost::shared_ptr<AssemblyMap>& locToGloMap);

            /// Set up the storage for the Schur complement or the top level
            /// of the multi-level Schur complement.
            void SetupTopLevel(
                    const boost::shared_ptr<AssemblyMap>& locToGloMap);

            /// Assemble the Schur complement matrix.
            void AssembleSchurComplement(
                                         const boost::shared_ptr<AssemblyMap>& locToGloMap,
                                         const MatrixStorage matStorage);

            ///
            void ConstructNextLevelCondensedSystem(
                    const boost::shared_ptr<AssemblyMap>& locToGloMap);
         };
    }
}

#endif
