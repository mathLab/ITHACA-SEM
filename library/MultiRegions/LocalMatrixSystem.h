///////////////////////////////////////////////////////////////////////////////
//
// File LocalMatrixSystem.h
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
// Description: LocalMatrixSystem definition
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_LOCALMATRIXSYSTEM_H
#define NEKTAR_LIB_MULTIREGIONS_LOCALMATRIXSYSTEM_H

#include "StdRegions/StdExpansion.h"
#include "MultiRegions/GlobalLinSysKey.h"

namespace Nektar
{
    namespace MultiRegions
    {
        class LocalMatrixSystem
        {
        public:
            virtual ~LocalMatrixSystem() {}

            virtual GlobalSysSolnType GetSystemType() = 0;

            Array<OneD, DNekScalBlkMatSharedPtr> GetLocalSystem();
            void DeallocateSystem(unsigned int pIndex);

        protected:
            Array<OneD, DNekScalBlkMatSharedPtr> m_sys;

            LocalMatrixSystem();
            LocalMatrixSystem(const GlobalLinSysKey &mkey,
                        boost::shared_ptr<StdRegions::StdExpansionVector>&,
                        const boost::shared_ptr<LocalToGlobalBaseMap>
                                                                &pLocToGloMap,
                        const Array<OneD, int>& pOffsets,
                        const map<int, RobinBCInfoSharedPtr>& pRobinBCInfo);
        };
        /// Pointer to a GlobalLinSys object.
        typedef boost::shared_ptr<LocalMatrixSystem> LocalMatrixSystemSharedPtr;


        class LocalMatrixSystemFull : public LocalMatrixSystem
        {
        public:
            LocalMatrixSystemFull(const GlobalLinSysKey &mkey,
                        boost::shared_ptr<StdRegions::StdExpansionVector>&,
                        const boost::shared_ptr<LocalToGlobalBaseMap>
                                                                &pLocToGloMap,
                        const Array<OneD, int>& pOffsets,
                        const map<int, RobinBCInfoSharedPtr>& pRobinBCInfo);
            virtual ~LocalMatrixSystemFull() {}
            virtual GlobalSysSolnType GetSystemType()
            {
            	return eDirectFullMatrix;
            }
        };
        /// Pointer to a GlobalLinSys object.
        typedef boost::shared_ptr<LocalMatrixSystemFull> LocalMatrixSystemFullSharedPtr;


        class LocalMatrixSystemStaticCond : public LocalMatrixSystem
        {
        public:
            LocalMatrixSystemStaticCond(
                        const boost::shared_ptr<LocalToGlobalBaseMap>
                                                                &pLocToGloMap);
            LocalMatrixSystemStaticCond(const GlobalLinSysKey &mkey,
                        boost::shared_ptr<StdRegions::StdExpansionVector>&,
                        const boost::shared_ptr<LocalToGlobalBaseMap>
                                                                &pLocToGloMap,
                        const Array<OneD, int>& pOffsets,
                        const map<int, RobinBCInfoSharedPtr>& pRobinBCInfo);
            virtual ~LocalMatrixSystemStaticCond() {}
            virtual GlobalSysSolnType GetSystemType()
            {
            	return eDirectStaticCond;
            }

        private:
            void Setup(const GlobalLinSysKey &mkey,
                    boost::shared_ptr<StdRegions::StdExpansionVector>&,
                    const boost::shared_ptr<LocalToGlobalBaseMap>
                                                            &pLocToGloMap,
                    const Array<OneD, int>& pOffsets,
                    const map<int, RobinBCInfoSharedPtr>& pRobinBCInfo);
            void SetupBnd(const GlobalLinSysKey &mkey,
                        boost::shared_ptr<StdRegions::StdExpansionVector>&,
                        const boost::shared_ptr<LocalToGlobalBaseMap>
                                                                &pLocToGloMap,
                        const Array<OneD, int>& pOffsets,
                        const map<int, RobinBCInfoSharedPtr>& pRobinBCInfo);
        };
        /// Pointer to a GlobalLinSys object.
        typedef boost::shared_ptr<LocalMatrixSystemStaticCond> LocalMatrixSystemStaticCondSharedPtr;

    }
}

#endif /* LOCALMATRIXSYSTEM_H_ */
