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
// Description: GlobalLinSys header 
//
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYS_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYS_H

#include <MultiRegions/GlobalLinSysKey.h>

namespace Nektar
{
    namespace MultiRegions
    {
        class LocalToGlobalC0ContMap;
	class ExpList;

	class GlobalLinSys
        {
        public:
            GlobalLinSys(const GlobalLinSysKey &mkey, 
                         const DNekLinSysSharedPtr linSys);

            GlobalLinSys(const GlobalLinSysKey &mkey, 
                         const DNekLinSysSharedPtr linsys,
                         const DNekScalBlkMatSharedPtr SchurCompl,
                         const DNekScalBlkMatSharedPtr BinvD,
                         const DNekScalBlkMatSharedPtr invDC,
                         const DNekScalBlkMatSharedPtr invD);


            const GlobalLinSysKey &GetKey(void) const
            {
                return m_linSysKey;
            }
            
            const DNekLinSysSharedPtr GetLinSys(void) const 
            {
                return m_linSys;
            }

            void Solve(const Array<OneD,const NekDouble> &in, 
                             Array<OneD,      NekDouble> &out);

            void Solve(const Array<OneD, const NekDouble> &in, 
                             Array<OneD,       NekDouble> &out,
                             LocalToGlobalC0ContMap &locToGloMap,
                             ExpList* exp,
                       const Array<OneD, const NekDouble> &dirForcing = NullNekDouble1DArray);

        private:
            GlobalLinSysKey                     m_linSysKey;
            DNekLinSysSharedPtr                 m_linSys;
            Array<OneD,DNekScalBlkMatSharedPtr> m_blkMatrices;
	};

        typedef boost::shared_ptr<GlobalLinSys> GlobalLinSysSharedPtr;

        typedef map<GlobalLinSysKey,GlobalLinSysSharedPtr> GlobalLinSysMap;
        typedef boost::shared_ptr<GlobalLinSysMap> GlobalLinSysMapShPtr; 

    } //end of namespace
} //end of namespace

#endif
