///////////////////////////////////////////////////////////////////////////////
//
// File GlobalLinSysIterativeSolver.h
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
// Description: GlobalLinSysIterativeSolver header: just to solve a linear sysytem input by a xml file
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSIterativeCG_H
#define NEKTAR_LIB_MULTIREGIONS_GLOBALLINSYSIterativeCG_H

#include <MultiRegions/GlobalLinSysIterative.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;

        /// A global linear system.
        class GlobalLinSysIterativeSolver : public GlobalLinSysIterative
        {
        public:
            /// Creates an instance of this class
            // static GlobalLinSysSharedPtr create(
            //         const LibUtilities::SessionReaderSharedPtr &session)
            // {
            //     return MemoryManager<GlobalLinSysIterativeSolver>
            //         ::AllocateSharedPtr(session);
            // }

            /// Name of class
            // static std::string className;

            // /// Constructor for full direct matrix solve.
            // MULTI_REGIONS_EXPORT GlobalLinSysIterativeSolver();


            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT GlobalLinSysIterativeSolver(
                    const GlobalLinSysKey &pLinSysKey,
                    const std::weak_ptr<ExpList> &pExpList,
                    const std::shared_ptr<AssemblyMap>
                                                           &pLocToGloMap);
                

            MULTI_REGIONS_EXPORT virtual ~GlobalLinSysIterativeSolver();

            void initializeLinSys(const LibUtilities::SessionReaderSharedPtr &session);

            inline void setdimension(const unsigned int ndimens)
            {
                m_nlinsys = ndimens;
                return;
            }

            inline int getdimension()
            {
                return m_nlinsys;
            }

            inline Array<OneD,NekDouble> &getrhs()
            {
                return m_rhs;
            }
        protected:
            
            // stores the dimension of the linear sysytem
            unsigned int    m_nlinsys;
            
            // stores the A of the linear system Ax = f
            Array<OneD, NekDouble> m_mat;

            // stores the f of the linear system Ax = f
            Array<OneD,       NekDouble> m_rhs;

        private:

            // Local to global map.
            std::shared_ptr<AssemblyMap>     m_locToGloMap;

            virtual void v_DoMatrixMultiply(
                    const Array<OneD, NekDouble>& pInput,
                          Array<OneD, NekDouble>& pOutput);


            virtual void v_Solve(
                    const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out,
                    const AssemblyMapSharedPtr &locToGloMap,
                    const Array<OneD, const NekDouble> &dirForcing
                                                        = NullNekDouble1DArray);
            virtual void v_UniqueMap();



        };
    }
}

#endif
