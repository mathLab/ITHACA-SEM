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
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <MultiRegions/GlobalLinSysKey.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class ExpList;
        class GlobalLinSys;

        /// Pointer to a GlobalLinSys object.
        typedef boost::shared_ptr<GlobalLinSys> GlobalLinSysSharedPtr;
        /// Mapping between GlobalLinSys objects and their associated keys.
        typedef map<GlobalLinSysKey,GlobalLinSysSharedPtr> GlobalLinSysMap;
        /// Pointer to a GlobalLinSys/key map.
        typedef boost::shared_ptr<GlobalLinSysMap> GlobalLinSysMapShPtr;

        /// Datatype of the NekFactory used to instantiate classes derived from
        /// the EquationSystem class.
        typedef LibUtilities::NekFactory< std::string, GlobalLinSys, 
            const GlobalLinSysKey&,
            const boost::shared_ptr<ExpList>&,
            const boost::shared_ptr<LocalToGlobalBaseMap>& > GlobalLinSysFactory;
        GlobalLinSysFactory& GetGlobalLinSysFactory();


        /// A global linear system.
        class GlobalLinSys
        {
        public:
            /// Default constructor
            MULTI_REGIONS_EXPORT
            GlobalLinSys();

            /// Constructor for full direct matrix solve.
            MULTI_REGIONS_EXPORT
            GlobalLinSys(const GlobalLinSysKey &pKey,
                    const boost::shared_ptr<ExpList> &pExpList,
                    const boost::shared_ptr<LocalToGlobalBaseMap>
                                                           &pLocToGloMap);

            /// Returns the key associated with the system.
            const inline GlobalLinSysKey &GetKey(void) const;

            /// Solve the linear system for given input and output vectors
            /// using a specified local to global map.
            MULTI_REGIONS_EXPORT
            inline void Solve(
                    const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out,
                    const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                    const Array<OneD, const NekDouble> &dirForcing
                                                = NullNekDouble1DArray);

        protected:
            /// Key associated with this linear system.
            GlobalLinSysKey                         m_linSysKey;
            /// Local Matrix System
            boost::shared_ptr<ExpList>              m_expList;

            /// Solve the linear system for given input and output vectors.
            inline void SolveLinearSystem(
                    const int pNumRows,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const int pNumDir = 0);

            DNekScalMatSharedPtr GetBlock(unsigned int n);
            DNekScalBlkMatSharedPtr GetStaticCondBlock(unsigned int n);

        private:
            /// Solve a linear system based on mapping.
            virtual void v_Solve(
                        const Array<OneD, const NekDouble> &in,
                              Array<OneD,       NekDouble> &out,
                        const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                        const Array<OneD, const NekDouble> &dirForcing
                                                = NullNekDouble1DArray) = 0;

            /// Solve a basic matrix system.
            virtual void v_SolveLinearSystem(
                    const int pNumRows,
                    const Array<OneD,const NekDouble> &pInput,
                          Array<OneD,      NekDouble> &pOutput,
                    const int pNumDir) = 0;
        };


        /**
         *
         */
        const inline GlobalLinSysKey &GlobalLinSys::GetKey(void) const
        {
            return m_linSysKey;
        }


        /**
         *
         */
        inline void GlobalLinSys::Solve(
                    const Array<OneD, const NekDouble> &in,
                          Array<OneD,       NekDouble> &out,
                    const LocalToGlobalBaseMapSharedPtr &locToGloMap,
                    const Array<OneD, const NekDouble> &dirForcing)
        {
            v_Solve(in,out,locToGloMap,dirForcing);
        }


        /**
         *
         */
        inline void GlobalLinSys::SolveLinearSystem(
                const int pNumRows,
                const Array<OneD,const NekDouble> &pInput,
                      Array<OneD,      NekDouble> &pOutput,
                const int pNumDir)
        {
            v_SolveLinearSystem(pNumRows, pInput, pOutput, pNumDir);
        }

    } //end of namespace
} //end of namespace

#endif
