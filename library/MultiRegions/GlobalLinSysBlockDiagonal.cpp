///////////////////////////////////////////////////////////////////////////////
//
// File GlobalVecLinSys.cpp
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
// Description: GlobalVecLinSysBlockDiagonal definition
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <MultiRegions/GlobalVecLinSysBlockDiagonal.h>


namespace Nektar
{
    namespace MultiRegions
    {

        /**
         * Registers the class with the Factory.
         */
        string GlobalVecLinSysDirectFull::className
                = GetGlobalVecLinSysFactory().RegisterCreatorFunction(
                    "BlockDiagonal",
                    GlobalLinSysDirectFull::create,
                    "BlockDiagonal.");
        
        /**
         * Given a block matrix, construct a global matrix system according to
         * a local to global mapping. #m_linSys is constructed by
         * AssembleFullMatrix().
         * @param   pkey        Associated linear system key.
         * @param   locToGloMap Local to global mapping.
         */
        GlobalVecLinSys::GlobalVecLinSys(const GlobalLinSysKey &pKey,
                const std::weak_ptr<VecExpList> &pExpList,
                const std::shared_ptr<AssemblyMap>
                                   &pLocToGloMap):
            m_linSysKey(pKey),
            m_expList(pExpList),
            m_verbose(true)
        {
            boost::ignore_unused(pLocToGloMap);
        }

        /**
         *
         */
        GlobalVecLinSysFactory& GetGlobalVecLinSysFactory()
        {
            static GlobalVecLinSysFactory instance;
            return instance;
        }

    } //end of namespace
} //end of namespace

