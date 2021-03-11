///////////////////////////////////////////////////////////////////////////////
//
// File PreconCfs.h
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
// Description: PreconCfs header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_PRECONCFS
#define NEKTAR_SOLVERS_COMPRESSIBLEFLOWSOLVER_PRECONCFS

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
    class PreconCfs
    {
    public:
        PreconCfs(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const LibUtilities::CommSharedPtr &vComm);

        virtual ~PreconCfs() {}

        inline void DoPreconCfs(
            const Array<OneD, NekDouble> &pInput,
            Array<OneD, NekDouble> &pOutput,
            const bool &flag);

        inline void BuildPreconCfs();

        inline void InitObject();

        inline bool UpdatePreconMatCheck();

    protected:
        // PreconCfsType                       m_preconType;
        LibUtilities::CommSharedPtr         m_comm;
        int                                 m_PreconMatFreezNumb;
    
        virtual void v_InitObject();

    private:

        void DoNullPrecon(
            const Array<OneD, NekDouble> &pInput,
            Array<OneD, NekDouble> &pOutput,
            const bool &flag);


        virtual void v_DoPreconCfs(
            const Array<OneD, NekDouble> &pInput,
            Array<OneD, NekDouble> &pOutput,
            const bool &flag);

        virtual void v_BuildPreconCfs();
    };
    typedef std::shared_ptr<PreconCfs>  PreconCfsSharedPtr;

    /**
     *
     */
    inline void PreconCfs::InitObject()
    {
        v_InitObject();
    }

    /**
     *
     */
    inline void PreconCfs::BuildPreconCfs()
    {
        v_BuildPreconCfs();
    }
}

#endif
