///////////////////////////////////////////////////////////////////////////////
//
// File:  PreconCfs.cpp
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
// Description:  PreconCfs definition
//
///////////////////////////////////////////////////////////////////////////////

#include <CompressibleFlowSolver/Preconditioner/PreconCfs.h>

using namespace std;

namespace Nektar
{
    PreconCfs::PreconCfs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const LibUtilities::CommSharedPtr &vComm)
    {
        m_Comm      = vComm;
        m_verbose   = false;
        m_root      = false;

        if (0 == m_Comm->GetRank())
        {
            m_root = true;
        }
        m_verbose = pSession->DefinesCmdLineArgument("verbose");

        m_spacedim = pFields[0]->GetGraph()->GetSpaceDimension();
        pSession->LoadParameter("PreconMatFreezNumb",     
            m_PreconMatFreezNumb, 200);
    }
    void PreconCfs::v_InitObject()
    {

    }

    void PreconCfs::DoNullPrecon(
            const Array<OneD, NekDouble> &pInput,
            Array<OneD, NekDouble> &pOutput,
            const bool &flag)
    {
        boost::ignore_unused(flag);
        Vmath::Vcopy(pInput.size(), pInput, 1, pOutput, 1);
    }

    /**
     *
     */
    void PreconCfs::DoPreconCfs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, NekDouble> &pInput,
        Array<OneD, NekDouble> &pOutput,
        const bool &flag)
    {
        ASSERTL0(pInput.size() == pOutput.size(), 
            "In and Out not the same size in DoNullPrecon");
        v_DoPreconCfs(pFields, pInput, pOutput, flag);
        m_PreconTimesCounter ++;
    }


    void PreconCfs::v_DoPreconCfs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, NekDouble> &pInput,
        Array<OneD, NekDouble> &pOutput,
        const bool &flag)
    {
        boost::ignore_unused(pFields, pInput, pOutput, flag);
        NEKERROR(ErrorUtil::efatal, "v_DoPreconCfs not defined");
    }

    void PreconCfs::v_BuildPreconCfs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, const Array<OneD, NekDouble>>   &intmp,
        const NekDouble                                   time,
        const NekDouble                                   lambda)
    {
        boost::ignore_unused(pFields, intmp, time, lambda);
        NEKERROR(ErrorUtil::efatal, "v_BuildPreconCfs not defined");
    }

    bool PreconCfs::UpdatePreconMatCheck(
            const Array<OneD, const NekDouble>  &res,
            const NekDouble                     dtLambda)
    {
        NEKERROR(ErrorUtil::efatal, "UpdatePreconMatCheck not defined");
    }

} // namespace Nektar
