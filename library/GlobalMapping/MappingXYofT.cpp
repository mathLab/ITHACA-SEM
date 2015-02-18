///////////////////////////////////////////////////////////////////////////////
//
// File: MappingXYofT.cpp
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
// Description: Mapping of the type X = x + f(t), Y = y + g(t)
//
///////////////////////////////////////////////////////////////////////////////

#include <GlobalMapping/MappingXYofT.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace GlobalMapping
{

    std::string MappingXYofT::className =
            GetMappingFactory().RegisterCreatorFunction("XYofT",
                    MappingXYofT::create, "X = x + f(t), Y = y +g(t)");

    /**
     *
     */
    MappingXYofT::MappingXYofT(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const Array<OneD, MultiRegions::ExpListSharedPtr>&   pFields)
        : MappingIdentity(pSession, pFields)
    {
    }


    /**
     *
     */
    void MappingXYofT::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const TiXmlElement                                *pMapping)
    {
        MappingIdentity::v_InitObject(pFields, pMapping);
        
        m_timeDependent    = true;
        
        ASSERTL0(m_nConvectiveFields>=2,
               "Mapping X = x + f(t), Y = y+g(t) needs 2 velocity components.");
    }

    void MappingXYofT::v_UpdateMapping(const NekDouble time)
    {
        std::string s_XFieldStr = m_session->GetVariable(0);
        std::string s_YFieldStr = m_session->GetVariable(1);
        // Evaluate x-functions
        EvaluateFunction(m_fields, m_session, s_XFieldStr, m_coords[0],
                m_funcName, time); 
        EvaluateFunction(m_fields, m_session, s_XFieldStr, m_coordsVel[0],
                m_velFuncName, time);
        
        // Evaluate y-functions 
        EvaluateFunction(m_fields, m_session, s_YFieldStr, m_coords[1],
                m_funcName, time);
        EvaluateFunction(m_fields, m_session, s_YFieldStr, m_coordsVel[1],
                m_velFuncName, time);
    }

}
}
