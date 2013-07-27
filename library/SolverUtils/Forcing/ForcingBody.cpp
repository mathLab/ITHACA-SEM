///////////////////////////////////////////////////////////////////////////////
//
// File: Forcing.cpp
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
// Description: Abstract base class for forcing terms.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingBody.h>

namespace Nektar
{
    namespace SolverUtils
    {
	
	std::string ForcingBody::className = GetForcingFactory().RegisterCreatorFunction("BodyForcing", ForcingBody::create, "Body Forcing");

	ForcingBody::ForcingBody()
	{
	}
        void ForcingBody::v_InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields,
                SpatialDomains::MeshGraphSharedPtr                pGraph)
        {
		m_Session = pSession; 
                v_ReadForceInfo(pSession,pFields,pGraph);
        }
	void ForcingBody::v_Apply(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray)
	{
		if (m_Session->DefinesFunction("BodyForce"))
		{
			for (int i=0; i<m_NumVariable;i++)
			{	
				Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1, m_Forcing[i], 1,outarray[i], 1);
			}
		}
        }
	void ForcingBody::v_ReadForceInfo(
		LibUtilities::SessionReaderSharedPtr              pSession,
		Array<OneD, MultiRegions::ExpListSharedPtr>       pFields,
                SpatialDomains::MeshGraphSharedPtr                pGraph)
	{
                std::string  m_SolverInfo = pSession->GetSolverInfo("SolverType");
                int nvariables  = pSession->GetVariables().size();
		
                if(m_SolverInfo == "VelocityCorrectionScheme")
                {
                   m_NumVariable = nvariables-1; // e.g. (u v w p) for 3D case
                }
                if(m_SolverInfo == "CoupledLinearisedNS")
                {
                   m_NumVariable = nvariables;  // e.g. (u v w)  for 3D case
                }
		
		if(pSession->DefinesFunction("BodyForce"))
		{		
			m_Forcing     = Array<OneD, Array<OneD, NekDouble> >(m_NumVariable);
			for(int i = 0; i < m_NumVariable; ++i)
                	{						
				m_Forcing[i]  = Array<OneD, NekDouble> (pFields[0]->GetTotPoints(),0.0);
			}
		}
		
	        std::string s_FieldStr;
                for(int i = 0; i < m_NumVariable; ++i)
                {
                      	s_FieldStr = pSession->GetVariable(i);
		     	if(pSession->DefinesFunction("BodyForce"))
		    	{
		   	 	EvaluateFunction(pFields, pSession, s_FieldStr, m_Forcing[i],"BodyForce");
		   	}
                }	
	}/// The end of reading sponge info.
	
    }
}
/// Hui XU  2013 Jul 26
