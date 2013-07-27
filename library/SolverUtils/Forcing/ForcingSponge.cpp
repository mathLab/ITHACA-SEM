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

#include <SolverUtils/Forcing/ForcingSponge.h>

namespace Nektar
{
    namespace SolverUtils
    {
	
	std::string ForcingSponge::className = GetForcingFactory().RegisterCreatorFunction("SpongeForcing", ForcingSponge::create, "Forcing Sponge");

	ForcingSponge::ForcingSponge()
	{
	}
        void ForcingSponge::v_InitObject(
                LibUtilities::SessionReaderSharedPtr              pSession,
                Array<OneD, MultiRegions::ExpListSharedPtr>       pFields,
                SpatialDomains::MeshGraphSharedPtr                pGraph)
        {
		m_Session = pSession; 
		v_ReadSpongeInfo(pSession,pFields,pGraph);
        }
	void ForcingSponge::v_Apply(
                const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                      Array<OneD, Array<OneD, NekDouble> >        &outarray)
	{
		if (m_Session->DefinesFunction("RefFields"))
		{
			for (int i=0; i<m_NumVariable;i++)
			{
				Vmath::Vsub(m_Forcing[i].num_elements(),inarray[i],1,m_Refflow[i],1,m_Forcing[i],1);
				Vmath::Vmul(m_Forcing[i].num_elements(), m_Sponge[i], 1, m_Forcing[i], 1, m_Forcing[i], 1);	
				Vmath::Vadd(outarray[i].num_elements(), m_Forcing[i], 1, outarray[i], 1,outarray[i], 1);
			}
		}
		else
		{
			for (int i=0; i<m_NumVariable;i++)
			{
				Vmath::Vmul(m_Forcing[i].num_elements(), m_Sponge[i], 1, inarray[i], 1, m_Forcing[i], 1);
				Vmath::Vadd(outarray[i].num_elements(),  m_Forcing[i], 1, outarray[i], 1,outarray[i], 1);
			}
		}
        }
	void ForcingSponge::EvaluateFunction(
		Array<OneD, MultiRegions::ExpListSharedPtr>       pFields,
		LibUtilities::SessionReaderSharedPtr              pSession,
		std::string 					  pFieldName,
            	Array<OneD, NekDouble>&                           pArray,
            	const std::string&                                pFunctionName,
                NekDouble                                         pTime)
	{
	
		ASSERTL0(pSession->DefinesFunction(pFunctionName),
                     "Function '" + pFunctionName + "' does not exist.");

            	unsigned int nq = pFields[0]->GetNpoints();
		if (pArray.num_elements() != nq)
		{
		      pArray = Array<OneD, NekDouble>(nq);
		}

		LibUtilities::FunctionType vType;
		vType = pSession->GetFunctionType(pFunctionName, pFieldName);
                if (vType == LibUtilities::eFunctionTypeExpression)
                {
		        Array<OneD,NekDouble> x0(nq);
		        Array<OneD,NekDouble> x1(nq);
		        Array<OneD,NekDouble> x2(nq);
		        
		        pFields[0]->GetCoords(x0,x1,x2);
		        LibUtilities::EquationSharedPtr ffunc
		            = pSession->GetFunction(pFunctionName, pFieldName);

		        ffunc->Evaluate(x0,x1,x2,pTime,pArray);
            	}
		else if (vType == LibUtilities::eFunctionTypeFile)
            	{
                	std::string filename
                	    = pSession->GetFunctionFilename(pFunctionName, pFieldName);
		        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
		        std::vector<std::vector<NekDouble> > FieldData;
		        Array<OneD, NekDouble> vCoeffs(pFields[0]->GetNcoeffs());
		        Vmath::Zero(vCoeffs.num_elements(),vCoeffs,1);
		        
		        LibUtilities::Import(filename,FieldDef,FieldData);
		        int idx = -1;
		        for(int i = 0; i < FieldDef.size(); ++i)
		        {
		            for(int j = 0; j < FieldDef[i]->m_fields.size(); ++j)
		            {
		                if (FieldDef[i]->m_fields[j] == pFieldName)
		                {
		                    idx = j;
		                }
		            }

		            if(idx >= 0 )
		            {
		                pFields[0]->ExtractDataToCoeffs(FieldDef[i], FieldData[i],
		                                                 FieldDef[i]->m_fields[idx],
		                                                 vCoeffs);
		            }
		            else
		            {
		                cout << "Field " + pFieldName + " not found." << endl;
		            }
		        }
		        pFields[0]->BwdTrans_IterPerExp(vCoeffs, pArray);
		    }
	}
	void ForcingSponge::v_ReadSpongeInfo(
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
		if(pSession->DefinesFunction("SpongeCoefficient"))
		{		
			m_Sponge      = Array<OneD, Array<OneD, NekDouble> >(m_NumVariable);
			m_Forcing     = Array<OneD, Array<OneD, NekDouble> >(m_NumVariable);
			for(int i = 0; i < m_NumVariable; ++i)
                	{						
	               		m_Sponge[i]   = Array<OneD, NekDouble> (pFields[0]->GetTotPoints(),0.0);
				m_Forcing[i]  = Array<OneD, NekDouble> (pFields[0]->GetTotPoints(),0.0);
			}
		}
		if (pSession->DefinesFunction("RefFields"))
		{
                	m_Refflow     = Array<OneD, Array<OneD, NekDouble> >(m_NumVariable);
			for(int i = 0; i < m_NumVariable; ++i)
                	{						
	               		m_Refflow[i]  = Array<OneD, NekDouble> (pFields[0]->GetTotPoints(),0.0);
			}
		}
                
	        std::string s_FieldStr;
                for(int i = 0; i < m_NumVariable; ++i)
                {
                      	s_FieldStr = pSession->GetVariable(i);
		     	if(pSession->DefinesFunction("SpongeCoefficient"))
		    	{
		   	 	EvaluateFunction(pFields, pSession, s_FieldStr, m_Sponge[i],"SpongeCoefficient");
		   	}
			if(pSession->DefinesFunction("RefFields"))
			{
				EvaluateFunction(pFields, pSession, s_FieldStr, m_Refflow[i],"RefFields");
			}
                }	
	}/// The end of reading sponge info.
	
    }
}
/// Hui XU  2013 Jul 21
