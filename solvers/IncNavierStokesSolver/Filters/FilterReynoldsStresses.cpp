///////////////////////////////////////////////////////////////////////////////
//
// File FilterReynoldsStresses.cpp
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
// Description: Generate information to be able to calculate the
// Reynolds streses
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/Filters/FilterReynoldsStresses.h>

namespace Nektar
{
    namespace SolverUtils
    {
      std::string FilterReynoldsStresses::className = GetFilterFactory().RegisterCreatorFunction("ReynoldsStresses", FilterReynoldsStresses::create);

        FilterReynoldsStresses::FilterReynoldsStresses(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::map<std::string, std::string> &pParams) :
            Filter(pSession)
        {
            if (pParams.find("OutputFile") == pParams.end())
            {
                m_outputFile = m_session->GetSessionName();
            }
            else
            {
                ASSERTL0(!(pParams.find("OutputFile")->second.empty()),
                         "Missing parameter 'OutputFile'.");
                m_outputFile = pParams.find("OutputFile")->second;
            }

            if(pParams.find("SampleFrequency") == pParams.end())
            {
                m_sampleFrequency = 1;
            }
            else
            {
                m_sampleFrequency = atoi(pParams.find("SampleFrequency")->second.c_str());
            }

            if(pParams.find("OutputFrequency") == pParams.end())
            {
                m_outputFrequency = m_session->GetParameter("NumSteps"); 
            }
            else
            {
                m_outputFrequency = atoi(pParams.find("OutputFrequency")->second.c_str());
            }

            m_numAverages = 0;
            m_index       = 0;
            m_outputIndex = 0;
            m_fld = MemoryManager<LibUtilities::FieldIO>::AllocateSharedPtr(pSession->GetComm());

        }

        FilterReynoldsStresses::~FilterReynoldsStresses()
        {
        }

        void FilterReynoldsStresses::v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {

	  int dim = pFields.num_elements()-1;
	  
	  if(dim == 2)
	  {
	    m_fields = Array<OneD,  Array<OneD, NekDouble> >(3);
	  }
	  else
	  {
	    m_fields = Array<OneD,  Array<OneD, NekDouble> >(6);
	  }

	  for(int n =0; n < m_fields.num_elements(); ++n)
	  {
	    m_fields[n] = Array<OneD, NekDouble>(pFields[0]->GetTotPoints(),0.0);
	  }
            m_avgFieldMetaData["InitialTime"] = boost::lexical_cast<std::string>(time);
        }

        void FilterReynoldsStresses::v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            m_index++;
            if (m_index % m_sampleFrequency > 0)
            {
                return;
            }
	    Array<OneD, NekDouble> tmp(m_fields[0].num_elements());
	    
	    // collect velocities squared. 
	    int dim = pFields.num_elements()-1;
            for(int n = 0; n < pFields.num_elements()-1; ++n)
            {
	      Vmath::Vmul(m_fields[0].num_elements(),pFields[n]->GetPhys(),1,
			  pFields[n]->GetPhys(),1, tmp,1);
	      
	      Vmath::Vadd(m_fields[0].num_elements(),
			  tmp,1,m_fields[n],1,m_fields[n],1);
            }
	    
	    if(dim == 2)
	    {
	      // get uv field
	      Vmath::Vmul(m_fields[0].num_elements(),pFields[0]->GetPhys(),1,
			  pFields[1]->GetPhys(),1, tmp,1);
	      
	      Vmath::Vadd(m_fields[0].num_elements(),
			  tmp,1,m_fields[2],1,m_fields[2],1);
	    }
	    else
	    {
	      // get uv field
	      Vmath::Vmul(m_fields[0].num_elements(),pFields[0]->GetPhys(),1,
			  pFields[1]->GetPhys(),1, tmp,1);
	      
	      Vmath::Vadd(m_fields[0].num_elements(),
			  tmp,1,m_fields[3],1,m_fields[3],1);

	      // get uw field
	      Vmath::Vmul(m_fields[0].num_elements(),pFields[0]->GetPhys(),1,
			  pFields[2]->GetPhys(),1, tmp,1);
	      
	      Vmath::Vadd(m_fields[0].num_elements(),
			  tmp,1,m_fields[4],1,m_fields[4],1);

	      // get vw field
	      Vmath::Vmul(m_fields[0].num_elements(),pFields[1]->GetPhys(),1,
			  pFields[2]->GetPhys(),1, tmp,1);
	      
	      Vmath::Vadd(m_fields[0].num_elements(),
			  tmp,1,m_fields[5],1,m_fields[5],1);
	    }
	    
            m_numAverages += 1;
            // update FinalTime here since this is when last field will be added. 
            m_avgFieldMetaData["FinalTime"] = boost::lexical_cast<std::string>(time);

            if (m_index % m_outputFrequency == 0)
            {
                OutputAvgField(pFields,++m_outputIndex);
            }
            
        }
        
        void FilterReynoldsStresses::v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
        {
            OutputAvgField(pFields);
        }

        void FilterReynoldsStresses::OutputAvgField(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, int dump)
        {
            for(int n = 0; n < m_fields.num_elements(); ++n)
            {
                Vmath::Smul(m_fields[n].num_elements(), 1.0/m_numAverages,m_fields[n],1,m_fields[n],1);
            }

            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                = pFields[0]->GetFieldDefinitions();
            std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
	    Array<OneD, std::string> Variables;
	    
	    if(pFields.num_elements()-1 == 2)
	    {
	      Variables = Array<OneD, std::string>(3);
	      Variables[0] = "uu";
	      Variables[1] = "vv";
	      Variables[2] = "uv";
	    }
	    else
	    {
	      Variables = Array<OneD, std::string>(6);
	      Variables[0] = "uu";
	      Variables[1] = "vv";
	      Variables[2] = "ww";
	      Variables[3] = "uv";
	      Variables[4] = "uw";
	      Variables[5] = "vw";

	    }

            
            Array<OneD, NekDouble>  fieldcoeffs(pFields[0]->GetNcoeffs());



            // copy Data into FieldData and set variable
            for(int j = 0; j < m_fields.num_elements(); ++j)
	    {                
	      // project data to coefficient space
	      pFields[0]->FwdTrans_IterPerExp(m_fields[j],fieldcoeffs);

	      for(int i = 0; i < FieldDef.size(); ++i)
              {
		// Could do a search here to find correct variable
		FieldDef[i]->m_fields.push_back(Variables[j]);
		pFields[0]->AppendFieldData(FieldDef[i], FieldData[i], fieldcoeffs);
	      }
            }

            m_avgFieldMetaData["NumberOfFieldDumps"] = boost::lexical_cast<std::string>(m_numAverages);

            std::stringstream outname; 
            if(dump == -1) // final dump
            {
                outname <<  m_outputFile << "_rey.fld";
            }
            else
            {
                outname << m_outputFile<< "_" << dump << "_rey.fld";
            }
            
            m_fld->Write(outname.str(),FieldDef,FieldData,m_avgFieldMetaData);

            if(dump != -1) // not final dump so rescale cummulative average
            {
                for(int n = 0; n < m_fields.num_elements(); ++n)
                {
                    Vmath::Smul(m_fields[n].num_elements(), 
				(NekDouble) m_numAverages,
                                m_fields[n],1,m_fields[n],1);
                }
            }
        }

        bool FilterReynoldsStresses::v_IsTimeDependent()
        {
            return true;
        }
    }
}
