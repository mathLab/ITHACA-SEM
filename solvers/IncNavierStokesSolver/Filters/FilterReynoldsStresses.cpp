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
        std::string FilterReynoldsStresses::className = GetFilterFactory().
            RegisterCreatorFunction("ReynoldsStresses",
                                    FilterReynoldsStresses::create);

        FilterReynoldsStresses::FilterReynoldsStresses(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::map<std::string, std::string> &pParams) :
            FilterAverageFields(pSession, pParams)
        {
        }

        FilterReynoldsStresses::~FilterReynoldsStresses()
        {
        }

        void FilterReynoldsStresses::v_Initialise(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble                                         &time)
        {
            FilterAverageFields::v_Initialise(pFields, time);

            int dim = pFields.num_elements()-1;
            int nExtraFields = dim == 2 ? 3 : 6;
            int origFields = pFields.num_elements();

            m_fields.resize(nExtraFields);
            m_avgFields.resize(m_avgFields.size() + nExtraFields);

            for(int n = 0; n < nExtraFields; ++n)
            {
                m_avgFields[n + origFields] = Array<OneD, NekDouble>(
                    pFields[0]->GetNcoeffs(), 0.0);
                m_fields[n] = Array<OneD, NekDouble>(
                    pFields[0]->GetTotPoints(), 0.0);
            }

            if (dim == 2)
            {
                m_variables.push_back("uu");
                m_variables.push_back("vv");
                m_variables.push_back("uv");
            }
            else if (dim == 3)
            {
                m_variables.push_back("uu");
                m_variables.push_back("vv");
                m_variables.push_back("ww");
                m_variables.push_back("uv");
                m_variables.push_back("vw");
                m_variables.push_back("vw");
            }
            else
            {
                ASSERTL0(false, "Unsupported dimension");
            }
        }

        void FilterReynoldsStresses::v_AddExtraFields(
            const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields,
            const NekDouble                                         &time)
        {
	    Array<OneD, NekDouble> tmp(m_fields[0].num_elements());
	    int dim = pFields.num_elements()-1;
            int extraFields = dim == 2 ? 3 : 6;
	    
	    // collect velocities squared. 
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

            // Forward transform and put into m_avgFields
            int o = pFields.num_elements();
            for (int i = 0; i < extraFields; ++i)
            {
                Array<OneD, NekDouble> tmp2(pFields[0]->GetNcoeffs());
                pFields[0]->FwdTrans_IterPerExp(m_fields[i], tmp2);
                Vmath::Vadd(m_avgFields[i+o].num_elements(),
                            tmp2, 1, m_avgFields[i+o], 1, m_avgFields[i+o], 1);
            }
        }
        
        bool FilterReynoldsStresses::v_IsTimeDependent()
        {
            return true;
        }
    }
}
