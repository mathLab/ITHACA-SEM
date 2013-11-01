///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingBody.cpp
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
// Description: Implementation of a forcing term to simulate a wavy body
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Forcing/ForcingWavyness.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace SolverUtils
{

    std::string ForcingWavyness::className = GetForcingFactory().
                                RegisterCreatorFunction("Wavyness",
                                                        ForcingBody::create,
                                                        "Wavyness Forcing");

    ForcingWavyness::ForcingWavyness(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : Forcing(pSession)
    {
    }

    void ForcingWavyness::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
            const unsigned int& pNumForcingFields,
            const TiXmlElement* pForce)
    {
		// Just 3D homogenous 1D problems can use this techinque
		ASSERTL0(pFields[0]->GetExpType()==MultiRegions::e3DH1D,"Wavyness implemented just for"
                                                                "3D Homogenous 1D expansions.");
		
        // Work in progress
    }

    void ForcingWavyness::v_Apply(const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                    const Array<OneD, Array<OneD, NekDouble> > &inarray,
                                    Array<OneD, Array<OneD, NekDouble> > &outarray,
                                    NekDouble time)
    {
		
		// Work in progress
        
		//calcualte the forcing components Ax,Ay,Az and put them in m_Forcing
		CalculateForcing(fields,inarray);
		
		// Apply forcing terms
        for (int i = 0; i < m_NumVariable; i++)
        {
            Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1,
                        m_Forcing[i], 1, outarray[i], 1);
        }
    }
				 
	void ForcingWavyness::CalculateForcing(const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                                             const Array<OneD, Array<OneD, NekDouble> > &inarray)
	{
        //////////////// Work in progress
		
	}
}
}
