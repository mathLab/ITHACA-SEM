///////////////////////////////////////////////////////////////////////////////
//
// File: StimulusCircle.cpp
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
// Description: Cell model base class.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/VmathArray.hpp>

#include <CardiacEPSolver/Stimuli/StimulusCircle.h>

#include <StdRegions/StdNodalTriExp.h>
//#include <LibUtilities/LinearAlgebra/Blas.hpp>

namespace Nektar
{
    /**
     * @class StimulusCircle
     */
    
    /**
     * Cell model base class constructor.
     */
    StimulusCircle::StimulusCircle(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const MultiRegions::ExpListSharedPtr& pField,
            const TiXmlElement* pXml)
            : Stimulus(pSession, pField, pXml)
    {

    }
    
    
    /**
     * Initialise the cell model. Allocate workspace and variable storage.
     */
    void StimulusCircle::Initialise()
    {

    }
    
    void StimulusCircle::v_Update(Array<OneD, Array<OneD, NekDouble> >&outarray,
                          const NekDouble time)
    {

    }

    void StimulusCircle::v_PrintSummary(std::ostream &out)
    {

    }

    void StimulusCircle::v_SetInitialConditions()
    {

    }

}
