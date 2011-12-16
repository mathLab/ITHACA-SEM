///////////////////////////////////////////////////////////////////////////////
//
// File InterpCoeff.cpp
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
// Description: Definition of Interpolation methods for Coefficients 
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/Foundations/InterpCoeff.h>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Basis.h>

namespace Nektar
{    
    namespace LibUtilities
    {     

      void InterpCoeff1D(const BasisKey &fbasis0, 
			 const Array<OneD, const NekDouble>& from,  
			 const BasisKey &tbasis0, 
			 Array<OneD, NekDouble> &to)
      { 
	ASSERTL0(fbasis0.GetNumModes() == tbasis0.GetNumModes(), "Number of modes must be the same for interpolating coefficients");
	if(fbasis0.GetBasisType() == tbasis0.GetBasisType()) //check to see if the same base
	  {
	    Vmath::Vcopy(fbasis0.GetNumModes(),from,1,to,1);
	  }
	else
	  {
	    // interpolate
	    
	    DNekMatSharedPtr ftB = BasisManager()[fbasis0]->GetI(tbasis0);
	    
	    NekVector<NekDouble> in(fbasis0.GetNumModes(),from,eWrapper);
	    NekVector<NekDouble> out(tbasis0.GetNumModes(),to,eWrapper);
	    
	    out  = (*ftB)*in;
	  }
      }
    }
}
