///////////////////////////////////////////////////////////////////////////////
//
// File BoussinesqEquations.h
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
// Description: Boussinesq equations definition in two-dimensions
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_BOUSSINESQEQUATIONS_BOUSSINESQEQUATIONS_H
#define NEKTAR_SOLVERS_BOUSSINESQEQUATIONS_BOUSSINESQEQUATIONS_H

#include <MultiRegions/DisContField2D.h>
#include <Auxiliary/ADRBase.h>
#include <ShallowWaterSolver/ShallowWaterEquations.h>

namespace Nektar
{     
    /**
     * 
     * 
     **/
    
  class BoussinesqEquations: public ShallowWaterEquations 
    {
    public:           

        /**
         * Default constructor. 
	 **/ 
        BoussinesqEquations();


        /**
         * Constructor.
	 **/
        BoussinesqEquations(string &fileStringName);
	
        void ODEforcing(const Array<OneD, const  Array<OneD, NekDouble> >&inarray, 
			Array<OneD, Array<OneD, NekDouble> >&outarray, NekDouble time);
	
        void ExplicitlyIntegrateAdvection(int nsteps);
	
        void Summary(std::ostream &out);

	enum DispersiveFluxType
        {           ///< flux not defined
	  eNotSet,  ///< averaged (or centred) flux
	  eAverage, ///< simple upwind flux
	};

	enum BoussinesqType
	{
	  ePeregrine,
	  eMadsenSorensen,
	  eFullyNonlinear,
	};
	
    protected:

    private: 
	NekDouble m_alpha_0; 
	NekDouble m_alpha_1;

	void NumericalFluxWaveCont(Array<OneD, Array<OneD, NekDouble> > &inarray,
				   Array<OneD, NekDouble> &numfluxX, 
				   Array<OneD, NekDouble> &numfluxY);
	
    };
    
    typedef boost::shared_ptr<BoussinesqEquations> BoussinesqEquationsSharedPtr;
    
} //end of namespace

#endif //NEKTAR_SOLVERS_BOUSSINESQEQUATIONS_BOUSSINESQEQUATIONS_H

