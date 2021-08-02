///////////////////////////////////////////////////////////////////////////////
//
// File CoupledLinearNS.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Coupled Stokes solver scheme header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_COUPLEDSTOKESSPARSE_H
#define NEKTAR_SOLVERS_COUPLEDSTOKESSPARSE_H

//#include <IncNavierStokesSolver/EquationSystems/CoupledLocalToGlobalC0ContMap.h>
//#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <IncNavierStokesSolver/EquationSystems/CoupledLinearNS.h>
//#include <MultiRegions/GlobalLinSys.h>
//#include <MultiRegions/ExpList3DHomogeneous1D.h>
//#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
//#include <MultiRegions/GlobalLinSysDirectStaticCond.h>

namespace Nektar
{     
    
    
    
    class CoupledLinearNS_sparse: public CoupledLinearNS
    {
	public:
        friend class MemoryManager<CoupledLinearNS_sparse>;
        
        /// Creates an instance of this class
        static SolverUtils::EquationSystemSharedPtr create(
            const LibUtilities::SessionReaderSharedPtr& pSession,
            const SpatialDomains::MeshGraphSharedPtr &pGraph)
        {
            SolverUtils::EquationSystemSharedPtr p = MemoryManager<CoupledLinearNS_sparse>::AllocateSharedPtr(pSession, pGraph);
            p->InitObject();
            return p;
        }

        /// Name of class
        static std::string className;

    protected:
        CoupledLinearNS_sparse(const LibUtilities::SessionReaderSharedPtr &pSesssion,
                        const SpatialDomains::MeshGraphSharedPtr &pGraph);

    private:
        /// Id to identify when single mode is mean mode (i.e. beta=0);
        bool m_zeroMode;



	};
    

    
} //end of namespace

#endif //COUPLEDSTOKESSCHEME_H
