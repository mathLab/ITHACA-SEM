///////////////////////////////////////////////////////////////////////////////
//
// File AdvectionTerm.h
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
// Description: Driver class for the stability solver
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DRIVERSTANDARD_H
#define NEKTAR_SOLVERUTILS_DRIVERSTANDARD_H

#include <SolverUtils/Driver.h>

namespace Nektar
{
    namespace SolverUtils
    {
        /// Base class for the development of solvers.
        class DriverStandard: public Driver
        {
        public:
            friend class MemoryManager<DriverStandard>;

            /// Creates an instance of this class
            static DriverSharedPtr create(
                const LibUtilities::SessionReaderSharedPtr& pSession,
                const SpatialDomains::MeshGraphSharedPtr& pGraph)
            {
                DriverSharedPtr p = MemoryManager<DriverStandard>
                    ::AllocateSharedPtr(pSession, pGraph);
                p->InitObject();
                return p;
            }
	
            ///Name of the class
            static std::string className;
	
        protected:
            /// Constructor
            SOLVER_UTILS_EXPORT DriverStandard(
                const LibUtilities::SessionReaderSharedPtr pSession,
                const SpatialDomains::MeshGraphSharedPtr pGraph);

            /// Destructor
            SOLVER_UTILS_EXPORT virtual ~DriverStandard();
        
            /// Second-stage initialisation
            SOLVER_UTILS_EXPORT virtual void v_InitObject(std::ostream &out = std::cout);

            /// Virtual function for solve implementation.
            SOLVER_UTILS_EXPORT virtual void v_Execute(std::ostream &out = std::cout);
		
            static std::string driverLookupId;
	};
    }	
} //end of namespace

#endif //NEKTAR_SOLVERUTILS_DRIVERSTANDARD_H

