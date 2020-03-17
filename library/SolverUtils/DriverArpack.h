///////////////////////////////////////////////////////////////////////////////
//
// File DriverArpack.h
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
// Description: Driver class for the stability solver using Arpack
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_DRIVERARPACK_H
#define NEKTAR_SOLVERUTILS_DRIVERARPACK_H

#include <LibUtilities/LinearAlgebra/Arpack.hpp>

#include <SolverUtils/DriverArnoldi.h>

namespace Nektar
{
namespace SolverUtils
{

/// Base class for the development of solvers.
class DriverArpack: public DriverArnoldi
{
public:
    friend class MemoryManager<DriverArpack>;

    /// Creates an instance of this class
    static DriverSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr& pSession,
        const SpatialDomains::MeshGraphSharedPtr& pGraph)
    {
        DriverSharedPtr p = MemoryManager<DriverArpack>::AllocateSharedPtr(
            pSession, pGraph);
        p->InitObject();
        return p;
    }

    ///Name of the class
    static std::string className;



protected:
    int m_maxn;//Maximum size of the problem
    int m_maxnev;//maximum number of eigenvalues requested
    int m_maxncv;//Largest number of basis vector used in Implicitly Restarted Arnoldi

    /// Constructor
    DriverArpack(const LibUtilities::SessionReaderSharedPtr pSession,
                 const SpatialDomains::MeshGraphSharedPtr pGraph);

    /// Destructor
    virtual ~DriverArpack();

    /// Virtual function for initialisation implementation.
    virtual void v_InitObject(std::ostream &out = std::cout);

    /// Virtual function for solve implementation.
    virtual void v_Execute(std::ostream &out = std::cout);

    static std::string arpackProblemTypeLookupIds[];
    static std::string arpackProblemTypeDefault;
    static std::string driverLookupId;

private:
    static std::string ArpackProblemTypeTrans[];
};

}
} //end of namespace

#endif //NEKTAR_SOLVERUTILS_DRIVERARPACK_H

