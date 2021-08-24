///////////////////////////////////////////////////////////////////////////////
//
// File Driver.h
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
// Description: Base class for Drivers.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef SOLVERUTILS_DRIVER_H
#define SOLVERUTILS_DRIVER_H

#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <SolverUtils/SolverUtils.hpp>
#include <SolverUtils/EquationSystem.h>

namespace Nektar
{
namespace SolverUtils
{

class Driver;

/// A shared pointer to a Driver object
typedef std::shared_ptr<Driver> DriverSharedPtr;

/// Datatype of the NekFactory used to instantiate classes derived from
/// the Driver class.
typedef LibUtilities::NekFactory<
    std::string, Driver,
    const LibUtilities::SessionReaderSharedPtr &,
    const SpatialDomains::MeshGraphSharedPtr &> DriverFactory;

SOLVER_UTILS_EXPORT DriverFactory& GetDriverFactory();

/// Base class for the development of solvers.
class Driver
{
public:
    /// Destructor
    virtual ~Driver();

    /// Initialise Object
    SOLVER_UTILS_EXPORT inline void InitObject(std::ostream &out = std::cout);

    /// Execute driver
    SOLVER_UTILS_EXPORT inline void Execute(std::ostream &out = std::cout);

    SOLVER_UTILS_EXPORT inline Array<OneD, EquationSystemSharedPtr> GetEqu();

    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> GetRealEvl(void);
    SOLVER_UTILS_EXPORT Array<OneD, NekDouble> GetImagEvl(void);


protected:
    /// Communication object
    LibUtilities::CommSharedPtr                 m_comm;

    /// Session reader object
    LibUtilities::SessionReaderSharedPtr        m_session;

    /// I the Coupling between SFD and arnoldi
    LibUtilities::SessionReaderSharedPtr        session_LinNS;

    /// MeshGraph object
    SpatialDomains::MeshGraphSharedPtr          m_graph;

    /// Equation system to solve
    Array<OneD, EquationSystemSharedPtr>        m_equ;

    ///number of equations
    int m_nequ;

    ///Evolution Operator
    enum EvolutionOperatorType m_EvolutionOperator;

    /// Initialises EquationSystem class members.
    Driver(const LibUtilities::SessionReaderSharedPtr pSession,
           const SpatialDomains::MeshGraphSharedPtr   pGraph);

    SOLVER_UTILS_EXPORT virtual void v_InitObject(std::ostream &out = std::cout);

    /// Virtual function for solve implementation.
    SOLVER_UTILS_EXPORT virtual void v_Execute(std::ostream &out = std::cout) = 0;


    SOLVER_UTILS_EXPORT virtual Array<OneD, NekDouble> v_GetRealEvl(void);
    SOLVER_UTILS_EXPORT virtual Array<OneD, NekDouble> v_GetImagEvl(void);


    static std::string evolutionOperatorLookupIds[];
    static std::string evolutionOperatorDef;
    static std::string driverDefault;

};

inline void Driver::InitObject(std::ostream &out)
{
    v_InitObject(out);
}

inline void Driver::Execute(std::ostream &out)
{
    v_Execute(out);
}

inline Array<OneD, EquationSystemSharedPtr> Driver::GetEqu()
{
    return m_equ;
}

inline Array<OneD, NekDouble> Driver::GetRealEvl()
{
    return v_GetRealEvl();
}

inline Array<OneD, NekDouble> Driver::GetImagEvl()
{
    return v_GetImagEvl();
}

}
} //end of namespace

#endif //NEKTAR_SOLVERS_AUXILIARY_ADRBASE_H

