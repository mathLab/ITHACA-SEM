#include <boost/algorithm/string.hpp>


#include <IncNavierStokesSolver/EquationSystems/CoupledLinearNS_sparse.h>
//include <IncNavierStokesSolver/EquationSystems/CoupledLinearNS.h>

#include <LibUtilities/BasicUtils/Timer.h>
#include <LocalRegions/MatrixKey.h>
#include <MultiRegions/GlobalLinSysDirectStaticCond.h>
#include <MultiRegions/ContField.h>

using namespace std;

namespace Nektar
{

    string CoupledLinearNS_sparse::className = SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction("CoupledLinearNS_sparse", CoupledLinearNS_sparse::create);


    /**
     *  @class CoupledLinearNS
     *
     * Set up expansion field for velocity and pressure, the local to
     * global mapping arrays and the basic memory definitions for
     * coupled matrix system
     */
    CoupledLinearNS_sparse::CoupledLinearNS_sparse(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
        : UnsteadySystem(pSession, pGraph),
          CoupledLinearNS(pSession, pGraph),
          m_zeroMode(false)
    {
    }














}
