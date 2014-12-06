/*
 * AdvectionSystem.cpp
 *
 *  Created on: 28 May 2014
 *      Author: cc
 */

#include <SolverUtils/AdvectionSystem.h>

namespace Nektar {
namespace SolverUtils {

/**
 *
 */
AdvectionSystem::AdvectionSystem(
        const LibUtilities::SessionReaderSharedPtr& pSession)
    : UnsteadySystem(pSession)
{
}


/**
 *
 */
AdvectionSystem::~AdvectionSystem()
{

}


/**
 *
 */
void AdvectionSystem::v_InitObject()
{
    UnsteadySystem::v_InitObject();
}

}
}
