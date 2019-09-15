///////////////////////////////////////////////////////////////////////////////
//
// File CommCwipi.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2018 Kilian Lackhove
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
// Description: MPI communication implementation
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Communication/CommCwipi.h>

#include <cwipi.h>

namespace Nektar
{
namespace LibUtilities
{

std::string CommCwipi::className = GetCommFactory().RegisterCreatorFunction(
    "CWIPI", CommCwipi::create, "Parallel communication using MPI with CWIPI.");

/**
 *
 */
CommCwipi::CommCwipi(int narg, char *arg[]) : CommMpi()
{
    int init = 0;
    MPI_Initialized(&init);
    ASSERTL0(!init, "MPI has already been initialised.");

    int retval = MPI_Init(&narg, &arg);
    if (retval != MPI_SUCCESS)
    {
        ASSERTL0(false, "Failed to initialise MPI");
    }

    std::string localName = "";
    for (int i = 0; i < narg; ++i)
    {
        if (!std::strcmp(arg[i], "--cwipi"))
        {
            localName = arg[i + 1];
        }
    }

    MPI_Comm localComm;
    cwipi_init(MPI_COMM_WORLD, localName.c_str(), &localComm);
    m_comm = localComm;

    MPI_Comm_size(m_comm, &m_size);
    MPI_Comm_rank(m_comm, &m_rank);

    m_type = "Parallel MPI with CWIPI";
}

/**
 *
 */
CommCwipi::~CommCwipi()
{
}

/**
 *
 */
void CommCwipi::v_Finalise()
{
    cwipi_finalize();
    CommMpi::v_Finalise();
}
}
}
