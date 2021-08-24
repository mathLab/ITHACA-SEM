///////////////////////////////////////////////////////////////////////////////
//
// File Metis.hpp
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
// Description: wrapper of functions around METIS routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASICUTILS_METIS_HPP
#define NEKTAR_LIB_UTILITIES_BASICUTILS_METIS_HPP

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include <metis.h>

namespace Metis
{
    inline static void PartGraphVKway(
            int&                              nVerts,
            int&                              nVertConds,
            Nektar::Array<Nektar::OneD, int>& xadj,
            Nektar::Array<Nektar::OneD, int>& adjcy,
            Nektar::Array<Nektar::OneD, int>& vertWgt,
            Nektar::Array<Nektar::OneD, int>& vertSize,
            Nektar::Array<Nektar::OneD, int>& edgeWgt,
            int&                              nparts,
            int&                              volume,
            Nektar::Array<Nektar::OneD, int>& part)
    {
        int *vwgt   = 0;
        int *vsize  = 0;
        int *adjwgt = 0;
        if (vertWgt.size() > 0)
        {
            vwgt = &vertWgt[0];
        }
        if (vertSize.size() > 0)
        {
            vsize = &vertSize[0];
        }
        if (edgeWgt.size() > 0)
        {
            adjwgt = &edgeWgt[0];
        }
        // number of balancing conditions (size of vertex multi-weight)
        int ncon = nVertConds;
        int options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        METIS_PartGraphKway(&nVerts, &ncon, &xadj[0], &adjcy[0], vwgt, vsize,
                            adjwgt, &nparts, 0, 0, options, &volume, &part[0]);
    }
}
#endif //NEKTAR_LIB_UTILITIES_BASICUTILS_METIS_HPP
