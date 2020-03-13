///////////////////////////////////////////////////////////////////////////////
//
// File Scotch.hpp
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
// Description: wrapper of functions around Scotch routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASICUTILS_SCOTCH_HPP
#define NEKTAR_LIB_UTILITIES_BASICUTILS_SCOTCH_HPP

namespace Scotch
{
    extern "C" {
        // Scotch version of old metis call
        void   METIS_PartGraphVKway(int *nVerts, int *xadj, int *adjcy,
                int *vertWgt, int *vertSize,
                int *wgtFlag, int *numflag, int *nparts,
                int *options, int *volume, int *part);

    }

    inline static void PartGraphVKway(
            int&                              nVerts,
            int&                              nVertConds,
            Nektar::Array<Nektar::OneD, int>& xadj,
            Nektar::Array<Nektar::OneD, int>& adjcy,
            Nektar::Array<Nektar::OneD, int>& vertWgt,
            Nektar::Array<Nektar::OneD, int>& vertSize,
            int&                              nparts,
            int&                              volume,
            Nektar::Array<Nektar::OneD, int>& part)
    {
        int wgtflag = 0;
        int *vwgt = 0;
        int *vsize = 0;
        if (vertWgt.num_elements() > 0)
        {
            wgtflag += 1;
            vwgt = &vertWgt[0];
        }
        if (vertSize.num_elements() > 0)
        {
            wgtflag += 2;
            vsize = &vertSize[0];
        }
        int numflag = 0;
        // number of balancing conditions (size of vertex multi-weight)
        int options[5];
        options[0] = 0;
        //METIS_SetDefaultOptions(options);
        METIS_PartGraphVKway(&nVerts, &xadj[0], &adjcy[0], vwgt, vsize,
                            &wgtflag, &numflag, &nparts, options, &volume,
                            &part[0]);
    }

}

