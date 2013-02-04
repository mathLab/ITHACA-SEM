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
// Description: wrapper of functions around METIS routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASICUTILS_METIS_HPP
#define NEKTAR_LIB_UTILITIES_BASICUTILS_METIS_HPP

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

namespace Metis
{
    extern "C"
    {
        // -- Sparse MAtrix Reordering (equivalent to onmetis)
        void   METIS_NodeND (int *nVerts, int *xadj, int *adjncy, int *numflag, int *options, 
                             int *perm, int *iperm);  

        void AS_METIS_NodeND(int *nVerts, int *xadj, int *adjncy, int *numflag, int *options, 
                             int *perm, int *iperm, int *map, int *mdswitch);

        void   METIS_PartMeshNodal(int *nElmts, int *nVerts, int *elType, int *numflag,
                             int *nparts, int *edgecut, int *epart, int *npart);

        void   METIS_PartGraphVKway(int *nVerts, int *xadj, int *adjcy, int *vertWgt, int *vertSize,
                             int *wgtFlag, int *numflag, int *nparts, int *options, int *volume, int *part);
    }

//#ifdef NEKTAR_USING_METIS
    inline static void onmetis(int *nVerts, int *xadj, int *adjncy, int *numflag, int *options, 
                  int *perm, int *iperm)
    {
        METIS_NodeND(nVerts,xadj,adjncy,numflag,options,perm,iperm);
    }

    inline static void as_onmetis(int *nVerts, int *xadj, int *adjncy, int *numflag, int *options, 
                           int *perm, int *iperm, int *map, int *mdswitch)
    {
        AS_METIS_NodeND(nVerts,xadj,adjncy,numflag,options,perm,iperm,map,mdswitch) ;   
    }

    inline static void onmetis(int nVerts, Nektar::Array<Nektar::OneD, int> xadj, Nektar::Array<Nektar::OneD, int> adjncy,
                        Nektar::Array<Nektar::OneD, int> perm,  Nektar::Array<Nektar::OneD, int> iperm)
    {
        ASSERTL1(xadj.num_elements() == nVerts+1,"Array xadj out of bounds");
        ASSERTL1(perm.num_elements() == nVerts,"Array perm out of bounds");
        ASSERTL1(iperm.num_elements() == nVerts,"Array iperm out of bounds");

        int numflag = 0;
        int options[8];
        options[0]=0;
        METIS_NodeND(&nVerts,&xadj[0],&adjncy[0],&numflag,options,&perm[0],&iperm[0]);
    }

    inline static void as_onmetis(int nVerts, Nektar::Array<Nektar::OneD, int> xadj, Nektar::Array<Nektar::OneD, int> adjncy,
                           Nektar::Array<Nektar::OneD, int> perm,  Nektar::Array<Nektar::OneD, int> iperm, Nektar::Array<Nektar::OneD, int> map,
                           int mdswitch)
    {
        ASSERTL1(xadj.num_elements() == nVerts+1,"Array xadj out of bounds");
        ASSERTL1(perm.num_elements() == nVerts,"Array perm out of bounds");
        ASSERTL1(iperm.num_elements() == nVerts,"Array iperm out of bounds");

        int numflag = 0;
        int options[8];
        options[0]=0;
        AS_METIS_NodeND(&nVerts,&xadj[0],&adjncy[0],&numflag,options,&perm[0],&iperm[0],&map[0],&mdswitch);
    }
   
    inline static void MeshPartition(int nElmts, int nVerts, Nektar::Array<Nektar::OneD, int>& mesh, int type, int nparts,
                            Nektar::Array<Nektar::OneD, int>& edgePart, Nektar::Array<Nektar::OneD, int>& nodePart)
    {
        int numflag = 0;
        METIS_PartMeshNodal(&nElmts, &nVerts, &mesh[0], &type, &numflag, &nparts, &edgePart[0], &nodePart[0]);
    }

    inline static void PartGraphVKway( int& nVerts,
                                Nektar::Array<Nektar::OneD, int>& xadj,
                                Nektar::Array<Nektar::OneD, int>& adjcy,
                                Nektar::Array<Nektar::OneD, int>& vertWgt,
                                Nektar::Array<Nektar::OneD, int>& vertSize,
                                int& nparts,
                                int& volume,
                                Nektar::Array<Nektar::OneD, int>& part)
    {
        int wgtflag = 0;
        int *wgts = 0;
        int *sizes = 0;
        if (vertWgt.num_elements() > 0)
        {
            wgtflag += 1;
            wgts = &vertWgt[0];
        }
        if (vertSize.num_elements() > 0)
        {
            wgtflag += 2;
            sizes = &vertSize[0];
        }
        int numflag = 0;
        int options[5];
        options[0]=0;
        METIS_PartGraphVKway(&nVerts, &xadj[0], &adjcy[0], wgts, sizes, &wgtflag,
                             &numflag, &nparts, options, &volume, &part[0]);
    }
//#endif //NEKTAR_USING_METIS
}
#endif //NEKTAR_LIB_UTILITIES_BASICUTILS_METIS_HPP
