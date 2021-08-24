///////////////////////////////////////////////////////////////////////////////
//
// File MeshPartitionScotch.cpp
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
// Description: Scotch partitioner interface
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SpatialDomains/MeshPartitionScotch.h>

namespace Nektar
{
namespace SpatialDomains
{

    std::string MeshPartitionScotch::className
        = GetMeshPartitionFactory().RegisterCreatorFunction(
            "Scotch",
            MeshPartitionScotch::create,
            "Partitioning using the Scotch library.");

    std::string MeshPartitionScotch::cmdSwitch
        = LibUtilities::SessionReader::RegisterCmdLineFlag(
            "use-scotch","","Use Scotch for mesh partitioning.");

    MeshPartitionScotch::MeshPartitionScotch(
        const LibUtilities::SessionReaderSharedPtr session,
        int                                        meshDim,
        std::map<int, MeshEntity>                  element,
        CompositeDescriptor                        compMap)
        : MeshPartition(session, meshDim, element, compMap)
    {
    }

    MeshPartitionScotch::~MeshPartitionScotch()
    {

    }

    void MeshPartitionScotch::PartitionGraphImpl(
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
        boost::ignore_unused(nVertConds, edgeWgt);

        int wgtflag = 0;
        int *vwgt = 0;
        int *vsize = 0;
        if (vertWgt.size() > 0)
        {
            wgtflag += 1;
            vwgt = &vertWgt[0];
        }
        if (vertSize.size() > 0)
        {
            wgtflag += 2;
            vsize = &vertSize[0];
        }
        int numflag = 0;

        PartGraphVKway(&nVerts, &xadj[0], &adjcy[0], vwgt, vsize,
                            &wgtflag, &numflag, &nparts, &volume,
                            &part[0]);
    }


    void MeshPartitionScotch::PartGraphVKway (
            const SCOTCH_Num * const    n,
            const SCOTCH_Num * const    xadj,
            const SCOTCH_Num * const    adjncy,
            const SCOTCH_Num * const    vwgt,
            const SCOTCH_Num * const    vsize,
            const SCOTCH_Num * const    wgtflag,
            const SCOTCH_Num * const    numflag,
            const SCOTCH_Num * const    nparts,
            SCOTCH_Num * const          volume,
            SCOTCH_Num * const          part)
    {
        SCOTCH_Num          baseval;
        const SCOTCH_Num *  vwgt2;
        const SCOTCH_Num *  vsize2;
        SCOTCH_Num          vsizval; /// Communication volume of current vertex
        SCOTCH_Num          vertnbr;
        SCOTCH_Num          vertnum;
        SCOTCH_Num          edgenum;
        const SCOTCH_Num *  edgetax;
        const SCOTCH_Num *  parttax;
        SCOTCH_Num *        nghbtab;
        SCOTCH_Num                  commvol;

        vsize2  = ((*wgtflag & 1) != 0) ? vsize : NULL;
        vwgt2   = ((*wgtflag & 2) != 0) ? vwgt  : NULL;
        baseval = *numflag;
        vertnbr = *n;
        edgetax = adjncy - baseval;

        // If no communication load   data provided
        if (vsize2 == NULL) {
            if (PartGraph2 (n, xadj, adjncy, vwgt2, NULL, numflag, nparts,
                            part, SCOTCH_STRATQUALITY, 0.01) != 0)
                return;
        }

        // Will have to turn communication volumes into edge loads
        else {
            const SCOTCH_Num *  vsiztax;
            SCOTCH_Num          edgenbr;
            SCOTCH_Num *        edlotax;
            int                 o;

            edgenbr = xadj[vertnbr] - baseval;
            if ((edlotax =
                 (SCOTCH_Num*) malloc (edgenbr * sizeof (SCOTCH_Num))) == NULL)
            {
                return;
            }

            edlotax -= baseval; // Base access to edlotax
            vsiztax  = vsize2 - baseval;

            // Un-based scan of vertex array xadj
            for (vertnum = 0, edgenum = baseval;
                    vertnum < vertnbr; vertnum ++) {
                SCOTCH_Num      vsizval; // Communication size of current vertex
                SCOTCH_Num      edgennd;

                vsizval = vsize2[vertnum];
                // Based traversal of edge array adjncy
                for (edgennd = xadj[vertnum + 1];
                     edgenum < edgennd;
                     edgenum ++) {

                    SCOTCH_Num  vertend; // Based end vertex number

                    vertend = edgetax[edgenum];
                    edlotax[edgenum] = vsizval + vsiztax[vertend];
                }
            }

            o = PartGraph2 (n, xadj, adjncy, vwgt2, edlotax + baseval, numflag,
                            nparts, part, SCOTCH_STRATQUALITY, 0.01);

            free (edlotax + baseval);

            if (o != 0)
                return;
        }

        if ((nghbtab =
             (SCOTCH_Num*) malloc (*nparts * sizeof (SCOTCH_Num))) == NULL)
        {
            return;
        }

        memset (nghbtab, ~0, *nparts * sizeof (SCOTCH_Num));

        parttax = part - baseval;
        vsizval = 1; // Assume no vertex communication sizes

        // Un-based scan of vertex array xadj
        for (vertnum = 0, edgenum = baseval, commvol = 0;
                vertnum < vertnbr; vertnum ++) {
            SCOTCH_Num      partval;
            SCOTCH_Num      edgennd;

            partval = part[vertnum];
            nghbtab[partval] = vertnum; // Do not count local neighbors in
                                        // communication volume
            if (vsize2 != NULL)
                vsizval = vsize2[vertnum];

            // Based traversal of edge array adjncy
            for (edgennd = xadj[vertnum + 1]; edgenum < edgennd; edgenum ++) {
                SCOTCH_Num  vertend; // Based end vertex number
                SCOTCH_Num  partend;

                vertend = edgetax[edgenum];
                partend = parttax[vertend];

                // If first neighbor in this part set part as accounted for
                if (nghbtab[partend] != vertnum) {
                    nghbtab[partend] = vertnum;
                    commvol += vsizval;
                }
            }
        }
        *volume = commvol;

        free (nghbtab);
    }

    int MeshPartitionScotch::PartGraph2 (
            const SCOTCH_Num * const    n,
            const SCOTCH_Num * const    xadj,
            const SCOTCH_Num * const    adjncy,
            const SCOTCH_Num * const    vwgt,
            const SCOTCH_Num * const    adjwgt,
            const SCOTCH_Num * const    numflag,
            const SCOTCH_Num * const    nparts,
            SCOTCH_Num * const          part,
            SCOTCH_Num                  flagval,
            double                      kbalval)
    {
        // Scotch graph object to interface with libScotch
        SCOTCH_Graph *      grafdat = SCOTCH_graphAlloc();
        SCOTCH_Strat        stradat;
        SCOTCH_Num          baseval;
        SCOTCH_Num          vertnbr;
        int                 o;

        SCOTCH_graphInit (grafdat);

        baseval = *numflag;
        vertnbr = *n;

        o = 1; // Assume something will go wrong
        if (SCOTCH_graphBuild (grafdat, baseval, vertnbr, xadj, xadj + 1,
                               vwgt, NULL, xadj[vertnbr] - baseval, adjncy,
                               adjwgt) == 0) {
            SCOTCH_stratInit          (&stradat);
            SCOTCH_stratGraphMapBuild (&stradat, flagval, *nparts, kbalval);
#ifdef SCOTCH_DEBUG_ALL
            // TRICK: next instruction called only if graph is consistent
            if (SCOTCH_graphCheck (grafdat) == 0)
#endif /* SCOTCH_DEBUG_ALL */
                o = SCOTCH_graphPart (grafdat, *nparts, &stradat, part);
            SCOTCH_stratExit (&stradat);
        }
        SCOTCH_graphExit (grafdat);

        if (o != 0)
            return (1);

        if (baseval != 0) { // MeTiS part array is based, Scotch is not
            SCOTCH_Num          vertnum;

            for (vertnum = 0; vertnum < vertnbr; vertnum ++)
                part[vertnum] += baseval;
        }

        return (0);
    }

}
}
