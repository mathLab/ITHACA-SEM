/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * ometis.c
 *
 * This file contains the top level routines for the multilevel recursive
 * bisection algorithm PMETIS.
 *
 * Started 7/24/97
 * George
 *
 * $Id: ometis.c 10513 2011-07-07 22:06:03Z karypis $
 *
 */

#include "metislib.h"

void AS_MlevelNestedDissection(ctrl_t *ctrl, graph_t *graph, idx_t *order, 
                               idx_t lastvtx, idx_t *map, idx_t mdswitch)
{
  idx_t i, j, nvtxs, nbnd;
  idx_t *label, *bndind;
  graph_t *lgraph, *rgraph;
  /* WGGao: add the parameters for separators */
  static int length = 0, level = 0, sex = 1;
  int first = 0, second = 0, third = 0;
  /* WGGao End */

  nvtxs = graph->nvtxs;

  MlevelNodeBisectionMultiple(ctrl, graph);

  IFSET(ctrl->dbglvl, METIS_DBG_SEPINFO, 
      printf("Nvtxs: %6"PRIDX", [%6"PRIDX" %6"PRIDX" %6"PRIDX"]\n", 
        graph->nvtxs, graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]));


  /* Order the nodes in the separator */
  nbnd   = graph->nbnd;
  bndind = graph->bndind;
  label  = graph->label;
  for (i=0; i<nbnd; i++) 
    order[label[bndind[i]]] = --lastvtx;

  SplitGraphOrder(ctrl, graph, &lgraph, &rgraph);

  /* WGGao: print out the partitioning structure after the first level
            of partition */  
  level++;
  length++;
  /* printf("  LEVEL %3d\t", level); */
  for ( i = 0; i < nvtxs; i++)
      if (graph->where[i] == 0) ++first;
      else if (graph->where[i] == 1) ++second;
      else ++third;
  /* printf("First %-12dSecond %-12dThird %-12d\n", first, second, third); */
  map[(length-1)*5]     = level;
  map[(length-1)*5 + 1] = sex;
  map[(length-1)*5 + 2] = first;
  map[(length-1)*5 + 3] = second;
  map[(length-1)*5 + 4] = third;
  /* WGGao end */

  /* Free the memory of the top level graph */
  FreeGraph(&graph);

  /* Recurse on lgraph first, as its lastvtx depends on rgraph->nvtxs, which
     will not be defined upon return from MlevelNestedDissection. */
  if (lgraph->nvtxs > mdswitch && lgraph->nedges > 0) {
    /* WGGao : left branch */
    sex = 1;
    /* WGGao End */
    AS_MlevelNestedDissection(ctrl, lgraph, order, lastvtx-rgraph->nvtxs, map, mdswitch);
  } else {
    MMDOrder(ctrl, lgraph, order, lastvtx-rgraph->nvtxs); 
    FreeGraph(&lgraph);
  }
  if (rgraph->nvtxs > mdswitch && rgraph->nedges > 0) {
    /* WGGao : right branch */
    sex = 2;
    /* WGGao End */
    AS_MlevelNestedDissection(ctrl, rgraph, order, lastvtx, map, mdswitch);
  } else {
    MMDOrder(ctrl, rgraph, order, lastvtx); 
    FreeGraph(&rgraph);
  }
  /* WGGao: back to last level */
  level--;
  if (level == 0) length = 0;
  /* WGGao End */
}

/*************************************************************************/
/*! This function is the entry point for the multilevel nested dissection 
    ordering code. At each bisection, a node-separator is computed using
    a node-based refinement approach.

    \param nvtxs is the number of vertices in the graph.
    \param xadj is of length nvtxs+1 marking the start of the adjancy 
           list of each vertex in adjncy.
    \param adjncy stores the adjacency lists of the vertices. The adjnacy
           list of a vertex should not contain the vertex itself.
    \param vwgt is an array of size nvtxs storing the weight of each 
           vertex. If vwgt is NULL, then the vertices are considered 
           to have unit weight.
    \param numflag is either 0 or 1 indicating that the numbering of 
           the vertices starts from 0 or 1, respectively.
    \param options is an array of size METIS_NOPTIONS used to pass 
           various options impacting the of the algorithm. A NULL
           value indicates use of default options.
    \param perm is an array of size nvtxs such that if A and A' are
           the original and permuted matrices, then A'[i] = A[perm[i]].
    \param iperm is an array of size nvtxs such that if A and A' are
           the original and permuted matrices, then A[i] = A'[iperm[i]].
*/
/*************************************************************************/
int AS_METIS_NodeND(idx_t *nvtxs, idx_t *xadj, idx_t *adjncy, idx_t *vwgt,
                    idx_t *options, idx_t *perm, idx_t *iperm, idx_t *map,
                    idx_t *mdswitch) 
{
  int sigrval=0, renumber=0;
  idx_t i, ii, j, l, nnvtxs=0;
  graph_t *graph=NULL;
  ctrl_t *ctrl;
  idx_t *cptr, *cind, *piperm;
  int numflag = 0;

  /* set up malloc cleaning code and signal catchers */
  if (!gk_malloc_init()) 
    return METIS_ERROR_MEMORY;

  gk_sigtrap();

  if ((sigrval = gk_sigcatch()) != 0) 
    goto SIGTHROW;


  /* set up the run time parameters */
  ctrl = SetupCtrl(METIS_OP_OMETIS, options, 1, 3, NULL, NULL);
  if (!ctrl) {
    gk_siguntrap();
    return METIS_ERROR_INPUT;
  }

  /* if required, change the numbering to 0 */
  if (ctrl->numflag == 1) {
    Change2CNumbering(*nvtxs, xadj, adjncy);
    renumber = 1;
  }

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, InitTimers(ctrl));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_startcputimer(ctrl->TotalTmr));

  /* WGGao: force-disable compression */
  ctrl->compress = 0;
  /* WGGao: end */

  /* prune the dense columns */
  if (ctrl->pfactor > 0.0) { 
    piperm = imalloc(*nvtxs, "OMETIS: piperm");

    graph = PruneGraph(ctrl, *nvtxs, xadj, adjncy, vwgt, piperm, ctrl->pfactor);
    if (graph == NULL) {
      /* if there was no prunning, cleanup the pfactor */
      gk_free((void **)&piperm, LTERM);
      ctrl->pfactor = 0.0;
    }
    else {
      nnvtxs = graph->nvtxs;
      ctrl->compress = 0;  /* disable compression if prunning took place */
    }
  }

  /* compress the graph; note that compression only happens if not prunning 
     has taken place. */
  if (ctrl->compress) { 
    cptr = imalloc(*nvtxs+1, "OMETIS: cptr");
    cind = imalloc(*nvtxs, "OMETIS: cind");

    graph = CompressGraph(ctrl, *nvtxs, xadj, adjncy, vwgt, cptr, cind);
    if (graph == NULL) {
      /* if there was no compression, cleanup the compress flag */
      gk_free((void **)&cptr, &cind, LTERM);
      ctrl->compress = 0; 
    }
    else {
      nnvtxs = graph->nvtxs;
      ctrl->cfactor = 1.0*(*nvtxs)/nnvtxs;
      if (ctrl->cfactor > 1.5 && ctrl->nseps == 1)
        ctrl->nseps = 2;
      //ctrl->nseps = (idx_t)(ctrl->cfactor*ctrl->nseps);
    }
  }

  /* if no prunning and no compression, setup the graph in the normal way. */
  if (ctrl->pfactor == 0.0 && ctrl->compress == 0) 
    graph = SetupGraph(ctrl, *nvtxs, 1, xadj, adjncy, vwgt, NULL, NULL);

  ASSERT(CheckGraph(graph, ctrl->numflag, 1));

  /* allocate workspace memory */
  AllocateWorkSpace(ctrl, graph);

  /* do the nested dissection ordering  */
  if (ctrl->ccorder) 
    MlevelNestedDissectionCC(ctrl, graph, iperm, graph->nvtxs);
  else
    AS_MlevelNestedDissection(ctrl, graph, iperm, graph->nvtxs, map, *mdswitch);


  if (ctrl->pfactor > 0.0) { /* Order any prunned vertices */
    icopy(nnvtxs, iperm, perm);  /* Use perm as an auxiliary array */
    for (i=0; i<nnvtxs; i++)
      iperm[piperm[i]] = perm[i];
    for (i=nnvtxs; i<*nvtxs; i++)
      iperm[piperm[i]] = i;

    gk_free((void **)&piperm, LTERM);
  }
  else if (ctrl->compress) { /* Uncompress the ordering */
    /* construct perm from iperm */
    for (i=0; i<nnvtxs; i++)
      perm[iperm[i]] = i; 
    for (l=ii=0; ii<nnvtxs; ii++) {
      i = perm[ii];
      for (j=cptr[i]; j<cptr[i+1]; j++)
        iperm[cind[j]] = l++;
    }

    gk_free((void **)&cptr, &cind, LTERM);
  }

  for (i=0; i<*nvtxs; i++)
    perm[iperm[i]] = i;

  IFSET(ctrl->dbglvl, METIS_DBG_TIME, gk_stopcputimer(ctrl->TotalTmr));
  IFSET(ctrl->dbglvl, METIS_DBG_TIME, PrintTimers(ctrl));

  /* clean up */
  FreeCtrl(&ctrl);

SIGTHROW:
  /* if required, change the numbering back to 1 */
  if (renumber)
    Change2FNumberingOrder(*nvtxs, xadj, adjncy, perm, iperm);

  gk_siguntrap();
  gk_malloc_cleanup(0);

  return metis_rcode(sigrval);
}
