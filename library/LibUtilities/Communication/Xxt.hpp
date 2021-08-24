///////////////////////////////////////////////////////////////////////////////
//
// File Xxt.hpp
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
// Description: wrapper of functions around XXT routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_COMMUNICATION_XXT_HPP
#define NEKTAR_LIB_UTILITIES_COMMUNICATION_XXT_HPP

#include <iostream>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/Communication/Comm.h>
#ifdef NEKTAR_USE_MPI
#include <LibUtilities/Communication/CommMpi.h>
#endif

namespace Xxt
{
using namespace Nektar;
#ifdef NEKTAR_USE_MPI
typedef MPI_Comm comm_ext;
typedef MPI_Request comm_req;
#else
typedef int comm_ext;
typedef int comm_req;
#endif

struct comm
{
    unsigned int id;
    unsigned int np;
    comm_ext c;
};

struct sparse_cholesky
{
    unsigned int n, *Lrp, *Lj;
    double *L, *D;
};

struct csr_mat
{
    unsigned int n, *Arp, *Aj;
    double *A;
};

struct crs_data
{

    /* communication */
    struct comm comm;
    unsigned int pcoord; /* coordinate in communication tree */
    unsigned plevels;    /* # of stages of communication */
    signed int *pother;  /* let p = pother[i], then during stage i of fan-in,
                            if p>=0, receive from p
                            if p< 0, send to (-p-1)
                    fan-out is just the reverse ...
                    on proc 0, pother is never negative
                    on others, pother is negative for the last stage only */
    comm_req *req;

    /* separators */
    unsigned nsep;          /* number of separators */
    unsigned int *sep_size; /* # of dofs on each separator,
                 ordered from the bottom to the top of the tree:
                 separator 0 is the bottom-most one (dofs not shared)
                 separator nsep-1 is the root of the tree */

    unsigned null_space;
    double *share_weight;

    /* vector sizes */
    unsigned int un; /* user's vector size */

    /* xxt_solve works with "condensed" vectors;
       same dofs as in user's vectors, but no duplicates and no Dirichlet
       nodes, and also ordered topologically (children before parents)
       according to the separator tree */

    unsigned int cn;      /* size of condensed vectors */
    signed int *perm_u2c; /* permutation from user vector to condensed
                        vector, p=perm_u2c[i]; xu[i] = p=-1 ? 0 : xc[p];*/
    unsigned int ln, sn;  /* xc[0 ... ln-1] are not shared
                            (ln=sep_size[0])
               xc[ln ... ln+sn-1] are shared
               ln+sn = cn                    */

    unsigned int xn; /* # of columns of x = sum_i(sep_size[i]) -
                         sep_size[0] */

    /* data */
    struct sparse_cholesky fac_A_ll;
    struct csr_mat A_sl;
    unsigned int *Xp;
    double *X; /* column i of X starts at X[Xp[i]] */

    /* execution buffers */
    double *vl, *vc, *vx, *combuf;
};

extern "C" {
struct crs_data *nektar_crs_setup(unsigned int n, const unsigned long *id,
                                  unsigned int nz, const unsigned int *Ai,
                                  const unsigned int *Aj, const double *A,
                                  unsigned int null_space,
                                  const struct comm *comm);
void nektar_crs_solve(double *x, struct crs_data *data, double *b);
void nektar_crs_stats(struct crs_data *data);
void nektar_crs_free(struct crs_data *data);
}

/**
 * @brief Initialise the matrix-solve.
 *
 * On each process an array of IDs for each global degree of freedom is
 * supplied which corresponds to a unique numbering of universal degrees of
 * freedom. Three vectors describing the matrix are also provided. The
 * parallel matrix solve is then set up.
 *
 * @param   pId         Array of integers providing universal IDs for each
 *                      global DOF on the process.
 * @param   pAi         Row indices of matrix entries
 * @param   pAj         Column indices of matrix entries
 * @param   pAr         Values of matrix entries
 * @param   pComm       Communication object used for inter-process
 *                      communication.
 * @returns crs_data structure
 */
static inline struct crs_data *Init(
    unsigned int pRank, const Nektar::Array<OneD, unsigned long> pId,
    const Nektar::Array<OneD, unsigned int> pAi,
    const Nektar::Array<OneD, unsigned int> pAj,
    const Nektar::Array<OneD, NekDouble> pAr,
    const LibUtilities::CommSharedPtr &pComm)
{
#ifdef NEKTAR_USE_MPI
    unsigned int nz = pAr.size();
    LibUtilities::CommMpiSharedPtr vCommMpi =
        std::dynamic_pointer_cast<LibUtilities::CommMpi>(pComm);
    ASSERTL1(vCommMpi, "Failed to cast MPI Comm object.");
    comm vComm;
    MPI_Comm_dup(vCommMpi->GetComm(), &vComm.c);
    vComm.id         = vCommMpi->GetRank();
    vComm.np         = vCommMpi->GetSize();
    crs_data *result = nektar_crs_setup(pRank, &pId[0], nz, &pAi[0], &pAj[0],
                                        &pAr[0], 0, &vComm);
    MPI_Comm_free(&vComm.c);
    return result;
#else
    return 0;
#endif
}

/**
 * @brief Solve the matrix system for a given input vector b.
 */
static inline void Solve(Nektar::Array<OneD, NekDouble> pX,
                         struct crs_data *pCrs,
                         Nektar::Array<OneD, NekDouble> pB)
{
#ifdef NEKTAR_USE_MPI
    if (!pCrs)
    {
        return;
    }
    nektar_crs_solve(&pX[0], pCrs, &pB[0]);
#endif
}

/**
 * @brief Deallocates the crs mapping data.
 */
static inline void Finalise(crs_data *pCrs)
{
#ifdef NEKTAR_USE_MPI
    int finalized;
    MPI_Finalized(&finalized);
    if (pCrs && !finalized)
    {
        nektar_crs_free(pCrs);
    }
#endif
}
}

#endif
