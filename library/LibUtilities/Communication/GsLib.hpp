///////////////////////////////////////////////////////////////////////////////
//
// File GsLib.hpp
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
// Description: wrapper of functions around GSLib routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_COMMUNICATION_GSLIB_HPP
#define NEKTAR_LIB_UTILITIES_COMMUNICATION_GSLIB_HPP

#include <iostream>
using namespace std;

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#ifdef NEKTAR_USE_MPI
#include <LibUtilities/Communication/CommMpi.h>
#endif
using namespace Nektar;

namespace Gs
{
    typedef enum { gs_double, gs_float, gs_int, gs_long, gs_dom_n } gs_dom;
    typedef enum { gs_add, gs_mul, gs_min, gs_max, gs_bpr, gs_op_n } gs_op;
    typedef enum { mode_plain, mode_vec, mode_many, mode_dry_run } gs_mode;

    typedef struct { void *ptr; size_t n,max; } array;
    typedef array buffer;
#ifdef NEKTAR_USE_MPI
    typedef MPI_Comm comm_ext;
    typedef MPI_Request comm_req;
#else
    typedef int comm_ext;
    typedef int comm_req;
#endif

    struct comm {
      unsigned int id;
      unsigned int np;
      comm_ext c;
    };

    typedef struct {
      unsigned int n;      /* number of messages */
      unsigned int *p;     /* message source/dest proc */
      unsigned int *size;  /* size of message */
      unsigned int total;  /* sum of message sizes */
    } pw_comm_data;

    typedef struct {
      pw_comm_data comm[2];
      const unsigned int *map[2];
      comm_req *req;
      unsigned int buffer_size;
    } pw_data;

    typedef struct {
      const unsigned int *scatter_map, *gather_map;
      unsigned int size_r, size_r1, size_r2;
      unsigned int size_sk, size_s, size_total;
      unsigned int p1, p2;
      unsigned int nrecvn;
    } cr_stage_data;

    typedef struct {
      cr_stage_data *stage[2];
      unsigned int nstages;
      unsigned int buffer_size, stage_buffer_size;
    } cr_data;

    typedef struct {
      const unsigned int *map_to_buf[2], *map_from_buf[2];
      unsigned int buffer_size;
    } allreduce_data;

    typedef void exec_fun(
      void *data, gs_mode mode, unsigned vn, gs_dom dom, gs_op op,
      unsigned transpose, const void *execdata, const struct comm *comm, char *buf);

    typedef struct {
      struct comm comm;
      const unsigned int *map_local[2]; /* 0=unflagged, 1=all */
      const unsigned int *flagged_primaries;
      pw_data *pwd;
      cr_data *crd;
      allreduce_data *ard;
      unsigned int buffer_size;
      void *execdata;
      exec_fun *exec;
    } gs_data;

    extern "C"
    {
        void nektar_gs(void *u, gs_dom dom, gs_op op, unsigned transpose,
                gs_data *gsh, buffer *buf);
        gs_data *nektar_gs_setup(const long *id, unsigned int n, const struct comm *comm);
        void nektar_gs_free(gs_data *gsh);
        void nektar_gs_unique(const long *id, unsigned int n, const struct comm *comm);
    }


    /**
     * @brief Initialise Gather-Scatter map.
     *
     * On each process an array of IDs for each global degree of freedom is
     * supplied which corresponds to a unique numbering of universal degrees of
     * freedom. This is used to initialise the GSLib mapping between process-
     * boundary degrees of freedom on different processes.
     * @param   pId         Array of integers providing universal IDs for each
     *                      global DOF on the process.
     * @param   pComm       Communication object used for inter-process
     *                      communication.
     * @return GSLib data structure containing mapping information.
     */
    static inline gs_data* Init (  const Nektar::Array<OneD, long> pId,
                            const LibUtilities::CommSharedPtr& pComm)
    {
#ifdef NEKTAR_USE_MPI
        if (pComm->GetSize() == 1)
        {
            return 0;
        }
        LibUtilities::CommMpiSharedPtr vCommMpi = boost::dynamic_pointer_cast<LibUtilities::CommMpi> (pComm);
        ASSERTL1(vCommMpi, "Failed to cast MPI Comm object.");
        comm vComm;
        MPI_Comm_dup(vCommMpi->GetComm(), &vComm.c);
        vComm.id = vCommMpi->GetRank();
        vComm.np = vCommMpi->GetSize();
        return nektar_gs_setup(&pId[0], pId.num_elements(), &vComm);
#else
        return 0;
#endif
    }


    /**
     * @brief Updates pId to negate all-but-one references to each universal ID.
     *
     * The array of universal IDs corresponding to the process-local DOF are
     * updated such that the ID of only one instance of each universal ID
     * remains positive. This allows the consistent formulation of universally
     * -distributed dot products, for which the contributions of each DOF must
     * be included only once.
     */
    static inline void Unique(const Nektar::Array<OneD, long> pId,
                              const LibUtilities::CommSharedPtr& pComm)
    {
#ifdef NEKTAR_USE_MPI
        if (pComm->GetSize() == 1)
        {
            return;
        }
        LibUtilities::CommMpiSharedPtr vCommMpi = boost::dynamic_pointer_cast<LibUtilities::CommMpi> (pComm);
        ASSERTL1(vCommMpi, "Failed to cast MPI Comm object.");
        comm vComm;
        vComm.c  = vCommMpi->GetComm();
        vComm.id = vCommMpi->GetRank();
        vComm.np = vCommMpi->GetSize();
        nektar_gs_unique(&pId[0], pId.num_elements(), &vComm);
#endif
    }


    /**
     * @brief Deallocates the GSLib mapping data.
     */
    static inline void Finalise (gs_data *pGsh)
    {
#ifdef NEKTAR_USE_MPI
        if (pGsh)
        {
            nektar_gs_free(pGsh);
        }
#endif
    }


    /**
     * @brief Performs a gather-scatter operation of the provided values.
     *
     * The
     */
    static inline void Gather(Nektar::Array<OneD, NekDouble> pU, gs_op pOp,
                       gs_data *pGsh, Nektar::Array<OneD, NekDouble> pBuffer
                                                        = NullNekDouble1DArray)
    {
#ifdef NEKTAR_USE_MPI
        if (!pGsh)
        {
            return;
        }
        if (pBuffer.num_elements() == 0)
        {
            nektar_gs(&pU[0], gs_double, pOp, false, pGsh, 0);
        }
        else
        {
            array buf;
            buf.ptr = &pBuffer[0];
            buf.n = pBuffer.num_elements();
            nektar_gs(&pU[0], gs_double, pOp, false, pGsh, &buf);
        }
#endif
    }

}

#endif
