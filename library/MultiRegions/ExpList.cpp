///////////////////////////////////////////////////////////////////////////////
//
// File ExpList.cpp
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
// Description: Expansion list definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/GlobalLinSys.h>

#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdTriExp.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdTetExp.h>
#include <StdRegions/StdPyrExp.h>
#include <StdRegions/StdPrismExp.h>
#include <StdRegions/StdHexExp.h>

#include <LocalRegions/MatrixKey.h>     // for MatrixKey
#include <LocalRegions/Expansion.h>     // for Expansion

#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>  // for AssemblyMapCG, etc
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>  // for AssemblyMapCG, etc
#include <MultiRegions/GlobalLinSysKey.h>  // for GlobalLinSysKey
#include <MultiRegions/GlobalMatrix.h>  // for GlobalMatrix, etc
#include <MultiRegions/GlobalMatrixKey.h>  // for GlobalMatrixKey

#include <LibUtilities/LinearAlgebra/SparseMatrixFwd.hpp>
#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

#include <Collections/CollectionOptimisation.h>
#include <Collections/Operator.h>

namespace Nektar
{
    namespace MultiRegions
    {
        /**
         * @class ExpList
         * All multi-elemental expansions \f$u^{\delta}(\boldsymbol{x})\f$ can
         * be considered as the assembly of the various elemental contributions.
         * On a discrete level, this yields,
         * \f[u^{\delta}(\boldsymbol{x}_i)=\sum_{e=1}^{{N_{\mathrm{el}}}}
         * \sum_{n=0}^{N^{e}_m-1}\hat{u}_n^e\phi_n^e(\boldsymbol{x}_i).\f]
         * where \f${N_{\mathrm{el}}}\f$ is the number of elements and
         * \f$N^{e}_m\f$ is the local elemental number of expansion modes.
         * As it is the lowest level class, it contains the definition of the
         * common data and common routines to all multi-elemental expansions.
         *
         * The class stores a vector of expansions, \a m_exp, (each derived from
         * StdRegions#StdExpansion) which define the constituent components of
         * the domain. The coefficients from these expansions are concatenated
         * in \a m_coeffs, while the expansion evaluated at the quadrature
         * points is stored in \a m_phys.
         */

        /**
         * Creates an empty expansion list. The expansion list will typically be
         * populated by a derived class (namely one of MultiRegions#ExpList1D,
         * MultiRegions#ExpList2D or MultiRegions#ExpList3D).
         */
        ExpList::ExpList():
            m_comm(),
            m_session(),
            m_graph(),
            m_ncoeffs(0),
            m_npoints(0),
            m_coeffs(),
            m_phys(),
            m_physState(false),
            m_exp(MemoryManager<LocalRegions::ExpansionVector>
                      ::AllocateSharedPtr()),
            m_coeff_offset(),
            m_phys_offset(),
            m_offset_elmt_id(),
            m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
            m_WaveSpace(false)
        {
            SetExpType(eNoType);
        }


        /**
         * Creates an empty expansion list. The expansion list will typically be
         * populated by a derived class (namely one of MultiRegions#ExpList1D,
         * MultiRegions#ExpList2D or MultiRegions#ExpList3D).
         */
        ExpList::ExpList(
                const LibUtilities::SessionReaderSharedPtr &pSession):
            m_comm(pSession->GetComm()),
            m_session(pSession),
            m_graph(),
            m_ncoeffs(0),
            m_npoints(0),
            m_coeffs(),
            m_phys(),
            m_physState(false),
            m_exp(MemoryManager<LocalRegions::ExpansionVector>
                      ::AllocateSharedPtr()),
            m_coeff_offset(),
            m_phys_offset(),
            m_offset_elmt_id(),
            m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
            m_WaveSpace(false)
        {
            SetExpType(eNoType);
        }


        /**
         * Creates an empty expansion list. The expansion list will typically be
         * populated by a derived class (namely one of MultiRegions#ExpList1D,
         * MultiRegions#ExpList2D or MultiRegions#ExpList3D).
         */
        ExpList::ExpList(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &pGraph):
            m_comm(pSession->GetComm()),
            m_session(pSession),
            m_graph(pGraph),
            m_ncoeffs(0),
            m_npoints(0),
            m_coeffs(),
            m_phys(),
            m_physState(false),
            m_exp(MemoryManager<LocalRegions::ExpansionVector>
                      ::AllocateSharedPtr()),
            m_coeff_offset(),
            m_phys_offset(),
            m_offset_elmt_id(),
            m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
            m_WaveSpace(false)
        {
            SetExpType(eNoType);
        }
        
        /**
         * Copies the eIds elements from an existing expansion list.
         * @param   in              Source expansion list.
         * @param   in              elements that will be in the new exp list.
         */
        ExpList::ExpList(const ExpList &in, 
                         const std::vector<unsigned int> &eIDs,
                         const bool DeclareCoeffPhysArrays):
            m_comm(in.m_comm),
            m_session(in.m_session),
            m_graph(in.m_graph),
            m_ncoeffs(0),
            m_npoints(0),
            m_coeffs(),
            m_phys(),
            m_physState(false),
            m_exp(MemoryManager<LocalRegions::ExpansionVector>
                      ::AllocateSharedPtr()),
            m_coeff_offset(),
            m_phys_offset(),
            m_offset_elmt_id(),
            m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
            m_WaveSpace(false)
        {
            SetExpType(eNoType);
            
            for (int i=0; i < eIDs.size(); ++i)
            {
                (*m_exp).push_back( (*(in.m_exp))[eIDs[i]]);
                m_ncoeffs += (*m_exp)[i]->GetNcoeffs();
                m_npoints += (*m_exp)[i]->GetTotPoints();
            }

            if(DeclareCoeffPhysArrays)
            {
                m_coeffs = Array<OneD, NekDouble>(m_ncoeffs, 0.0);
                m_phys   = Array<OneD, NekDouble>(m_npoints, 0.0);
            }
        }


        /**
         * Copies an existing expansion list.
         * @param   in              Source expansion list.
         */
        ExpList::ExpList(const ExpList &in, const bool DeclareCoeffPhysArrays):
            m_comm(in.m_comm),
            m_session(in.m_session),
            m_graph(in.m_graph),
            m_ncoeffs(in.m_ncoeffs),
            m_npoints(in.m_npoints),
            m_physState(false),
            m_exp(in.m_exp),
            m_collections(in.m_collections),
            m_coll_coeff_offset(in.m_coll_coeff_offset),
            m_coll_phys_offset(in.m_coll_phys_offset),
            m_coeff_offset(in.m_coeff_offset),
            m_phys_offset(in.m_phys_offset),
            m_offset_elmt_id(in.m_offset_elmt_id),
            m_globalOptParam(in.m_globalOptParam),
            m_blockMat(in.m_blockMat),
            m_WaveSpace(false)
        {
            SetExpType(eNoType);

            if(DeclareCoeffPhysArrays)
            {
                m_coeffs = Array<OneD, NekDouble>(m_ncoeffs, 0.0);
                m_phys   = Array<OneD, NekDouble>(m_npoints, 0.0);
            }
        }

        /**
         *
         */
        ExpansionType ExpList::GetExpType(void)
        {
            return m_expType;
        }

        /**
         *
         */
        void ExpList::SetExpType(ExpansionType Type)
        {
            m_expType = Type;
        }

        ExpList::~ExpList()
        {
        }

        /**
         * The integration is evaluated locally, that is
         * \f[\int
         *    f(\boldsymbol{x})d\boldsymbol{x}=\sum_{e=1}^{{N_{\mathrm{el}}}}
         * \left\{\int_{\Omega_e}f(\boldsymbol{x})d\boldsymbol{x}\right\},  \f]
         * where the integration over the separate elements is done by the
         * function StdRegions#StdExpansion#Integral, which discretely
         * evaluates the integral using Gaussian quadrature.
         *
         * Note that the array #m_phys should be filled with the values of the
         * function \f$f(\boldsymbol{x})\f$ at the quadrature points
         * \f$\boldsymbol{x}_i\f$.
         *
         * @return  The value of the discretely evaluated integral
         *          \f$\int f(\boldsymbol{x})d\boldsymbol{x}\f$.
         */
        NekDouble ExpList::PhysIntegral()
        {
            ASSERTL1(m_physState == true,
                     "local physical space is not true ");

            return PhysIntegral(m_phys);
        }


        /**
         * The integration is evaluated locally, that is
         * \f[\int
         *    f(\boldsymbol{x})d\boldsymbol{x}=\sum_{e=1}^{{N_{\mathrm{el}}}}
         * \left\{\int_{\Omega_e}f(\boldsymbol{x})d\boldsymbol{x}\right\},  \f]
         * where the integration over the separate elements is done by the
         * function StdRegions#StdExpansion#Integral, which discretely
         * evaluates the integral using Gaussian quadrature.
         *
         * @param   inarray         An array of size \f$Q_{\mathrm{tot}}\f$
         *                          containing the values of the function
         *                          \f$f(\boldsymbol{x})\f$ at the quadrature
         *                          points \f$\boldsymbol{x}_i\f$.
         * @return  The value of the discretely evaluated integral
         *          \f$\int f(\boldsymbol{x})d\boldsymbol{x}\f$.
         */
        NekDouble ExpList::PhysIntegral(
                                const Array<OneD, const NekDouble> &inarray)
        {
            int       i;
            NekDouble sum = 0.0;

            for(i = 0; i < (*m_exp).size(); ++i)
            {
                sum += (*m_exp)[i]->Integral(inarray + m_phys_offset[i]);
            }

            return sum;
        }


        /**
         * Retrieves the block matrix specified by \a bkey, and computes
         * \f$ y=Mx \f$.
         * @param   gkey        GlobalMatrixKey specifying the block matrix to
         *                      use in the matrix-vector multiply.
         * @param   inarray     Input vector \f$ x \f$.
         * @param   outarray    Output vector \f$ y \f$.
         */
        void ExpList::MultiplyByBlockMatrix(
                                const GlobalMatrixKey             &gkey,
                                const Array<OneD,const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray)
        {
            // Retrieve the block matrix using the given key.
            const DNekScalBlkMatSharedPtr& blockmat = GetBlockMatrix(gkey);
            int nrows = blockmat->GetRows();
            int ncols = blockmat->GetColumns();

            // Create NekVectors from the given data arrays
            NekVector<NekDouble> in (ncols,inarray, eWrapper);
            NekVector<      NekDouble> out(nrows,outarray,eWrapper);

            // Perform matrix-vector multiply.
            out = (*blockmat)*in;
        }


        /**
         * The operation is evaluated locally for every element by the function
         * StdRegions#StdExpansion#IProductWRTBase.
         *
         * @param   inarray         An array of size \f$Q_{\mathrm{tot}}\f$
         *                          containing the values of the function
         *                          \f$f(\boldsymbol{x})\f$ at the quadrature
         *                          points \f$\boldsymbol{x}_i\f$.
         * @param   outarray        An array of size \f$N_{\mathrm{eof}}\f$
         *                          used to store the result.
         */
        void ExpList::v_IProductWRTBase_IterPerExp(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray)
        {
            Array<OneD,NekDouble>  tmp;
            for (int i = 0; i < m_collections.size(); ++i)
            {

                m_collections[i].ApplyOperator(Collections::eIProductWRTBase,
                                               inarray + m_coll_phys_offset[i],
                                               tmp = outarray + m_coll_coeff_offset[i]);
            }
        }

        /**
         * The operation is evaluated locally for every element by the function
         * StdRegions#StdExpansion#IProductWRTDerivBase.
         *
         * @param   dir             {0,1} is the direction in which the
         *                          derivative of the basis should be taken
         * @param   inarray         An array of size \f$Q_{\mathrm{tot}}\f$
         *                          containing the values of the function
         *                          \f$f(\boldsymbol{x})\f$ at the quadrature
         *                          points \f$\boldsymbol{x}_i\f$.
         * @param   outarray        An array of size \f$N_{\mathrm{eof}}\f$
         *                          used to store the result.
         */
        void ExpList::IProductWRTDerivBase(const int dir,
                                           const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD, NekDouble> &outarray)
        {
            int    i;

            Array<OneD,NekDouble> e_outarray;

            for(i = 0; i < (*m_exp).size(); ++i)
            {
                (*m_exp)[i]->IProductWRTDerivBase(dir,inarray+m_phys_offset[i],
                                                  e_outarray = outarray+m_coeff_offset[i]);
            }
        }


        /**
         * The operation is evaluated locally for every element by the function
         * StdRegions#StdExpansion#IProductWRTDerivBase.
         *
         * @param   inarray         An array of arrays of size \f$Q_{\mathrm{tot}}\f$
         *                          containing the values of the function
         *                          \f$f(\boldsymbol{x})\f$ at the quadrature
         *                          points \f$\boldsymbol{x}_i\f$ in dir directions.
         * @param   outarray        An array of size \f$N_{\mathrm{eof}}\f$
         *                          used to store the result.
         */
        void ExpList::IProductWRTDerivBase(const Array<OneD, const Array<OneD, NekDouble> > &inarray,
                                           Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, NekDouble> tmp0,tmp1,tmp2;
            // assume coord dimension defines the size of Deriv Base
            int dim = GetCoordim(0);

            ASSERTL1(inarray.num_elements() >= dim,"inarray is not of sufficient dimension");

            switch(dim)
            {
            case 1:
                for (int i = 0; i < m_collections.size(); ++i)
                {
                    m_collections[i].ApplyOperator(
                                                   Collections::eIProductWRTDerivBase,
                                                   inarray[0] + m_coll_phys_offset[i],
                                                   tmp0 = outarray + m_coll_coeff_offset[i]);
                }
                break;
            case 2:
                for (int i = 0; i < m_collections.size(); ++i)
                {
                    m_collections[i].ApplyOperator(
                                                   Collections::eIProductWRTDerivBase,
                                                   inarray[0] + m_coll_phys_offset[i],
                                                   tmp0 = inarray[1] + m_coll_phys_offset[i],
                                                   tmp1 = outarray + m_coll_coeff_offset[i]);
                }
                break;
            case 3:
                for (int i = 0; i < m_collections.size(); ++i)
                {
                    m_collections[i].ApplyOperator(
                                                   Collections::eIProductWRTDerivBase,
                                                   inarray[0] + m_coll_phys_offset[i],
                                                   tmp0 = inarray[1] + m_coll_phys_offset[i],
                                                   tmp1 = inarray[2] + m_coll_phys_offset[i],
                                                   tmp2 = outarray + m_coll_coeff_offset[i]);
                }
                break;
            default:
                ASSERTL0(false,"Dimension of inarray not correct");
                break;
            }
        }
        /**
         * Given a function \f$f(\boldsymbol{x})\f$ evaluated at
         * the quadrature points, this function calculates the
         * derivatives \f$\frac{d}{dx_1}\f$, \f$\frac{d}{dx_2}\f$
         * and \f$\frac{d}{dx_3}\f$ of the function
         * \f$f(\boldsymbol{x})\f$ at the same quadrature
         * points. The local distribution of the quadrature points
         * allows an elemental evaluation of the derivative. This
         * is done by a call to the function
         * StdRegions#StdExpansion#PhysDeriv.
         *
         * @param   inarray         An array of size \f$Q_{\mathrm{tot}}\f$
         *                          containing the values of the function
         *                          \f$f(\boldsymbol{x})\f$ at the quadrature
         *                          points \f$\boldsymbol{x}_i\f$.
         * @param   out_d0          The discrete evaluation of the
         *                          derivative\f$\frac{d}{dx_1}\f$ will
         *                          be stored in this array of size
         *                          \f$Q_{\mathrm{tot}}\f$.
         * @param   out_d1          The discrete evaluation of the
         *                          derivative\f$\frac{d}{dx_2}\f$ will be
         *                          stored in this array of size
         *                          \f$Q_{\mathrm{tot}}\f$. Note that if no
         *                          memory is allocated for \a out_d1, the
         *                          derivative \f$\frac{d}{dx_2}\f$ will not be
         *                          calculated.
         * @param   out_d2          The discrete evaluation of the
         *                          derivative\f$\frac{d}{dx_3}\f$ will be
         *                          stored in this array of size
         *                          \f$Q_{\mathrm{tot}}\f$. Note that if no
         *                          memory is allocated for \a out_d2, the
         *                          derivative \f$\frac{d}{dx_3}\f$ will not be
         *                          calculated.
         */
        void ExpList::v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD, NekDouble> &out_d0,
                                  Array<OneD, NekDouble> &out_d1,
                                  Array<OneD, NekDouble> &out_d2)
        {
            Array<OneD, NekDouble> e_out_d0;
            Array<OneD, NekDouble> e_out_d1;
            Array<OneD, NekDouble> e_out_d2;
            for (int i = 0; i < m_collections.size(); ++i)
            {
                int offset = m_coll_phys_offset[i];
                e_out_d0 = out_d0  + offset;
                e_out_d1 = out_d1  + offset;
                e_out_d2 = out_d2  + offset;

                m_collections[i].ApplyOperator(Collections::ePhysDeriv,
                                               inarray + offset,
                                               e_out_d0,e_out_d1, e_out_d2);

            }
        }

        void ExpList::v_PhysDeriv(const int dir,
                                  const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD, NekDouble> &out_d)
        {
            Direction edir = DirCartesianMap[dir];
            v_PhysDeriv(edir, inarray,out_d);
        }

        void ExpList::v_PhysDeriv(Direction edir, const Array<OneD, const NekDouble> &inarray,
                Array<OneD, NekDouble> &out_d)
        {
            int i;
            if(edir==MultiRegions::eS)
            {
                Array<OneD, NekDouble> e_out_ds;
                for(i=0; i<(*m_exp).size(); ++i)
                {
                    e_out_ds = out_d + m_phys_offset[i];
                    (*m_exp)[i]->PhysDeriv_s(inarray+m_phys_offset[i],e_out_ds);
                }
            }
            else if(edir==MultiRegions::eN)
            {
                Array<OneD, NekDouble > e_out_dn;
                for(i=0; i<(*m_exp).size(); i++)
                {
                    e_out_dn = out_d +m_phys_offset[i];
                    (*m_exp)[i]->PhysDeriv_n(inarray+m_phys_offset[i],e_out_dn);
                }
            }
            else
            {
                // convert enum into int
                int intdir= (int)edir;
                Array<OneD, NekDouble> e_out_d;
                for(i= 0; i < (*m_exp).size(); ++i)
                {
                    e_out_d = out_d + m_phys_offset[i];
                    (*m_exp)[i]->PhysDeriv(intdir, inarray+m_phys_offset[i], e_out_d);
                }

            }
        }


        /**
         * The coefficients of the function to be acted upon
         * should be contained in the \param inarray. The
         * resulting coefficients are stored in \param outarray
         *
         * @param   inarray         An array of size \f$N_{\mathrm{eof}}\f$
         *                          containing the inner product.
         */
        void ExpList::MultiplyByElmtInvMass(
                                            const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD, NekDouble> &outarray)
        {
            GlobalMatrixKey mkey(StdRegions::eInvMass);
            const DNekScalBlkMatSharedPtr& InvMass = GetBlockMatrix(mkey);

            // Inverse mass matrix
            NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);
            if(inarray.get() == outarray.get())
            {
                NekVector<NekDouble> in(m_ncoeffs,inarray); // copy data
                out = (*InvMass)*in;
            }
            else
            {
                NekVector<NekDouble> in(m_ncoeffs,inarray,eWrapper);
                out = (*InvMass)*in;
            }
        }

        /**
         * Given a function \f$u(\boldsymbol{x})\f$ defined at the
         * quadrature points, this function determines the
         * transformed elemental coefficients \f$\hat{u}_n^e\f$
         * employing a discrete elemental Galerkin projection from
         * physical space to coefficient space. For each element,
         * the operation is evaluated locally by the function
         * StdRegions#StdExpansion#IproductWRTBase followed by a
         * call to #MultiRegions#MultiplyByElmtInvMass.
         *
         * @param   inarray         An array of size \f$Q_{\mathrm{tot}}\f$
         *                          containing the values of the function
         *                          \f$f(\boldsymbol{x})\f$ at the quadrature
         *                          points \f$\boldsymbol{x}_i\f$.
         * @param   outarray        The resulting coefficients
         *                          \f$\hat{u}_n^e\f$ will be stored in this
         *                          array of size \f$N_{\mathrm{eof}}\f$.
         */
        void ExpList::v_FwdTrans_IterPerExp(const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD, NekDouble> &outarray)
        {
            Array<OneD,NekDouble> f(m_ncoeffs);

            IProductWRTBase_IterPerExp(inarray,f);
            MultiplyByElmtInvMass(f,outarray);

        }

        void ExpList::FwdTrans_BndConstrained(
                                              const Array<OneD, const NekDouble>& inarray,
                                              Array<OneD, NekDouble> &outarray)
        {
            int i;

            Array<OneD,NekDouble> e_outarray;

            for(i= 0; i < (*m_exp).size(); ++i)
            {
                (*m_exp)[i]->FwdTrans_BndConstrained(inarray+m_phys_offset[i],
                                                     e_outarray = outarray+m_coeff_offset[i]);
            }
        }

        /**
         * This function smooth a field after some calculaitons which have
         * been done elementally.
         *
         * @param   field     An array containing the field in physical space
         *
         */
        void ExpList::v_SmoothField(Array<OneD, NekDouble> &field)
        {
            // Do nothing unless the method is implemented in the appropriate
            // class, i.e. ContField1D,ContField2D, etc.

            // So far it has been implemented just for ContField2D and
            // ContField3DHomogeneous1D

            // Block in case users try the smoothing with a modal expansion.
            // Maybe a different techique for the smoothing require
            // implementation for modal basis.

            ASSERTL0((*m_exp)[0]->GetBasisType(0)
                     == LibUtilities::eGLL_Lagrange ||
                     (*m_exp)[0]->GetBasisType(0)
                     == LibUtilities::eGauss_Lagrange,
                     "Smoothing is currently not allowed unless you are using "
                     "a nodal base for efficiency reasons. The implemented "
                     "smoothing technique requires the mass matrix inversion "
                     "which is trivial just for GLL_LAGRANGE_SEM and "
                     "GAUSS_LAGRANGE_SEMexpansions.");
        }


        /**
         * This function assembles the block diagonal matrix
         * \f$\underline{\boldsymbol{M}}^e\f$, which is the
         * concatenation of the local matrices
         * \f$\boldsymbol{M}^e\f$ of the type \a mtype, that is
         *
         * \f[
         * \underline{\boldsymbol{M}}^e = \left[
         * \begin{array}{cccc}
         * \boldsymbol{M}^1 & 0 & \hspace{3mm}0 \hspace{3mm}& 0 \\
         *  0 & \boldsymbol{M}^2 & 0 & 0 \\
         *  0 &  0 & \ddots &  0 \\
         *  0 &  0 & 0 & \boldsymbol{M}^{N_{\mathrm{el}}} \end{array}\right].\f]
         *
         * @param   mtype           the type of matrix to be assembled
         * @param   scalar          an optional parameter
         * @param   constant        an optional parameter
         */
        const DNekScalBlkMatSharedPtr ExpList::GenBlockMatrix(
                                                              const GlobalMatrixKey &gkey)
        {
            int i,cnt1;
            int n_exp = 0;
            DNekScalMatSharedPtr    loc_mat;
            DNekScalBlkMatSharedPtr BlkMatrix;
            map<int,int> elmt_id;
            LibUtilities::ShapeType ShapeType = gkey.GetShapeType();

            if(ShapeType != LibUtilities::eNoShapeType)
            {
                for(i = 0 ; i < (*m_exp).size(); ++i)
                {
                    if((*m_exp)[m_offset_elmt_id[i]]->DetShapeType()
                       == ShapeType)
                    {
                        elmt_id[n_exp++] = m_offset_elmt_id[i];
                    }
                }
            }
            else
            {
                n_exp = (*m_exp).size();
                for(i = 0; i < n_exp; ++i)
                {
                    elmt_id[i] = m_offset_elmt_id[i];
                }
            }

            Array<OneD,unsigned int> nrows(n_exp);
            Array<OneD,unsigned int> ncols(n_exp);

            switch(gkey.GetMatrixType())
            {
            case StdRegions::eBwdTrans:
                {
                    // set up an array of integers for block matrix construction
                    for(i = 0; i < n_exp; ++i)
                    {
                        nrows[i] = (*m_exp)[elmt_id.find(i)->second]->GetTotPoints();
                        ncols[i] = (*m_exp)[elmt_id.find(i)->second]->GetNcoeffs();
                    }
                }
                break;
            case StdRegions::eIProductWRTBase:
                {
                    // set up an array of integers for block matrix construction
                    for(i = 0; i < n_exp; ++i)
                    {
                        nrows[i] = (*m_exp)[elmt_id.find(i)->second]->GetNcoeffs();
                        ncols[i] = (*m_exp)[elmt_id.find(i)->second]->GetTotPoints();
                    }
                }
                break;
            case StdRegions::eMass:
            case StdRegions::eInvMass:
            case StdRegions::eHelmholtz:
            case StdRegions::eLaplacian:
            case StdRegions::eInvHybridDGHelmholtz:
                {
                    // set up an array of integers for block matrix construction
                    for(i = 0; i < n_exp; ++i)
                    {
                        nrows[i] = (*m_exp)[elmt_id.find(i)->second]->GetNcoeffs();
                        ncols[i] = (*m_exp)[elmt_id.find(i)->second]->GetNcoeffs();
                    }
                }
                break;

            case StdRegions::eHybridDGLamToU:
                {
                    // set up an array of integers for block matrix construction
                    for(i = 0; i < n_exp; ++i)
                    {
                        nrows[i] = (*m_exp)[elmt_id.find(i)->second]->GetNcoeffs();
                        ncols[i] = (*m_exp)[elmt_id.find(i)->second]->NumDGBndryCoeffs();
                    }
                }
                break;

            case StdRegions::eHybridDGHelmBndLam:
                {
                    // set up an array of integers for block matrix construction
                    for(i = 0; i < n_exp; ++i)
                    {
                        nrows[i] = (*m_exp)[elmt_id.find(i)->second]->NumDGBndryCoeffs();
                        ncols[i] = (*m_exp)[elmt_id.find(i)->second]->NumDGBndryCoeffs();
                    }
                }
                break;

            default:
                {
                    NEKERROR(ErrorUtil::efatal,
                             "Global Matrix creation not defined for this type "
                             "of matrix");
                }
            }

            MatrixStorage blkmatStorage = eDIAGONAL;
            BlkMatrix = MemoryManager<DNekScalBlkMat>
                ::AllocateSharedPtr(nrows,ncols,blkmatStorage);

            int nvarcoeffs = gkey.GetNVarCoeffs();
            int eid;
            Array<OneD, NekDouble> varcoeffs_wk;

            for(i = cnt1 = 0; i < n_exp; ++i)
            {
                // need to be initialised with zero size for non variable coefficient case
                StdRegions::VarCoeffMap varcoeffs;

                eid = elmt_id[i];
                if(nvarcoeffs>0)
                {
                    StdRegions::VarCoeffMap::const_iterator x;
                    for (x = gkey.GetVarCoeffs().begin(); x != gkey.GetVarCoeffs().end(); ++x)
                    {
                        varcoeffs[x->first] = x->second + m_phys_offset[eid];
                    }
                }

                LocalRegions::MatrixKey matkey(gkey.GetMatrixType(),
                                               (*m_exp)[eid]->DetShapeType(),
                                               *(*m_exp)[eid],
                                               gkey.GetConstFactors(),
                                               varcoeffs );

                loc_mat = boost::dynamic_pointer_cast<LocalRegions::Expansion>((*m_exp)[elmt_id.find(i)->second])->GetLocMatrix(matkey);
                BlkMatrix->SetBlock(i,i,loc_mat);
            }

            return BlkMatrix;
        }

        const DNekScalBlkMatSharedPtr& ExpList::GetBlockMatrix(
                                                               const GlobalMatrixKey &gkey)
        {
            BlockMatrixMap::iterator matrixIter = m_blockMat->find(gkey);

            if(matrixIter == m_blockMat->end())
            {
                return ((*m_blockMat)[gkey] = GenBlockMatrix(gkey));
            }
            else
            {
                return matrixIter->second;
            }
        }

        void ExpList::GeneralMatrixOp_IterPerExp(
                                                 const GlobalMatrixKey             &gkey,
                                                 const Array<OneD,const NekDouble> &inarray,
                                                 Array<OneD,      NekDouble> &outarray)
        {
            const Array<OneD, const bool>  doBlockMatOp
                        = m_globalOptParam->DoBlockMatOp(gkey.GetMatrixType());
            const Array<OneD, const int> num_elmts
                        = m_globalOptParam->GetShapeNumElements();

            Array<OneD,NekDouble> tmp_outarray;
            int cnt = 0;
            int eid;
            for(int n = 0; n < num_elmts.num_elements(); ++n)
            {
                if(doBlockMatOp[n])
                {
                    const LibUtilities::ShapeType vType
                                    = m_globalOptParam->GetShapeList()[n];
                    const MultiRegions::GlobalMatrixKey vKey(gkey, vType);
                    if (cnt < m_offset_elmt_id.num_elements())
                    {
                        eid = m_offset_elmt_id[cnt];
                        MultiplyByBlockMatrix(vKey,inarray + m_coeff_offset[eid],
                                              tmp_outarray = outarray + m_coeff_offset[eid]);
                        cnt += num_elmts[n];
                    }
                }
                else
                {
                    int i;
                    int nvarcoeffs = gkey.GetNVarCoeffs();

                    for(i= 0; i < num_elmts[n]; ++i)
                    {
                        // need to be initialised with zero size for non variable coefficient case
                        StdRegions::VarCoeffMap varcoeffs;

                        eid = m_offset_elmt_id[cnt++];
                        if(nvarcoeffs>0)
                        {
                            StdRegions::VarCoeffMap::const_iterator x;
                            for (x = gkey.GetVarCoeffs().begin(); x != gkey.GetVarCoeffs().end(); ++x)
                            {
                                varcoeffs[x->first] = x->second + m_phys_offset[eid];
                            }
                        }

                        StdRegions::StdMatrixKey mkey(gkey.GetMatrixType(),
                                                      (*m_exp)[eid]->DetShapeType(),
                                                      *((*m_exp)[eid]),
                                                      gkey.GetConstFactors(),varcoeffs);

                        (*m_exp)[eid]->GeneralMatrixOp(inarray + m_coeff_offset[eid],
                                                       tmp_outarray = outarray+m_coeff_offset[eid],
                                                       mkey);
                    }
                }
            }
        }

        /**
         * Retrieves local matrices from each expansion in the expansion list
         * and combines them together to generate a global matrix system.
         * @param   mkey        Matrix key for the matrix to be generated.
         * @param   locToGloMap Local to global mapping.
         * @returns Shared pointer to the generated global matrix.
         */
        GlobalMatrixSharedPtr ExpList::GenGlobalMatrix(
                                                       const GlobalMatrixKey &mkey,
                                                       const AssemblyMapCGSharedPtr &locToGloMap)
        {
            int i,j,n,gid1,gid2,cntdim1,cntdim2;
            NekDouble sign1,sign2;
            DNekScalMatSharedPtr loc_mat;

            unsigned int glob_rows;
            unsigned int glob_cols;
            unsigned int loc_rows;
            unsigned int loc_cols;

            bool assembleFirstDim;
            bool assembleSecondDim;

            switch(mkey.GetMatrixType())
            {
            case StdRegions::eBwdTrans:
                {
                    glob_rows = m_npoints;
                    glob_cols = locToGloMap->GetNumGlobalCoeffs();

                    assembleFirstDim  = false;
                    assembleSecondDim = true;
                }
                break;
            case StdRegions::eIProductWRTBase:
                {
                    glob_rows = locToGloMap->GetNumGlobalCoeffs();
                    glob_cols = m_npoints;

                    assembleFirstDim  = true;
                    assembleSecondDim = false;
                }
                break;
            case StdRegions::eMass:
            case StdRegions::eHelmholtz:
            case StdRegions::eLaplacian:
            case StdRegions::eHybridDGHelmBndLam:
                {
                    glob_rows = locToGloMap->GetNumGlobalCoeffs();
                    glob_cols = locToGloMap->GetNumGlobalCoeffs();

                    assembleFirstDim  = true;
                    assembleSecondDim = true;
                }
                break;
            default:
                {
                    NEKERROR(ErrorUtil::efatal,
                             "Global Matrix creation not defined for this type "
                             "of matrix");
                }
            }

            COOMatType spcoomat;
            CoordType  coord;

            int nvarcoeffs = mkey.GetNVarCoeffs();
            int eid;

            // fill global matrix
            for(n = cntdim1 = cntdim2 = 0; n < (*m_exp).size(); ++n)
            {
                // need to be initialised with zero size for non variable coefficient case
                StdRegions::VarCoeffMap varcoeffs;

                eid = m_offset_elmt_id[n];
                if(nvarcoeffs>0)
                {
                    StdRegions::VarCoeffMap::const_iterator x;
                    for (x = mkey.GetVarCoeffs().begin(); x != mkey.GetVarCoeffs().end(); ++x)
                    {
                        varcoeffs[x->first] = x->second + m_phys_offset[eid];
                    }
                }

                LocalRegions::MatrixKey matkey(mkey.GetMatrixType(),
                                              (*m_exp)[eid]->DetShapeType(),
                                              *((*m_exp)[eid]),
                                              mkey.GetConstFactors(),varcoeffs);

                loc_mat = boost::dynamic_pointer_cast<LocalRegions::Expansion>((*m_exp)[m_offset_elmt_id[n]])->GetLocMatrix(matkey);

                loc_rows = loc_mat->GetRows();
                loc_cols = loc_mat->GetColumns();

                for(i = 0; i < loc_rows; ++i)
                {
                    if(assembleFirstDim)
                    {
                        gid1  = locToGloMap->GetLocalToGlobalMap (cntdim1 + i);
                        sign1 = locToGloMap->GetLocalToGlobalSign(cntdim1 + i);
                    }
                    else
                    {
                        gid1  = cntdim1 + i;
                        sign1 = 1.0;
                    }

                    for(j = 0; j < loc_cols; ++j)
                    {
                        if(assembleSecondDim)
                        {
                            gid2  = locToGloMap
                                ->GetLocalToGlobalMap(cntdim2 + j);
                            sign2 = locToGloMap
                                ->GetLocalToGlobalSign(cntdim2 + j);
                        }
                        else
                        {
                            gid2  = cntdim2 + j;
                            sign2 = 1.0;
                        }

                        // sparse matrix fill
                        coord = make_pair(gid1,gid2);
                        if( spcoomat.count(coord) == 0 )
                        {
                            spcoomat[coord] = sign1*sign2*(*loc_mat)(i,j);
                        }
                        else
                        {
                            spcoomat[coord] += sign1*sign2*(*loc_mat)(i,j);
                        }
                    }
                }
                cntdim1 += loc_rows;
                cntdim2 += loc_cols;
            }

            return MemoryManager<GlobalMatrix>
                ::AllocateSharedPtr(m_session,glob_rows,glob_cols,spcoomat);
        }


        DNekMatSharedPtr ExpList::GenGlobalMatrixFull(const GlobalLinSysKey &mkey, const AssemblyMapCGSharedPtr &locToGloMap)
        {
            int i,j,n,gid1,gid2,loc_lda,eid;
            NekDouble sign1,sign2,value;
            DNekScalMatSharedPtr loc_mat;

            int totDofs     = locToGloMap->GetNumGlobalCoeffs();
            int NumDirBCs   = locToGloMap->GetNumGlobalDirBndCoeffs();

            unsigned int rows = totDofs - NumDirBCs;
            unsigned int cols = totDofs - NumDirBCs;
            NekDouble zero = 0.0;

            DNekMatSharedPtr Gmat;
            int bwidth = locToGloMap->GetFullSystemBandWidth();

            int nvarcoeffs = mkey.GetNVarCoeffs();
            MatrixStorage matStorage;

            map<int, RobinBCInfoSharedPtr> RobinBCInfo = GetRobinBCInfo();

            switch(mkey.GetMatrixType())
            {
                // case for all symmetric matices
            case StdRegions::eHelmholtz:
            case StdRegions::eLaplacian:
                if( (2*(bwidth+1)) < rows)
                {
                    matStorage = ePOSITIVE_DEFINITE_SYMMETRIC_BANDED;
                    Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero,matStorage,bwidth,bwidth);
                }
                else
                {
                    matStorage = ePOSITIVE_DEFINITE_SYMMETRIC;
                    Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero,matStorage);
                }

                break;
            default: // Assume general matrix - currently only set up for full invert
                {
                    matStorage = eFULL;
                    Gmat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,zero,matStorage);
                }
            }

            // fill global symmetric matrix
            for(n = 0; n < (*m_exp).size(); ++n)
            {
                // need to be initialised with zero size for non variable coefficient case
                StdRegions::VarCoeffMap varcoeffs;

                eid = m_offset_elmt_id[n];
                if(nvarcoeffs>0)
                {
                    StdRegions::VarCoeffMap::const_iterator x;
                    for (x = mkey.GetVarCoeffs().begin(); x != mkey.GetVarCoeffs().end(); ++x)
                    {
                        varcoeffs[x->first] = x->second + m_phys_offset[eid];
                    }
                }

                LocalRegions::MatrixKey matkey(mkey.GetMatrixType(),
                                              (*m_exp)[eid]->DetShapeType(),
                                              *((*m_exp)[eid]),
                                              mkey.GetConstFactors(),varcoeffs);

                loc_mat = boost::dynamic_pointer_cast<LocalRegions::Expansion>((*m_exp)[n])->GetLocMatrix(matkey);


                if(RobinBCInfo.count(n) != 0) // add robin mass matrix
                {
                    RobinBCInfoSharedPtr rBC;

                    // declare local matrix from scaled matrix.
                    int rows = loc_mat->GetRows();
                    int cols = loc_mat->GetColumns();
                    const NekDouble *dat = loc_mat->GetRawPtr();
                    DNekMatSharedPtr new_mat = MemoryManager<DNekMat>::AllocateSharedPtr(rows,cols,dat);
                    Blas::Dscal(rows*cols,loc_mat->Scale(),new_mat->GetRawPtr(),1);

                    // add local matrix contribution
                    for(rBC = RobinBCInfo.find(n)->second;rBC; rBC = rBC->next)
                    {
                        (*m_exp)[n]->AddRobinMassMatrix(rBC->m_robinID,rBC->m_robinPrimitiveCoeffs,new_mat);
                    }

                    NekDouble one = 1.0;
                    // redeclare loc_mat to point to new_mat plus the scalar.
                    loc_mat = MemoryManager<DNekScalMat>::AllocateSharedPtr(one,new_mat);
                }

                loc_lda = loc_mat->GetColumns();

                for(i = 0; i < loc_lda; ++i)
                {
                    gid1 = locToGloMap->GetLocalToGlobalMap(m_coeff_offset[n] + i) - NumDirBCs;
                    sign1 =  locToGloMap->GetLocalToGlobalSign(m_coeff_offset[n] + i);
                    if(gid1 >= 0)
                    {
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2 = locToGloMap->GetLocalToGlobalMap(m_coeff_offset[n] + j) - NumDirBCs;
                            sign2 = locToGloMap->GetLocalToGlobalSign(m_coeff_offset[n] + j);
                            if(gid2 >= 0)
                            {
                                // When global matrix is symmetric,
                                // only add the value for the upper
                                // triangular part in order to avoid
                                // entries to be entered twice
                                if((matStorage == eFULL)||(gid2 >= gid1))
                                {
                                    value = Gmat->GetValue(gid1,gid2) + sign1*sign2*(*loc_mat)(i,j);
                                    Gmat->SetValue(gid1,gid2,value);
                                }
                            }
                        }
                    }
                }
            }

            return Gmat;
        }


        /**
         * Consider a linear system
         * \f$\boldsymbol{M\hat{u}}_g=\boldsymbol{f}\f$ to be solved. Dependent
         * on the solution method, this function constructs
         * - <b>The full linear system</b><BR>
         *   A call to the function #GenGlobalLinSysFullDirect
         * - <b>The statically condensed linear system</b><BR>
         *   A call to the function #GenGlobalLinSysStaticCond
         *
         * @param   mkey            A key which uniquely defines the global
         *                          matrix to be constructed.
         * @param   locToGloMap     Contains the mapping array and required
         *                          information for the transformation from
         *                          local to global degrees of freedom.
         * @return  (A shared pointer to) the global linear system in
         *          required format.
         */
        GlobalLinSysSharedPtr ExpList::GenGlobalLinSys(
                    const GlobalLinSysKey &mkey,
                    const AssemblyMapCGSharedPtr &locToGloMap)
        {
            GlobalLinSysSharedPtr returnlinsys;
            boost::shared_ptr<ExpList> vExpList = GetSharedThisPtr();

            MultiRegions::GlobalSysSolnType vType = mkey.GetGlobalSysSolnType();

            if (vType >= eSIZE_GlobalSysSolnType)
            {
                ASSERTL0(false,"Matrix solution type not defined");
            }
            std::string vSolnType = MultiRegions::GlobalSysSolnTypeMap[vType];

            return GetGlobalLinSysFactory().CreateInstance( vSolnType, mkey,
                                                        vExpList,  locToGloMap);
        }

        GlobalLinSysSharedPtr ExpList::GenGlobalBndLinSys(
                    const GlobalLinSysKey     &mkey,
                    const AssemblyMapSharedPtr &locToGloMap)
        {
            boost::shared_ptr<ExpList> vExpList = GetSharedThisPtr();
            const map<int,RobinBCInfoSharedPtr> vRobinBCInfo = GetRobinBCInfo();

            MultiRegions::GlobalSysSolnType vType = mkey.GetGlobalSysSolnType();

            if (vType >= eSIZE_GlobalSysSolnType)
            {
                ASSERTL0(false,"Matrix solution type not defined");
            }
            std::string vSolnType = MultiRegions::GlobalSysSolnTypeMap[vType];

            return GetGlobalLinSysFactory().CreateInstance(vSolnType,mkey,
                                                        vExpList,locToGloMap);
        }


        /**
         * Given the elemental coefficients \f$\hat{u}_n^e\f$ of
         * an expansion, this function evaluates the spectral/hp
         * expansion \f$u^{\delta}(\boldsymbol{x})\f$ at the
         * quadrature points \f$\boldsymbol{x}_i\f$. The operation
         * is evaluated locally by the elemental function
         * StdRegions#StdExpansion#BwdTrans.
         *
         * @param   inarray         An array of size \f$N_{\mathrm{eof}}\f$
         *                          containing the local coefficients
         *                          \f$\hat{u}_n^e\f$.
         * @param   outarray        The resulting physical values at the
         *                          quadrature points
         *                          \f$u^{\delta}(\boldsymbol{x}_i)\f$
         *                          will be stored in this array of size
         *                          \f$Q_{\mathrm{tot}}\f$.
         */
        void ExpList::v_BwdTrans_IterPerExp(const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, NekDouble> tmp;
            for (int i = 0; i < m_collections.size(); ++i)
            {
                m_collections[i].ApplyOperator(Collections::eBwdTrans,
                                               inarray + m_coll_coeff_offset[i],
                                               tmp = outarray + m_coll_phys_offset[i]);
            }
        }

        LocalRegions::ExpansionSharedPtr& ExpList::GetExp(
                    const Array<OneD, const NekDouble> &gloCoord)
        {
            Array<OneD, NekDouble> stdCoord(GetCoordim(0),0.0);
            for (int i = 0; i < (*m_exp).size(); ++i)
            {
                if ((*m_exp)[i]->GetGeom()->ContainsPoint(gloCoord))
                {
                    return (*m_exp)[i];
                }
            }
            ASSERTL0(false, "Cannot find element for this point.");
            return (*m_exp)[0]; // avoid warnings
        }


        /**
         * @todo need a smarter search here that first just looks at bounding
         * vertices - suggest first seeing if point is within 10% of
         * region defined by vertices. The do point search.
         */
        int ExpList::GetExpIndex(
                                 const Array<OneD, const NekDouble> &gloCoord,
                                 NekDouble tol,
                                 bool returnNearestElmt)
        {
            Array<OneD, NekDouble> Lcoords(gloCoord.num_elements());

            return GetExpIndex(gloCoord,Lcoords,tol,returnNearestElmt);
        }


        int ExpList::GetExpIndex(const Array<OneD, const NekDouble> &gloCoords,
                                 Array<OneD, NekDouble> &locCoords,
                                 NekDouble tol,
                                 bool returnNearestElmt)
        {
            NekDouble nearpt = 1e6;

            if (GetNumElmts() == 0)
            {
                return -1;
            }
            std::vector<std::pair<int,NekDouble> > elmtIdDist;

            // Manifold case (point may match multiple elements)
            if (GetExp(0)->GetCoordim() > GetExp(0)->GetShapeDimension())
            {
                SpatialDomains::PointGeomSharedPtr v;
                SpatialDomains::PointGeom w;
                NekDouble dist = 0.0;

                // Scan all elements and store those which may contain the point
                for (int i = 0; i < (*m_exp).size(); ++i)
                {
                    if ((*m_exp)[i]->GetGeom()->ContainsPoint(gloCoords,
                                                              locCoords,
                                                              tol, nearpt))
                    {
                        w.SetX(gloCoords[0]);
                        w.SetY(gloCoords[1]);
                        w.SetZ(gloCoords[2]);

                        // Find closest vertex
                        for (int j = 0; j < (*m_exp)[i]->GetNverts(); ++j) {
                            v = m_graph->GetVertex(
                                            (*m_exp)[i]->GetGeom()->GetVid(j));
                            if (j == 0 || dist > v->dist(w))
                            {
                                dist = v->dist(w);
                            }
                        }
                        elmtIdDist.push_back(
                                    std::pair<int, NekDouble>(i, dist));
                    }
                }

                // Find nearest element
                if (!elmtIdDist.empty())
                {
                    int         min_id = elmtIdDist[0].first;
                    NekDouble   min_d  = elmtIdDist[0].second;

                    for (int i = 1; i < elmtIdDist.size(); ++i)
                    {
                        if (elmtIdDist[i].second < min_d) {
                            min_id = elmtIdDist[i].first;
                            min_d = elmtIdDist[i].second;
                        }
                    }

                    // retrieve local coordinate of point
                    (*m_exp)[min_id]->GetGeom()->GetLocCoords(gloCoords,
                                                              locCoords);
                    return min_id;
                }
                else
                {
                    return -1;
                }
            }
            // non-embedded mesh (point can only match one element)
            else
            {
                static int start = 0;
                int min_id  = 0;
                NekDouble nearpt_min = 1e6;
                Array<OneD, NekDouble> savLocCoords(locCoords.num_elements());

                // restart search from last found value
                for (int i = start; i < (*m_exp).size(); ++i)
                {
                    if ((*m_exp)[i]->GetGeom()->ContainsPoint(gloCoords, 
                                                              locCoords,
                                                              tol, nearpt))
                    {
                        start = i;
                        return i;
                    }
                    else
                    {
                        if(nearpt < nearpt_min)
                        {
                            min_id    = i;
                            nearpt_min = nearpt;
                            Vmath::Vcopy(locCoords.num_elements(),locCoords,1,savLocCoords,1);
                        }
                    }
                }

                for (int i = 0; i < start; ++i)
                {
                    if ((*m_exp)[i]->GetGeom()->ContainsPoint(gloCoords, 
                                                              locCoords,
                                                              tol, nearpt))
                    {
                        start = i;
                        return i;
                    }
                    else
                    {
                        if(nearpt < nearpt_min)
                        {
                            min_id    = i;
                            nearpt_min = nearpt;
                            Vmath::Vcopy(locCoords.num_elements(),
                                         locCoords,1,savLocCoords,1);
                        }
                    }
                }

                std::string msg = "Failed to find point within element to tolerance of "
                    + boost::lexical_cast<std::string>(tol)
                    + " using local point ("
                    + boost::lexical_cast<std::string>(locCoords[0]) +","
                    + boost::lexical_cast<std::string>(locCoords[1]) +","
                    + boost::lexical_cast<std::string>(locCoords[1]) 
                    + ") in element: "
                    + boost::lexical_cast<std::string>(min_id);
                WARNINGL1(false,msg.c_str());

                if(returnNearestElmt)
                {
                    Vmath::Vcopy(locCoords.num_elements(),savLocCoords,1,locCoords,1);
                    return min_id;
                }
                else
                {
                    return -1;
                }

            }
        }


        /**
         * Configures geometric info, such as tangent direction, on each
         * expansion.
         * @param   graph2D         Mesh
         */
        void ExpList::ApplyGeomInfo()
        {

        }

        /**
         * @brief Reset geometry information, metrics, matrix managers and
         * geometry information.
         *
         * This routine clears all matrix managers and resets all geometry
         * information, which allows the geometry information to be dynamically
         * updated as the solver is run.
         */
        void ExpList::v_Reset()
        {
            // Reset matrix managers.
            LibUtilities::NekManager<LocalRegions::MatrixKey,
                DNekScalMat, LocalRegions::MatrixKey::opLess>::ClearManager();
            LibUtilities::NekManager<LocalRegions::MatrixKey,
                DNekScalBlkMat, LocalRegions::MatrixKey::opLess>::ClearManager();

            // Loop over all elements and reset geometry information.
            for (int i = 0; i < m_exp->size(); ++i)
            {
                (*m_exp)[i]->GetGeom()->Reset(m_graph->GetCurvedEdges(),
                                              m_graph->GetCurvedFaces());
            }

            // Loop over all elements and rebuild geometric factors.
            for (int i = 0; i < m_exp->size(); ++i)
            {
                (*m_exp)[i]->Reset();
            }
        }

        /**
         * Write Tecplot Files Header
         * @param   outfile Output file name.
         * @param   var                 variables names
         */
        void ExpList::v_WriteTecplotHeader(std::ostream &outfile,
                                           std::string    var)
        {
            if (GetNumElmts() == 0)
            {
                return;
            }

            int coordim  = GetExp(0)->GetCoordim();
            char vars[3] = { 'x', 'y', 'z' };

            if (m_expType == e3DH1D)
            {
                coordim += 1;
            }
            else if (m_expType == e3DH2D)
            {
                coordim += 2;
            }

            outfile << "Variables = x";
            for (int i = 1; i < coordim; ++i)
            {
                outfile << ", " << vars[i];
            }

            if (var.size() > 0)
            {
                outfile << ", " << var;
            }

            outfile << std::endl << std::endl;
        }

        /**
         * Write Tecplot Files Zone
         * @param   outfile    Output file name.
         * @param   expansion  Expansion that is considered
         */
        void ExpList::v_WriteTecplotZone(std::ostream &outfile, int expansion)
        {
            int i, j;
            int coordim = GetCoordim(0);
            int nPoints = GetTotPoints();
            int nBases  = (*m_exp)[0]->GetNumBases();
            int numBlocks = 0;

            Array<OneD, Array<OneD, NekDouble> > coords(3);

            if (expansion == -1)
            {
                nPoints = GetTotPoints();

                coords[0] = Array<OneD, NekDouble>(nPoints);
                coords[1] = Array<OneD, NekDouble>(nPoints);
                coords[2] = Array<OneD, NekDouble>(nPoints);

                GetCoords(coords[0], coords[1], coords[2]);

                for (i = 0; i < m_exp->size(); ++i)
                {
                    int numInt = 1;

                    for (j = 0; j < nBases; ++j)
                    {
                        numInt *= (*m_exp)[i]->GetNumPoints(j)-1;
                    }

                    numBlocks += numInt;
                }
            }
            else
            {
                nPoints = (*m_exp)[expansion]->GetTotPoints();

                coords[0] = Array<OneD, NekDouble>(nPoints);
                coords[1] = Array<OneD, NekDouble>(nPoints);
                coords[2] = Array<OneD, NekDouble>(nPoints);

                (*m_exp)[expansion]->GetCoords(coords[0], coords[1], coords[2]);

                numBlocks = 1;
                for (j = 0; j < nBases; ++j)
                {
                    numBlocks *= (*m_exp)[expansion]->GetNumPoints(j)-1;
                }
            }

            if (m_expType == e3DH1D)
            {
                nBases += 1;
                coordim += 1;
                int nPlanes = GetZIDs().num_elements();
                NekDouble tmp = numBlocks * (nPlanes-1.0) / nPlanes;
                numBlocks = (int)tmp;
            }
            else if (m_expType == e3DH2D)
            {
                nBases    += 2;
                coordim += 1;
            }

            outfile << "Zone, N=" << nPoints << ", E="
                    << numBlocks << ", F=FEBlock" ;

            switch(nBases)
            {
                case 2:
                    outfile << ", ET=QUADRILATERAL" << std::endl;
                    break;
                case 3:
                    outfile << ", ET=BRICK" << std::endl;
                    break;
                default:
                    ASSERTL0(false,"Not set up for this type of output");
                    break;
            }

            // Write out coordinates
            for (j = 0; j < coordim; ++j)
            {
                for (i = 0; i < nPoints; ++i)
                {
                    outfile << coords[j][i] << " ";
                    if (i % 1000 == 0 && i)
                    {
                        outfile << std::endl;
                    }
                }
                outfile << std::endl;
            }
        }

        void ExpList::v_WriteTecplotConnectivity(std::ostream &outfile,
                                                 int expansion)
        {
            int i,j,k,l;
            int nbase = (*m_exp)[0]->GetNumBases();
            int cnt = 0;

            boost::shared_ptr<LocalRegions::ExpansionVector> exp = m_exp;

            if (expansion != -1)
            {
                exp = boost::shared_ptr<LocalRegions::ExpansionVector>(
                    new LocalRegions::ExpansionVector(1));
                (*exp)[0] = (*m_exp)[expansion];
            }

            if (nbase == 2)
            {
                for(i = 0; i < (*exp).size(); ++i)
                {
                    const int np0 = (*exp)[i]->GetNumPoints(0);
                    const int np1 = (*exp)[i]->GetNumPoints(1);

                    for(j = 1; j < np1; ++j)
                    {
                        for(k = 1; k < np0; ++k)
                        {
                            outfile << cnt + (j-1)*np0 + k   << " ";
                            outfile << cnt + (j-1)*np0 + k+1 << " ";
                            outfile << cnt +  j   *np0 + k+1 << " ";
                            outfile << cnt +  j   *np0 + k   << endl;
                        }
                    }

                    cnt += np0*np1;
                }
            }
            else if (nbase == 3)
            {
                for(i = 0; i < (*exp).size(); ++i)
                {
                    const int np0 = (*exp)[i]->GetNumPoints(0);
                    const int np1 = (*exp)[i]->GetNumPoints(1);
                    const int np2 = (*exp)[i]->GetNumPoints(2);
                    const int np01 = np0*np1;

                    for(j = 1; j < np2; ++j)
                    {
                        for(k = 1; k < np1; ++k)
                        {
                            for(l = 1; l < np0; ++l)
                            {
                                outfile << cnt + (j-1)*np01 + (k-1)*np0 + l   << " ";
                                outfile << cnt + (j-1)*np01 + (k-1)*np0 + l+1 << " ";
                                outfile << cnt + (j-1)*np01 +  k   *np0 + l+1 << " ";
                                outfile << cnt + (j-1)*np01 +  k   *np0 + l   << " ";
                                outfile << cnt +  j   *np01 + (k-1)*np0 + l   << " ";
                                outfile << cnt +  j   *np01 + (k-1)*np0 + l+1 << " ";
                                outfile << cnt +  j   *np01 +  k   *np0 + l+1 << " ";
                                outfile << cnt +  j   *np01 +  k   *np0 + l   << endl;
                            }
                        }
                    }
                    cnt += np0*np1*np2;
                }
            }
            else
            {
                ASSERTL0(false,"Not set up for this dimension");
            }
        }

        /**
         * Write Tecplot Files Field
         * @param   outfile    Output file name.
         * @param   expansion  Expansion that is considered
         */
        void ExpList::v_WriteTecplotField(std::ostream &outfile, int expansion)
        {
            if (expansion == -1)
            {
                int totpoints = GetTotPoints();
                if(m_physState == false)
                {
                    BwdTrans(m_coeffs,m_phys);
                }

                for(int i = 0; i < totpoints; ++i)
                {
                    outfile << m_phys[i] << " ";
                    if(i % 1000 == 0 && i)
                    {
                        outfile << std::endl;
                    }
                }
                outfile << std::endl;

            }
            else
            {
                int nPoints = (*m_exp)[expansion]->GetTotPoints();

                for (int i = 0; i < nPoints; ++i)
                {
                    outfile << m_phys[i + m_phys_offset[expansion]] << " ";
                }

                outfile << std::endl;
            }
        }

        void ExpList::WriteVtkHeader(std::ostream &outfile)
        {
            outfile << "<?xml version=\"1.0\"?>" << endl;
            outfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
                    << "byte_order=\"LittleEndian\">" << endl;
            outfile << "  <UnstructuredGrid>" << endl;
        }

        void ExpList::WriteVtkFooter(std::ostream &outfile)
        {
            outfile << "  </UnstructuredGrid>" << endl;
            outfile << "</VTKFile>" << endl;
        }

        void ExpList::v_WriteVtkPieceHeader(std::ostream &outfile, int expansion, int istrip)
        {
            ASSERTL0(false, "Routine not implemented for this expansion.");
        }

        void ExpList::WriteVtkPieceFooter(std::ostream &outfile, int expansion)
        {
            outfile << "      </PointData>" << endl;
            outfile << "    </Piece>" << endl;
        }

        void ExpList::v_WriteVtkPieceData(std::ostream &outfile, int expansion,
                                        std::string var)
        {
            int i;
            int nq = (*m_exp)[expansion]->GetTotPoints();

            // printing the fields of that zone
            outfile << "        <DataArray type=\"Float64\" Name=\""
                    << var << "\">" << endl;
            outfile << "          ";
            const Array<OneD, NekDouble> phys = m_phys + m_phys_offset[expansion];
            for(i = 0; i < nq; ++i)
            {
                outfile << (fabs(phys[i]) < NekConstants::kNekZeroTol ? 0 : phys[i]) << " ";
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
        }

        /**
         * Given a spectral/hp approximation
         * \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the quadrature points
         * (which should be contained in #m_phys), this function calculates the
         * \f$L_\infty\f$ error of this approximation with respect to an exact
         * solution. The local distribution of the quadrature points allows an
         * elemental evaluation of this operation through the functions
         * StdRegions#StdExpansion#Linf.
         *
         * The exact solution, also evaluated at the quadrature
         * points, should be contained in the variable #m_phys of
         * the ExpList object \a Sol.
         *
         * @param   soln            A 1D array, containing the discrete
         *                          evaluation of the exact solution at the
         *                          quadrature points in its array #m_phys.
         * @return  The \f$L_\infty\f$ error of the approximation.
         */
        NekDouble ExpList::Linf(
            const Array<OneD, const NekDouble> &inarray,
            const Array<OneD, const NekDouble> &soln)
        {
            NekDouble err = 0.0;

            if (soln == NullNekDouble1DArray)
            {
                err = Vmath::Vmax(m_npoints, inarray, 1);
            }
            else
            {
                for (int i = 0; i < m_npoints; ++i)
                {
                    err = max(err, abs(inarray[i] - soln[i]));
                }
            }

            m_comm->GetRowComm()->AllReduce(err, LibUtilities::ReduceMax);

            return err;
        }

        /**
         * Given a spectral/hp approximation \f$u^{\delta}(\boldsymbol{x})\f$
         * evaluated at the quadrature points (which should be contained in
         * #m_phys), this function calculates the \f$L_2\f$ error of this
         * approximation with respect to an exact solution. The local
         * distribution of the quadrature points allows an elemental evaluation
         * of this operation through the functions StdRegions#StdExpansion#L2.
         *
         * The exact solution, also evaluated at the quadrature points, should
         * be contained in the variable #m_phys of the ExpList object \a Sol.
         *
         * @param   Sol             An ExpList, containing the discrete
         *                          evaluation of the exact solution at the
         *                          quadrature points in its array #m_phys.
         * @return  The \f$L_2\f$ error of the approximation.
         */
        NekDouble ExpList::v_L2(
            const Array<OneD, const NekDouble> &inarray,
            const Array<OneD, const NekDouble> &soln)
        {
            NekDouble err = 0.0, errl2;
            int    i;

            if (soln == NullNekDouble1DArray)
            {
                for (i = 0; i < (*m_exp).size(); ++i)
                {
                    errl2 = (*m_exp)[i]->L2(inarray + m_phys_offset[i]);
                    err += errl2*errl2;
                }
            }
            else
            {
                for (i = 0; i < (*m_exp).size(); ++i)
                {
                    errl2 = (*m_exp)[i]->L2(inarray + m_phys_offset[i],
                                            soln    + m_phys_offset[i]);
                    err += errl2*errl2;
                }
            }

            m_comm->GetRowComm()->AllReduce(err, LibUtilities::ReduceSum);

            return sqrt(err);
        }

        NekDouble ExpList::v_Integral(const Array<OneD, const NekDouble> &inarray)
        {
            NekDouble err = 0.0;
            int       i   = 0;

            for (i = 0; i < (*m_exp).size(); ++i)
            {
                err += (*m_exp)[m_offset_elmt_id[i]]->Integral(inarray + m_phys_offset[i]);
            }
            m_comm->GetRowComm()->AllReduce(err, LibUtilities::ReduceSum);

            return err;
        }

        Array<OneD, const NekDouble> ExpList::v_HomogeneousEnergy (void)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            Array<OneD, NekDouble> NoEnergy(1,0.0);
            return NoEnergy;
        }

        LibUtilities::TranspositionSharedPtr ExpList::v_GetTransposition(void)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            LibUtilities::TranspositionSharedPtr trans;
            return trans;
        }

        NekDouble ExpList::v_GetHomoLen(void)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            NekDouble len = 0.0;
            return len;
        }

        Array<OneD, const unsigned int> ExpList::v_GetZIDs(void)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            Array<OneD, unsigned int> NoModes(1);
            return NoModes;
        }

        Array<OneD, const unsigned int> ExpList::v_GetYIDs(void)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            Array<OneD, unsigned int> NoModes(1);
            return NoModes;
        }


        void ExpList::v_PhysInterp1DScaled(const NekDouble scale, const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_PhysGalerkinProjection1DScaled(const NekDouble scale, const Array<OneD, NekDouble> &inarray, Array<OneD, NekDouble> &outarray)        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::ExtractFileBCs(
            const std::string               &fileName,
            const std::string               &varName,
            const boost::shared_ptr<ExpList> locExpList)
        {
            string varString = fileName.substr(0, fileName.find_last_of("."));
            int j, k, len = varString.length();
            varString = varString.substr(len-1, len);

            std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
            std::vector<std::vector<NekDouble> > FieldData;

            LibUtilities::FieldIO f(m_session->GetComm());
            f.Import(fileName, FieldDef, FieldData);

            bool found = false;
            for (j = 0; j < FieldDef.size(); ++j)
            {
                for (k = 0; k < FieldDef[j]->m_fields.size(); ++k)
                {
                    if (FieldDef[j]->m_fields[k] == varName)
                    {
                        // Copy FieldData into locExpList
                        locExpList->ExtractDataToCoeffs(
                            FieldDef[j], FieldData[j],
                            FieldDef[j]->m_fields[k],
                            locExpList->UpdateCoeffs());
                        found = true;
                    }
                }
            }

            ASSERTL0(found, "Could not find variable '"+varName+
                            "' in file boundary condition "+fileName);
            locExpList->BwdTrans_IterPerExp(
                locExpList->GetCoeffs(),
                locExpList->UpdatePhys());
        }

        /**
         * Given a spectral/hp approximation
         * \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the quadrature points
         * (which should be contained in #m_phys), this function calculates the
         * \f$H^1_2\f$ error of this approximation with respect to an exact
         * solution. The local distribution of the quadrature points allows an
         * elemental evaluation of this operation through the functions
         * StdRegions#StdExpansion#H1.
         *
         * The exact solution, also evaluated at the quadrature points, should
         * be contained in the variable #m_phys of the ExpList object \a Sol.
         *
         * @param   soln        An 1D array, containing the discrete evaluation
         *                      of the exact solution at the quadrature points.
         *
         * @return  The \f$H^1_2\f$ error of the approximation.
         */
        NekDouble ExpList::H1(
            const Array<OneD, const NekDouble> &inarray,
            const Array<OneD, const NekDouble> &soln)
        {
            NekDouble err = 0.0, errh1;
            int    i;

            for (i = 0; i < (*m_exp).size(); ++i)
            {
                errh1 = (*m_exp)[i]->H1(inarray + m_phys_offset[i],
                                        soln    + m_phys_offset[i]);
                err += errh1*errh1;
            }

            m_comm->GetRowComm()->AllReduce(err, LibUtilities::ReduceSum);

            return sqrt(err);
        }

        void  ExpList::GeneralGetFieldDefinitions(std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef,
                                                  int NumHomoDir,
                                                  int NumHomoStrip,
                                                  Array<OneD, LibUtilities::BasisSharedPtr> &HomoBasis,
                                                  std::vector<NekDouble> &HomoLen,
                                                  std::vector<unsigned int> &HomoZIDs,
                                                  std::vector<unsigned int> &HomoYIDs)
        {
            int startenum = (int) LibUtilities::eSegment;
            int endenum   = (int) LibUtilities::eHexahedron;
            int s         = 0;
            LibUtilities::ShapeType shape;

            ASSERTL1(NumHomoDir == HomoBasis.num_elements(),"Homogeneous basis is not the same length as NumHomoDir");
            ASSERTL1(NumHomoDir == HomoLen.size(),"Homogeneous length vector is not the same length as NumHomDir");

            // count number of shapes
            switch((*m_exp)[0]->GetShapeDimension())
            {
            case 1:
                startenum = (int) LibUtilities::eSegment;
                endenum   = (int) LibUtilities::eSegment;
                break;
            case 2:
                startenum = (int) LibUtilities::eTriangle;
                endenum   = (int) LibUtilities::eQuadrilateral;
                break;
            case 3:
                startenum = (int) LibUtilities::eTetrahedron;
                endenum   = (int) LibUtilities::eHexahedron;
                break;
            }

            for(s = startenum; s <= endenum; ++s)
            {
                std::vector<unsigned int>             elementIDs;
                std::vector<LibUtilities::BasisType>  basis;
                std::vector<unsigned int>             numModes;
                std::vector<std::string>              fields;

                bool first    = true;
                bool UniOrder = true;
                int n;

                shape = (LibUtilities::ShapeType) s;

                for(int i = 0; i < (*m_exp).size(); ++i)
                {
                    if((*m_exp)[i]->GetGeom()->GetShapeType() == shape)
                    {
                        elementIDs.push_back((*m_exp)[i]->GetGeom()->GetGlobalID());
                        if(first)
                        {
                            for(int j = 0; j < (*m_exp)[i]->GetNumBases(); ++j)
                            {
                                basis.push_back((*m_exp)[i]->GetBasis(j)->GetBasisType());
                                numModes.push_back((*m_exp)[i]->GetBasis(j)->GetNumModes());
                            }

                            // add homogeneous direction details if defined
                            for(n = 0 ; n < NumHomoDir; ++n)
                            {
                                basis.push_back(HomoBasis[n]->GetBasisType());
                                numModes.push_back(HomoBasis[n]->GetNumModes());
                            }

                            first = false;
                        }
                        else
                        {
                            ASSERTL0((*m_exp)[i]->GetBasis(0)->GetBasisType() == basis[0],"Routine is not set up for multiple bases definitions");

                            for(int j = 0; j < (*m_exp)[i]->GetNumBases(); ++j)
                            {
                                numModes.push_back((*m_exp)[i]->GetBasis(j)->GetNumModes());
                                if(numModes[j] != (*m_exp)[i]->GetBasis(j)->GetNumModes())
                                {
                                    UniOrder = false;
                                }
                            }
                            // add homogeneous direction details if defined
                            for(n = 0 ; n < NumHomoDir; ++n)
                            {
                                numModes.push_back(HomoBasis[n]->GetNumModes());
                            }
                        }
                    }
                }


                if(elementIDs.size() > 0)
                {
                    for(int i = 0; i < NumHomoStrip; ++i)
                    {
                        LibUtilities::FieldDefinitionsSharedPtr fdef  =
                            MemoryManager<LibUtilities::FieldDefinitions>::
                                AllocateSharedPtr(shape, elementIDs, basis,
                                                  UniOrder, numModes,fields,
                                                  NumHomoDir, HomoLen, HomoZIDs,
                                                  HomoYIDs);
                        fielddef.push_back(fdef);
                    }
                }
            }
        }


        //
        // Virtual functions
        //
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> ExpList::v_GetFieldDefinitions()
        {
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> returnval;
            v_GetFieldDefinitions(returnval);
            return returnval;
        }

        void  ExpList::v_GetFieldDefinitions(std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef)
        {
            GeneralGetFieldDefinitions(fielddef);
        }

        //Append the element data listed in elements
        //fielddef->m_ElementIDs onto fielddata
        void ExpList::v_AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata)
        {
            v_AppendFieldData(fielddef,fielddata,m_coeffs);
        }

        void ExpList::v_AppendFieldData(LibUtilities::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, Array<OneD, NekDouble> &coeffs)
        {
            int i;
            // Determine mapping from element ids to location in
            // expansion list
            // Determine mapping from element ids to location in
            // expansion list
            map<int, int> ElmtID_to_ExpID;
            for(i = 0; i < (*m_exp).size(); ++i)
            {
                ElmtID_to_ExpID[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
            }

            for(i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                int eid     = ElmtID_to_ExpID[fielddef->m_elementIDs[i]];
                int datalen = (*m_exp)[eid]->GetNcoeffs();
                fielddata.insert(fielddata.end(),&coeffs[m_coeff_offset[eid]],&coeffs[m_coeff_offset[eid]]+datalen);
            }

        }

        /// Extract the data in fielddata into the coeffs
        void ExpList::ExtractDataToCoeffs(
                                   LibUtilities::FieldDefinitionsSharedPtr &fielddef,
                                   std::vector<NekDouble> &fielddata,
                                   std::string &field,
                                   Array<OneD, NekDouble> &coeffs)
        {
            v_ExtractDataToCoeffs(fielddef,fielddata,field,coeffs);
        }

        void ExpList::ExtractCoeffsToCoeffs(const boost::shared_ptr<ExpList> &fromExpList, const Array<OneD, const NekDouble> &fromCoeffs, Array<OneD, NekDouble> &toCoeffs)
        {
            v_ExtractCoeffsToCoeffs(fromExpList,fromCoeffs,toCoeffs);
        }

        /**
         * @brief Extract data from raw field data into expansion list.
         *
         * @param fielddef   Field definitions.
         * @param fielddata  Data for associated field.
         * @param field      Field variable name.
         * @param coeffs     Resulting coefficient array.
         */
        void ExpList::v_ExtractDataToCoeffs(
            LibUtilities::FieldDefinitionsSharedPtr   &fielddef,
            std::vector<NekDouble>                    &fielddata,
            std::string                               &field,
            Array<OneD, NekDouble>                    &coeffs)
        {
            int i, expId;
            int offset       = 0;
            int modes_offset = 0;
            int datalen      = fielddata.size()/fielddef->m_fields.size();

            // Find data location according to field definition
            for(i = 0; i < fielddef->m_fields.size(); ++i)
            {
                if(fielddef->m_fields[i] == field)
                {
                    break;
                }
                offset += datalen;
            }

            ASSERTL0(i != fielddef->m_fields.size(),
                     "Field (" + field + ") not found in file.");

            // Determine mapping from element ids to location in expansion list
            map<int, int> elmtToExpId;

            // Loop in reverse order so that in case where using a Homogeneous
            // expansion it sets geometry ids to first part of m_exp
            // list. Otherwise will set to second (complex) expansion
            for(i = (*m_exp).size()-1; i >= 0; --i)
            {
                elmtToExpId[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
            }

            for (i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                // Reset modes_offset in the case where all expansions of
                // the same order.
                if (fielddef->m_uniOrder == true)
                {
                    modes_offset = 0;
                }

                datalen = LibUtilities::GetNumberOfCoefficients(fielddef->m_shapeType,
                                                                fielddef->m_numModes, modes_offset);

                const int elmtId = fielddef->m_elementIDs[i];
                if (elmtToExpId.count(elmtId) == 0)
                {
                    offset += datalen;
                    modes_offset += (*m_exp)[0]->GetNumBases();
                    continue;
                }

                expId   = elmtToExpId[elmtId];

                if (datalen == (*m_exp)[expId]->GetNcoeffs())
                {
                    Vmath::Vcopy(datalen, &fielddata[offset], 1,
                                 &coeffs[m_coeff_offset[expId]], 1);
                }
                else
                {
                    (*m_exp)[expId]->ExtractDataToCoeffs(
                                                         &fielddata[offset], fielddef->m_numModes,
                                                         modes_offset, &coeffs[m_coeff_offset[expId]]);
                }

                offset += datalen;
                modes_offset += (*m_exp)[0]->GetNumBases();
            }

            return;
        }

        void ExpList::v_ExtractCoeffsToCoeffs(const boost::shared_ptr<ExpList> &fromExpList, const Array<OneD, const NekDouble> &fromCoeffs, Array<OneD, NekDouble> &toCoeffs)
        {
            int i;
            int offset = 0;

            for(i = 0; i < (*m_exp).size(); ++i)
            {
                std::vector<unsigned int> nummodes;
                int eid = m_offset_elmt_id[i];
                for(int j= 0; j < fromExpList->GetExp(eid)->GetNumBases(); ++j)
                {
                    nummodes.push_back(fromExpList->GetExp(eid)->GetBasisNumModes(j));
                }

                (*m_exp)[eid]->ExtractDataToCoeffs(&fromCoeffs[offset], nummodes,0,
                                                   &toCoeffs[m_coeff_offset[eid]]);

                offset += fromExpList->GetExp(eid)->GetNcoeffs();
            }
        }


        const Array<OneD,const boost::shared_ptr<ExpList> >
                                        &ExpList::v_GetBndCondExpansions(void)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            static Array<OneD,const boost::shared_ptr<ExpList> > result;
            return result;
        }

        boost::shared_ptr<ExpList>  &ExpList::v_UpdateBndCondExpansion(int i)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            static boost::shared_ptr<ExpList> result;
            return result;
        }

        void ExpList::v_Upwind(
            const Array<OneD, const Array<OneD,       NekDouble> > &Vec,
            const Array<OneD,                   const NekDouble>   &Fwd,
            const Array<OneD,                   const NekDouble>   &Bwd,
                  Array<OneD,                         NekDouble>   &Upwind)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_Upwind(
            const Array<OneD, const NekDouble> &Vn,
            const Array<OneD, const NekDouble> &Fwd,
            const Array<OneD, const NekDouble> &Bwd,
                  Array<OneD,       NekDouble> &Upwind)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        boost::shared_ptr<ExpList> &ExpList::v_GetTrace()
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            static boost::shared_ptr<ExpList> returnVal;
            return returnVal;
        }

        boost::shared_ptr<AssemblyMapDG> &ExpList::v_GetTraceMap()
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            static boost::shared_ptr<AssemblyMapDG> result;
            return result;
        }

        const Array<OneD, const int> &ExpList::v_GetTraceBndMap()
        {
            return GetTraceMap()->GetBndCondTraceToGlobalTraceMap();
        }

        void ExpList::v_GetNormals(
            Array<OneD, Array<OneD, NekDouble> > &normals)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_AddTraceIntegral(
                                const Array<OneD, const NekDouble> &Fx,
                                const Array<OneD, const NekDouble> &Fy,
                                      Array<OneD, NekDouble> &outarray)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_AddTraceIntegral(
                                const Array<OneD, const NekDouble> &Fn,
                                      Array<OneD, NekDouble> &outarray)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_AddFwdBwdTraceIntegral(
                                const Array<OneD, const NekDouble> &Fwd,
                                const Array<OneD, const NekDouble> &Bwd,
                                      Array<OneD, NekDouble> &outarray)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd,
                                           Array<OneD,NekDouble> &Bwd)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_GetFwdBwdTracePhys(
                                const Array<OneD,const NekDouble>  &field,
                                      Array<OneD,NekDouble> &Fwd,
                                      Array<OneD,NekDouble> &Bwd)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_ExtractTracePhys(Array<OneD,NekDouble> &outarray)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_ExtractTracePhys(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_MultiplyByInvMassMatrix(
                                const Array<OneD,const NekDouble> &inarray,
                                Array<OneD,      NekDouble> &outarray,
                                CoeffState coeffstate)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_HelmSolve(
                const Array<OneD, const NekDouble> &inarray,
                      Array<OneD,       NekDouble> &outarray,
                const FlagList &flags,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff,
                const Array<OneD, const NekDouble> &dirForcing)
        {
            ASSERTL0(false, "HelmSolve not implemented.");
        }

        void ExpList::v_LinearAdvectionDiffusionReactionSolve(
                       const Array<OneD, Array<OneD, NekDouble> > &velocity,
                       const Array<OneD, const NekDouble> &inarray,
                       Array<OneD, NekDouble> &outarray,
                       const NekDouble lambda,
                       CoeffState coeffstate,
                       const Array<OneD, const NekDouble>&  dirForcing)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_LinearAdvectionReactionSolve(
                       const Array<OneD, Array<OneD, NekDouble> > &velocity,
                       const Array<OneD, const NekDouble> &inarray,
                       Array<OneD, NekDouble> &outarray,
                       const NekDouble lambda,
                       CoeffState coeffstate,
                       const Array<OneD, const NekDouble>&  dirForcing)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD, NekDouble> &outarray,
                                            CoeffState coeffstate,
                                            bool Shuff,
                                            bool UnShuff)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD, NekDouble> &outarray,
                                            CoeffState coeffstate,
                                            bool Shuff,
                                            bool UnShuff)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_DealiasedProd(const Array<OneD, NekDouble> &inarray1,const Array<OneD, NekDouble> &inarray2,Array<OneD, NekDouble> &outarray,CoeffState coeffstate)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_GetBCValues(Array<OneD, NekDouble> &BndVals,
                                    const Array<OneD, NekDouble> &TotField,
                                    int BndID)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_NormVectorIProductWRTBase(Array<OneD, const NekDouble> &V1,
                                                  Array<OneD, const NekDouble> &V2,
                                                  Array<OneD, NekDouble> &outarray,
                                                  int BndID)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }
        
        void ExpList::v_NormVectorIProductWRTBase(
                                Array<OneD, Array<OneD, NekDouble> > &V,
                                Array<OneD, NekDouble> &outarray)
        {
            Array< OneD, NekDouble> tmp;
            switch (GetCoordim(0))
            {
                case 1:
                {
                    for(int i = 0; i < GetExpSize(); ++i)
                    {
                        (*m_exp)[i]->NormVectorIProductWRTBase(
                                        V[0] + GetPhys_Offset(i), 
                                        tmp = outarray + GetCoeff_Offset(i));
                    }
                }
                break;
                case 2:
                {
                    for(int i = 0; i < GetExpSize(); ++i)
                    {
                        (*m_exp)[i]->NormVectorIProductWRTBase(
                                        V[0] + GetPhys_Offset(i),
                                        V[1] + GetPhys_Offset(i),
                                        tmp = outarray + GetCoeff_Offset(i));
                    }
                }
                break;
                case 3:
                {
                    for(int i = 0; i < GetExpSize(); ++i)
                    {
                        (*m_exp)[i]->NormVectorIProductWRTBase(
                                        V[0] + GetPhys_Offset(i),
                                        V[1] + GetPhys_Offset(i),
                                        V[2] + GetPhys_Offset(i),
                                        tmp = outarray + GetCoeff_Offset(i));
                    }
                }
                break;
                default:
                    ASSERTL0(false,"Dimension not supported");
                    break;
            }
        }

        void ExpList::v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        /**
         */
        void ExpList::v_FillBndCondFromField()
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_LocalToGlobal(void)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_GlobalToLocal(void)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }


        void ExpList::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,       NekDouble> &outarray,
                                 CoeffState coeffstate)
        {
            v_BwdTrans_IterPerExp(inarray,outarray);
        }

        void ExpList::v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                 Array<OneD,       NekDouble> &outarray,
                                 CoeffState coeffstate)
        {
            v_FwdTrans_IterPerExp(inarray,outarray);
        }

        void ExpList::v_IProductWRTBase(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray,
                                CoeffState coeffstate)
        {
            Array<OneD,NekDouble>  tmp;
            for (int i = 0; i < m_collections.size(); ++i)
            {

                m_collections[i].ApplyOperator(Collections::eIProductWRTBase,
                                               inarray + m_coll_phys_offset[i],
                                               tmp = outarray + m_coll_coeff_offset[i]);
            }
        }

        void ExpList::v_GeneralMatrixOp(
                                        const GlobalMatrixKey             &gkey,
                                        const Array<OneD,const NekDouble> &inarray,
                                        Array<OneD,      NekDouble> &outarray,
                                        CoeffState coeffstate)
        {
            GeneralMatrixOp_IterPerExp(gkey,inarray,outarray);
        }

        /**
         * The operation is evaluated locally by the elemental
         * function StdRegions#StdExpansion#GetCoords.
         *
         * @param   coord_0         After calculation, the \f$x_1\f$ coordinate
         *                          will be stored in this array.
         * @param   coord_1         After calculation, the \f$x_2\f$ coordinate
         *                          will be stored in this array.
         * @param   coord_2         After calculation, the \f$x_3\f$ coordinate
         *                          will be stored in this array.
         */
        void ExpList::v_GetCoords(Array<OneD, NekDouble> &coord_0,
                                  Array<OneD, NekDouble> &coord_1,
                                  Array<OneD, NekDouble> &coord_2)
        {
            if (GetNumElmts() == 0)
            {
                return;
            }

            int    i;
            Array<OneD, NekDouble> e_coord_0;
            Array<OneD, NekDouble> e_coord_1;
            Array<OneD, NekDouble> e_coord_2;

            switch(GetExp(0)->GetCoordim())
            {
            case 1:
                for(i= 0; i < (*m_exp).size(); ++i)
                {
                    e_coord_0 = coord_0 + m_phys_offset[i];
                    (*m_exp)[i]->GetCoords(e_coord_0);
                }
                break;
            case 2:
                ASSERTL0(coord_1.num_elements() != 0,
                         "output coord_1 is not defined");

                for(i= 0; i < (*m_exp).size(); ++i)
                {
                    e_coord_0 = coord_0 + m_phys_offset[i];
                    e_coord_1 = coord_1 + m_phys_offset[i];
                    (*m_exp)[i]->GetCoords(e_coord_0,e_coord_1);
                }
                break;
            case 3:
                ASSERTL0(coord_1.num_elements() != 0,
                         "output coord_1 is not defined");
                ASSERTL0(coord_2.num_elements() != 0,
                         "output coord_2 is not defined");

                for(i= 0; i < (*m_exp).size(); ++i)
                {
                    e_coord_0 = coord_0 + m_phys_offset[i];
                    e_coord_1 = coord_1 + m_phys_offset[i];
                    e_coord_2 = coord_2 + m_phys_offset[i];
                    (*m_exp)[i]->GetCoords(e_coord_0,e_coord_1,e_coord_2);
                }
                break;
            }
        }

        /**
         */
        void ExpList::v_SetUpPhysNormals()
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        /**
         */
        void ExpList::v_GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                            Array<OneD,int> &EdgeID)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }
        
        /**
         */
        void ExpList::v_GetBndElmtExpansion(int i,
                            boost::shared_ptr<ExpList> &result)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }
        
        /**
         */
        void ExpList::v_ExtractElmtToBndPhys(int i,
                            Array<OneD, NekDouble> &element,
                            Array<OneD, NekDouble> &boundary)
        {
            int n, cnt;
            Array<OneD, NekDouble> tmp1, tmp2;
            StdRegions::StdExpansionSharedPtr elmt;
            
            Array<OneD, int> ElmtID,EdgeID;
            GetBoundaryToElmtMap(ElmtID,EdgeID);
            
            // Initialise result
            boundary = Array<OneD, NekDouble> 
                            (GetBndCondExpansions()[i]->GetTotPoints(), 0.0);
            
            // Skip other boundary regions
            for (cnt = n = 0; n < i; ++n)
            {
                cnt += GetBndCondExpansions()[n]->GetExpSize();
            }
            
            int offsetBnd;
            int offsetElmt = 0;
            for (n = 0; n < GetBndCondExpansions()[i]->GetExpSize(); ++n)
            {
                offsetBnd = GetBndCondExpansions()[i]->GetPhys_Offset(n);
                
                elmt   = GetExp(ElmtID[cnt+n]);
                elmt->GetTracePhysVals(EdgeID[cnt+n],
                                      GetBndCondExpansions()[i]->GetExp(n),
                                      tmp1 = element + offsetElmt,
                                      tmp2 = boundary + offsetBnd);
                
                offsetElmt += elmt->GetTotPoints();
            }
        }
        
        /**
         */
        void ExpList::v_ExtractPhysToBndElmt(int i,
                            const Array<OneD, const NekDouble> &phys,
                            Array<OneD, NekDouble> &bndElmt)
        {
            int n, cnt, nq;
            Array<OneD, NekDouble> tmp1, tmp2;
            
            Array<OneD, int> ElmtID,EdgeID;
            GetBoundaryToElmtMap(ElmtID,EdgeID);
            
            // Skip other boundary regions
            for (cnt = n = 0; n < i; ++n)
            {
                cnt += GetBndCondExpansions()[n]->GetExpSize();
            }
            
            // Count number of points
            int npoints = 0;
            for (n = 0; n < GetBndCondExpansions()[i]->GetExpSize(); ++n)
            {
                npoints += GetExp(ElmtID[cnt+n])->GetTotPoints();
            }
            
            // Initialise result
            bndElmt = Array<OneD, NekDouble> (npoints, 0.0);
            
            // Extract data
            int offsetPhys;
            int offsetElmt = 0;
            for (n = 0; n < GetBndCondExpansions()[i]->GetExpSize(); ++n)
            {
                nq = GetExp(ElmtID[cnt+n])->GetTotPoints();
                offsetPhys = GetPhys_Offset(ElmtID[cnt+n]);
                Vmath::Vcopy(nq, tmp1 = phys    + offsetPhys, 1,
                                 tmp2 = bndElmt + offsetElmt, 1);
                offsetElmt += nq;
            }
        }
        
        /**
         */
        void ExpList::v_GetBoundaryNormals(int i,
                        Array<OneD, Array<OneD, NekDouble> > &normals)
        {
            int j, n, cnt, nq;
            int coordim = GetCoordim(0);
            Array<OneD, NekDouble> tmp;
            StdRegions::StdExpansionSharedPtr elmt;
            
            Array<OneD, int> ElmtID,EdgeID;
            GetBoundaryToElmtMap(ElmtID,EdgeID);
            
            // Initialise result
            normals = Array<OneD, Array<OneD, NekDouble> > (coordim);
            for (j = 0; j < coordim; ++j)
            {
                normals[j] = Array<OneD, NekDouble> ( 
                                GetBndCondExpansions()[i]->GetTotPoints(), 0.0);
            }
            
            // Skip other boundary regions
            for (cnt = n = 0; n < i; ++n)
            {
                cnt += GetBndCondExpansions()[n]->GetExpSize();
            }
            
            int offset;
            for (n = 0; n < GetBndCondExpansions()[i]->GetExpSize(); ++n)
            {
                offset = GetBndCondExpansions()[i]->GetPhys_Offset(n);
                nq = GetBndCondExpansions()[i]->GetExp(n)->GetTotPoints();
                
                elmt   = GetExp(ElmtID[cnt+n]);
                const Array<OneD, const Array<OneD, NekDouble> > normalsElmt
                            = elmt->GetSurfaceNormal(EdgeID[cnt+n]);
                // Copy to result
                for (j = 0; j < coordim; ++j)
                {
                    Vmath::Vcopy(nq, normalsElmt[j], 1,
                                     tmp = normals[j] + offset, 1);
                }
            }
        }

        /**
         */
        void ExpList::v_ReadGlobalOptimizationParameters()
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        /**
         */
        const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
                                            &ExpList::v_GetBndConditions(void)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            static Array<OneD, const SpatialDomains::BoundaryConditionShPtr>
                                                                        result;
            return result;
        }

        /**
         */
        Array<OneD,SpatialDomains::BoundaryConditionShPtr> &ExpList::v_UpdateBndConditions()
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            static Array<OneD, SpatialDomains::BoundaryConditionShPtr> result;
            return result;
        }

        /**
         */
        void ExpList::v_EvaluateBoundaryConditions(
            const NekDouble time,
            const std::string varName,
            const NekDouble x2_in,
            const NekDouble x3_in)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        /**
         */
        map<int, RobinBCInfoSharedPtr> ExpList::v_GetRobinBCInfo(void)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            static map<int,RobinBCInfoSharedPtr> result;
            return result;
        }

        /**
         */
        void ExpList::v_GetPeriodicEntities(
            PeriodicMap &periodicVerts,
            PeriodicMap &periodicEdges,
            PeriodicMap &periodicFaces)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        SpatialDomains::BoundaryConditionShPtr ExpList::GetBoundaryCondition(
            const SpatialDomains::BoundaryConditionCollection& collection,
            unsigned int regionId,
            const std::string& variable)
        {
            SpatialDomains::BoundaryConditionCollection::const_iterator collectionIter = collection.find(regionId);
            ASSERTL1(collectionIter != collection.end(), "Unable to locate collection "+boost::lexical_cast<string>(regionId));
            const SpatialDomains::BoundaryConditionMapShPtr boundaryConditionMap = (*collectionIter).second;
            SpatialDomains::BoundaryConditionMap::const_iterator conditionMapIter = boundaryConditionMap->find(variable);
            ASSERTL1(conditionMapIter != boundaryConditionMap->end(), "Unable to locate condition map.");
            const SpatialDomains::BoundaryConditionShPtr boundaryCondition = (*conditionMapIter).second;
            return boundaryCondition;
        }

        ExpListSharedPtr &ExpList::v_GetPlane(int n)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            return NullExpListSharedPtr;
        }


        StdRegions::StdExpansionSharedPtr GetStdExp(StdRegions::StdExpansionSharedPtr exp)
        {

            StdRegions::StdExpansionSharedPtr stdExp;

            switch(exp->DetShapeType())
            {
            case LibUtilities::eSegment:
                stdExp = MemoryManager<StdRegions::StdSegExp>
                    ::AllocateSharedPtr(exp->GetBasis(0)->GetBasisKey());
                break;
            case LibUtilities::eTriangle:
                {
                    StdRegions::StdNodalTriExpSharedPtr nexp;
                    if((nexp = exp->as<StdRegions::StdNodalTriExp>()))
                    {
                        stdExp = MemoryManager<StdRegions::StdNodalTriExp>
                            ::AllocateSharedPtr(exp->GetBasis(0)->GetBasisKey(),
                                                exp->GetBasis(1)->GetBasisKey(),
                                                nexp->GetNodalPointsKey().GetPointsType());
                    }
                    else
                    {
                        stdExp = MemoryManager<StdRegions::StdTriExp>
                            ::AllocateSharedPtr(exp->GetBasis(0)->GetBasisKey(),
                                                exp->GetBasis(1)->GetBasisKey());
                    }
                }
                break;
            case LibUtilities::eQuadrilateral:
                stdExp = MemoryManager<StdRegions::StdQuadExp>
                    ::AllocateSharedPtr(exp->GetBasis(0)->GetBasisKey(),
                                        exp->GetBasis(1)->GetBasisKey());
                break;
            case LibUtilities::eTetrahedron:
                    stdExp = MemoryManager<StdRegions::StdTetExp>
                        ::AllocateSharedPtr(exp->GetBasis(0)->GetBasisKey(),
                                            exp->GetBasis(1)->GetBasisKey(),
                                            exp->GetBasis(2)->GetBasisKey());
                    break;
            case LibUtilities::ePyramid:
                stdExp = MemoryManager<StdRegions::StdPyrExp>
                    ::AllocateSharedPtr(exp->GetBasis(0)->GetBasisKey(),
                                        exp->GetBasis(1)->GetBasisKey(),
                                        exp->GetBasis(2)->GetBasisKey());
                break;
            case LibUtilities::ePrism:
                stdExp = MemoryManager<StdRegions::StdPrismExp>
                    ::AllocateSharedPtr(exp->GetBasis(0)->GetBasisKey(),
                                        exp->GetBasis(1)->GetBasisKey(),
                                        exp->GetBasis(2)->GetBasisKey());
                break;
            case LibUtilities::eHexahedron:
                    stdExp = MemoryManager<StdRegions::StdHexExp>
                        ::AllocateSharedPtr(exp->GetBasis(0)->GetBasisKey(),
                                            exp->GetBasis(1)->GetBasisKey(),
                                            exp->GetBasis(2)->GetBasisKey());
                    break;
            default:
                ASSERTL0(false,"Shape type not setup");
                break;
            }

            return stdExp;
        }

        /**
         * @brief Construct collections of elements containing a single element
         * type and polynomial order from the list of expansions.
         */
        void ExpList::CreateCollections(Collections::ImplementationType ImpType)
        {
            map<LibUtilities::ShapeType,
                vector<std::pair<LocalRegions::ExpansionSharedPtr,int> > > collections;
            map<LibUtilities::ShapeType,
                vector<std::pair<LocalRegions::ExpansionSharedPtr,int> > >::iterator it;

            // Figure out optimisation parameters if provided in
            // session file or default given
            Collections::CollectionOptimisation colOpt(m_session, ImpType);
            ImpType = colOpt.GetDefaultImplementationType();

            bool autotuning = colOpt.IsUsingAutotuning();
            bool verbose    = (m_session->DefinesCmdLineArgument("verbose")) &&
                              (m_comm->GetRank() == 0);
            int  collmax    = (colOpt.GetMaxCollectionSize() > 0
                                        ? colOpt.GetMaxCollectionSize()
                                        : 2*m_exp->size());

            // clear vectors in case previously called
            m_collections.clear();
            m_coll_coeff_offset.clear();
            m_coll_phys_offset.clear();

            // Loop over expansions, and create collections for each element type
            for (int i = 0; i < m_exp->size(); ++i)
            {
                collections[(*m_exp)[i]->DetShapeType()].push_back(
                    std::pair<LocalRegions::ExpansionSharedPtr,int> ((*m_exp)[i],i));
            }

            for (it = collections.begin(); it != collections.end(); ++it)
            {
                LocalRegions::ExpansionSharedPtr exp = it->second[0].first;

                Collections::OperatorImpMap impTypes = colOpt.GetOperatorImpMap(exp);
                vector<StdRegions::StdExpansionSharedPtr> collExp;

                int prevCoeffOffset     = m_coeff_offset[it->second[0].second];
                int prevPhysOffset      = m_phys_offset [it->second[0].second];
                int collcnt;

                m_coll_coeff_offset.push_back(prevCoeffOffset);
                m_coll_phys_offset .push_back(prevPhysOffset);

                if(it->second.size() == 1) // single element case
                {
                    collExp.push_back(it->second[0].first);

                    // if no Imp Type provided and No settign in xml file.
                    // reset impTypes using timings
                    if(autotuning)
                    {
                        impTypes = colOpt.SetWithTimings(collExp,
                                                         impTypes, verbose);
                    }

                    Collections::Collection tmp(collExp, impTypes);
                    m_collections.push_back(tmp);
                }
                else
                {
                    // set up first geometry
                    collExp.push_back(it->second[0].first);
                    int prevnCoeff = it->second[0].first->GetNcoeffs();
                    int prevnPhys  = it->second[0].first->GetTotPoints();
                    collcnt = 1;

                    for (int i = 1; i < it->second.size(); ++i)
                    {
                        int nCoeffs     = it->second[i].first->GetNcoeffs();
                        int nPhys       = it->second[i].first->GetTotPoints();
                        int coeffOffset = m_coeff_offset[it->second[i].second];
                        int physOffset  = m_phys_offset [it->second[i].second];

                        // check to see if next elmt is different or
                        // collmax reached and if so end collection
                        // and start new one
                        if(prevCoeffOffset + nCoeffs != coeffOffset ||
                           prevnCoeff != nCoeffs ||
                           prevPhysOffset + nPhys != physOffset ||
                           prevnPhys != nPhys || collcnt >= collmax)
                        {

                            // if no Imp Type provided and No
                            // settign in xml file. reset
                            // impTypes using timings
                            if(autotuning)
                            {
                                impTypes = colOpt.SetWithTimings(collExp,
                                                                 impTypes,
                                                                 verbose);
                            }

                            Collections::Collection tmp(collExp, impTypes);
                            m_collections.push_back(tmp);


                            // start new geom list
                            collExp.clear();

                            m_coll_coeff_offset.push_back(coeffOffset);
                            m_coll_phys_offset .push_back(physOffset);
                            collExp.push_back(it->second[i].first);
                            collcnt = 1;
                        }
                        else // add to list of collections
                        {
                            collExp.push_back(it->second[i].first);
                            collcnt++;
                        }

                        // if end of list finish up collection
                        if (i == it->second.size() - 1)
                        {
                            // if no Imp Type provided and No
                            // settign in xml file.
                            if(autotuning)
                            {
                                impTypes = colOpt.SetWithTimings(collExp,
                                                                 impTypes,verbose);
                            }

                            Collections::Collection tmp(collExp, impTypes);
                            m_collections.push_back(tmp);
                            collExp.clear();
                            collcnt = 0;

                        }

                        prevCoeffOffset = coeffOffset;
                        prevPhysOffset  = physOffset;
                        prevnCoeff      = nCoeffs;
                        prevnPhys       = nPhys;
                    }
                }
            }
        }
    } //end of namespace
} //end of namespace

