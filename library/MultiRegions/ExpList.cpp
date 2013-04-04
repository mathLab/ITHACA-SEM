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

#include <LocalRegions/MatrixKey.h>     // for MatrixKey
#include <LocalRegions/Expansion.h>     // for Expansion

#include <MultiRegions/AssemblyMap/AssemblyMapCG.h>  // for AssemblyMapCG, etc
#include <MultiRegions/GlobalLinSysKey.h>  // for GlobalLinSysKey
#include <MultiRegions/GlobalMatrix.h>  // for GlobalMatrix, etc
#include <MultiRegions/GlobalMatrixKey.h>  // for GlobalMatrixKey

#include <LibUtilities/LinearAlgebra/NekTypeDefs.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>


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
            m_exp(MemoryManager<StdRegions::StdExpansionVector>
                      ::AllocateSharedPtr()),
            m_coeff_offset(),
            m_phys_offset(),
            m_offset_elmt_id(),
            m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
            m_WaveSpace(false)
        {
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
            m_exp(MemoryManager<StdRegions::StdExpansionVector>
                      ::AllocateSharedPtr()),
            m_coeff_offset(),
            m_phys_offset(),
            m_offset_elmt_id(),
            m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
            m_WaveSpace(false)
        {
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
            m_exp(MemoryManager<StdRegions::StdExpansionVector>
                      ::AllocateSharedPtr()),
            m_coeff_offset(),
            m_phys_offset(),
            m_offset_elmt_id(),
            m_blockMat(MemoryManager<BlockMatrixMap>::AllocateSharedPtr()),
            m_WaveSpace(false)
        {
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
            m_coeff_offset(in.m_coeff_offset),
            m_phys_offset(in.m_phys_offset),
            m_offset_elmt_id(in.m_offset_elmt_id),
            m_globalOptParam(in.m_globalOptParam),
            m_blockMat(in.m_blockMat),
            m_WaveSpace(false)
        {
            if(DeclareCoeffPhysArrays)
            {
                m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
                m_phys   = Array<OneD, NekDouble>(m_npoints);
            }
        }
        
        
        //boost::shared_ptr<ExpList> do_clone(void) const = 0; {}
		
        /**
         * For each element, copy the coefficients from \a m_coeffs into their
         * respective element expansion from \a m_exp.
         */
        void ExpList::PutCoeffsInToElmtExp()
        {
            int i;
            int order_e;

            for(i = 0; i < (*m_exp).size(); ++i)
            {
                order_e = (*m_exp)[i]->GetNcoeffs();
                Vmath::Vcopy(order_e,&m_coeffs[m_coeff_offset[i]], 1,
                                         &((*m_exp)[i]->UpdateCoeffs())[0],1);
            }
        }


        /**
         * Copy the coefficients associated with element \a eid from \a m_coeffs
         * to the corresponding element expansion object from \a m_exp.
         * @param   eid         Index of element for which copy is performed.
         */
        void ExpList::PutCoeffsInToElmtExp(int eid)
        {
            int order_e;
            int cnt = 0;

            order_e = (*m_exp)[eid]->GetNcoeffs();
            cnt = m_coeff_offset[eid];
            Vmath::Vcopy(order_e,&m_coeffs[cnt], 1,
                         &((*m_exp)[eid]->UpdateCoeffs())[0],1);
        }


        /**
         * Coefficients from each local expansion are copied into the
         * concatenated list of coefficients for all elements.
         */
        void ExpList::PutElmtExpInToCoeffs(void)
        {
            int i;
            int order_e;

            for(i = 0; i < (*m_exp).size(); ++i)
            {
                order_e = (*m_exp)[i]->GetNcoeffs();
                Vmath::Vcopy(order_e, &((*m_exp)[i]->UpdateCoeffs())[0],1,
                             &m_coeffs[m_coeff_offset[i]],1);
            }
        }


        /**
         * Coefficients for a single element are copied from the associated
         * element expansion to the concatenated coefficient array.
         * @param   eid         Index of element to copy.
         */
        void ExpList::PutElmtExpInToCoeffs(int eid)
        {
            int order_e;
            int cnt = 0;

            order_e = (*m_exp)[eid]->GetNcoeffs();
            cnt = m_coeff_offset[eid];

            Vmath::Vcopy(order_e, &((*m_exp)[eid]->UpdateCoeffs())[0],1,
                             &m_coeffs[cnt],1);
        }


        /**
         * The local expansion objects are populated with the physical
         * evaluation at the quadrature points stored in the \a m_phys storage.
         */
        void ExpList::PutPhysInToElmtExp()
        {
            PutPhysInToElmtExp(m_phys);
        }


        /**
         * The local expansion objects are populated with the supplied physical
         * evaluations at the quadrature points. The layout and order of the
         * supplied data is assumed to conform to the expansion list.
         * @param   in          Physical quadrature data.
         */
        void ExpList::PutPhysInToElmtExp(Array<OneD,const NekDouble> &in)
        {
            int i;
            int npoints_e;

            for(i = 0; i < (*m_exp).size(); ++i)
            {
                npoints_e = (*m_exp)[i]->GetTotPoints();
                Vmath::Vcopy(npoints_e, &in[m_phys_offset[i]],1,
                                        &((*m_exp)[i]->UpdatePhys())[0],1);
            }
        }


        /**
         * The physical evaluations at the quadrature points from the expansion
         * objects are concatenated and stored in \a out.
         * @param   out         Storage for physical values.
         */
        void ExpList::PutElmtExpInToPhys(Array<OneD,NekDouble> &out)
        {
            int i;
            int npoints_e;

            for(i = 0; i < (*m_exp).size(); ++i)
            {
                npoints_e = (*m_exp)[i]->GetTotPoints();
                Vmath::Vcopy(npoints_e, &((*m_exp)[i]->GetPhys())[0],1,
                             &out[m_phys_offset[i]],1);
            }
        }


        /**
         * The physical evaluations at the quadrature points in the element
         * expansion \a eid are copied to \a out.
         * @param   out         Storage for physical values.
         */
        void ExpList::PutElmtExpInToPhys(int eid, Array<OneD,NekDouble> &out)
        {
            int npoints_e;
            int cnt = m_phys_offset[eid];

            npoints_e = (*m_exp)[eid]->GetTotPoints();
            Vmath::Vcopy(npoints_e, &((*m_exp)[eid]->GetPhys())[0],1,
                         &out[cnt],1);
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
            // get optimisation information about performing block
            // matrix multiplies
            const Array<OneD, const bool>  doBlockMatOp
                = m_globalOptParam->DoBlockMatOp(StdRegions::eIProductWRTBase);
            const Array<OneD, LibUtilities::ShapeType> shape = m_globalOptParam->GetShapeList();
            const Array<OneD, const int> num_elmts = m_globalOptParam->GetShapeNumElements();

            Array<OneD,NekDouble> tmp_outarray;
            int cnt = 0,eid;

            for(int n = 0; n < shape.num_elements(); ++n)
            {
                if(doBlockMatOp[n])
                {
                    if(num_elmts[n])
                    {
                        GlobalMatrixKey mkey(StdRegions::eIProductWRTBase,
                                             shape[n]);
                        eid = m_offset_elmt_id[cnt];
                        MultiplyByBlockMatrix(mkey,inarray + m_phys_offset[eid],
                                              tmp_outarray = outarray + m_coeff_offset[eid]);
                        cnt += num_elmts[n];
                    }
                }
                else
                {
                    int    i;
                    for(i = 0; i < num_elmts[n]; ++i)
                    {
                        eid = m_offset_elmt_id[cnt++];
                        (*m_exp)[eid]->IProductWRTBase(inarray+m_phys_offset[eid],
                                                       tmp_outarray = outarray+m_coeff_offset[eid]);
                    }
                }
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
            int  i;
            Array<OneD, NekDouble> e_out_d0;
            Array<OneD, NekDouble> e_out_d1;
            Array<OneD, NekDouble> e_out_d2;

            for(i= 0; i < (*m_exp).size(); ++i)
            {
                e_out_d0 = out_d0 + m_phys_offset[i];
                if(out_d1.num_elements())
                {
                    e_out_d1 = out_d1 + m_phys_offset[i];
                }

                if(out_d2.num_elements())
                {
                    e_out_d2 = out_d2 + m_phys_offset[i];
                }
                (*m_exp)[i]->PhysDeriv(inarray+m_phys_offset[i],e_out_d0,e_out_d1,e_out_d2);
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

            map< pair< int,  int>, NekDouble > spcoomat;
            pair<int,int> coord;

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
                ::AllocateSharedPtr(glob_rows,glob_cols,spcoomat);
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
            // get optimisation information about performing block
            // matrix multiplies
            const Array<OneD, const bool>  doBlockMatOp
                = m_globalOptParam->DoBlockMatOp(StdRegions::eBwdTrans);
            const Array<OneD, LibUtilities::ShapeType> shape = m_globalOptParam->GetShapeList();
            const Array<OneD, const int> num_elmts = m_globalOptParam->GetShapeNumElements();

            Array<OneD,NekDouble> tmp_outarray;
            int cnt = 0,eid;

            for(int n = 0; n < num_elmts.num_elements(); ++n)
            {
                if(doBlockMatOp[n])
                {
                    if(num_elmts[n])
                    {
                        GlobalMatrixKey mkey(StdRegions::eBwdTrans, shape[n]);
                        eid = m_offset_elmt_id[cnt];
                        MultiplyByBlockMatrix(mkey,inarray + m_coeff_offset[eid],
                                              tmp_outarray = outarray + m_phys_offset[eid]);
                        cnt += num_elmts[n];
                    }
                }
                else
                {
                    int  i;

                    for(i= 0; i < num_elmts[n]; ++i)
                    {
                        eid = m_offset_elmt_id[cnt++];
                        (*m_exp)[eid]->BwdTrans(inarray + m_coeff_offset[eid],
                                   tmp_outarray = outarray+m_phys_offset[eid]);
                    }
                }
            }
        }

        StdRegions::StdExpansionSharedPtr& ExpList::GetExp(
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
                                 NekDouble tol)
        {
            static int start = 0;
            // start search at previous element or 0 
            for (int i = start; i < (*m_exp).size(); ++i)
            {
                if ((*m_exp)[i]->GetGeom()->ContainsPoint(gloCoord,tol))
                {
                    start = i;
                    return i;
                }
            }

            for (int i = 0; i < start; ++i)
            {
                if ((*m_exp)[i]->GetGeom()->ContainsPoint(gloCoord,tol))
                {
                    start = i;
                    return i;
                }
            }
            return -1;
        }


        /**
         * The operation is evaluated locally by the elemental
         * function StdRegions#StdExpansion#GetSurfaceNormal.
         */
        void ExpList::GetSurfaceNormal(Array<OneD, NekDouble> &SurfaceNormal,
                                const int k)
        {
            int i;
            Array<OneD, NekDouble> normals;

            for(i = 0; i < (*m_exp).size(); ++i)
            {
                //e_SN = SurfaceNormal + m_phys_offset[i];
                normals = (*m_exp)[i]->GetSurfaceNormal()[k];
                Vmath::Vcopy(normals.num_elements(), &normals[0], 1, &SurfaceNormal[0] + m_phys_offset[i], 1);
            }
        }

        void ExpList::GetTangents(
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &tangents)
        {
            int i,j,k,e_npoints,offset;
            Array<OneD,Array<OneD, NekDouble> > loctangent;

            // Assume whole array is of same coordinate dimension
            int coordim = (*m_exp)[0]->GetGeom()->GetCoordim();

            ASSERTL0(tangents.num_elements() > 0,
                     "Must have storage for at least one tangent");
            ASSERTL1(tangents[0].num_elements() >= coordim,
                     "Output vector does not have sufficient dimensions to "
                     "match coordim");

            // Process each expansion.
            for(i = 0; i < m_exp->size(); ++i)
            {
                // Get the number of points and normals for this expansion.
                e_npoints  = (*m_exp)[i]->GetTotPoints();
                offset = m_phys_offset[i];

                for (j = 0; j < tangents.num_elements(); ++j)
                {
                    loctangent = (*m_exp)[i]->GetMetricInfo()->GetTangent(j);
                    // Get the physical data offset for this expansion.

                    for (k = 0; k < coordim; ++k)
                    {
                        Vmath::Vcopy(e_npoints, &(loctangent[k][0]), 1,
                                                &(tangents[j][k][offset]), 1);
                    }
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
            std::string dir = "TangentX";
            Array<OneD,NekDouble> coords(2);

            m_session->LoadGeometricInfo("TangentDir",dir,"TangentX");
            m_session->LoadGeometricInfo("TangentCentreX",coords[0],0.0);
            m_session->LoadGeometricInfo("TangentCentreY",coords[1],0.0);

            // Apply geometric info to each expansion.
            for (int i = 0; i < m_exp->size(); ++i)
            {
                (*m_exp)[i]->GetMetricInfo()->SetTangentOrientation(dir);
                (*m_exp)[i]->GetMetricInfo()->SetTangentCircularCentre(coords);
            }
        }


        /**
         * The coordinates of the quadrature points, together with
         * the content of the array #m_phys, are written to the
         * file \a out.
         *
         * @param   out             The file to which the solution should be
         *                          written.
         */
        void ExpList::WriteToFile(std::ofstream &out, OutputFormat format,
                                std::string var)
        {
            if(format==eTecplot)
            {
                int i;

                Array<OneD, const NekDouble> phys = m_phys;

                if(m_physState == false)
                {
                    BwdTrans(m_coeffs,m_phys);
                }

                (*m_exp)[0]->SetPhys(phys+m_phys_offset[0]);
                (*m_exp)[0]->WriteToFile(out,eTecplot,true,var);

                for(i= 1; i < (*m_exp).size(); ++i)
                {
                    (*m_exp)[i]->SetPhys(phys+m_phys_offset[i]);
                    (*m_exp)[i]->WriteToFile(out,eTecplot,false,var);
                }
            }
            else if(format==eGnuplot)
            {
                int i;

                Array<OneD, const NekDouble> phys = m_phys;

                if(m_physState == false)
                {
                    BwdTrans(m_coeffs,m_phys);
                }

                (*m_exp)[0]->SetPhys(phys+m_phys_offset[0]);
                (*m_exp)[0]->WriteToFile(out,eGnuplot,true,var);

                for(i= 1; i < (*m_exp).size(); ++i)
                {
                    (*m_exp)[i]->SetPhys(phys+m_phys_offset[i]);
                    (*m_exp)[i]->WriteToFile(out,eTecplot,false,var);
                }
            }
            else if(format==eGmsh)
            {

                out<<"View.MaxRecursionLevel = 4;"<<endl;
                out<<"View.TargetError = 0.00;"<<endl;
                out<<"View.AdaptVisualizationGrid = 1;"<<endl;

                int i,j,k;
                int nElementalCoeffs =  (*m_exp)[0]->GetBasisNumModes(0);
                int nDumpCoeffs =  nElementalCoeffs*nElementalCoeffs;
                Array<TwoD, int> exponentMap(nDumpCoeffs,3,0);
                int cnt = 0;
                for(i = 0; i < nElementalCoeffs; i++)
                {
                    for(j = 0; j < nElementalCoeffs; j++)
                    {
                        exponentMap[cnt][0] = j;
                        exponentMap[cnt++][1] = i;
                    }
                }

                PutCoeffsInToElmtExp();
                bool dumpNewView = true;
                bool closeView = false;
                for(i= 0; i < (*m_exp).size(); ++i)
                {
                    if(nElementalCoeffs != (*m_exp)[i]->GetBasisNumModes(0))
                    {
                        ASSERTL0(false,"Not all elements have the same number "
                                       "of expansions, this will probably lead "
                                       "to a corrupt Gmsh-output file.")
                    }

                    if(i>0)
                    {
                        if ( ((*m_exp)[i]->DetShapeType())
                                        !=((*m_exp)[i-1]->DetShapeType()) )
                        {
                            dumpNewView = true;
                        }
                        else
                        {
                            dumpNewView = false;
                        }
                    }
                    if(i<(*m_exp).size()-1)
                    {
                        if ( ((*m_exp)[i]->DetShapeType())
                                        !=((*m_exp)[i+1]->DetShapeType()) )
                        {
                            closeView = true;
                        }
                        else
                        {
                            closeView = false;
                        }
                    }
                    else
                    {
                            closeView = true;
                    }

                    if(dumpNewView)
                    {
                        out<<"View \" \" {"<<endl;
                    }

                    (*m_exp)[i]->WriteToFile(out,eGmsh,false);

                    if(closeView)
                    {
                        out<<"INTERPOLATION_SCHEME"<<endl;
                        out<<"{"<<endl;
                        for(k=0; k < nDumpCoeffs; k++)
                        {
                            out<<"{";
                            for(j = 0; j < nDumpCoeffs; j++)
                            {
                                if(k==j)
                                {
                                    out<<"1.00";
                                }
                                else
                                {
                                    out<<"0.00";
                                }
                                if(j < nDumpCoeffs - 1)
                                {
                                    out<<", ";
                                }
                            }
                            if(k < nDumpCoeffs - 1)
                            {
                                out<<"},"<<endl;
                            }
                            else
                            {
                                out<<"}"<<endl<<"}"<<endl;
                            }
                        }

                        out<<"{"<<endl;
                        for(k=0; k < nDumpCoeffs; k++)
                        {
                            out<<"{";
                            for(j = 0; j < 3; j++)
                            {
                                out<<exponentMap[k][j];
                                if(j < 2)
                                {
                                    out<<", ";
                                }
                            }
                            if(k < nDumpCoeffs - 1)
                            {
                                out<<"},"<<endl;
                            }
                            else
                            {
                                out<<"}"<<endl<<"};"<<endl;
                            }
                        }
                        out<<"};"<<endl;
                    }
                }
                out<<"Combine ElementsFromAllViews;"<<endl;
                out<<"View.Name = \"\";"<<endl;
            }
            else
            {
                ASSERTL0(false, "Output routine not implemented for requested "
                                "type of output");
            }
        }


        /**
         * Write Tecplot Files Header
         * @param   outfile Output file name.
         * @param   var                 variables names
         */
        void ExpList::v_WriteTecplotHeader(std::ofstream &outfile,
                        std::string var)
        {

            int coordim  = GetExp(0)->GetCoordim();
            outfile << "Variables = x";

            if(coordim == 2)
            {
                outfile << ", y";
            }
            else if (coordim == 3)
            {
                outfile << ", y, z";
            }
            outfile << ", "<< var << std::endl << std::endl;
        }


        /**
         * Write Tecplot Files Zone
         * @param   outfile    Output file name.
         * @param   expansion  Expansion that is considered
         */
        void ExpList::v_WriteTecplotZone(std::ofstream &outfile, int expansion)
        {
            (*m_exp)[expansion]->WriteTecplotZone(outfile);
        }

        /**
         * Write Tecplot Files Field
         * @param   outfile    Output file name.
         * @param   expansion  Expansion that is considered
         */
        void ExpList::v_WriteTecplotField(std::ofstream &outfile, int expansion)
        {
            (*m_exp)[expansion]->SetPhys(m_phys+m_phys_offset[expansion]);
            (*m_exp)[expansion]->WriteTecplotField(outfile);
        }


        void ExpList::WriteVtkHeader(std::ofstream &outfile)
        {
            outfile << "<?xml version=\"1.0\"?>" << endl;
            outfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
                    << "byte_order=\"LittleEndian\">" << endl;
            outfile << "  <UnstructuredGrid>" << endl;
        }

        void ExpList::WriteVtkFooter(std::ofstream &outfile)
        {
            outfile << "  </UnstructuredGrid>" << endl;
            outfile << "</VTKFile>" << endl;
        }

        void ExpList::v_WriteVtkPieceHeader(std::ofstream &outfile, int expansion)
        {
            ASSERTL0(false, "Routine not implemented for this expansion.");
        }

        void ExpList::WriteVtkPieceFooter(std::ofstream &outfile, int expansion)
        {
            outfile << "      </PointData>" << endl;
            outfile << "    </Piece>" << endl;
        }

        void ExpList::v_WriteVtkPieceData(std::ofstream &outfile, int expansion,
                                        std::string var)
        {
            int i;
            int nq = (*m_exp)[expansion]->GetTotPoints();

            // printing the fields of that zone
            outfile << "        <DataArray type=\"Float32\" Name=\""
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

        void ExpList::ReadFromFile(std::ifstream &in, OutputFormat format)
        {
            if(format==eTecplot)
            {
                int i,npts;
                Array<OneD, NekDouble> phys = m_phys;

                npts = (*m_exp)[0]->GetTotPoints();
                (*m_exp)[0]->ReadFromFile(in,eTecplot,true);
                Vmath::Vcopy(npts,&(*m_exp)[0]->GetPhys()[0],1,&phys[m_phys_offset[0]],1);

                for(i= 1; i < (*m_exp).size(); ++i)
                {
                    npts = (*m_exp)[i]->GetTotPoints();
                    (*m_exp)[i]->ReadFromFile(in,eTecplot,false);
                    Vmath::Vcopy(npts,&((*m_exp)[i]->GetPhys())[0],1,
                                 &phys[m_phys_offset[i]],1);
                }
                FwdTrans(m_phys,m_coeffs);
            }
            else
            {
                ASSERTL0(false, "Output routine not implemented for requested "
                                "type of output");
            }
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
        NekDouble  ExpList::Linf(const Array<OneD, const NekDouble> &soln)
        {
            NekDouble err = 0.0;
            int       i;

            for(i= 0; i < (*m_exp).size(); ++i)
            {
                // set up physical solution in local element
                (*m_exp)[i]->SetPhys(m_phys+m_phys_offset[i]);
                err  = std::max(err,(*m_exp)[i]->Linf(soln + m_phys_offset[i]));
            }
            m_comm->GetRowComm()->AllReduce(err, LibUtilities::ReduceMax);

            return err;
        }

        /**
         * Given a spectral/hp approximation
         * \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
         * quadrature points (which should be contained in #m_phys),
         * this function calculates the \f$L_\infty\f$ error of this
         * approximation. The local distribution of the quadrature
         * points allows an elemental evaluation of this operation
         * through the functions StdRegions#StdExpansion#Linf.
         *s
         * The exact solution, also evaluated at the quadrature
         * points, should be contained in the variable #m_phys of
         * the ExpList object \a Sol.
         *
         * @return  The \f$L_\infty\f$ error of the approximation.
         */
        NekDouble  ExpList::Linf(void)
        {
            NekDouble err = 0.0;
            int       i;

            for(i= 0; i < (*m_exp).size(); ++i)
            {
                // set up physical solution in local element
                (*m_exp)[i]->SetPhys(m_phys+m_phys_offset[i]);
                err  = std::max(err,(*m_exp)[i]->Linf());
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
        NekDouble ExpList::v_L2(const Array<OneD, const NekDouble> &soln)
        {
            NekDouble err = 0.0,errl2;
            int    i;

            for(i= 0; i < (*m_exp).size(); ++i)
            {				
                // set up physical solution in local element
                (*m_exp)[i]->SetPhys(m_phys+m_phys_offset[i]);
                errl2 = (*m_exp)[i]->L2(soln+m_phys_offset[i]);
                err += errl2*errl2;
            }
            m_comm->GetRowComm()->AllReduce(err, LibUtilities::ReduceSum);

            return sqrt(err);
        }


        /**
         * Given a spectral/hp approximation
         * \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
         * quadrature points (which should be contained in #m_phys),
         * this function calculates the \f$L_2\f$ measure of this
         * approximation. The local distribution of the quadrature
         * points allows an elemental evaluation of this operation
         * through the functions StdRegions#StdExpansion#L2.
         *
         * The exact solution, also evaluated at the quadrature points, should
         * be contained in the variable #m_phys of the ExpList object \a Sol.
         *
         * @param   soln            A 1D array, containing the discrete
         *                          evaluation of the exact solution at the
         *                          quadrature points.
         * @return  The \f$L_2\f$ error of the approximation.
         */
        NekDouble ExpList::v_L2(void)
        {
            NekDouble err = 0.0,errl2;
            int    i;

            for(i= 0; i < (*m_exp).size(); ++i)
            {
                // set up physical solution in local element
                (*m_exp)[i]->SetPhys(m_phys+m_phys_offset[i]);
                errl2 = (*m_exp)[i]->L2();
                err += errl2*errl2;
            }
            m_comm->GetRowComm()->AllReduce(err, LibUtilities::ReduceSum);
            
            return sqrt(err);
        }
		

        NekDouble ExpList::v_Integral(const Array<OneD, const NekDouble> &inarray)
        {
            NekDouble err = 0.0;
            int       i   = 0;

            for(i = 0; i < (*m_exp).size(); ++i)
            {
                err += (*m_exp)[m_offset_elmt_id[i]]->Integral(inarray+m_phys_offset[i]);
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
        NekDouble ExpList::H1(const Array<OneD, const NekDouble> &soln)
        {

            NekDouble err = 0.0,errh1;
            int    i;

            for(i= 0; i < (*m_exp).size(); ++i)
            {
                // set up physical solution in local element
                (*m_exp)[i]->SetPhys(m_phys+m_phys_offset[i]);
                errh1 =  (*m_exp)[i]->H1(soln+m_phys_offset[i]);
                err  += errh1*errh1;
            }
            m_comm->AllReduce(err, LibUtilities::ReduceSum);
            return sqrt(err);
        }

        void  ExpList::GeneralGetFieldDefinitions(std::vector<LibUtilities::FieldDefinitionsSharedPtr> &fielddef, 
												  int NumHomoDir, 
												  Array<OneD, LibUtilities::BasisSharedPtr> &HomoBasis, 
												  std::vector<NekDouble> &HomoLen,
												  std::vector<unsigned int> &HomoZIDs,
												  std::vector<unsigned int> &HomoYIDs)
        {
            int startenum, endenum, s;

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
                LibUtilities::ShapeType               shape;
                std::vector<unsigned int>             elementIDs;
                std::vector<LibUtilities::BasisType>  basis;
                std::vector<unsigned int>             numModes;
                std::vector<std::string>              fields;

                bool first    = true;
                bool UniOrder = true;
                int n;

                for(int i = 0; i < (*m_exp).size(); ++i)
                {
                    if((*m_exp)[i]->GetGeom()->GetShapeType() == (LibUtilities::ShapeType) s)
                    {
                        elementIDs.push_back((*m_exp)[i]->GetGeom()->GetGlobalID());
                        if(first)
                        {
                            shape = (LibUtilities::ShapeType) s;
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
                    LibUtilities::FieldDefinitionsSharedPtr fdef  = MemoryManager<LibUtilities::FieldDefinitions>::AllocateSharedPtr(shape, elementIDs, basis, UniOrder, numModes,fields, NumHomoDir, HomoLen, HomoZIDs, HomoYIDs);
                    fielddef.push_back(fdef);
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
            int i, cnt, expId;
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

            // If no session is set, we use the non-parallel version of this
            // routine. This is used for reading BCs from files - note therefore
            // that this will probably not work in parallel.
            if (!m_session)
            {
                for (cnt = i = 0; i < fielddef->m_elementIDs.size(); ++i)
                {
                    const int elmtId = fielddef->m_elementIDs[i];
                    if (elmtToExpId.count(elmtId) == 0)
                    {
                        continue;
                    }
                    
                    expId   = elmtToExpId[elmtId];
                    datalen = (*m_exp)[expId]->CalcNumberOfCoefficients(
                        fielddef->m_numModes, modes_offset);

                    // Reset modes_offset in the case where all expansions of
                    // the same order.
                    if (fielddef->m_uniOrder == true)
                    {
                        modes_offset = 0;
                    }

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
                }

                return;
            }

            // Determine rank and number of processors.
            LibUtilities::CommSharedPtr vComm =
                m_session->GetComm()->GetRowComm();
            int n = vComm->GetSize();
            int p = vComm->GetRank();

            // Determine number of elements in fielddef located on this process.
            for(cnt = i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                if (elmtToExpId.count(fielddef->m_elementIDs[i]) == 0)
                {
                    continue;
                }
                ++cnt;
            }

            // Exchange this information between processors.
            Array<OneD, int> numEls(n, 0);
            numEls[p] = cnt;
            vComm->AllReduce(numEls, LibUtilities::ReduceSum);
            int totEls = Vmath::Vsum(n, numEls, 1);
                
            Array<OneD, int> elOffsets(n, 0);
            elOffsets[0] = 0;
            for (i = 1; i < n; ++i)
            {
                elOffsets[i] = elOffsets[i-1] + numEls[i-1];
            }
            
            // Storage holding number of coefficients per element and their
            // global IDs.
            Array<OneD, int> coeffsPerEl  (totEls, 0);
            Array<OneD, int> elmtGlobalIds(totEls, 0);
                
            // Determine number of coefficients in each local (to this
            // partition) element and store in the arrays above.
            for(cnt = i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                const int elmtId = fielddef->m_elementIDs[i];
                
                if (elmtToExpId.count(elmtId) == 0)
                {
                    continue;
                }

                expId   = elmtToExpId[elmtId];
                datalen = (*m_exp)[expId]->CalcNumberOfCoefficients(
                    fielddef->m_numModes, modes_offset);

                if(fielddef->m_uniOrder == true)
                {
                    modes_offset = 0;
                }

                elmtGlobalIds[cnt + elOffsets[p]] = fielddef->m_elementIDs[i];
                coeffsPerEl  [cnt + elOffsets[p]] = datalen;
                cnt++;
            }

            // Exchange this information so that each processor knows about all
            // coefficients per element and the corresponding global ID.
            vComm->AllReduce(coeffsPerEl,   LibUtilities::ReduceSum);
            vComm->AllReduce(elmtGlobalIds, LibUtilities::ReduceSum);

            // Map taking element global ID to number of coefficients for that
            // element.
            map<int,int> coeffsElmtMap;

            for (i = 0; i < totEls; ++i)
            {
                // Ensure that global IDs are mapped precisely once.
                ASSERTL0(coeffsElmtMap.count(elmtGlobalIds[i]) == 0,
                         "Error in communicating global ids for field "+
                         field + "!");
                coeffsElmtMap[elmtGlobalIds[i]] = coeffsPerEl[i];
            }
                
            for (i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                const int elmtId = fielddef->m_elementIDs[i];

                if (elmtToExpId.count(elmtId) == 0)
                {
                    ASSERTL1(coeffsElmtMap.count(elmtId) == 1,
                             "Couldn't find element!");
                    offset += coeffsElmtMap[elmtId];
                    continue;
                }
                    
                expId   = elmtToExpId  [elmtId];
                datalen = coeffsElmtMap[elmtId];

                if(fielddef->m_uniOrder == true)
                {
                    modes_offset = 0;
                }

                if(datalen == (*m_exp)[expId]->GetNcoeffs())
                {
                    // Copy data if it is the same length as expansion.
                    Vmath::Vcopy(datalen, &fielddata[offset], 1, 
                                 &coeffs[m_coeff_offset[expId]], 1);
                }
                else
                {
                    // unpack data to new order
                    (*m_exp)[expId]->ExtractDataToCoeffs(
                        &fielddata[offset], fielddef->m_numModes,
                        modes_offset, &coeffs[m_coeff_offset[expId]]);
                }
                    
                offset += datalen;
            }                
        }

        void ExpList::v_ExtractCoeffsToCoeffs(const boost::shared_ptr<ExpList> &fromExpList, const Array<OneD, const NekDouble> &fromCoeffs, Array<OneD, NekDouble> &toCoeffs)
        {     	
            int i;
            int offset = 0;

            // check if the same and if so just copy over coeffs
            if(fromExpList->GetNcoeffs() == m_ncoeffs)
            {
                Vmath::Vcopy(m_ncoeffs,fromCoeffs,1,toCoeffs,1);
            }
            else
            {
                std::vector<unsigned int> nummodes;
                for(i = 0; i < (*m_exp).size(); ++i)
                {
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

        boost::shared_ptr<ExpList> &ExpList::v_GetTrace(int i)
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
	
        void ExpList::v_ImposeDirichletConditions(Array<OneD,NekDouble>& outarray)
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
            v_IProductWRTBase_IterPerExp(inarray,outarray);
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
        void ExpList::v_GetCoords(NekDouble &x,NekDouble &y,NekDouble &z)
        {
			ASSERTL0(false,
                     "This method is not defined or valid for this class type");
		}
		
		/**
         */
        void ExpList::v_GetCoord(Array<OneD, NekDouble> &coords)
        {
			ASSERTL0(false,
                     "This method is not defined or valid for this class type");
		}
		
		/**
         */
		void ExpList::v_SetCoeff(NekDouble val)
        {
			ASSERTL0(false,
                     "This method is not defined or valid for this class type");
		}
		
		/**
         */
		void ExpList::v_SetPhys(NekDouble val)
        {
			ASSERTL0(false,
                     "This method is not defined or valid for this class type");
		}
		
		/**
         */
		const SpatialDomains::VertexComponentSharedPtr &ExpList::v_GetGeom(void) const
		{
			ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            static SpatialDomains::VertexComponentSharedPtr result;
            return result;
		}
		
		/**
         */
		const SpatialDomains::VertexComponentSharedPtr &ExpList::v_GetVertex(void) const
		{
			ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            static SpatialDomains::VertexComponentSharedPtr result;
            return result;
		}
		
		/**
         */
        void ExpList::v_SetUpPhysNormals()
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        void ExpList::v_SetUpPhysTangents(
                    const StdRegions::StdExpansionVector &locexp)
        {
            ASSERTL0(false,
                      "This method is not defined or valid for this class type");
        }

        /**
         */
        void ExpList::v_SetUpTangents()
        {
            for (int i = 0; i < (*m_exp).size(); ++i) {
                (*m_exp)[i]->GetMetricInfo()->SetUpTangents();
            }
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
            static Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
                                                                        result;
            return result;
        }

		/**
         */
        Array<OneD,SpatialDomains::BoundaryConditionShPtr> &ExpList::v_UpdateBndConditions()
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
            static Array<OneD,SpatialDomains::BoundaryConditionShPtr>
                                                                        result;
            return result;
        }

		/**
         */
        void ExpList::v_EvaluateBoundaryConditions(const NekDouble time, const NekDouble x2_in, const NekDouble x3_in)
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
        void ExpList::v_GetPeriodicEdges(
            vector<map<int,int> > &periodicVertices,
            map<int,int>          &periodicEdges)
        {
            ASSERTL0(false,
                     "This method is not defined or valid for this class type");
        }

        SpatialDomains::BoundaryConditionShPtr ExpList::GetBoundaryCondition(const SpatialDomains::BoundaryConditionCollection& collection,
                                                                             unsigned int index, const std::string& variable)
        {
            SpatialDomains::BoundaryConditionCollection::const_iterator collectionIter = collection.find(index);
            ASSERTL1(collectionIter != collection.end(), "Unable to locate collection.");
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
    } //end of namespace
} //end of namespace

