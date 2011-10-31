///////////////////////////////////////////////////////////////////////////////
//
// File SegExp.h
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
// Description: Header file for SegExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef SEGEXP_H
#define SEGEXP_H

#include <LocalRegions/LocalRegions.hpp>

#include <StdRegions/StdSegExp.h>
#include <SpatialDomains/SegGeom.h>

#include <SpatialDomains/GeomFactors1D.h>

#include <LocalRegions/MatrixKey.h>

#include <LocalRegions/Expansion1D.h>
#include <LocalRegions/LocalRegionsDeclspec.h>

#include <fstream>

namespace Nektar
{
    namespace LocalRegions
    {

        class SegExp: virtual public StdRegions::StdSegExp, virtual public Expansion1D
        {

        public:
            /// Constructor using BasisKey class for quadrature points and
            /// order definition.
            LOCAL_REGIONS_EXPORT SegExp(const LibUtilities::BasisKey &Ba,
                   const SpatialDomains::Geometry1DSharedPtr &geom);

            /// Copy Constructor
            LOCAL_REGIONS_EXPORT SegExp(const SegExp &S);

            /// Destructor
            LOCAL_REGIONS_EXPORT ~SegExp();

            /// Returns the locations of the quadrature points in up to
            /// three-dimensions.
            LOCAL_REGIONS_EXPORT virtual void v_GetCoords(Array<OneD,NekDouble> &coords_1,
                           Array<OneD,NekDouble> &coords_2
                                        = NullNekDouble1DArray,
                           Array<OneD,NekDouble> &coords_3
                                        = NullNekDouble1DArray);

            ///
            LOCAL_REGIONS_EXPORT virtual void v_GetCoord(const Array<OneD, const NekDouble>& Lcoords,
                                Array<OneD,NekDouble> &coords);

            virtual int v_GetCoordim()
            {
                return m_geom->GetCoordim();
            }

            /// Returns a pointer to the GeomFactors object describing the
            /// metric information for the segment.
            virtual const SpatialDomains::GeomFactorsSharedPtr& v_GetMetricInfo() const
            {
                return m_metricinfo;
            }

            /// Returns a pointer to a Geometry object describing the
            /// geometry of the segment.
            virtual const SpatialDomains::GeometrySharedPtr v_GetGeom() const
            {
                return m_geom;
            }

            /// Returns a pointer to a Geometry1D object describing the
            /// geometry of the segment.
            virtual const SpatialDomains::Geometry1DSharedPtr& v_GetGeom1D() const
            {
                return m_geom;
            }

            /// Writes out the physical space data to file.
            LOCAL_REGIONS_EXPORT virtual void v_WriteToFile(std::ofstream &outfile, 
                            OutputFormat format,
                            const bool dumpVar = true,
                            std::string var = "v");


            //----------------------------
            // Integration Methods
            //----------------------------

            /// \brief Integrate the physical point list \a inarray over region
            LOCAL_REGIONS_EXPORT virtual NekDouble v_Integral(const Array<OneD, const NekDouble>& inarray);

            /** \brief  Inner product of \a inarray over region with respect to
                the expansion basis (this)->_Base[0] and return in \a outarray

                Wrapper call to SegExp::IProduct_WRT_B

                Input:\n

                - \a inarray: array of function evaluated at the physical
                collocation points

                Output:\n

                - \a outarray: array of inner product with respect to each
                basis over region
            */
            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray)
            {
                v_IProductWRTBase(m_base[0]->GetBdata(),inarray,outarray,1);
            }

            LOCAL_REGIONS_EXPORT virtual void v_IProductWRTDerivBase(const int dir,
                                      const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> & outarray);

            //-----------------------------
            // Differentiation Methods
            //-----------------------------


            /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
                physical quadrature points given by \a inarray and return in \a
                outarray. */
            LOCAL_REGIONS_EXPORT virtual void v_PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                           Array<OneD,NekDouble> &out_d0,
                           Array<OneD,NekDouble> &out_d1 = NullNekDouble1DArray,
                           Array<OneD,NekDouble> &out_d2 = NullNekDouble1DArray);

            LOCAL_REGIONS_EXPORT virtual void v_PhysDeriv(const int dir,
                           const Array<OneD, const NekDouble>& inarray,
                           Array<OneD, NekDouble> &outarray);

	   /**
	   *\brief Evaluate the derivative along a line: \f$ d/ds=\frac{spacedim}{||tangent||}d/d{\xi}  \f$.
	   * The derivative is calculated performing
	   *the product \f$ du/d{s}=\nabla u \cdot tangent \f$.
	   *\param inarray function to derive
	   *\param out_ds result of the derivative operation 
	   **/
           LOCAL_REGIONS_EXPORT virtual void v_PhysDeriv_s(const Array<OneD, const NekDouble>& inarray, 
            	    Array<OneD, NekDouble> &out_ds);

	   /**
	   *\brief Evaluate the derivative normal to a line: \f$ d/dn=\frac{spacedim}{||normal||}d/d{\xi}  \f$.
	   * The derivative is calculated performing
	   *the product \f$ du/d{s}=\nabla u \cdot normal \f$.
	   *\param inarray function to derive
	   *\param out_dn result of the derivative operation 
	   **/     
	   LOCAL_REGIONS_EXPORT virtual void v_PhysDeriv_n(const Array<OneD, const NekDouble>& inarray,
	   	    Array<OneD, NekDouble>& out_dn);
            //----------------------------
            // Evaluations Methods
            //---------------------------

            LOCAL_REGIONS_EXPORT virtual void v_SetCoeffsToOrientation(StdRegions::EdgeOrientation dir);

            LOCAL_REGIONS_EXPORT virtual void v_SetCoeffsToOrientation(StdRegions::EdgeOrientation dir,
                                        Array<OneD, const NekDouble> &inarray,
                                        Array<OneD, NekDouble> &outarray);



            /** \brief Forward transform from physical quadrature space
                stored in \a inarray and evaluate the expansion coefficients and
                store in \a (this)->_coeffs  */
            LOCAL_REGIONS_EXPORT virtual void v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                          Array<OneD,NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT virtual void v_FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD, NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& coord);

            LOCAL_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& coord, const Array<OneD, const NekDouble> & physvals);

            LOCAL_REGIONS_EXPORT virtual void v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT virtual void v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const double lambda);

        protected:

            virtual DNekMatSharedPtr v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey);


            /**
               \brief  Inner product of \a inarray over region with respect to
               expansion basis \a base and return in \a outarray

               Calculate \f$ I[p] = \int^{1}_{-1} \phi_p(\xi_1) u(\xi_1) d\xi_1
               = \sum_{i=0}^{nq-1} \phi_p(\xi_{1i}) u(\xi_{1i}) w_i \f$ where
               \f$ outarray[p] = I[p], inarray[i] = u(\xi_{1i}), base[p*nq+i] =
               \phi_p(\xi_{1i}) \f$.

               Inputs: \n

               - \a base: an array definiing the local basis for the inner
               product usually passed from Basis->get_bdata() or
               Basis->get_Dbdata()
               - \a inarray: physical point array of function to be integrated
               \f$ u(\xi_1) \f$
               - \a coll_check: Flag to identify when a Basis->collocation()
               call should be performed to see if this is a GLL_Lagrange basis
               with a collocation property. (should be set to 0 if taking the
               inner product with respect to the derivative of basis)

               Output: \n

               - \a outarray: array of coefficients representing the inner
               product of function with ever  mode in the exapnsion

            **/
            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble>& base,
                                 const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 int coll_check);

            SpatialDomains::Geometry1DSharedPtr m_geom;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;

            virtual void v_SetUpPhysTangents(const StdRegions::StdExpansionSharedPtr &exp2D, const int edge);

            virtual const Array<OneD, const NekDouble>& v_GetPhysNormals(void)
            {
                NEKERROR(ErrorUtil::efatal, "Got to SegExp");
                return NullNekDouble1DArray;
            }

            virtual SpatialDomains::GeomType v_MetricInfoType()
            {
                return m_metricinfo->GetGtype();
            }

            virtual void v_NormVectorIProductWRTBase(
                    const Array<OneD, const NekDouble>   &Fx,
                    const Array<OneD, const NekDouble> &Fy,
                    Array< OneD, NekDouble> &outarray,
                    bool NegateNormal = false);

            /** \brief Virtual function to evaluate the discrete \f$ L_\infty\f$
                error \f$ |\epsilon|_\infty = \max |u - u_{exact}|\f$ where \f$
                u_{exact}\f$ is given by the array \a sol.

                The full function is defined in StdExpansion::Linf

                Input:

                - \a _phys: takes the physical value space array as
                approximate solution

                - \a sol: array of solution function  at physical quadrature points

                output:

                - returns the \f$ L_\infty \f$ error as a NekDouble.
            */
            virtual NekDouble v_Linf(const Array<OneD, const NekDouble>& sol)
            {
                return Linf(sol);
            }

            /** \brief Virtual function to evaluate the \f$ L_\infty \f$ norm of
                the function defined at the physical points \a (this)->_phys.

                The full function is defined in StdExpansion::Linf

                Input:

                - \a _phys: uses the physical value space array as discrete
                function to be evaulated.

                output:

                - returns the \f$ L_\infty \f$  as a NekDouble.
            */
            virtual NekDouble v_Linf()
            {
                return Linf();
            }

            /** \brief Virtual function to evaluate the \f$ L_2\f$, \f$ |
                \epsilon |_{2} = \left [ \int^1_{-1} [u - u_{exact}]^2 dx
                \right]^{1/2} d\xi_1 \f$ where \f$ u_{exact}\f$ is given by the
                array sol.

                The full function is defined in StdExpansion::L2

                Input:

                - \a _phys: takes the physical value space array as
                approximate solution
                - \a sol: array of solution function  at physical quadrature points

                output:

                - returns the \f$ L_2 \f$ error as a NekDouble.
            */
            virtual NekDouble v_L2(const Array<OneD, const NekDouble>& sol)
            {
                return StdExpansion::L2(sol);
            }

            /** \brief Virtual function to evaluate the \f$ L_2\f$ norm of the
                function defined at the physical points \a (this)->_phys.

                The full function is defined in StdExpansion::L2

                Input:

                - \a _phys: uses the physical value space array as discrete
                function to be evaulated.

                output:

                - returns the \f$ L_2 \f$  as a NekDouble.
            */
            virtual NekDouble v_L2()
            {
                return StdExpansion::L2();
            }

            virtual DNekScalMatSharedPtr& v_GetLocMatrix(const MatrixKey &mkey)
            {
                return m_matrixManager[mkey];
            }


            virtual DNekScalMatSharedPtr& v_GetLocMatrix(const StdRegions::MatrixType mtype, NekDouble lambdaval, NekDouble tau)
            {
                MatrixKey mkey(mtype,DetExpansionType(),*this,lambdaval,tau);
                return m_matrixManager[mkey];
            }

            virtual DNekScalBlkMatSharedPtr& v_GetLocStaticCondMatrix(const MatrixKey &mkey)
            {
                return m_staticCondMatrixManager[mkey];
            }

            virtual void v_ExtractDataToCoeffs(const std::vector<NekDouble> &data, 
                                               const int offset, 
                                               const std::vector<unsigned int > &nummodes, 
                                               const int nmode_offset,
                                               Array<OneD, NekDouble> &coeffs);

            virtual DNekMatSharedPtr v_GenMatrix(const StdRegions::StdMatrixKey &mkey);


        private:
            SegExp();

            LOCAL_REGIONS_EXPORT void ReverseCoeffsAndSign(const Array<OneD,NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray);

            /** \brief Inverse Matrix Product */
            LOCAL_REGIONS_EXPORT void MultiplyByElmtInvMass(
                            const Array<OneD, const NekDouble>& inarray,
                                  Array<OneD,NekDouble> &outarray);

            DNekScalMatSharedPtr     CreateMatrix(const MatrixKey &mkey);

            DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(const MatrixKey &mkey);


            virtual int v_GetNumPoints(const int dir) const
            {
                return GetNumPoints(dir);
            }


            virtual int v_GetNcoeffs(void) const
            {
                return m_ncoeffs;
            }

            virtual const LibUtilities::BasisSharedPtr& v_GetBasis(int dir) const
            {
                return GetBasis(dir);
            }


            virtual int v_NumBndryCoeffs() const
            {
                return 2;
            }


            virtual int v_NumDGBndryCoeffs() const
            {
                return 2;
            }

        };

        // type defines for use of SegExp in a boost vector
        typedef boost::shared_ptr<SegExp>      SegExpSharedPtr;
        typedef std::vector< SegExpSharedPtr > SegExpVector;
        typedef std::vector< SegExpSharedPtr >::iterator SegExpVectorIter;
    } //end of namespace
} //end of namespace

#endif // SEGEXP_H

