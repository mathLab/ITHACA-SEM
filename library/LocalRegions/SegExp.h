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

        class SegExp: public StdRegions::StdSegExp, public Expansion1D
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

            /// Return Shape of region, using  ShapeType enum list. i.e. Segment
            StdRegions::ExpansionType DetExpansionType() const
            {
                return StdRegions::eSegment;
            }

            /// Returns the locations of the quadrature points in up to
            /// three-dimensions.
            LOCAL_REGIONS_EXPORT void GetCoords(Array<OneD,NekDouble> &coords_1,
                           Array<OneD,NekDouble> &coords_2
                                        = NullNekDouble1DArray,
                           Array<OneD,NekDouble> &coords_3
                                        = NullNekDouble1DArray);

            ///
            LOCAL_REGIONS_EXPORT void GetCoord(const Array<OneD, const NekDouble>& Lcoords,
                                Array<OneD,NekDouble> &coords);


            /// Returns a pointer to the GeomFactors object describing the
            /// metric information for the segment.
            const SpatialDomains::GeomFactorsSharedPtr& GetMetricInfo() const
            {
                return m_metricinfo;
            }

            /// Returns a pointer to a Geometry object describing the
            /// geometry of the segment.
            const SpatialDomains::GeometrySharedPtr GetGeom() const
            {
                return m_geom;
            }

            /// Returns a pointer to a Geometry1D object describing the
            /// geometry of the segment.
            const SpatialDomains::Geometry1DSharedPtr& GetGeom1D() const
            {
                return m_geom;
            }

            /// Writes out the physical space data to file.
            LOCAL_REGIONS_EXPORT void WriteToFile(std::ofstream &outfile, 
                            OutputFormat format, 
                            const bool dumpVar = true, 
                            std::string var = "v");


            //----------------------------
            // Integration Methods
            //----------------------------

            /// \brief Integrate the physical point list \a inarray over region
            LOCAL_REGIONS_EXPORT NekDouble Integral(const Array<OneD, const NekDouble>& inarray);

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
            void IProductWRTBase(const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray)
            {
                IProductWRTBase(m_base[0]->GetBdata(),inarray,outarray,1);
            }

            LOCAL_REGIONS_EXPORT void IProductWRTDerivBase(const int dir,
                                      const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> & outarray);

            //-----------------------------
            // Differentiation Methods
            //-----------------------------


            /** \brief Evaluate the derivative \f$ d/d{\xi_1} \f$ at the
                physical quadrature points given by \a inarray and return in \a
                outarray. */
            LOCAL_REGIONS_EXPORT void PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                           Array<OneD,NekDouble> &out_d0,
                           Array<OneD,NekDouble> &out_d1 = NullNekDouble1DArray,
                           Array<OneD,NekDouble> &out_d2 = NullNekDouble1DArray);

            LOCAL_REGIONS_EXPORT void PhysDeriv(const int dir,
                           const Array<OneD, const NekDouble>& inarray,
                           Array<OneD, NekDouble> &outarray);

	   /**
	   *\brief Evaluate the derivative along a line: \f$ d/ds=\frac{spacedim}{||tangent||}d/d{\xi}  \f$.
	   * The derivative is calculated performing
	   *the product \f$ du/d{s}=\nabla u \cdot tangent \f$.
	   *\param inarray function to derive
	   *\param out_ds result of the derivative operation 
	   **/
           LOCAL_REGIONS_EXPORT void PhysDeriv_s(const Array<OneD, const NekDouble>& inarray, 
            	    Array<OneD, NekDouble> &out_ds);

	   /**
	   *\brief Evaluate the derivative normal to a line: \f$ d/dn=\frac{spacedim}{||normal||}d/d{\xi}  \f$.
	   * The derivative is calculated performing
	   *the product \f$ du/d{s}=\nabla u \cdot normal \f$.
	   *\param inarray function to derive
	   *\param out_dn result of the derivative operation 
	   **/     
	   LOCAL_REGIONS_EXPORT void PhysDeriv_n(const Array<OneD, const NekDouble>& inarray,
	   	    Array<OneD, NekDouble>& out_dn);
            //----------------------------
            // Evaluations Methods
            //---------------------------

            LOCAL_REGIONS_EXPORT void SetCoeffsToOrientation(StdRegions::EdgeOrientation dir);

            LOCAL_REGIONS_EXPORT void SetCoeffsToOrientation(StdRegions::EdgeOrientation dir,
                                        Array<OneD, const NekDouble> &inarray,
                                        Array<OneD, NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT void ReverseCoeffsAndSign(const Array<OneD,NekDouble> &inarray,
                                      Array<OneD,NekDouble> &outarray);

            /** \brief Inverse Matrix Product */
            LOCAL_REGIONS_EXPORT void MultiplyByElmtInvMass(
                            const Array<OneD, const NekDouble>& inarray,
                                  Array<OneD,NekDouble> &outarray);

            /** \brief Forward transform from physical quadrature space
                stored in \a inarray and evaluate the expansion coefficients and
                store in \a (this)->_coeffs  */
            LOCAL_REGIONS_EXPORT void FwdTrans(const Array<OneD, const NekDouble>& inarray,
                          Array<OneD,NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT void FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray,
                                         Array<OneD, NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT NekDouble PhysEvaluate(const Array<OneD, const NekDouble>& coord);

            LOCAL_REGIONS_EXPORT void LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray);

            LOCAL_REGIONS_EXPORT void HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                   Array<OneD,NekDouble> &outarray,
                                   const double lambda);

        protected:

            DNekMatSharedPtr GenMatrix(const StdRegions::StdMatrixKey &mkey);
            DNekMatSharedPtr CreateStdMatrix(const StdRegions::StdMatrixKey &mkey);
            DNekScalMatSharedPtr     CreateMatrix(const MatrixKey &mkey);
            DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(const MatrixKey &mkey);


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
            void IProductWRTBase(const Array<OneD, const NekDouble>& base,
                                 const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray,
                                 int coll_check);

            SpatialDomains::Geometry1DSharedPtr m_geom;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;

            virtual void v_SetUpPhysNormals(const StdRegions::StdExpansionSharedPtr &exp2D, const int edge);


	    virtual void v_SetUpPhysTangents(const StdRegions::StdExpansionSharedPtr &exp2D, const int edge);
            
            virtual StdRegions::ExpansionType v_DetExpansionType() const
            {
                return DetExpansionType();
            }

            virtual const Array<OneD, const NekDouble>& v_GetPhysNormals(void)
            {
                NEKERROR(ErrorUtil::efatal, "Got to SegExp");
                return NullNekDouble1DArray;
            }

            virtual int v_GetVertexMap(int vert)
            {
                return StdSegExp::GetVertexMap(vert);
            }


            virtual const SpatialDomains::GeomFactorsSharedPtr& v_GetMetricInfo() const
            {
                return GetMetricInfo();
            }

            virtual const SpatialDomains::GeometrySharedPtr v_GetGeom() const
            {
                return GetGeom();
            }


            virtual const SpatialDomains::Geometry1DSharedPtr& v_GetGeom1D() const
            {
                return GetGeom1D();
            }

            virtual void v_GetCoords(Array<OneD,NekDouble> &coords_0,
                                     Array<OneD,NekDouble> &coords_1 = NullNekDouble1DArray,
                                     Array<OneD,NekDouble> &coords_2 = NullNekDouble1DArray)
            {
                GetCoords(coords_0, coords_1, coords_2);
            }

            virtual void v_GetCoord(const Array<OneD, const NekDouble>& lcoord,
                                    Array<OneD,NekDouble> &coord)
            {
                GetCoord(lcoord, coord);
            }

            virtual int v_GetCoordim()
            {
                return m_geom->GetCoordim();
            }

            virtual void v_WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v")
            {
                WriteToFile(outfile,format,dumpVar,var);
            }

            virtual SpatialDomains::GeomType v_MetricInfoType()
            {
                return m_metricinfo->GetGtype();
            }

            virtual void v_SetCoeffsToOrientation(StdRegions::EdgeOrientation dir,
                                                  Array<OneD, const NekDouble> &inarray,
                                                  Array<OneD, NekDouble> &outarray)
            {
                SetCoeffsToOrientation(dir,inarray,outarray);
            }

            virtual void v_SetCoeffsToOrientation(StdRegions::EdgeOrientation dir)
            {
                SetCoeffsToOrientation(dir);
            }


            /// \brief Virtual call to integrate the physical point list \a inarray
            /// over region (see SegExp::Integral)
            virtual NekDouble v_Integral(const Array<OneD, const NekDouble>& inarray )
            {
                return Integral(inarray);
            }

            /** \brief Virtual call to SegExp::IProduct_WRT_B */
            virtual void v_IProductWRTBase(const Array<OneD, const NekDouble>& inarray,
                                           Array<OneD,NekDouble> &outarray)
            {
                IProductWRTBase(inarray,outarray);
            }

            virtual void v_IProductWRTDerivBase (const int dir,
                                                 const Array<OneD, const NekDouble> &inarray,
                                                 Array<OneD, NekDouble> &outarray)
            {
                IProductWRTDerivBase(dir,inarray,outarray);
            }

            virtual void v_NormVectorIProductWRTBase(
                    const Array<OneD, const NekDouble>   &Fx, 
                    const Array<OneD, const NekDouble> &Fy, 
                    Array< OneD, NekDouble> &outarray, 
                    bool NegateNormal = false);
                    
            virtual void v_PhysDeriv(const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD,NekDouble> &out_d0,
                                     Array<OneD,NekDouble> &out_d1 = NullNekDouble1DArray,
                                     Array<OneD,NekDouble> &out_d2 = NullNekDouble1DArray)
            {
                PhysDeriv(inarray, out_d0);
            }

            virtual void v_PhysDeriv(const int dir,
                                     const Array<OneD, const NekDouble>& inarray,
                                     Array<OneD, NekDouble> &outarray)
            {
                PhysDeriv(dir,inarray,outarray);
            }

	    virtual void v_PhysDeriv_s(const Array<OneD, const NekDouble>& inarray,
	    	    		       Array<OneD,NekDouble> &out_ds)
	    {
	        PhysDeriv_s(inarray, out_ds);
	    }

            virtual void v_PhysDeriv_n(const Array<OneD, const NekDouble>& inarray,
            	                       Array<OneD, NekDouble> &out_dn)
            {
            	PhysDeriv_n(inarray,out_dn);
            }
	    
            /// Virtual call to SegExp::FwdTrans
            virtual void v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                                    Array<OneD,NekDouble> &outarray)
            {
                FwdTrans(inarray,outarray);
            }

            virtual void v_FwdTrans_BndConstrained(const Array<OneD, const NekDouble>& inarray,
                                                   Array<OneD, NekDouble> &outarray)
            {
                FwdTrans_BndConstrained(inarray, outarray);
            }

            /// Virtual call to SegExp::Evaluate
            virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& coords)
            {
                return PhysEvaluate(coords);
            }

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

            virtual DNekMatSharedPtr v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
            {
                return CreateStdMatrix(mkey);
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

            virtual void v_LaplacianMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray)
            {
                LaplacianMatrixOp(inarray,outarray);
            }

            virtual void v_HelmholtzMatrixOp(const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,NekDouble> &outarray,
                                             const double lambda)
            {
                HelmholtzMatrixOp(inarray,outarray,lambda);
            }
            
            virtual void v_ExtractDataToCoeffs(const std::vector<NekDouble> &data, 
                                               const int offset, 
                                               const std::vector<unsigned int > &nummodes, 
                                               const int nmode_offset,
                                               Array<OneD, NekDouble> &coeffs);


        private:
            SegExp();


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

            virtual bool v_IsBoundaryInteriorExpansion()
            {
                return StdSegExp::IsBoundaryInteriorExpansion();
            }


            virtual int v_NumBndryCoeffs() const
            {
                return 2;
            }


            virtual int v_NumDGBndryCoeffs() const
            {
                return 2;
            }

            /// Virtual call to TriExp::BwdTrans
            virtual void v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                                    Array<OneD, NekDouble> &outarray)
            {
                BwdTrans(inarray,outarray);
            }

            virtual void v_BwdTrans_SumFac(const Array<OneD, const NekDouble> &inarray,
                                           Array<OneD, NekDouble> &outarray)
            {
                BwdTrans(inarray,outarray);
            }

            virtual void v_AddHDGHelmholtzTraceTerms(const NekDouble tau,
                                                     const Array<OneD, const NekDouble> &inarray,
                                                     Array<OneD,NekDouble> &outarray)
            {
                Expansion1D::AddHDGHelmholtzTraceTerms(tau,inarray,outarray);
            }


            virtual void v_AddRobinMassMatrix(const int edgeid, const Array<OneD, const NekDouble > &primCoeffs, DNekMatSharedPtr &inoutmat)
            {
                Expansion1D::AddRobinMassMatrix(edgeid,primCoeffs, inoutmat);
            }

            virtual void v_AddRobinEdgeContribution(const int edgeid, const Array<OneD, const NekDouble > &primCoeffs, Array<OneD, NekDouble> &coeffs)
            {
                Expansion1D::AddRobinEdgeContribution(edgeid,primCoeffs, coeffs);
            }

            virtual void v_GetBoundaryMap(Array<OneD, unsigned int> &maparray)
            {
                StdSegExp::GetBoundaryMap(maparray);
            }


            virtual DNekMatSharedPtr v_GenMatrix(const StdRegions::StdMatrixKey &mkey)
            {
                return GenMatrix(mkey);
            }


        };

        // type defines for use of SegExp in a boost vector
        typedef boost::shared_ptr<SegExp>      SegExpSharedPtr;
        typedef std::vector< SegExpSharedPtr > SegExpVector;
        typedef std::vector< SegExpSharedPtr >::iterator SegExpVectorIter;
    } //end of namespace
} //end of namespace

#endif // SEGEXP_H

//
// $Log: SegExp.h,v $
// Revision 1.46  2009/12/15 18:09:02  cantwell
// Split GeomFactors into 1D, 2D and 3D
// Added generation of tangential basis into GeomFactors
// Updated ADR2DManifold solver to use GeomFactors for tangents
// Added <GEOMINFO> XML session section support in MeshGraph
// Fixed const-correctness in VmathArray
// Cleaned up LocalRegions code to generate GeomFactors
// Removed GenSegExp
// Temporary fix to SubStructuredGraph
// Documentation for GlobalLinSys and GlobalMatrix classes
//
// Revision 1.45  2009/10/06 09:45:54  cbiotto
// Adding MultiplyByElmtInvMass
//
// Revision 1.44  2009/04/27 21:34:07  sherwin
// Updated WriteToField
//
// Revision 1.43  2009/04/20 16:12:28  sherwin
// Updates related to output format and optimising DG solver
//
// Revision 1.42  2009/03/04 14:17:38  pvos
// Removed all methods that take and Expansion as argument
//
// Revision 1.41  2008/08/14 22:12:57  sherwin
// Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
//
// Revision 1.40  2008/07/29 22:25:35  sherwin
// general update for DG Advection including separation of GetGeom() into GetGeom1D,2D,3D()
//
// Revision 1.39  2008/07/19 21:15:38  sherwin
// Removed MapTo function, made orientation anticlockwise, changed enum from BndSys to BndLam
//
// Revision 1.38  2008/07/12 17:27:07  sherwin
// Update for AddBoundaryInt and moved various members to be private rather than protected
//
// Revision 1.37  2008/07/04 10:19:05  pvos
// Some updates
//
// Revision 1.36  2008/07/02 14:09:18  pvos
// Implementation of HelmholtzMatOp and LapMatOp on shape level
//
// Revision 1.35  2008/05/30 00:33:48  delisi
// Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
//
// Revision 1.34  2008/05/29 21:33:37  pvos
// Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
//
// Revision 1.33  2008/05/14 18:06:50  sherwin
// mods to fix Seggeom to Geometry1D casting
//
// Revision 1.32  2008/05/10 18:27:33  sherwin
// Modifications necessary for QuadExp Unified DG Solver
//
// Revision 1.31  2008/04/06 05:59:05  bnelson
// Changed ConstArray to Array<const>
//
// Revision 1.30  2008/04/02 22:19:26  pvos
// Update for 2D local to global mapping
//
// Revision 1.29  2008/03/12 15:24:29  pvos
// Clean up of the code
//
// Revision 1.28  2008/02/28 10:04:11  sherwin
// Modes for UDG codes
//
// Revision 1.27  2008/01/21 19:59:32  sherwin
// Updated to take SegGeoms instead of EdgeComponents
//
// Revision 1.26  2007/12/17 13:04:30  sherwin
// Modified GenMatrix to take a StdMatrixKey and removed m_constant from MatrixKey
//
// Revision 1.25  2007/11/20 16:28:46  sherwin
// Added terms for UDG Helmholtz solver
//
// Revision 1.24  2007/10/03 11:37:50  sherwin
// Updates relating to static condensation implementation
//
// Revision 1.23  2007/08/10 03:38:15  jfrazier
// Updated with new rev of NekManager.
//
// Revision 1.22  2007/07/28 05:09:33  sherwin
// Fixed version with updated MemoryManager
//
// Revision 1.21  2007/07/22 23:04:19  bnelson
// Backed out Nektar::ptr.
//
// Revision 1.20  2007/07/20 00:45:51  bnelson
// Replaced boost::shared_ptr with Nektar::ptr
//
// Revision 1.19  2007/07/16 18:28:42  sherwin
// Modification to introduce non-zero Dirichlet boundary conditions into the Helmholtz1D Demo
//
// Revision 1.18  2007/07/13 09:02:23  sherwin
// Mods for Helmholtz solver
//
// Revision 1.17  2007/07/11 19:26:04  sherwin
// update for new Manager structure
//
// Revision 1.16  2007/07/10 17:17:26  sherwin
// Introduced Scaled Matrices into the MatrixManager
//
// Revision 1.15  2007/05/31 19:13:12  pvos
// Updated NodalTriExp + LocalRegions/Project2D + some other modifications
//
// Revision 1.14  2007/05/28 16:15:00  sherwin
// Updated files in MultiRegions to make 1D demos work
//
// Revision 1.13  2007/05/28 08:35:25  sherwin
// Updated for localregions up to Project1D
//
// Revision 1.12  2007/05/27 16:10:29  bnelson
// Update to new Array type.
//
// Revision 1.11  2007/04/26 15:00:16  sherwin
// SJS compiling working version using SHaredArrays
//
// Revision 1.10  2007/04/08 03:33:31  jfrazier
// Minor reformatting and fixing SharedArray usage.
//
// Revision 1.9  2007/03/25 15:48:21  sherwin
// UPdate LocalRegions to take new NekDouble and shared_array formats. Added new Demos
//
// Revision 1.8  2007/03/20 09:13:38  kirby
// new geomfactor routines; update for metricinfo; update style
//
// Revision 1.7  2007/03/14 21:24:07  sherwin
// Update for working version of MultiRegions up to ExpList1D
//
// Revision 1.6  2007/03/02 12:01:55  sherwin
// Update for working version of LocalRegions/Project1D
//
// Revision 1.5  2007/01/15 21:12:26  sherwin
// First definition
//
// Revision 1.4  2006/06/13 18:05:01  sherwin
// Modifications to make MultiRegions demo ProjectLoc2D execute properly.
//
// Revision 1.3  2006/05/30 14:00:04  sherwin
// Updates to make MultiRegions and its Demos work
//
// Revision 1.2  2006/05/29 17:05:49  sherwin
// Modified to put shared_ptr around geom definitions
//
// Revision 1.1  2006/05/04 18:58:46  kirby
// *** empty log message ***
//
// Revision 1.33  2006/03/13 18:20:33  sherwin
//
// Fixed error in ResetGmat
//
// Revision 1.32  2006/03/12 21:59:48  sherwin
//
// compiling version of LocalRegions
//
//
