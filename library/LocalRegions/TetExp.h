///////////////////////////////////////////////////////////////////////////////
//
// File TetExp.h
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

#ifndef TETEXP_H
#define TETEXP_H

#include <LocalRegions/LocalRegions.hpp>

#include <StdRegions/StdTetExp.h>
#include <SpatialDomains/TetGeom.h>

#include <SpatialDomains/GeomFactors3D.h>
#include <LocalRegions/MatrixKey.h>

#include <LocalRegions/Expansion3D.h>

namespace Nektar
{
    namespace LocalRegions
    {

        class TetExp: public StdRegions::StdTetExp, public Expansion3D
        {

        public:
            /// Constructor using BasisKey class for quadrature points and
            /// order definition.
            TetExp( const LibUtilities::BasisKey &Ba,
                            const LibUtilities::BasisKey &Bb,
                            const LibUtilities::BasisKey &Bc,
                            const SpatialDomains::TetGeomSharedPtr &geom);

            /// Copy Constructor
            TetExp(const TetExp &T);

            /// Destructor
            ~TetExp();

            void GetCoords(Array<OneD,NekDouble> &coords_0,
                Array<OneD,NekDouble> &coords_1,
                Array<OneD,NekDouble> &coords_2);

            void GetCoord(const Array<OneD, const NekDouble> &Lcoords, 
                Array<OneD,NekDouble> &coords);

        protected:
            /// Calculate the inner product of inarray with respect to the
            /// basis B=m_base0*m_base1*m_base2 and put into outarray:
            virtual void v_IProductWRTBase(
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> & outarray);

            /// Forward transform from physical quadrature space stored in
            /// \a inarray and evaluate the expansion coefficients and store
            /// in \a (this)->_coeffs  */
            virtual void v_FwdTrans(const Array<OneD, const NekDouble> & inarray,Array<OneD,NekDouble> &outarray);


            //----------------------------
            // Integration Methods
            //----------------------------
            /// Integrate the physical point list \a inarray over region
            NekDouble v_Integral(const Array<OneD, const NekDouble> &inarray);


            //-----------------------------
            // Differentiation Methods
            //-----------------------------
            /// Differentiate \a inarray in the three coordinate directions.
            virtual void v_PhysDeriv(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD, NekDouble> &out_d0,
                            Array<OneD, NekDouble> &out_d1,
                            Array<OneD, NekDouble> &out_d2);

            virtual NekDouble v_PhysEvaluate(
                            const Array<OneD, const NekDouble> &coords);

            /// Get the x,y,z coordinates of each quadrature point.
            virtual void v_GetCoords(Array<OneD,NekDouble> &coords_0,
                            Array<OneD,NekDouble> &coords_1,
                            Array<OneD,NekDouble> &coords_2)
            {
                GetCoords(coords_0, coords_1, coords_2);
            }

            virtual void v_GetCoord(
                            const Array<OneD, const NekDouble> &Lcoords,
                            Array<OneD,NekDouble> &coords)
            {
                GetCoord(Lcoords, coords);
            }

            virtual void v_WriteToFile( std::ofstream &outfile,
                            OutputFormat format,
                            const bool dumpVar = true,
                            std::string var = "v");

            virtual const SpatialDomains::GeometrySharedPtr v_GetGeom() const;

            virtual const SpatialDomains::Geometry3DSharedPtr& v_GetGeom3D()
                                                                        const;

            DNekMatSharedPtr GenMatrix(
                            const StdRegions::StdMatrixKey &mkey);

            DNekMatSharedPtr CreateStdMatrix(
                            const StdRegions::StdMatrixKey &mkey);

            DNekScalMatSharedPtr  CreateMatrix(
                            const MatrixKey &mkey);

            DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(
                            const MatrixKey &mkey);

            virtual void v_HelmholtzMatrixOp(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdRegions::StdMatrixKey &mkey);

            virtual void v_HelmholtzMatrixOp_MatFree(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdRegions::StdMatrixKey &mkey);

            virtual void v_LaplacianMatrixOp(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdRegions::StdMatrixKey &mkey);

            virtual void v_LaplacianMatrixOp(const int k1, const int k2,
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdRegions::StdMatrixKey &mkey);

            virtual void v_LaplacianMatrixOp_MatFree(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdRegions::StdMatrixKey &mkey);

        private:
            SpatialDomains::Geometry3DSharedPtr m_geom;
            SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

            LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
            LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;

            TetExp();

            void GeneralMatrixOp_MatOp(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdRegions::StdMatrixKey &mkey);

            void MultiplyByQuadratureMetric(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            /// Return Shape of region, using  ShapeType enum list.
            virtual StdRegions::ExpansionType v_DetExpansionType() const;

            virtual const SpatialDomains::GeomFactorsSharedPtr&
                                                        v_GetMetricInfo() const;

            virtual int v_GetCoordim();

            virtual NekDouble v_Linf(const Array<OneD, const NekDouble> &sol);

            virtual NekDouble v_Linf();

            virtual NekDouble v_L2(const Array<OneD, const NekDouble> &sol);

            virtual NekDouble v_L2();

            virtual DNekMatSharedPtr v_CreateStdMatrix(
                            const StdRegions::StdMatrixKey &mkey);

            virtual DNekScalMatSharedPtr& v_GetLocMatrix(
                            const MatrixKey &mkey);

            virtual DNekScalBlkMatSharedPtr& v_GetLocStaticCondMatrix(
                            const MatrixKey &mkey);
        };

        // type defines for use of TetExp in a boost vector
        typedef boost::shared_ptr<TetExp> TetExpSharedPtr;
        typedef std::vector< TetExpSharedPtr > TetExpVector;
        typedef std::vector< TetExpSharedPtr >::iterator TetExpVectorIter;

    } //end of namespace
} //end of namespace

#endif // TETEXP_H

/**
 *    $Log: TetExp.h,v $
 *    Revision 1.24  2010/02/26 13:52:45  cantwell
 *    Tested and fixed where necessary Hex/Tet projection and differentiation in
 *      StdRegions, and LocalRegions for regular and deformed (where applicable).
 *    Added SpatialData and SpatialParameters classes for managing spatiall-varying
 *      data.
 *    Added TimingGeneralMatrixOp3D for timing operations on 3D geometries along
 *      with some associated input meshes.
 *    Added 3D std and loc projection demos for tet and hex.
 *    Added 3D std and loc regression tests for tet and hex.
 *    Fixed bugs in regression tests in relation to reading OK files.
 *    Extended Elemental and Global optimisation parameters for 3D expansions.
 *    Added GNUPlot output format option.
 *    Updated ADR2DManifoldSolver to use spatially varying data.
 *    Added Barkley model to ADR2DManifoldSolver.
 *    Added 3D support to FldToVtk and XmlToVtk.
 *    Renamed History.{h,cpp} to HistoryPoints.{h,cpp}
 *
 *    Revision 1.23  2009/12/15 18:09:02  cantwell
 *    Split GeomFactors into 1D, 2D and 3D
 *    Added generation of tangential basis into GeomFactors
 *    Updated ADR2DManifold solver to use GeomFactors for tangents
 *    Added <GEOMINFO> XML session section support in MeshGraph
 *    Fixed const-correctness in VmathArray
 *    Cleaned up LocalRegions code to generate GeomFactors
 *    Removed GenSegExp
 *    Temporary fix to SubStructuredGraph
 *    Documentation for GlobalLinSys and GlobalMatrix classes
 *
 *    Revision 1.22  2009/04/27 21:34:07  sherwin
 *    Updated WriteToField
 *
 *    Revision 1.21  2009/04/20 16:12:28  sherwin
 *    Updates related to output format and optimising DG solver
 *
 *    Revision 1.20  2008/08/14 22:12:57  sherwin
 *    Introduced Expansion classes and used them to define HDG routines, has required quite a number of virtual functions to be added
 *
 *    Revision 1.19  2008/07/29 22:25:35  sherwin
 *    general update for DG Advection including separation of GetGeom() into GetGeom1D,2D,3D()
 *
 *    Revision 1.18  2008/07/12 17:27:07  sherwin
 *    Update for AddBoundaryInt and moved various members to be private rather than protected
 *
 *    Revision 1.17  2008/07/04 10:19:05  pvos
 *    Some updates
 *
 *    Revision 1.16  2008/06/06 23:25:43  ehan
 *    Added doxygen documentation
 *
 *    Revision 1.15  2008/05/30 00:33:48  delisi
 *    Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 *    Revision 1.14  2008/05/29 21:33:37  pvos
 *    Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 *    Revision 1.13  2008/05/10 18:27:33  sherwin
 *    Modifications necessary for QuadExp Unified DG Solver
 *
 *    Revision 1.12  2008/04/06 05:59:05  bnelson
 *    Changed ConstArray to Array<const>
 *
 *    Revision 1.11  2008/03/18 16:18:13  pvos
 *    Fixed some compiler complaints
 *
 *    Revision 1.10  2008/03/17 10:36:28  pvos
 *    Clean up of the code
 *
 *    Revision 1.9  2008/03/12 15:24:29  pvos
 *    Clean up of the code
 *
 *    Revision 1.8  2008/02/16 05:53:01  ehan
 *    Added PhysDeriv and virtual functions.
 *
 *    Revision 1.7  2008/02/05 00:41:37  ehan
 *    Added initial tetrahedral expansion.
 *
 *    Revision 1.6  2007/07/22 23:04:19  bnelson
 *    Backed out Nektar::ptr.
 *
 *    Revision 1.5  2007/07/20 00:45:51  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.4  2007/01/15 21:12:26  sherwin
 *    First definition
 *
 *    Revision 1.3  2006/06/13 18:05:01  sherwin
 *    Modifications to make MultiRegions demo ProjectLoc2D execute properly.
 *
 *    Revision 1.2  2006/05/30 14:00:04  sherwin
 *    Updates to make MultiRegions and its Demos work
 *
 *    Revision 1.1  2006/05/04 18:58:46  kirby
 *    *** empty log message ***
 *
 *    Revision 1.9  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.8  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/
