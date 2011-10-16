///////////////////////////////////////////////////////////////////////////////
//
// File StdTetExp.h
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
// Description: Header field for tetrahedral routines built upon
// StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_STDREGIONS_STDTETEXP_H
#define NEKTAR_LIB_STDREGIONS_STDTETEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion3D.h>
#include <StdRegions/StdMatrixKey.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
    namespace StdRegions
    {

        namespace StdTetData
        {
            // Adds up the number of cells in a truncated Nc by Nc by Nc pyramid,
            // where the longest Na rows and longest Nb columns are kept.
            // Example: (Na, Nb, Nc) = (3, 4, 5); The number of coefficients is the
            // sum of the elements of the following matrix:
            //     |5  4  3  2  0|
            //     |4  3  2  0   |
            //     |3  2  0      |
            //     |0  0         |
            //     |0            |
            // Sum = 28 = number of tet coefficients
            inline int getNumberOfCoefficients( int Na, int Nb, int Nc )
            {
                int nCoef = 0;
                for( int a = 0; a < Na; ++a )
                {
                    for( int b = 0; b < Nb - a; ++b )
                    {
                        for( int c = 0; c < Nc - a - b; ++c )
                        {
                            ++nCoef;
                        }
                    }
                }
                return nCoef;
            }
        }

        class StdTetExp: virtual public StdExpansion3D
        {

        public:

            STD_REGIONS_EXPORT StdTetExp();
            /** \brief Constructor using BasisKey class for quadrature points and order definition
             */
            STD_REGIONS_EXPORT StdTetExp(const  LibUtilities::BasisKey &Ba, const  LibUtilities::BasisKey &Bb, const  LibUtilities::BasisKey &Bc);

            /** \brief Constructor using BasisKey class for quadrature points and order definition where m_coeffs and m_phys are all set.
             */
            STD_REGIONS_EXPORT StdTetExp(const  LibUtilities::BasisKey &Ba, const  LibUtilities::BasisKey &Bb, const  LibUtilities::BasisKey &Bc, double *coeffs, double *phys);

            /** \brief Copy Constructor */
            STD_REGIONS_EXPORT StdTetExp(const StdTetExp &T);

            /** \brief Destructor */
            STD_REGIONS_EXPORT ~StdTetExp();

            /** \brief Return Shape of region, using  ShapeType enum list.
             *  i.e. Tetrahedron
             */
            ExpansionType DetExpansionType() const
            {
                return eTetrahedron;
            }

            STD_REGIONS_EXPORT void TripleTensorProduct(
                                const Array<OneD, const NekDouble>& fx,
                                const Array<OneD, const NekDouble>& gy,
                                const Array<OneD, const NekDouble>& hz,
                                const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> & outarray );

            STD_REGIONS_EXPORT NekDouble TripleInnerProduct(
                                const Array<OneD, const NekDouble>& fxyz,
                                const Array<OneD, const NekDouble>& wx,
                                const Array<OneD, const NekDouble>& wy,
                                const Array<OneD, const NekDouble>& wz );

            STD_REGIONS_EXPORT NekDouble Integral3D(const Array<OneD, const NekDouble>& inarray,
                                const Array<OneD, const NekDouble>& w0,
                                const Array<OneD, const NekDouble>& w1,
                                const Array<OneD, const NekDouble>&w2);

            /// Integrate the physical point list \a inarray over tetrahedral
            /// region and return the value
            STD_REGIONS_EXPORT NekDouble Integral(const Array<OneD, const NekDouble>& inarray);


            STD_REGIONS_EXPORT void FillMode(const int mode, Array<OneD, NekDouble> &outarray);


            /// Backward tranform for tetrahedral elements
//            void BwdTrans(const Array<OneD, const NekDouble>& inarray,
//                          Array<OneD, NekDouble> &outarray);

            /// Forward transform from physical quadrature space stored in \a
            /// inarray and evaluate the expansion coefficients and store in \a
            /// outarray.
//            void FwdTrans(const Array<OneD, const NekDouble>& inarray,
//                          Array<OneD, NekDouble> &outarray);

            /** \brief Single Point Evaluation */
            STD_REGIONS_EXPORT NekDouble PhysEvaluate3D(const Array<OneD, const NekDouble>& coords);

            STD_REGIONS_EXPORT NekDouble PhysEvaluate3D(const Array<OneD, const NekDouble>& coords,  const Array<OneD, const NekDouble> & physvals);


            //----------------------------------
            // Generate Matrix Routine
            //----------------------------------

            DNekMatSharedPtr GenMatrix(const StdMatrixKey &mkey)
            {
                return StdExpansion::CreateGeneralMatrix(mkey);
            }

        protected:

            /// Compute backward transform of modes to quadrature points.
            STD_REGIONS_EXPORT virtual void v_BwdTrans(const Array<OneD, const NekDouble>& inarray,
                                    Array<OneD, NekDouble> &outarray );

            /// Performs the sum factorisation form of the BwdTrans operation.
            STD_REGIONS_EXPORT void BwdTrans_SumFacKernel(
                    const Array<OneD, const NekDouble>& base0,
                    const Array<OneD, const NekDouble>& base1,
                    const Array<OneD, const NekDouble>& base2,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray,
                          Array<OneD, NekDouble> &wsp,
                    bool doCheckCollDir0,
                    bool doCheckCollDir1,
                    bool doCheckCollDir2);

            /// Perform forward transform of quadrature data to coefficients.
            STD_REGIONS_EXPORT virtual void v_FwdTrans(const Array<OneD, const NekDouble>& inarray,
                                    Array<OneD, NekDouble> &outarray );

            /// Inner product of \a inarray over region with respect to the
            /// expansion basis m_base[0]->GetBdata(),m_base[1]->GetBdata(),
            /// m_base[2]->GetBdata() and return in \a outarray.
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> & outarray);

            /// Fundamental Tet sum-factorisation implementation of IProduct.
            STD_REGIONS_EXPORT void IProductWRTBase_SumFacKernel (
                    const Array<OneD, const NekDouble>& base0,
                    const Array<OneD, const NekDouble>& base1,
                    const Array<OneD, const NekDouble>& base2,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray,
                          Array<OneD, NekDouble> &wsp,
                    bool doCheckCollDir0,
                    bool doCheckCollDir1,
                    bool doCheckCollDir2);

            /// Compute inner product with respect to derivative basis
            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase(const int dir,
                    const Array<OneD, const NekDouble>& inarray,
                    Array<OneD, NekDouble> & outarray);

            /// Compute inner product with respect to derivative basis using
            /// sum-factorisation technique.
            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase_SumFac(const int dir,
                                             const Array<OneD, const NekDouble>& inarray,
                                             Array<OneD, NekDouble> & outarray);
            /// Compute inner product with respect to derivative basis using
            /// local matrix operation.
            STD_REGIONS_EXPORT virtual void v_IProductWRTDerivBase_MatOp(const int dir,
                                            const Array<OneD, const NekDouble>& inarray,
                                            Array<OneD, NekDouble> & outarray);


            STD_REGIONS_EXPORT virtual void v_PhysDeriv(
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &out_dx,
                            Array<OneD, NekDouble> &out_dy,
                            Array<OneD, NekDouble> &out_dz );

            /// Calculate the derivative of the physical points in a single
            /// direction.
            STD_REGIONS_EXPORT virtual void v_PhysDeriv(const int dir,
                           const Array<OneD, const NekDouble>& inarray,
                           Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_PhysDirectionalDeriv(
                            const Array<OneD, const NekDouble>& inarray,
                            const Array<OneD, const NekDouble>& direction,
                            Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_StdPhysDeriv(
                            const Array<OneD, const NekDouble>& inarray,
                            Array<OneD, NekDouble> &out_d0,
                            Array<OneD, NekDouble> &out_d1,
                            Array<OneD, NekDouble> &out_d2);

            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& coords);

            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(const Array<OneD, const NekDouble>& coords,  const Array<OneD, const NekDouble> & physvals);

        private:
            STD_REGIONS_EXPORT virtual void v_BwdTrans_SumFac(const Array<OneD, const NekDouble>& inarray,
                                 Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_IProductWRTBase_MatOp (
                                const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> & outarray);

            STD_REGIONS_EXPORT virtual void v_IProductWRTBase_SumFac (
                                const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> & outarray);

            STD_REGIONS_EXPORT void MultiplyByQuadratureMetric(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            /// The type of expansion.
            STD_REGIONS_EXPORT virtual ExpansionType v_DetExpansionType() const;

            STD_REGIONS_EXPORT virtual bool v_IsBoundaryInteriorExpansion();

            /// Number of vertices in a tetrahedron.
            STD_REGIONS_EXPORT virtual int  v_GetNverts() const;

            /// Number of edges in a tetrahedron.
            STD_REGIONS_EXPORT virtual int  v_GetNedges() const;

            /// Number of faces in a tetrahedron.
            STD_REGIONS_EXPORT virtual int  v_GetNfaces() const;

            /// Number of boundary coefficients in this tetrahedron.
            STD_REGIONS_EXPORT virtual int  v_NumBndryCoeffs() const;

            /// Number of coefficients on a given face.
            STD_REGIONS_EXPORT virtual int  v_GetFaceNcoeffs(const int i) const;

            /// Number of coefficients interior to a given face.
            STD_REGIONS_EXPORT virtual int  v_GetFaceIntNcoeffs(const int i) const;

            /// Number of coefficients in a given edge.
            STD_REGIONS_EXPORT virtual int  v_GetEdgeNcoeffs(const int i) const;

            STD_REGIONS_EXPORT virtual int  v_CalcNumberOfCoefficients(const std::vector<unsigned int> &nummodes, int &modes_offset);

            STD_REGIONS_EXPORT virtual void v_GetBoundaryMap(Array<OneD, unsigned int>& outarray);

            STD_REGIONS_EXPORT virtual void v_GetInteriorMap(Array<OneD, unsigned int>& outarray);

            STD_REGIONS_EXPORT virtual int  v_GetVertexMap(int localVertexId);

            STD_REGIONS_EXPORT virtual void v_GetEdgeInteriorMap(
                            const int eid, const EdgeOrientation edgeOrient,
                            Array<OneD, unsigned int> &maparray,
                            Array<OneD, int> &signarray);

            STD_REGIONS_EXPORT virtual void v_GetFaceInteriorMap(const int fid, const FaceOrientation faceOrient,
                                              Array<OneD, unsigned int> &maparray,
                                              Array<OneD, int> &signarray);

            STD_REGIONS_EXPORT virtual void v_GetFaceToElementMap(const int fid,
                            const FaceOrientation faceOrient,
                            Array<OneD, unsigned int> &maparray,
                            Array<OneD, int>& signarray);

            STD_REGIONS_EXPORT int GetMode(const int i, const int j, const int k);

            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_GenMatrix(const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_CreateStdMatrix(const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual LibUtilities::BasisType v_GetEdgeBasisType(const int i) const;

            STD_REGIONS_EXPORT virtual void v_GetCoords( Array<OneD, NekDouble> &coords_x,
                            Array<OneD, NekDouble> &coords_y,
                            Array<OneD, NekDouble> &coords_z);

            STD_REGIONS_EXPORT virtual NekDouble v_Integral(
                            const Array<OneD, const NekDouble>& inarray );

            STD_REGIONS_EXPORT virtual void v_FillMode(
                            const int mode, Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT virtual void v_WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar = true, std::string var = "v");

            STD_REGIONS_EXPORT virtual void v_WriteCoeffsToFile(std::ofstream &outfile);

        };

        typedef boost::shared_ptr<StdTetExp> StdTetExpSharedPtr;

    } //end of namespace
} //end of namespace

#endif //STDTETEXP_H

/**
 * $Log: StdTetExp.h,v $
 * Revision 1.30  2010/02/26 13:52:46  cantwell
 * Tested and fixed where necessary Hex/Tet projection and differentiation in
 *   StdRegions, and LocalRegions for regular and deformed (where applicable).
 * Added SpatialData and SpatialParameters classes for managing spatiall-varying
 *   data.
 * Added TimingGeneralMatrixOp3D for timing operations on 3D geometries along
 *   with some associated input meshes.
 * Added 3D std and loc projection demos for tet and hex.
 * Added 3D std and loc regression tests for tet and hex.
 * Fixed bugs in regression tests in relation to reading OK files.
 * Extended Elemental and Global optimisation parameters for 3D expansions.
 * Added GNUPlot output format option.
 * Updated ADR2DManifoldSolver to use spatially varying data.
 * Added Barkley model to ADR2DManifoldSolver.
 * Added 3D support to FldToVtk and XmlToVtk.
 * Renamed History.{h,cpp} to HistoryPoints.{h,cpp}
 *
 * Revision 1.29  2009/12/15 18:09:02  cantwell
 * Split GeomFactors into 1D, 2D and 3D
 * Added generation of tangential basis into GeomFactors
 * Updated ADR2DManifold solver to use GeomFactors for tangents
 * Added <GEOMINFO> XML session section support in MeshGraph
 * Fixed const-correctness in VmathArray
 * Cleaned up LocalRegions code to generate GeomFactors
 * Removed GenSegExp
 * Temporary fix to SubStructuredGraph
 * Documentation for GlobalLinSys and GlobalMatrix classes
 *
 * Revision 1.28  2009/11/10 19:02:20  sehunchun
 * *** empty log message ***
 *
 * Revision 1.27  2009/04/27 21:32:45  sherwin
 * Updated WriteToField method
 *
 * Revision 1.26  2009/04/20 16:11:47  sherwin
 * Mods to handle output and optimise DG work
 *
 * Revision 1.25  2009/01/01 02:34:46  ehan
 * Added virtual functions.
 *
 * Revision 1.24  2008/11/17 09:01:58  ehan
 * Implemented GetFaceNcoeffs
 *
 * Revision 1.23  2008/09/17 13:46:06  pvos
 * Added LocalToGlobalC0ContMap for 3D expansions
 *
 * Revision 1.22  2008/07/19 21:12:54  sherwin
 * Removed MapTo function and made orientation convention anticlockwise in UDG routines
 *
 * Revision 1.21  2008/07/04 10:18:41  pvos
 * Some updates
 *
 * Revision 1.20  2008/06/16 22:46:51  ehan
 * Populated the function GetFaceToElementMap(..)
 *
 * Revision 1.19  2008/06/06 23:23:05  ehan
 * Added doxygen documentation
 *
 * Revision 1.18  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.17  2008/05/29 21:36:25  pvos
 * Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 * Revision 1.16  2008/05/15 22:42:29  ehan
 * Added WriteToFile() function and its virtual function
 *
 * Revision 1.15  2008/05/15 04:15:52  ehan
 * Added virtual function v_CreatStdMatrix()
 *
 * Revision 1.14  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.13  2008/03/25 08:40:20  ehan
 * Added GetEdgeNcoeffs() and GetEdgeBasisType().
 *
 * Revision 1.12  2008/03/17 10:37:58  pvos
 * Clean up of the code
 *
 * Revision 1.11  2008/01/08 22:48:41  ehan
 * Fixed the call signature of a shadowed virtual function: Added a const qualifier to the passed parameter StdMatrixKey in the virtual function v_GenMatrix().  This enables Nektar to generate the correct standard mass matrix at initialization time.
 *
 * Revision 1.10  2007/12/17 13:03:51  sherwin
 * Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
 *
 * Revision 1.9  2007/10/28 18:32:34  bnelson
 * Fixed visual c++ compile error.
 *
 * Revision 1.8  2007/10/15 20:40:07  ehan
 * Completed Basis, Backward, and Forward transformation
 *
 * Revision 1.7  2007/07/20 02:16:55  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.6  2007/07/10 21:05:17  kirby
 * even more fixes
 *
 * Revision 1.5  2007/07/10 20:41:52  kirby
 * more fixes
 *
 * Revision 1.4  2007/01/17 16:05:41  pvos
 * updated doxygen documentation
 *
 * Revision 1.3  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
 *
 * Revision 1.2  2006/07/02 17:16:19  sherwin
 *
 * Modifications to make MultiRegions work for a connected domain in 2D (Tris)
 *
 * Revision 1.1  2006/05/04 18:58:33  kirby
 * *** empty log message ***
 *
 * Revision 1.24  2006/03/12 21:59:48  sherwin
 *
 * compiling version of LocalRegions
 *
 * Revision 1.23  2006/03/06 17:12:46  sherwin
 *
 * Updated to properly execute all current StdRegions Demos.
 *
 * Revision 1.22  2006/03/04 20:26:55  bnelson
 * Added comments after #endif.
 *
 * Revision 1.21  2006/03/01 08:25:05  sherwin
 *
 * First compiling version of StdRegions
 *
 **/

