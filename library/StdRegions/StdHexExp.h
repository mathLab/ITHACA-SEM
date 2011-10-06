///////////////////////////////////////////////////////////////////////////////
//
// File StdHexExp.h
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
// Description: Hex routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_STDREGSION_STDHEXEXP_H
#define NEKTAR_LIBS_STDREGSION_STDHEXEXP_H

#include <StdRegions/StdRegions.hpp>
#include <StdRegions/StdExpansion3D.h>
#include <StdRegions/StdRegionsDeclspec.h>

namespace Nektar
{
    namespace StdRegions
    {
        /// Class representing a hexehedral element in reference space.
        class StdHexExp: virtual public StdExpansion3D
        {

        public:
            /// Default constructor.
            STD_REGIONS_EXPORT StdHexExp();

            /// Constructor using BasisKey class for quadrature points and
            /// order definition.
            STD_REGIONS_EXPORT StdHexExp(  const  LibUtilities::BasisKey &Ba,
                        const  LibUtilities::BasisKey &Bb,
                        const  LibUtilities::BasisKey &Bc);

            /// Constructor using BasisKey class for quadrature points and
            /// order definition where m_coeffs and m_phys are all set.
            STD_REGIONS_EXPORT StdHexExp(  const  LibUtilities::BasisKey &Ba,
                        const  LibUtilities::BasisKey &Bb,
                        const  LibUtilities::BasisKey &Bc,
                        double *coeffs,
                        double *phys);

            /// Copy constructor.
            STD_REGIONS_EXPORT StdHexExp(const StdHexExp &T);

            /// Destructor.
            STD_REGIONS_EXPORT ~StdHexExp();


            /// Fill outarray with mode \a mode of expansion
            STD_REGIONS_EXPORT void FillMode(const int mode, Array<OneD, NekDouble> &outarray);


            //-------------------------------
            // Integration Methods
            //-------------------------------
            STD_REGIONS_EXPORT void TripleTensorProduct(const Array<OneD, const NekDouble>& fx,
                                const Array<OneD, const NekDouble>& gy,
                                const Array<OneD, const NekDouble>& hz,
                                const Array<OneD, const NekDouble>& inarray,
                                      Array<OneD, NekDouble> & outarray );

            STD_REGIONS_EXPORT NekDouble TripleInnerProduct(
                                const Array<OneD, const NekDouble>& fxyz,
                                const Array<OneD, const NekDouble>& wx,
                                const Array<OneD, const NekDouble>& wy,
                                const Array<OneD, const NekDouble>& wz);

            STD_REGIONS_EXPORT NekDouble Integral3D(const Array<OneD, const NekDouble>& inarray,
                                 const Array<OneD, const NekDouble>& wx,
                                 const Array<OneD, const NekDouble>& wy,
                                 const Array<OneD, const NekDouble>& wz);

            /// Integrate the physical point list \a inarray over hexahedral
            /// region and return the value
            STD_REGIONS_EXPORT NekDouble Integral(const Array<OneD, const NekDouble>& inarray);


            //----------------------------
            // Differentiation Methods
            //----------------------------

            /// Calculate the deritive of the physical points
            STD_REGIONS_EXPORT void PhysDeriv(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &out_d0,
                          Array<OneD, NekDouble> &out_d1,
                          Array<OneD, NekDouble> &out_d2);

            STD_REGIONS_EXPORT void StdPhysDeriv(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &out_d0,
                          Array<OneD, NekDouble> &out_d1,
                          Array<OneD, NekDouble> &out_d2);

            STD_REGIONS_EXPORT void WriteToFile(std::ofstream &outfile,
                        OutputFormat format,
                        const bool dumpVar = true,
                        std::string var = "v");

            STD_REGIONS_EXPORT void WriteCoeffsToFile(std::ofstream &outfile);

            //----------------------------------
            // Local Matrix Routines
            //----------------------------------
            STD_REGIONS_EXPORT DNekMatSharedPtr GenMassMatrix();

            STD_REGIONS_EXPORT DNekMatSharedPtr GenLaplacianMatrix();

            STD_REGIONS_EXPORT DNekMatSharedPtr GenLaplacianMatrix(const int i, const int j);

            STD_REGIONS_EXPORT DNekMatSharedPtr GenWeakDerivMatrix(const int i);

            STD_REGIONS_EXPORT DNekMatSharedPtr GenNBasisTransMatrix();

            STD_REGIONS_EXPORT DNekMatSharedPtr GenBwdTransMatrix();


            DNekMatSharedPtr GenMatrix(const StdMatrixKey &mkey)
            {
                return StdExpansion::CreateGeneralMatrix(mkey);
            }

        protected:
            /// Backward transformation is evaluated at the quadrature points
            /// \f$ u^{\delta} (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{m(pqr)}
            /// \hat u_{pqr} \phi_{pqr} (\xi_{1i}, \xi_{2j}, \xi_{3k})\f$
            STD_REGIONS_EXPORT virtual void v_BwdTrans(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

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

            /// Forward transform from physical quadrature space stored in \a
            /// inarray and evaluate the expansion coefficients and store in
            /// #m_coeffs
            STD_REGIONS_EXPORT virtual void v_FwdTrans(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            /// Inner product of \a inarray over region with respect to the
            /// expansion basis m_base[0]->GetBdata(),m_base[1]->GetBdata(),
            /// m_base[2]->GetBdata() and return in \a outarray.
            /**
             * Wrapper call to StdHexExp::IProductWRTBase
             * Input:\n
             * - \a inarray: array of function evaluated at the physical
             *   collocation points
             * Output:\n
             * - \a outarray: array of inner product with respect to each basis
             *   over region
             */
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            /// Calculate the inner product of inarray with respect to the
            /// basis B = base0 * base1 * base2.
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase(
                    const Array<OneD, const NekDouble>& bx,
                    const Array<OneD, const NekDouble>& by,
                    const Array<OneD, const NekDouble>& bz,
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> & outarray,
                    int coll_check);

            // the arguments doCheckCollDir0 and doCheckCollDir1 allow you to specify whether
            // to check if the basis has collocation properties (i.e. for the classical spectral 
            // element basis, In this case the 1D 'B' matrix is equal to the identity matrix
            // which can be exploited to speed up the calculations).
            // However, as this routine also allows to pass the matrix 'DB' (derivative of the basis),
            // the collocation property cannot always be used. Therefor follow this rule:
            // if base0 == m_base[0]->GetBdata() --> set doCheckCollDir0 == true;
            //    base1 == m_base[1]->GetBdata() --> set doCheckCollDir1 == true;
            //    base0 == m_base[0]->GetDbdata() --> set doCheckCollDir0 == false;
            //    base1 == m_base[1]->GetDbdata() --> set doCheckCollDir1 == false;
            STD_REGIONS_EXPORT void IProductWRTBase_SumFacKernel(
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
            STD_REGIONS_EXPORT void IProductWRTDerivBase_SumFac(const int dir, 
                                             const Array<OneD, const NekDouble>& inarray, 
                                             Array<OneD, NekDouble> & outarray);
            /// Compute inner product with respect to derivative basis using
            /// local matrix operation.            
            STD_REGIONS_EXPORT void IProductWRTDerivBase_MatOp(const int dir, 
                                            const Array<OneD, const NekDouble>& inarray, 
                                            Array<OneD, NekDouble> & outarray);

            STD_REGIONS_EXPORT virtual NekDouble v_PhysEvaluate(Array<OneD,
                                const NekDouble>& Lcoords);

            STD_REGIONS_EXPORT virtual void v_MassMatrixOp(
                            const Array<OneD, const NekDouble> &inarray, 
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_LaplacianMatrixOp(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_LaplacianMatrixOp(const int k1, const int k2, 
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey);

            STD_REGIONS_EXPORT virtual void v_WeakDerivMatrixOp(const int i,
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey);
            
            STD_REGIONS_EXPORT virtual void v_HelmholtzMatrixOp(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey);
            
            STD_REGIONS_EXPORT virtual void v_LaplacianMatrixOp_MatFree(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey);
            
            STD_REGIONS_EXPORT virtual void v_HelmholtzMatrixOp_MatFree(
                            const Array<OneD, const NekDouble> &inarray,
                            Array<OneD,NekDouble> &outarray,
                            const StdMatrixKey &mkey);


        private:
            STD_REGIONS_EXPORT virtual void v_BwdTrans_SumFac(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);
                          
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase_MatOp(
                    const Array<OneD, const NekDouble>& inarray, 
                          Array<OneD, NekDouble> &outarray);
                          
            STD_REGIONS_EXPORT virtual void v_IProductWRTBase_SumFac(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);
                          
            STD_REGIONS_EXPORT void MultiplyByQuadratureMetric(
                    const Array<OneD, const NekDouble>& inarray,
                          Array<OneD, NekDouble> &outarray);

            STD_REGIONS_EXPORT void GeneralMatrixOp_MatOp(const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,NekDouble> &outarray,
                                       const StdMatrixKey &mkey);    

                                                    
            //---------------------------
            // Helper functions
            //---------------------------
            STD_REGIONS_EXPORT virtual int  v_GetNverts() const;
            STD_REGIONS_EXPORT virtual int  v_GetNedges() const;
            STD_REGIONS_EXPORT virtual int  v_GetNfaces() const;
            STD_REGIONS_EXPORT virtual ExpansionType v_DetExpansionType() const;
            STD_REGIONS_EXPORT virtual int  v_NumBndryCoeffs() const;
            STD_REGIONS_EXPORT virtual int  v_GetFaceNcoeffs(const int i) const;
            STD_REGIONS_EXPORT virtual int  v_GetFaceIntNcoeffs(const int i) const;
            STD_REGIONS_EXPORT virtual int  v_CalcNumberOfCoefficients(const std::vector<unsigned int> &nummodes, int &modes_offset);

            //--------------------------
            // Mappings
            //--------------------------
            STD_REGIONS_EXPORT virtual void v_GetBoundaryMap(Array<OneD, unsigned int>& outarray);
            STD_REGIONS_EXPORT virtual void v_GetInteriorMap(Array<OneD, unsigned int>& outarray);
            STD_REGIONS_EXPORT virtual int  v_GetVertexMap(int localVertexId);
            STD_REGIONS_EXPORT virtual void v_GetEdgeInteriorMap(const int eid,
                                const EdgeOrientation edgeOrient,
                                Array<OneD, unsigned int> &maparray,
                                Array<OneD, int> &signarray);
            STD_REGIONS_EXPORT virtual void v_GetFaceInteriorMap(const int fid,
                                const FaceOrientation faceOrient,
                                Array<OneD, unsigned int> &maparray,
                                Array<OneD, int>& signarray);
            STD_REGIONS_EXPORT virtual void v_GetFaceToElementMap(const int fid,
                                const FaceOrientation faceOrient,
                                Array<OneD, unsigned int> &maparray,
                                Array<OneD, int>& signarray);
                                
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_GenMatrix(const StdMatrixKey &mkey);
            STD_REGIONS_EXPORT virtual DNekMatSharedPtr v_CreateStdMatrix(
                                const StdMatrixKey &mkey);
            STD_REGIONS_EXPORT virtual LibUtilities::BasisType v_GetEdgeBasisType(
                                const int i) const;
            STD_REGIONS_EXPORT virtual void v_GetCoords(Array<OneD, NekDouble> &coords_x,
                                     Array<OneD, NekDouble> &coords_y,
                                     Array<OneD, NekDouble> &coords_z);
            STD_REGIONS_EXPORT virtual NekDouble v_Integral(const Array<OneD,
                                const NekDouble>& inarray );
                                
            STD_REGIONS_EXPORT virtual void v_FillMode(const int mode,
                                Array<OneD, NekDouble> &outarray);
            STD_REGIONS_EXPORT virtual void v_PhysDeriv(const Array<OneD,
                                const NekDouble>& inarray,
                                Array<OneD, NekDouble> &out_d0,
                                Array<OneD, NekDouble> &out_d1,
                                Array<OneD, NekDouble> &out_d2);
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

            STD_REGIONS_EXPORT virtual int  v_GetEdgeNcoeffs(const int i) const;
            STD_REGIONS_EXPORT virtual void v_WriteToFile(std::ofstream &outfile,
                                OutputFormat format,
                                const bool dumpVar = true,
                                std::string var = "v");
            STD_REGIONS_EXPORT virtual void v_WriteCoeffsToFile(std::ofstream &outfile);
        };

        typedef boost::shared_ptr<StdHexExp> StdHexExpSharedPtr;


        inline void StdHexExp::StdPhysDeriv(
                const Array<OneD, const NekDouble>& inarray,
                      Array<OneD, NekDouble> &out_d0,
                      Array<OneD, NekDouble> &out_d1,
                      Array<OneD, NekDouble> &out_d2)
        {
            PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }
    } //end of namespace
} //end of namespace

#endif //STDHEXEXP_H

/**
 * $Log: StdHexExp.h,v $
 * Revision 1.32  2010/03/07 14:45:22  cantwell
 * Added support for solving Helmholtz on Hexes
 *
 * Revision 1.31  2010/02/26 13:52:46  cantwell
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
 * Revision 1.30  2009/12/15 18:09:02  cantwell
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
 * Revision 1.29  2009/11/10 19:02:20  sehunchun
 * *** empty log message ***
 *
 * Revision 1.28  2009/04/27 21:32:45  sherwin
 * Updated WriteToField method
 *
 * Revision 1.27  2008/09/17 13:46:06  pvos
 * Added LocalToGlobalC0ContMap for 3D expansions
 *
 * Revision 1.26  2008/09/15 13:18:08  pvos
 * Added more hexahedron mapping routines
 *
 * Revision 1.25  2008/09/12 11:26:39  pvos
 * Updates for mappings in 3D
 *
 * Revision 1.24  2008/07/04 10:18:40  pvos
 * Some updates
 *
 * Revision 1.23  2008/06/16 22:45:51  ehan
 * Populated the function GetFaceToElementMap(..)
 *
 * Revision 1.22  2008/06/05 15:06:06  pvos
 * Added documentation
 *
 * Revision 1.21  2008/05/30 00:33:49  delisi
 * Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
 *
 * Revision 1.20  2008/05/29 21:36:25  pvos
 * Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
 *
 * Revision 1.19  2008/05/15 22:40:05  ehan
 * Clean up the codes
 *
 * Revision 1.18  2008/05/15 04:14:37  ehan
 * Added virtual function v_CreatStdMatrix()
 *
 * Revision 1.17  2008/04/06 06:04:15  bnelson
 * Changed ConstArray to Array<const>
 *
 * Revision 1.16  2008/03/25 08:39:45  ehan
 * Added GetEdgeNcoeffs() and GetEdgeBasisType().
 *
 * Revision 1.15  2008/03/17 10:37:32  pvos
 * Clean up of the code
 *
 * Revision 1.14  2008/01/20 06:09:37  bnelson
 * Fixed visual c++ compile errors.
 *
 * Revision 1.13  2008/01/08 22:30:43  ehan
 * Clean up the codes.
 *
 * Revision 1.12  2007/12/17 13:03:51  sherwin
 * Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
 *
 * Revision 1.11  2007/12/01 00:52:32  ehan
 * Completed implementing and testing following functions:
 * Integral, IProductWRTBase, PhysDeriv. BwdTrans, FwdTrans, and PhysEvaluate.
 *
 * Revision 1.10  2007/07/20 02:16:54  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.9  2007/07/10 21:05:16  kirby
 * even more fixes
 *
 * Revision 1.7  2007/01/17 16:36:58  pvos
 * updating doxygen documentation
 *
 * Revision 1.6  2007/01/17 16:05:40  pvos
 * updated doxygen documentation
 *
 * Revision 1.5  2006/12/10 19:00:54  sherwin
 * Modifications to handle nodal expansions
 *
 * Revision 1.4  2006/07/02 17:16:18  sherwin
 *
 * Modifications to make MultiRegions work for a connected domain in 2D (Tris)
 *
 * Revision 1.3  2006/06/01 14:13:36  kirby
 * *** empty log message ***
 *
 * Revision 1.2  2006/05/23 15:08:19  jfrazier
 * Minor cleanup to correct compile warnings.
 *
 * Revision 1.1  2006/05/04 18:58:31  kirby
 * *** empty log message ***
 *
 * Revision 1.30  2006/03/12 14:20:44  sherwin
 *
 * First compiling version of SpatialDomains and associated modifications
 *
 * Revision 1.29  2006/03/06 17:12:45  sherwin
 *
 * Updated to properly execute all current StdRegions Demos.
 *
 * Revision 1.28  2006/03/04 20:26:54  bnelson
 * Added comments after #endif.
 *
 * Revision 1.27  2006/03/01 08:25:04  sherwin
 *
 * First compiling version of StdRegions
 *
 * Revision 1.26  2006/02/27 23:47:23  sherwin
 *
 * Standard coding update upto compilation of StdHexExp.cpp
 *
 *
 **/



