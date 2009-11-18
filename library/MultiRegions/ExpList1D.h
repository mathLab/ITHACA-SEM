///////////////////////////////////////////////////////////////////////////////
//
// File ExpList1D.h
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
// Description: Expansion list 1D definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_EXPLIST1D_H
#define NEKTAR_LIB_MULTIREGIONS_EXPLIST1D_H

#include <vector>
#include <fstream>

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/GenSegExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/PointExp.h>
#include <SpatialDomains/MeshGraph1D.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <LibUtilities/Kernel/kernel.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations for typedefs
        class ExpList1D;

        /// Shared pointer to an ExpList1D object.
        typedef boost::shared_ptr<ExpList1D>      ExpList1DSharedPtr;
        /// Vector of pointers to ExpList1D objects.
        typedef std::vector<ExpList1DSharedPtr>   ExpList1DVector;
        /// Iterator for the vector of ExpList1D pointers.
        typedef std::vector<ExpList1DSharedPtr>::iterator ExpList1DVectorIter;

        /// This class is the abstraction of a one-dimensional multi-elemental
        /// expansions which is merely a collection of local expansions.
        class ExpList1D: public ExpList
        {
        public:
            /// The default constructor.
            ExpList1D();

            /// The copy constructor.
            ExpList1D(const ExpList1D &In);

            /// Construct an ExpList1D from a given graph.
            ExpList1D(const LibUtilities::BasisKey &Ba,
                      const SpatialDomains::MeshGraph1D &graph1D,
                      bool UseGenSegExp = false);

            /// This constructor sets up a list of local expansions based on an
            /// input mesh.
            ExpList1D(SpatialDomains::MeshGraph1D &graph1D,
                      bool UseGenSegExp = false);

            /// Specialised constructor for Neumann boundary conditions in
            /// DisContField2D and ContField2D.
            ExpList1D(const SpatialDomains::CompositeVector &domain,
                      SpatialDomains::MeshGraph2D &graph2D,
                      bool UseGenSegExp = false);

            /// Specialised constructor for trace expansions.
            ExpList1D(const Array<OneD,const ExpList1DSharedPtr> &bndConstraint,
                      const Array<OneD,const SpatialDomains
                                            ::BoundaryConditionShPtr>  &bndCond,
                      const StdRegions::StdExpansionVector &locexp,
                      SpatialDomains::MeshGraph2D &graph2D,
                      const map<int,int> &periodicEdges,
                      bool UseGenSegExp = false);

            /// Destructor.
            ~ExpList1D();

            /// Populates the list of boundary condition expansions.
            void SetBoundaryConditionExpansion(
                                const SpatialDomains::MeshGraph1D &graph1D,
                                      SpatialDomains::BoundaryConditions &bcs,
                                const std::string variable,
                                Array<OneD, LocalRegions::PointExpSharedPtr>
                                                            &bndCondExpansions,
                                Array<OneD, SpatialDomains
                                    ::BoundaryConditionShPtr> &bndConditions);

            /// Generate a associative map of periodic vertices in a mesh.
            void GetPeriodicVertices(
                                const SpatialDomains::MeshGraph1D &graph1D,
                                      SpatialDomains::BoundaryConditions &bcs,
                                const std::string variable,
                                      map<int,int>& periodicVertices);

            // Evaluates boundary conditions.
            void EvaluateBoundaryConditions(
                                const NekDouble time,
                                Array<OneD, LocalRegions::PointExpSharedPtr>
                                                            &bndCondExpansions,
                                Array<OneD, SpatialDomains
                                    ::BoundaryConditionShPtr> &bndConditions);


            /// Performs the post-processing on a specified element.
            void PostProcess(   LibUtilities::KernelSharedPtr kernel,
                                Array<OneD,NekDouble> &inarray,
                                Array<OneD,NekDouble> &outarray,
                                NekDouble h,
                                int elmId = 0);

            /// Evaluates the global spectral/hp expansion at some arbitray set
            /// of points.
            void PeriodicEval(  Array<OneD,NekDouble> &inarray1,
                                Array<OneD,NekDouble> &inarray2,
                                NekDouble h, int nmodes,
                                Array<OneD,NekDouble> &outarray);

            /// Set up the normals on each expansion.
            void SetUpPhysNormals(const StdRegions::StdExpansionVector &locexp);

            // direction =  1: Upwind
            // direction = -1: Downwind

            /// Upwind the \a Fwd and \a Bwd states based on the velocity field
            /// given by \a Vec.
            void Upwind(const Array<OneD, const Array<OneD, NekDouble> > &Vec,
                        const Array<OneD, const NekDouble> &Fwd,
                        const Array<OneD, const NekDouble> &Bwd,
                        Array<OneD, NekDouble> &Upwind,
                        int direction=1);

            /// Upwind the \a Fwd and \a Bwd states based on the one-
            /// dimensional normal velocity field given by \a Vn.
            void Upwind(const Array<OneD, const NekDouble> &Vn,
                        const Array<OneD, const NekDouble> &Fwd,
                        const Array<OneD, const NekDouble> &Bwd,
                        Array<OneD, NekDouble> &Upwind,
                        int direction=1);

            /// Populate \a normals with the normals of all expansions.
            void GetNormals(Array<OneD, Array<OneD, NekDouble> > &normals);

        protected:
            /// Flag to indicate if GenSegExp's are being used.
            bool m_UseGenSegExp;

        private:
            /// Set up the normals on each expansion.
            virtual void v_SetUpPhysNormals(
                                const StdRegions::StdExpansionVector &locexp);
        };

        /// Empty ExpList1D object.
        const static Array<OneD, ExpList1DSharedPtr>
                                NullExpList1DSharedPtrArray;

    } //end of namespace
} //end of namespace

#endif//NEKTAR_LIB_MULTIREGIONS_EXPLIST1D_H

/**
 * Revision 1.23  2008/09/16 13:36:06  pvos
 * Restructured the LocalToGlobalMap classes
 *
 * Revision 1.22  2008/08/14 22:15:51  sherwin
 * Added LocalToglobalMap and DGMap and depracted LocalToGlobalBndryMap1D,2D. Made DisContField classes compatible with updated ContField formats
 *
 * Revision 1.21  2008/07/29 22:27:33  sherwin
 * Updates for DG solvers, including using GenSegExp, fixed forcing function on UDG HelmSolve and started to tidy up the mapping arrays to be 1D rather than 2D
 *
 * Revision 1.20  2008/06/23 14:21:01  pvos
 * updates for 1D ExpLists
 *
 * Revision 1.19  2008/05/10 18:27:33  sherwin
 * Modifications necessary for QuadExp Unified DG Solver
 *
 * Revision 1.18  2007/12/06 22:52:30  pvos
 * 2D Helmholtz solver updates
 *
 * Revision 1.17  2007/09/25 14:25:29  pvos
 * Update for helmholtz1D with different expansion orders
 *
 * Revision 1.16  2007/09/03 19:58:31  jfrazier
 * Formatting.
 *
 * Revision 1.15  2007/07/22 23:04:20  bnelson
 * Backed out Nektar::ptr.
 *
 * Revision 1.14  2007/07/20 02:04:12  bnelson
 * Replaced boost::shared_ptr with Nektar::ptr
 *
 * Revision 1.13  2007/07/10 08:54:30  pvos
 * Updated ContField1D constructor
 *
 * Revision 1.12  2007/07/06 18:39:34  pvos
 * ContField1D constructor updates
 *
 * Revision 1.11  2007/06/05 16:36:55  pvos
 * Updated Explist2D ContExpList2D and corresponding demo-codes
 *
 **/
