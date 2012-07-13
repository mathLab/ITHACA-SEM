///////////////////////////////////////////////////////////////////////////////
//
// File ExpList2D.h
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
// Description: Expansion list 2D header definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef EXPLIST2D_H
#define EXPLIST2D_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <vector>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <SpatialDomains/MeshGraph3D.h>
#include <SpatialDomains/Conditions.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declaration for typedefs
        class ExpList2D;

        /// Shared pointer to an ExpList2D object.
        typedef boost::shared_ptr<ExpList2D>      ExpList2DSharedPtr;
        /// Vector of pointers to ExpList2D objects.
        typedef std::vector< ExpList2DSharedPtr > ExpList2DVector;
        /// Iterator for the vector of ExpList2D pointers.
        typedef std::vector< ExpList2DSharedPtr >::iterator ExpList2DVectorIter;

        typedef pair<int, StdRegions::Orientation> PeriodicFace;
        
        /// Abstraction of a two-dimensional multi-elemental expansion which
        /// is merely a collection of local expansions.
        class ExpList2D: public ExpList
        {
        public:
            /// Default constructor.
            MULTI_REGIONS_EXPORT ExpList2D();

            /// Copy constructor.
            MULTI_REGIONS_EXPORT ExpList2D(  const ExpList2D &In,
                        const bool DeclareCoeffPhysArrays = true);

            /// Sets up a list of local expansions based on an input mesh.
            MULTI_REGIONS_EXPORT ExpList2D(
                      const LibUtilities::SessionReaderSharedPtr &pSession,
                      const SpatialDomains::MeshGraphSharedPtr &graph2D,
                      const bool DelcareCoeffPhysArrays = true,
                      const std::string &var = "DefaultVar");

            /// Sets up a list of local expansions based on an expansion  Mapr
            MULTI_REGIONS_EXPORT ExpList2D(
                      const LibUtilities::SessionReaderSharedPtr &pSession,
                      const SpatialDomains::ExpansionMap &expansions,
                      const bool DeclareCoeffPhysArrays = true);
            
            /// Sets up a list of local expansions based on an input mesh
            /// and separately defined basiskeys
            MULTI_REGIONS_EXPORT ExpList2D(
                      const LibUtilities::SessionReaderSharedPtr &pSession,
                      const LibUtilities::BasisKey &TriBa,
                      const LibUtilities::BasisKey &TriBb,
                      const LibUtilities::BasisKey &QuadBa,
                      const LibUtilities::BasisKey &QuadBb,
                      const SpatialDomains::MeshGraphSharedPtr &graph2D,
                      const LibUtilities::PointsType
                      TriNb = LibUtilities::SIZE_PointsType);

			//Specialized constructor for trace expansions.
            MULTI_REGIONS_EXPORT ExpList2D(
                      const Array<OneD,const ExpListSharedPtr> &bndConstraint,
                      const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>  &bndCond,
                      const StdRegions::StdExpansionVector &locexp,
                      const SpatialDomains::MeshGraphSharedPtr &graph3D,
                      const map<int,PeriodicFace> &periodicFaces,
                      const bool DeclareCoeffPhysArrays = true);

            /// Specialised constructor for Neumann boundary conditions in
            /// DisContField3D and ContField3D.
            MULTI_REGIONS_EXPORT ExpList2D(  
                        const LibUtilities::SessionReaderSharedPtr &pSession,
                        const SpatialDomains::CompositeMap &domain,
                        const SpatialDomains::MeshGraphSharedPtr &graph3D);
            
            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~ExpList2D();

        protected:
            /// Upwind the \a Fwd and \a Bwd states based on the one-
            /// dimensional normal velocity field given by \a Vn.
            MULTI_REGIONS_EXPORT void v_Upwind(
                const Array<OneD, const NekDouble> &Vn,
                const Array<OneD, const NekDouble> &Fwd,
                const Array<OneD, const NekDouble> &Bwd,
                      Array<OneD,       NekDouble> &Upwind);
            
            void v_GetNormals(
                Array<OneD, Array<OneD, NekDouble> > &normals);

        private:
            /// Definition of the total number of degrees of freedom and
            /// quadrature points and offsets to access datax
            void SetCoeffPhysOffsets(void);

            /// Set up the normals on each expansion.
            virtual void v_SetUpPhysNormals();

            virtual void v_ReadGlobalOptimizationParameters();

            virtual void v_WriteVtkPieceHeader(std::ofstream &outfile, int expansion);

        };

        /// Empty ExpList2D object.
        const static Array<OneD, ExpList2DSharedPtr>
                                                NullExpList2DSharedPtrArray;
    } //end of namespace
} //end of namespace

#endif//EXPLIST2D_H

/**
* $Log: ExpList2D.h,v $
* Revision 1.24  2009/12/15 18:09:03  cantwell
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
* Revision 1.23  2009/11/25 14:51:29  pvos
* Updates for added Timings directory
*
* Revision 1.22  2009/11/19 23:30:36  cantwell
* Documentation for ExpList2D and GlobalMatrixKey
* Updated doxygen pages.
*
* Revision 1.21  2009/11/02 19:15:43  cantwell
* Moved ContField1D to inherit from DisContField1D.
* Moved ContField3D to inherit from DisContField3D.
* Incorporated GenExpList1D functionality into ExpList1D.
* Tidied up and added documentation to various classes.
* Moved Namespace documentation and introductions to separate files along with
* doxygen configuration.
* Added option to use system ZLIB library instead of libboost_zlib on UNIX.
* Added extra search paths to FindMetis.cmake and FindNektar++.cmake.
* Updated Linux compiling instructions.
* Updated regDemo to use Helmholtz2D-g when built as debug.
*
* Revision 1.20  2009/09/06 22:28:45  sherwin
* Updates for Navier-Stokes solver
*
* Revision 1.19  2009/05/10 23:17:12  sherwin
* Updated mainly to handle doubly periodic meshes which required modification to vertex handling from a numbering perspective
*
* Revision 1.18  2009/03/04 14:17:38  pvos
* Removed all methods that take and Expansion as argument
*
* Revision 1.17  2008/10/04 20:04:26  sherwin
* Modifications for solver access
*
* Revision 1.16  2008/09/23 18:21:00  pvos
* Updates for working ProjectContField3D demo
*
* Revision 1.15  2008/09/17 13:46:40  pvos
* Added AssemblyMapCG for 3D expansions
*
* Revision 1.14  2008/08/14 22:15:51  sherwin
* Added LocalToglobalMap and DGMap and depracted LocalToGlobalBndryMap1D,2D. Made DisContField classes compatible with updated ContField formats
*
* Revision 1.13  2008/06/05 15:06:58  pvos
* Added documentation
*
* Revision 1.12  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
*
* Revision 1.11  2007/07/22 23:04:21  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.10  2007/07/20 02:04:12  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.9  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/

