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
#include <MultiRegions/ExpList.h>
#include <SpatialDomains/Conditions.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declaration for typedefs
        class ExpList2D;

        /// Shared pointer to an ExpList2D object.
        typedef std::shared_ptr<ExpList2D>      ExpList2DSharedPtr;
        /// Vector of pointers to ExpList2D objects.
        typedef std::vector< ExpList2DSharedPtr > ExpList2DVector;

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
            
            /// Constructor copying only elements defined in eIds.
            MULTI_REGIONS_EXPORT ExpList2D(  const ExpList2D &In,
                const std::vector<unsigned int> &eIDs,
                const bool DeclareCoeffPhysArrays = true,
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);

            /// Sets up a list of local expansions based on an input mesh.
            MULTI_REGIONS_EXPORT ExpList2D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const bool DelcareCoeffPhysArrays = true,
                const std::string &var = "DefaultVar",
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);

            /// Sets up a list of local expansions based on an expansion  Map
            MULTI_REGIONS_EXPORT ExpList2D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::ExpansionMap &expansions,
                const bool DeclareCoeffPhysArrays = true,
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);
            
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
                TriNb = LibUtilities::SIZE_PointsType,
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);

			//Specialized constructor for trace expansions.
            MULTI_REGIONS_EXPORT ExpList2D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const Array<OneD,const ExpListSharedPtr> &bndConstraint,
                const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
                                                                       &bndCond,
                const LocalRegions::ExpansionVector &locexp,
                const SpatialDomains::MeshGraphSharedPtr &graph3D,
                const PeriodicMap &periodicFaces,
                const bool DeclareCoeffPhysArrays = true,
                const std::string variable = "DefaultVar",
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);

            /// Specialised constructor for Neumann boundary conditions in
            /// DisContField3D and ContField3D.
            MULTI_REGIONS_EXPORT ExpList2D(  
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::CompositeMap &domain,
                const SpatialDomains::MeshGraphSharedPtr &graph3D,
                const std::string variable = "DefaultVar",
                const LibUtilities::CommSharedPtr comm =
                                                  LibUtilities::CommSharedPtr(),
                const Collections::ImplementationType ImpType
                                                     = Collections::eNoImpType);

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

            virtual void v_GetElmtNormalLength(
                Array<OneD, NekDouble>  &lengthsFwd,
                Array<OneD, NekDouble>  &lengthsBwd);
            
        private:
            /// Set up the normals on each expansion.
            virtual void v_SetUpPhysNormals();

            virtual void v_WriteVtkPieceHeader(
                        std::ostream &outfile, int expansion, int istrip);

            virtual void v_PhysInterp1DScaled(
                const NekDouble scale,
                const Array<OneD, NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray);

            virtual void v_PhysGalerkinProjection1DScaled(
                const NekDouble scale,
                const Array<OneD, NekDouble> &inarray,
                      Array<OneD, NekDouble> &outarray);
        };

        /// Empty ExpList2D object.
        const static Array<OneD, ExpList2DSharedPtr>
            NullExpList2DSharedPtrArray;
    } //end of namespace
} //end of namespace

#endif//EXPLIST2D_H

