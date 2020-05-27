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

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <vector>
#include <MultiRegions/ExpList.h>
#include <LocalRegions/SegExp.h>
#include <LibUtilities/Kernel/kernel.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations for typedefs
        class ExpList1D;

        /// Shared pointer to an ExpList1D object.
        typedef std::shared_ptr<ExpList1D>      ExpList1DSharedPtr;
        /// Vector of pointers to ExpList1D objects.
        typedef std::vector<ExpList1DSharedPtr>   ExpList1DVector;

        /// This class is the abstraction of a one-dimensional multi-elemental
        /// expansions which is merely a collection of local expansions.
        class ExpList1D: public ExpList
        {
        public:
            /// The default constructor.
            MULTI_REGIONS_EXPORT ExpList1D();

            /// The copy constructor.
            MULTI_REGIONS_EXPORT ExpList1D(
                const ExpList1D &In,
                const bool DeclareCoeffPhysArrays = true);
                
            
            /// Constructor copying only elements defined in eIds.
            MULTI_REGIONS_EXPORT ExpList1D(  const ExpList1D &In,
                const std::vector<unsigned int> &eIDs,
                const bool DeclareCoeffPhysArrays = true,
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);

            /// Construct an ExpList1D from a given graph.
            MULTI_REGIONS_EXPORT ExpList1D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const LibUtilities::BasisKey &Ba,
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);

            /// This constructor sets up a list of local expansions based on an
            /// input graph1D.
            MULTI_REGIONS_EXPORT ExpList1D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const bool DeclareCoeffPhysArrays = true,
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);


            /// This constructor sets up a list of local expansions based on an
            /// input compositeMap
            MULTI_REGIONS_EXPORT ExpList1D(
                      const LibUtilities::SessionReaderSharedPtr &pSession,
                      const SpatialDomains::MeshGraphSharedPtr &graph1D,
                      const SpatialDomains::CompositeMap &domain,
                      const bool DeclareCoeffPhysArrays = true,
                      const std::string var = "DefaultVar",
                      bool SetToOneSpaceDimension = false,
                      const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);



            /// Specialised constructor for Neumann boundary conditions in
            /// DisContField2D and ContField2D.
            MULTI_REGIONS_EXPORT ExpList1D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::CompositeMap &domain,
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const bool DeclareCoeffPhysArrays = true,
                const std::string variable = "DefaultVar",
                const LibUtilities::CommSharedPtr comm
                                                = LibUtilities::CommSharedPtr(),
                const Collections::ImplementationType ImpType
                                                     = Collections::eNoImpType);

            MULTI_REGIONS_EXPORT ExpList1D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const SpatialDomains::CompositeMap &domain,
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                int i,
                const bool DeclareCoeffPhysArrays = true,
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);  
           
            /// Specialised constructor for trace expansions.
            MULTI_REGIONS_EXPORT ExpList1D(
                const LibUtilities::SessionReaderSharedPtr &pSession,
                const Array<OneD,const ExpListSharedPtr> &bndConstraint,
                const Array<OneD,const SpatialDomains
                                            ::BoundaryConditionShPtr>  &bndCond,
                const LocalRegions::ExpansionVector &locexp,
                const SpatialDomains::MeshGraphSharedPtr &graph2D,
                const PeriodicMap &periodicEdges,
                const bool DeclareCoeffPhysArrays = true,
                const std::string variable = "DefaultVar",
                const Collections::ImplementationType ImpType
                                             = Collections::eNoImpType);
                

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~ExpList1D();

            /// Performs the post-processing on a specified element.
            MULTI_REGIONS_EXPORT void PostProcess(
                LibUtilities::KernelSharedPtr kernel,
                Array<OneD,NekDouble> &inarray,
                Array<OneD,NekDouble> &outarray,
                NekDouble h,
                int elmId = 0);

            /// Evaluates the global spectral/hp expansion at some arbitray set
            /// of points.
            MULTI_REGIONS_EXPORT void PeriodicEval(
                Array<OneD,NekDouble> &inarray1,
                Array<OneD,NekDouble> &inarray2,
                NekDouble h, int nmodes,
                Array<OneD,NekDouble> &outarray);

            /// Set up the normals on each expansion.
            //MULTI_REGIONS_EXPORT void SetUpPhysNormals();
            //const StdRegions::StdExpansionVector &locexp);

        MULTI_REGIONS_EXPORT void ParNormalSign(Array<OneD, NekDouble> &normsign);

        protected:
            /// Upwind the \a Fwd and \a Bwd states based on the velocity
            /// field given by \a Vec.
            void v_Upwind(
                const Array<OneD, const Array<OneD,       NekDouble> > &Vec,
                const Array<OneD,                   const NekDouble>   &Fwd,
                const Array<OneD,                   const NekDouble>   &Bwd,
                      Array<OneD,                         NekDouble>   &Upwind);

            /// Upwind the \a Fwd and \a Bwd states based on the one-
            /// dimensional normal velocity field given by \a Vn.
            void v_Upwind(
                const Array<OneD, const NekDouble> &Vn,
                const Array<OneD, const NekDouble> &Fwd,
                const Array<OneD, const NekDouble> &Bwd,
                      Array<OneD,       NekDouble> &Upwind);

            /// Populate \a normals with the normals of all expansions.
            void v_GetNormals(Array<OneD, Array<OneD, NekDouble> > &normals);

            virtual void v_GetElmtNormalLength(
                Array<OneD, NekDouble>  &lengthsFwd,
                Array<OneD, NekDouble>  &lengthsBwd);

        private:
            /// Set up the normals on each expansion.
            virtual void v_SetUpPhysNormals();
            /// const StdRegions::StdExpansionVector &locexp);

            virtual void v_WriteVtkPieceHeader(std::ostream &outfile, int expansion, int istrip);

        };

        /// Empty ExpList1D object.
        const static Array<OneD, ExpList1DSharedPtr>
                                NullExpList1DSharedPtrArray;

    } //end of namespace
} //end of namespace

#endif//NEKTAR_LIB_MULTIREGIONS_EXPLIST1D_H

