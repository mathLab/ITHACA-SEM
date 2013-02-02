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

#include <MultiRegions/MultiRegionsDeclspec.h>
#include <vector>
#include <MultiRegions/ExpList.h>
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
            MULTI_REGIONS_EXPORT ExpList1D();

            /// The copy constructor.
            MULTI_REGIONS_EXPORT ExpList1D(const ExpList1D &In, const bool DeclareCoeffPhysArrays = true);

            /// Construct an ExpList1D from a given graph.
            MULTI_REGIONS_EXPORT ExpList1D(
                      const LibUtilities::SessionReaderSharedPtr &pSession,
                      const LibUtilities::BasisKey &Ba,
                      const SpatialDomains::MeshGraphSharedPtr &graph1D);

            /// This constructor sets up a list of local expansions based on an
            /// input mesh.
            MULTI_REGIONS_EXPORT ExpList1D(
                      const LibUtilities::SessionReaderSharedPtr &pSession,
                      const SpatialDomains::MeshGraphSharedPtr &graph1D,
                      const bool DeclareCoeffPhysArrays = true);

            /// Specialised constructor for Neumann boundary conditions in
            /// DisContField2D and ContField2D.
            MULTI_REGIONS_EXPORT ExpList1D(
                      const SpatialDomains::CompositeMap &domain,
                      const SpatialDomains::MeshGraphSharedPtr &graph2D,
                      const bool DeclareCoeffPhysArrays = true,
                      const std::string variable = "DefaultVar");
			
			
            MULTI_REGIONS_EXPORT ExpList1D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                           const SpatialDomains::CompositeMap &domain,
                                           const SpatialDomains::MeshGraphSharedPtr &graph1D,
                                           int i,
                                           const bool DeclareCoeffPhysArrays = true);
            
            /// Specialised constructor for trace expansions.
            MULTI_REGIONS_EXPORT ExpList1D(
                      const Array<OneD,const ExpListSharedPtr> &bndConstraint,
                      const Array<OneD,const SpatialDomains
                                            ::BoundaryConditionShPtr>  &bndCond,
                      const StdRegions::StdExpansionVector &locexp,
                      const SpatialDomains::MeshGraphSharedPtr &graph2D,
                      const map<int,int> &periodicEdges,
                      const bool DeclareCoeffPhysArrays = true,
                      const std::string variable = "DefaultVar");
            

            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~ExpList1D();

            /// Performs the post-processing on a specified element.
            MULTI_REGIONS_EXPORT void PostProcess(   LibUtilities::KernelSharedPtr kernel,
                                Array<OneD,NekDouble> &inarray,
                                Array<OneD,NekDouble> &outarray,
                                NekDouble h,
                                int elmId = 0);

            /// Evaluates the global spectral/hp expansion at some arbitray set
            /// of points.
            MULTI_REGIONS_EXPORT void PeriodicEval(  Array<OneD,NekDouble> &inarray1,
                                Array<OneD,NekDouble> &inarray2,
                                NekDouble h, int nmodes,
                                Array<OneD,NekDouble> &outarray);

            /// Set up the normals on each expansion.
            //MULTI_REGIONS_EXPORT void SetUpPhysNormals();//const StdRegions::StdExpansionVector &locexp);

	    /// Set up the tangents on each expansion
	    MULTI_REGIONS_EXPORT void SetUpPhysTangents(const StdRegions::StdExpansionVector &locexp);
	    
            // direction =  1: Upwind
            // direction = -1: Downwind


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

        private:
            /// Definition of the total number of degrees of freedom and
            /// quadrature points. Sets up the storage for \a m_coeff and \a
            ///  m_phys.
            void SetCoeffPhysOffsets(void);
            
            virtual void v_ReadGlobalOptimizationParameters();
            
            /// Set up the normals on each expansion.
            virtual void v_SetUpPhysNormals();
            /// const StdRegions::StdExpansionVector &locexp);

            /// Set up the tangents on each expansion.
            virtual void v_SetUpPhysTangents(
            	    		const StdRegions::StdExpansionVector &locexp);
            
            virtual void v_WriteVtkPieceHeader(std::ofstream &outfile, int expansion);

            int m_firstIntEl;
            
            Array<OneD, NekDouble> m_normSign;
        };

        /// Empty ExpList1D object.
        const static Array<OneD, ExpList1DSharedPtr>
                                NullExpList1DSharedPtrArray;

    } //end of namespace
} //end of namespace

#endif//NEKTAR_LIB_MULTIREGIONS_EXPLIST1D_H

