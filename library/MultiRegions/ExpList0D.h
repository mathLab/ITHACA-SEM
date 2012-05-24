///////////////////////////////////////////////////////////////////////////////
//
// File ExpList0D.h
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
// Description: Expansion list 0D definition. This is not really a class descibing an expansion, but
// just a utilty class to manage boundary conditions for the 1D case and the homogenous cases.
// It basically represents a 0-dimensional expansion, or better a collection of points which are
// generally located on the boundaries. It is a wrap around LocalRegion::PointExp
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_MULTIREGIONS_EXPLIST0D_H
#define NEKTAR_LIB_MULTIREGIONS_EXPLIST0D_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <vector>
#include <fstream>

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ExpList.h>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/TriExp.h>
#include <LocalRegions/QuadExp.h>
#include <LocalRegions/PointExp.h>
#include <SpatialDomains/MeshGraph1D.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <LibUtilities/Kernel/kernel.h>
#include <SpatialDomains/MeshComponents.h>
//#include <LocalRegions/Expansion0D.h>


namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations for typedefs
        class ExpList0D;

        /// Shared pointer to an ExpList0D object.
        typedef boost::shared_ptr<ExpList0D>      ExpList0DSharedPtr;
        /// Vector of pointers to ExpList0D objects.
        typedef std::vector<ExpList0DSharedPtr>   ExpList0DVector;
        /// Iterator for the vector of ExpList0D pointers.
        typedef std::vector<ExpList0DSharedPtr>::iterator ExpList0DVectorIter;

        /// This class is the abstraction of a collection of
        /// zero-dimensional expansions which is merely a collection of points/values.
        class ExpList0D: public ExpList
        {
        public:
			
			/// The defualt constructor.
            MULTI_REGIONS_EXPORT ExpList0D();

            /// The copy constructor.
            MULTI_REGIONS_EXPORT ExpList0D(const ExpList0D &In, bool DeclareCoeffPhysArrays);

            // wrap around LocalRegion::PointExp
            MULTI_REGIONS_EXPORT ExpList0D(const SpatialDomains::VertexComponentSharedPtr &m_geom);

			/// Specialised constructor for trace expansions (croth)
            MULTI_REGIONS_EXPORT ExpList0D(
                const Array<OneD,const ExpListSharedPtr> &bndConstraint,
                const Array<OneD,const SpatialDomains
                            ::BoundaryConditionShPtr>  &bndCond,
                const StdRegions::StdExpansionVector &locexp,
                const SpatialDomains::MeshGraphSharedPtr &graph1D,
                const map<int,int> &periodicVertices,
                const bool DeclareCoeffPhysArrays = true);
            
            /// Destructor.
            MULTI_REGIONS_EXPORT virtual ~ExpList0D();
			
            LocalRegions::PointExpSharedPtr m_point;
            
            
        protected:
            virtual void v_Upwind(
                const Array<OneD, const NekDouble> &Vn,
                const Array<OneD, const NekDouble> &Fwd,
                const Array<OneD, const NekDouble> &Bwd,
                      Array<OneD,       NekDouble> &Upwind);
            
            virtual void v_GetNormals(
                Array<OneD, Array<OneD, NekDouble> > &normals);

            virtual void v_GetCoords(NekDouble &x, NekDouble &y, NekDouble &z);
            
            virtual void v_GetCoord(Array<OneD,NekDouble> &coords);
            
            virtual void v_SetCoeff(NekDouble val);
            
            virtual void v_SetPhys(NekDouble val);
            
            virtual const SpatialDomains::VertexComponentSharedPtr &v_GetGeom(void) const;
            
            virtual const SpatialDomains::VertexComponentSharedPtr &v_GetVertex(void) const;
            
        private:
            void SetCoeffPhysOffsets(void);
        };

        /// Empty ExpList0D object.
        const static Array<OneD, ExpList0DSharedPtr>
                                NullExpList0DSharedPtrArray;

    } //end of namespace
} //end of namespace

#endif//NEKTAR_LIB_MULTIREGIONS_EXPLIST0D_H

/**
 * $Log: v $
 *
 **/
