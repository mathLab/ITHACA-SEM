///////////////////////////////////////////////////////////////////////////////
//
// File ExpList0D.cpp
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
// Description: Expansion list 0D definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList0D.h>

namespace Nektar
{
    namespace MultiRegions
    {

		/**
         * Default constructor ExpList0D object.
         */
        ExpList0D::ExpList0D():
		ExpList()
        {
        }

        /**
         * Creates an identical copy of another ExpList0D object.
         */
        ExpList0D::ExpList0D(const ExpList0D &In, bool DeclareCoeffPhysArrays):
		ExpList(In,DeclareCoeffPhysArrays)
        {
        }
		
        ExpList0D::ExpList0D(const SpatialDomains::VertexComponentSharedPtr &m_geom):
            ExpList()
        {
            m_point = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(m_geom);
            
            m_ncoeffs = 1;
            m_npoints = 1;
            
            // Set up m_coeffs, m_phys.
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
        }

        /**
         * Store expansions for the trace space expansions used in
         * DisContField1D.
         *
         * @param   bndConstraint   Array of ExpList1D objects each containing a
         *                      1D spectral/hp element expansion on a single
         *                      boundary region.
         * @param   bndCond     Array of BoundaryCondition objects which contain
         *                      information about the boundary conditions on the
         *                      different boundary regions.
         * @param   locexp      Complete domain expansion list.
         * @param   graph1D     1D mesh corresponding to the expansion list.
         * @param   periodicVertices   List of periodic Vertices.
         * @param   UseGenSegExp If true, create general segment expansions
         *                      instead of just normal segment expansions.
         */
        ExpList0D::ExpList0D(
            const Array<OneD, const ExpListSharedPtr> &bndConstraint,
            const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> 
                                                      &bndCond,
            const StdRegions::StdExpansionVector      &locexp,
            const SpatialDomains::MeshGraphSharedPtr  &graph1D,
            const map<int,int>                        &periodicVertices,
            const bool                                 DeclareCoeffPhysArrays)
            : ExpList()
        {
            int i, j, id, elmtid=0;
            map<int,int> EdgeDone;
            map<int,int> NormalSet;

            SpatialDomains::VertexComponentSharedPtr PointGeom;
            LocalRegions::PointExpSharedPtr Point;
			
            // First loop over boundary conditions to renumber Dirichlet boundaries
            for(i = 0; i < bndCond.num_elements(); ++i)
            {
                if(bndCond[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    for(j = 0; j < bndConstraint[i]->GetExpSize(); ++j)
                    {
						PointGeom = bndConstraint[i]->GetVertex();
                        Point = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(PointGeom);
                        
						EdgeDone[PointGeom->GetVid()] = elmtid;
						
                        Point->SetElmtId(elmtid++);
                        (*m_exp).push_back(Point);
                    }
                }
            }
			
            // loop over all other edges and fill out other connectivities
            for(i = 0; i < locexp.size(); ++i)
            {
                for(j = 0; j < 2; ++j)
                {
                    PointGeom = (locexp[i]->GetGeom1D())->GetVertex(j);
					id = PointGeom->GetVid();
					
                    if(EdgeDone.count(id)==0)
                    {						
                        Point = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(PointGeom);
                        EdgeDone[id] = elmtid;
						
                        //if (periodicVertices.count(id) > 0)
                        //{
						//   EdgeDone[periodicVertices.find(id)->second] = elmtid;
                        //}
						
                        Point->SetElmtId(elmtid++);
						(*m_exp).push_back(Point);
                    }
                    /*else // variable modes/points
                    {
                        LibUtilities::BasisKey EdgeBkey
						= locexp[i]->DetEdgeBasisKey(j);
						
                        if((*m_exp)[EdgeDone[id]]->GetNumPoints(0) >= EdgeBkey.GetNumPoints()
						   && (*m_exp)[EdgeDone[id]]->GetBasisNumModes(0) >= EdgeBkey.GetNumModes())
                        {
                        }
                        else if((*m_exp)[EdgeDone[id]]->GetNumPoints(0) <= EdgeBkey.GetNumPoints()
								&& (*m_exp)[EdgeDone[id]]->GetBasisNumModes(0) <= EdgeBkey.GetNumModes())
                        {
                            Seg = MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(EdgeBkey, SegGeom);
                            Seg->SetElmtId(EdgeDone[id]);
                            (*m_exp)[EdgeDone[id]] = Seg;
                            NormalSet.erase(id);
                        }
                        else
                        {
                            ASSERTL0(false,
									 "inappropriate number of points/modes (max "
									 "num of points is not set with max order)");
                        }
                    }*/
			 }
            }
		 
		 
			
            // Setup Default optimisation information.
            int nel = GetExpSize();
            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>::AllocateSharedPtr(nel);
			
            // Set up offset information and array sizes
            SetCoeffPhysOffsets();
			
            // Set up m_coeffs, m_phys.
            if(DeclareCoeffPhysArrays)
            {
                m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
                m_phys   = Array<OneD, NekDouble>(m_npoints);
            }
        }
		
		
		void ExpList0D::SetCoeffPhysOffsets()
        {
            int i;
			
            // Set up offset information and array sizes
            m_coeff_offset   = Array<OneD,int>(m_exp->size());
            m_phys_offset    = Array<OneD,int>(m_exp->size());
            m_offset_elmt_id = Array<OneD,int>(m_exp->size());
			
            m_ncoeffs = m_npoints = 0;
			
            for(i = 0; i < m_exp->size(); ++i)
            {
                m_coeff_offset[i]   = m_ncoeffs;
                m_phys_offset [i]   = m_npoints;
                m_offset_elmt_id[i] = i;
                m_ncoeffs += (*m_exp)[i]->GetNcoeffs();
                m_npoints += (*m_exp)[i]->GetTotPoints();
            }
        }

        /**
         *
         */
        ExpList0D::~ExpList0D()
        {
        }
		
        void ExpList0D::v_GetCoords(NekDouble &x, NekDouble &y, NekDouble &z)
        {
            m_point->GetCoords(x,y,z);
        }
		
        void ExpList0D::v_GetCoord(Array<OneD,NekDouble> &coords)
        {
            m_point->GetCoords(coords);
        }
		
        void ExpList0D::v_SetCoeff(NekDouble val)
        {
            m_coeffs[0] = val;
        }
		
        void ExpList0D::v_SetPhys(NekDouble val)
        {
            m_phys[0] = val;
        }
		
        const SpatialDomains::VertexComponentSharedPtr &ExpList0D::v_GetGeom(void) const
        {
            return m_point->GetGeom();
        }
		
        const SpatialDomains::VertexComponentSharedPtr &ExpList0D::v_GetVertex(void) const
        {
            return m_point->GetVertex();
        }
		
        /**
         * For each local element, copy the normals stored in the element list
         * into the array \a normals.
         * @param   normals     Multidimensional array in which to copy normals
         *                      to. Must have dimension equal to or larger than
         *                      the spatial dimension of the elements.
         */
        void ExpList0D::v_GetNormals(Array<OneD, Array<OneD, NekDouble> > &normals)
        {
			int i,j,k,e_npoints,offset;
            Array<OneD,Array<OneD,NekDouble> > locnormals;
			
            // Assume whole array is of same coordinate dimension
			int coordim = normals.num_elements();
			
            ASSERTL1(normals.num_elements() >= coordim,
                     "Output vector does not have sufficient dimensions to "
                     "match coordim");
			
            // Process each expansion.
            for(i = 0; i < m_exp->size(); ++i)
            {
                LocalRegions::Expansion0DSharedPtr loc_exp = boost::dynamic_pointer_cast<LocalRegions::Expansion0D>((*m_exp)[i]);

                LocalRegions::Expansion1DSharedPtr loc_elmt = loc_exp->GetLeftAdjacentElementExp();
                
                // Get the number of points and normals for this expansion.
                e_npoints  = 1;
                locnormals = loc_elmt->GetVertexNormal(loc_exp->GetLeftAdjacentElementVertex());
				
                // Get the physical data offset for this expansion.
                offset = m_phys_offset[i];
				
                // Process each point in the expansion.
                for(j = 0; j < e_npoints; ++j)
                {
                    // Process each spatial dimension and copy the values into
                    // the output array.
                    for(k = 0; k < coordim; ++k)
                    {
                        normals[k][offset] = locnormals[k][0];
                    }
                }
            }
            //negate first normal inwards facing into domain
            for (k=0; k<coordim; k++)
            {
                normals[k][0] = locnormals[k][0];
            }
	
        }
		
        /**
         * One-dimensional upwind.
         * 
         * @param   Vn          Velocity field.
         * @param   Fwd         Left state.
         * @param   Bwd         Right state.
         * @param   Upwind      Output vector.
         */
        void ExpList0D::v_Upwind(const Array<OneD, const NekDouble> &Vn,
                                 const Array<OneD, const NekDouble> &Fwd,
                                 const Array<OneD, const NekDouble> &Bwd,
                                       Array<OneD,       NekDouble> &Upwind)
        {
            // Process each point in the expansion.
            for(int j = 0; j < Fwd.num_elements(); ++j)
            {
                // Upwind based on one-dimensional velocity.
                if(Vn[j] > 0.0)
                {
                    Upwind[j] = Fwd[j];
                }
                else
                {
                    Upwind[j] = Bwd[j];
                }
            }
        }
    } //end of namespace
} //end of namespace
