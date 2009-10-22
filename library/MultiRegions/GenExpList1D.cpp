///////////////////////////////////////////////////////////////////////////////
//
// File GenExpList1D.cpp
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
// Description: Generalised Expansion list 1D definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/GenExpList1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
    
        GenExpList1D::GenExpList1D():
            ExpList1D()
        {
        }
        
        GenExpList1D::~GenExpList1D()
        {
        }
        
        GenExpList1D::GenExpList1D(const GenExpList1D &In):
            ExpList1D(In)
        {
        }

        // specialised constructor for Neumann boundary conditions in
        // DisContField2D and ContField2d
        GenExpList1D::GenExpList1D(const SpatialDomains::CompositeVector &domain, SpatialDomains::MeshGraph2D &graph2D)
        {
            int i,j, nel,cnt,id=0;
            SpatialDomains::Composite comp;
            SpatialDomains::SegGeomSharedPtr SegmentGeom;
            LocalRegions::SegExpSharedPtr seg;
            
            nel = 0;
            for(i = 0; i < domain.size(); ++i)
            {
                nel += (domain[i])->size();
            }

            m_coeff_offset = Array<OneD,int>(nel,0);
            m_phys_offset  = Array<OneD,int>(nel,0);

            cnt = 0;
            for(i = 0; i < domain.size(); ++i)
            {
                comp = domain[i];
                
                for(j = 0; j < comp->size(); ++j)
                {
                    if(SegmentGeom = boost::dynamic_pointer_cast<SpatialDomains::SegGeom>((*comp)[j]))
                    {
                        LibUtilities::BasisKey bkey = graph2D.GetEdgeBasisKey(SegmentGeom);
                        seg = MemoryManager<LocalRegions::GenSegExp>::AllocateSharedPtr(bkey, SegmentGeom);
                        
                        seg->SetElmtId(id++);
                        (*m_exp).push_back(seg);  

                        m_coeff_offset[cnt] = m_ncoeffs; 
                        m_phys_offset[cnt]  = m_npoints; cnt++;
                        m_ncoeffs += bkey.GetNumModes();
                        m_npoints += bkey.GetNumPoints();
                     }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a SegGeom failed");
                    }  
                }
                
            }

            ExpList::SetCoeffPhys();
        }

        // Specialised constructor for DisContField2D for trace expansion

        GenExpList1D::GenExpList1D(const Array<OneD,const MultiRegions::ExpList1DSharedPtr> &bndConstraint,  
                                   const Array<OneD, const SpatialDomains::BoundaryConditionShPtr>  &bndCond, 
                                   const StdRegions::StdExpansionVector &locexp, 
                                   SpatialDomains::MeshGraph2D &graph2D,
                                   const map<int,int> &periodicEdges)
        {
            int i,j,id, elmtid=0, edgeid;
            // Could replace these with maps
            map<int, int> EdgeDone;
            map<int, int> NormalSet;

            StdRegions::EdgeOrientation      orient;
            LocalRegions::GenSegExpSharedPtr Seg;
            SpatialDomains::ElementEdgeVectorSharedPtr con_elmt;
            SpatialDomains::Geometry2DSharedPtr AdjGeom;
            SpatialDomains::Geometry1DSharedPtr SegGeom;

            // First loop over boundary conditions to renumber
            // Dirichlet boundaries to be first 

            for(i = 0; i < bndCond.num_elements(); ++i)
	      {
                if(bndCond[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    for(j = 0; j < bndConstraint[i]->GetExpSize(); ++j)
                    {
                        LibUtilities::BasisKey bkey = bndConstraint[i]->GetExp(j)->GetBasis(0)->GetBasisKey();
                        SegGeom  = bndConstraint[i]->GetExp(j)->GetGeom1D();  
                        Seg = MemoryManager<LocalRegions::GenSegExp>::AllocateSharedPtr(bkey, SegGeom);
                        con_elmt = graph2D.GetElementsFromEdge(SegGeom);
                        
                        EdgeDone[SegGeom->GetEid()] = elmtid;
                        
                        Seg->SetElmtId(elmtid++);
                        (*m_exp).push_back(Seg);
                    }
                }
            }

            // loop over all other edges and fill out other connectivities	    
            for(i = 0; i < locexp.size(); ++i)
            {
                for(j = 0; j < locexp[i]->GetNedges(); ++j)
                {   
                    const SpatialDomains::Geometry1DSharedPtr& SegGeom = (locexp[i]->GetGeom2D())->GetEdge(j);
                    
                    id = SegGeom->GetEid();

                    if(EdgeDone.count(id)==0)
		      {
			
                        LibUtilities::BasisKey EdgeBkey = locexp[i]->DetEdgeBasisKey(j);
			
                        Seg = MemoryManager<LocalRegions::GenSegExp>::AllocateSharedPtr(EdgeBkey, SegGeom);

                        EdgeDone[id] = elmtid;
			
                        // set periodic edge 
                        if(periodicEdges.count(id) > 0)
			  {
                            EdgeDone[periodicEdges.find(id)->second] = elmtid;
			  }
			
                        Seg->SetElmtId(elmtid++);
                        
                        (*m_exp).push_back(Seg);
		      }
		    else // variable modes/points
		      {
			LibUtilities::BasisKey EdgeBkey = locexp[i]->DetEdgeBasisKey(j);
			
			if((*m_exp)[EdgeDone[id]]->GetNumPoints(0) >= EdgeBkey.GetNumPoints() && (*m_exp)[EdgeDone[id]]->GetBasisNumModes(0) >= EdgeBkey.GetNumModes())
			  {
			  }
			else if((*m_exp)[EdgeDone[id]]->GetNumPoints(0) <= EdgeBkey.GetNumPoints() && (*m_exp)[EdgeDone[id]]->GetBasisNumModes(0) <= EdgeBkey.GetNumModes())
			  {
			    Seg = MemoryManager<LocalRegions::GenSegExp>::AllocateSharedPtr(EdgeBkey, SegGeom);
			    Seg->SetElmtId(EdgeDone[id]);
			    (*m_exp)[EdgeDone[id]] = Seg;
			    NormalSet.erase(id);          
			  }
			else
			  {
			    ASSERTL0(false,"inappropriate number of points/modes (max num of points is not set with max order)")
			  }
		      }
		    
                    if(NormalSet.count(id) == 0)
                    {
                        Seg = boost::dynamic_pointer_cast<LocalRegions::GenSegExp>((*m_exp)[EdgeDone.find(id)->second]);
                        
                        // Set up normals at all Segment Quadrature points
                        Seg->SetUpPhysNormals(locexp[i],j);
                        
                        NormalSet[id] = 1;
                    }
                }
            }
            ExpList::SetCoeffPhys();
        }

        void GenExpList1D::SetUpPhysNormals(const StdRegions::StdExpansionVector &locexp)
        {
            map<int, int> EdgeGID;
            int i,j,cnt,n,id; 

            // setup map of all global ids along boundary
            for(cnt = i = 0; i < (*m_exp).size(); ++i)
            {
                id =  (*m_exp)[i]->GetGeom1D()->GetEid();
                EdgeGID[id] = cnt++;
            }
            
            // Loop over elements and find edges that match;
            for(cnt = n = 0; n < locexp.size(); ++n)
            {
                for(i = 0; i < locexp[n]->GetNedges(); ++i)
                {
                    id = locexp[n]->GetGeom2D()->GetEid(i);
                    
                    if(EdgeGID.count(id) > 0)
                    {
                        (*m_exp)[EdgeGID.find(id)->second]->SetUpPhysNormals(locexp[n],i);
                    }
                }
            }                
        }        

        // Upwind the left and right states given by the Arrays Fwd
        // and Bwd using the vector quantity Vec and ouput the
        // upwinded value in the array upwind
        void GenExpList1D::Upwind(const Array<OneD, const Array<OneD, NekDouble> > &Vec,
                                  const Array<OneD, const NekDouble> &Fwd, 
                                  const Array<OneD, const NekDouble> &Bwd, 
                                  Array<OneD, NekDouble> &Upwind,
                                  int direction)
        {
            int i,j,k,e_npoints,offset;
            Array<OneD,NekDouble> normals; 
            NekDouble Vn;

            // Assume whole array is of same coordimate dimention
            int coordim = (*m_exp)[0]->GetGeom1D()->GetCoordim();
            
            ASSERTL1(Vec.num_elements() >= coordim,
              "Input vector does not have sufficient dimensions to match coordim");

            for(i = 0; i < m_exp->size(); ++i)
            {
                e_npoints = (*m_exp)[i]->GetNumPoints(0);
		 
                normals   = (*m_exp)[i]->GetPhysNormals();
                
                offset = m_phys_offset[i];

                for(j = 0; j < e_npoints; ++j)
                {
                    // Calculate normal velocity
                    Vn = 0.0;
                    for(k = 0; k < coordim; ++k)
                    {
                        Vn += Vec[k][offset+j]*normals[k*e_npoints + j];
                    }

                    // Upwind
                    if(Vn > 0.0)
                    {
                        Upwind[offset + j] = Fwd[offset + j];
                    }
                    else
                    {
                        Upwind[offset + j] = Bwd[offset + j];
                    }
                }
            }            
        }

        void GenExpList1D::Upwind(const Array<OneD, const NekDouble> &Vn, 
                                  const Array<OneD, const NekDouble> &Fwd, 
                                  const Array<OneD, const NekDouble> &Bwd, 
                                  Array<OneD, NekDouble> &Upwind,
                                  int direction)
        {
            int i,j,k,e_npoints,offset;
            Array<OneD,NekDouble> normals; 

            // Assume whole array is of same coordimate dimention
            int coordim = (*m_exp)[0]->GetGeom1D()->GetCoordim();
            
            for(i = 0; i < m_exp->size(); ++i)
            {
                e_npoints = (*m_exp)[i]->GetNumPoints(0);
                offset = m_phys_offset[i];

                for(j = 0; j < e_npoints; ++j)
                {
                    // Upwind
                    if(Vn[offset + j] > 0.0)
                    {
                        Upwind[offset + j] = Fwd[offset + j];
                    }
                    else
                    {
                        Upwind[offset + j] = Bwd[offset + j];
                    }
                }
            }            
        }
		

        void GenExpList1D::GetNormals(Array<OneD, Array<OneD, NekDouble> > &normals) 
        {
            
            int i,j,k,e_npoints,offset;
            Array<OneD,NekDouble> locnormals; 
            
            // Assume whole array is of same coordinate dimension
            int coordim = (*m_exp)[0]->GetGeom1D()->GetCoordim();
            
            ASSERTL1(normals.num_elements() >= coordim,
                     "Output vector does not have sufficient dimensions to match coordim");
            
            for(i = 0; i < m_exp->size(); ++i)
            {
                e_npoints = (*m_exp)[i]->GetNumPoints(0);
                locnormals   = (*m_exp)[i]->GetPhysNormals();
                
                offset = m_phys_offset[i];
                
                for(j = 0; j < e_npoints; ++j)
                {
                    for(k = 0; k < coordim; ++k)
                    {
                        normals[k][offset+j] = locnormals[k*e_npoints + j];
                    }
                }
            }
        }
        
    } //end of namespace
} //end of namespace

/**
* $Log: GenExpList1D.cpp,v $
* Revision 1.13  2009/09/08 14:56:47  rcantao
* - Fixed an extra class qualifier inside *.h
* - Fixed missing parameter direction on Upwind methods
*
* Revision 1.12  2009/09/06 22:28:45  sherwin
* Updates for Navier-Stokes solver
*
* Revision 1.11  2009/07/09 09:01:49  sehunchun
* Upwind function is modified in a faster form
*
* Revision 1.10  2009/04/20 16:14:06  sherwin
* Updates for optimising bandwidth of DG solver and allowing write import on explist
*
* Revision 1.9  2009/02/28 21:28:40  sehunchun
*  Now upwind has "forward" direction and "backward" direction. Default is forward and no changes are necessary for previous file.
*
* Revision 1.8  2008/10/29 22:46:35  sherwin
* Updates for const correctness and a few other bits
*
* Revision 1.7  2008/10/09 21:47:36  ehan
* Fixed error from the function Upwind().
*
* Revision 1.6  2008/10/04 19:54:28  sherwin
* Added upwind method using only the normal flux
*
* Revision 1.5  2008/09/09 15:06:03  sherwin
* Modifications related to curved elements.
*
* Revision 1.4  2008/08/22 09:42:32  pvos
* Updates for Claes' Shallow Water and Euler solver
*
* Revision 1.3  2008/08/14 22:15:51  sherwin
* Added LocalToglobalMap and DGMap and depracted LocalToGlobalBndryMap1D,2D. Made DisContField classes compatible with updated ContField formats
*
* Revision 1.2  2008/07/31 11:17:13  sherwin
* Changed GetEdgeBasis with DetEdgeBasisKey
*
* Revision 1.1  2008/07/29 22:26:35  sherwin
* Generalised 1D Segment list which includes a normal direction at physical points
*
**/
