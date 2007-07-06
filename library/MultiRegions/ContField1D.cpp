///////////////////////////////////////////////////////////////////////////////
//
// File ContField1D.cpp
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
// Description: Field definition for 1D domain with boundary conditions
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContField1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
	
	ContField1D::ContField1D(void)
	{
	}

        ContField1D::ContField1D(const LibUtilities::BasisKey &Ba, 
                                 const SpatialDomains::Domain &domain1D,
                                 SpatialDomains::BoundaryConditions &bcs):
            ContExpList1D(Ba,domain1D)
        {
	    int i,nbnd;
            LocalRegions::PointExpSharedPtr  p_exp;
            int NumDirichet = 0;

            SpatialDomains::BoundaryRegionCollectionType    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollectionType &bconditions = bcs.GetBoundaryConditions();

            nbnd = bregions.size();

            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {
                if(  ((*(bconditions[i]))["u"])->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    SpatialDomains::VertexComponentSharedPtr vert;
                    
                    if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[0])[0]))
                    {
                        p_exp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                        m_bndConstraint.push_back(p_exp);
                        m_bndTypes.push_back(SpatialDomains::eDirichlet);
                        NumDirichet++;
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a vertex faiiled");
                    }
                }
            }

            for(i = 0; i < nbnd; ++i)
            {
                if( ((*(bconditions[i]))["u"])->GetBoundaryConditionType() != SpatialDomains::eDirichlet)
                {
                   SpatialDomains:: VertexComponentSharedPtr vert;

                    if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[0])[0]))
                    {
                        p_exp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                        m_bndConstraint.push_back(p_exp);
                        m_bndTypes.push_back(SpatialDomains::eDirichlet);
                        NumDirichet++;
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a vertex faiiled");
                    }
                }
            }
            

            // Need to reset numbering according to Dirichlet Bondary Condition
            //m_locToGloMap->ResetMapping(NumDirichlet,bcs);
            // Need to know how to get global vertex id 
            // I note that boundary Regions is a composite vector and so belive we can use this to get the desired quantities. 

        }

        ContField1D::~ContField1D()
	{
	}

	void ContField1D::FwdTrans(const ExpList &In)
	{
            int i, NumDirBcs;
            int NumDirBCs = m_locToGloMap->GetNumDirichletBCs();
            Array<OneD,NekDouble> tmp;

	    IProductWRTBase(In);

	    if(!(m_mass.get()))
	    {
		GenMassMatrixLinSys();
	    }

            for(i = 0; i < NumDirBcs; ++i)
            {
                m_contCoeffs[i] = m_bndConstraint[i]->GetValue();
            }
            
            // need to modify right hand side vector here to take
            // account of Dirichlet contribution to element


	    DNekVec v(m_contNcoeffs-NumDirBcs, tmp = m_contCoeffs+NumDirBcs,eWrapper);
	    m_mass->Solve(v,v);
	    m_transState = eContinuous;
	    m_physState = false;
	}

	void ContField1D::GenMassMatrixLinSys(void)
	{
	    if(!(m_mass.get()))
	    {
                int NumDirBCs = m_locToGloMap->GetNumDirichletBCs();
		int   i,j,n,gid1,gid2,loc_lda,cnt;
		DNekMatSharedPtr loc_mass;
		StdRegions::StdExpansionVectorIter def;

		DNekMatSharedPtr Gmass = MemoryManager<DNekMat>::AllocateSharedPtr(m_contNcoeffs - NumDirBCs,m_contNcoeffs - NumDirBCs);
		
		// fill global matrix 
		for(n = cnt = 0; n < (*m_exp).size(); ++n)
		{
		    loc_mass = (*m_exp)[n]->GetLocMatrix(StdRegions::eMassMatrix);
		    loc_lda = loc_mass->GetColumns();
		    
		    for(i = 0; i < loc_lda; ++i)
		    {
			gid1 = m_locToGloMap->GetMap(cnt + i);
                        if(gid1 >= NumDirBCs)
                        {
                            for(j = 0; j < loc_lda; ++j)
                            {
                                gid2 = m_locToGloMap->GetMap(cnt + j);
                                if(gid2 >= NumDirBCs)
                                {
                                    (*Gmass)(gid1-NumDirBCs,gid2-NumDirBCs) 
                                        += (*loc_mass)(i,j);
                                }
                            }		
                        }
                    }
                    cnt += (*m_exp)[n]->GetNcoeffs();
		}
                m_mass = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmass);
	    }
	}

    } // end of namespace
} //end of namespace
