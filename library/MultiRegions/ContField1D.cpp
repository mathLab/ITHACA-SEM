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

        ContField1D::ContField1D(const ContField1D &In):
            ContExpList1D(In)
        {
        }

        ContField1D::ContField1D(const LibUtilities::BasisKey &Ba, 
                                 const SpatialDomains::Composite &cmps,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const int bc_loc):
            ContExpList1D(Ba,cmps)  
        {
            
            GenerateField1D(bcs,bcs.GetVariable(bc_loc));
        }

        ContField1D::ContField1D(const LibUtilities::BasisKey &Ba, 
                                 const SpatialDomains::Composite &cmps,
                                 SpatialDomains::BoundaryConditions &bcs, 
                                 const std::string variable):
            ContExpList1D(Ba,cmps)
        {
            GenerateField1D(bcs,variable);
        }

        
        void ContField1D::GenerateField1D(SpatialDomains::BoundaryConditions &bcs, 
                                     const std::string variable)
        {
	    int i,nbnd;
            LocalRegions::PointExpSharedPtr  p_exp;
            int NumDirichlet = 0;

            SpatialDomains::BoundaryRegionCollection    &bregions = bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryConditionCollection &bconditions = bcs.GetBoundaryConditions();

            nbnd = bregions.size();

            // list Dirichlet boundaries first
            for(i = 0; i < nbnd; ++i)
            {
                SpatialDomains::BoundaryConditionShPtr LocBCond = (*(bconditions[i]))[variable];
                
                if(LocBCond->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    SpatialDomains::VertexComponentSharedPtr vert;
                    
                    if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[0])[0]))
                    {
                        Array<OneD,NekDouble> coords(3,0.0);

                        p_exp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                        vert->GetCoords(coords);
                        p_exp->SetValue(boost::static_pointer_cast<SpatialDomains::DirichletBoundaryCondition>(LocBCond)->m_DirichletCondition.Evaluate(coords[0],coords[1],coords[2]));


                        m_bndConstraint.push_back(p_exp);
                        m_bndTypes.push_back(SpatialDomains::eDirichlet);
                        NumDirichlet++;
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a vertex failed");
                    }
                }
            }

            for(i = 0; i < nbnd; ++i)
            {
                SpatialDomains::BoundaryConditionShPtr LocBCond = (*(bconditions[i]))[variable];
                SpatialDomains::BoundaryConditionType BCtype  = LocBCond->GetBoundaryConditionType();

                if(BCtype != SpatialDomains::eDirichlet)
                {
                    SpatialDomains:: VertexComponentSharedPtr vert;
                    
                    if(vert = boost::dynamic_pointer_cast<SpatialDomains::VertexComponent>((*(*bregions[i])[0])[0]))
                    {
                        Array<OneD,NekDouble> coords(3,0.0);
                        p_exp = MemoryManager<LocalRegions::PointExp>::AllocateSharedPtr(vert);
                        vert->GetCoords(coords);

                        if(BCtype == SpatialDomains::eNeumann)
                        {
                            p_exp->SetValue(boost::static_pointer_cast<SpatialDomains::NeumannBoundaryCondition>(LocBCond)->m_NeumannCondition.Evaluate(coords[0],coords[1],coords[2]));
                        }
                        else // Robin boundary condition
                        {
                            boost::shared_ptr<SpatialDomains::RobinBoundaryCondition> robinBC  = boost::static_pointer_cast<SpatialDomains::RobinBoundaryCondition>(LocBCond);
                            p_exp->SetValue(robinBC->m_a.Evaluate(coords[0],coords[1],coords[2])/robinBC->m_b.Evaluate(coords[0],coords[1],coords[2]));

                        }

                        m_bndConstraint.push_back(p_exp);
                        m_bndTypes.push_back(LocBCond->GetBoundaryConditionType());
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a vertex failed");
                    }
                }
            }
            

            // Need to reset numbering according to Dirichlet Bondary Condition
            m_locToGloMap->ResetMapping(NumDirichlet,bcs,variable);
            // Need to know how to get global vertex id 
            // I note that boundary Regions is a composite vector and so belive we can use this to get the desired quantities. 

        }

        ContField1D::~ContField1D()
	{
	}

	void ContField1D::FwdTrans(const ExpList &In)
	{
            GlobalLinSysKey key(StdRegions::eMass);
            GlobalSolve(key,In);

	    m_transState = eContinuous;
	    m_physState = false;
	}

        // Solve the helmholtz problem assuming that m_contCoeff vector 
        // contains an intial estimate for solution
	void ContField1D::HelmSolve(const ExpList &In, NekDouble lambda)
        {
            GlobalLinSysKey key(StdRegions::eHelmholtz,lambda);
            GlobalSolve(key,In);
        }


	void ContField1D::GlobalSolve(const GlobalLinSysKey &key, const ExpList &Rhs)
	{
            int i;
            int NumDirBcs = m_locToGloMap->GetNumDirichletBCs();
            Array<OneD,NekDouble> sln;
            Array<OneD,NekDouble> init(m_contNcoeffs,0.0);
            Array<OneD,NekDouble> Dir_fce(m_contNcoeffs,0.0);

            //assume m_contCoeffs contains initial estimate
            // Set BCs in m_contCoeffs
            Blas::Dcopy(m_contNcoeffs,&m_contCoeffs[0],1,&init[0],1);
            for(i = 0; i < NumDirBcs; ++i)
            {
                init[i] = m_bndConstraint[i]->GetValue();
            }

            GeneralMatrixOp(key.GetLinSysType(), init, Dir_fce,
                            key.GetScaleFactor());

            // Set up forcing function
	    IProductWRTBase(Rhs);
            Vmath::Neg(m_contNcoeffs,&m_contCoeffs[0],1);

            // Forcing function with Dirichlet conditions 
            Vmath::Vsub(m_contNcoeffs,&m_contCoeffs[0],1,
                        &Dir_fce[0],1,&m_contCoeffs[0],1);

            if(m_contNcoeffs - NumDirBcs > 0)
            {
                GlobalLinSysSharedPtr LinSys = GetGlobalLinSys(key);

                sln = m_contCoeffs+NumDirBcs;
                DNekVec v(m_contNcoeffs-NumDirBcs,sln,eWrapper);
                LinSys->GetLinSys()->Solve(v,v);
            }

            // Recover solution by addinig intial conditons
            Vmath::Zero(NumDirBcs,&m_contCoeffs[0],1);
            Vmath::Vadd(m_contNcoeffs,&init[0],1,&m_contCoeffs[0],1,
                        &m_contCoeffs[0],1);

	    m_transState = eContinuous;
	    m_physState = false;
	}

        GlobalLinSysSharedPtr ContField1D::GetGlobalLinSys(const GlobalLinSysKey &mkey) 
       {
            GlobalLinSysSharedPtr glo_matrix;
            GlobalLinSysMap::iterator matrixIter = m_globalMat->find(mkey);
           
            if(matrixIter == m_globalMat->end())
            {
                glo_matrix = GenGlobalLinSys(mkey,
                                             m_locToGloMap->GetNumDirichletBCs());
                (*m_globalMat)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }

    } // end of namespace
} //end of namespace
