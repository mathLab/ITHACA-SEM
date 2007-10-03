///////////////////////////////////////////////////////////////////////////////
//
// File ContExpList2D.cpp
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
// Description: Continuous Expansion list definition in 2D
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContExpList2D.h>

namespace Nektar
{
    namespace MultiRegions
    {
    
    ContExpList2D::ContExpList2D()
    {
    }
    
    ContExpList2D::~ContExpList2D()
    {
    }

        ContExpList2D::ContExpList2D(const ContExpList2D &In):
            ExpList2D(In),
            m_contNcoeffs(In.m_contNcoeffs),
            m_locToGloMap(In.m_locToGloMap),
            m_mass(In.m_mass)
            
        {
            m_contCoeffs = Array<OneD,NekDouble>(m_contNcoeffs);
        }

        ContExpList2D::ContExpList2D(const LibUtilities::BasisKey &TriBa, 
                                     const LibUtilities::BasisKey &TriBb,
                                     const LibUtilities::BasisKey &QuadBa, 
                                     const LibUtilities::BasisKey &QuadBb, 
                                     const SpatialDomains::MeshGraph2D &graph2D,
                                     const LibUtilities::PointsType TriNb):
        ExpList2D(TriBa,TriBb,QuadBa,QuadBb,graph2D,TriNb)
    {
        
        ASSERTL1((TriBa.GetBasisType() == LibUtilities::eModified_A)&&
             (TriNb == LibUtilities::SIZE_PointsType),
             "Expansion not of an boundary-interior type");
        
        ASSERTL1((TriBb.GetBasisType() == LibUtilities::eModified_B)&&
             (TriNb == LibUtilities::SIZE_PointsType),
             "Expansion not of an boundary-interior type");

        ASSERTL1((QuadBa.GetBasisType() == LibUtilities::eModified_A)
             ||(QuadBa.GetBasisType() == LibUtilities::eGLL_Lagrange),
             "Expansion not of an boundary-interior type");

        ASSERTL1((QuadBb.GetBasisType() == LibUtilities::eModified_A)
             ||(QuadBb.GetBasisType() == LibUtilities::eGLL_Lagrange),
             "Expansion not of an boundary-interior type");
        
        if(TriBa.GetBasisType() == LibUtilities::eModified_A)
        {
        ASSERTL1((QuadBb.GetBasisType() == LibUtilities::eModified_A),
             "Quad and Tri Expansions are not of the same type");
        }
               
        // setup mapping array 
        m_locToGloMap = MemoryManager<LocalToGlobalMap2D>::AllocateSharedPtr(m_ncoeffs,*m_exp,graph2D);
 
        m_contNcoeffs = m_locToGloMap->GetTotGloLen();
        m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs);
    }

        void ContExpList2D::IProductWRTBase(const ExpList &In)
    {
        ExpList2D::IProductWRTBase(In);
        Assemble();
        m_transState = eLocalCont;
    }

    void ContExpList2D::FwdTrans(const ExpList &In)
    {
        IProductWRTBase(In);

        if(!(m_mass.get()))
        {
            GenMassMatrixLinSys();
        }        
        
        DNekVec v(m_contNcoeffs,m_contCoeffs,eWrapper);
        m_mass->Solve(v,v);
        m_transState = eContinuous;
        m_physState = false;
    }

    void ContExpList2D::BwdTrans(const ExpList &In)
    {
        
        if(m_transState == eContinuous)
        {
        ContToLocal();
        }
        
        ExpList2D::BwdTrans(In);
    }
    

    void ContExpList2D::GenMassMatrixLinSys(void)
    {
        if(!(m_mass.get()))
        {
        int   i,j,n,cnt,gid1,gid2,loc_lda;
        NekDouble sign1,sign2;
        DNekScalMatSharedPtr loc_mass;
        StdRegions::StdExpansionVectorIter def;

        DNekMatSharedPtr Gmass = MemoryManager<DNekMat>::AllocateSharedPtr(m_contNcoeffs,m_contNcoeffs);
    
        // fill global matrix 
        for(n = cnt = 0; n < (*m_exp).size(); ++n)
        {
                    loc_mass = (*m_exp)[n]->GetLocMatrix(StdRegions::eMass);
                    loc_lda = loc_mass->GetColumns();
            
                    for(i = 0; i < loc_lda; ++i)
                    {
                        gid1  =  m_locToGloMap->GetMap(cnt + i);
                        sign1 =  m_locToGloMap->GetSign(cnt + i);
            
                        for(j = 0; j < loc_lda; ++j)
                        {
                            gid2  = m_locToGloMap->GetMap(cnt + j);    
                            sign2 = m_locToGloMap->GetSign(cnt + j);
                            
                            (*Gmass)(gid1,gid2) +=  sign1*sign2*(*loc_mass)(i,j);
                        }
                    }
                    cnt += (*m_exp)[n]->GetNcoeffs();            
        }
                m_mass = MemoryManager<DNekLinSys>::AllocateSharedPtr(Gmass);
        }
    }
    } //end of namespace
} //end of namespace

/**
* $Log: ContExpList2D.cpp,v $
* Revision 1.7  2007/07/20 02:04:10  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.6  2007/07/13 16:48:46  pvos
* Another HelmHoltz update (homogeneous dir BC multi-elemental solver does work)
*
* Revision 1.5  2007/06/07 15:54:19  pvos
* Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
* Also made corrections to various ASSERTL2 calls
*
* Revision 1.4  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/
