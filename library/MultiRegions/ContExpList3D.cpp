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
// Description: Continuous Expansion list definition in 3D
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ContExpList3D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        
        ContExpList3D::ContExpList3D():
            ExpList3D(),
            m_locToGloMap(),
            m_contNcoeffs(0),
            m_contCoeffs(),
            m_globalLinSys()
        {
        }
                
        ContExpList3D::ContExpList3D(const ContExpList3D &In):
            ExpList3D(In),
            m_locToGloMap(In.m_locToGloMap),
            m_contNcoeffs(In.m_contNcoeffs),
            m_contCoeffs(m_contNcoeffs,0.0),
            m_globalLinSys(In.m_globalLinSys)            
        {
        }

        
        ContExpList3D::~ContExpList3D()
        {
        }

        ContExpList3D::ContExpList3D(const LibUtilities::BasisKey &Ba,
                                     const LibUtilities::BasisKey &Bb,
                                     const LibUtilities::BasisKey &Bc,
                                     const SpatialDomains::MeshGraph3D &graph3D,
                                     const LibUtilities::PointsType TetNb,
                                     const bool constructMap):
            ExpList3D(Ba,Bb,Bc,graph3D,TetNb),
            m_locToGloMap(),
            m_contNcoeffs(),
            m_contCoeffs(),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {
            int Tet, Hex, Prism, Pyramid; //TODO: check this.
            if(Tet ){ // Tetrahedron
            ASSERTL1((Ba.GetBasisType() == LibUtilities::eModified_A)&&
                     (TetNb == LibUtilities::SIZE_PointsType),
                     "Expansion not of an boundary-interior type");
            
            ASSERTL1((Bb.GetBasisType() == LibUtilities::eModified_B)&&
                     (TetNb == LibUtilities::SIZE_PointsType),
                     "Expansion not of an boundary-interior type");

            ASSERTL1((Bc.GetBasisType() == LibUtilities::eModified_C)&&
                     (TetNb == LibUtilities::SIZE_PointsType),
                     "Expansion not of an boundary-interior type");
            
            } else if (Hex) // Hexahedron
            {
            ASSERTL1((Ba.GetBasisType() == LibUtilities::eModified_A)
                     ||(Ba.GetBasisType() == LibUtilities::eGLL_Lagrange),
                     "Expansion not of an boundary-interior type");
            
            ASSERTL1((Bb.GetBasisType() == LibUtilities::eModified_A)
                     ||(Bb.GetBasisType() == LibUtilities::eGLL_Lagrange),
                     "Expansion not of an boundary-interior type");
                     
             ASSERTL1((Bc.GetBasisType() == LibUtilities::eModified_A)
                     ||(Bc.GetBasisType() == LibUtilities::eGLL_Lagrange),
                     "Expansion not of an boundary-interior type");
            }
            else if (Prism) // Prism
            {
            ASSERTL1((Ba.GetBasisType() == LibUtilities::eModified_A)
                     ||(Ba.GetBasisType() == LibUtilities::eGLL_Lagrange),
                     "Expansion not of an boundary-interior type");
            
            ASSERTL1((Bb.GetBasisType() == LibUtilities::eModified_A)
                     ||(Bb.GetBasisType() == LibUtilities::eGLL_Lagrange),
                     "Expansion not of an boundary-interior type");
                     
             ASSERTL1((Bc.GetBasisType() == LibUtilities::eModified_B)
                     &&(TetNb == LibUtilities::SIZE_PointsType),
                     "Expansion not of an boundary-interior type");
            }
             else if (Pyramid) // Pyramid
            {
            ASSERTL1((Ba.GetBasisType() == LibUtilities::eModified_A)
                     ||(Ba.GetBasisType() == LibUtilities::eGLL_Lagrange),
                     "Expansion not of an boundary-interior type");
            
            ASSERTL1((Bb.GetBasisType() == LibUtilities::eModified_A)
                     ||(Bb.GetBasisType() == LibUtilities::eGLL_Lagrange),
                     "Expansion not of an boundary-interior type");
                     
             ASSERTL1((Bc.GetBasisType() == LibUtilities::eModified_C)
                     &&(TetNb == LibUtilities::SIZE_PointsType),
                     "Expansion not of an boundary-interior type");
            }
            
                        
            // setup mapping array 
            if(constructMap)
            { 
                m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp);
                m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
                m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
            }
        }

        ContExpList3D::ContExpList3D(SpatialDomains::MeshGraph3D &graph3D, const bool constructMap):
            ExpList3D(graph3D),
            m_locToGloMap(),
            m_contNcoeffs(),
            m_contCoeffs(),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {                        
            // setup mapping array 
            if(constructMap)
            {
                m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp);
                m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
                m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
            }
        }

        void ContExpList3D::IProductWRTBase(const Array<OneD, const NekDouble> &inarray,  
                                                  Array<OneD,       NekDouble> &outarray,
                                            bool  UseContCoeffs)
        {
            if(UseContCoeffs)
            {
                Array<OneD, NekDouble> wsp(m_ncoeffs);
                IProductWRTBase_IterPerExp(inarray,wsp);
                Assemble(wsp,outarray);
            }
            else
            {
                IProductWRTBase_IterPerExp(inarray,outarray);
            }
        }

        void ContExpList3D::GeneralMatrixOp(const GlobalMatrixKey             &gkey,
                                            const Array<OneD,const NekDouble> &inarray, 
                                                  Array<OneD,      NekDouble> &outarray,
                                            bool  UseContCoeffs)
            
        {
            if(UseContCoeffs)
            {
                Array<OneD,NekDouble> tmp1(2*m_ncoeffs);
                Array<OneD,NekDouble> tmp2(tmp1+m_ncoeffs);
                GlobalToLocal(inarray,tmp1);
                GeneralMatrixOp_IterPerExp(gkey,tmp1,tmp2);
                Assemble(tmp2,outarray);
            }
            else
            {
                GeneralMatrixOp_IterPerExp(gkey,inarray,outarray);
            }
        }

        void ContExpList3D::FwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                           Array<OneD,       NekDouble> &outarray,
                                     bool  UseContCoeffs)
        {
            GlobalLinSysSharedPtr mass_matrix;
            GlobalLinSysKey key(StdRegions::eMass, m_locToGloMap);
            GlobalLinSysMap::iterator matrixIter = m_globalLinSys->find(key);
            
            if(matrixIter == m_globalLinSys->end())
            {
                mass_matrix = GenGlobalLinSys(key,m_locToGloMap);
                (*m_globalLinSys)[key] = mass_matrix;
            }
            else
            {
                mass_matrix = matrixIter->second;
            }
            
            IProductWRTBase(inarray,outarray,true);

            if(UseContCoeffs)
            {                
                mass_matrix->Solve(outarray,outarray,*m_locToGloMap,this);
            }
            else
            {
                Array<OneD, NekDouble> wsp(m_contNcoeffs);
                mass_matrix->Solve(outarray,wsp,*m_locToGloMap,this);
                GlobalToLocal(wsp,outarray);
            }

        }
        
        void ContExpList3D::BwdTrans(const Array<OneD, const NekDouble> &inarray, 
                                           Array<OneD,       NekDouble> &outarray,
                                     bool  UseContCoeffs)
        {
            if(UseContCoeffs)
            {
                Array<OneD, NekDouble> wsp(m_ncoeffs);
                GlobalToLocal(inarray,wsp);
                BwdTrans_IterPerExp(wsp,outarray);
            }
            else
            {
                BwdTrans_IterPerExp(inarray,outarray);
            }
        }
        
           
    } //end of namespace
} //end of namespace


