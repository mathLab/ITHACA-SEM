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
            m_globalMat()
        {
        }
                
        ContExpList3D::ContExpList3D(const ContExpList3D &In):
            ExpList3D(In),
            m_locToGloMap(In.m_locToGloMap),
            m_contNcoeffs(In.m_contNcoeffs),
            m_contCoeffs(m_contNcoeffs,0.0),
            m_globalMat(In.m_globalMat)            
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
            m_globalMat(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
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
            m_globalMat(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
        {                        
            // setup mapping array 
            if(constructMap)
            {
                m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp);
                m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
                m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
            }
        }

        void ContExpList3D::IProductWRTBase(const ExpList &In)
        {
            IProductWRTBase(In.GetPhys(),m_contCoeffs, m_coeffs);
            m_transState = eLocalCont;
        }


        // Note inarray is assumed to be of size m_contExpList2D;
        void ContExpList3D::IProductWRTBase(const Array<OneD, const NekDouble> &inarray,  Array<OneD, NekDouble> &outarray, Array<OneD, NekDouble> &wksp)
        {
            Array<OneD, NekDouble> tmp;
            
            if(wksp == NullNekDouble1DArray)
            {
                tmp =  Array<OneD, NekDouble>(m_ncoeffs);
            }
            else
            {
                tmp = wksp;
            }

            IProductWRTBase_IterPerExp(inarray,tmp);
            
            Assemble(tmp,outarray);
        }

        void ContExpList3D::GeneralMatrixOp(const GlobalLinSysKey &gkey,
                                            const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD, NekDouble> &outarray)
            
        {
                Array<OneD,NekDouble> tmp(m_ncoeffs);
                GlobalToLocal(inarray,tmp);
                ExpList3D::GeneralMatrixOp(gkey,tmp,tmp);
                Assemble(tmp,outarray);
        }

        void ContExpList3D::FwdTrans(const ExpList &In)
        {
            FwdTrans(In.GetPhys(),m_contCoeffs);
            
	    m_transState = eContinuous;
	    m_physState = false;
        }

        void ContExpList3D::FwdTrans(const Array<OneD, const NekDouble> &inarray, 
                             Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(inarray,outarray);

            GlobalLinSysSharedPtr mass_matrix;
            GlobalLinSysKey key(StdRegions::eMass);
            GlobalLinSysMap::iterator matrixIter = m_globalMat->find(key);
           
            if(matrixIter == m_globalMat->end())
            {
                mass_matrix = GenGlobalLinSys(key,m_locToGloMap);
                (*m_globalMat)[key] = mass_matrix;
            }
            else
            {
                mass_matrix = matrixIter->second;
            }
            mass_matrix->Solve(outarray,outarray,*m_locToGloMap);

        }

        void ContExpList3D::BwdTrans(const ExpList &In)
        {
            Array<OneD, NekDouble> tmp;

            if(m_transState == eContinuous)
            {
                GlobalToLocal(In.GetContCoeffs(),m_coeffs);
            }         
            else
            {
                tmp = In.GetCoeffs();
            }

            BwdTrans_IterPerExp(In.GetCoeffs(),m_phys);

	    m_physState = true;
        }

        void ContExpList3D::BwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, NekDouble> tmp(m_ncoeffs);
            
            GlobalToLocal(inarray,tmp);

            BwdTrans_IterPerExp(tmp,outarray);
        }
        
           
    } //end of namespace
} //end of namespace


