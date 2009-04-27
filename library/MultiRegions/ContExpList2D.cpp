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
        
        ContExpList2D::ContExpList2D():
            ExpList2D(),
            m_locToGloMap(),
            m_contNcoeffs(0),
            m_contCoeffs(),
            m_globalMat(),
            m_globalLinSys()
        {
        }
        
        ContExpList2D::~ContExpList2D()
        {
        }
        
        ContExpList2D::ContExpList2D(const ContExpList2D &In):
            ExpList2D(In),
            m_locToGloMap(In.m_locToGloMap),
            m_contNcoeffs(In.m_contNcoeffs),
            m_contCoeffs(m_contNcoeffs,0.0),
            m_globalMat(In.m_globalMat),
            m_globalLinSys(In.m_globalLinSys)            
        {
        }

        ContExpList2D::ContExpList2D(const LibUtilities::BasisKey &TriBa, 
                                     const LibUtilities::BasisKey &TriBb,
                                     const LibUtilities::BasisKey &QuadBa, 
                                     const LibUtilities::BasisKey &QuadBb, 
                                     const SpatialDomains::MeshGraph2D &graph2D,
                                     const LibUtilities::PointsType TriNb,
                                     const bool constructMap):
            ExpList2D(TriBa,TriBb,QuadBa,QuadBb,graph2D,TriNb),
            m_locToGloMap(),
            m_contNcoeffs(),
            m_contCoeffs(),
            m_globalMat(MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
            m_globalLinSys(MemoryManager<GlobalLinSysMap>::AllocateSharedPtr())
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
            if(constructMap)
            {
                m_locToGloMap = MemoryManager<LocalToGlobalC0ContMap>::AllocateSharedPtr(m_ncoeffs,*m_exp);
                m_contNcoeffs = m_locToGloMap->GetNumGlobalCoeffs();
                m_contCoeffs  = Array<OneD,NekDouble>(m_contNcoeffs,0.0);
            }
        }

        ContExpList2D::ContExpList2D(SpatialDomains::MeshGraph2D &graph2D,
                                     const bool constructMap):
            ExpList2D(graph2D),
            m_locToGloMap(),
            m_contNcoeffs(),
            m_contCoeffs(),
            m_globalMat(MemoryManager<GlobalMatrixMap>::AllocateSharedPtr()),
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

        void ContExpList2D::FwdTrans(const Array<OneD, const NekDouble> &inarray, 
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
                mass_matrix->Solve(outarray,outarray,*m_locToGloMap);
            }
            else
            {
                Array<OneD, NekDouble> wsp(m_contNcoeffs);
                mass_matrix->Solve(outarray,wsp,*m_locToGloMap);
                GlobalToLocal(wsp,outarray);
            }

        }

        GlobalMatrixSharedPtr ContExpList2D::GetGlobalMatrix(const GlobalMatrixKey &mkey) 
        {
            ASSERTL1(mkey.LocToGloMapIsDefined(),
                     "To use method must have a LocalToGlobalBaseMap "
                     "attached to key");

            GlobalMatrixSharedPtr glo_matrix;
            GlobalMatrixMap::iterator matrixIter = m_globalMat->find(mkey);

            if(matrixIter == m_globalMat->end())
            {
                glo_matrix = GenGlobalMatrix(mkey,m_locToGloMap);
                (*m_globalMat)[mkey] = glo_matrix;
            }
            else
            {
                glo_matrix = matrixIter->second;
            }

            return glo_matrix;
        }
   
    } //end of namespace
} //end of namespace

/**
* $Log: ContExpList2D.cpp,v $
* Revision 1.19  2009/03/04 14:17:38  pvos
* Removed all methods that take and Expansion as argument
*
* Revision 1.18  2009/02/08 09:06:19  sherwin
* Added updated definition of GlobalLinSysKey to include localtoglobalbasemap
*
* Revision 1.17  2009/01/06 21:05:56  sherwin
* Added virtual function calls for BwdTrans, FwdTrans and IProductWRTBase from the class ExpList. Introduced _IterPerExp versions of these methods in ExpList.cpp§
*
* Revision 1.16  2008/09/16 13:36:05  pvos
* Restructured the LocalToGlobalMap classes
*
* Revision 1.15  2008/08/26 02:42:51  ehan
* Many modifications in order to be consistent with changes in 2D and 1D.
*
* Revision 1.14  2008/05/07 16:05:55  pvos
* Mapping + Manager updates
*
* Revision 1.13  2008/04/06 06:00:06  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.12  2008/03/18 14:14:13  pvos
* Update for nodal triangular helmholtz solver
*
* Revision 1.11  2008/03/12 15:25:45  pvos
* Clean up of the code
*
* Revision 1.10  2007/12/17 13:05:04  sherwin
* Made files compatible with modifications in StdMatrixKey which now holds constants
*
* Revision 1.9  2007/12/06 22:52:29  pvos
* 2D Helmholtz solver updates
*
* Revision 1.8  2007/10/03 11:37:50  sherwin
* Updates relating to static condensation implementation
*
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
