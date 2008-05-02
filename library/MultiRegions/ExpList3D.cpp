///////////////////////////////////////////////////////////////////////////////
//
// File ExpList3D.cpp
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
// Description: Expansion list 3D definition
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList3D.h>

namespace Nektar
{
  namespace MultiRegions
  {

    ExpList3D::ExpList3D(): ExpList()
    {
    }
        
    ExpList3D::ExpList3D(const ExpList3D &In): ExpList(In)
    {
    }

    ExpList3D::~ExpList3D()
    {
    }

    namespace
        {
            // Adds up the number of cells in a truncated Nc by Nc by Nc pyramid, 
            // where the longest Na rows and longest Nb columns are kept.
            // Example: (Na, Nb, Nc) = (3, 4, 5); The number of coefficients is the 
            // sum of the elements of the following matrix:
            //     |5  4  3  2  0|
            //     |4  3  2  0   |
            //     |3  2  0      |
            //     |0  0         |
            //     |0            |
            // Sum = 28 = number of tet coefficients
            inline int getNumberOfCoefficients( int Na, int Nb, int Nc ) 
            {
                int nCoef = 0;
                for( int a = 0; a < Na; ++a )
                {
                    for( int b = 0; b < Nb - a; ++b )
                    {
                        for( int c = 0; c < Nc - a - b; ++c )
                        {
                            ++nCoef;
                        }
                    }
                }
                return nCoef;
            }
        }

     ExpList3D::ExpList3D(const LibUtilities::BasisKey &Ba,
                          const LibUtilities::BasisKey &Bb,
                          const LibUtilities::BasisKey &Bc,
                          const SpatialDomains::MeshGraph3D &graph3D,
                          const LibUtilities::PointsType TetNb):
            ExpList()
        {

            LocalRegions::TetExpSharedPtr tet;
            LocalRegions::HexExpSharedPtr hex;
            LocalRegions::PrismExpSharedPtr prism;
            LocalRegions::PyrExpSharedPtr pyramid;

            const SpatialDomains::ExpansionVector &expansions = graph3D.GetExpansions();
                        
            for(int i = 0; i < expansions.size(); ++i)
            {
                SpatialDomains::TetGeomSharedPtr TetGeom;
                SpatialDomains::HexGeomSharedPtr HexGeom;
                SpatialDomains::PrismGeomSharedPtr PrismGeom;
                SpatialDomains::PyrGeomSharedPtr PyrGeom;
                
                if(TetGeom = boost::dynamic_pointer_cast<SpatialDomains::TetGeom>(expansions[i]->m_GeomShPtr)) // Tetrahedron
                {
                    if(TetNb < LibUtilities::SIZE_PointsType)
                    {
//                         Ntet = MemoryManager<LocalRegions::NodalTetExp>::AllocateSharedPtr(TetBa,TetBb,TetBc,TetNb,TetGeom);
//                         (*m_exp).push_back(Ntet);
                    }
                    else
                    {
                        tet = MemoryManager<LocalRegions::TetExp>::AllocateSharedPtr(Ba,Bb,Bc,TetGeom);
                        (*m_exp).push_back(tet);
                    }

                       m_ncoeffs += getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes());
                        
                       m_npoints += Ba.GetNumPoints()*Bb.GetNumPoints()*Bc.GetNumPoints();
                }
                else if(PrismGeom = boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(expansions[i]->m_GeomShPtr)) // Prism
                {
                      prism = MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(Ba,Bb,Bc,PrismGeom);
                      (*m_exp).push_back(prism);

                      m_ncoeffs += getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes());
                      m_npoints +=  Ba.GetNumPoints()*Bb.GetNumPoints()*Bc.GetNumPoints();
                
                }
                else if(PyrGeom = boost::dynamic_pointer_cast<SpatialDomains::PyrGeom>(expansions[i]->m_GeomShPtr)) // Pyramid
                {
                     pyramid = MemoryManager<LocalRegions::PyrExp>::AllocateSharedPtr(Ba,Bb,Bc,PyrGeom);
                     (*m_exp).push_back(pyramid);

                      m_ncoeffs += getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes());
                      m_npoints +=  Ba.GetNumPoints()*Bb.GetNumPoints()*Bc.GetNumPoints();

                }
                else if(HexGeom = boost::dynamic_pointer_cast<SpatialDomains::HexGeom>(expansions[i]->m_GeomShPtr)) // Hexahedron
                {
                    hex = MemoryManager<LocalRegions::HexExp>::AllocateSharedPtr(Ba,Bb,Bc, HexGeom);
                    (*m_exp).push_back(hex);
                    
                    m_ncoeffs += Ba.GetNumModes()*Bb.GetNumModes()*Bc.GetNumModes();
                    m_npoints += Ba.GetNumPoints()*Bb.GetNumPoints()*Bc.GetNumPoints();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry3D failed");
                }  
                
            }
            
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
        }


        ExpList3D::ExpList3D(SpatialDomains::MeshGraph3D &graph3D):ExpList()
        {

            LocalRegions::TetExpSharedPtr tet;
            LocalRegions::HexExpSharedPtr hex;
            LocalRegions::PrismExpSharedPtr prism;
            LocalRegions::PyrExpSharedPtr pyramid;
            LibUtilities::PointsType TetNb;

            const SpatialDomains::ExpansionVector &expansions = graph3D.GetExpansions();
            
            for(int i = 0; i < expansions.size(); ++i)
            {
                SpatialDomains::TetGeomSharedPtr TetGeom;
                SpatialDomains::HexGeomSharedPtr HexGeom;
                SpatialDomains::PrismGeomSharedPtr PrismGeom;
                SpatialDomains::PyrGeomSharedPtr PyrGeom;
                
                if(TetGeom = boost::dynamic_pointer_cast<SpatialDomains::TetGeom>(expansions[i]->m_GeomShPtr))
                {
                    LibUtilities::BasisKey TetBa = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey TetBb = graph3D.GetBasisKey(expansions[i],1);
                    LibUtilities::BasisKey TetBc = graph3D.GetBasisKey(expansions[i],2);
                    
                    if(expansions[i]->m_ExpansionType == SpatialDomains::eNodal)
                    {
//                         TriNb = LibUtilities::eNodalTriElec;
//                         Ntri = MemoryManager<LocalRegions::NodalTriExp>::AllocateSharedPtr(TriBa,TriBb,TriNb,TriangleGeom);
//                         (*m_exp).push_back(Ntri);
                    }
                    else
                    {
                        tet = MemoryManager<LocalRegions::TetExp>::AllocateSharedPtr(TetBa,TetBb,TetBc,TetGeom);
                        (*m_exp).push_back(tet);
                    }
                        
                    m_ncoeffs += getNumberOfCoefficients(TetBa.GetNumModes(), TetBb.GetNumModes(), TetBc.GetNumModes());
                    m_npoints += TetBa.GetNumPoints()*TetBb.GetNumPoints()*TetBc.GetNumPoints();
                }
                else if(PrismGeom = boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(expansions[i]->m_GeomShPtr))
                {
                    LibUtilities::BasisKey PrismBa = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey PrismBb = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey PrismBc = graph3D.GetBasisKey(expansions[i],1);

                    prism = MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(PrismBa,PrismBb,PrismBc,PrismGeom);
                    (*m_exp).push_back(prism);

                    m_ncoeffs += getNumberOfCoefficients(PrismBa.GetNumModes(), PrismBb.GetNumModes(), PrismBc.GetNumModes());
                    m_npoints += PrismBa.GetNumPoints()*PrismBb.GetNumPoints()*PrismBc.GetNumPoints();
                }
                else if(PyrGeom = boost::dynamic_pointer_cast<SpatialDomains::PyrGeom>(expansions[i]->m_GeomShPtr))
                {
                    LibUtilities::BasisKey PyrBa = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey PyrBb = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey PyrBc = graph3D.GetBasisKey(expansions[i],2);

                    pyramid = MemoryManager<LocalRegions::PyrExp>::AllocateSharedPtr(PyrBa,PyrBb,PyrBc,PyrGeom);
                    (*m_exp).push_back(pyramid);

                    m_ncoeffs += getNumberOfCoefficients(PyrBa.GetNumModes(), PyrBb.GetNumModes(), PyrBc.GetNumModes());
                    m_npoints += PyrBa.GetNumPoints()*PyrBb.GetNumPoints()*PyrBc.GetNumPoints();
                }
                else if(HexGeom = boost::dynamic_pointer_cast<SpatialDomains::HexGeom>(expansions[i]->m_GeomShPtr))
                {
                    LibUtilities::BasisKey HexBa = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey HexBb = graph3D.GetBasisKey(expansions[i],0);
                    LibUtilities::BasisKey HexBc = graph3D.GetBasisKey(expansions[i],0);
                    
                    hex = MemoryManager<LocalRegions::HexExp>::AllocateSharedPtr(HexBa,HexBb,HexBc,HexGeom);
                    (*m_exp).push_back(hex);
                    
                    m_ncoeffs += HexBa.GetNumModes()*HexBb.GetNumModes()*HexBc.GetNumModes();
                    m_npoints += HexBa.GetNumPoints()*HexBb.GetNumPoints()*HexBc.GetNumPoints();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry3D failed");
                }  
                
            }            
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
        }


  } //end of namespace
} //end of namespace

