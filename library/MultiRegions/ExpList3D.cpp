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

#include <iomanip>
#include <MultiRegions/ExpList3D.h>

#include <LocalRegions/HexExp.h>
#include <LocalRegions/PrismExp.h>
#include <LocalRegions/PyrExp.h>
#include <LocalRegions/TetExp.h>

#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/PhysGalerkinProject.h>

using namespace std;

namespace Nektar
{
    namespace MultiRegions
    {

        ExpList3D::ExpList3D(): ExpList(e3D)
        {
        }

        ExpList3D::ExpList3D(const ExpList3D &In): ExpList(In)
        {
        }
        
        ExpList3D::ExpList3D(const ExpList3D &In,
                const std::vector<unsigned int> &eIDs,
                const bool DeclareCoeffPhysArrays,
                const Collections::ImplementationType ImpType):
            ExpList(In, eIDs, DeclareCoeffPhysArrays,ImpType)
        {
        }

        ExpList3D::~ExpList3D()
        {
        }


        ExpList3D::ExpList3D(const LibUtilities::SessionReaderSharedPtr &pSession,
                             const LibUtilities::BasisKey &TBa,
                             const LibUtilities::BasisKey &TBb,
                             const LibUtilities::BasisKey &TBc,
                             const LibUtilities::BasisKey &HBa,
                             const LibUtilities::BasisKey &HBb,
                             const LibUtilities::BasisKey &HBc,
                             const SpatialDomains::MeshGraphSharedPtr &graph3D,
                             const LibUtilities::PointsType TetNb,
                             const Collections::ImplementationType ImpType):
            ExpList(e3D, pSession,graph3D)
        {

            LocalRegions::TetExpSharedPtr tet;
            LocalRegions::HexExpSharedPtr hex;
            LocalRegions::PrismExpSharedPtr prism;
            LocalRegions::PyrExpSharedPtr pyramid;

            const SpatialDomains::ExpansionMap &expansions = graph3D->GetExpansions();

            for (auto &expIt : expansions)
            {
                SpatialDomains::TetGeomSharedPtr TetGeom;
                SpatialDomains::HexGeomSharedPtr HexGeom;
                SpatialDomains::PrismGeomSharedPtr PrismGeom;
                SpatialDomains::PyrGeomSharedPtr PyrGeom;

                if((TetGeom = std::dynamic_pointer_cast<SpatialDomains::TetGeom>(expIt.second->m_geomShPtr)))
                {
                    if(TetNb < LibUtilities::SIZE_PointsType)
                    {
//                         Ntet = MemoryManager<LocalRegions::NodalTetExp>::AllocateSharedPtr(TetBa,TetBb,TetBc,TetNb,TetGeom);
//                         (*m_exp).push_back(Ntet);
                    }
                    else
                    {
                        tet = MemoryManager<LocalRegions::TetExp>::AllocateSharedPtr(TBa,TBb,TBc,TetGeom);
                        (*m_exp).push_back(tet);
                    }

                    m_ncoeffs += LibUtilities::StdTetData::getNumberOfCoefficients(TBa.GetNumModes(), TBb.GetNumModes(), TBc.GetNumModes());
                    
                    m_npoints += TBa.GetNumPoints()*TBb.GetNumPoints()*TBc.GetNumPoints();
                }
/*
                else if((PrismGeom = std::dynamic_pointer_cast<SpatialDomains::PrismGeom>(expansions[i]->m_geomShPtr)))
                {
                      prism = MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(Ba,Bb,Bc,PrismGeom);
                      (*m_exp).push_back(prism);

                      m_ncoeffs += StdRegions::StdPrismData::getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes());
                      m_npoints +=  Ba.GetNumPoints()*Bb.GetNumPoints()*Bc.GetNumPoints();

                }
                else if((PyrGeom = std::dynamic_pointer_cast<SpatialDomains::PyrGeom>(expansions[i]->m_geomShPtr)))
                {
                     pyramid = MemoryManager<LocalRegions::PyrExp>::AllocateSharedPtr(Ba,Bb,Bc,PyrGeom);
                     (*m_exp).push_back(pyramid);

                     m_ncoeffs += StdRegions::StdPyrData::getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes());
                      m_npoints +=  Ba.GetNumPoints()*Bb.GetNumPoints()*Bc.GetNumPoints();

                }
*/
                else if((HexGeom = std::dynamic_pointer_cast<SpatialDomains::HexGeom>(expIt.second->m_geomShPtr)))
                {
                    hex = MemoryManager<LocalRegions::HexExp>::AllocateSharedPtr(HBa,HBb,HBc, HexGeom);
                    (*m_exp).push_back(hex);

                    m_ncoeffs += HBa.GetNumModes()*HBb.GetNumModes()*HBc.GetNumModes();
                    m_npoints += HBa.GetNumPoints()*HBb.GetNumPoints()*HBc.GetNumPoints();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry3D failed");
                }

            }

            SetCoeffPhysOffsets();

            // Set up m_coeffs, m_phys.
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs,0.0);
            m_phys   = Array<OneD, NekDouble>(m_npoints,0.0);

            CreateCollections(ImpType);
        }

        /**
         * Given a mesh \a graph3D, containing information about the domain and
         * the spectral/hp element expansion, this constructor fills the list
         * of local expansions \texttt{m_exp} with the proper expansions,
         * calculates the total number of quadrature points
         * \f$\boldsymbol{x}_i\f$ and the local expansion coefficients
         * \f$\hat{u}^e_n\f$ and allocates memory for the arrays #m_coeffs and
         * #m_phys.
         *
         * @param   graph3D     A mesh, containing information about the domain
         *                      and the spectral/hp element expansion.
         */
        ExpList3D::ExpList3D(const LibUtilities::SessionReaderSharedPtr &pSession,
                             const SpatialDomains::MeshGraphSharedPtr &graph3D,
                             const std::string  &variable,
                             const Collections::ImplementationType ImpType):
            ExpList(e3D, pSession,graph3D)
        {
            int elmtid = 0;
            LocalRegions::TetExpSharedPtr   tet;
            LocalRegions::HexExpSharedPtr   hex;
            LocalRegions::PrismExpSharedPtr prism;
            LocalRegions::PyrExpSharedPtr   pyramid;

            const SpatialDomains::ExpansionMap &expansions
                                        = graph3D->GetExpansions(variable);

            for (auto &expIt : expansions)
            {
                SpatialDomains::TetGeomSharedPtr   TetGeom;
                SpatialDomains::HexGeomSharedPtr   HexGeom;
                SpatialDomains::PrismGeomSharedPtr PrismGeom;
                SpatialDomains::PyrGeomSharedPtr   PyrGeom;

                if((TetGeom = std::dynamic_pointer_cast<
                        SpatialDomains::TetGeom>(expIt.second->m_geomShPtr)))
                {
                    LibUtilities::BasisKey TetBa
                                        = expIt.second->m_basisKeyVector[0];
                    LibUtilities::BasisKey TetBb
                                        = expIt.second->m_basisKeyVector[1];
                    LibUtilities::BasisKey TetBc
                                        = expIt.second->m_basisKeyVector[2];

                    if(TetBa.GetBasisType() == LibUtilities::eGLL_Lagrange ||
                       TetBa.GetBasisType() == LibUtilities::eGauss_Lagrange)
                    {
                      ASSERTL0(false,"LocalRegions::NodalTetExp is not "
                                     "implemented yet");
                    }
                    else
                    {
                        tet = MemoryManager<LocalRegions::TetExp>
                                        ::AllocateSharedPtr(TetBa,TetBb,TetBc,
                                                            TetGeom);
                        tet->SetElmtId(elmtid++);
                        (*m_exp).push_back(tet);
                    }
                }
                else if((PrismGeom = std::dynamic_pointer_cast<SpatialDomains
                             ::PrismGeom>(expIt.second->m_geomShPtr)))
                {
                    LibUtilities::BasisKey PrismBa
                                        = expIt.second->m_basisKeyVector[0];
                    LibUtilities::BasisKey PrismBb
                                        = expIt.second->m_basisKeyVector[1];
                    LibUtilities::BasisKey PrismBc
                                        = expIt.second->m_basisKeyVector[2];

                    prism = MemoryManager<LocalRegions::PrismExp>
                                        ::AllocateSharedPtr(PrismBa,PrismBb,
                                                            PrismBc,PrismGeom);
                    prism->SetElmtId(elmtid++);
                    (*m_exp).push_back(prism);
                }
                else if((PyrGeom = std::dynamic_pointer_cast<
                         SpatialDomains::PyrGeom>(expIt.second->m_geomShPtr)))
                {
                    LibUtilities::BasisKey PyrBa
                                        = expIt.second->m_basisKeyVector[0];
                    LibUtilities::BasisKey PyrBb
                                        = expIt.second->m_basisKeyVector[1];
                    LibUtilities::BasisKey PyrBc
                                        = expIt.second->m_basisKeyVector[2];

                    pyramid = MemoryManager<LocalRegions::PyrExp>
                                        ::AllocateSharedPtr(PyrBa,PyrBb,PyrBc,
                                                            PyrGeom);
                    pyramid->SetElmtId(elmtid++);
                    (*m_exp).push_back(pyramid);
                }
                else if((HexGeom = std::dynamic_pointer_cast<
                         SpatialDomains::HexGeom>(expIt.second->m_geomShPtr)))
                {
                    LibUtilities::BasisKey HexBa
                                        = expIt.second->m_basisKeyVector[0];
                    LibUtilities::BasisKey HexBb
                                        = expIt.second->m_basisKeyVector[1];
                    LibUtilities::BasisKey HexBc
                                        = expIt.second->m_basisKeyVector[2];

                    hex = MemoryManager<LocalRegions::HexExp>
                                        ::AllocateSharedPtr(HexBa,HexBb,HexBc,
                                                            HexGeom);
                    hex->SetElmtId(elmtid++);
                    (*m_exp).push_back(hex);
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry3D "
                                   "failed");
                }

            }

            SetCoeffPhysOffsets();

            // Set up m_coeffs, m_phys.
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs,0.0);
            m_phys   = Array<OneD, NekDouble>(m_npoints,0.0);

            CreateCollections(ImpType);
        }

        /**
         * Given an expansion vector \a expansions, containing
         * information about the domain and the spectral/hp element
         * expansion, this constructor fills the list of local
         * expansions \texttt{m_exp} with the proper expansions,
         * calculates the total number of quadrature points
         * \f$\boldsymbol{x}_i\f$ and the local expansion coefficients
         * \f$\hat{u}^e_n\f$ and allocates memory for the arrays
         * #m_coeffs and #m_phys.
         *
         * @param expansions An expansion vector, containing
         *                   information about the domain and the
         *                   spectral/hp element expansion.
         */
        ExpList3D::ExpList3D(const SpatialDomains::ExpansionMap &expansions,
                             const Collections::ImplementationType ImpType):
            ExpList(e3D)
        {
            int elmtid = 0;
            LocalRegions::TetExpSharedPtr tet;
            LocalRegions::HexExpSharedPtr hex;
            LocalRegions::PrismExpSharedPtr prism;
            LocalRegions::PyrExpSharedPtr pyramid;


            for(int i = 0; i < expansions.size(); ++i)
            {
                SpatialDomains::TetGeomSharedPtr TetGeom;
                SpatialDomains::HexGeomSharedPtr HexGeom;
                SpatialDomains::PrismGeomSharedPtr PrismGeom;
                SpatialDomains::PyrGeomSharedPtr PyrGeom;

                auto expmap = expansions.find(i);
                ASSERTL1(expmap != expansions.end(), "Unable to find expansion.");
                const SpatialDomains::ExpansionShPtr exp = expmap->second;

                if((TetGeom = std::dynamic_pointer_cast<
                        SpatialDomains::TetGeom>(exp->m_geomShPtr)))
                {
                    LibUtilities::BasisKey TetBa
                                        = exp->m_basisKeyVector[0];
                    LibUtilities::BasisKey TetBb
                                        = exp->m_basisKeyVector[1];
                    LibUtilities::BasisKey TetBc
                                        = exp->m_basisKeyVector[2];

                    if(TetBa.GetBasisType() == LibUtilities::eGLL_Lagrange ||
                       TetBa.GetBasisType() == LibUtilities::eGauss_Lagrange)
                    {
                      ASSERTL0(false,"LocalRegions::NodalTetExp is not "
                                     "implemented yet");
                    }
                    else
                    {
                        tet = MemoryManager<LocalRegions::TetExp>
                                        ::AllocateSharedPtr(TetBa,TetBb,TetBc,
                                                            TetGeom);
                        tet->SetElmtId(elmtid++);
                        (*m_exp).push_back(tet);
                    }
                }
                else if((PrismGeom = std::dynamic_pointer_cast<
                         SpatialDomains::PrismGeom>(exp->m_geomShPtr)))
                {
                    LibUtilities::BasisKey PrismBa
                                        = exp->m_basisKeyVector[0];
                    LibUtilities::BasisKey PrismBb
                                        = exp->m_basisKeyVector[1];
                    LibUtilities::BasisKey PrismBc
                                        = exp->m_basisKeyVector[2];

                    prism = MemoryManager<LocalRegions::PrismExp>
                                        ::AllocateSharedPtr(PrismBa,PrismBb,
                                                            PrismBc,PrismGeom);
                    prism->SetElmtId(elmtid++);
                    (*m_exp).push_back(prism);
                }
                else if((PyrGeom = std::dynamic_pointer_cast<
                         SpatialDomains::PyrGeom>(exp->m_geomShPtr)))
                {
                    LibUtilities::BasisKey PyrBa
                                        = exp->m_basisKeyVector[0];
                    LibUtilities::BasisKey PyrBb
                                        = exp->m_basisKeyVector[1];
                    LibUtilities::BasisKey PyrBc
                                        = exp->m_basisKeyVector[2];

                    pyramid = MemoryManager<LocalRegions::PyrExp>
                                        ::AllocateSharedPtr(PyrBa,PyrBb,PyrBc,
                                                            PyrGeom);
                    pyramid->SetElmtId(elmtid++);
                    (*m_exp).push_back(pyramid);
                }
                else if((HexGeom = std::dynamic_pointer_cast<
                         SpatialDomains::HexGeom>(exp->m_geomShPtr)))
                {
                    LibUtilities::BasisKey HexBa
                                        = exp->m_basisKeyVector[0];
                    LibUtilities::BasisKey HexBb
                                        = exp->m_basisKeyVector[1];
                    LibUtilities::BasisKey HexBc
                                        = exp->m_basisKeyVector[2];

                    hex = MemoryManager<LocalRegions::HexExp>
                                        ::AllocateSharedPtr(HexBa,HexBb,HexBc,
                                                            HexGeom);
                    hex->SetElmtId(elmtid++);
                    (*m_exp).push_back(hex);
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry3D "
                                   "failed");
                }

            }

            SetCoeffPhysOffsets();

            // Set up m_coeffs, m_phys.
            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs,0.0);
            m_phys   = Array<OneD, NekDouble>(m_npoints,0.0);

            CreateCollections(ImpType);
        }

  } //end of namespace
} //end of namespace

