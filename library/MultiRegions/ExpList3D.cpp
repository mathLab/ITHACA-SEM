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

#include <boost/core/ignore_unused.hpp>

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

        ExpList3D::ExpList3D(): ExpList()
        {
            SetExpType(e3D);
        }

        ExpList3D::ExpList3D(const ExpList3D &In): ExpList(In)
        {
            SetExpType(e3D);
        }
        
        ExpList3D::ExpList3D(const ExpList3D &In,
                const std::vector<unsigned int> &eIDs,
                const bool DeclareCoeffPhysArrays,
                const Collections::ImplementationType ImpType):
            ExpList(In, eIDs, DeclareCoeffPhysArrays)
        {
            SetExpType(e3D);
            
            SetCoeffPhysOffsets();

            if (DeclareCoeffPhysArrays)
            {
                // Set up m_coeffs, m_phys.
                m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
                m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};
             }

            CreateCollections(ImpType);
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
            ExpList(pSession,graph3D)
        {
            SetExpType(e3D);

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
            m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
            m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};

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
            ExpList(pSession,graph3D)
        {
            SetExpType(e3D);

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
            m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
            m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};

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
            ExpList()
        {
            SetExpType(e3D);

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
            m_coeffs = Array<OneD, NekDouble> {size_t(m_ncoeffs), 0.0};
            m_phys   = Array<OneD, NekDouble> {size_t(m_npoints), 0.0};

            CreateCollections(ImpType);
        }

        void ExpList3D::v_WriteVtkPieceHeader(std::ostream &outfile, int expansion, int istrip)
        {
            boost::ignore_unused(istrip);

            int i,j,k;
            int nquad0 = (*m_exp)[expansion]->GetNumPoints(0);
            int nquad1 = (*m_exp)[expansion]->GetNumPoints(1);
            int nquad2 = (*m_exp)[expansion]->GetNumPoints(2);
            int ntot = nquad0*nquad1*nquad2;
            int ntotminus = (nquad0-1)*(nquad1-1)*(nquad2-1);

            Array<OneD,NekDouble> coords[3];
            coords[0] = Array<OneD,NekDouble>(ntot);
            coords[1] = Array<OneD,NekDouble>(ntot);
            coords[2] = Array<OneD,NekDouble>(ntot);
            (*m_exp)[expansion]->GetCoords(coords[0],coords[1],coords[2]);

            outfile << "    <Piece NumberOfPoints=\""
                    << ntot << "\" NumberOfCells=\""
                    << ntotminus << "\">" << endl;
            outfile << "      <Points>" << endl;
            outfile << "        <DataArray type=\"Float64\" "
                    << "NumberOfComponents=\"3\" format=\"ascii\">" << endl;
            outfile << "          ";
            for (i = 0; i < ntot; ++i)
            {
                for (j = 0; j < 3; ++j)
                {
                    outfile << setprecision(8) << scientific 
                            << (float)coords[j][i] << " ";
                }
                outfile << endl;
            }

            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "      </Points>" << endl;
            outfile << "      <Cells>" << endl;
            outfile << "        <DataArray type=\"Int32\" "
                    << "Name=\"connectivity\" format=\"ascii\">" << endl;
            for (i = 0; i < nquad0-1; ++i)
            {
                for (j = 0; j < nquad1-1; ++j)
                {
                    for (k = 0; k < nquad2-1; ++k)
                    {
                        outfile << k*nquad0*nquad1 + j*nquad0 + i << " "
                                << k*nquad0*nquad1 + j*nquad0 + i + 1 << " "
                                << k*nquad0*nquad1 + (j+1)*nquad0 + i + 1 << " "
                                << k*nquad0*nquad1 + (j+1)*nquad0 + i << " "
                                << (k+1)*nquad0*nquad1 + j*nquad0 + i << " "
                                << (k+1)*nquad0*nquad1 + j*nquad0 + i + 1 << " "
                                << (k+1)*nquad0*nquad1 + (j+1)*nquad0 + i + 1 << " "
                                << (k+1)*nquad0*nquad1 + (j+1)*nquad0 + i << " " << endl;
                    }
                }
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "        <DataArray type=\"Int32\" "
                    << "Name=\"offsets\" format=\"ascii\">" << endl;
            for (i = 0; i < ntotminus; ++i)
            {
                outfile << i*8+8 << " ";
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "        <DataArray type=\"UInt8\" "
                    << "Name=\"types\" format=\"ascii\">" << endl;
            for (i = 0; i < ntotminus; ++i)
            {
                outfile << "12 ";
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
            outfile << "      </Cells>" << endl;
            outfile << "      <PointData>" << endl;
        }

        void ExpList3D::v_SetUpPhysNormals()
        {
            int i, j;
            for (i = 0; i < m_exp->size(); ++i)
            {
                for (j = 0; j < (*m_exp)[i]->GetNfaces(); ++j)
                {
                    (*m_exp)[i]->ComputeFaceNormal(j);
                }
            }
        }

        void ExpList3D::v_PhysInterp1DScaled(const NekDouble scale, 
                                  const Array<OneD, NekDouble> &inarray, 
                                  Array<OneD, NekDouble> &outarray)
        {
            int cnt,cnt1;

            cnt = cnt1 = 0;
            for(int i = 0; i < GetExpSize(); ++i)
            {
                // get new points key
                int pt0 = (*m_exp)[i]->GetNumPoints(0);
                int pt1 = (*m_exp)[i]->GetNumPoints(1);
                int pt2 = (*m_exp)[i]->GetNumPoints(2);
                int npt0 = (int) pt0*scale;
                int npt1 = (int) pt1*scale;
                int npt2 = (int) pt2*scale;
                
                LibUtilities::PointsKey newPointsKey0(npt0,(*m_exp)[i]->GetPointsType(0));
                LibUtilities::PointsKey newPointsKey1(npt1,(*m_exp)[i]->GetPointsType(1));
                LibUtilities::PointsKey newPointsKey2(npt2,(*m_exp)[i]->GetPointsType(2));

                // Interpolate points; 
                LibUtilities::Interp3D((*m_exp)[i]->GetBasis(0)->GetPointsKey(),
                                       (*m_exp)[i]->GetBasis(1)->GetPointsKey(),
                                       (*m_exp)[i]->GetBasis(2)->GetPointsKey(),
                                       &inarray[cnt], newPointsKey0,
                                       newPointsKey1, newPointsKey2,
                                       &outarray[cnt1]);

                cnt  += pt0*pt1*pt2;
                cnt1 += npt0*npt1*npt2;
            }
        }
        
        void ExpList3D::v_PhysGalerkinProjection1DScaled(const NekDouble scale, 
                                           const Array<OneD, NekDouble> &inarray,
                                           Array<OneD, NekDouble> &outarray)
        {
            int cnt,cnt1;

            cnt = cnt1 = 0;
            for(int i = 0; i < GetExpSize(); ++i)
            {
                // get new points key
                int pt0 = (*m_exp)[i]->GetNumPoints(0);
                int pt1 = (*m_exp)[i]->GetNumPoints(1);
                int pt2 = (*m_exp)[i]->GetNumPoints(2);
                int npt0 = (int) pt0*scale;
                int npt1 = (int) pt1*scale;
                int npt2 = (int) pt2*scale;
                
                LibUtilities::PointsKey newPointsKey0(npt0,(*m_exp)[i]->GetPointsType(0));
                LibUtilities::PointsKey newPointsKey1(npt1,(*m_exp)[i]->GetPointsType(1));
                LibUtilities::PointsKey newPointsKey2(npt2,(*m_exp)[i]->GetPointsType(2));

                // Project points; 
                LibUtilities::PhysGalerkinProject3D(newPointsKey0, 
                                                    newPointsKey1,
                                                    newPointsKey2,
                                                    &inarray[cnt],
                                       (*m_exp)[i]->GetBasis(0)->GetPointsKey(),
                                       (*m_exp)[i]->GetBasis(1)->GetPointsKey(),
                                       (*m_exp)[i]->GetBasis(2)->GetPointsKey(),
                                       &outarray[cnt1]);
                
                cnt  += npt0*npt1*npt2;
                cnt1 += pt0*pt1*pt2;
            }

        }
  } //end of namespace
} //end of namespace

