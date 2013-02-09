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


        ExpList3D::ExpList3D(const LibUtilities::SessionReaderSharedPtr &pSession,
                             const LibUtilities::BasisKey &TBa,
                             const LibUtilities::BasisKey &TBb,
                             const LibUtilities::BasisKey &TBc,
                             const LibUtilities::BasisKey &HBa,
                             const LibUtilities::BasisKey &HBb,
                             const LibUtilities::BasisKey &HBc,
                             const SpatialDomains::MeshGraphSharedPtr &graph3D,
                             const LibUtilities::PointsType TetNb):
            ExpList(pSession,graph3D)
        {

            LocalRegions::TetExpSharedPtr tet;
            LocalRegions::HexExpSharedPtr hex;
            LocalRegions::PrismExpSharedPtr prism;
            LocalRegions::PyrExpSharedPtr pyramid;

            const SpatialDomains::ExpansionMap &expansions = graph3D->GetExpansions();

            SpatialDomains::ExpansionMap::const_iterator expIt;
            for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
            {
                SpatialDomains::TetGeomSharedPtr TetGeom;
                SpatialDomains::HexGeomSharedPtr HexGeom;
                SpatialDomains::PrismGeomSharedPtr PrismGeom;
                SpatialDomains::PyrGeomSharedPtr PyrGeom;

                if((TetGeom = boost::dynamic_pointer_cast<SpatialDomains::TetGeom>(expIt->second->m_geomShPtr)))
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

                    m_ncoeffs += StdRegions::StdTetData::getNumberOfCoefficients(TBa.GetNumModes(), TBb.GetNumModes(), TBc.GetNumModes());
                    
                    m_npoints += TBa.GetNumPoints()*TBb.GetNumPoints()*TBc.GetNumPoints();
                }
/*
                else if((PrismGeom = boost::dynamic_pointer_cast<SpatialDomains::PrismGeom>(expansions[i]->m_geomShPtr)))
                {
                      prism = MemoryManager<LocalRegions::PrismExp>::AllocateSharedPtr(Ba,Bb,Bc,PrismGeom);
                      (*m_exp).push_back(prism);

                      m_ncoeffs += StdRegions::StdPrismData::getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes());
                      m_npoints +=  Ba.GetNumPoints()*Bb.GetNumPoints()*Bc.GetNumPoints();

                }
                else if((PyrGeom = boost::dynamic_pointer_cast<SpatialDomains::PyrGeom>(expansions[i]->m_geomShPtr)))
                {
                     pyramid = MemoryManager<LocalRegions::PyrExp>::AllocateSharedPtr(Ba,Bb,Bc,PyrGeom);
                     (*m_exp).push_back(pyramid);

                     m_ncoeffs += StdRegions::StdPyrData::getNumberOfCoefficients(Ba.GetNumModes(), Bb.GetNumModes(), Bc.GetNumModes());
                      m_npoints +=  Ba.GetNumPoints()*Bb.GetNumPoints()*Bc.GetNumPoints();

                }
*/
                else if((HexGeom = boost::dynamic_pointer_cast<SpatialDomains::HexGeom>(expIt->second->m_geomShPtr)))
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

            // Setup Default optimisation information.
            int nel = GetExpSize();
            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                ::AllocateSharedPtr(nel);

            SetCoeffPhys();

            ReadGlobalOptimizationParameters();
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
                const SpatialDomains::MeshGraphSharedPtr &graph3D) :
            ExpList(pSession,graph3D)
        {
            LocalRegions::TetExpSharedPtr   tet;
            LocalRegions::HexExpSharedPtr   hex;
            LocalRegions::PrismExpSharedPtr prism;
            LocalRegions::PyrExpSharedPtr   pyramid;

            const SpatialDomains::ExpansionMap &expansions
                                        = graph3D->GetExpansions();

            SpatialDomains::ExpansionMap::const_iterator expIt;
            for (expIt = expansions.begin(); expIt != expansions.end(); ++expIt)
            {
                SpatialDomains::TetGeomSharedPtr   TetGeom;
                SpatialDomains::HexGeomSharedPtr   HexGeom;
                SpatialDomains::PrismGeomSharedPtr PrismGeom;
                SpatialDomains::PyrGeomSharedPtr   PyrGeom;

                if((TetGeom = boost::dynamic_pointer_cast<
                        SpatialDomains::TetGeom>(expIt->second->m_geomShPtr)))
                {
                    LibUtilities::BasisKey TetBa
                                        = expIt->second->m_basisKeyVector[0];
                    LibUtilities::BasisKey TetBb
                                        = expIt->second->m_basisKeyVector[1];
                    LibUtilities::BasisKey TetBc
                                        = expIt->second->m_basisKeyVector[2];

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
                        (*m_exp).push_back(tet);
                    }
                }
                else if((PrismGeom = boost::dynamic_pointer_cast<SpatialDomains
                             ::PrismGeom>(expIt->second->m_geomShPtr)))
                {
                    LibUtilities::BasisKey PrismBa
                                        = expIt->second->m_basisKeyVector[0];
                    LibUtilities::BasisKey PrismBb
                                        = expIt->second->m_basisKeyVector[1];
                    LibUtilities::BasisKey PrismBc
                                        = expIt->second->m_basisKeyVector[2];

                    prism = MemoryManager<LocalRegions::PrismExp>
                                        ::AllocateSharedPtr(PrismBa,PrismBb,
                                                            PrismBc,PrismGeom);
                    (*m_exp).push_back(prism);
                }
                else if((PyrGeom = boost::dynamic_pointer_cast<
                         SpatialDomains::PyrGeom>(expIt->second->m_geomShPtr)))
                {
                    LibUtilities::BasisKey PyrBa
                                        = expIt->second->m_basisKeyVector[0];
                    LibUtilities::BasisKey PyrBb
                                        = expIt->second->m_basisKeyVector[1];
                    LibUtilities::BasisKey PyrBc
                                        = expIt->second->m_basisKeyVector[2];

                    pyramid = MemoryManager<LocalRegions::PyrExp>
                                        ::AllocateSharedPtr(PyrBa,PyrBb,PyrBc,
                                                            PyrGeom);
                    (*m_exp).push_back(pyramid);
                }
                else if((HexGeom = boost::dynamic_pointer_cast<
                         SpatialDomains::HexGeom>(expIt->second->m_geomShPtr)))
                {
                    LibUtilities::BasisKey HexBa
                                        = expIt->second->m_basisKeyVector[0];
                    LibUtilities::BasisKey HexBb
                                        = expIt->second->m_basisKeyVector[1];
                    LibUtilities::BasisKey HexBc
                                        = expIt->second->m_basisKeyVector[2];

                    hex = MemoryManager<LocalRegions::HexExp>
                                        ::AllocateSharedPtr(HexBa,HexBb,HexBc,
                                                            HexGeom);
                    (*m_exp).push_back(hex);
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry3D "
                                   "failed");
                }

            }

            // Setup Default optimisation information.
            int nel = GetExpSize();
            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                ::AllocateSharedPtr(nel);

            SetCoeffPhys();
            ReadGlobalOptimizationParameters();
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
        ExpList3D::ExpList3D(const SpatialDomains::ExpansionMap &expansions):
            ExpList()
        {
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

                SpatialDomains::ExpansionMap::const_iterator expmap = expansions.find(i);
                ASSERTL1(expmap != expansions.end(), "Unable to find expansion.");
                const SpatialDomains::ExpansionShPtr exp = expmap->second;

                if((TetGeom = boost::dynamic_pointer_cast<
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
                        (*m_exp).push_back(tet);
                    }
                }
                else if((PrismGeom = boost::dynamic_pointer_cast<
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
                    (*m_exp).push_back(prism);
                }
                else if((PyrGeom = boost::dynamic_pointer_cast<
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
                    (*m_exp).push_back(pyramid);
                }
                else if((HexGeom = boost::dynamic_pointer_cast<
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
                    (*m_exp).push_back(hex);
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a proper Geometry3D "
                                   "failed");
                }

            }

            // Setup Default optimisation information.
            int nel = GetExpSize();
            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                ::AllocateSharedPtr(nel);

            SetCoeffPhys();
        }

        /**
         * Set up the storage for the concatenated list of
         * coefficients and physical evaluations at the quadrature
         * points. Each expansion (local element) is processed in turn
         * to determine the number of coefficients and physical data
         * points it contributes to the domain. Three arrays,
         * #m_coeff_offset, #m_phys_offset and #m_offset_elmt_id, are
         * also initialised and updated to store the data offsets of
         * each element in the #m_coeffs and #m_phys arrays, and the
         * element id that each consecutive block is associated
         * respectively.
         */
        void ExpList3D::SetCoeffPhys()
        {
            int i;

            // Set up offset information and array sizes
            m_coeff_offset   = Array<OneD,int>(m_exp->size());
            m_phys_offset    = Array<OneD,int>(m_exp->size());
            m_offset_elmt_id = Array<OneD,int>(m_exp->size());

            m_ncoeffs = m_npoints = 0;

            for(i = 0; i < m_exp->size(); ++i)
            {
                m_coeff_offset[i]   = m_ncoeffs;
                m_phys_offset [i]   = m_npoints;
                m_offset_elmt_id[i] = i;
                m_ncoeffs += (*m_exp)[i]->GetNcoeffs();
                m_npoints += (*m_exp)[i]->GetTotPoints();
            }

            m_coeffs = Array<OneD, NekDouble>(m_ncoeffs);
            m_phys   = Array<OneD, NekDouble>(m_npoints);
        }

        void ExpList3D::v_ReadGlobalOptimizationParameters()
        {
            Array<OneD, int> NumShape(4,0);

            for(int i = 0; i < GetExpSize(); ++i)
            {
                switch ((*m_exp)[i]->DetExpansionType())
                {
                    case StdRegions::eTetrahedron:  NumShape[0]++; break;
                    case StdRegions::ePyramid:      NumShape[1]++; break;
                    case StdRegions::ePrism:        NumShape[2]++; break;
                    case StdRegions::eHexahedron:   NumShape[3]++; break;
                    default:
                        ASSERTL0(false, "Unknown expansion type.");
                        break;
                }
            }

            int three = 3;
            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                ::AllocateSharedPtr(m_session,three,NumShape);
        }

        void ExpList3D::v_WriteVtkPieceHeader(std::ofstream &outfile, int expansion)
        {
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
            outfile << "        <DataArray type=\"Float32\" "
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

  } //end of namespace
} //end of namespace

