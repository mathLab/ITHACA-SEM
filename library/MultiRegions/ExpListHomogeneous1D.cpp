///////////////////////////////////////////////////////////////////////////////
//
// File ExpListHomogeneous1D.cpp
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
// Description: An ExpList which is homogeneous in 1-direction
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpListHomogeneous1D.h>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declaration for typedefs
        ExpListHomogeneous1D::ExpListHomogeneous1D():
            ExpList(),
            m_homogeneousBasis(LibUtilities::NullBasisSharedPtr),
            m_lhom(1),
            m_homogeneous1DBlockMat(MemoryManager<Homo1DBlockMatrixMap>::AllocateSharedPtr())
        {
        }

        ExpListHomogeneous1D::ExpListHomogeneous1D(const LibUtilities::BasisKey &HomoBasis, const NekDouble lhom):
            ExpList(),
            m_lhom(lhom),
            m_homogeneous1DBlockMat(MemoryManager<Homo1DBlockMatrixMap>::AllocateSharedPtr())
        {
            ASSERTL2(m_homogeneousBasis != LibUtilities::NullBasisKey,
                     "Homogeneous Basis is not a null basis");
            m_homogeneousBasis = LibUtilities::BasisManager()[HomoBasis];

            int nzplanes = m_homogeneousBasis->GetNumPoints();

            m_planes = Array<OneD,ExpListSharedPtr>(nzplanes);
        }


        /**
         * @param   In          ExpListHomogeneous1D object to copy.
         */
        ExpListHomogeneous1D::ExpListHomogeneous1D(const ExpListHomogeneous1D &In):
            ExpList(In,false),
            m_homogeneousBasis(In.m_homogeneousBasis),
            m_homogeneous1DBlockMat(In.m_homogeneous1DBlockMat),
            m_lhom(In.m_lhom)
        {
            m_planes = Array<OneD, ExpListSharedPtr>(In.m_planes.num_elements());
        }

        /**
         * Destructor
         */
        ExpListHomogeneous1D::~ExpListHomogeneous1D()
        {
        }

        void ExpListHomogeneous1D::v_FwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            int nplanes = m_planes.num_elements();

            for(int n = 0; n < nplanes; ++n)
            {
                m_planes[n]->FwdTrans(inarray+cnt, tmparray = outarray + cnt1,
                                      UseContCoeffs);
                cnt   += m_planes[n]->GetTotPoints();

                if(UseContCoeffs)
                {
                    cnt1  += m_planes[n]->GetContNcoeffs();
                }
                else
                {
                    cnt1  += m_planes[n]->GetNcoeffs();
                }
            }

            Homogeneous1DFwdTrans(outarray,outarray,UseContCoeffs);
        }

        void ExpListHomogeneous1D::v_BwdTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            int nplanes = m_planes.num_elements();

            for(int n = 0; n < nplanes; ++n)
            {
                m_planes[n]->BwdTrans(inarray+cnt, tmparray = outarray + cnt1,
                                      UseContCoeffs);
                if(UseContCoeffs)
                {
                    cnt    += m_planes[n]->GetContNcoeffs();
                }
                else
                {
                    cnt    += m_planes[n]->GetNcoeffs();
                }
                cnt1   += m_planes[n]->GetTotPoints();
            }

                Homogeneous1DBwdTrans(outarray,outarray);
        }


        void ExpListHomogeneous1D::v_IProductWRTBase(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool UseContCoeffs)
        {
            int cnt = 0, cnt1 = 0;
            Array<OneD, NekDouble> tmparray;
            int nplanes = m_planes.num_elements();

            for(int n = 0; n < nplanes; ++n)
            {
                m_planes[n]->IProductWRTBase(inarray+cnt, tmparray = outarray + cnt1,UseContCoeffs);

                if(UseContCoeffs)
                {
                    cnt    += m_planes[n]->GetContNcoeffs();
                }
                else
                {
                    cnt    += m_planes[n]->GetNcoeffs();
                }
                cnt1   += m_planes[n]->GetTotPoints();
            }
        }

        void ExpListHomogeneous1D::Homogeneous1DTrans(const Array<OneD, const NekDouble> &inarray, Array<OneD, NekDouble> &outarray, bool IsForwards, bool UseContCoeffs)
        {
            DNekBlkMatSharedPtr blkmat;

            if(inarray.num_elements() == m_npoints) //transform phys space
            {
                if(IsForwards)
                {
                    blkmat = GetHomogeneous1DBlockMatrix(eForwardsPhysSpace);
                }
                else
                {
                    blkmat = GetHomogeneous1DBlockMatrix(eBackwardsPhysSpace);
                }
            }
            else
            {
                if(IsForwards)
                {
                    blkmat = GetHomogeneous1DBlockMatrix(eForwardsCoeffSpace,UseContCoeffs);
                }
                else
                {
                    blkmat = GetHomogeneous1DBlockMatrix(eBackwardsCoeffSpace,UseContCoeffs);
                }
            }

            int nrows = blkmat->GetRows();
            int ncols = blkmat->GetColumns();

            Array<OneD, NekDouble> sortedinarray(ncols);
            Array<OneD, NekDouble> sortedoutarray(nrows);


            ShuffleIntoHomogeneous1DClosePacked(inarray,sortedinarray,!IsForwards);

            // Create NekVectors from the given data arrays
            NekVector<const NekDouble> in (ncols,sortedinarray,eWrapper);
            NekVector<      NekDouble> out(nrows,sortedoutarray,eWrapper);

            // Perform matrix-vector multiply.
            out = (*blkmat)*in;

            UnshuffleFromHomogeneous1DClosePacked(sortedoutarray,outarray,IsForwards);
        }

        void ExpListHomogeneous1D::ShuffleIntoHomogeneous1DClosePacked(
                              const Array<OneD, const NekDouble> &inarray,
                              Array<OneD, NekDouble> &outarray,
                              bool UseNumModes)
        {
            int i, pts_per_plane;
            int n = inarray.num_elements();
            int packed_len;

            pts_per_plane = n/m_planes.num_elements();

            if(UseNumModes)
            {
                packed_len = m_homogeneousBasis->GetNumModes();
            }
            else
            {
                packed_len = m_planes.num_elements();
            }

            ASSERTL1(&inarray[0] != &outarray[0],"Inarray and outarray cannot be the same");

            for(i = 0; i < packed_len; ++i)
            {
                Vmath::Vcopy(pts_per_plane,&(inarray[i*pts_per_plane]),1,
                             &(outarray[i]),packed_len);
            }
        }

        void ExpListHomogeneous1D::UnshuffleFromHomogeneous1DClosePacked(
                              const Array<OneD, const NekDouble> &inarray,
                              Array<OneD, NekDouble> &outarray,
                              bool UseNumModes)
        {
            int i,pts_per_plane;
            int n = inarray.num_elements();
            int packed_len;

            // use length of inarray to determine data storage type
            // (i.e.modal or physical).
            pts_per_plane = n/m_planes.num_elements();

            if(UseNumModes)
            {
                packed_len = m_homogeneousBasis->GetNumModes();
            }
            else
            {
                packed_len = m_planes.num_elements();
            }

            ASSERTL1(&inarray[0] != &outarray[0],"Inarray and outarray cannot be the same");


            for(i = 0; i < packed_len; ++i)
            {
                Vmath::Vcopy(pts_per_plane,&(inarray[i]),packed_len,
                             &(outarray[i*pts_per_plane]),1);
            }
        }

        DNekBlkMatSharedPtr ExpListHomogeneous1D::GetHomogeneous1DBlockMatrix(Homogeneous1DMatType mattype, bool UseContCoeffs) const
        {
            Homo1DBlockMatrixMap::iterator matrixIter = m_homogeneous1DBlockMat->find(mattype);

            if(matrixIter == m_homogeneous1DBlockMat->end())
            {
                return ((*m_homogeneous1DBlockMat)[mattype] =
                        GenHomogeneous1DBlockMatrix(mattype,UseContCoeffs));
            }
            else
            {
                return matrixIter->second;
            }
        }


        DNekBlkMatSharedPtr ExpListHomogeneous1D::GenHomogeneous1DBlockMatrix(Homogeneous1DMatType mattype, bool UseContCoeffs) const
        {
            int i;
            int n_exp = 0;
            DNekMatSharedPtr    loc_mat;
            DNekBlkMatSharedPtr BlkMatrix;

            if((mattype == eForwardsCoeffSpace)
               ||(mattype == eBackwardsCoeffSpace)) // will operate on m_coeffs
            {
                if(UseContCoeffs)
                {
                    n_exp = m_planes[0]->GetContNcoeffs();
                }
                else
                {
                    n_exp = m_planes[0]->GetNcoeffs();
                }
            }
            else
            {
                n_exp = m_planes[0]->GetTotPoints(); // will operatore on m_phys
            }

            Array<OneD,unsigned int> nrows(n_exp);
            Array<OneD,unsigned int> ncols(n_exp);

            if((mattype == eForwardsCoeffSpace)||(mattype == eForwardsPhysSpace))
            {
                nrows = Array<OneD, unsigned int>(n_exp,m_homogeneousBasis->GetNumModes());
                ncols = Array<OneD, unsigned int>(n_exp,m_planes.num_elements());
            }
            else
            {
                nrows = Array<OneD, unsigned int>(n_exp,m_planes.num_elements());
                ncols = Array<OneD, unsigned int>(n_exp,m_homogeneousBasis->GetNumModes());
            }

            MatrixStorage blkmatStorage = eDIAGONAL;
            BlkMatrix = MemoryManager<DNekBlkMat>
                ::AllocateSharedPtr(nrows,ncols,blkmatStorage);

            StdRegions::StdSegExp StdSeg(m_homogeneousBasis->GetBasisKey());

            if((mattype == eForwardsCoeffSpace)||(mattype == eForwardsPhysSpace))
            {
                StdRegions::StdMatrixKey matkey(StdRegions::eFwdTrans,
                                                StdSeg.DetExpansionType(),
                                                StdSeg);

                loc_mat = StdSeg.GetStdMatrix(matkey);
            }
            else
            {
                StdRegions::StdMatrixKey matkey(StdRegions::eBwdTrans,
                                                StdSeg.DetExpansionType(),
                                                StdSeg);

                loc_mat = StdSeg.GetStdMatrix(matkey);
            }

            // set up array of block matrices.
            for(i = 0; i < n_exp; ++i)
            {
                BlkMatrix->SetBlock(i,i,loc_mat);
            }

            return BlkMatrix;
        }

        std::vector<SpatialDomains::FieldDefinitionsSharedPtr> ExpListHomogeneous1D::v_GetFieldDefinitions()
        {
            std::vector<SpatialDomains::FieldDefinitionsSharedPtr> returnval;
            // Set up Homogeneous length details.
            Array<OneD,LibUtilities::BasisSharedPtr> HomoBasis(1,m_homogeneousBasis);
            std::vector<NekDouble> HomoLen;
            HomoLen.push_back(m_lhom);

            m_planes[0]->GeneralGetFieldDefinitions(returnval, 1, HomoBasis, HomoLen);
            return returnval;
        }

        void  ExpListHomogeneous1D::v_GetFieldDefinitions(std::vector<SpatialDomains::FieldDefinitionsSharedPtr> &fielddef)
        {
            // Set up Homogeneous length details.
            Array<OneD,LibUtilities::BasisSharedPtr> HomoBasis(1,m_homogeneousBasis);
            std::vector<NekDouble> HomoLen;
            HomoLen.push_back(m_lhom);

             // enforce NumHomoDir == 1 by direct call
            m_planes[0]->GeneralGetFieldDefinitions(fielddef,1, HomoBasis,HomoLen);
        }


        void ExpListHomogeneous1D::v_AppendFieldData(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata)
        {
            int i,n;
            int nzmodes = m_homogeneousBasis->GetNumModes();
            int ncoeffs_per_plane = m_planes[0]->GetNcoeffs();

            // Determine mapping from element ids to location in
            // expansion list
            map<int, int> ElmtID_to_ExpID;
            for(i = 0; i < m_planes[0]->GetExpSize(); ++i)
            {
                ElmtID_to_ExpID[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
            }

            for(i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                int eid     = ElmtID_to_ExpID[fielddef->m_elementIDs[i]];
                int datalen = (*m_exp)[eid]->GetNcoeffs();

                for(n = 0; n < nzmodes; ++n)
                {
                    fielddata.insert(fielddata.end(),&m_coeffs[m_coeff_offset[eid]+n*ncoeffs_per_plane],&m_coeffs[m_coeff_offset[eid]+n*ncoeffs_per_plane]+datalen);
                }
            }
        }

        //Extract the data in fielddata into the m_coeff list
        void ExpListHomogeneous1D::v_ExtractDataToCoeffs(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, std::string &field)
        {
            int i,n;
            int offset = 0;
            int nzmodes = m_homogeneousBasis->GetNumModes();
            int datalen = fielddata.size()/fielddef->m_fields.size();
            int ncoeffs_per_plane = m_planes[0]->GetNcoeffs();

            // Find data location according to field definition
            for(i = 0; i < fielddef->m_fields.size(); ++i)
            {
                if(fielddef->m_fields[i] == field)
                {
                    break;
                }
                offset += datalen;
            }

            ASSERTL0(i!= fielddef->m_fields.size(),"Field not found in data file");

            // Determine mapping from element ids to location in
            // expansion list
            map<int, int> ElmtID_to_ExpID;
            for(i = 0; i < m_planes[0]->GetExpSize(); ++i)
            {
                ElmtID_to_ExpID[(*m_exp)[i]->GetGeom()->GetGlobalID()] = i;
            }

            for(i = 0; i < fielddef->m_elementIDs.size(); ++i)
            {
                int eid = ElmtID_to_ExpID[fielddef->m_elementIDs[i]];
                int datalen = (*m_exp)[eid]->GetNcoeffs();
                for(n = 0; n < nzmodes; ++n)
                {
                    Vmath::Vcopy(datalen,&fielddata[offset],1,&m_coeffs[m_coeff_offset[eid] + n*ncoeffs_per_plane],1);
                    offset += datalen;
                }
            }
        }


        /**
         * Write Tecplot Files Header
         * @param   outfile Output file name.
         * @param   var                 variables names
         */
        void ExpListHomogeneous1D::v_WriteTecplotHeader(std::ofstream &outfile, std::string var)
        {
            if(GetExp(0)->GetCoordim() == 1)
            {
                outfile << "Variables = x, y";
            }
            else
            {
                outfile << "Variables = x, y, z";
            }
            outfile << ", "<< var << std::endl << std::endl;
        }

        /**
         * Write Tecplot Files Field
         * @param   outfile    Output file name.
         * @param   expansion  Expansion that is considered
         */
        void ExpListHomogeneous1D::v_WriteTecplotField(std::ofstream &outfile, int expansion)
        {
            int npoints_per_plane = m_planes[0]->GetTotPoints();

            for(int n = 0; n < m_planes.num_elements(); ++n)
            {
                (*m_exp)[expansion]->SetPhys(m_phys+m_phys_offset[expansion]+
                                             n*npoints_per_plane);
                (*m_exp)[expansion]->WriteTecplotField(outfile);
            }
        }

        void ExpListHomogeneous1D::v_WriteVtkPieceData(std::ofstream &outfile, int expansion,
                                        std::string var)
        {
            int i;
            int nq = (*m_exp)[expansion]->GetTotPoints();
            int npoints_per_plane = m_planes[0]->GetTotPoints();

            // printing the fields of that zone
            outfile << "        <DataArray type=\"Float32\" Name=\""
                    << var << "\">" << endl;
            outfile << "          ";
            for (int n = 0; n < m_planes.num_elements(); ++n)
            {
                const Array<OneD, NekDouble> phys = m_phys + m_phys_offset[expansion] + n*npoints_per_plane;
                for(i = 0; i < nq; ++i)
                {
                    outfile << (fabs(phys[i]) < NekConstants::kNekZeroTol ? 0 : phys[i]) << " ";
                }
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
        }


   } //end of namespace
} //end of namespace


/**
* $Log: v $
*
**/

