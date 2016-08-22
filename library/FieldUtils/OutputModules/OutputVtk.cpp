////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputVtk.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: VTK file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "OutputVtk.h"
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <boost/format.hpp>
#include <iomanip>

namespace Nektar
{
namespace FieldUtils
{

ModuleKey OutputVtk::m_className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eOutputModule, "vtu"), OutputVtk::create, "Writes a VTU file.");

OutputVtk::OutputVtk(FieldSharedPtr f) : OutputModule(f)
{
    m_requireEquiSpaced = true;
}

OutputVtk::~OutputVtk()
{
}

void OutputVtk::Process(po::variables_map &vm)
{
    LibUtilities::PtsFieldSharedPtr fPts = m_f->m_fieldPts;

    // Do nothing if no expansion defined
    if (fPts == LibUtilities::NullPtsField && !m_f->m_exp.size())
    {
        return;
    }

    int i, j;
    if (m_f->m_verbose)
    {
        if (m_f->m_comm->TreatAsRankZero())
        {
            cout << "OutputVtk: Writing file..." << endl;
        }
    }

    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();
    string path;

    // amend for parallel output if required
    if (m_f->m_session->GetComm()->GetSize() != 1)
    {
        int dot      = filename.find_last_of('.');
        string ext   = filename.substr(dot, filename.length() - dot);
        string start = filename.substr(0, dot);
        path         = start + "_vtu";

        boost::format pad("P%1$07d.vtu");
        pad % m_f->m_session->GetComm()->GetRank();
        filename = pad.str();

        fs::path poutfile(filename.c_str());
        fs::path specPath(path.c_str());

        if (m_f->m_comm->TreatAsRankZero())
        {
            try
            {
                fs::create_directory(specPath);
            }
            catch (fs::filesystem_error &e)
            {
                ASSERTL0(false, "Filesystem error: " + string(e.what()));
            }
            cout << "Writing files to directory: " << specPath << endl;
        }

        fs::path fulloutname = specPath / poutfile;
        filename             = LibUtilities::PortablePath(fulloutname);
        m_f->m_comm->Block();
    }
    else
    {
        fs::path specPath(filename.c_str());
        cout << "Writing: " << specPath << endl;
        filename = LibUtilities::PortablePath(specPath);
    }

    // Write solution.
    ofstream outfile(filename.c_str());
    m_f->m_exp[0]->WriteVtkHeader(outfile);
    int nfields = 0;
    int dim     = 0;

    vector<string> fieldname;
    if (fPts ==
        LibUtilities::NullPtsField) // standard output in collapsed coordinates
    {
        int nstrips;
        if (m_f->m_fielddef.size() == 0)
        {
            nfields = 0;
        }
        else
        {
            nfields = m_f->m_fielddef[0]->m_fields.size();
        }
        m_f->m_session->LoadParameter("Strip_Z", nstrips, 1);

        // Homogeneous strip variant
        for (int s = 0; s < nstrips; ++s)
        {
            // For each field write out field data for each expansion.
            for (i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
            {
                m_f->m_exp[0]->WriteVtkPieceHeader(outfile, i, s);

                // For this expansion write out each field.
                for (j = 0; j < nfields; ++j)
                {
                    m_f->m_exp[s * nfields + j]->WriteVtkPieceData(
                        outfile, i, m_f->m_fielddef[0]->m_fields[j]);
                }
                m_f->m_exp[0]->WriteVtkPieceFooter(outfile, i);
            }
        }

        if (m_f->m_exp[0]->GetNumElmts() == 0)
        {
            WriteEmptyVtkPiece(outfile);
        }
        // save field names for parallel output
        for (i = 0; i < nfields; ++i)
        {
            fieldname.push_back(m_f->m_fielddef[0]->m_fields[i]);
        }
    }
    else // write out data stored in fPts (for example if equispaced output is
         // called).
    {
        int i = 0;
        int j = 0;

        dim = fPts->GetDim();

        int nvert, vtktype;
        switch (fPts->GetPtsType())
        {
            case LibUtilities::ePtsFile:
            case LibUtilities::ePtsLine:
            {
                ASSERTL0(false,
                    "VTK output needs setting up for ePtsFile and ePtsLine");
                break;
            }
            case LibUtilities::ePtsPlane:
            {
                ASSERTL0(false, "VTK output needs settig up for PtsPlane");
                break;
            }
            case LibUtilities::ePtsBox:
            {
                ASSERTL0(false,"VTK output needs settig up for PtsBox");
                break;
            }
            case LibUtilities::ePtsSegBlock:
            {
                nvert = 2;
                vtktype = 3;
                break;
            }
            case LibUtilities::ePtsTriBlock:
            {
                nvert   = 3;
                vtktype = 5;
                break;
            }
            case LibUtilities::ePtsTetBlock:
            {
                nvert   = 4;
                vtktype = 10;
                break;
            }
            default:
                ASSERTL0(false, "ptsType not supported yet.");
        }

        vector<Array<OneD, int> > ptsConn;
        fPts->GetConnectivity(ptsConn);

        nfields = fPts->GetNFields();

        int nPts      = fPts->GetNpoints();
        int numBlocks = 0;
        for (i = 0; i < ptsConn.size(); ++i)
        {
            numBlocks += ptsConn[i].num_elements() / nvert;
        }

        // write out pieces of data.
        outfile << "    <Piece NumberOfPoints=\"" << nPts
                << "\" NumberOfCells=\"" << numBlocks << "\">" << endl;
        outfile << "      <Points>" << endl;
        outfile << "        <DataArray type=\"Float64\" "
                << "NumberOfComponents=\"" << 3 << "\" format=\"ascii\">"
                << endl;
        for (i = 0; i < nPts; ++i)
        {
            for (j = 0; j < dim; ++j)
            {
                outfile << "          " << setprecision(8) << scientific
                        << fPts->GetPointVal(j, i) << " ";
            }
            for (j = dim; j < 3;
                 ++j) // pack to 3D since paraview does not seem to handle 2D
            {
                outfile << "          0.000000";
            }
            outfile << endl;
        }
        outfile << "        </DataArray>" << endl;
        outfile << "      </Points>" << endl;
        outfile << "      <Cells>" << endl;
        outfile << "        <DataArray type=\"Int32\" "
                << "Name=\"connectivity\" format=\"ascii\">" << endl;

        // dump connectivity data if it exists
        outfile << "          ";
        int cnt = 1;
        for (i = 0; i < ptsConn.size(); ++i)
        {
            for (j = 0; j < ptsConn[i].num_elements(); ++j)
            {
                outfile << ptsConn[i][j] << " ";
                if ((!(cnt % nvert)) && cnt)
                {
                    outfile << std::endl;
                    outfile << "          ";
                }
                cnt++;
            }
        }
        outfile << "        </DataArray>" << endl;
        outfile << "        <DataArray type=\"Int32\" "
                << "Name=\"offsets\" format=\"ascii\">" << endl;

        outfile << "          ";
        for (i = 0; i < numBlocks; ++i)
        {
            outfile << i * nvert + nvert << " ";
        }
        outfile << endl;
        outfile << "        </DataArray>" << endl;
        outfile << "        <DataArray type=\"UInt8\" "
                << "Name=\"types\" format=\"ascii\">" << endl;
        outfile << "          ";
        for (i = 0; i < numBlocks; ++i)
        {
            outfile << vtktype << " ";
        }
        outfile << endl;
        outfile << "        </DataArray>" << endl;
        outfile << "      </Cells>" << endl;
        outfile << "      <PointData>" << endl;

        // printing the fields
        for (j = 0; j < nfields; ++j)
        {
            fieldname.push_back(fPts->GetFieldName(j));
            outfile << "        <DataArray type=\"Float64\" Name=\""
                    << fPts->GetFieldName(j) << "\">" << endl;
            outfile << "          ";
            for (i = 0; i < fPts->GetNpoints(); ++i)
            {
                outfile << fPts->GetPointVal(dim + j, i) << " ";
            }
            outfile << endl;
            outfile << "        </DataArray>" << endl;
        }

        outfile << "      </PointData>" << endl;
        outfile << "    </Piece>" << endl;
    }

    m_f->m_exp[0]->WriteVtkFooter(outfile);
    cout << "Written file: " << filename << endl;

    // output parallel outline info if necessary
    if (m_f->m_comm->GetRank() == 0)
    {
        ASSERTL1(fieldname.size() == nfields,
                 "fieldname not the same size as nfields");
        int nprocs = m_f->m_comm->GetSize();
        if (nprocs != 1)
        {
            filename    = m_config["outfile"].as<string>();
            int dot     = filename.find_last_of('.');
            string body = filename.substr(0, dot);
            filename    = body + ".pvtu";

            ofstream outfile(filename.c_str());

            outfile << "<?xml version=\"1.0\"?>" << endl;
            outfile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
                    << "byte_order=\"LittleEndian\">" << endl;
            outfile << "<PUnstructuredGrid GhostLevel=\"0\">" << endl;
            outfile << "<PPoints> " << endl;
            outfile << "<PDataArray type=\"Float64\" NumberOfComponents=\"" << 3
                    << "\"/> " << endl;
            outfile << "</PPoints>" << endl;
            outfile << "<PCells>" << endl;
            outfile << "<PDataArray type=\"Int32\" Name=\"connectivity\" "
                       "NumberOfComponents=\"1\"/>"
                    << endl;
            outfile << "<PDataArray type=\"Int32\" Name=\"offsets\"      "
                       "NumberOfComponents=\"1\"/>"
                    << endl;
            outfile << "<PDataArray type=\"UInt8\" Name=\"types\"        "
                       "NumberOfComponents=\"1\"/>"
                    << endl;
            outfile << "</PCells>" << endl;
            outfile << "<PPointData Scalars=\"Material\">" << endl;
            for (int i = 0; i < nfields; ++i)
            {
                outfile << "<PDataArray type=\"Float64\" Name=\""
                        << fieldname[i] << "\"/>" << endl;
            }
            outfile << "</PPointData>" << endl;

            for (int i = 0; i < nprocs; ++i)
            {
                boost::format pad("P%1$07d.vtu");
                pad % i;
                outfile << "<Piece Source=\"" << path << "/" << pad.str()
                        << "\"/>" << endl;
            }
            outfile << "</PUnstructuredGrid>" << endl;
            outfile << "</VTKFile>" << endl;
            cout << "Written file: " << filename << endl;
        }
    }
}

void OutputVtk::WriteEmptyVtkPiece(std::ofstream &outfile)
{
    // write out empty piece of data.
    outfile << "    <Piece NumberOfPoints=\"" << 0 << "\" NumberOfCells=\"" << 0
            << "\">" << endl;
    outfile << "      <Points>" << endl;
    outfile << "        <DataArray type=\"Float64\" "
            << "NumberOfComponents=\"" << 3 << "\" format=\"ascii\">" << endl;
    outfile << "        </DataArray>" << endl;
    outfile << "      </Points>" << endl;
    outfile << "      <Cells>" << endl;
    outfile << "        <DataArray type=\"Int32\" "
            << "Name=\"connectivity\" format=\"ascii\">" << endl;
    outfile << "        </DataArray>" << endl;
    outfile << "        <DataArray type=\"Int32\" "
            << "Name=\"offsets\" format=\"ascii\">" << endl;

    outfile << "          ";
    outfile << endl;
    outfile << "        </DataArray>" << endl;
    outfile << "        <DataArray type=\"UInt8\" "
            << "Name=\"types\" format=\"ascii\">" << endl;
    outfile << "          ";
    outfile << endl;
    outfile << "        </DataArray>" << endl;
    outfile << "      </Cells>" << endl;
    outfile << "      <PointData>" << endl;

    outfile << "      </PointData>" << endl;
    outfile << "    </Piece>" << endl;
}
}
}
