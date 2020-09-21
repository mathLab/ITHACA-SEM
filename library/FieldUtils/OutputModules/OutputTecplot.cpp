///////////////////////////////////////////////////////////////////////
//
//  File: OutputTecplot.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
//  Description: Dat file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <set>
#include <string>
using namespace std;

#include <boost/core/ignore_unused.hpp>

#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/PtsIO.h>

#include "OutputTecplot.h"

namespace Nektar
{
namespace FieldUtils
{

std::string TecplotZoneTypeMap[] = {
    "ORDERED",
    "LINESEG",
    "TRIANGLE",
    "QUADRILATERAL",
    "TETRAHEDRON",
    "BRICK",
    "POLYGON",
    "POLYHEDRON"
};

ModuleKey OutputTecplot::m_className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eOutputModule, "dat"),
        OutputTecplot::create,
        "Writes a Tecplot file.");
ModuleKey OutputTecplotBinary::m_className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eOutputModule, "plt"),
        OutputTecplotBinary::create,
        "Writes a Tecplot file in binary plt format.");

OutputTecplot::OutputTecplot(FieldSharedPtr f) : OutputFileBase(f),
                                                 m_binary(false),
                                                 m_oneOutputFile(false)
{
    m_requireEquiSpaced = true;
    m_config["double"] =
        ConfigOption(true, "0", "Write double-precision format data:"
                                "more accurate but more disk space"
                                " required");
}

OutputTecplot::~OutputTecplot()
{
}

void OutputTecplot::Process(po::variables_map &vm)
{

    if(m_config["writemultiplefiles"].as<bool>())
    {
        m_oneOutputFile = false;
    }
    else
    {
        m_oneOutputFile = (m_f->m_comm->GetSize()> 1);
    }

    OutputFileBase::Process(vm);
}

/**
 * @brief Helper function to write binary data to stream.
 */
template<typename T> void WriteStream(std::ostream &outfile, T data)
{
    T tmp = data;
    outfile.write(reinterpret_cast<char *>(&tmp), sizeof(T));
}

/**
 * @brief Specialisation of WriteStream to support writing std::string.
 *
 * Tecplot binary formats represent all strings by writing out their characters
 * as 32-bit integers, followed by a 32-bit integer null (0) character to denote
 * the end of the string.
 */
template<> void WriteStream(std::ostream &outfile, std::string data)
{
    // Convert string to array of int32_t
    for (std::string::size_type i = 0; i < data.size(); ++i)
    {
        char strChar = data[i];
        NekInt32 strCharInt = strChar;
        WriteStream(outfile, strCharInt);
    }

    // Now dump out zero character to terminate
    WriteStream(outfile, 0);
}

/**
 * @brief Specialisation of WriteStream to support writing Nektar::Array
 * datatype.
 */
template<typename T> void WriteStream(std::ostream &outfile,
                                      Array<OneD, T> data)
{
    outfile.write(reinterpret_cast<char *>(&data[0]),
                  data.size() * sizeof(T));
}

/**
 * @brief Specialisation of WriteStream to support writing std::vector datatype.
 */
template<typename T> void WriteStream(std::ostream  &outfile,
                                      std::vector<T> data)
{
    outfile.write(reinterpret_cast<char *>(&data[0]),
                  data.size() * sizeof(T));
}

void OutputTecplot::OutputFromPts(po::variables_map &vm)
{
    LibUtilities::PtsFieldSharedPtr fPts = m_f->m_fieldPts;

    // do not output if zone is empty
    if (fPts->GetNpoints() == 0)
    {
        return;
    }

    int rank   = m_f->m_comm->GetRank();
    m_numBlocks = 0;

    m_coordim = fPts->GetDim();

    // Grab connectivity information.
    fPts->GetConnectivity(m_conn);

    switch (fPts->GetPtsType())
    {
        case LibUtilities::ePtsFile:
            m_numPoints.resize(1);
            m_numPoints[0] = fPts->GetNpoints();
            m_f->m_comm->AllReduce(m_numPoints[0], LibUtilities::ReduceSum);
            m_zoneType     = eOrdered;
            break;
        case LibUtilities::ePtsLine:
            m_numPoints.resize(1);
            m_numPoints[0] = fPts->GetPointsPerEdge(0);
            m_zoneType     = eOrdered;
            break;
        case LibUtilities::ePtsPlane:
            m_numPoints.resize(2);
            m_numPoints[0] = fPts->GetPointsPerEdge(0);
            m_numPoints[1] = fPts->GetPointsPerEdge(1);
            m_zoneType     = eOrdered;
            break;
        case LibUtilities::ePtsBox:
            m_numPoints.resize(3);
            m_numPoints[0] = fPts->GetPointsPerEdge(0);
            m_numPoints[1] = fPts->GetPointsPerEdge(1);
            m_numPoints[2] = fPts->GetPointsPerEdge(2);
            m_zoneType     = eOrdered;
            break;
        case LibUtilities::ePtsTriBlock:
        {
            m_zoneType = eFETriangle;
            for (int i = 0; i < m_conn.size(); ++i)
            {
                m_numBlocks += m_conn[i].size() / 3;
            }
            break;
        }
        case LibUtilities::ePtsTetBlock:
        {
            m_zoneType = eFETetrahedron;
            for (int i = 0; i < m_conn.size(); ++i)
            {
                m_numBlocks += m_conn[i].size() / 4;
            }
            break;
        }
        default:
            ASSERTL0(false, "This points type is not supported yet.");
    }

    // Get fields and coordinates
    m_fields =
        Array<OneD, Array<OneD, NekDouble> >(m_f->m_variables.size()+m_coordim);

    // We can just grab everything from points. This should be a
    // reference, not a copy.
    fPts->GetPts(m_fields);

    // Only write header if we're root or FE block; binary files always
    // write header
    m_writeHeader =
        (m_zoneType != eOrdered || rank == 0) || m_binary;

    WriteTecplotFile(vm);
}

void OutputTecplot::OutputFromExp(po::variables_map &vm)
{
    m_numBlocks = 0;
    m_writeHeader = true;

    // Calculate number of FE blocks
    m_numBlocks = GetNumTecplotBlocks();

    // Calculate coordinate dimension
    int nBases = m_f->m_exp[0]->GetExp(0)->GetNumBases();

    m_coordim = m_f->m_exp[0]->GetExp(0)->GetCoordim();
    int totpoints = m_f->m_exp[0]->GetTotPoints();

    if (m_f->m_numHomogeneousDir > 0)
    {
        int nPlanes = m_f->m_exp[0]->GetZIDs().size();
        if (nPlanes == 1) // halfMode case
        {
            // do nothing
        }
        else
        {
            // If Fourier points, output extra plane to fill domain
            if (m_f->m_exp[0]->GetExpType() == MultiRegions::e3DH1D)
            {
                nPlanes += 1;
                totpoints += m_f->m_exp[0]->GetPlane(0)->GetTotPoints();
            }
            nBases += m_f->m_numHomogeneousDir;
            m_coordim += m_f->m_numHomogeneousDir;
            NekDouble tmp = m_numBlocks * (nPlanes - 1);
            m_numBlocks   = (int)tmp;
        }
    }

    m_zoneType = (TecplotZoneType)(2*(nBases-1) + 1);

    // Calculate connectivity
    CalculateConnectivity();

    // Set up storage for output fields
    m_fields =
        Array<OneD, Array<OneD, NekDouble> >(m_f->m_variables.size()+m_coordim);

    // Get coordinates
    for (int i = 0; i < m_coordim; ++i)
    {
        m_fields[i] = Array<OneD, NekDouble>(totpoints);
    }

    if (m_coordim == 1)
    {
        m_f->m_exp[0]->GetCoords(m_fields[0]);
    }
    else if (m_coordim == 2)
    {
        m_f->m_exp[0]->GetCoords(m_fields[0], m_fields[1]);
    }
    else
    {
        m_f->m_exp[0]->GetCoords(m_fields[0], m_fields[1], m_fields[2]);
    }

    if (m_f->m_exp[0]->GetExpType() == MultiRegions::e3DH1D)
    {
        // Copy values
        for (int i = 0; i < m_f->m_variables.size(); ++i)
        {
            m_fields[i + m_coordim] = Array<OneD, NekDouble>(totpoints);
            Vmath::Vcopy( m_f->m_exp[0]->GetTotPoints(),
                          m_f->m_exp[i]->UpdatePhys(), 1,
                          m_fields[i + m_coordim], 1);
        }
    }
    else
    {
        // Add references to m_fields
        for (int i = 0; i < m_f->m_variables.size(); ++i)
        {
            m_fields[i + m_coordim] = m_f->m_exp[i]->UpdatePhys();
        }
    }

    // If Fourier, fill extra plane with data
    if (m_f->m_exp[0]->GetExpType() == MultiRegions::e3DH1D)
    {
        int points_on_plane = m_f->m_exp[0]->GetPlane(0)->GetTotPoints();
        const int offset = totpoints-points_on_plane;
        NekDouble z = m_fields[m_coordim-1][totpoints-2*points_on_plane] +
                      (m_fields[m_coordim-1][points_on_plane] - m_fields[m_coordim-1][0]);
        // x and y
        Array<OneD, NekDouble> tmp = m_fields[0] + offset;
        Vmath::Vcopy(points_on_plane, m_fields[0], 1, tmp, 1 );
        tmp = m_fields[1] + offset;
        Vmath::Vcopy(points_on_plane, m_fields[1], 1, tmp, 1 );
        // z coordinate
        tmp = m_fields[2] + offset;
        Vmath::Vcopy(points_on_plane, m_fields[2], 1, tmp, 1 );
        Vmath::Sadd(points_on_plane, z, m_fields[2], 1, tmp, 1 );

        //variables
        for (int i = 0; i < m_f->m_variables.size(); ++i)
        {
            tmp = m_fields[i + m_coordim] + offset;
            Vmath::Vcopy(points_on_plane, m_fields[i + m_coordim], 1, tmp, 1 );
        }
    }

    WriteTecplotFile(vm);
}

void OutputTecplot::OutputFromData(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    NEKERROR(ErrorUtil::efatal,
             "OutputTecplot can't write using only FieldData.");
}

fs::path OutputTecplot::GetPath(std::string &filename,
                                    po::variables_map &vm)
{
    boost::ignore_unused(vm);

    int nprocs = m_f->m_comm->GetSize();
    string       returnstr(filename);

    // Amend for parallel output if required
    if (nprocs != 1 && !m_oneOutputFile)
    {
        int rank   = m_f->m_comm->GetRank();
        int dot       = filename.find_last_of('.');
        string ext    = filename.substr(dot, filename.length() - dot);
        string procId = "_P" + boost::lexical_cast<std::string>(rank);
        string start  = filename.substr(0, dot);
        returnstr     = start + procId + ext;
    }
    return   fs::path(returnstr);
}

fs::path OutputTecplot::GetFullOutName(std::string &filename,
                                        po::variables_map &vm)
{
    return   GetPath(filename, vm);
}

void OutputTecplot::WriteTecplotFile(po::variables_map &vm)
{
    // Variable names
    std::string coordVars[] = { "x", "y", "z" };
    vector<string> variables = m_f->m_variables;
    variables.insert(variables.begin(), coordVars, coordVars + m_coordim);

    int nprocs = m_f->m_comm->GetSize();
    int rank   = m_f->m_comm->GetRank();


    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();
    string outFile  = LibUtilities::PortablePath(GetFullOutName(filename, vm));
    // Open output file
    ofstream outfile;
    if ((m_oneOutputFile && rank == 0) || !m_oneOutputFile)
    {
        outfile.open(outFile.c_str(), m_binary ? ios::binary : ios::out);
    }

    if (m_oneOutputFile)
    {
        // Reduce on number of blocks and number of points.
        m_f->m_comm->AllReduce(m_numBlocks, LibUtilities::ReduceSum);

        // Root process needs to know how much data everyone else has for
        // writing in parallel.
        m_rankFieldSizes       = Array<OneD, int>(nprocs, 0);
        m_rankConnSizes        = Array<OneD, int>(nprocs, 0);
        m_rankFieldSizes[rank] = m_fields[0].size();

        m_totConn = 0;
        for (int i = 0; i < m_conn.size(); ++i)
        {
            m_totConn += m_conn[i].size();
        }

        m_rankConnSizes[rank] = m_totConn;

        m_f->m_comm->AllReduce(m_rankFieldSizes, LibUtilities::ReduceSum);
        m_f->m_comm->AllReduce(m_rankConnSizes, LibUtilities::ReduceSum);
    }

    if (m_writeHeader)
    {
        WriteTecplotHeader(outfile, variables);
    }

    // Write zone data.
    WriteTecplotZone(outfile);

    // If we're a FE block format, write connectivity (m_conn will be empty for
    // point data).
    WriteTecplotConnectivity(outfile);

    if ((m_oneOutputFile && rank == 0) || !m_oneOutputFile)
    {
        cout << "Written file: " << GetFullOutName(filename,vm) << endl;
    }
}

/**
 * @brief Write Tecplot files header
 *
 * @param   outfile   Output file name
 * @param   var       Variables names
 */
void OutputTecplot::WriteTecplotHeader(std::ofstream &outfile,
                                       std::vector<std::string> &var)
{
    outfile << "Variables = " << var[0];

    for (int i = 1; i < var.size(); ++i)
    {
        outfile << ", " << var[i];
    }

    outfile << std::endl << std::endl;
}

/**
 * @brief Write Tecplot files header in binary format
 *
 * @param   outfile   Output file name
 * @param   var       Variables names
 */
void OutputTecplotBinary::WriteTecplotHeader(std::ofstream &outfile,
                                             std::vector<std::string> &var)
{
    if (m_oneOutputFile && m_f->m_comm->GetRank() > 0)
    {
        return;
    }

    // Version number
    outfile << "#!TDV112";

    // Int value of 1 for endian check
    WriteStream(outfile, 1);

    // We'll probably write a full solution field
    WriteStream(outfile, 0);

    // Title
    std::string title = "";
    WriteStream(outfile, title);

    // Number of variables
    WriteStream(outfile, (int)var.size());

    for (int i = 0; i < var.size(); ++i)
    {
        WriteStream(outfile, var[i]);
    }
}


/**
 * Write Tecplot zone output in ASCII
 *
 * @param   outfile    Output file name.
 * @param   expansion  Expansion that is considered
 */
void OutputTecplot::WriteTecplotZone(std::ofstream &outfile)
{
    bool useDoubles = m_config["double"].as<bool>();

    if (useDoubles)
    {
        int precision = std::numeric_limits<double>::max_digits10;
        outfile << std::setprecision(precision);

    }

    // Write either points or finite element block
    if (m_zoneType != eOrdered)
    {
        if ((m_oneOutputFile && m_f->m_comm->GetRank() == 0) || !m_oneOutputFile)
        {
            // Number of points in zone
            int nPoints = m_oneOutputFile ?
                Vmath::Vsum(m_f->m_comm->GetSize(), m_rankFieldSizes, 1) :
                m_fields[0].size();

            outfile << "Zone, N=" << nPoints << ", E="
                    << m_numBlocks << ", F=FEBlock, ET="
                    << TecplotZoneTypeMap[m_zoneType] << std::endl;
        }


        if (m_oneOutputFile && m_f->m_comm->GetRank() == 0)
        {
            for (int j = 0; j < m_fields.size(); ++j)
            {
                for (int i = 0; i < m_fields[j].size(); ++i)
                {
                    if ((!(i % 1000)) && i)
                    {
                        outfile << std::endl;
                    }
                    outfile << m_fields[j][i] << " ";
                }

                for (int n = 1; n < m_f->m_comm->GetSize(); ++n)
                {
                    if(m_rankFieldSizes[n])
                    {
                        Array<OneD, NekDouble> tmp(m_rankFieldSizes[n]);
                        m_f->m_comm->Recv(n, tmp);

                        for (int i = 0; i < m_rankFieldSizes[n]; ++i)
                        {
                            if ((!(i % 1000)) && i)
                            {
                                outfile << std::endl;
                            }
                            outfile << tmp[i] << " ";
                        }
                    }
                }
                outfile << std::endl;
            }
        }
        else if (m_oneOutputFile && m_f->m_comm->GetRank() > 0)
        {
            if(m_fields[0].size())
            {
                for (int i = 0; i < m_fields.size(); ++i)
                {
                    m_f->m_comm->Send(0, m_fields[i]);
                }
            }
        }
        else
        {
            // Write out coordinates and field data: ordered by field
            // and then its data.
            for (int j = 0; j < m_fields.size(); ++j)
            {
                for (int i = 0; i < m_fields[j].size(); ++i)
                {
                    if ((!(i % 1000)) && i)
                    {
                        outfile << std::endl;
                    }
                    outfile << m_fields[j][i] << " ";
                }
                outfile << std::endl;
            }
        }
    }
    else
    {
        if((m_oneOutputFile && m_f->m_comm->GetRank() == 0) || !m_oneOutputFile)
        {
            std::string dirs[] = { "I", "J", "K" };
            outfile << "Zone";
            for (int i = 0; i < m_numPoints.size(); ++i)
            {
                outfile << ", " << dirs[i] << "=" << m_numPoints[i];
            }
            outfile << ", F=POINT" << std::endl;
        }

        if (m_oneOutputFile && m_f->m_comm->GetRank() == 0)
        {
            Array<OneD, NekDouble> tmp(m_fields.size());
            for (int i = 0; i < m_fields[0].size(); ++i)
            {
                for (int j = 0; j < m_fields.size(); ++j)
                {
                    outfile << setw(12) << m_fields[j][i] << " ";
                }
                outfile << std::endl;
            }

            for (int n = 1; n < m_f->m_comm->GetSize(); ++n)
            {
                for (int i = 0; i < m_rankFieldSizes[n]; ++i)
                {
                    m_f->m_comm->Recv(n, tmp);
                    for (int j = 0; j < m_fields.size(); ++j)
                    {
                        outfile << setw(12) << tmp[j] << " ";
                    }
                    outfile << std::endl;
                }
            }
        }
        else if (m_oneOutputFile && m_f->m_comm->GetRank() > 0)
        {
            Array<OneD, NekDouble> tmp(m_fields.size());
            for (int i = 0; i < m_fields[0].size(); ++i)
            {
                for (int j = 0; j < m_fields.size(); ++j)
                {
                    tmp[j] = m_fields[j][i];
                }
                m_f->m_comm->Send(0, tmp);
            }
        }
        else
        {
            // Write out coordinates and field data: ordered by each
            // point then each field.
            for (int i = 0; i < m_fields[0].size(); ++i)
            {
                for (int j = 0; j < m_fields.size(); ++j)
                {
                    outfile << setw(12) << m_fields[j][i] << " ";
                }
                outfile << std::endl;
            }
        }
    }
}

/**
 * @brief Write either double-precision or single-precision output of field
 * data.
 *
 * @param   outfile    Output file name.
 * @param   expansion  Expansion that is considered
 */
void OutputTecplotBinary::WriteDoubleOrFloat(std::ofstream          &outfile,
                                             Array<OneD, NekDouble> &data)
{
    // Data format: either double or single depending on user options
    bool useDoubles = m_config["double"].as<bool>();

    if (useDoubles)
    {
        // For doubles, we can just write data.
        WriteStream(outfile, data);
    }
    else
    {
        // For single precision, needs typecast first.
        int nPts = data.size();
        vector<float> tmp(data.size());
        std::copy(&data[0], &data[0] + nPts, &tmp[0]);
        WriteStream(outfile, tmp);
    }
}

/**
 * Write Tecplot zone output in binary
 *
 * @param   outfile    Output file name.
 * @param   expansion  Expansion that is considered
 */
void OutputTecplotBinary::WriteTecplotZone(std::ofstream &outfile)
{
    Array<OneD, NekDouble> fieldMin(m_fields.size());
    Array<OneD, NekDouble> fieldMax(m_fields.size());

    // Data format: either double or single depending on user options
    bool useDoubles = m_config["double"].as<bool>();

    if ((m_oneOutputFile && m_f->m_comm->GetRank() == 0) || !m_oneOutputFile)
    {
        // Don't bother naming zone
        WriteStream(outfile, 299.0f); // Zone marker

        // Write same name as preplot
        int rank   = m_f->m_comm->GetRank();
        string zonename = "ZONE " + boost::lexical_cast<string>(rank);
        WriteStream(outfile, zonename);

        WriteStream(outfile, -1); // No parent zone
        WriteStream(outfile, -1); // No strand ID
        WriteStream(outfile, 0.0); // Solution time
        WriteStream(outfile, -1); // Unused, set to -1

        // Zone type: 1 = lineseg, 3 = quad, 5 = brick
        WriteStream(outfile, (int)m_zoneType);

        WriteStream(outfile, 0); // Data at nodes
        WriteStream(outfile, 0); // No 1-1 face neighbours
        WriteStream(outfile, 0); // No user-defined connections

        if (m_zoneType == eOrdered)
        {
            for (int i = 0; i < m_numPoints.size(); ++i)
            {
                WriteStream(outfile, m_numPoints[i]);
            }

            for (int i = m_numPoints.size(); i < 3; ++i)
            {
                WriteStream(outfile, 0);
            }
        }
        else
        {
            // Number of points in zone
            int nPoints = m_oneOutputFile ?
                Vmath::Vsum(m_f->m_comm->GetSize(), m_rankFieldSizes, 1) :
                m_fields[0].size();

            WriteStream(outfile, nPoints); // Total number of points
            WriteStream(outfile, m_numBlocks); // Number of blocks
            WriteStream(outfile, 0); // Unused
            WriteStream(outfile, 0); // Unused
            WriteStream(outfile, 0); // Unused
        }

        WriteStream(outfile, 0); // No auxiliary data names

        // Finalise header
        WriteStream(outfile, 357.0f);

        // Now start to write data section so that we can dump geometry
        // information

        // Data marker
        WriteStream(outfile, 299.0f);

        for (int j = 0; j < m_fields.size(); ++j)
        {
            WriteStream(outfile, useDoubles ? 2 : 1);
        }

        // No passive variables or variable sharing, no zone connectivity
        // sharing (we only dump one zone)
        WriteStream(outfile, 0);
        WriteStream(outfile, 0);
        WriteStream(outfile, -1);
    }

    for (int i = 0; i < m_fields.size(); ++i)
    {
        fieldMin[i] = Vmath::Vmin(m_fields[i].size(), m_fields[i], 1);
        fieldMax[i] = Vmath::Vmax(m_fields[i].size(), m_fields[i], 1);
    }

    m_f->m_comm->AllReduce(fieldMin, LibUtilities::ReduceMin);
    m_f->m_comm->AllReduce(fieldMax, LibUtilities::ReduceMax);

    // Write out min/max of field data
    if ((m_oneOutputFile && m_f->m_comm->GetRank() == 0) || !m_oneOutputFile)
    {
        for (int i = 0; i < m_fields.size(); ++i)
        {
            WriteStream(outfile, fieldMin[i]);
            WriteStream(outfile, fieldMax[i]);
        }
    }

    if (m_oneOutputFile && m_f->m_comm->GetRank() == 0)
    {
        for (int i = 0; i < m_fields.size(); ++i)
        {
            WriteDoubleOrFloat(outfile, m_fields[i]);

            for (int n = 1; n < m_f->m_comm->GetSize(); ++n)
            {
                Array<OneD, NekDouble> tmp(m_rankFieldSizes[n]);
                m_f->m_comm->Recv(n, tmp);
                WriteDoubleOrFloat(outfile, tmp);
            }
        }
    }
    else if (m_oneOutputFile && m_f->m_comm->GetRank() > 0)
    {
        for (int i = 0; i < m_fields.size(); ++i)
        {
            m_f->m_comm->Send(0, m_fields[i]);
        }
    }
    else
    {
        for (int i = 0; i < m_fields.size(); ++i)
        {
            WriteDoubleOrFloat(outfile, m_fields[i]);
        }
    }
}

/**
 * @brief Write Tecplot connectivity information (ASCII)
 *
 * @param   outfile    Output file
 */
void OutputTecplot::WriteTecplotConnectivity(std::ofstream &outfile)
{
    // Ordered data have no connectivity information.
    if (m_zoneType == eOrdered)
    {
        return;
    }

    if (m_oneOutputFile && m_f->m_comm->GetRank() > 0)
    {
        // Need to amalgamate connectivity information
        if (m_totConn)
        {
            Array<OneD, int> conn(m_totConn);
            for (int i = 0, cnt = 0; i < m_conn.size(); ++i)
            {
                if(m_conn[i].size())
                {
                    Vmath::Vcopy(m_conn[i].size(), &m_conn[i][0], 1,
                                 &conn[cnt], 1);
                    cnt += m_conn[i].size();
                }
            }
            m_f->m_comm->Send(0, conn);
        }
    }
    else
    {
        int cnt = 1;
        for (int i = 0; i < m_conn.size(); ++i)
        {
            const int nConn = m_conn[i].size();
            for (int j = 0; j < nConn; ++j,++cnt)
            {
                outfile << m_conn[i][j] + 1 << " ";
                if (!(cnt % 1000))
                {
                    outfile << std::endl;
                }
            }
        }
        outfile << endl;

        if (m_oneOutputFile && m_f->m_comm->GetRank() == 0)
        {
            int offset = m_rankFieldSizes[0];

            for (int n = 1; n < m_f->m_comm->GetSize(); ++n)
            {
                if(m_rankConnSizes[n])
                {
                    Array<OneD, int> conn(m_rankConnSizes[n]);
                    m_f->m_comm->Recv(n, conn);
                    for (int j = 0; j < conn.size(); ++j)
                    {
                        outfile << conn[j] + offset + 1 << " ";
                        if ((!(j % 1000)) && j)
                        {
                            outfile << std::endl;
                        }
                    }
                }
                offset += m_rankFieldSizes[n];
            }
        }
    }
}

void OutputTecplotBinary::WriteTecplotConnectivity(std::ofstream &outfile)
{
    if (m_oneOutputFile && m_f->m_comm->GetRank() > 0)
    {
        // Need to amalgamate connectivity information
        Array<OneD, int> conn(m_totConn);
        for (int i = 0, cnt = 0; i < m_conn.size(); ++i)
        {
            Vmath::Vcopy(m_conn[i].size(), &m_conn[i][0], 1,
                         &conn[cnt], 1);
            cnt += m_conn[i].size();
        }
        m_f->m_comm->Send(0, conn);
    }
    else
    {
        for (int i = 0; i < m_conn.size(); ++i)
        {
            WriteStream(outfile, m_conn[i]);
        }

        if (m_oneOutputFile && m_f->m_comm->GetRank() == 0)
        {
            int offset = m_rankFieldSizes[0];

            for (int n = 1; n < m_f->m_comm->GetSize(); ++n)
            {
                Array<OneD, int> conn(m_rankConnSizes[n]);
                m_f->m_comm->Recv(n, conn);

                for (int j = 0; j < conn.size(); ++j)
                {
                    conn[j] += offset;
                }

                WriteStream(outfile, conn);
                offset += m_rankFieldSizes[n];
            }
        }
    }
}

/**
 * @brief Calculate number of Tecplot blocks.
 *
 * @param   outfile    Output file
 */
int OutputTecplot::GetNumTecplotBlocks()
{
    int returnval = 0;

    if (m_f->m_exp[0]->GetExp(0)->GetNumBases() == 1)
    {
        for (int i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
        {
            returnval += (m_f->m_exp[0]->GetExp(i)->GetNumPoints(0) - 1);
        }
    }
    else if (m_f->m_exp[0]->GetExp(0)->GetNumBases() == 2)
    {
        for (int i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
        {
            returnval += (m_f->m_exp[0]->GetExp(i)->GetNumPoints(0) - 1) *
                         (m_f->m_exp[0]->GetExp(i)->GetNumPoints(1) - 1);
        }
    }
    else
    {
        for (int i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
        {
            returnval += (m_f->m_exp[0]->GetExp(i)->GetNumPoints(0) - 1) *
                         (m_f->m_exp[0]->GetExp(i)->GetNumPoints(1) - 1) *
                         (m_f->m_exp[0]->GetExp(i)->GetNumPoints(2) - 1);
        }
    }

    return returnval;
}

/**
 * @brief Calculate connectivity information for each expansion dimension.
 *
 * @param   outfile    Output file
 */
void OutputTecplot::CalculateConnectivity()
{
    int i, j, k, l;
    int nbase = m_f->m_exp[0]->GetExp(0)->GetNumBases();
    int cnt   = 0;

    m_conn.resize(m_f->m_exp[0]->GetNumElmts());

    for (i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
    {
        cnt = m_f->m_exp[0]->GetPhys_Offset(i);

        if (nbase == 1)
        {
            int cnt2    = 0;
            int np0     = m_f->m_exp[0]->GetExp(i)->GetNumPoints(0);
            int nPlanes = 1;

            if (m_f->m_exp[0]->GetExpType() == MultiRegions::e2DH1D)
            {
                nPlanes = m_f->m_exp[0]->GetZIDs().size();

                if (nPlanes > 1)
                {
                    int totPoints = m_f->m_exp[0]->GetPlane(0)->GetTotPoints();

                    Array<OneD, int> conn(4 * (np0 - 1) * (nPlanes - 1));
                    for (int n = 1; n < nPlanes; ++n)
                    {
                        for (k = 1; k < np0; ++k)
                        {
                            conn[cnt2++] = cnt + (n - 1) * totPoints + k;
                            conn[cnt2++] = cnt + (n - 1) * totPoints + k - 1;
                            conn[cnt2++] = cnt +  n      * totPoints + k - 1;
                            conn[cnt2++] = cnt +  n      * totPoints + k;
                        }
                    }
                    m_conn[i] = conn;
                }
            }

            if (nPlanes == 1)
            {
                Array<OneD, int> conn(2 * (np0 - 1));

                for (k = 1; k < np0; ++k)
                {
                    conn[cnt2++] = cnt + k;
                    conn[cnt2++] = cnt + k - 1;
                }

                m_conn[i] = conn;
            }
        }
        else if (nbase == 2)
        {
            int np0       = m_f->m_exp[0]->GetExp(i)->GetNumPoints(0);
            int np1       = m_f->m_exp[0]->GetExp(i)->GetNumPoints(1);
            int totPoints = m_f->m_exp[0]->GetTotPoints();
            int nPlanes   = 1;
            int cnt2      = 0;

            if (m_f->m_exp[0]->GetExpType() == MultiRegions::e3DH1D)
            {
                nPlanes = m_f->m_exp[0]->GetZIDs().size();

                // default to 2D case for HalfMode when nPlanes = 1
                if (nPlanes > 1)
                {
                    // If Fourier points, output extra plane to fill domain
                    nPlanes += 1;
                    totPoints = m_f->m_exp[0]->GetPlane(0)->GetTotPoints();

                    Array<OneD, int>
                            conn(8 * (np1 - 1) * (np0 - 1) * (nPlanes - 1));

                    for (int n = 1; n < nPlanes; ++n)
                    {
                        for (j = 1; j < np1; ++j)
                        {
                            for (k = 1; k < np0; ++k)
                            {
                                conn[cnt2++] = cnt + (n - 1) * totPoints +
                                    (j - 1) * np0 + k - 1;
                                conn[cnt2++] = cnt + (n - 1) * totPoints +
                                    (j - 1) * np0 + k;
                                conn[cnt2++] = cnt + (n - 1) * totPoints +
                                    j * np0 + k;
                                conn[cnt2++] = cnt + (n - 1) * totPoints +
                                    j * np0 + k - 1;
                                conn[cnt2++] = cnt + n * totPoints +
                                    (j - 1) * np0 + k - 1;
                                conn[cnt2++] = cnt + n * totPoints +
                                    (j - 1) * np0 + k;
                                conn[cnt2++] = cnt + n * totPoints +
                                    j * np0 + k;
                                conn[cnt2++] = cnt + n * totPoints +
                                    j * np0 + k - 1;
                            }
                        }
                    }
                    m_conn[i] = conn;
                }
            }

            if (nPlanes == 1)
            {
                Array<OneD, int> conn(4 * (np0 - 1) * (np1 - 1));
                for (j = 1; j < np1; ++j)
                {
                    for (k = 1; k < np0; ++k)
                    {
                        conn[cnt2++] = cnt + (j - 1) * np0 + k - 1;
                        conn[cnt2++] = cnt + (j - 1) * np0 + k;
                        conn[cnt2++] = cnt + j * np0 + k;
                        conn[cnt2++] = cnt + j * np0 + k - 1;
                    }
                }
                m_conn[i] = conn;
            }
        }
        else if (nbase == 3)
        {
            int np0  = m_f->m_exp[0]->GetExp(i)->GetNumPoints(0);
            int np1  = m_f->m_exp[0]->GetExp(i)->GetNumPoints(1);
            int np2  = m_f->m_exp[0]->GetExp(i)->GetNumPoints(2);
            int cnt2 = 0;

            Array<OneD, int> conn(8 * (np0 - 1) * (np1 - 1) * (np2 - 1));

            for (j = 1; j < np2; ++j)
            {
                for (k = 1; k < np1; ++k)
                {
                    for (l = 1; l < np0; ++l)
                    {
                        conn[cnt2++] =
                            cnt + (j - 1) * np0 * np1 + (k - 1) * np0 + l - 1;
                        conn[cnt2++] =
                            cnt + (j - 1) * np0 * np1 + (k - 1) * np0 + l;
                        conn[cnt2++] =
                            cnt + (j - 1) * np0 * np1 + k * np0 + l;
                        conn[cnt2++] =
                            cnt + (j - 1) * np0 * np1 + k * np0 + l - 1;
                        conn[cnt2++] =
                            cnt + j * np0 * np1 + (k - 1) * np0 + l - 1;
                        conn[cnt2++] =
                            cnt + j * np0 * np1 + (k - 1) * np0 + l;
                        conn[cnt2++] =
                            cnt + j * np0 * np1 + k * np0 + l;
                        conn[cnt2++] =
                            cnt + j * np0 * np1 + k * np0 + l - 1;
                    }
                }
            }

            m_conn[i] = conn;
        }
        else
        {
            ASSERTL0(false, "Not set up for this dimension");
        }

    }
}
}
}
