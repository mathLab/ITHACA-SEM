////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputTecplot.h
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
//  Description: Tecplot output module
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_OUTPUTTECPLOT
#define FIELDUTILS_OUTPUTTECPLOT

#include "OutputFileBase.h"
#include <tinyxml.h>

namespace Nektar
{
namespace FieldUtils
{

enum TecplotZoneType{
    eOrdered = 0,
    eFELineSeg,
    eFETriangle,
    eFEQuadrilateral,
    eFETetrahedron,
    eFEBrick,
    eFEPolygon,
    eFEPolyhedron
};

/**
 * @brief Tecplot output class.
 */
class OutputTecplot : public OutputFileBase
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<OutputTecplot>::AllocateSharedPtr(f);
    }

    static ModuleKey m_className;
    OutputTecplot(FieldSharedPtr f);
    virtual ~OutputTecplot();


    virtual void Process(po::variables_map &vm);

    virtual std::string GetModuleName()
    {
        return "OutputTecplot";
    }

protected:
    /// Write from pts to output file.
    virtual void OutputFromPts(po::variables_map &vm);

    /// Write from m_exp to output file.
    virtual void OutputFromExp(po::variables_map &vm);

    /// Write from data to output file.
    virtual void OutputFromData(po::variables_map &vm);

    virtual fs::path GetPath(std::string &filename,
                                    po::variables_map &vm);

    virtual fs::path GetFullOutName(std::string &filename,
                                    po::variables_map &vm);

    /// True if writing binary field output
    bool            m_binary;
    /// True if writing a single output file
    bool            m_oneOutputFile;
    /// True if writing header
    bool            m_writeHeader;
    /// Tecplot zone type of output
    TecplotZoneType m_zoneType;
    /// Number of points per block in Tecplot file
    vector<int>     m_numPoints;
    /// Number of blocks in Tecplot file
    int             m_numBlocks;
    /// Coordinate dimension of output
    int             m_coordim;
    /// Total number of connectivity entries
    int             m_totConn;
    /// Connectivty for each block: one per element
    vector<Array<OneD, int> > m_conn;
    /// Each rank's field sizes
    Array<OneD, int> m_rankFieldSizes;
    /// Each rank's connectivity sizes
    Array<OneD, int> m_rankConnSizes;
    /// Field data to output
    Array<OneD, Array<OneD, NekDouble> > m_fields;

    virtual void WriteTecplotHeader(std::ofstream &outfile,
                                    std::vector<std::string> &var);
    virtual void WriteTecplotZone(std::ofstream &outfile);
    virtual void WriteTecplotConnectivity(std::ofstream &outfile);

    void WriteTecplotFile(po::variables_map &vm);

    int GetNumTecplotBlocks();
    void CalculateConnectivity();

};

/**
 * @brief Tecplot output class, specifically for binary field output.
 */
class OutputTecplotBinary : public OutputTecplot
{
public:
    /// Creates an instance of this class
    static std::shared_ptr<Module> create(FieldSharedPtr f)
    {
        return MemoryManager<OutputTecplotBinary>::AllocateSharedPtr(f);
    }

    static ModuleKey m_className;
    OutputTecplotBinary(FieldSharedPtr f) : OutputTecplot(f)
    {
        m_binary = true;
    }

    virtual ~OutputTecplotBinary()
    {
    }

protected:
    void WriteDoubleOrFloat(std::ofstream          &outfile,
                            Array<OneD, NekDouble> &data);
    virtual void WriteTecplotHeader(std::ofstream &outfile,
                                    std::vector<std::string> &var);
    virtual void WriteTecplotZone(std::ofstream &outfile);
    virtual void WriteTecplotConnectivity(std::ofstream &outfile);
};

}
}

#endif
