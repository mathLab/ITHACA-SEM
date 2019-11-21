////////////////////////////////////////////////////////////////////////////////
//
//  File: InputNek5000.cpp
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
//  Description: Reads a Nek5000 checkpoint file.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string.hpp>

#include "InputNek5000.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey InputNek5000::m_className[1] = {
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "fld5000"), InputNek5000::create,
        "Reads Nek5000 field file.")
};

/**
 * @brief Set up InputNek5000 object.
 *
 */
InputNek5000::InputNek5000(FieldSharedPtr f) : InputModule(f)
{
    m_allowedFiles.insert("fld5000");
}

/**
 *
 */
InputNek5000::~InputNek5000()
{
}

/**
 * @brief Process Nek5000 input file.
 *
 * This routine reads a binary-format Nek5000 field file, loads the data into
 * memory and populates the field definitions to match the data format. Nek5000
 * is a classic nodal-Lagrangian spectral element code at a single polynomial
 * order, meaning that the field data are set up according to this structure.
 *
 * This module is adapted from the VisIt visualisation software, which supports
 * a number of Nek5000 inputs.
 */
void InputNek5000::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    ifstream file(m_config["infile"].as<string>().c_str(), ios::binary);

    // Header: 132 bytes for binary.
    vector<char> data(132);
    file.read(&data[0], 132);

    // Check header: should be the four characters #std
    string check(&data[0], 4);
    string header(&data[4], 128);

    ASSERTL0(check == "#std", "Unable to read file");

    // Determine whether we need to byte-swap data: 4-byte float at byte 80
    // should be 6.54321.
    bool byteSwap = false;
    float test;

    file.read((char *)(&test), 4);
    if (test > 6.5 && test < 6.6)
    {
        byteSwap = false;
    }
    else
    {
        swap_endian(test);
        ASSERTL0(test > 6.5 && test < 6.6,
                 "Unable to determine endian-ness of input file");
        byteSwap = true;
    }

    stringstream ss;
    ss.str(header);

    int nBytes, nBlocksXYZ[3], nBlocks, nTotBlocks, dir, nDirs, nCycle, nDim;
    NekDouble time;

    // Read header information (this is written as ASCII)
    string remain;
    ss >> nBytes >> nBlocksXYZ[0] >> nBlocksXYZ[1] >> nBlocksXYZ[2]
       >> nBlocks >> nTotBlocks >> time >> nCycle >> dir >> nDirs >> remain;
    boost::trim(remain);

    nDim = nBlocksXYZ[2] == 1 ? 2 : 3;

    // Print some basic information for input if in verbose mode.
    if (m_f->m_verbose)
    {
        cout << "Found header information:" << endl;
        cout << " -- " << (byteSwap ? "" : "do not ") << "need to swap endian"
             << endl;
        cout << " -- Data byte size       : " << nBytes << endl;
        cout << " -- Number of xyz blocks : " << nBlocksXYZ[0] << "x"
             << nBlocksXYZ[1] << "x" << nBlocksXYZ[2] << endl;
        cout << " -- Blocks in file/total : " << nBlocks << "/" << nTotBlocks
             << endl;
        cout << " -- Simulation time      : " << time << endl;
        cout << " -- Number of cycles     : " << nCycle << endl;
        cout << " -- Number of dirs       : " << dir << "/" << nDirs << endl;
        cout << " -- Remaining header     : " << remain << endl;
    }

    // Major limitation: we don't read out of multiple directories
    ASSERTL0(nDirs == 1, "Number of directories must be one");

    // We also don't read partial files.
    ASSERTL0(nBlocks == nTotBlocks, "Partial field output not supported");

    // We don't support non-double data
    ASSERTL0(nBytes == 8, "Data file must contain double-precision data");

    // Set up a field definition
    LibUtilities::FieldDefinitionsSharedPtr fielddef = MemoryManager<
        LibUtilities::FieldDefinitions>::AllocateSharedPtr();
    fielddef->m_shapeType         = LibUtilities::eHexahedron;
    fielddef->m_numHomogeneousDir = 0;
    fielddef->m_homoStrips        = false;
    fielddef->m_pointsDef         = false;
    fielddef->m_uniOrder          = true;
    fielddef->m_numPointsDef      = false;

    for (int i = 0; i < nDim; ++i)
    {
        fielddef->m_basis.push_back(LibUtilities::eGLL_Lagrange);
        fielddef->m_numModes.push_back(nBlocksXYZ[i]);
    }

    // Read element IDs
    NekUInt32 maxID = 0, minID = numeric_limits<NekUInt32>::max();
    for (NekUInt32 i = 0; i < nBlocks; ++i)
    {
        NekUInt32 blockNum;
        file.read((char *)&blockNum, 4);
        if (byteSwap)
        {
            swap_endian(blockNum);
        }
        fielddef->m_elementIDs.push_back(blockNum-1);

        maxID = maxID > blockNum ? maxID : blockNum;
        minID = minID < blockNum ? minID : blockNum;
    }

    // Determine how many fields we have to extract
    size_t blockSize = nBlocksXYZ[0] * nBlocksXYZ[1] * nBlocksXYZ[2];
    size_t dataSize = blockSize * nBlocks;

    for (string::size_type i = 0; i < remain.size(); ++i)
    {
        switch (remain[i])
        {
            case 'U':
                fielddef->m_fields.push_back("u");
                fielddef->m_fields.push_back("v");
                if (nDim == 3)
                {
                    fielddef->m_fields.push_back("w");
                }
                break;
            case 'P':
                fielddef->m_fields.push_back("p");
                break;
            case 'T':
                fielddef->m_fields.push_back("T");
                break;
            case '1':
            case '2':
            case '3':
                fielddef->m_fields.push_back(string("scalar") + remain[i]);
                break;
            case ' ':
                continue;
            default:
                cerr << "Field contains unknown variable: "
                     << remain[i] << endl;
                abort();
        }
    }

    m_f->m_data.resize(1);
    m_f->m_data[0].resize(fielddef->m_fields.size() * dataSize);

    // Now reprocess since different fields need different logic
    for (size_t i = 0, cnt = 0; i < remain.size(); ++i)
    {
        switch (remain[i])
        {
            case 'U':
            {
                size_t cntVel[3] = {
                    cnt, cnt + dataSize, cnt + 2*dataSize
                };

                for (size_t j = 0; j < nBlocks; ++j)
                {
                    for (size_t k = 0; k < nDim; ++k)
                    {
                        file.read(
                            (char *)&m_f->m_data[0][cntVel[k]],
                            blockSize * sizeof(NekDouble));
                        cntVel[k] += blockSize;
                    }
                }

                cnt += nDim * dataSize;
                break;
            }
            case 'P':
            {
                file.read(
                    (char *)&m_f->m_data[0][cnt],
                    dataSize * sizeof(NekDouble));
                cnt += dataSize;
                break;
            }
            case '1':
            case '2':
            case '3':
            {
                file.read(
                    (char *)&m_f->m_data[0][cnt],
                    dataSize * sizeof(NekDouble));
                cnt += dataSize;
                break;
            }
            case ' ':
                continue;
        }
    }

    m_f->m_fielddef.push_back(fielddef);

    // save field names
    m_f->m_variables = m_f->m_fielddef[0]->m_fields;
}
}
}
