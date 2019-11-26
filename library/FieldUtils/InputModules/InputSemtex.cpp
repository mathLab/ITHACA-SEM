////////////////////////////////////////////////////////////////////////////////
//
//  File: InputSemtex.cpp
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
//  Description: Reads a Semtex checkpoint file.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

#include <boost/core/ignore_unused.hpp>
#include <boost/algorithm/string.hpp>

#include <LibUtilities/BasicUtils/CompressData.h>

#include "InputSemtex.h"

namespace Nektar
{
namespace FieldUtils
{

ModuleKey InputSemtex::m_className[1] = {
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "fldsem"), InputSemtex::create,
        "Reads Semtex field file.")
};

/**
 * @brief Set up InputSemtex object.
 *
 */
InputSemtex::InputSemtex(FieldSharedPtr f) : InputModule(f)
{
    m_allowedFiles.insert("fldsem");
}

/**
 *
 */
InputSemtex::~InputSemtex()
{
}

/**
 * @brief Process Semtex input file.
 *
 * This routine reads a binary-format Semtex field file, loads the data into
 * memory and populates the field definitions to match the data format. Semtex
 * is a classic nodal-Lagrangian spectral element code at a single polynomial
 * order, meaning that the field data are set up according to this structure.
 */
void InputSemtex::Process(po::variables_map &vm)
{
    boost::ignore_unused(vm);

    // Variables to be read from session file
    string sessionName, date, fields, endian;
    int nr, ns, nz, nelmt, step;
    NekDouble time, dt, kinvis, beta;

    ifstream file(m_config["infile"].as<string>().c_str(), ios::binary);

    // -- Read header information.
    char buf[25];
    string line;

    // Session name
    file.read(buf, 25);
    sessionName = string(buf, 25);
    boost::trim(sessionName);
    getline(file, line);
    m_f->m_fieldMetaDataMap["SessionName0"] = sessionName;

    // Date
    file.read(buf, 25);
    date = string(buf, 25);
    boost::trim(date);
    getline(file, line);

    // nP, nZ, nElmt
    file >> nr >> ns >> nz >> nelmt;
    getline(file, line);

    // Step
    file >> step;
    getline(file, line);

    // Time
    file >> time;
    getline(file, line);
    m_f->m_fieldMetaDataMap["Time"] = boost::lexical_cast<string>(time);

    // Timestep
    file >> dt;
    getline(file, line);
    m_f->m_fieldMetaDataMap["TimeStep"] = boost::lexical_cast<string>(dt);

    // Viscosity
    file >> kinvis;
    getline(file, line);
    m_f->m_fieldMetaDataMap["Kinvis"] = boost::lexical_cast<string>(kinvis);

    // Beta
    file >> beta;
    getline(file, line);

    // Fields
    file.read(buf, 25);
    fields = string(buf, 25);
    boost::trim(fields);
    getline(file, line);

    // Endian-ness
    LibUtilities::EndianType systemEndian =  LibUtilities::Endianness();
    std::string endianSearch;
    if (systemEndian == LibUtilities::eEndianBig)
    {
        endianSearch = "big";
    }
    else if (systemEndian == LibUtilities::eEndianLittle)
    {
        endianSearch = "little";
    }
    else
    {
        ASSERTL0(false, "Only little- or big-endian systems are supported");
    }

    file.read(buf, 25);
    endian = string(buf, 25);
    bool byteSwap = endian.find(endianSearch) == string::npos;
    getline(file, line);

    // Print some basic information for input if in verbose mode.
    if (m_f->m_verbose)
    {
        cout << "Found header information:" << endl;
        cout << " -- From session         : " << sessionName << endl;
        cout << " -- File generated       : " << date        << endl;
        cout << " -- Polynomial order     : " << nr-1        << endl;
        cout << " -- Number of planes     : " << nz          << endl;
        cout << " -- Number of elements   : " << nelmt       << endl;
        cout << " -- Simulation time      : " << time        << endl;
        cout << " -- Timestep             : " << dt          << endl;
        cout << " -- Viscosity            : " << kinvis      << endl;
        cout << " -- Fields               : " << fields
             << " (" << fields.size() << " total)" << endl;

        if (nz > 1)
        {
            cout << " -- Homogeneous length   : " << 2*M_PI/beta << endl;
        }

        cout << " -- " << (byteSwap ? "" : "do not ") << "need to swap endian"
             << endl;
    }

    ASSERTL0(nr == ns, "Semtex reader assumes values of nr and ns are equal");

    // Set up a field definition
    LibUtilities::FieldDefinitionsSharedPtr fielddef = MemoryManager<
        LibUtilities::FieldDefinitions>::AllocateSharedPtr();
    fielddef->m_shapeType         = LibUtilities::eQuadrilateral;
    fielddef->m_homoStrips        = false;
    fielddef->m_pointsDef         = false;
    fielddef->m_uniOrder          = true;
    fielddef->m_numPointsDef      = false;

    // Set up basis
    fielddef->m_basis.push_back(LibUtilities::eGLL_Lagrange);
    fielddef->m_basis.push_back(LibUtilities::eGLL_Lagrange);
    fielddef->m_numModes.push_back(nr);
    fielddef->m_numModes.push_back(nr);

    // Set up elements
    fielddef->m_elementIDs.resize(nelmt);
    for (int i = 0; i < nelmt; ++i)
    {
        fielddef->m_elementIDs[i] = i;
    }

    // Deal with homogeneous direction.
    if (nz > 1)
    {
        fielddef->m_numHomogeneousDir = 1;
        fielddef->m_homogeneousLengths.push_back(2 * M_PI / beta);
        fielddef->m_numModes.push_back(nz);
        fielddef->m_basis.push_back(LibUtilities::eFourier);

        for (int i = 0; i < nz; ++i)
        {
            fielddef->m_homogeneousZIDs.push_back(i);
        }
    }
    else
    {
        fielddef->m_numHomogeneousDir = 0;
    }

    for (string::size_type i = 0; i < fields.size(); ++i)
    {
        fielddef->m_fields.push_back(string(&fields[i], 1));
    }

    // Size of data to read.
    size_t elmtSize  = nr * ns;
    size_t planeSize = elmtSize * nelmt;
    size_t fieldSize = planeSize * nz;
    size_t dataSize = fieldSize * fields.size();

    // Allocate our storage.
    m_f->m_data.resize(1);
    m_f->m_data[0].resize(dataSize);

    // Temporary storage for one plane of data.
    vector<NekDouble> tmp(planeSize);
    size_t offset = nz * nr * ns;

    // Now reorder data; Semtex ordering traverses memory fastest over planes,
    // whereas Nektar++ expects it over elements
    for (int i = 0; i < fields.size(); ++i)
    {
        NekDouble *data = &m_f->m_data[0][i * fieldSize];
        for (int j = 0; j < nz; ++j)
        {
            size_t elSizeJ = j * elmtSize;
            file.read((char *)&tmp[0], planeSize * sizeof(NekDouble));

            if (byteSwap)
            {
                swap_endian(tmp);
            }

            double* x = &tmp[0];
            for (int k = 0; k < nelmt; ++k)
            {
                std::copy(x, x + elmtSize, data + k * offset + elSizeJ);
                x += elmtSize;

            }
        }
    }

    m_f->m_fielddef.push_back(fielddef);

    // save field names
    m_f->m_variables = m_f->m_fielddef[0]->m_fields;
}

}
}
