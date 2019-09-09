////////////////////////////////////////////////////////////////////////////////
//
//  File: Module.cpp
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
//  Description: Abstract input/output modules.
//
////////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include "Module.h"

using namespace std;

namespace Nektar
{
namespace FieldUtils
{
/**
 * Returns an instance of the module factory, held as a singleton.
 */
ModuleFactory &GetModuleFactory()
{
    static ModuleFactory instance;
    return instance;
}

/**
 * Prints a given module key to a stream.
 */
std::ostream &operator<<(std::ostream &os, const ModuleKey &rhs)
{
    return os << ModuleTypeMap[rhs.first] << ": " << rhs.second;
}

InputModule::InputModule(FieldSharedPtr m) : Module(m)
{
    m_config["infile"] = ConfigOption(false, "", "Input filename.");
}

OutputModule::OutputModule(FieldSharedPtr m) : Module(m)
{
    m_config["outfile"] = ConfigOption(false, "", "Output filename.");
}

void InputModule::AddFile(string fileType, string fileName)
{
    // Check to see if this file type is allowed
    if (m_allowedFiles.count(fileType) == 0)
    {
        cerr << "File type " << fileType << " not supported for this "
             << "module." << endl;
    }

    m_f->m_inputfiles[fileType].push_back(fileName);
}
/**
 * @brief Open a file for output.
 */
void OutputModule::OpenStream()
{
    string fname = m_config["outfile"].as<string>();
    m_fldFile.open(fname.c_str());
    if (!m_fldFile.good())
    {
        cerr << "Error opening file: " << fname << endl;
        abort();
    }
}

/**
 * @brief Register a configuration option with a module.
 */
void Module::RegisterConfig(string key, string val)
{
    auto it = m_config.find(key);
    if (it == m_config.end())
    {
        cerr << "WARNING: Unrecognised config option " << key
             << ", proceeding anyway." << endl;
    }

    it->second.m_beenSet = true;

    if (it->second.m_isBool && val=="")
    {
        it->second.m_value = "1";
    }
    else
    {
        it->second.m_value = val;
    }
}

/**
 * @brief Print out all configuration options for a module.
 */
void Module::PrintConfig()
{
    if (m_config.size() == 0)
    {
        cerr << "No configuration options for this module." << endl;
        return;
    }

    for (auto &it : m_config)
    {
        cerr << setw(10) << it.first << ": " << it.second.m_desc << endl;
    }
}

/**
 * @brief Sets default configuration options for those which have not
 * been set.
 */
void Module::SetDefaults()
{
    for (auto &it : m_config)
    {
        if (!it.second.m_beenSet)
        {
            it.second.m_value = it.second.m_defValue;
        }
    }
}

/**
 * @brief Tries to guess the format of the input file.
 */
string InputModule::GuessFormat(string filename)
{
    // Read first 64 bytes of data, assuming input is this long.
    ifstream inFile(filename.c_str(), ios::binary);
    vector<char> data(64, 0);
    inFile.read(&data[0], 64);

    string check(&data[0], 64);

    // Nek5000 format: first four characters are: #std
    if (check.compare(0, 4, "#std") == 0)
    {
        inFile.close();
        return "fld5000";
    }

    // Semtex format: first line should contain the string "Session" at
    // character 27.
    if (check.compare(26, 7, "Session") == 0)
    {
        inFile.close();
        return "fldsem";
    }

    // Otherwise don't really know -- try to guess from file extension.
    inFile.close();
    return "";
}

/**
 * @brief Print a brief summary of information.
 */
void InputModule::PrintSummary()
{
    cout << "Field size = " << m_f->m_data[0].size() * sizeof(NekDouble)
         << endl;
}
}
}
