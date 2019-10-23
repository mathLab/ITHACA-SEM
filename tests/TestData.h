///////////////////////////////////////////////////////////////////////////////
//
// File: Tester.cpp
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
// Description: Encapsulation of test XML file.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_TESTER_TESTDATA
#define NEKTAR_TESTER_TESTDATA

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>

#include <tinyxml.h>

namespace fs = boost::filesystem;
namespace po = boost::program_options;

namespace Nektar
{
    struct DependentFile
    {
        std::string m_description;
        std::string m_filename;
    };

    struct Command
    {
        fs::path     m_executable;
        std::string  m_parameters;
        unsigned int m_processes;
        bool         m_pythonTest;
    };

    class TestData
    {
    public:
        TestData(const fs::path& pFilename, po::variables_map& pVm);
        TestData(const TestData& pSrc);

        const std::string& GetDescription() const;
        const Command &GetCommand(unsigned int pId) const;
        unsigned int GetNumCommands() const;

        std::string GetMetricType(unsigned int pId) const;
        unsigned int GetNumMetrics() const;
        TiXmlElement* GetMetric(unsigned int pId);
        unsigned int GetMetricId(unsigned int pId);

        DependentFile GetDependentFile(unsigned int pId) const;
        unsigned int GetNumDependentFiles() const;

        void SaveFile();

    private:
        po::variables_map               m_cmdoptions;
        std::string                     m_description;
        std::vector<Command>            m_commands;
        TiXmlDocument*                  m_doc;
        std::vector<TiXmlElement*>      m_metrics;
        std::vector<DependentFile>      m_files;

        void Parse(TiXmlDocument* pDoc);
        Command ParseCommand(TiXmlElement *pElmt) const;
    };
}

#endif
