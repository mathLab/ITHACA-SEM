///////////////////////////////////////////////////////////////////////////////
//
// File: SessionReader.cpp
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
// Description: Python wrapper for SessionReader.
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Python/NekPyConfig.hpp>

using namespace Nektar::LibUtilities;

/**
 * @brief Thin wrapper around SessionReader to provide a nicer Pythonic
 * interface.
 *
 * This allows us to do, for example
 *
 *     session = SessionReader.CreateInstance(sys.argv)
 *
 * which is more natural in Python.
 */
SessionReaderSharedPtr SessionReader_CreateInstance(py::list &ns)
{
    int i, argc = py::len(ns), bufSize = 0;
    char **argv = new char *[argc+1], *p;

    // Create argc, argv to give to the session reader. Note that this needs to
    // be a contiguous block in memory, otherwise MPI (specifically OpenMPI)
    // will likely segfault.
    for (i = 0; i < argc; ++i)
    {
        std::string tmp = py::extract<std::string>(ns[i]);
        bufSize += tmp.size() + 1;
    }

    std::vector<char> buf(bufSize);
    for (i = 0, p = &buf[0]; i < argc; ++i)
    {
        std::string tmp = py::extract<std::string>(ns[i]);
        std::copy(tmp.begin(), tmp.end(), p);
        p[tmp.size()] = '\0';
        argv[i] = p;
        p += tmp.size()+1;
    }

    // Also make sure we set argv[argc] = NULL otherwise OpenMPI will also
    // segfault.
    argv[argc] = NULL;

    // Create session reader.
    SessionReaderSharedPtr sr = SessionReader::CreateInstance(argc, argv);

    // Clean up.
    delete [] argv;

    return sr;
}

/**
 * @brief SessionReader exports.
 *
 * Currently wrapped functions:
 *   - SessionReader::CreateInstance for creating objects
 *   - SessionReader::GetSessionName to return the session name
 *   - SessionReader::Finalise to deal with finalising things
 */
void export_SessionReader()
{
    py::class_<SessionReader,
           std::shared_ptr<SessionReader>,
           boost::noncopyable>(
               "SessionReader", py::no_init)

        .def("CreateInstance", SessionReader_CreateInstance)
        .staticmethod("CreateInstance")

        .def("GetSessionName", &SessionReader::GetSessionName,
             py::return_value_policy<py::copy_const_reference>())

        .def("Finalise", &SessionReader::Finalise)

        .def("DefinesParameter", &SessionReader::DefinesParameter)
        .def("GetParameter", &SessionReader::GetParameter,
             py::return_value_policy<py::return_by_value>())

        .def("GetVariable", &SessionReader::GetVariable,
             py::return_value_policy<py::copy_const_reference>())

        .def("GetComm", &SessionReader::GetComm)

        ;
}
