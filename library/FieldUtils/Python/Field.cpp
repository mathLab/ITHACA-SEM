///////////////////////////////////////////////////////////////////////////////
//
// File: Field.cpp
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
// Description: Python wrapper for Field class.
//
///////////////////////////////////////////////////////////////////////////////

#include <FieldUtils/Field.hpp>
#include <FieldUtils/FieldConvertComm.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>

using namespace Nektar;
using namespace Nektar::FieldUtils;

// Converts Python's sys.argv into C++'s argv.
char **ConvertCommandLine(py::list &py_argv)
{
    int i, argc = py::len(py_argv), bufSize = 0;
    char **argv = new char *[argc + 1], *p;

    for (i = 0; i < argc; ++i)
    {
        std::string tmp = py::extract<std::string>(py_argv[i]);
        bufSize += tmp.size() + 1;
    }

    std::vector<char> buf(bufSize);
    for (i = 0, p = &buf[0]; i < argc; ++i)
    {
        std::string tmp = py::extract<std::string>(py_argv[i]);
        std::copy(tmp.begin(), tmp.end(), p);
        p[tmp.size()] = '\0';
        argv[i]       = p;
        p += tmp.size() + 1;
    }

    argv[argc] = NULL;
    return argv;
}

// Called at the start of every loop over nparts.
// Args: FieldSharedPtr "f",
// Python's sys.argv stored as a Python list called "py_argv",
// and the partition number stored as an integer called "part".
// Function: clears data stored in f, defines f->m_partComm.
void NewPartition(FieldSharedPtr f, py::list &py_argv, int part)
{
    std::cout << std::endl << "Processing partition: " << part << std::endl;
    f->ClearField();
    int argc      = py::len(py_argv);
    char **argv   = ConvertCommandLine(py_argv);
    f->m_partComm = std::shared_ptr<FieldConvertComm>(
        new FieldConvertComm(argc, argv, f->m_nParts, part));
}

// Wrapper around the Field constructor
// Args: FieldConvert command line arguments.
// Function: constructs a FieldSharedPtr "f", defines f->m_defComm if
// nparts specified, and stores a line arguments in f->m_vm.
// Returns: f
FieldSharedPtr Field_Init(py::list &py_argv, int nparts = 0,
                          int output_points = 0, int output_points_hom_z = 0,
                          bool error = false, bool forceoutput = false,
                          std::string domain = "", bool noequispaced = false,
                          int npz = 0, std::string onlyshape = "",
                          int part_only = 0, int part_only_overlapping = 0,
                          bool useSessionVariables = false,
                          bool useSessionExpansion = false,
                          bool verbose             = false)
{
    // Construct shared pointer to a Field object.
    std::shared_ptr<Field> f = MemoryManager<Field>::AllocateSharedPtr();

    // Get argc and argv from the Python command line.
    int argc    = py::len(py_argv);
    char **argv = ConvertCommandLine(py_argv);

    // Define the MPI Communicator.
    f->m_comm =
        LibUtilities::GetCommFactory().CreateInstance("Serial", argc, argv);

    if (nparts)
    {
        f->m_vm.insert(
            std::make_pair("nparts", po::variable_value(nparts, false)));
        if (nparts > 1)
        {
            f->m_nParts  = nparts;
            f->m_defComm = f->m_comm;
        }
    }

    // Populate m_vm.
    if (output_points)
    {
        f->m_vm.insert(std::make_pair(
            "output-points", po::variable_value(output_points, false)));
    }

    if (output_points_hom_z)
    {
        f->m_vm.insert(
            std::make_pair("output-points-hom-z",
                           po::variable_value(output_points_hom_z, false)));
    }

    if (error)
    {
        f->m_vm.insert(std::make_pair("error", po::variable_value()));
    }

    if (forceoutput)
    {
        f->m_vm.insert(std::make_pair("forceoutput", po::variable_value()));
    }

    if (domain.size())
    {
        f->m_vm.insert(
            std::make_pair("range", po::variable_value(domain, false)));
    }

    if (noequispaced)
    {
        f->m_vm.insert(std::make_pair("noequispaced", po::variable_value()));
    }

    if (npz)
    {
        f->m_vm.insert(std::make_pair("npz", po::variable_value(npz, false)));
    }

    if (onlyshape.size())
    {
        f->m_vm.insert(
            std::make_pair("onlyshape", po::variable_value(onlyshape, false)));
    }

    if (part_only)
    {
        f->m_vm.insert(
            std::make_pair("part-only", po::variable_value(part_only, false)));
    }

    if (part_only_overlapping)
    {
        f->m_vm.insert(
            std::make_pair("part-only-overlapping",
                           po::variable_value(part_only_overlapping, false)));
    }

    if (useSessionVariables)
    {
        f->m_vm.insert(
            std::make_pair("useSessionVariables", po::variable_value()));
    }

    if (useSessionExpansion)
    {
        f->m_vm.insert(
            std::make_pair("useSessionExpansion", po::variable_value()));
    }

    if (verbose)
    {
        f->m_vm.insert(std::make_pair("verbose", po::variable_value()));
    }

    return f;
}

void export_Field()
{
    py::class_<Field, std::shared_ptr<Field>>("Field", py::no_init)
        .def("ClearField", &Field::ClearField)
        .def("NewPartition", &NewPartition)
        .def("__init__",
             py::make_constructor(
                 &Field_Init, py::default_call_policies(),
                 (py::arg("py_argv"), py::arg("nparts") = 0,
                  py::arg("output_points")       = 0,
                  py::arg("output_points_hom_z") = 0, py::arg("error") = false,
                  py::arg("forceoutput") = false, py::arg("domain") = "",
                  py::arg("noequispaced") = false, py::arg("npz") = 0,
                  py::arg("onlyshape") = "", py::arg("part_only") = 0,
                  py::arg("part_only_overlapping") = 0,
                  py::arg("useSessionVariables")   = false,
                  py::arg("useSessionExpansion")   = false,
                  py::arg("verbose")               = false)));
}
