#include <LibUtilities/BasicUtils/SessionReader.h>
#include <NekPyConfig.hpp>

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

        ;
}
