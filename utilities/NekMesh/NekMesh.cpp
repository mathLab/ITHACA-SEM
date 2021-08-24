////////////////////////////////////////////////////////////////////////////////
//
//  File: NekMesh.cpp
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
//  Description: Mesh conversion utility.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <boost/algorithm/string.hpp>
#include <LibUtilities/BasicConst/GitRevision.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <boost/program_options.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/format.hpp>

#include <NekMesh/Module/Module.h>

using namespace std;
using namespace Nektar::NekMesh;

namespace po = boost::program_options;
namespace ip = boost::asio::ip;

int main(int argc, char* argv[])
{
    po::options_description desc("Available options");
    desc.add_options()
        ("help,h",         "Produce this help message.")
        ("modules-list,l", "Print the list of available modules.")
        ("modules-opt,p",  po::value<string>(),
             "Print options for a module.")
        ("module,m",       po::value<vector<string> >(),
             "Specify modules which are to be used.")
        ("verbose,v",      "Enable verbose output.");

    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-file",   po::value<vector<string> >(), "Input filename");

    po::options_description cmdline_options;
    cmdline_options.add(hidden).add(desc);

    po::options_description visible("Allowed options");
    visible.add(desc);

    po::positional_options_description p;
    p.add("input-file", -1);

    po::variables_map vm;

    try
    {
        po::store(po::command_line_parser(argc, argv).
                  options(cmdline_options).positional(p).run(), vm);
        po::notify(vm);
    }
    catch (const std::exception& e)
    {
        cerr << e.what() << endl;
        cerr << desc;
        return 1;
    }

    // Create a logger.
    auto logOutput = std::make_shared<StreamOutput>(std::cout);
    Logger log(logOutput, vm.count("verbose") ? VERBOSE : INFO);

    // Print available modules.
    if (vm.count("modules-list"))
    {
        GetModuleFactory().PrintAvailableClasses(std::cerr);
        return 1;
    }

    if (vm.count("modules-opt"))
    {
        vector<string> tmp1;
        boost::split(tmp1, vm["modules-opt"].as<string>(), boost::is_any_of(":"));

        if (tmp1.size() != 2)
        {
            cerr << "ERROR: To specify a module, use one of in, out or proc "
                 << "together with the filename; for example in:vtk." << endl;
            return 1;
        }

        if (tmp1[0] != "in" && tmp1[0] != "out" && tmp1[0] != "proc")
        {
            cerr << "ERROR: Invalid module type " << tmp1[0] << endl;
            return 1;
        }

        ModuleType t;

        if (tmp1[0] == "in")
        {
            t = eInputModule;
        }
        else if (tmp1[0] == "out")
        {
            t = eOutputModule;
        }
        else if (tmp1[0] == "proc")
        {
            t = eProcessModule;
        }

        MeshSharedPtr m = std::shared_ptr<Mesh>(new Mesh());
        ModuleSharedPtr mod = GetModuleFactory().CreateInstance(
            ModuleKey(t, tmp1[1]), m);
        cerr << "Options for module " << tmp1[1] << ":" << endl;
        mod->SetLogger(log);
        mod->PrintConfig();
        return 1;
    }

    if (vm.count("help") || vm.count("input-file") != 1) {
        cerr << "Usage: NekMesh [options] inputfile.ext1 outputfile.ext2"
             << endl;
        cout << desc;
        return 1;
    }

    vector<string> inout = vm["input-file"].as<vector<string> >();

    if (inout.size() < 2)
    {
        cerr << "ERROR: You must specify an input and output file." << endl;
        return 1;
    }

    /*
     * Process list of modules. Each element of the vector of module strings can
     * be in the following form:
     *
     * modname:arg1=a:arg2=b:arg3=c:arg4:arg5=asd
     *
     * where the only required argument is 'modname', specifing the name of the
     * module to load.
     */

    MeshSharedPtr mesh = std::shared_ptr<Mesh>(new Mesh());

    // Add provenance information to mesh.
    stringstream ss;
    for(int i = 1; i < argc; i++)
    {
        ss << argv[i] << " ";
    }
    mesh->m_metadata["NekMeshCommandLine"] = ss.str();

    vector<ModuleSharedPtr> modules;
    vector<string>          modcmds;

    if (vm.count("module"))
    {
        modcmds = vm["module"].as<vector<string> >();
    }

    // Add input and output modules to beginning and end of this vector.
    modcmds.insert   (modcmds.begin(), inout[0]);
    modcmds.push_back(inout[1]);

    // Keep track of maximum string length for nicer verbose output.
    size_t maxModName = 0;

    for (int i = 0; i < modcmds.size(); ++i)
    {
        // First split each command by the colon separator.
        vector<string> tmp1;
        ModuleKey module;
        int offset = 1;

        boost::split(tmp1, modcmds[i], boost::is_any_of(":"));

        if (i == 0 || i == modcmds.size() - 1)
        {
            module.first = (i == 0 ? eInputModule : eOutputModule);

            // If no colon detected, automatically detect mesh type from
            // file extension. Otherwise override and use tmp1[1] as the
            // module to load. This also allows us to pass options to
            // input/output modules. So, for example, to override
            // filename.xml to be read as vtk, you use:
            //
            // filename.xml:vtk:opt1=arg1:opt2=arg2
            if (tmp1.size() == 1)
            {
                int    dot    = tmp1[0].find_last_of('.') + 1;
                string ext    = tmp1[0].substr(dot, tmp1[0].length() - dot);

                if (ext == "gz")
                {
                    string tmp = tmp1[0].substr(0, tmp1[0].length() - 3);
                    dot = tmp.find_last_of('.') + 1;
                    ext = tmp.substr(dot, tmp.length() - dot);
                }

                module.second = ext;
                tmp1.push_back(string(i == 0 ? "infile=" : "outfile=")+tmp1[0]);
            }
            else
            {
                module.second = tmp1[1];
                tmp1.push_back(string(i == 0 ? "infile=" : "outfile=")+tmp1[0]);
                offset++;
            }
        }
        else
        {
            module.first  = eProcessModule;
            module.second = tmp1[0];
        }

        // Create module.
        ModuleSharedPtr mod = GetModuleFactory().CreateInstance(module,mesh);
        mod->SetLogger(log);
        modules.push_back(mod);

        // Set options for this module.
        for (int j = offset; j < tmp1.size(); ++j)
        {
            vector<string> tmp2;
            boost::split(tmp2, tmp1[j], boost::is_any_of("="));

            if (tmp2.size() == 1)
            {
                mod->RegisterConfig(tmp2[0]);
            }
            else if (tmp2.size() == 2)
            {
                mod->RegisterConfig(tmp2[0], tmp2[1]);
            }
            else
            {
                cerr << "ERROR: Invalid module configuration: format is "
                     << "either :arg or :arg=val" << endl;
                return 1;
            }
        }

        // Ensure configuration options have been set.
        mod->SetDefaults();

        // Track maximum module name length.
        std::string modName = mod->GetModuleName();
        maxModName = std::max(maxModName, modName.length());
    }

    log.SetPrefixLen(maxModName);

    // Run mesh process.
    for (int i = 0; i < modules.size(); ++i)
    {
        Nektar::LibUtilities::Timer t;
        t.Start();
        try
        {
            modules[i]->GetLogger().SetPrefixLen(maxModName);
            modules[i]->Process();
        }
        catch (NekMeshError &e)
        {
            return 1;
        }
        t.Stop();

        log.SetPrefix(modules[i]->GetModuleName());
        log(VERBOSE) << "  - Elapsed time: "
                     << t.TimePerTest(1) << "s." << std::endl;
    }

    return 0;
}
