////////////////////////////////////////////////////////////////////////////////
//
//  File: FieldConvert.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and laimitations under
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
//  Description: Field conversion utility.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include "Module.h"

using namespace std;
using namespace Nektar;
using namespace Nektar::Utilities;

int main(int argc, char* argv[])
{
    Timer     timer;
    timer.Start();
    po::options_description desc("Available options");
    desc.add_options()
        ("help,h",
                "Produce this help message.")
        ("modules-list,l",
                "Print the list of available modules.")
        ("output-points,n", po::value<int>(),
                "Output at n equipspaced points along the collapsed coordinates (for .dat, .vtk).")
        ("output-points-hom-z", po::value<int>(),
                "Number of planes in the z-direction for output of Homogeneous 1D expansion(for .dat, .vtk).")
        ("error,e",
                "Write error of fields for regression checking")
        ("forceoutput,f",
                "Force the output to be written without any checks")
        ("range,r", po::value<string>(),
                "Define output range i.e. (-r xmin,xmax,ymin,ymax,zmin,zmax) "
                "in which any vertex is contained.")
        ("noequispaced","Do not use equispaced output. Currently stops the output-points option")
        ("nprocs", po::value<int>(),
                "Used to define nprocs if running serial problem to mimic "
                "parallel run.")
        ("onlyshape", po::value<string>(),
                 "Only use element with defined shape type i.e. -onlyshape "
                 " Tetrahedron")
        ("part-only", po::value<int>(),
                "Partition into specfiied npart partitions and exit")
        ("part-only-overlapping", po::value<int>(),
                "Partition into specfiied npart overlapping partitions and exit")
        ("procid", po::value<int>(),
                "Process as single procid of a partition of size nproc "
                "(-nproc must be specified).")
        ("modules-opt,p", po::value<string>(),
                "Print options for a module.")
        ("module,m", po::value<vector<string> >(),
                "Specify modules which are to be used.")
        ("shared-filesystem,s", "Using shared filesystem.")
        ("useSessionVariables",
                "Use variables defined in session for output")
        ("verbose,v",
                "Enable verbose mode.");

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

    // Print available modules.
    if (vm.count("modules-list"))
    {
        GetModuleFactory().PrintAvailableClasses(std::cerr);
        return 1;
    }

    if (vm.count("modules-opt"))
    {
        vector<string> tmp1;
        boost::split(tmp1, vm["modules-opt"].as<string>(),
                     boost::is_any_of(":"));

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

        FieldSharedPtr f = boost::shared_ptr<Field>(new Field());
        ModuleSharedPtr mod = GetModuleFactory().CreateInstance(
            ModuleKey(t, tmp1[1]), f);
        cerr << "Options for module " << tmp1[1] << ":" << endl;
        mod->PrintConfig();
        return 1;
    }

    if (vm.count("help") || vm.count("input-file") != 1) {
        cerr << "Usage: FieldConvert [options] inputfile.ext1 outputfile.ext2"
             << endl;
        cout << desc;
        cout << endl;
        cout << "Example Usage: \n" << endl;
        cout << "\t FieldConvert -m vorticity file.xml file.fld file_vort.fld "
             << endl;
        cout << "(This will add vorticity to file file.fld and put it in a "
                "new file file_vort.fld) " << endl;
        cout << endl;
        cout << "\t FieldConvert file.xml file_vort.fld file_vort.dat " << endl;
        cout << "(process file_vort.fld and make a tecplot output "
                "file_vort.dat) " << endl;

        return 1;
    }

    ASSERTL0(vm.count("input-file"),
             "Must specify input(s) and/or output file.");
    vector<string> inout = vm["input-file"].as<vector<string> >();


    /*
     * Process list of modules. Each element of the vector of module
     * strings can be in the following form:
     *
     * modname:arg1=a:arg2=b:arg3=c:arg4:arg5=asd
     *
     * where the only required argument is 'modname', specifing the
     * name of the module to load.
     */

    FieldSharedPtr f = boost::shared_ptr<Field>(new Field());
    if (LibUtilities::GetCommFactory().ModuleExists("ParallelMPI"))
    {
        if(vm.count("procid"))
        {
            int nprocs, rank;

            ASSERTL0(vm.count("nprocs"),
                     "Must specify --nprocs when using --procid option");
            nprocs = vm["nprocs"].as<int>();
            rank   = vm["procid"].as<int>();

            f->m_comm = boost::shared_ptr<FieldConvertComm>(
                                new FieldConvertComm(argc, argv, nprocs,rank));

            // Set forceoutput option. Otherwise only procid 0 will write file
            vm.insert(std::make_pair("forceoutput", po::variable_value()));
        }
        else
        {
            f->m_comm = LibUtilities::GetCommFactory().CreateInstance(
                                                    "ParallelMPI", argc, argv);
        }
    }
    else
    {
        f->m_comm = LibUtilities::GetCommFactory().CreateInstance(
                                                    "Serial", argc, argv);

    }

    vector<ModuleSharedPtr> modules;
    vector<string>          modcmds;


    if (vm.count("verbose"))
    {
        f->m_verbose = true;
    }

    if (vm.count("module"))
    {
        modcmds = vm["module"].as<vector<string> >();
    }

    // Add input and output modules to beginning and end of this vector.
    modcmds.insert(modcmds.begin(), inout.begin(), inout.end()-1);
    modcmds.push_back(*(inout.end()-1));

    int nInput = inout.size()-1;

    // For special case of part-only or part-only-overlapping options
    // only require a single input file and so reset the nInputs to be
    // of size inout.size(). Since the code will exit before reaching
    // any output module this seems to work as expected
    if(vm.count("part-only")||vm.count("part-only-overlapping"))
    {
        nInput = inout.size();
    }

    InputModuleSharedPtr inputModule;

    for (int i = 0; i < modcmds.size(); ++i)
    {
        // First split each command by the colon separator.
        vector<string> tmp1;
        ModuleKey module;
        int offset = 1;

        boost::split(tmp1, modcmds[i], boost::is_any_of(":"));

        if (i < nInput || i == modcmds.size() - 1)
        {
            module.first = (i < nInput ? eInputModule : eOutputModule);

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

                if(ext == "gz")
                {
                    string tmp2 = tmp1[0].substr(0,dot-1);
                    dot = tmp2.find_last_of('.') + 1;
                    ext = tmp1[0].substr(dot,tmp1[0].length()-dot);
                }

                module.second = ext;
                tmp1.push_back(string(i < nInput ? "infile=" : "outfile=")
                               +tmp1[0]);
            }
            else
            {
                module.second = tmp1[1];
                tmp1.push_back(string(i < nInput ? "infile=" : "outfile=")
                               +tmp1[0]);
                offset++;
            }
        }
        else
        {
            module.first  = eProcessModule;
            module.second = tmp1[0];
        }

        // Create module.
        ModuleSharedPtr mod;
        mod = GetModuleFactory().CreateInstance(module, f);
        modules.push_back(mod);

        if (i < nInput)
        {
            inputModule = boost::dynamic_pointer_cast<InputModule>(mod);
            inputModule->AddFile(module.second, tmp1[0]);
        }

        // Set options for this module.
        for (int j = offset; j < tmp1.size(); ++j)
        {
            vector<string> tmp2;
            boost::split(tmp2, tmp1[j], boost::is_any_of("="));

            if (tmp2.size() == 1)
            {
                mod->RegisterConfig(tmp2[0], "1");
            }
            else if (tmp2.size() == 2)
            {
                mod->RegisterConfig(tmp2[0], tmp2[1]);
            }
            else
            {
                cerr << "ERROR: Invalid module configuration: format is "
                     << "either :arg or :arg=val" << endl;
                abort();
            }
        }

        // Ensure configuration options have been set.
        mod->SetDefaults();
    }

    // If any output module has to reset points then set intput modules to match
   if(vm.count("noequispaced"))
    {
        for (int i = 0; i < modules.size(); ++i)
        {
            modules[i]->SetRequireEquiSpaced(false);
        }
    }
    else
    {
        bool RequiresEquiSpaced = false;
        for (int i = 0; i < modules.size(); ++i)
        {
            if(modules[i]->GetRequireEquiSpaced())
            {
                RequiresEquiSpaced = true;
            }
        }
        if (RequiresEquiSpaced)
        {
            for (int i = 0; i < modules.size(); ++i)
            {
                modules[i]->SetRequireEquiSpaced(true);
            }
        }
    }
    
    // Run field process.
    for (int i = 0; i < modules.size(); ++i)
    {
        modules[i]->Process(vm);
        cout.flush();
    }

    if(f->m_verbose)
    {
        if(f->m_comm->GetRank() == 0)
        {
            timer.Stop();
            NekDouble cpuTime = timer.TimePerTest(1);
            
            stringstream ss;
            ss << cpuTime << "s";
            cout << "Total CPU Time: " << setw(8) << left
                 << ss.str() << endl;
            cpuTime = 0.0;
        }
        
    }
    return 0;
}
