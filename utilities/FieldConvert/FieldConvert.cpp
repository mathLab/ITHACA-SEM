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
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <FieldUtils/Module.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::FieldUtils;

void CheckModules(vector<ModuleSharedPtr> &modules);

void PrintExecutionSequence(vector<ModuleSharedPtr> &modules);

void RunModule(ModuleSharedPtr module, po::variables_map &vm, bool verbose);

int main(int argc, char* argv[])
{
    LibUtilities::Timer    timer;
    timer.Start();

    po::options_description desc("Available options");
    desc.add_options()
        ("help,h",
            "Produce this help message.")
        ("modules-list,l",
            "Print the list of available modules.")
        ("output-points,n", po::value<int>(),
            "Output at n equipspaced points along the "
            "collapsed coordinates (for .dat, .vtu).")
        ("output-points-hom-z", po::value<int>(),
            "Number of planes in the z-direction for output of "
            "Homogeneous 1D expansion(for .dat, .vtu).")
        ("error,e",
            "Write error of fields for regression checking")
        ("forceoutput,f",
            "Force the output to be written without any checks")
        ("range,r", po::value<string>(),
            "Define output range i.e. (-r xmin,xmax,ymin,ymax,zmin,zmax) "
            "in which any vertex is contained.")
        ("noequispaced",
            "Do not use equispaced output.")
        ("nparts", po::value<int>(),
            "Define nparts if running serial problem to mimic "
            "parallel run with many partitions.")
        ("npz", po::value<int>(),
            "Used to define number of partitions in z for Homogeneous1D "
            "expansions for parallel runs.")
        ("onlyshape", po::value<string>(),
            "Only use element with defined shape type i.e. -onlyshape "
            " Tetrahedron")
        ("part-only", po::value<int>(),
            "Partition into specified npart partitions and exit")
        ("part-only-overlapping", po::value<int>(),
            "Partition into specified npart overlapping partitions and exit")
        ("modules-opt,p", po::value<string>(),
            "Print options for a module.")
        ("module,m", po::value<vector<string> >(),
            "Specify modules which are to be used.")
        ("useSessionVariables",
            "Use variables defined in session for output")
        ("useSessionExpansion",
            "Use expansion defined in session.")
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

#ifdef NEKTAR_DISABLE_BACKUPS
    vm.insert(std::make_pair("forceoutput", po::variable_value()));
#endif

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

        FieldSharedPtr f = std::shared_ptr<Field>(new Field());
        ModuleSharedPtr mod = GetModuleFactory().CreateInstance(
            ModuleKey(t, tmp1[1]), f);
        cerr << "Options for module " << tmp1[1] << ":" << endl;
        mod->PrintConfig();
        return 1;
    }

    if (vm.count("help") || vm.count("input-file") != 1)
    {
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
    FieldSharedPtr f = std::shared_ptr<Field>(new Field());
    int nParts = 1;
    int MPInprocs = 1;
    int MPIrank   = 0;
    LibUtilities::CommSharedPtr MPIComm;

    if (LibUtilities::GetCommFactory().ModuleExists("ParallelMPI"))
    {
        // get hold of parallel communicator first
        MPIComm = LibUtilities::GetCommFactory().CreateInstance(
                                                    "ParallelMPI", argc, argv);

        if(vm.count("nparts"))
        {
            //work out number of processors to run in serial over partitions
            MPInprocs = MPIComm->GetSize();
            MPIrank   = MPIComm->GetRank();

            nParts = vm["nparts"].as<int>();

            f->m_comm = LibUtilities::GetCommFactory().CreateInstance(
                                                    "Serial", argc, argv);
        }
        else
        {
            f->m_comm = MPIComm;
        }
    }
    else
    {
        if(vm.count("nparts"))
        {
            nParts = vm["nparts"].as<int>();
        }

        f->m_comm = LibUtilities::GetCommFactory().CreateInstance(
                                                    "Serial", argc, argv);
    }

    vector<ModuleSharedPtr> modules;
    vector<string>          modcmds;
    ModuleKey               module;
    ModuleSharedPtr         mod;

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
    string outfilename;

    for (int i = 0; i < modcmds.size(); ++i)
    {
        // First split each command by the colon separator.
        vector<string> tmp1;
        int offset = 1;

        boost::split(tmp1, modcmds[i], boost::is_any_of(":"));

        if (i < nInput || i == modcmds.size() - 1)
        {
            //assume all modules are input unless last, or specified to be :out
            module.first = (i < nInput ? eInputModule : eOutputModule);
            if (tmp1.size() > 1 && tmp1.back()=="out")
            {
                module.first = eOutputModule;
                tmp1.pop_back();
            }

            // If no colon detected, automatically detect mesh type from
            // file extension. Otherwise override and use tmp1[1] as the
            // module to load. This also allows us to pass options to
            // input/output modules. So, for example, to override
            // filename.xml to be read as vtk, you use:
            //
            // filename.xml:vtk:opt1=arg1:opt2=arg2
            if (tmp1.size() == 1)
            {
                // First, let's try to guess the input format if we're dealing
                // with input files.
                string guess;

                if (module.first == eInputModule)
                {
                    guess = InputModule::GuessFormat(tmp1[0]);
                }

                // Found file type.
                if (guess != "")
                {
                    if (f->m_verbose)
                    {
                        cout << "Using input module " << guess << " for: "
                             << tmp1[0] << endl;
                    }

                    module.second = guess;
                    tmp1.push_back(string("infile="+tmp1[0]));
                }
                else
                {
                    int    dot = tmp1[0].find_last_of('.') + 1;
                    string ext = tmp1[0].substr(dot, tmp1[0].length() - dot);

                    // Remove trailing separator from extension to allow
                    //    folder inputs using file.fld/
                    if(ext.back() == fs::path::preferred_separator)
                    {
                        ext.pop_back();
                    }

                    if(ext == "gz")
                    {
                        string tmp2 = tmp1[0].substr(0,dot-1);
                        dot = tmp2.find_last_of('.') + 1;
                        ext = tmp1[0].substr(dot,tmp1[0].length()-dot);
                    }

                    module.second = ext;
                    tmp1.push_back(string(module.first == eInputModule ? "infile=" :
                                          "outfile=")  +tmp1[0]);
                }
            }
            else
            {
                module.second = tmp1[1];
                tmp1.push_back(string(module.first == eInputModule ? "infile=" : "outfile=")
                               + tmp1[0]);
                offset++;
            }
        }
        else
        {
            module.first  = eProcessModule;
            module.second = tmp1[0];
        }

        // Create module.
        mod = GetModuleFactory().CreateInstance(module, f);
        modules.push_back(mod);

        if (module.first == eInputModule)
        {
            inputModule = std::dynamic_pointer_cast<InputModule>(mod);
            inputModule->AddFile(module.second, tmp1[0]);
        }
        else if(module.first == eOutputModule)
        {
            outfilename = tmp1[0];
            if(nParts > 1)
            {
                // if nParts is specified then ensure output modules
                // write out mutipile files
                mod->RegisterConfig("writemultiplefiles");
            }
        }

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
                abort();
            }
        }

        // Ensure configuration options have been set.
        mod->SetDefaults();
    }

    // Include dummy module to create m_exp
    module.first  = eProcessModule;
    module.second = string("createExp");
    mod = GetModuleFactory().CreateInstance(module, f);
    modules.push_back(mod);
    mod->SetDefaults();

    // Include equispacedoutput module if needed
    Array< OneD, int>  modulesCount(SIZE_ModulePriority,0);
    for (int i = 0; i < modules.size(); ++i)
    {
        ++modulesCount[modules[i]->GetModulePriority()];
    }
    if( modulesCount[eModifyPts] != 0 &&
        modulesCount[eCreatePts] == 0 &&
        modulesCount[eConvertExpToPts] == 0)
    {
        module.first  = eProcessModule;
        module.second = string("equispacedoutput");
        mod = GetModuleFactory().CreateInstance(module, f);
        modules.push_back(mod);
        mod->SetDefaults();
    }

    // Check if modules provided are compatible
    CheckModules(modules);
    // Can't have ContField with range option (because of boundaries)
    if (vm.count("range") && f->m_declareExpansionAsContField)
    {
        ASSERTL0(false, "Can't use range option with module requiring "
            "a continuous expansion.");
    }

    bool verbose = (f->m_verbose && f->m_comm->TreatAsRankZero());
    if(verbose)
    {
        PrintExecutionSequence(modules);
    }


    // Loop on partitions if required
    LibUtilities::CommSharedPtr defComm = f->m_comm;
    LibUtilities::CommSharedPtr partComm;
    for(int p = MPIrank; p < nParts; p += MPInprocs)
    {
        // write out which partition is being processed and defined a
        // new serial communicator
        if(nParts > 1)
        {
            cout << endl << "Processing partition: " << p << endl;
            
            int rank = p;
            f->ClearField();
            partComm = std::shared_ptr<FieldConvertComm>(
                             new FieldConvertComm(argc, argv, nParts,rank));
        }
        
        // Run field process.
        for (int n = 0; n < SIZE_ModulePriority; ++n)
        {
            ModulePriority priority = static_cast<ModulePriority>(n);

            if(nParts > 1)
            {
                if(((priority == eCreateGraph)||(priority == eOutput)))
                {
                    f->m_comm = partComm;
                }
                else
                {
                    f->m_comm = defComm;
                }
            }

            for (int i = 0; i < modules.size(); ++i)
            {
                if(modules[i]->GetModulePriority() == priority)
                {
                    RunModule(modules[i], vm, verbose);
                }
            }
        }
    }

    // write out Info file if required.
    if (nParts > 1)
    {
        int i;
        // check to see if we have created a fld file.
        for (i = 0; i < modules.size(); ++i)
        {
            if (boost::iequals(modules[i]->GetModuleName(), "OutputFld"))
            {
                break;
            }
        }

        if (i != modules.size())
        {
            if (MPInprocs > 1)
            {
                MPIComm->Block();
            }

            if (MPIrank == 0)
            {
                module.first  = eOutputModule;
                module.second = string("info");
                mod           = GetModuleFactory().CreateInstance(module, f);

                mod->RegisterConfig("nparts",
                                    boost::lexical_cast<string>(nParts));
                mod->SetDefaults();

                if (f->m_writeBndFld)
                {
                    // find ending of output file and insert _b1, _b2
                    int dot = outfilename.find_last_of('.') + 1;
                    string ext =
                        outfilename.substr(dot, outfilename.length() - dot);
                    string name = outfilename.substr(0, dot - 1);

                    for (int b = 0; b < f->m_bndRegionsToWrite.size(); ++b)
                    {
                        string outfilenew = name + "_b" +
                                            boost::lexical_cast<string>(
                                                f->m_bndRegionsToWrite[b]) +
                                            "." + ext;
                        mod->RegisterConfig("outfile", outfilenew);
                        RunModule(mod, vm, verbose);
                    }
                }
                else
                {
                    mod->RegisterConfig("outfile", outfilename);
                    RunModule(mod, vm, verbose);
                }
            }
        }
    }

    if(verbose)
    {
        timer.Stop();
        NekDouble cpuTime = timer.TimePerTest(1);

        stringstream ss;
        ss << cpuTime << "s";
        cout << "Total CPU Time: " << setw(8) << left
             << ss.str() << endl;
    }

    if(MPInprocs > 1)
    {
        MPIComm->Block();
        MPIComm->Finalise();
    }

    return 0;
}

// This function checks validity conditions for the list of modules provided
void CheckModules(vector<ModuleSharedPtr> &modules)
{
    // Count number of modules by priority
    Array< OneD, int>  modulesCount(SIZE_ModulePriority,0);
    for (int i = 0; i < modules.size(); ++i)
    {
        ++modulesCount[modules[i]->GetModulePriority()];
    }

    // Modules of type eModifyFieldData require a eCreateFieldData module
    if( modulesCount[eModifyFieldData] != 0 &&
        modulesCount[eCreateFieldData] == 0)
    {
        stringstream ss;
        ss << "Module(s): ";
        for (int i = 0; i < modules.size(); ++i)
        {
            if(modules[i]->GetModulePriority() == eModifyFieldData)
            {
                ss << modules[i]->GetModuleName()<<" ";
            }
        }
        ss << "require fld input.";
        ASSERTL0(false, ss.str());
    }

    // Modules of type eFillExp require eCreateGraph without eCreateFieldData
    if( modulesCount[eFillExp] != 0)
    {
        if( modulesCount[eCreateGraph]       == 0 ||
            modulesCount[eCreateFieldData]   != 0)
        {
            stringstream ss;
            ss << "Module(s): ";
            for (int i = 0; i < modules.size(); ++i)
            {
                if(modules[i]->GetModulePriority() == eFillExp)
                {
                    ss << modules[i]->GetModuleName()<<" ";
                }
            }
            ss << "require xml input without fld input.";
            ASSERTL0(false, ss.str());
        }
    }

    // Modules of type eModifyExp and eBndExtraction
    //      require a eCreateGraph module
    if( (modulesCount[eModifyExp] != 0 || modulesCount[eBndExtraction] != 0) &&
        modulesCount[eCreateGraph] == 0)
    {
        stringstream ss;
        ss << "Module(s): ";
        for (int i = 0; i < modules.size(); ++i)
        {
            if(modules[i]->GetModulePriority() == eModifyExp ||
               modules[i]->GetModulePriority() == eBndExtraction)
            {
                ss << modules[i]->GetModuleName()<<" ";
            }
        }
        ss << "require xml input.";
        ASSERTL0(false, ss.str());
    }

    // Modules of type eCreatePts should not be used with xml or fld inputs
    if( modulesCount[eCreatePts] != 0)
    {
        if(modulesCount[eCreateGraph]!=0 || modulesCount[eCreateFieldData]!=0)
        {
            stringstream ss;
            ss << "Module(s): ";
            for (int i = 0; i < modules.size(); ++i)
            {
                if(modules[i]->GetModulePriority() == eCreatePts)
                {
                    ss << modules[i]->GetModuleName()<<" ";
                }
            }
            ss << "should not use xml or fld inputs.";
            ASSERTL0(false, ss.str());
        }
    }

    // Modules of type eConvertExpToPts require eCreateGraph, but are not
    //    compatible with eBndExtraction
    if( modulesCount[eConvertExpToPts] != 0)
    {
        if( modulesCount[eCreateGraph] == 0)
        {
            stringstream ss;
            ss << "Module(s): ";
            for (int i = 0; i < modules.size(); ++i)
            {
                if(modules[i]->GetModulePriority() == eConvertExpToPts)
                {
                    ss << modules[i]->GetModuleName()<<" ";
                }
            }
            ss << "require xml input.";
            ASSERTL0(false, ss.str());
        }
        if( modulesCount[eBndExtraction] != 0)
        {
            stringstream ss;
            ss << "Module(s): ";
            for (int i = 0; i < modules.size(); ++i)
            {
                if(modules[i]->GetModulePriority() == eBndExtraction)
                {
                    ss << modules[i]->GetModuleName()<<" ";
                }
            }
            ss << "is not compatible with module(s): ";
            for (int i = 0; i < modules.size(); ++i)
            {
                if(modules[i]->GetModulePriority() == eConvertExpToPts)
                {
                    ss << modules[i]->GetModuleName()<<" ";
                }
            }
            ss << ".";
            ASSERTL0(false, ss.str());
        }
    }
}

void PrintExecutionSequence(vector<ModuleSharedPtr> &modules)
{
    bool first = true;
    cout << "Execution sequence:" << endl;
    for (int n = 0; n < SIZE_ModulePriority; ++n)
    {
        ModulePriority priority = static_cast<ModulePriority>(n);
        for (int i = 0; i < modules.size(); ++i)
        {
            if(modules[i]->GetModulePriority() == priority)
            {
                if(first)
                {
                    cout << "\t"  << modules[i]->GetModuleName();
                    first = false;
                }
                else
                {
                    cout << " -> " << modules[i]->GetModuleName();
                }
            }
        }
    }
    cout << endl;
}

void RunModule(ModuleSharedPtr module, po::variables_map &vm, bool verbose)
{
    LibUtilities::Timer moduleTimer;
    if(verbose)
    {
        moduleTimer.Start();

        cout << module->GetModuleName() << ": "
             << module->GetModuleDescription() << endl;
    }
    module->Process(vm);
    cout.flush();
    if(verbose)
    {
        moduleTimer.Stop();
        NekDouble cpuTime = moduleTimer.TimePerTest(1);

        stringstream ss;
        ss << cpuTime << "s";
        cout << module->GetModuleName()
             << " CPU Time: " << setw(8) << left
             << ss.str() << endl;
    }
}
