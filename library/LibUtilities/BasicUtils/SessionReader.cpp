///////////////////////////////////////////////////////////////////////////////
//
// File SessionReader.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Session reader
//
///////////////////////////////////////////////////////////////////////////////

#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicConst/GitRevision.h>

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/algorithm/string.hpp>
#include <tinyxml.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/MeshPartition.h>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/FileSystem.h>

#include <boost/program_options.hpp>
#include <boost/format.hpp>

#ifndef NEKTAR_VERSION
#define NEKTAR_VERSION "Unknown"
#endif

namespace po = boost::program_options;
namespace io = boost::iostreams;

namespace Nektar
{
    namespace LibUtilities
    {
        /**
         * @class SessionReader
         *
         * This class provides an interface to Nektar++-specific content in a
         * supplied XML document. It also initialises a Nektar++ session
         * including setting up communication for parallel execution and where
         * necessary partitioning the supplied mesh for running across multiple
         * processes.
         *
         * A session should be initialised at the beginning of a user's
         * application by passing the command-line arguments. This not only
         * allows the SessionReader to extract the name of the XML document to
         * load containing Nektar++ session information, but also supplies the
         * MPI arguments necessary for setting up parallel communication. The
         * SessionReader should be initialised using the #CreateInstance
         * function:
         * @code
         * LibUtilities::SessionReaderSharedPtr vSession
         *          = LibUtilities::SessionReader::CreateInstance(argc, argv);
         * @endcode
         * The instance \c vSession can now be passed to other key Nektar++
         * components during their construction.
         * @note At the end of the user application, it is important to call the
         * #Finalise routine in order to finalise any MPI communication and
         * correctly free resources.
         *
         * The SessionReader class provides streamlined, validated access to
         * session parameters, solver information and functions defined within a
         * Nektar++ XML document. The available routines and their usage is
         * documented below.
         *
         * In the case of solver information properties, the classes to which
         * these parameters are pertinent may register with the SessionReader
         * class the set of valid values for a given property. Such values may
         * also be associated with an enumeration value for more transparent use
         * of the property values in code.
         */

        /**
         * This map of maps stores the list of valid string values for a number
         * of solver information parameters. The top level map connects
         * different parameter names to their list of possible values. The list
         * of possible values is also a map, mapping a valid string to a
         * corresponding enum value.
         *
         * This list is populated through the #RegisterEnumValue static member
         * function which is called statically from various classes to register
         * the valid values for solver info parameters associated with them. The
         * map is therefore fully populated before the SessionReader class is
         * instantiated and a file is read in and parsed.
         */
        EnumMapList& SessionReader::GetSolverInfoEnums()
        {
            static EnumMapList solverInfoEnums;
            return solverInfoEnums;
        }


        /**
         * List of default values for solver information parameters to be used
         * in the case of them not being provided.
         *
         * This list is populated through the #RegisterDefaultSolverInfo static
         * member variable which is called statically from various classes to
         * register the default value for a given parameter.
         */
        SolverInfoMap& SessionReader::GetSolverInfoDefaults()
        {
            static SolverInfoMap solverInfoMap;
            return solverInfoMap;
        }


        /**
         * List of values for GlobalSysSoln parameters to be used to override
         * details given in SolverInfo
         *
         * This list is populated by ReadGlobalSysSoln if the
         * GLOBALSYSSOLNINFO section is defined in the input file.
         * This List allows for details to define for the Global Sys
         * solver for each variable. 
         */
        GloSysSolnInfoList& SessionReader::GetGloSysSolnList()
        {
            static GloSysSolnInfoList gloSysSolnInfoList;
            return gloSysSolnInfoList;
        }

        /**
         * Lists the possible command-line argument which can be specified for
         * this executable.
         *
         * This list is populated through the #RegisterCmdLineArgument static
         * member function which is called statically from various classes to
         * register command-line arguments they need.
         */
        CmdLineArgMap& SessionReader::GetCmdLineArgMap()
        {
            static CmdLineArgMap cmdLineArguments;
            return cmdLineArguments;
        }


        /**
         * This constructor parses the command-line arguments given to the user
         * application to set up any MPI communication, read supplied XML
         * session files, and partition meshes where necessary.
         *
         * @param   argc        Number of command-line arguments
         * @param   argv        Array of command-line arguments
         */
        SessionReader::SessionReader(int argc, char *argv[])
        {
            m_xmlDoc    = 0;
            m_filenames = ParseCommandLineArguments(argc, argv);

            ASSERTL0(m_filenames.size() > 0, "No session file(s) given.");

            m_sessionName = ParseSessionName(m_filenames);

            // Create communicator
            CreateComm(argc, argv);

            // If running in parallel change the default global sys solution
            // type.
            if (m_comm->GetSize() > 1)
            {
                GetSolverInfoDefaults()["GLOBALSYSSOLN"] = 
                    "IterativeStaticCond";
            }
        }


        /**
         *
         */
        SessionReader::SessionReader(
            int                             argc, 
            char                           *argv[], 
            const std::vector<std::string> &pFilenames, 
            const CommSharedPtr            &pComm)
        {
            ASSERTL0(pFilenames.size() > 0, "No filenames specified.");

            ParseCommandLineArguments(argc, argv);
            m_xmlDoc      = 0;
            m_filenames   = pFilenames;

            m_sessionName = ParseSessionName(m_filenames);

            // Create communicator
            if (!pComm.get())
            {
                CreateComm(argc, argv);
            }
            else
            {
                m_comm = pComm;

                if (m_comm->GetSize() > 1)
                {
                    GetSolverInfoDefaults()["GLOBALSYSSOLN"] = 
                        "IterativeStaticCond";
                }
            }

            // If running in parallel change the default global sys solution
            // type.
            if (m_comm->GetSize() > 1)
            {
                GetSolverInfoDefaults()["GLOBALSYSSOLN"] = 
                    "IterativeStaticCond";
            }
        }


        /**
         *
         */
        SessionReader::~SessionReader()
        {
            delete m_xmlDoc;
        }


        /**
         * Performs the main initialisation of the object. The XML file provided
         * on the command-line is loaded and any mesh partitioning is done. The
         * resulting process-specific XML file (containing the process's
         * geometry partition) is then reloaded and parsed.
         */
        void SessionReader::InitSession()
        {
            m_exprEvaluator.SetRandomSeed((m_comm->GetRank() + 1) * time(NULL));

            // Split up the communicator
            PartitionComm();

            // Partition mesh
            PartitionMesh();

            // Parse the XML data in #m_xmlDoc
            ParseDocument();

            // Override SOLVERINFO and parameters with any specified on the
            // command line.
            CmdLineOverride();

            // In verbose mode, print out parameters and solver info sections
            if (m_verbose && m_comm)
            {
                if (m_comm->GetRank() == 0 && m_parameters.size() > 0)
                {
                    cout << "Parameters:" << endl;
                    ParameterMap::iterator x;
                    for (x = m_parameters.begin(); x != m_parameters.end(); ++x)
                    {
                        cout << "\t" << x->first << " = " << x->second << endl;
                    }
                    cout << endl;
                }

                if (m_comm->GetRank() == 0 && m_solverInfo.size() > 0)
                {
                    cout << "Solver Info:" << endl;
                    SolverInfoMap::iterator x;
                    for (x = m_solverInfo.begin(); x != m_solverInfo.end(); ++x)
                    {
                        cout << "\t" << x->first << " = " << x->second << endl;
                    }
                    cout << endl;
                }
            }
        }


        /**
         * @brief Parses the command-line arguments for known options and
         * filenames.
         */
        std::vector<std::string> SessionReader::ParseCommandLineArguments(
            int argc, char *argv[])
        {
            // List the publically visible options (listed using --help)
            po::options_description desc("Allowed options");
            desc.add_options()
                ("verbose,v",    "be verbose")
                ("version,V",    "print version information")
                ("help,h",       "print this help message")
                ("solverinfo,I", po::value<vector<std::string> >(), 
                                 "override a SOLVERINFO property")
                ("parameter,P",  po::value<vector<std::string> >(),
                                 "override a parameter")
                ("shared-filesystem,s", "Using shared filesystem.")
                ("npx",          po::value<int>(),
                                 "number of procs in X-dir")
                ("npy",          po::value<int>(),
                                 "number of procs in Y-dir")
                ("npz",          po::value<int>(),
                                 "number of procs in Z-dir")
                ("nsz",          po::value<int>(),
                                 "number of slices in Z-dir")
                ("part-only",    po::value<int>(),
                                 "only partition mesh into N partitions.")
                ("part-info",    "Output partition information")
            ;
            
            CmdLineArgMap::const_iterator cmdIt;
            for (cmdIt  = GetCmdLineArgMap().begin();
                 cmdIt != GetCmdLineArgMap().end(); ++cmdIt)
            {
                std::string names = cmdIt->first;
                if (cmdIt->second.shortName != "")
                {
                    names += "," + cmdIt->second.shortName;
                }
                if (cmdIt->second.isFlag)
                {
                    desc.add_options()
                        (names.c_str(), cmdIt->second.description.c_str())
                    ;
                }
                else
                {
                    desc.add_options()
                        (names.c_str(), po::value<std::string>(),
                         cmdIt->second.description.c_str())
                    ;
                }
            }

            // List hidden options (e.g. session file arguments are not actually
            // specified using the input-file option by the user).
            po::options_description hidden("Hidden options");
            hidden.add_options()
                    ("input-file", po::value< vector<string> >(), 
                                   "input filename")
            ;

            // Combine all options for the parser
            po::options_description all("All options");
            all.add(desc).add(hidden);

            // Session file is a positional option
            po::positional_options_description p;
            p.add("input-file", -1);

            // Parse the command-line options
            po::parsed_options parsed = po::command_line_parser(argc, argv).
                                                options(all).
                                                positional(p).
                                                allow_unregistered().
                                                run();

            // Extract known options to map and update
            po::store(parsed, m_cmdLineOptions);
            po::notify(m_cmdLineOptions);

            // Help message
            if (m_cmdLineOptions.count("help"))
            {
                cout << desc;
                exit(0);
            }

            // Version information
            if (m_cmdLineOptions.count("version"))
            {
                cout << "Nektar++ version " << NEKTAR_VERSION;

                if (NekConstants::kGitSha1 != "GITDIR-NOTFOUND")
                {
                    string sha1(NekConstants::kGitSha1);
                    string branch(NekConstants::kGitBranch);
                    boost::replace_all(branch, "refs/heads/", "");

                    cout << " (git changeset " << sha1.substr(0, 8) << ", ";

                    if (branch == "")
                    {
                        cout << "detached head";
                    }
                    else
                    {
                        cout << "head " << branch;
                    }

                    cout << ")";
                }

                cout << endl;
                exit(0);
            }

            // Enable verbose mode
            if (m_cmdLineOptions.count("verbose"))
            {
                m_verbose = true;
            }
            else
            {
                m_verbose = false;
            }
            
            // Print a warning for unknown options
            std::vector< po::basic_option<char> >::iterator x;
            for (x = parsed.options.begin(); x != parsed.options.end(); ++x)
            {
                if (x->unregistered)
                {
                    cout << "Warning: Unknown option: " << x->string_key 
                         << endl;
                }
            }

            // Return the vector of filename(s) given as positional options
            if (m_cmdLineOptions.count("input-file"))
            {
                return m_cmdLineOptions["input-file"].as<
                    std::vector<std::string> >();
            }
            else
            {
                return std::vector<std::string>();
            }
        }


        /**
         *
         */
        std::string SessionReader::ParseSessionName(
                std::vector<std::string> &filenames)
        {
            ASSERTL0(!filenames.empty(),
                     "At least one filename expected.");

            std::string retval = "";

            // First input file defines the session name
            std::string fname = filenames[0];

            // If loading a pre-partitioned mesh, remove _xml extension
            if (fname.size() > 4 &&
                fname.substr(fname.size() - 4, 4) == "_xml")
            {
                retval = fname.substr(0, fname.find_last_of("_"));
            }
            // otherwise remove the .xml extension
            else if (fname.size() > 4 &&
                fname.substr(fname.size() - 4, 4) == ".xml")
            {
                retval = fname.substr(0, fname.find_last_of("."));
            }
            // If compressed .xml.gz, remove both extensions
            else if (fname.size() > 7 &&
                fname.substr(fname.size() - 7, 7) == ".xml.gz")
            {
                retval = fname.substr(0, fname.find_last_of("."));
                retval = retval.substr(0, retval.find_last_of("."));
            }

            return retval;
        }


        /**
         *
         */
        TiXmlDocument& SessionReader::GetDocument()
        {
            ASSERTL1(m_xmlDoc, "XML Document not defined.");
            return *m_xmlDoc;
        }


        /**
         * The single parameter specifies a path to the requested element in a
         * similar format to the filesystem path. Given the following XML:
         * @code
         * <NEKTAR>
         *   <CONDITIONS>
         *     <PARAMETERS>
         *     ...
         *     </PARAMETERS>
         *   </CONDITIONS>
         * </NEKTAR>
         * @endcode
         * the PARAMETERS element would be retrieved by requesting the path:
         * @code
         * Nektar/Conditions/Parameters
         * @endcode
         * @note Paths are case-insensitive.
         *
         * @param   pPath       Path to requested element.
         *
         * @return Direct pointer to requested XML Element.
         */
        TiXmlElement* SessionReader::GetElement(const string& pPath)
        {
            std::string vPath = boost::to_upper_copy(pPath);
            std::vector<std::string> st;
            boost::split(st, vPath, boost::is_any_of("\\/ "));
            ASSERTL0(st.size() > 0, "No path given in XML element request.");

            TiXmlElement* vReturn = m_xmlDoc->FirstChildElement(st[0].c_str());
            ASSERTL0(vReturn, std::string("Cannot find element '")
                              + st[0] + std::string("'."));
            for (int i = 1; i < st.size(); ++i)
            {
                vReturn = vReturn->FirstChildElement(st[i].c_str());
                ASSERTL0(vReturn, std::string("Cannot find element '")
                                  + st[i] + std::string("'."));
            }
            return vReturn;
        }


        /**
         *
         */
        bool SessionReader::DefinesElement(const std::string &pPath) const
        {
            std::string vPath = boost::to_upper_copy(pPath);
            std::vector<std::string> st;
            boost::split(st, vPath, boost::is_any_of("\\/ "));
            ASSERTL0(st.size() > 0, "No path given in XML element request.");

            TiXmlElement* vReturn = m_xmlDoc->FirstChildElement(st[0].c_str());
            ASSERTL0(vReturn, std::string("Cannot find element '")
                              + st[0] + std::string("'."));
            for (int i = 1; i < st.size(); ++i)
            {
                vReturn = vReturn->FirstChildElement(st[i].c_str());
                if (!vReturn) return false;
            }
            return true;
        }


        /**
         *
         */
        const std::vector<std::string>& SessionReader::GetFilenames() const
        {
            return m_filenames;
        }


        /**
         *
         */
        const std::string& SessionReader::GetSessionName() const
        {
            return m_sessionName;
        }


        /**
         * Output is of the form [sessionName]_P[idx] where idx is the rank
         * of the process.
         */
        const std::string SessionReader::GetSessionNameRank() const
        {
            std::string  dirname = m_sessionName + "_xml"; 
            fs::path     pdirname(dirname);
            
            std::string vFilename = "P" + boost::lexical_cast<std::string>(m_comm->GetRowComm()->GetRank());
            fs::path    pFilename(vFilename);            

            fs::path fullpath = pdirname / pFilename;

            return PortablePath(fullpath);
        }

        /**
         *
         */
        CommSharedPtr& SessionReader::GetComm()
        {
            return m_comm;
        }


        /**
         * This routine finalises any parallel communication.
         *
         * @note This routine should be called at the very end of a users
         * application.
         */
        void SessionReader::Finalise()
        {
            m_comm->Finalise();
        }


        /**
         *
         */
        bool SessionReader::DefinesParameter(const std::string& pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ParameterMap::const_iterator paramIter = m_parameters.find(vName);
            return (paramIter != m_parameters.end());
        }


        /**
         * If the parameter is not defined, termination occurs. Therefore, the
         * parameters existence should be tested for using #DefinesParameter
         * before calling this function.
         *
         * @param   pName       The name of a floating-point parameter.
         * @returns The value of the floating-point parameter.
         */
        const NekDouble& SessionReader::GetParameter(
            const std::string& pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ParameterMap::const_iterator paramIter = m_parameters.find(vName);

            ASSERTL0(paramIter != m_parameters.end(),
                     "Unable to find requested parameter: " + pName);

            return paramIter->second;
        }


        /**
         *
         */
        void SessionReader::LoadParameter(
            const std::string &pName, int &pVar) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ParameterMap::const_iterator paramIter = m_parameters.find(vName);
            ASSERTL0(paramIter != m_parameters.end(), "Required parameter '" + 
                     pName + "' not specified in session.");
            pVar = (int)floor(paramIter->second);
        }


        /**
         *
         */
        void SessionReader::LoadParameter(
            const std::string &pName, int &pVar, const int &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ParameterMap::const_iterator paramIter = m_parameters.find(vName);
            if(paramIter != m_parameters.end())
            {
                pVar = (int)floor(paramIter->second);
            }
            else
            {
                pVar  = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::LoadParameter(
            const std::string &pName, NekDouble& pVar) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ParameterMap::const_iterator paramIter = m_parameters.find(vName);
            ASSERTL0(paramIter != m_parameters.end(), "Required parameter '" + 
                     pName + "' not specified in session.");
            pVar = paramIter->second;
        }


        /**
         *
         */
        void SessionReader::LoadParameter(
            const std::string &pName, 
                  NekDouble   &pVar, 
            const NekDouble   &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            ParameterMap::const_iterator paramIter = m_parameters.find(vName);
            if(paramIter != m_parameters.end())
            {
                pVar = paramIter->second;
            }
            else
            {
                pVar = pDefault;
            }
        }



        /**
         *
         */
        void SessionReader::SetParameter(const std::string &pName, int &pVar) 
        {
            std::string vName = boost::to_upper_copy(pName);
            m_parameters[vName] = pVar;
        }


        /**
         *
         */
        void SessionReader::SetParameter(
            const std::string &pName, NekDouble& pVar) 
        {
            std::string vName = boost::to_upper_copy(pName);
            m_parameters[vName] = pVar;
        }



        /**
         *
         */
        bool SessionReader::DefinesSolverInfo(const std::string &pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            SolverInfoMap::const_iterator infoIter = m_solverInfo.find(vName);
            return (infoIter != m_solverInfo.end());
        }


        /**
         *
         */
        const std::string& SessionReader::GetSolverInfo(
            const std::string &pProperty) const
        {
            std::string vProperty = boost::to_upper_copy(pProperty);
            SolverInfoMap::const_iterator iter = m_solverInfo.find(vProperty);

            ASSERTL1(iter != m_solverInfo.end(),
                     "Unable to find requested property: " + pProperty);

            return iter->second;
        }

        /**
         *
         */
        void SessionReader::SetSolverInfo(
            const std::string &pProperty, const std::string &pValue) 
        {
            std::string vProperty = boost::to_upper_copy(pProperty);
            SolverInfoMap::iterator iter = m_solverInfo.find(vProperty);

            ASSERTL1(iter != m_solverInfo.end(),
                     "Unable to find requested property: " + pProperty);

            iter->second = pValue;
        }

        /**
         *
         */
        void SessionReader::LoadSolverInfo(
            const std::string &pName, 
                  std::string &pVar, 
            const std::string &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            SolverInfoMap::const_iterator infoIter = m_solverInfo.find(vName);
            if(infoIter != m_solverInfo.end())
            {
                pVar = infoIter->second;
            }
            else
            {
                pVar = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::MatchSolverInfo(
            const std::string &pName,
            const std::string &pTrueVal,
                  bool        &pVar,
            const bool        &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            SolverInfoMap::const_iterator infoIter = m_solverInfo.find(vName);
            if(infoIter != m_solverInfo.end())
            {
                pVar = boost::iequals(infoIter->second, pTrueVal);
            }
            else
            {
                pVar = pDefault;
            }
        }


        /**
         *
         */
        bool SessionReader::MatchSolverInfo(
            const std::string &pName,
            const std::string &pTrueVal) const
        {
            if (DefinesSolverInfo(pName))
            {
                std::string vName = boost::to_upper_copy(pName);
                SolverInfoMap::const_iterator iter = m_solverInfo.find(vName);
                if(iter != m_solverInfo.end())
                {
                    return true;
                }
            }
            return false;
        }


        /**
         *
         */
        bool SessionReader::DefinesGlobalSysSolnInfo(const std::string &pVariable, 
                                                     const std::string &pProperty) const
        {

            GloSysSolnInfoList::const_iterator iter =
                    GetGloSysSolnList().find(pVariable);
            if(iter == GetGloSysSolnList().end())
            {
                return false;
            }

            std::string vProperty = boost::to_upper_copy(pProperty);
            
            GloSysInfoMap::const_iterator iter1 = iter->second.find(vProperty);
            if(iter1 == iter->second.end())
            {
                return false;
            }
            
            return true;
        }

        
        /**
         *
         */
        const std::string &SessionReader::GetGlobalSysSolnInfo(const std::string &pVariable, const std::string &pProperty) const
        {
            GloSysSolnInfoList::const_iterator iter; 

            ASSERTL0( (iter = GetGloSysSolnList().find(pVariable)) !=
                              GetGloSysSolnList().end(),
                      "Failed to find variable in GlobalSysSolnInfoList");

            std::string vProperty = boost::to_upper_copy(pProperty);
            GloSysInfoMap::const_iterator iter1; 

            ASSERTL0( (iter1 = iter->second.find(vProperty)) != iter->second.end(),
                      "Failed to find property: " + vProperty + " in GlobalSysSolnInfoList");
            
            return iter1->second;
        }
        
        /**
         *
         */
        bool SessionReader::DefinesGeometricInfo(const std::string &pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            GeometricInfoMap::const_iterator iter = m_geometricInfo.find(vName);
            return (iter != m_geometricInfo.end());
        }


        /**
         *
         */
        void SessionReader::LoadGeometricInfo(
            const std::string &pName,
                  std::string &pVar,
            const std::string &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            GeometricInfoMap::const_iterator iter = m_geometricInfo.find(vName);
            if(iter != m_geometricInfo.end())
            {
                pVar = iter->second;
            }
            else
            {
                pVar = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::LoadGeometricInfo(
            const std::string &pName,
                  bool        &pVar,
            const bool        &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            GeometricInfoMap::const_iterator iter = m_geometricInfo.find(vName);
            if(iter != m_geometricInfo.end())
            {
                if (iter->second == "TRUE")
                {
                    pVar = true;
                }
                else
                {
                    pVar = false;
                }
            }
            else
            {
                pVar = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::LoadGeometricInfo(
            const std::string &pName,
                  NekDouble   &pVar,
            const NekDouble   &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            GeometricInfoMap::const_iterator iter = m_geometricInfo.find(vName);
            if(iter != m_geometricInfo.end())
            {
                pVar = std::atoi(iter->second.c_str());
            }
            else
            {
                pVar = pDefault;
            }
        }


        /**
         *
         */
        void SessionReader::MatchGeometricInfo(
            const std::string &pName,
            const std::string &pTrueVal,
                  bool        &pVar,
            const bool        &pDefault) const
        {
            std::string vName = boost::to_upper_copy(pName);
            GeometricInfoMap::const_iterator iter = m_geometricInfo.find(vName);
            if(iter != m_geometricInfo.end())
            {
                pVar = boost::iequals(iter->second, pTrueVal);
            }
            else
            {
                pVar  = pDefault;
            }
        }


        /**
         *
         */
        const std::string& SessionReader::GetVariable(
            const unsigned int &idx) const
        {
            ASSERTL0(idx < m_variables.size(), "Variable index out of range.");
            return m_variables[idx];
        }



        /**
         *
         */
        void SessionReader::SetVariable(const unsigned int &idx, 
                                        std::string newname) 
        {
            ASSERTL0(idx < m_variables.size(), "Variable index out of range.");
            m_variables[idx] = newname;
        }


        /**
         *
         */
        std::vector<std::string> SessionReader::GetVariables() const
        {
            return m_variables;
        }


        /**
         *
         */
        bool SessionReader::DefinesFunction(const std::string &pName) const
        {
            FunctionMap::const_iterator it1;
            std::string vName = boost::to_upper_copy(pName);

            if ((it1 = m_functions.find(vName)) != m_functions.end())
            {
                return true;
            }
            return false;
        }


        /**
         *
         */
        bool SessionReader::DefinesFunction(
            const std::string &pName,
            const std::string &pVariable,
            const int pDomain) const
        {
            FunctionMap::const_iterator it1;
            FunctionVariableMap::const_iterator it2;
            std::string vName = boost::to_upper_copy(pName);

            // Check function exists
            if ((it1 = m_functions.find(vName))     != m_functions.end())
            {
                pair<std::string, int> key(pVariable,pDomain);
                pair<std::string, int> defkey("*",pDomain);
                bool varExists =
                    (it2 = it1->second.find(key)) != it1->second.end() ||
                    (it2 = it1->second.find(defkey)) != it1->second.end();
                return varExists;
            }
            return false;
        }


        /**
         *
         */
        EquationSharedPtr SessionReader::GetFunction(
            const std::string &pName,
            const std::string &pVariable,
            const int pDomain) const
        {
            FunctionMap::const_iterator it1;
            FunctionVariableMap::const_iterator it2, it3;
            std::string vName = boost::to_upper_copy(pName);

            ASSERTL0((it1 = m_functions.find(vName)) != m_functions.end(),
                     std::string("No such function '") + pName
                     + std::string("' has been defined in the session file."));

            // Check for specific and wildcard definitions
            pair<std::string,int> key(pVariable,pDomain);
            pair<std::string,int> defkey("*",pDomain);
            bool specific = (it2 = it1->second.find(key)) !=
                it1->second.end();
            bool wildcard = (it3 = it1->second.find(defkey)) !=
                it1->second.end();

            // Check function is defined somewhere
            ASSERTL0(specific || wildcard,
                     "No such variable " + pVariable
                     + " in domain " + boost::lexical_cast<string>(pDomain) 
                     + " defined for function " + pName
                     + " in session file.");

            // If not specific, must be wildcard
            if (!specific)
            {
                it2 = it3;
            }

            ASSERTL0((it2->second.m_type == eFunctionTypeExpression),
                    std::string("Function is defined by a file."));
            return it2->second.m_expression;
        }


        /**
         *
         */
        EquationSharedPtr SessionReader::GetFunction(
            const std::string  &pName,
            const unsigned int &pVar,
            const int pDomain) const
        {
            ASSERTL0(pVar < m_variables.size(), "Variable index out of range.");
            return GetFunction(pName, m_variables[pVar],pDomain);
        }


        /**
         *
         */
        enum FunctionType SessionReader::GetFunctionType(
            const std::string &pName,
            const std::string &pVariable,
            const int pDomain) const
        {
            FunctionMap::const_iterator it1;
            FunctionVariableMap::const_iterator it2, it3;
            std::string vName = boost::to_upper_copy(pName);

            it1 = m_functions.find(vName);
            ASSERTL0 (it1 != m_functions.end(),
                      std::string("Function '") + pName
                      + std::string("' not found."));

            // Check for specific and wildcard definitions
            pair<std::string,int> key(pVariable,pDomain);
            pair<std::string,int> defkey("*",pDomain);
            bool specific = (it2 = it1->second.find(key)) !=
                            it1->second.end();
            bool wildcard = (it3 = it1->second.find(defkey)) !=
                            it1->second.end();

            // Check function is defined somewhere
            ASSERTL0(specific || wildcard,
                     "No such variable " + pVariable
                     + " in domain " + boost::lexical_cast<string>(pDomain) 
                     + " defined for function " + pName
                     + " in session file.");

            // If not specific, must be wildcard
            if (!specific)
            {
                it2 = it3;
            }

            return it2->second.m_type;
        }


        /**
         *
         */
        enum FunctionType SessionReader::GetFunctionType(
            const std::string  &pName,
            const unsigned int &pVar,
            const int pDomain) const
        {
            ASSERTL0(pVar < m_variables.size(), "Variable index out of range.");
            return GetFunctionType(pName, m_variables[pVar],pDomain);
        }


        /**
         *
         */
        std::string SessionReader::GetFunctionFilename(
            const std::string &pName, 
            const std::string &pVariable,
            const int pDomain) const
        {
            FunctionMap::const_iterator it1;
            FunctionVariableMap::const_iterator it2, it3;
            std::string vName = boost::to_upper_copy(pName);

            it1 = m_functions.find(vName);
            ASSERTL0 (it1 != m_functions.end(),
                      std::string("Function '") + pName
                      + std::string("' not found."));

            // Check for specific and wildcard definitions
            pair<std::string,int> key(pVariable,pDomain);
            pair<std::string,int> defkey("*",pDomain);
            bool specific = (it2 = it1->second.find(key)) !=
                            it1->second.end();
            bool wildcard = (it3 = it1->second.find(defkey)) !=
                            it1->second.end();

            // Check function is defined somewhere
            ASSERTL0(specific || wildcard,
                     "No such variable " + pVariable
                     + " in domain " + boost::lexical_cast<string>(pDomain) 
                     + " defined for function " + pName
                     + " in session file.");
            
            // If not specific, must be wildcard
            if (!specific)
            {
                it2 = it3;
            }

            return it2->second.m_filename;
        }


        /**
         *
         */
        std::string SessionReader::GetFunctionFilename(
            const std::string  &pName, 
            const unsigned int &pVar,
            const int pDomain) const
        {
            ASSERTL0(pVar < m_variables.size(), "Variable index out of range.");
            return GetFunctionFilename(pName, m_variables[pVar],pDomain);
        }


        /**
         *
         */
        std::string SessionReader::GetFunctionFilenameVariable(
            const std::string &pName,
            const std::string &pVariable,
            const int pDomain) const
        {
            FunctionMap::const_iterator it1;
            FunctionVariableMap::const_iterator it2, it3;
            std::string vName = boost::to_upper_copy(pName);

            it1 = m_functions.find(vName);
            ASSERTL0 (it1 != m_functions.end(),
                      std::string("Function '") + pName
                      + std::string("' not found."));

            // Check for specific and wildcard definitions
            pair<std::string,int> key(pVariable,pDomain);
            pair<std::string,int> defkey("*",pDomain);
            bool specific = (it2 = it1->second.find(key)) !=
                            it1->second.end();
            bool wildcard = (it3 = it1->second.find(defkey)) !=
                            it1->second.end();

            // Check function is defined somewhere
            ASSERTL0(specific || wildcard,
                     "No such variable " + pVariable
                     + " in domain " + boost::lexical_cast<string>(pDomain)
                     + " defined for function " + pName
                     + " in session file.");

            // If not specific, must be wildcard
            if (!specific)
            {
                it2 = it3;
            }

            return it2->second.m_fileVariable;
        }


        /**
         *
         */
        AnalyticExpressionEvaluator& SessionReader::GetExpressionEvaluator()
        {
            return m_exprEvaluator;
        }


        /**
         *
         */
        bool SessionReader::DefinesTag(const std::string &pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            TagMap::const_iterator vTagIterator = m_tags.find(vName);
            return (vTagIterator != m_tags.end());
        }


        /**
         *
         */
        void SessionReader::SetTag(
            const std::string &pName, 
            const std::string &pValue)
        {
            std::string vName = boost::to_upper_copy(pName);
            m_tags[vName] = pValue;
        }


        /**
         *
         */
        const std::string &SessionReader::GetTag(const std::string& pName) const
        {
            std::string vName = boost::to_upper_copy(pName);
            TagMap::const_iterator vTagIterator = m_tags.find(vName);
            ASSERTL0(vTagIterator != m_tags.end(),
                     "Requested tag does not exist.");
            return vTagIterator->second;
        }


        /**
         *
         */
        const FilterMap &SessionReader::GetFilters() const
        {
            return m_filters;
        }


        /**
         *
         */
        bool SessionReader::DefinesCmdLineArgument(
            const std::string& pName) const
        {
            return (m_cmdLineOptions.find(pName) != m_cmdLineOptions.end());
        }


        /**
         *
         */
        void SessionReader::SubstituteExpressions(std::string& pExpr)
        {
            ExpressionMap::iterator exprIter;
            for (exprIter  = m_expressions.begin(); 
                 exprIter != m_expressions.end(); ++exprIter)
            {
                boost::replace_all(pExpr, exprIter->first, exprIter->second);
            }
        }

        CompositeOrdering SessionReader::GetCompositeOrdering() const
        {
            return m_compOrder;
        }

        BndRegionOrdering SessionReader::GetBndRegionOrdering() const
        {
            return m_bndRegOrder;
        }

        /**
         *
         */
        void SessionReader::LoadDoc(
            const std::string &pFilename,
            TiXmlDocument* pDoc) const
        {
            if (pFilename.size() > 3 &&
                pFilename.substr(pFilename.size() - 3, 3) == ".gz")
            {
                ifstream file(pFilename.c_str(),
                              ios_base::in | ios_base::binary);
                ASSERTL0(file.good(), "Unable to open file: " + pFilename);
                stringstream ss;
                io::filtering_streambuf<io::input> in;
                in.push(io::gzip_decompressor());
                in.push(file);
                try
                {
                    io::copy(in, ss);
                    ss >> (*pDoc);
                }
                catch (io::gzip_error& e)
                {
                    ASSERTL0(false,
                             "Error: File '" + pFilename + "' is corrupt.");
                }
            }
            else if (pFilename.size() > 4 &&
                    pFilename.substr(pFilename.size() - 4, 4) == "_xml")
            {
                fs::path    pdirname(pFilename);
                boost::format pad("P%1$07d.xml");
                pad % m_comm->GetRank();
                fs::path    pRankFilename(pad.str());
                fs::path fullpath = pdirname / pRankFilename;

                ifstream file(PortablePath(fullpath).c_str());
                ASSERTL0(file.good(), "Unable to open file: " + fullpath.string());
                file >> (*pDoc);
            }
            else
            {
                ifstream file(pFilename.c_str());
                ASSERTL0(file.good(), "Unable to open file: " + pFilename);
                file >> (*pDoc);
            }
        }

        /**
         *
         */
        TiXmlDocument *SessionReader::MergeDoc(
            const std::vector<std::string> &pFilenames) const
        {
            ASSERTL0(pFilenames.size() > 0, "No filenames for merging.");

            // Read the first document
            TiXmlDocument *vMainDoc = new TiXmlDocument;
            LoadDoc(pFilenames[0], vMainDoc);

            TiXmlHandle vMainHandle(vMainDoc);
            TiXmlElement* vMainNektar = 
                vMainHandle.FirstChildElement("NEKTAR").Element();

            // Read all subsequent XML documents.
            // For each element within the NEKTAR tag, use it to replace the
            // version already present in the loaded XML data.
            for (int i = 1; i < pFilenames.size(); ++i)
            {
                if((pFilenames[i].compare(pFilenames[i].size()-3,3,"xml") == 0)
                   ||(pFilenames[i].compare(pFilenames[i].size()-6,6,"xml.gz") == 0))
                {
                    TiXmlDocument* vTempDoc = new TiXmlDocument;
                    LoadDoc(pFilenames[i], vTempDoc);
                    
                    TiXmlHandle docHandle(vTempDoc);
                    TiXmlElement* vTempNektar;
                    vTempNektar = docHandle.FirstChildElement("NEKTAR").Element();
                    ASSERTL0(vTempNektar, "Unable to find NEKTAR tag in file.");
                    TiXmlElement* p = vTempNektar->FirstChildElement();
                    
                    while (p)
                    {
                        TiXmlElement *vMainEntry = 
                            vMainNektar->FirstChildElement(p->Value());
                        TiXmlElement *q = new TiXmlElement(*p);
                        if (vMainEntry)
                        {
                            vMainNektar->RemoveChild(vMainEntry);
                        }
                        vMainNektar->LinkEndChild(q);
                        p = p->NextSiblingElement();
                    }
                    
                    delete vTempDoc;
                }
            }
            return vMainDoc;
        }


        /**
         *
         */
        void SessionReader::ParseDocument()
        {
            // Check we actually have a document loaded.
            ASSERTL0(m_xmlDoc, "No XML document loaded.");

            // Look for all data in CONDITIONS block.
            TiXmlHandle docHandle(m_xmlDoc);
            TiXmlElement* e;
            e = docHandle.FirstChildElement("NEKTAR").
                FirstChildElement("CONDITIONS").Element();

            // Read the various sections of the CONDITIONS block
            ReadParameters        (e);
            ReadSolverInfo        (e);
            ReadGlobalSysSolnInfo (e);
            ReadExpressions       (e);
            ReadVariables         (e);
            ReadFunctions         (e);

            e = docHandle.FirstChildElement("NEKTAR").
                FirstChildElement("FILTERS").Element();

            ReadFilters(e);
        }


        /**
         *
         */
        void SessionReader::CreateComm(
            int               &argc, 
            char*              argv[])
        {
            if (argc == 0)
            {
                m_comm = GetCommFactory().CreateInstance("Serial", 0, 0);
            }
            else
            {
                string vCommModule("Serial");
                if (GetCommFactory().ModuleExists("ParallelMPI"))
                {
                    vCommModule = "ParallelMPI";
                }

                m_comm = GetCommFactory().CreateInstance(vCommModule,argc,argv);
            }
        }


        /**
         *
         */
        void SessionReader::PartitionMesh()
        {
            ASSERTL0(m_comm.get(), "Communication not initialised.");

            // Get row of comm, or the whole comm if not split
            CommSharedPtr vCommMesh = m_comm->GetRowComm();
            const bool isRoot = (m_comm->GetRank() == 0);

            // Delete any existing loaded mesh
            if (m_xmlDoc)
            {
                delete m_xmlDoc;
            }

            // Load file for root process only (since this is always needed)
            // and determine if the provided geometry has already been
            // partitioned. This will be the case if the user provides the
            // directory of mesh partitions as an input. Partitioned geometries
            // have the attribute
            //    PARTITION=X
            // where X is the number of the partition (and should match the
            // process rank). The result is shared with all other processes.
            int isPartitioned = 0;
            if (isRoot)
            {
                m_xmlDoc = MergeDoc(m_filenames);
                if (DefinesElement("Nektar/Geometry"))
                {
                    if (GetElement("Nektar/Geometry")->Attribute("PARTITION"))
                    {
                        cout << "Using pre-partitioned mesh." << endl;
                        isPartitioned = 1;
                    }
                }
            }
            GetComm()->AllReduce(isPartitioned, LibUtilities::ReduceMax);

            // If the mesh is already partitioned, we are done. Remaining
            // processes must load their partitions.
            if (isPartitioned) {
                if (!isRoot)
                {
                    m_xmlDoc = MergeDoc(m_filenames);
                }
                return;
            }

            // Default partitioner to use is Metis. Use Scotch as default
            // if it is installed. Override default with command-line flags
            // if they are set.
            string vPartitionerName = "Metis";
            if (GetMeshPartitionFactory().ModuleExists("Scotch"))
            {
                vPartitionerName = "Scotch";
            }
            if (DefinesCmdLineArgument("use-metis"))
            {
                vPartitionerName = "Metis";
            }
            if (DefinesCmdLineArgument("use-scotch"))
            {
                vPartitionerName = "Scotch";
            }

            // Mesh has not been partitioned so do partitioning if required.
            // Note in the serial case nothing is done as we have already loaded
            // the mesh.
            if (DefinesCmdLineArgument("part-only"))
            {
                // Perform partitioning of the mesh only. For this we insist
                // the code is run in serial (parallel execution is pointless).
                ASSERTL0(GetComm()->GetSize() == 1,
                        "The 'part-only' option should be used in serial.");

                // Number of partitions is specified by the parameter.
                int nParts = GetCmdLineArgument<int>("part-only");
                SessionReaderSharedPtr vSession     = GetSharedThisPtr();
                MeshPartitionSharedPtr vPartitioner = 
                        GetMeshPartitionFactory().CreateInstance(
                                            vPartitionerName, vSession);
                vPartitioner->PartitionMesh(nParts, true);
                vPartitioner->WriteAllPartitions(vSession);
                vPartitioner->GetCompositeOrdering(m_compOrder);
                vPartitioner->GetBndRegionOrdering(m_bndRegOrder);

                if (isRoot && DefinesCmdLineArgument("part-info"))
                {
                    vPartitioner->PrintPartInfo(std::cout);
                }

                Finalise();
                exit(0);
            }
            else if (vCommMesh->GetSize() > 1)
            {
                SessionReaderSharedPtr vSession     = GetSharedThisPtr();
                int nParts = vCommMesh->GetSize();
                if (DefinesCmdLineArgument("shared-filesystem"))
                {
                    CommSharedPtr vComm = GetComm();
                    vector<unsigned int> keys, vals;
                    int i;

                    if (vComm->GetRank() == 0)
                    {
                        m_xmlDoc = MergeDoc(m_filenames);

                        MeshPartitionSharedPtr vPartitioner =
                                GetMeshPartitionFactory().CreateInstance(
                                                    vPartitionerName, vSession);
                        vPartitioner->PartitionMesh(nParts, true);
                        vPartitioner->WriteAllPartitions(vSession);
                        vPartitioner->GetCompositeOrdering(m_compOrder);
                        vPartitioner->GetBndRegionOrdering(m_bndRegOrder);

                        // Communicate orderings to the other processors.

                        // First send sizes of the orderings and boundary
                        // regions to allocate storage on the remote end.
                        keys.resize(2);
                        keys[0] = m_compOrder.size();
                        keys[1] = m_bndRegOrder.size();

                        for (i = 1; i < vComm->GetSize(); ++i)
                        {
                            vComm->Send(i, keys);
                        }

                        // Construct the keys and sizes of values for composite
                        // ordering
                        CompositeOrdering::iterator cIt;
                        keys.resize(m_compOrder.size());
                        vals.resize(m_compOrder.size());

                        for (cIt  = m_compOrder.begin(), i = 0;
                             cIt != m_compOrder.end(); ++cIt, ++i)
                        {
                            keys[i] = cIt->first;
                            vals[i] = cIt->second.size();
                        }

                        // Send across data.
                        for (i = 1; i < vComm->GetSize(); ++i)
                        {
                            vComm->Send(i, keys);
                            vComm->Send(i, vals);

                            for (cIt  = m_compOrder.begin();
                                 cIt != m_compOrder.end(); ++cIt)
                            {
                                vComm->Send(i, cIt->second);
                            }
                        }

                        // Construct the keys and sizes of values for composite
                        // ordering
                        BndRegionOrdering::iterator bIt;
                        keys.resize(m_bndRegOrder.size());
                        vals.resize(m_bndRegOrder.size());

                        for (bIt  = m_bndRegOrder.begin(), i = 0;
                             bIt != m_bndRegOrder.end(); ++bIt, ++i)
                        {
                            keys[i] = bIt->first;
                            vals[i] = bIt->second.size();
                        }

                        // Send across data.
                        for (i = 1; i < vComm->GetSize(); ++i)
                        {
                            vComm->Send(i, keys);
                            vComm->Send(i, vals);

                            for (bIt  = m_bndRegOrder.begin();
                                 bIt != m_bndRegOrder.end(); ++bIt)
                            {
                                vComm->Send(i, bIt->second);
                            }
                        }

                        if (DefinesCmdLineArgument("part-info"))
                        {
                            vPartitioner->PrintPartInfo(std::cout);
                        }
                    }
                    else
                    {
                        keys.resize(2);
                        vComm->Recv(0, keys);

                        int cmpSize = keys[0];
                        int bndSize = keys[1];

                        keys.resize(cmpSize);
                        vals.resize(cmpSize);
                        vComm->Recv(0, keys);
                        vComm->Recv(0, vals);

                        for (int i = 0; i < keys.size(); ++i)
                        {
                            vector<unsigned int> tmp(vals[i]);
                            vComm->Recv(0, tmp);
                            m_compOrder[keys[i]] = tmp;
                        }

                        keys.resize(bndSize);
                        vals.resize(bndSize);
                        vComm->Recv(0, keys);
                        vComm->Recv(0, vals);

                        for (int i = 0; i < keys.size(); ++i)
                        {
                            vector<unsigned int> tmp(vals[i]);
                            vComm->Recv(0, tmp);
                            m_bndRegOrder[keys[i]] = tmp;
                        }
                    }
                }
                else
                {
                    // Need to load mesh on non-root processes.
                    if (!isRoot)
                    {
                        m_xmlDoc = MergeDoc(m_filenames);
                    }

                    // Partitioner now operates in parallel
                    // Each process receives partitioning over interconnect
                    // and writes its own session file to the working directory.
                    MeshPartitionSharedPtr vPartitioner =
                                GetMeshPartitionFactory().CreateInstance(
                                                    vPartitionerName, vSession);
                    vPartitioner->PartitionMesh(nParts, false);
                    vPartitioner->WriteLocalPartition(vSession);
                    vPartitioner->GetCompositeOrdering(m_compOrder);
                    vPartitioner->GetBndRegionOrdering(m_bndRegOrder);

                    if (DefinesCmdLineArgument("part-info") && isRoot)
                    {
                        vPartitioner->PrintPartInfo(std::cout);
                    }
                }
                m_comm->Block();

                std::string  dirname = GetSessionName() + "_xml";
                fs::path    pdirname(dirname);
                boost::format pad("P%1$07d.xml");
                pad % m_comm->GetRowComm()->GetRank();
                fs::path    pFilename(pad.str());
                fs::path fullpath = pdirname / pFilename;

                std::string vFilename = PortablePath(fullpath);

                if (m_xmlDoc)
                {
                    delete m_xmlDoc;
                }
                m_xmlDoc = new TiXmlDocument(vFilename);

                ASSERTL0(m_xmlDoc, "Failed to create XML document object.");

                bool loadOkay = m_xmlDoc->LoadFile(vFilename);
                ASSERTL0(loadOkay, "Unable to load file: " + vFilename       +
                         ". Check XML standards compliance. Error on line: " +
                         boost::lexical_cast<std::string>(m_xmlDoc->Row()));
            }
            else
            {
                m_xmlDoc = MergeDoc(m_filenames);
            }
        }


        /**
         * Splits the processes into a cartesian grid and creates communicators
         * for each row and column of the grid. The grid is defined by the
         * PROC_X parameter which, if specified, gives the number of processes
         * spanned by the Fourier direction. PROC_X must exactly divide the
         * total number of processes or an error is thrown.
         */
        void SessionReader::PartitionComm()
        {
            if (m_comm->GetSize() > 1)
            {
                int nProcZ  = 1;
                int nProcY  = 1;
                int nProcX  = 1;
                int nStripZ = 1;
                if (DefinesCmdLineArgument("npx")) {
                    nProcX = GetCmdLineArgument<int>("npx");
                }
                if (DefinesCmdLineArgument("npy")) {
                    nProcY = GetCmdLineArgument<int>("npy");
                }
                if (DefinesCmdLineArgument("npz")) {
                    nProcZ = GetCmdLineArgument<int>("npz");
                }
                if (DefinesCmdLineArgument("nsz")) {
                    nStripZ = GetCmdLineArgument<int>("nsz");
                }
                ASSERTL0(m_comm->GetSize() % (nProcZ*nProcY*nProcX) == 0,
                         "Cannot exactly partition using PROC_Z value.");
                ASSERTL0(nProcZ % nProcY == 0,
                         "Cannot exactly partition using PROC_Y value.");
                ASSERTL0(nProcY % nProcX == 0,
                         "Cannot exactly partition using PROC_X value.");

                // Number of processes associated with the spectral method
                int nProcSm  = nProcZ * nProcY * nProcX;

                // Number of processes associated with the spectral element
                // method.
                int nProcSem = m_comm->GetSize() / nProcSm;

                m_comm->SplitComm(nProcSm,nProcSem);
                m_comm->GetColumnComm()->SplitComm(nProcZ/nStripZ,nStripZ);
                m_comm->GetColumnComm()->GetColumnComm()->SplitComm(
                                            (nProcY*nProcX),nProcZ/nStripZ);
                m_comm->GetColumnComm()->GetColumnComm()->GetColumnComm()
                                            ->SplitComm(nProcX,nProcY);
            }
        }


        /**
         *
         */
        void SessionReader::ReadParameters(TiXmlElement *conditions)
        {
            m_parameters.clear();

            if (!conditions)
            {
                return;
            }

            TiXmlElement *parametersElement = conditions->FirstChildElement(
                "PARAMETERS");

            // See if we have parameters defined.  They are optional so we go on
            // if not.
            if (parametersElement)
            {
                TiXmlElement *parameter = 
                    parametersElement->FirstChildElement("P");

                ParameterMap caseSensitiveParameters;

                // Multiple nodes will only occur if there is a comment in
                // between definitions.
                while (parameter)
                {
                    stringstream tagcontent;
                    tagcontent << *parameter;
                    TiXmlNode *node = parameter->FirstChild();

                    while (node && node->Type() != TiXmlNode::TINYXML_TEXT)
                    {
                        node = node->NextSibling();
                    }

                    if (node)
                    {
                        // Format is "paramName = value"
                        std::string line = node->ToText()->Value(), lhs, rhs;

                        try {
                            ParseEquals(line, lhs, rhs);
                        }
                        catch (...)
                        {
                            ASSERTL0(false, "Syntax error in parameter "
                                     "expression '" + line 
                                     + "' in XML element: \n\t'"
                                     + tagcontent.str() + "'");
                        }

                        // We want the list of parameters to have their RHS
                        // evaluated, so we use the expression evaluator to do
                        // the dirty work.
                        if (!lhs.empty() && !rhs.empty())
                        {
                            NekDouble value=0.0;
                            try
                            {
                                LibUtilities::Equation expession(
                                    GetSharedThisPtr(), rhs);
                                value = expession.Evaluate();
                            }
                            catch (const std::runtime_error &)
                            {
                                ASSERTL0(false, 
                                         "Error evaluating parameter expression"
                                         " '" + rhs + "' in XML element: \n\t'"
                                         + tagcontent.str() + "'");
                            }
                            m_exprEvaluator.SetParameter(lhs, value);
                            caseSensitiveParameters[lhs] = value;
                            boost::to_upper(lhs);
                            m_parameters[lhs] = value;
                        }
                    }

                    parameter = parameter->NextSiblingElement();
                }
            }
        }


        /**
         *
         */
        void SessionReader::ReadSolverInfo(TiXmlElement *conditions)
        {
            m_solverInfo.clear();
            m_solverInfo = GetSolverInfoDefaults();

            if (!conditions)
            {
                return;
            }

            TiXmlElement *solverInfoElement = 
                conditions->FirstChildElement("SOLVERINFO");

            if (solverInfoElement)
            {
                TiXmlElement *solverInfo = 
                    solverInfoElement->FirstChildElement("I");

                while (solverInfo)
                {
                    std::stringstream tagcontent;
                    tagcontent << *solverInfo;
                    // read the property name
                    ASSERTL0(solverInfo->Attribute("PROPERTY"),
                             "Missing PROPERTY attribute in solver info "
                             "XML element: \n\t'" + tagcontent.str() + "'");
                    std::string solverProperty = 
                        solverInfo->Attribute("PROPERTY");
                    ASSERTL0(!solverProperty.empty(),
                             "PROPERTY attribute must be non-empty in XML "
                             "element: \n\t'" + tagcontent.str() + "'");

                    // make sure that solver property is capitalised
                    std::string solverPropertyUpper =
                        boost::to_upper_copy(solverProperty);

                    // read the value
                    ASSERTL0(solverInfo->Attribute("VALUE"),
                            "Missing VALUE attribute in solver info "
                            "XML element: \n\t'" + tagcontent.str() + "'");
                    std::string solverValue    = solverInfo->Attribute("VALUE");
                    ASSERTL0(!solverValue.empty(),
                             "VALUE attribute must be non-empty in XML "
                             "element: \n\t'" + tagcontent.str() + "'");

                    EnumMapList::const_iterator propIt = 
                        GetSolverInfoEnums().find(solverPropertyUpper);
                    if (propIt != GetSolverInfoEnums().end())
                    {
                        EnumMap::const_iterator valIt = 
                            propIt->second.find(solverValue);
                        ASSERTL0(valIt != propIt->second.end(),
                                 "Value '" + solverValue + "' is not valid for "
                                 "property '" + solverProperty + "'");
                    }

                    // Set Variable
                    m_solverInfo[solverPropertyUpper] = solverValue;
                    solverInfo = solverInfo->NextSiblingElement("I");
                }
            }
            
            if (m_comm && m_comm->GetRowComm()->GetSize() > 1)
            {
                ASSERTL0(
                    m_solverInfo["GLOBALSYSSOLN"] == "IterativeFull"       ||
                    m_solverInfo["GLOBALSYSSOLN"] == "IterativeStaticCond" ||
                    m_solverInfo["GLOBALSYSSOLN"] == 
                        "IterativeMultiLevelStaticCond"                    ||
                    m_solverInfo["GLOBALSYSSOLN"] == "XxtFull"             ||
                    m_solverInfo["GLOBALSYSSOLN"] == "XxtStaticCond"       ||
                    m_solverInfo["GLOBALSYSSOLN"] ==
                        "XxtMultiLevelStaticCond"                          ||
                    m_solverInfo["GLOBALSYSSOLN"] == "PETScFull"           ||
                    m_solverInfo["GLOBALSYSSOLN"] == "PETScStaticCond"     ||
                    m_solverInfo["GLOBALSYSSOLN"] ==
                        "PETScMultiLevelStaticCond",
                    "A parallel solver must be used when run in parallel.");
            }
        }



        /**
         *
         */
        void SessionReader::ReadGlobalSysSolnInfo(TiXmlElement *conditions)
        {
            GetGloSysSolnList().clear();

            if (!conditions)
            {
                return;
            }

            TiXmlElement *GlobalSys =
                            conditions->FirstChildElement("GLOBALSYSSOLNINFO");

            if(!GlobalSys)
            {
                return;
            }

            TiXmlElement *VarInfo   = GlobalSys->FirstChildElement("V");

            while (VarInfo)
            {
                std::stringstream tagcontent;
                tagcontent << *VarInfo;
                ASSERTL0(VarInfo->Attribute("VAR"),
                         "Missing VAR attribute in GobalSysSolnInfo XML "
                         "element: \n\t'" + tagcontent.str() + "'");

                std::string VarList = VarInfo->Attribute("VAR");
                ASSERTL0(!VarList.empty(),
                         "VAR attribute must be non-empty in XML element:\n\t'"
                         + tagcontent.str() + "'");

                // generate a list of variables.
                std::vector<std::string> varStrings;
                bool valid = ParseUtils::GenerateOrderedStringVector(
                                                VarList.c_str(),varStrings);

                ASSERTL0(valid,"Unable to process list of variable in XML "
                               "element \n\t'" + tagcontent.str() + "'");

                if(varStrings.size())
                {
                    TiXmlElement *SysSolnInfo = VarInfo->FirstChildElement("I");

                    while (SysSolnInfo)
                    {
                        tagcontent.clear();
                        tagcontent << *SysSolnInfo;
                        // read the property name
                        ASSERTL0(SysSolnInfo->Attribute("PROPERTY"),
                                 "Missing PROPERTY attribute in "
                                 "GlobalSysSolnInfo for variable(s) '"
                                 + VarList + "' in XML element: \n\t'" 
                                 + tagcontent.str() + "'");

                        std::string SysSolnProperty =
                            SysSolnInfo->Attribute("PROPERTY");

                        ASSERTL0(!SysSolnProperty.empty(),
                                 "GlobalSysSolnIno properties must have a "
                                 "non-empty name for variable(s) : '"
                                 + VarList + "' in XML element: \n\t'"
                                 + tagcontent.str() + "'");

                        // make sure that solver property is capitalised
                        std::string SysSolnPropertyUpper =
                                        boost::to_upper_copy(SysSolnProperty);

                        // read the value
                        ASSERTL0(SysSolnInfo->Attribute("VALUE"),
                                 "Missing VALUE attribute in GlobalSysSolnInfo "
                                 "for variable(s) '" + VarList
                                 + "' in XML element: \n\t"
                                 + tagcontent.str() + "'");

                        std::string SysSolnValue =
                                            SysSolnInfo->Attribute("VALUE");
                        ASSERTL0(!SysSolnValue.empty(),
                                 "GlobalSysSolnInfo properties must have a "
                                 "non-empty value for variable(s) '"
                                 + VarList + "' in XML element: \n\t'"
                                 + tagcontent.str() + "'");

                        // Store values under variable map.
                        for(int i = 0; i < varStrings.size(); ++i)
                        {
                            GloSysSolnInfoList::iterator x;
                            if ((x = GetGloSysSolnList().find(varStrings[i])) ==
                                    GetGloSysSolnList().end())
                            {
                                (GetGloSysSolnList()[varStrings[i]])[
                                        SysSolnPropertyUpper] = SysSolnValue;
                            }
                            else
                            {
                                x->second[SysSolnPropertyUpper] = SysSolnValue;
                            }
                        }

                        SysSolnInfo = SysSolnInfo->NextSiblingElement("I");
                    }
                    VarInfo = VarInfo->NextSiblingElement("V");
                }
            }

            if (m_verbose && GetGloSysSolnList().size() > 0 && m_comm)
            {
                if(m_comm->GetRank() == 0)
                {
                    cout << "GlobalSysSoln Info:" << endl;

                    GloSysSolnInfoList::iterator x;
                    for (x = GetGloSysSolnList().begin();
                         x != GetGloSysSolnList().end();
                         ++x)
                    {
                        cout << "\t Variable: " << x->first <<  endl;

                        GloSysInfoMap::iterator y;
                        for (y = x->second.begin(); y != x->second.end(); ++y)
                        {
                            cout << "\t\t " << y->first  << " = " << y->second
                                 << endl;
                        }
                    }
                    cout << endl;
                }
            }
        }


        /**
         *
         */
        void SessionReader::ReadExpressions(TiXmlElement *conditions)
        {
            m_expressions.clear();

            if (!conditions)
            {
                return;
            }

            TiXmlElement *expressionsElement = 
                conditions->FirstChildElement("EXPRESSIONS");

            if (expressionsElement)
            {
                TiXmlElement *expr = expressionsElement->FirstChildElement("E");

                while (expr)
                {
                    stringstream tagcontent;
                    tagcontent << *expr;
                    ASSERTL0(expr->Attribute("NAME"),
                             "Missing NAME attribute in expression "
                             "definition: \n\t'" + tagcontent.str() + "'");
                    std::string nameString = expr->Attribute("NAME");
                    ASSERTL0(!nameString.empty(),
                             "Expressions must have a non-empty name: \n\t'"
                             + tagcontent.str() + "'");

                    ASSERTL0(expr->Attribute("VALUE"),
                             "Missing VALUE attribute in expression "
                             "definition: \n\t'" + tagcontent.str() + "'");
                    std::string valString = expr->Attribute("VALUE");
                    ASSERTL0(!valString.empty(),
                             "Expressions must have a non-empty value: \n\t'"
                            + tagcontent.str() + "'");

                    ExpressionMap::iterator exprIter
                                            = m_expressions.find(nameString);
                    ASSERTL0(exprIter == m_expressions.end(),
                             std::string("Expression '") + nameString
                             + std::string("' already specified."));

                    m_expressions[nameString] = valString;
                    expr = expr->NextSiblingElement("E");
                }
            }
        }


        /**
         *
         */
        void SessionReader::ReadVariables(TiXmlElement *conditions)
        {
            m_variables.clear();

            if (!conditions)
            {
                return;
            }

            TiXmlElement *variablesElement = 
                conditions->FirstChildElement("VARIABLES");

            // See if we have parameters defined. They are optional so we go on
            // if not.
            if (variablesElement)
            {
                TiXmlElement *varElement = 
                    variablesElement->FirstChildElement("V");

                // Sequential counter for the composite numbers.
                int nextVariableNumber = -1;

                while (varElement)
                {
                    stringstream tagcontent;
                    tagcontent << *varElement;

                    /// All elements are of the form: "<V ID="#"> name = value
                    /// </V>", with ? being the element type.
                    nextVariableNumber++;

                    int i;
                    int err = varElement->QueryIntAttribute("ID", &i);
                    ASSERTL0(err == TIXML_SUCCESS,
                             "Variables must have a unique ID number attribute "
                             "in XML element: \n\t'" + tagcontent.str() + "'");
                    ASSERTL0(i == nextVariableNumber,
                             "ID numbers for variables must begin with zero and"
                             " be sequential in XML element: \n\t'"
                             + tagcontent.str() + "'");

                    TiXmlNode* varChild = varElement->FirstChild();
                    // This is primarily to skip comments that may be present.
                    // Comments appear as nodes just like elements.  We are
                    // specifically looking for text in the body of the
                    // definition.
                    while(varChild && varChild->Type() != TiXmlNode::TINYXML_TEXT)
                    {
                        varChild = varChild->NextSibling();
                    }

                    ASSERTL0(varChild,
                             "Unable to read variable definition body for "
                             "variable with ID "
                             + boost::lexical_cast<string>(i) 
                             + " in XML element: \n\t'"
                             + tagcontent.str() + "'");
                    std::string variableName = varChild->ToText()->ValueStr();

                    std::istringstream variableStrm(variableName);
                    variableStrm >> variableName;

                    ASSERTL0(std::find(m_variables.begin(), m_variables.end(), 
                                       variableName) == m_variables.end(),
                             "Variable with ID "
                             + boost::lexical_cast<string>(i) 
                             + " in XML element \n\t'" + tagcontent.str()
                             + "'\nhas already been defined.");

                    m_variables.push_back(variableName);

                    varElement = varElement->NextSiblingElement("V");
                }

                ASSERTL0(nextVariableNumber > -1,
                         "Number of variables must be greater than zero.");
            }
        }


        /**
         *
         */
        void SessionReader::ReadFunctions(TiXmlElement *conditions)
        {
            m_functions.clear();

            if (!conditions)
            {
                return;
            }

            // Scan through conditions section looking for functions.
            TiXmlElement *function = conditions->FirstChildElement("FUNCTION");
            while (function)
            {
                stringstream tagcontent;
                tagcontent << *function;

                // Every function must have a NAME attribute
                ASSERTL0(function->Attribute("NAME"),
                         "Functions must have a NAME attribute defined in XML "
                         "element: \n\t'" + tagcontent.str() + "'");
                std::string functionStr = function->Attribute("NAME");
                ASSERTL0(!functionStr.empty(),
                         "Functions must have a non-empty name in XML "
                         "element: \n\t'" + tagcontent.str() + "'");

                // Store function names in uppercase to remain case-insensitive.
                boost::to_upper(functionStr);

                // Retrieve first entry (variable, or file)
                TiXmlElement *variable  = function->FirstChildElement();

                // Create new function structure with default type of none.
                FunctionVariableMap functionVarMap;

                // Process all entries in the function block
                while (variable)
                {
                    FunctionVariableDefinition funcDef;
                    std::string conditionType = variable->Value();

                    // If no var is specified, assume wildcard
                    std::string variableStr;
                    if (!variable->Attribute("VAR"))
                    {
                        variableStr = "*";
                    }
                    else
                    {
                        variableStr = variable->Attribute("VAR");
                    }

                    // Parse list of variables
                    std::vector<std::string> variableList;
                    ParseUtils::GenerateOrderedStringVector(variableStr.c_str(),
                                                            variableList);

                    // If no domain string put to 0
                    std::string domainStr;
                    if (!variable->Attribute("DOMAIN"))
                    {
                        domainStr = "0";
                    }
                    else
                    {
                        domainStr = variable->Attribute("DOMAIN");
                    }

                    // Parse list of variables
                    std::vector<std::string> varSplit;
                    std::vector<unsigned int> domainList;
                    ParseUtils::GenerateSeqVector(domainStr.c_str(), domainList);

                    // Expressions are denoted by E
                    if (conditionType == "E")
                    {
                        funcDef.m_type = eFunctionTypeExpression;

                        // Expression must have a VALUE.
                        ASSERTL0(variable->Attribute("VALUE"),
                                 "Attribute VALUE expected for function '"
                                 + functionStr + "'.");
                        std::string fcnStr = variable->Attribute("VALUE");

                        ASSERTL0(!fcnStr.empty(),
                                 (std::string("Expression for var: ")
                                 + variableStr
                                 + std::string(" must be specified.")).c_str());

                        SubstituteExpressions(fcnStr);

                        // set expression
                        funcDef.m_expression = MemoryManager<Equation>
                            ::AllocateSharedPtr(GetSharedThisPtr(),fcnStr);
                    }

                    // Files are denoted by F
                    else if (conditionType == "F")
                    {
                        if (variable->Attribute("TIMEDEPENDENT") &&
                            boost::lexical_cast<bool>(variable->Attribute("TIMEDEPENDENT")))
                        {
                            funcDef.m_type = eFunctionTypeTransientFile;
                        }
                        else
                        {
                            funcDef.m_type = eFunctionTypeFile;
                        }

                        // File must have a FILE.
                        ASSERTL0(variable->Attribute("FILE"),
                                 "Attribute FILE expected for function '"
                                 + functionStr + "'.");
                        std::string filenameStr = variable->Attribute("FILE");

                        ASSERTL0(!filenameStr.empty(),
                                 "A filename must be specified for the FILE "
                                 "attribute of function '" + functionStr
                                 + "'.");

                        std::vector<std::string> fSplit;
                        boost::split(fSplit, filenameStr, boost::is_any_of(":"));

                        ASSERTL0(fSplit.size() == 1 || fSplit.size() == 2,
                                 "Incorrect filename specification in function "
                                 + functionStr + "'. "
                                 "Specify variables inside file as: "
                                 "filename:var1,var2");

                        // set the filename
                        funcDef.m_filename = fSplit[0];

                        if (fSplit.size() == 2)
                        {
                            ASSERTL0(variableList[0] != "*",
                                     "Filename variable mapping not valid "
                                     "when using * as a variable inside "
                                     "function '" + functionStr + "'.");

                            boost::split(
                                varSplit, fSplit[1], boost::is_any_of(","));
                            ASSERTL0(varSplit.size() == variableList.size(),
                                     "Filename variables should contain the "
                                     "same number of variables defined in "
                                     "VAR in function " + functionStr + "'.");
                        }
                    }

                    // Nothing else supported so throw an error
                    else
                    {
                        stringstream tagcontent;
                        tagcontent << *variable;

                        ASSERTL0(false,
                                "Identifier " + conditionType + " in function "
                                + std::string(function->Attribute("NAME"))
                                + " is not recognised in XML element: \n\t'"
                                + tagcontent.str() + "'");
                    }


                    
                    // Add variables to function
                    for (unsigned int i = 0; i < variableList.size(); ++i)
                    {
                        for(unsigned int j = 0; j < domainList.size(); ++j)
                        {
                            // Check it has not already been defined
                            pair<std::string,int> key(variableList[i],domainList[j]);
                            FunctionVariableMap::iterator fcnsIter
                                = functionVarMap.find(key);
                            ASSERTL0(fcnsIter == functionVarMap.end(),
                                     "Error setting expression '" + variableList[i]
                                     + " in domain " 
                                     + boost::lexical_cast<std::string>(domainList[j]) 
                                     + "' in function '" + functionStr + "'. "
                                     "Expression has already been defined.");

                            if (varSplit.size() > 0)
                            {
                                FunctionVariableDefinition funcDef2 = funcDef;
                                funcDef2.m_fileVariable = varSplit[i];
                                functionVarMap[key] = funcDef2;
                            }
                            else
                            {
                                functionVarMap[key] = funcDef;
                            }
                        }
                    }
                    
                    variable = variable->NextSiblingElement();
                }
                // Add function definition to map
                m_functions[functionStr] = functionVarMap;
                function = function->NextSiblingElement("FUNCTION");
            }
        }


        /**
         *
         */
        void SessionReader::ReadFilters(TiXmlElement *filters)
        {
            if (!filters)
            {
                return;
            }

            m_filters.clear();

            TiXmlElement *filter = filters->FirstChildElement("FILTER");
            while (filter)
            {
                ASSERTL0(filter->Attribute("TYPE"),
                        "Missing attribute 'TYPE' for filter.");
                std::string typeStr = filter->Attribute("TYPE");

                std::map<std::string, std::string> vParams;

                TiXmlElement *param = filter->FirstChildElement("PARAM");
                while (param)
                {
                    ASSERTL0(param->Attribute("NAME"),
                            "Missing attribute 'NAME' for parameter in filter "
                            + typeStr + "'.");
                    std::string nameStr = param->Attribute("NAME");

                    ASSERTL0(param->GetText(), "Empty value string for param.");
                    std::string valueStr = param->GetText();

                    vParams[nameStr] = valueStr;

                    param = param->NextSiblingElement("PARAM");
                }

                m_filters.push_back(
                    std::pair<std::string, FilterParams>(typeStr, vParams));

                filter = filter->NextSiblingElement("FILTER");
            }
        }
        
        void SessionReader::ParseEquals(
            const std::string &line,
                  std::string &lhs,
                  std::string &rhs)
        {
            /// Pull out lhs and rhs and eliminate any spaces.
            int beg = line.find_first_not_of(" ");
            int end = line.find_first_of("=");
            // Check for no parameter name
            if (beg == end) throw 1;
            // Check for no parameter value
            if (end != line.find_last_of("=")) throw 1;
            // Check for no equals sign
            if (end == std::string::npos) throw 1;
            
            lhs = line.substr(line.find_first_not_of(" "),
                              end-beg);
            lhs = lhs .substr(0, lhs.find_last_not_of(" ")+1);
            rhs = line.substr(line.find_last_of("=")+1);
            rhs = rhs .substr(rhs.find_first_not_of(" "));
            rhs = rhs .substr(0, rhs.find_last_not_of(" ")+1);
        }
        
        /**
         *
         */
        void SessionReader::CmdLineOverride()
        {
            // Parse solver info overrides
            if (m_cmdLineOptions.count("solverinfo"))
            {
                std::vector<std::string> solverInfoList = 
                    m_cmdLineOptions["solverinfo"].as<
                        std::vector<std::string> >();
                
                for (int i = 0; i < solverInfoList.size(); ++i)
                {
                    std::string lhs, rhs;

                    try
                    {
                        ParseEquals(solverInfoList[i], lhs, rhs);
                    } 
                    catch (...)
                    {
                        ASSERTL0(false, "Parse error with command line "
                                 "option: "+solverInfoList[i]);
                    }

                    std::string lhsUpper = boost::to_upper_copy(lhs);
                    m_solverInfo[lhsUpper] = rhs;
                }
            }
            
            if (m_cmdLineOptions.count("parameter"))
            {
                std::vector<std::string> parametersList = 
                    m_cmdLineOptions["parameter"].as<
                        std::vector<std::string> >();
                
                for (int i = 0; i < parametersList.size(); ++i)
                {
                    std::string lhs, rhs;

                    try
                    {
                        ParseEquals(parametersList[i], lhs, rhs);
                    } 
                    catch (...)
                    {
                        ASSERTL0(false, "Parse error with command line "
                                 "option: "+parametersList[i]);
                    }

                    std::string lhsUpper = boost::to_upper_copy(lhs);
                    
                    try
                    {
                        m_parameters[lhsUpper] = 
                            boost::lexical_cast<NekDouble>(rhs);
                    }
                    catch (...)
                    {
                        ASSERTL0(false, "Unable to convert string: "+rhs+
                                 "to double value.");
                    }
                }
            }
        }

        void SessionReader::SetUpXmlDoc(void)
        {
            m_xmlDoc = MergeDoc(m_filenames);
        }
    }
}
