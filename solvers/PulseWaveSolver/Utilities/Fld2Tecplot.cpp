#include <cstdio>
#include <cstdlib>

#include <SolverUtils/Driver.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

#include <PulseWaveSolver/EquationSystems/PulseWaveSystem.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

static std::string cvar = LibUtilities::SessionReader::RegisterCmdLineFlag(
                "CharateristicVariables","c","Output characteristic variables");

static std::string SetToOneD = LibUtilities::SessionReader::RegisterCmdLineArgument("SetToOneSpaceDimension","1","Redefine mesh to be aligned to x-axis");

int main(int argc, char *argv[])
{
    if((argc < 3)||(argc > 4))
    {
        fprintf(stderr,"Usage: ./Fld2Tecplot [-c] file.xml file.fld\n");
        exit(1);
    }

    
    LibUtilities::SessionReaderSharedPtr session;
    SpatialDomains::MeshGraphSharedPtr graph;
    string vDriverModule;
    DriverSharedPtr drv;  


    try
    {
    
        // Define new input with extra argument to intialisae -OneD=false
        int newargc = argc+1;
        char **newargv = new char*[newargc];

        newargv[0] = argv[0];
        newargv[1] = new char[31];
        strcpy(newargv[1], "--SetToOneSpaceDimension=false");

        for(int i = 1; i < argc; ++i)
        {
            newargv[i+1] = argv[i];
        }
        
        // Create session reader and MeshGraph.
        session = LibUtilities::SessionReader::CreateInstance(newargc, newargv);
        graph = SpatialDomains::MeshGraph::Read(session);
        delete[] newargv;
        
        bool CalcCharacteristicVariables = false;

        if(session->DefinesCmdLineArgument(cvar))
        {
            CalcCharacteristicVariables = true;
        }

        
        // Create driver
        session->LoadSolverInfo("Driver", vDriverModule, "Standard");
        drv = GetDriverFactory().CreateInstance(vDriverModule, session, graph);
        
        EquationSystemSharedPtr EqSys = drv->GetEqu()[0];
        
        PulseWaveSystemSharedPtr PulseWave;
        if(!(PulseWave = std::dynamic_pointer_cast
             <PulseWaveSystem>(EqSys)))
        {
            ASSERTL0(false,"Failed to dynamically cast to PulseWaveSystemOutput");
        }
        
        std::string fname(argv[argc-1]);
        Array<OneD, MultiRegions::ExpListSharedPtr> Vessels;
        
        int ndomains = PulseWave->GetNdomains();

        PulseWave->ImportFldToMultiDomains(fname,Vessels = PulseWave->UpdateVessels(),
                                           ndomains);
        int fdot = fname.find_last_of('.');
        
        if (fdot != std::string::npos)
        {
            string ending = fname.substr(fdot);

            // If .chk or .fld we exchange the extension in the output file.
            // For all other files (e.g. .bse) we append the extension to avoid
            // conflicts.
            if (ending == ".chk" || ending == ".fld")
            {
                fname = fname.substr(0,fdot);
            }
        }
        
        fname = fname + ".dat";
        
        ofstream outfile(fname.c_str());
        int nvariables = session->GetVariables().size();
        std::string var = "";
        int j;
        for(j = 0; j < nvariables-1; ++j)
        {
            var += session->GetVariable(j) +  ", ";
        }
        var += session->GetVariable(j);
        
        if(CalcCharacteristicVariables)
        {
            var += ", Char1, Char2";
        }

        Vessels[0]->WriteTecplotHeader(outfile,var);
        
        for(int n = 0; n < ndomains; ++n)
        {
            Vessels[n*nvariables]->WriteTecplotZone(outfile);
            for(int j = 0; j < nvariables; ++j)
            {
                Vessels[n*nvariables+j]->WriteTecplotField(outfile);
            }

            if(CalcCharacteristicVariables)
            {
                PulseWave->CalcCharacteristicVariables(n*nvariables);
                
                for(int j = 0; j < nvariables; ++j)
                {
                    Vessels[n*nvariables+j]->WriteTecplotField(outfile);
                }
            }
            Vessels[n*nvariables]->WriteTecplotConnectivity(outfile);
        }
    }
    
    catch (const std::runtime_error&)
    {
        return 1;
    }
    catch (const std::string& eStr)
    {
        cout << "Error: " << eStr << endl;
    }

    return 0;
}        
