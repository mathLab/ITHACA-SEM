///////////////////////////////////////////////////////////////////////////////
//
// File VortexWaveInteractionSolver.cpp
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
// Description: Vortex Wave Interaction solver
//
///////////////////////////////////////////////////////////////////////////////

#include <Auxiliary/Driver.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    
    if(argc != 2)
    {
        cerr << "\n \t Usage: VortexWaveInteractionSolver  input \n" << endl;
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr IncNSsession, AdvDiffsession, LinNSsession;
    string vDriverModule;
    DriverSharedPtr LinNSdrv;
    int pargc = 2;
    std::string meshfile(argv[argc-1]);
    meshfile += ".xml";

    // system call strings; 
    string FldToRst("cp -f ");
    FldToRst += argv[argc-1]; 
    FldToRst += ".fld ";
    FldToRst += argv[argc-1]; 
    FldToRst += ".rst";

    string FldToBase("cp -f ");
    FldToBase += argv[argc-1]; 
    FldToBase += ".fld ";
    FldToBase += argv[argc-1]; 
    FldToBase += "-Base.fld";

    string FldToStreak("cp -f ");
    FldToStreak += argv[argc-1]; 
    FldToStreak += ".fld ";
    FldToStreak += argv[argc-1]; 
    FldToStreak += "_streak.fld";

    try
    {
        // Initialise NS solver 
        std::string IncCondFile(argv[argc-1]); 
        IncCondFile += "_IncNSCond.xml"; 
        std::vector<std::string> IncNSFilenames;
        IncNSFilenames.push_back(meshfile);   
        IncNSFilenames.push_back(IncCondFile);
        
        // Initialise AdvDiff solver 
        std::string AdvDiffCondFile(argv[argc-1]); 
        AdvDiffCondFile += "_AdvDiffCond.xml"; 
        std::vector<std::string> AdvDiffFilenames;
        AdvDiffFilenames.push_back(meshfile);   
        AdvDiffFilenames.push_back(AdvDiffCondFile);


        // Initialise LinNS solver 
        std::string LinNSCondFile(argv[argc-1]); 
        LinNSCondFile += "_LinNSCond.xml"; 
        std::vector<std::string> LinNSFilenames;
        LinNSFilenames.push_back(meshfile);   
        LinNSFilenames.push_back(LinNSCondFile);
        
        // Create IncompressibleNavierStokesSolver session reader.
        IncNSsession = LibUtilities::SessionReader::CreateInstance(argc, argv, IncNSFilenames);
        // Setup Incompresible solver and execute 
        std::string vEquation = IncNSsession->GetSolverInfo("SolverType");
        EquationSystemSharedPtr NSeqn = GetEquationSystemFactory().CreateInstance(vEquation,IncNSsession);

        NSeqn->PrintSummary(cout);
        NSeqn->DoInitialise();
        NSeqn->DoSolve();
        NSeqn->Output();

        // Copy .fld file to .rst and base.fld
        cout << "Executing " <<  FldToRst << endl;
        system(FldToRst.c_str());
        cout << "Executing " <<  FldToBase << endl;
        system(FldToBase.c_str());


        // Create AdvDiffusion session reader.
        AdvDiffsession = LibUtilities::SessionReader::CreateInstance(argc, argv, AdvDiffFilenames);
        // Setup and execute Advection Diffusion solver 
        vEquation = AdvDiffsession->GetSolverInfo("EqType");
        EquationSystemSharedPtr AdvDiffeqn = GetEquationSystemFactory().CreateInstance(vEquation,AdvDiffsession);

        AdvDiffeqn->PrintSummary(cout);
        AdvDiffeqn->DoInitialise();
        AdvDiffeqn->DoSolve();
        AdvDiffeqn->Output();
        
        cout << "Executing " <<  FldToStreak << endl;
        system(FldToStreak.c_str());

        // Create Linearised NS stability session reader.
        LinNSsession = LibUtilities::SessionReader::CreateInstance(argc, argv, LinNSFilenames);
        // Create driver
        LinNSsession->LoadSolverInfo("Driver", vDriverModule, "ModifiedArnoldi");
        LinNSdrv = GetDriverFactory().CreateInstance(vDriverModule, LinNSsession);        
        /// Do NavierStokes Session 
        LinNSdrv->Execute();

        // Finalise communications
        IncNSsession->Finalise();
        AdvDiffsession->Finalise();
        LinNSsession->Finalise();
        
    }
    catch (const std::runtime_error& e)
    {
        return 1;
    }
    catch (const std::string& eStr)
    {
        cout << "Error: " << eStr << endl;
    }
    
    return 0;
}
