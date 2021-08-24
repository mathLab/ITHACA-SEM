#include <cstdio>
#include <cstdlib>

#include <SolverUtils/Driver.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <IncNavierStokesSolver/AdvectionTerms/NavierStokesAdvection.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        fprintf(stderr,"Usage: ./Aliasing file.xml \n");
        fprintf(stderr,"\t Method will read intiial conditions section of .xml file for input \n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr session;
    SpatialDomains::MeshGraphSharedPtr graph;
    string vDriverModule;
    DriverSharedPtr drv;  
    try
    {
        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);
        
        // Create MeshGraph.
        graph = SpatialDomains::MeshGraph::Read(session);

        // Create driver
        session->LoadSolverInfo("Driver", vDriverModule, "Standard");
        drv = GetDriverFactory().CreateInstance(vDriverModule, session, graph);

        EquationSystemSharedPtr EqSys = drv->GetEqu()[0];
        IncNavierStokesSharedPtr IncNav = EqSys->as<IncNavierStokes>();
        
        IncNav->SetInitialConditions(0.0,false);
        Array<OneD, MultiRegions::ExpListSharedPtr> fields = IncNav->UpdateFields();
        
        int i;
        int nConvectiveFields = IncNav->GetNConvectiveFields();
        int nphys = fields[0]->GetTotPoints();
        Array<OneD, Array<OneD, NekDouble> > VelFields(nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> > NonLinear(nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> > NonLinearDealiased(nConvectiveFields);
        
        for(i = 0; i < nConvectiveFields; ++i)
        {
            VelFields[i] = fields[i]->UpdatePhys();
            NonLinear[i] = Array<OneD, NekDouble> (nphys);
            NonLinearDealiased[i] = Array<OneD, NekDouble> (nphys);
        }

        std::shared_ptr<NavierStokesAdvection> A
            = std::dynamic_pointer_cast<NavierStokesAdvection>(IncNav->GetAdvObject());

        if (!A)
        {
            cout << "Must use non-linear Navier-Stokes advection" << endl;
            exit(-1);
        }

        // calculate non-linear terms without dealiasing
        A->SetSpecHPDealiasing(false);
        A->Advect(nConvectiveFields, fields,
                                            VelFields, VelFields, 
                                            NonLinear, 0.0);


        // calculate non-linear terms with dealiasing
        A->SetSpecHPDealiasing(true);
        A->Advect(nConvectiveFields, fields,
                                            VelFields, VelFields, 
                                            NonLinearDealiased, 0.0);

        // Evaulate Difference and put into fields;
        for(i = 0; i < nConvectiveFields; ++i)
        {
            Vmath::Vsub(nphys,NonLinearDealiased[i],1,NonLinear[i],1,NonLinear[i],1);
            fields[i]->FwdTrans_IterPerExp(NonLinear[i],fields[i]->UpdateCoeffs());
            // Need to reset varibale name for output
            string name = "NL_Aliasing_"+session->GetVariable(i);
            session->SetVariable(i,name.c_str());
        }


        // Reset session name for output file
        std::string outname = IncNav->GetSessionName();
        
        outname += "_NonLinear_Aliasing";
        IncNav->ResetSessionName(outname);
        IncNav->Output();

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

