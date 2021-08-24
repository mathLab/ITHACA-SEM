#include <cstdio>
#include <cstdlib>

#include <SolverUtils/Driver.h>
#include <LibUtilities/BasicUtils/SessionReader.h>

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>

using namespace std;
using namespace Nektar;
using namespace Nektar::SolverUtils;

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        fprintf(stderr,"Usage: ./CflStep file.xml \n");
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
        Array<OneD, NekDouble> cfl = IncNav->GetElmtCFLVals();

        // Reset Pressure field with CFL values
        Array<OneD, MultiRegions::ExpListSharedPtr> fields = IncNav->UpdateFields();
        int i,n,nquad,cnt;
        int nfields = fields.size();
        int nexp = fields[0]->GetExpSize();

        int elmtid = Vmath::Imax(nexp,cfl,1);

        cout << "Max CFL: "<< cfl[elmtid] << " In element " << elmtid << endl;


        for(n = 0; n < nfields; ++n)
        {
            if(session->GetVariable(n) == "p")
            {
                break;
            }
        }

        ASSERTL0(n != nfields, "Could not find field named p in m_fields");

        Array<OneD, NekDouble> phys = fields[n]->UpdatePhys();

        cnt = 0;
        for(i = 0; i < fields[n]->GetExpSize(); ++i)
        {
            nquad = fields[n]->GetExp(i)->GetTotPoints();
            Vmath::Fill(nquad,cfl[i],&phys[cnt],1);
            cnt += nquad;
        }

        fields[n]->FwdTrans_IterPerExp(fields[n]->GetPhys(),fields[n]->UpdateCoeffs());

        // Need to reset varibale name for output
        session->SetVariable(n,"CFL");

        // Reset session name for output file
        std::string outname = IncNav->GetSessionName();

        outname += "_CFLStep";
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

