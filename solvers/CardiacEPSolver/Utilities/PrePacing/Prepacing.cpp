#include <SpatialDomains/MeshComponents.h>
#include <MultiRegions/ExpList0D.h>
#include <CardiacEPSolver/CellModels/CellModel.h>
#include <CardiacEPSolver/Stimuli/Stimulus.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    SpatialDomains::VertexComponentSharedPtr    vPoint;
    MultiRegions::ExpListSharedPtr              vExp;
    LibUtilities::SessionReaderSharedPtr        vSession;
    std::string                                 vCellModel;
    CellModelSharedPtr                          vCell;
    std::vector<StimulusSharedPtr>              vStimulus;
    Array<OneD, Array<OneD, NekDouble> >        vWsp(1);
    Array<OneD, Array<OneD, NekDouble> >        vSol(1);
    NekDouble                                   vDeltaT;
    NekDouble                                   vTime;
    unsigned int                                nSteps;

    // Create a session reader to read pacing parameters
    vSession = LibUtilities::SessionReader::CreateInstance(argc, argv);

    try
    {
        // Construct a field consisting of a single vertex
        vPoint = MemoryManager<SpatialDomains::VertexComponent>
                        ::AllocateSharedPtr(3, 0, 0.0, 0.0, 0.0);
        vExp = MemoryManager<MultiRegions::ExpList0D>
                        ::AllocateSharedPtr(vPoint);

        // Get cell model name and create it
        vSession->LoadSolverInfo("CELLMODEL", vCellModel, "");
        ASSERTL0(vCellModel != "", "Cell Model not specified.");

        vCell = GetCellModelFactory().CreateInstance(
                                            vCellModel, vSession, vExp);
        vCell->Initialise();

        // Load the stimuli
        vStimulus = Stimulus::LoadStimuli(vSession, vExp);

        // Set up solution arrays, workspace and read in parameters
        vSol[0] = Array<OneD, NekDouble>(1, 0.0);
        vWsp[0] = Array<OneD, NekDouble>(1, 0.0);
        vDeltaT = vSession->GetParameter("TimeStep");
        vTime   = 0.0;
        nSteps  = vSession->GetParameter("NumSteps");

        vSol[0][0] = -81.0;

        // Time integrate cell model
        for (unsigned int i = 0; i < nSteps; ++i)
        {
            // Compute J_ion
            vCell->TimeIntegrate(vSol, vWsp, vTime);

            // Add stimuli J_stim
            for (unsigned int i = 0; i < vStimulus.size(); ++i)
            {
                vStimulus[i]->Update(vWsp, vTime);
            }

            // Time-step with forward Euler
            Vmath::Svtvp(1, vDeltaT, vWsp[0], 1, vSol[0], 1, vSol[0], 1);

            // Increment time
            vTime += vDeltaT;

            // Output current solution to stdout
            cout << vTime << "   " << vSol[0][0] << endl;
        }
    }
    catch (...)
    {
        cerr << "An error occured" << endl;
    }

    return 0;
}
