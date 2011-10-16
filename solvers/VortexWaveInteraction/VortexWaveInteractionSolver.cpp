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

void CalcNonLinearForce(EquationSystemSharedPtr &eqn);

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
#if 1
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
#endif

        // Initialise LinNS solver 
        std::string LinNSCondFile(argv[argc-1]); 
        LinNSCondFile += "_LinNSCond.xml"; 
        std::vector<std::string> LinNSFilenames;
        LinNSFilenames.push_back(meshfile);   
        LinNSFilenames.push_back(LinNSCondFile);
        

        // Create Linearised NS stability session reader.
        LinNSsession = LibUtilities::SessionReader::CreateInstance(argc, argv, LinNSFilenames);
        // Create driver
        LinNSsession->LoadSolverInfo("Driver", vDriverModule, "ModifiedArnoldi");
        LinNSdrv = GetDriverFactory().CreateInstance(vDriverModule, LinNSsession);        
        /// Do linearised NavierStokes Session  with Modified Arnoldi
        LinNSdrv->Execute();

        // Calculate Non-linear wave terms 
        CalcNonLinearForce(LinNSdrv->GetEqu()[0]);

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
}

void CalcNonLinearForce(EquationSystemSharedPtr &eqn)
{
    Array<OneD, MultiRegions::ExpListSharedPtr> fields = eqn->UpdateFields();
    MultiRegions::ExpListSharedPtr pres = eqn->GetPressure();

    int npts    = fields[0]->GetPlane(0)->GetNpoints();
    int ncoeffs = fields[0]->GetPlane(0)->GetNcoeffs();
    Array<OneD, NekDouble> val(npts), der1(2*npts);
    Array<OneD, NekDouble> der2 = der1 + npts; 
    Array<OneD, Array<OneD, NekDouble> >  outfield(2);

    outfield[0] = Array<OneD, NekDouble> (2*ncoeffs); 
    outfield[1] = outfield[0] + ncoeffs;
    
    // determine inverse of area normalised field. 
    pres->GetPlane(0)->BwdTrans(pres->GetPlane(0)->GetCoeffs(),
                                pres->GetPlane(0)->UpdatePhys());
    pres->GetPlane(0)->BwdTrans(pres->GetPlane(1)->GetCoeffs(),
                                pres->GetPlane(1)->UpdatePhys());
    NekDouble norm = eqn->GetPressure()->L2();
    Vmath::Fill(2*npts,1.0,der1,1);
    norm = sqrt(fields[0]->PhysIntegral(der1)/(norm*norm));
    
    // Get hold of arrays. 
    fields[0]->GetPlane(0)->BwdTrans(fields[0]->GetPlane(0)->GetCoeffs(),fields[0]->GetPlane(0)->UpdatePhys());
    Array<OneD, NekDouble> u_real = fields[0]->GetPlane(0)->UpdatePhys();
    Vmath::Smul(npts,norm,u_real,1,u_real,1);
    fields[0]->GetPlane(1)->BwdTrans(fields[0]->GetPlane(1)->GetCoeffs(),fields[0]->GetPlane(1)->UpdatePhys());
    Array<OneD, NekDouble> u_imag = fields[0]->GetPlane(1)->UpdatePhys();
    Vmath::Smul(npts,norm,u_imag,1,u_imag,1);
    fields[1]->GetPlane(0)->BwdTrans(fields[1]->GetPlane(0)->GetCoeffs(),fields[1]->GetPlane(0)->UpdatePhys());
    Array<OneD, NekDouble> v_real = fields[1]->GetPlane(0)->UpdatePhys(); 
    Vmath::Smul(npts,norm,v_real,1,v_real,1);
    fields[1]->GetPlane(1)->BwdTrans(fields[1]->GetPlane(1)->GetCoeffs(),fields[1]->GetPlane(1)->UpdatePhys());
    Array<OneD, NekDouble> v_imag = fields[1]->GetPlane(1)->UpdatePhys();
    Vmath::Smul(npts,norm,v_imag,1,v_imag,1);
    
    // Calculate non-linear terms for x and y directions
    // d/dx(u u* + u* u)
    Vmath::Vmul (npts,u_real,1,u_real,1,val,1);
    Vmath::Vvtvp(npts,u_imag,1,u_imag,1,val,1,val,1);
    Vmath::Smul (npts,2.0,val,1,val,1);
    fields[0]->GetPlane(0)->PhysDeriv(0,val,der1);
    
    // d/dy(v u* + v* u)
    Vmath::Vmul (npts,u_real,1,v_real,1,val,1);
    Vmath::Vvtvp(npts,u_imag,1,v_imag,1,val,1,val,1);
    Vmath::Smul (npts,2.0,val,1,val,1);
    fields[0]->GetPlane(0)->PhysDeriv(1,val,der2);
    
    Vmath::Vadd(npts,der1,1,der2,1,der1,1);
    
    fields[0]->GetPlane(0)->FwdTrans_IterPerExp(der1,outfield[0]);
    Vmath::Neg(ncoeffs,outfield[0],1);
    
    // d/dx(u v* + u* v)
    fields[0]->GetPlane(0)->PhysDeriv(0,val,der1);
    
    // d/dy(v v* + v* v)
    Vmath::Vmul(npts,v_real,1,v_real,1,val,1);
    Vmath::Vvtvp(npts,v_imag,1,v_imag,1,val,1,val,1);
    Vmath::Smul (npts,2.0,val,1,val,1);
    fields[0]->GetPlane(0)->PhysDeriv(1,val,der2);
    
    fields[0]->GetPlane(0)->FwdTrans_IterPerExp(der1,outfield[1]);
    Vmath::Neg(ncoeffs,outfield[1],1);

    // Symmetrise forcing
    //-> Get coordinates 
    Array<OneD, NekDouble> coord(2);
    Array<OneD, NekDouble> coord_x(npts);
    Array<OneD, NekDouble> coord_y(npts);
    
    //-> Impose symmetry (x -> -x + Lx/2, y-> -y) on coordinates
    fields[0]->GetPlane(0)->GetCoords(coord_x,coord_y);
    NekDouble xmax = Vmath::Vmax(npts,coord_x,1);
    Vmath::Neg(npts,coord_x,1);
    Vmath::Sadd(npts,xmax,coord_x,1,coord_x,1);
    Vmath::Neg(npts,coord_y,1);

    int i, physoffset;

    //-> Obtain list of expansion element ids for each point. 
    Array<OneD, int> Eid(npts);
    // This search may not be necessary every iteration
    for(i = 0; i < npts; ++i)
    {
        coord[0] = coord_x[i];
        coord[1] = coord_y[i];

        // Note this will not quite be symmetric. 
        Eid[i] = fields[0]->GetPlane(0)->GetExpIndex(coord,1e-6);
    }

    // Interpolate field 0 
    fields[0]->GetPlane(0)->BwdTrans_IterPerExp(outfield[0],der1);
    for(i = 0; i < npts; ++i)
    {
        physoffset = fields[0]->GetPlane(0)->GetPhys_Offset(Eid[i]);
        coord[0] = coord_x[i];
        coord[1] = coord_y[i];
        der2[i] = fields[0]->GetPlane(0)->GetExp(Eid[i])->PhysEvaluate(coord,
                                                          der1 + physoffset);
    }
    //-> Average field 0 
    Vmath::Vsub(npts,der1,1,der2,1,der2,1);
    Vmath::Smul(npts,0.5,der2,1,der2,1);
    fields[0]->GetPlane(0)->FwdTrans_IterPerExp(der2, outfield[0]);
    
    //-> Interpoloate field 1
    fields[0]->GetPlane(0)->BwdTrans_IterPerExp(outfield[1],der1);
    for(i = 0; i < npts; ++i)
    {
        physoffset = fields[0]->GetPlane(0)->GetPhys_Offset(Eid[i]);
        coord[0] = coord_x[i];
        coord[1] = coord_y[i];
        der2[i] = fields[0]->GetPlane(0)->GetExp(Eid[i])->PhysEvaluate(coord,
                                                         der1 + physoffset);
    }

    //-> Average field 1
    Vmath::Vsub(npts,der1,1,der2,1,der2,1);
    Vmath::Smul(npts,0.5,der2,1,der2,1);
    fields[0]->GetPlane(0)->FwdTrans_IterPerExp(der2, outfield[1]);
    
    // dump output
    Array<OneD, std::string> variables(2);
    variables[0] = "u";   variables[1] = "v";

    std::string outname = eqn->GetSessionName().substr(0,eqn->GetSessionName().find_last_of('.')) + ".vwi";

    eqn->WriteFld(outname, fields[0]->GetPlane(0), outfield, variables);
}
