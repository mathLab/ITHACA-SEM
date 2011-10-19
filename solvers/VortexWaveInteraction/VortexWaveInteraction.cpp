///////////////////////////////////////////////////////////////////////////////
//
// File VortexWaveInteraction.cpp
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
// Description: Vortex Wave Interaction class
//
///////////////////////////////////////////////////////////////////////////////
#include <VortexWaveInteraction/VortexWaveInteraction.h>

namespace Nektar
{
    
    VortexWaveInteraction::VortexWaveInteraction(int argc, char * argv[])
    {
        m_sessionName = argv[argc-1];
        string meshfile = m_sessionName + ".xml";

        LibUtilities::SessionReaderSharedPtr session; 

        std::string VWICondFile(argv[argc-1]); 
        VWICondFile += "_VWI.xml"; 
        std::vector<std::string> VWIFilenames;
        VWIFilenames.push_back(meshfile);   
        VWIFilenames.push_back(VWICondFile);

        // Create Incompressible NavierStokesSolver session reader.
        m_sessionVWI = LibUtilities::SessionReader::CreateInstance(argc, argv, VWIFilenames);
        
        
        if(m_sessionVWI->DefinesParameter("StartIteration"))
        {
            m_iterStart = m_sessionVWI->GetParameter("StartIteration");
        }
        else
        {
            m_iterStart = 0; 
        }

        if(m_sessionVWI->DefinesParameter("EndIteration"))
        {
            m_iterEnd = m_sessionVWI->GetParameter("EndIteration");
        }
        else
        {
            m_iterEnd = 1; 
        }

        m_rho   = m_sessionVWI->GetParameter("Rho");
        m_alpha = m_sessionVWI->GetParameter("Alpha");
        
        // Initialise NS Roll solver 
        std::string IncCondFile(argv[argc-1]); 
        IncCondFile += "_IncNSCond.xml"; 
        std::vector<std::string> IncNSFilenames;
        IncNSFilenames.push_back(meshfile);   
        IncNSFilenames.push_back(IncCondFile);

        // Create Incompressible NavierStokesSolver session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv, IncNSFilenames);
        std::string vEquation = session->GetSolverInfo("SolverType");
        m_solverRoll = GetEquationSystemFactory().CreateInstance(vEquation, session);

        // Create AdvDiff Streak solver 
        std::string AdvDiffCondFile(argv[argc-1]); 
        AdvDiffCondFile += "_AdvDiffCond.xml"; 
        std::vector<std::string> AdvDiffFilenames;
        AdvDiffFilenames.push_back(meshfile);   
        AdvDiffFilenames.push_back(AdvDiffCondFile);

        // Create AdvDiffusion session reader.
        m_sessionStreak = LibUtilities::SessionReader::CreateInstance(argc, argv, AdvDiffFilenames);

        // Initialise LinNS solver 
        std::string LinNSCondFile(argv[argc-1]); 
        LinNSCondFile += "_LinNSCond.xml"; 
        std::vector<std::string> LinNSFilenames;
        LinNSFilenames.push_back(meshfile);   
        LinNSFilenames.push_back(LinNSCondFile);
        
        // Create Linearised NS stability session reader.
        m_sessionWave = LibUtilities::SessionReader::CreateInstance(argc, argv, LinNSFilenames);

        // Set the initial beta value in stability to be equal to VWI file
        std::string LZstr("LZ");
        NekDouble LZ = 2*M_PI/m_alpha;
        cout << "Setting LZ in Linearised solver to " << LZ << endl;
        m_sessionWave->LoadParameter(LZstr,LZ);


        // Check for restart
        bool restart;
        m_sessionVWI->MatchSolverInfo("RestartIteration","True",restart,false);
        if(restart)
        {
            string nstr =  boost::lexical_cast<std::string>(m_iterStart);
            cout << "Restarting from iteration " << m_iterStart << endl;
            std::string rstfile = "cp -f Save/" + m_sessionName + ".rst." + nstr + " " + m_sessionName + ".rst"; 
            system(rstfile.c_str());
            std::string vwifile = "cp -f Save/" + m_sessionName + ".vwi." + nstr + " " + m_sessionName + ".vwi"; 
            system(vwifile.c_str());
        }

        // system call strings; 
        m_fldToRst    = "cp -f " + m_sessionName + ".fld " + m_sessionName + ".rst";
        m_fldToBase   = "cp -f " + m_sessionName + ".fld " + m_sessionName + "-Base.fld";
        m_fldToStreak = "cp -f " + m_sessionName + ".fld " + m_sessionName + "_streak.fld";
    }


    VortexWaveInteraction::~VortexWaveInteraction()
    {
        m_solverRoll->GetSession()->Finalise();
        m_sessionStreak->Finalise();
        m_sessionWave->Finalise();
    }
    
    void VortexWaveInteraction::ExecuteRoll(void)
    {
        // Execute Roll 
        cout << "Executing Roll solver" << endl;
        m_solverRoll->DoInitialise();
        m_solverRoll->DoSolve();
        m_solverRoll->Output();
        
        // Copy .fld file to .rst and base.fld
        cout << "Executing " <<  m_fldToRst << endl;
        system(m_fldToRst.c_str());
        cout << "Executing " <<  m_fldToBase << endl;
        system(m_fldToBase.c_str());
    }


    void VortexWaveInteraction::ExecuteStreak(void)
    {
        // Setup and execute Advection Diffusion solver 
        string vEquation = m_sessionStreak->GetSolverInfo("EqType");
        EquationSystemSharedPtr solverStreak = GetEquationSystemFactory().CreateInstance(vEquation,m_sessionStreak);

        cout << "Executing Streak Solver" << endl;
        solverStreak->DoInitialise();
        solverStreak->DoSolve();
        solverStreak->Output();
        
        cout << "Executing " <<  m_fldToStreak << endl;
        system(m_fldToStreak.c_str());
    }

    void VortexWaveInteraction::ExecuteWaveAndForce(void)
    {
        // Create driver
        std::string vDriverModule;
        m_sessionWave->LoadSolverInfo("Driver", vDriverModule, "ModifiedArnoldi");
        cout << "Setting up linearised NS sovler" << endl;
        DriverSharedPtr solverWave = GetDriverFactory().CreateInstance(vDriverModule, m_sessionWave);  

        /// Do linearised NavierStokes Session  with Modified Arnoldi
        cout << "Executing wave solution " << endl;
        solverWave->Execute();

        m_leading_real_evl = solverWave->GetRealEvl()[0];
        m_leading_imag_evl = solverWave->GetImagEvl()[0];
        
        cout << "Calculating Nonlinear Wave Forcing" << endl;
        CalcNonLinearWaveForce(solverWave->GetEqu()[0]);
    }

    void VortexWaveInteraction::CalcNonLinearWaveForce(EquationSystemSharedPtr &eqn)
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
    
    void VortexWaveInteraction::SaveFile(string fileend, string dir, int n)
    {
        static bool init = true; 

        if(init == true)
        {
            // make directory and presume will fail if it already exists
            string mkdir = "mkdir " + dir;
            system(mkdir.c_str());
            init = false;
        }
        
        string cpfile   = m_sessionName + fileend;
        string savefile = dir + "/" + cpfile + "." + boost::lexical_cast<std::string>(n);
        string syscall  = "cp -f "  + cpfile + " " + savefile; 

        system(syscall.c_str());
    }

    void VortexWaveInteraction::AppendEvlToFile(string file, int n)
    {
        FILE *fp;
        fp = fopen(file.c_str(),"a");

        NekDouble invmag = 1.0/(m_leading_real_evl*m_leading_real_evl + 
                                m_leading_imag_evl*m_leading_imag_evl);
        fprintf(fp, "%d: %lf %16.12le  %16.12le\n",n, m_alpha, -m_leading_real_evl*invmag,m_leading_imag_evl*invmag);
        fclose(fp);
    }

}
