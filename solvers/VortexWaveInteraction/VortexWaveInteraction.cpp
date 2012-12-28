//////////////////////////////////////////////////////////////////////////////
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
#include <MultiRegions/GlobalLinSysKey.h>
#include <MultiRegions/ExpList1D.h>
#include <SolverUtils/Driver.h>

namespace Nektar
{

    VortexWaveInteraction::VortexWaveInteraction(int argc, char * argv[]):
        m_nOuterIterations(0)
    {

        int storesize;// number of previous iterations to store

        m_sessionName = argv[argc-1];
        string meshfile = m_sessionName + ".xml";

        LibUtilities::SessionReaderSharedPtr session; 

        std::string VWICondFile(argv[argc-1]); 
        VWICondFile += "_VWI.xml"; 
        std::vector<std::string> VWIFilenames;
        VWIFilenames.push_back(meshfile);   
        VWIFilenames.push_back(VWICondFile);

        // Create Incompressible NavierStokesSolver session reader.
        m_sessionVWI = LibUtilities::SessionReader::CreateInstance(argc, argv, 
                                                                   VWIFilenames);
        
        m_sessionVWI->LoadParameter("AlphaStep",              m_alphaStep,0.05);
        m_sessionVWI->LoadParameter("OuterIterationStoreSize",storesize,10);

        //set a low tol for interfaceVWI
	m_sessionVWI->LoadParameter("EigenvalueRelativeTol",  m_eigRelTol,1e-3);
       
        
        m_sessionVWI->LoadParameter("NeutralPointTolerance",  m_neutralPointTol,1e-4);
        m_sessionVWI->LoadParameter("MaxOuterIterations",     m_maxOuterIterations,100);
        
        m_sessionVWI->LoadParameter("StartIteration",m_iterStart, 0);
        m_sessionVWI->LoadParameter("EndIteration",  m_iterEnd, 0);

        m_sessionVWI->LoadParameter("WaveForceMagStep",   m_waveForceMagStep,0.01);
        m_sessionVWI->LoadParameter("DAlphaDWaveForceMag", m_dAlphaDWaveForceMag,0.0);
        m_sessionVWI->LoadParameter("MaxWaveForceMagIter",m_maxWaveForceMagIter,1);
        m_sessionVWI->LoadParameter("RollForceScale",     m_rollForceScale,1.0);
        
        if(m_sessionVWI->DefinesSolverInfo("DeltaFcnApprox"))
        {
            m_deltaFcnApprox = true;
            m_sessionVWI->LoadParameter("DeltaFcnDecay", m_deltaFcnDecay,1.0/500);
        }
        else
        {
            m_deltaFcnApprox = false;
            m_deltaFcnDecay = 0.0;
        }

        if(m_sessionVWI->DefinesSolverInfo("LinfPressureNorm"))
        {
            m_useLinfPressureNorm  = true;
        }
        else
        {
            m_useLinfPressureNorm  = false;
        }

        if( m_sessionVWI->DefinesSolverInfo("INTERFACE")  )
        {
            m_iterinterface = true;
        }
        else
        {
            m_iterinterface = false;
        }

        if(m_sessionVWI->DefinesSolverInfo("MoveMeshToCriticalLayer"))
        {
            m_moveMeshToCriticalLayer = true;
        }
        else
        {
            m_moveMeshToCriticalLayer = false;
        }

        m_alpha = Array<OneD, NekDouble> (storesize);
        m_alpha[0]         = m_sessionVWI->GetParameter("Alpha");
        m_waveForceMag     = Array<OneD, NekDouble> (storesize);
        m_waveForceMag[0]  = m_sessionVWI->GetParameter("WaveForceMag");
        
        m_leading_real_evl = Array<OneD, NekDouble> (storesize);
        m_leading_imag_evl = Array<OneD, NekDouble> (storesize);
        
        if(m_sessionVWI->DefinesParameter("Relaxation"))
        {
            m_vwiRelaxation = m_sessionVWI->GetParameter("Relaxation");
            // fix minimum number of iterations to be number of
            // iterations required to make contribution of innitial
            // forcing to 0.1
            m_minInnerIterations = (int) (log(0.1)/log(m_vwiRelaxation)); 
        }
        else
        {
            m_vwiRelaxation = 0.0;
            m_minInnerIterations = 1;
        }
        
        // Initialise NS Roll solver 
        std::string IncCondFile(argv[argc-1]); 
        IncCondFile += "_IncNSCond.xml"; 
        std::vector<std::string> IncNSFilenames;
        IncNSFilenames.push_back(meshfile);   
        IncNSFilenames.push_back(IncCondFile);

        // Create Incompressible NavierStokesSolver session reader.
        m_sessionRoll = LibUtilities::SessionReader::CreateInstance(argc, argv, IncNSFilenames, 
                                                                    m_sessionVWI->GetComm());
        std::string vEquation = m_sessionRoll->GetSolverInfo("SolverType");
        m_solverRoll = GetEquationSystemFactory().CreateInstance(vEquation, m_sessionRoll);
        m_solverRoll->PrintSummary(cout);


        if(m_sessionRoll->DefinesSolverInfo("INTERFACE"))
	{
            m_sessionVWI->LoadParameter("EigenvalueRelativeTol",  m_eigRelTol,1e-2);
	}

        int ncoeffs = m_solverRoll->UpdateFields()[0]->GetNcoeffs();
        m_vwiForcing = Array<OneD, Array<OneD, NekDouble> > (4);
        m_vwiForcing[0] = Array<OneD, NekDouble> (4*ncoeffs);
        for(int i = 1; i < 4; ++i)
        {
            m_vwiForcing[i] = m_vwiForcing[i-1] + ncoeffs;
        }

        // Fill forcing into m_vwiForcing
        Vmath::Vcopy(ncoeffs,m_solverRoll->UpdateForces()[0]->GetCoeffs(),1,m_vwiForcing[0],1);
        Vmath::Vcopy(ncoeffs,m_solverRoll->UpdateForces()[1]->GetCoeffs(),1,m_vwiForcing[1],1);
        

        // Create AdvDiff Streak solver 
        std::string AdvDiffCondFile(argv[argc-1]); 
        AdvDiffCondFile += "_AdvDiffCond.xml"; 
        std::vector<std::string> AdvDiffFilenames;
        AdvDiffFilenames.push_back(meshfile);   
        AdvDiffFilenames.push_back(AdvDiffCondFile);
        
        // Create AdvDiffusion session reader.
        m_sessionStreak = LibUtilities::SessionReader::CreateInstance(argc, argv, AdvDiffFilenames, m_sessionVWI->GetComm());

        // Initialise LinNS solver 
        std::string LinNSCondFile(argv[argc-1]); 
        LinNSCondFile += "_LinNSCond.xml"; 
        std::vector<std::string> LinNSFilenames;
        LinNSFilenames.push_back(meshfile);   
        LinNSFilenames.push_back(LinNSCondFile);
        
        // Create Linearised NS stability session reader.
        m_sessionWave = LibUtilities::SessionReader::CreateInstance(argc, argv, LinNSFilenames, m_sessionVWI->GetComm());

        // Set the initial beta value in stability to be equal to VWI file
        std::string LZstr("LZ");
        NekDouble LZ = 2*M_PI/m_alpha[0];
        cout << "Setting LZ in Linearised solver to " << LZ << endl;
        m_sessionWave->SetParameter(LZstr,LZ);

        // Check for iteration type 
        if(m_sessionVWI->DefinesSolverInfo("VWIIterationType"))
        {
            std::string IterationTypeStr = m_sessionVWI->GetSolverInfo("VWIIterationType");
            for(int i = 0; i < (int) eVWIIterationTypeSize; ++i)
            {
                if(m_solverRoll->NoCaseStringCompare(VWIIterationTypeMap[i],IterationTypeStr) == 0 )
                {
                    m_VWIIterationType = (VWIIterationType)i; 
                    break;
                }
            }
        }
        else
        {
            m_VWIIterationType = eFixedAlphaWaveForcing;
        }



        // Check for restart
        bool restart;
        m_sessionVWI->MatchSolverInfo("RestartIteration","True",restart,false);
        if(restart)
        {
            switch(m_VWIIterationType)
            {
            case  eFixedAlpha:
                break;
            case  eFixedWaveForcing:
                {
                    FILE *fp;
                    // Check for OuterIter.his file to read
                    if((fp = fopen("OuterIter.his","r")))
                    {
                        char buf[BUFSIZ];
                        std::vector<NekDouble> Alpha, Growth, Phase;
                        NekDouble alpha,growth,phase;
                        while(fgets(buf,BUFSIZ,fp))
                        {
                            sscanf(buf,"%*d:%lf%lf%lf",&alpha,&growth,&phase);
                            Alpha.push_back(alpha);
                            Growth.push_back(growth);
                            Phase.push_back(phase);
                        }

                        m_nOuterIterations = Alpha.size();
                        
                        int nvals = std::min(m_nOuterIterations,(int)m_alpha.num_elements());
                        
                        for(int i = 0; i < nvals; ++i)
                        {
                            m_alpha[nvals-1-i]            = Alpha[m_nOuterIterations-nvals+i];
                            m_leading_real_evl[nvals-1-i] = Growth[m_nOuterIterations-nvals+i];
                            m_leading_imag_evl[nvals-1-i] = Phase [m_nOuterIterations-nvals+i];
                        }

                        UpdateAlpha(m_nOuterIterations++);
                    }
                    else
                    {
                        cout << " No File OuterIter.his to restart from" << endl;
                    }
                }
                break;
            case eFixedAlphaWaveForcing:
                {
                    string nstr =  boost::lexical_cast<std::string>(m_iterStart);
                    cout << "Restarting from iteration " << m_iterStart << endl;
                    std::string rstfile = "cp -f Save/" + m_sessionName + ".rst." + nstr + " " + m_sessionName + ".rst"; 
                    cout << "      " << rstfile << endl;
                    if(system(rstfile.c_str()))
                    {
                        ASSERTL0(false,rstfile.c_str());
                    }
                    std::string vwifile = "cp -f Save/" + m_sessionName + ".vwi." + nstr + " " + m_sessionName + ".vwi"; 
                    cout << "      " << vwifile << endl;
                    if(system(vwifile.c_str()))
                    {
                        ASSERTL0(false,vwifile.c_str());
                    }
                }
                break;
            default:
                ASSERTL0(false,"Unknown VWIITerationType in restart");
            }
        }

        // Check for ConveredSoln to update DAlphaDWaveForce
        {
            FILE *fp;
            if((fp = fopen("ConvergedSolns","r")))
            {
                char buf[BUFSIZ];
                std::vector<NekDouble> WaveForce, Alpha;
                NekDouble waveforce,alpha;
                while(fgets(buf,BUFSIZ,fp))
                {
                    sscanf(buf,"%*d:%lf%lf",&waveforce,&alpha);
                    WaveForce.push_back(waveforce);
                    Alpha.push_back(alpha);
                }

                if(Alpha.size() > 1)
                {
                    int min_i = 0;
                    NekDouble min_alph = fabs(m_alpha[0]-Alpha[min_i]);
                    // find nearest point 
                    for(int i = 1; i < Alpha.size(); ++i)
                    {
                        if(fabs(m_alpha[0]-Alpha[min_i]) < min_alph)
                        {
                            min_i = i;
                            min_alph = fabs(m_alpha[0]-Alpha[min_i]);
                        }
                    }

                    // find next nearest point 
                    int min_j = (min_i == 0)? 1:0;
                    min_alph = fabs(m_alpha[0]-Alpha[min_j]);
                    for(int i = 0; i < Alpha.size(); ++i)
                    {
                        if(i != min_i)
                        {
                            if(fabs(m_alpha[0]-Alpha[min_j]) < min_alph)
                            {
                                min_j = i;
                                min_alph = fabs(m_alpha[0]-Alpha[min_j]);
                            }
                        }
                    }
                    
                    if(fabs(Alpha[min_i] - Alpha[min_j]) > 1e-4)
                    {
                        m_dAlphaDWaveForceMag = (Alpha[min_i]-Alpha[min_j])/(WaveForce[min_i]-WaveForce[min_j]);
                    }
                }
            }
        }
    }


    VortexWaveInteraction::~VortexWaveInteraction()
    {
        m_sessionVWI->Finalise();
    }
    
    void VortexWaveInteraction::ExecuteRoll(void)
    {
        //set up the equation system to update the mesh
        if(m_sessionRoll->DefinesSolverInfo("INTERFACE"))
        {
            string vEquation = m_sessionRoll->GetSolverInfo("solvertype");
            EquationSystemSharedPtr solverRoll = GetEquationSystemFactory().CreateInstance(vEquation,m_sessionRoll);
            // The forcing terms are inserted as N bcs
            // Execute Roll 
            cout << "Executing Roll solver" << endl;
            solverRoll->DoInitialise();
            solverRoll->DoSolve();
            solverRoll->Output();
            m_rollField = solverRoll->UpdateFields();
            for(int g=0; g< solverRoll->GetNvariables(); ++g)
            {
                NekDouble vL2Error = solverRoll->L2Error(g,false);
                NekDouble vLinfError = solverRoll->LinfError(g);
                cout << "L 2 error (variable " << solverRoll->GetVariable(g) << ") : " << vL2Error << endl;
                cout << "L inf error (variable " << solverRoll->GetVariable(g) << ") : " << vLinfError << endl;
            }  
        }
        else
        {
            if(m_moveMeshToCriticalLayer)
            {
                string vEquation = m_sessionRoll->GetSolverInfo("solvertype");
                EquationSystemSharedPtr solverRoll = GetEquationSystemFactory().CreateInstance(vEquation,m_sessionRoll);
            }
            else
            {
                static int init = 1;
                if(init)
                {
                    // Read vwi file
                    std::string forcefile
                        = m_sessionRoll->GetFunctionFilename("BodyForce", 0);
                    
                    m_solverRoll->DoInitialise();
                    
                    if(forcefile != "")
                    {
                        m_solverRoll->ImportFld(forcefile,m_solverRoll->UpdateForces());
                        int ncoeffs = m_solverRoll->UpdateForces()[0]->GetNcoeffs();
                        
                        // Scale forcing
                        int npoints = m_solverRoll->UpdateForces()[0]->GetNpoints();
                        for(int i = 0; i < m_solverRoll->UpdateForces().num_elements(); ++i)
                        {
                            Vmath::Smul(npoints,m_rollForceScale,m_solverRoll->UpdateForces()[i]->UpdatePhys(),1,m_solverRoll->UpdateForces()[i]->UpdatePhys(),1);
                            
                            Vmath::Vcopy(ncoeffs,m_solverRoll->UpdateForces()[i]->GetCoeffs(),1,m_vwiForcing[2+i],1);

                        }
                    }
                    
                    init = 0;
                }
                else // use internal definition of forcing in m_vwiForcing
                {
                    // Scale forcing
                    int npoints = m_solverRoll->UpdateForces()[0]->GetNpoints();
                    Array<OneD, NekDouble> physForce(npoints);
                    for(int i = 0; i < m_solverRoll->UpdateForces().num_elements(); ++i)
                    {
                        m_solverRoll->UpdateForces()[i]->BwdTrans(m_vwiForcing[i],physForce);
                        Vmath::Smul(npoints,m_rollForceScale,physForce,1,m_solverRoll->UpdateForces()[i]->UpdatePhys(),1);
                    }
                    
                    // Shift m_vwiForcing for new restart in case of relaxation 
                    int ncoeffs = m_solverRoll->UpdateForces()[0]->GetNcoeffs();
                    Vmath::Vcopy(ncoeffs,m_vwiForcing[0],1,m_vwiForcing[2],1);
                    Vmath::Vcopy(ncoeffs,m_vwiForcing[1],1,m_vwiForcing[3],1);
                }
            }

            // Execute Roll 
            cout << "Executing Roll solver" << endl;
            m_solverRoll->DoSolve();
            m_solverRoll->Output();
            m_rollField = m_solverRoll->UpdateFields();
            for(int g=0; g< m_solverRoll->GetNvariables(); ++g)
            {
                NekDouble vL2Error = m_solverRoll->L2Error(g,false);
                NekDouble vLinfError = m_solverRoll->LinfError(g);
                cout << "L 2 error (variable " << m_solverRoll->GetVariable(g) << ") : " << vL2Error << endl;
                cout << "L inf error (variable " << m_solverRoll->GetVariable(g) << ") : " << vLinfError << endl;
            }  


        }
        
        // Copy .fld file to .rst and base.fld
        cout << "Executing cp -f session.fld session.rst" << endl;
        CopyFile(".fld",".rst");
        
        // Write out data into base flow with variable Vx,Vy
        cout << "Writing data to session-Base.fld" << endl;
        
        Array<OneD, std::string> variables(2);
        variables[0] = "Vx";   variables[1] = "Vy";
        Array<OneD, Array<OneD, NekDouble> > outfield(2);
        outfield[0]  = m_solverRoll->UpdateFields()[0]->UpdateCoeffs(); 
        outfield[1]  = m_solverRoll->UpdateFields()[1]->UpdateCoeffs(); 
        std::string outname = m_sessionName  + "-Base.fld";
        m_solverRoll->WriteFld(outname, m_solverRoll->UpdateFields()[0], 
                               outfield, variables);
    }


    void VortexWaveInteraction::ExecuteStreak(void)
    {
        // Create driver
#if 1
        std::string vDriverModule;
        m_sessionStreak->LoadSolverInfo("Driver", vDriverModule, "Standard");
        
        DriverSharedPtr solverStreak = GetDriverFactory().CreateInstance(vDriverModule, m_sessionStreak); 
        solverStreak->Execute();

        m_streakField = solverStreak->GetEqu()[0]->UpdateFields();
#else        
        // Setup and execute Advection Diffusion solver 
        string vEquation = m_sessionStreak->GetSolverInfo("EqType");
        EquationSystemSharedPtr solverStreak = GetEquationSystemFactory().CreateInstance(vEquation,m_sessionStreak);

        cout << "Executing Streak Solver" << endl;
        solverStreak->DoInitialise();
        solverStreak->DoSolve();
        solverStreak->Output();

        m_streakField = solverStreak->UpdateFields();

        if(m_sessionVWI->DefinesSolverInfo("INTERFACE"))
        {
            for(int g=0; g< solverStreak->GetNvariables(); ++g)
            {
                NekDouble vL2Error = solverStreak->L2Error(g,false);
                NekDouble vLinfError = solverStreak->LinfError(g);
                cout << "L 2 error (variable " << solverStreak->GetVariable(g) << ") : " << vL2Error << endl;
                cout << "L inf error (variable " << solverStreak->GetVariable(g) << ") : " << vLinfError << endl;
            }
        }
#endif

        cout << "Executing cp -f session.fld session_streak.fld" << endl;
        CopyFile(".fld","_streak.fld");
    }

    void VortexWaveInteraction::ExecuteWave(void)
    {

        // Set the initial beta value in stability to be equal to VWI file
        std::string LZstr("LZ");
        NekDouble LZ = 2*M_PI/m_alpha[0];
        cout << "Setting LZ in Linearised solver to " << LZ << endl;
        m_sessionWave->SetParameter(LZstr,LZ);

        // Create driver
        std::string vDriverModule;
        m_sessionWave->LoadSolverInfo("Driver", vDriverModule, "ModifiedArnoldi");
        cout << "Setting up linearised NS sovler" << endl;
        DriverSharedPtr solverWave = GetDriverFactory().CreateInstance(vDriverModule, m_sessionWave);  

        /// Do linearised NavierStokes Session  with Modified Arnoldi
        cout << "Executing wave solution " << endl;
        solverWave->Execute();

        // Copy file to a rst location for next restart
        cout << "Executing cp -f session_eig_0 session_eig_0.rst" << endl;
        CopyFile("_eig_0","_eig_0.rst");

        // Store data relevant to other operations 
        m_leading_real_evl[0] = solverWave->GetRealEvl()[0];
        m_leading_imag_evl[0] = solverWave->GetImagEvl()[0];

        // note this will only be true for modified Arnoldi
        NekDouble realShift = 0.0;
        if(m_sessionWave->DefinesParameter("RealShift"))
        {
            bool defineshift;
            // only use shift in modifiedArnoldi solver since
            // implicitly handled in Arpack.
            m_sessionWave->MatchSolverInfo("Driver","ModifiedArnoldi",defineshift,true);
            if(defineshift)
            {
                realShift = m_sessionWave->GetParameter("RealShift");
            }

        }

        // Set up leading eigenvalue from inverse 
        NekDouble invmag = 1.0/(m_leading_real_evl[0]*m_leading_real_evl[0] + 
                                m_leading_imag_evl[0]*m_leading_imag_evl[0]);
        m_leading_real_evl[0] *= -invmag;
        m_leading_real_evl[0] += realShift;
        m_leading_imag_evl[0] *= invmag;


        m_waveVelocities = solverWave->GetEqu()[0]->UpdateFields();
        m_wavePressure   = solverWave->GetEqu()[0]->GetPressure();


        if( m_sessionVWI->DefinesSolverInfo("INTERFACE")  )
        {
            cout << "Growth =" <<m_leading_real_evl[0]<<endl; 
            cout << "Phase =" <<m_leading_imag_evl[0]<<endl; 
        }
        
    }

    void VortexWaveInteraction::CalcNonLinearWaveForce(void)
    {
        if(m_sessionRoll->DefinesSolverInfo("INTERFACE") )
        {
            static int cnt=0;
            string wavefile  = m_sessionName +".fld"; 
            string movedmesh = m_sessionName + "_advPost_moved.xml";
            string filestreak   = m_sessionName + "_streak.fld";
            char c[16]="";
            sprintf(c,"%d",cnt);               
            char c_alpha[16]="";
            sprintf(c_alpha,"%f",m_alpha[0]);    
            string syscall;
            if( m_sessionVWI->GetSolverInfo("INTERFACE")=="phase" )
            {
                string filePost = m_sessionName + "_advPost.xml";
                syscall = "../../utilities/PostProcessing/Extras/FldCalcBCs  "
                    + filePost +"     "+
                    "meshhalf_pos_Spen_stability_moved.fld  meshhalf_pos_Spen_advPost_moved.fld "
                    +c_alpha +"  > data_alpha0";
                cout<<syscall.c_str()<<endl;
                if(system(syscall.c_str()))
                {
                    ASSERTL0(false,syscall.c_str());
                }
             
                syscall = "cp -f meshhalf_pos_Spen_stability_moved_u_5.bc  "+m_sessionName+"_u_5.bc";  
                cout<<syscall.c_str()<<endl;
                if(system(syscall.c_str()))
                {
                    ASSERTL0(false,syscall.c_str());
                }
                syscall = "cp -f meshhalf_pos_Spen_stability_moved_v_5.bc  "+m_sessionName+"_v_5.bc";  
                cout<<syscall.c_str()<<endl;
                if(system(syscall.c_str()))
                {
                    ASSERTL0(false,syscall.c_str());
                }
            }
            else
            {
                syscall =  "../../utilities/PostProcessing/Extras/FldCalcBCs  "
                    + movedmesh + "  " + wavefile + "  " + filestreak + "   "+c_alpha +"  >  datasub_"+c;
                cout<<syscall.c_str()<<endl;
                if(system(syscall.c_str()))
                {
                    ASSERTL0(false,syscall.c_str());
                }
            }

             
             
            //relaxation for different alpha values? does it make sense?

            //save the wave
            string wave_subalp = m_sessionName + "_wave_subalp_"+c+".fld";
            syscall = "cp -f " + wavefile + "  " + wave_subalp;
            cout<<syscall.c_str()<<endl;
            if(system(syscall.c_str()))
            {
                ASSERTL0(false,syscall.c_str());
            }
            //FileRelaxation(3);
            cnt++;
     
        }
        else
        {
            int npts    = m_waveVelocities[0]->GetPlane(0)->GetNpoints();
            int ncoeffs = m_waveVelocities[0]->GetPlane(0)->GetNcoeffs();
            Array<OneD, NekDouble> val(npts), der1(2*npts);
            Array<OneD, NekDouble> der2 = der1 + npts; 
            Array<OneD, NekDouble> streak;
            static int projectfield = -1;

            if(m_deltaFcnApprox)
            {
                streak = Array<OneD, NekDouble> (npts);
                m_streakField[0]->BwdTrans(m_streakField[0]->GetCoeffs(), streak);
            }
 
            // Set project field to be first field that has a Neumann
            // boundary since this not impose any condition on the vertical boundaries
            // Othersise set to zero. 
            if(projectfield == -1)
            {
                Array<OneD, const SpatialDomains::BoundaryConditionShPtr > BndConds;
                 
                for(int i = 0; i < m_waveVelocities.num_elements(); ++i)
                {
                    BndConds = m_waveVelocities[i]->GetBndConditions();
                    for(int j = 0; j < BndConds.num_elements(); ++j)
                    {
                        if(BndConds[j]->GetBoundaryConditionType() == SpatialDomains::eNeumann)
                        {
                            projectfield = i;
                            break;
                        }
                    }
                    if(projectfield != -1)
                    {
                        break;
                    }
                }
                if(projectfield == -1)
                {
                    cout << "using first field to project non-linear forcing which imposes a Dirichlet condition" << endl;
                    projectfield = 0;
                }
            }
        
            // Shift m_vwiForcing in case of relaxation 
            Vmath::Vcopy(ncoeffs,m_vwiForcing[0],1,m_vwiForcing[2],1);
            Vmath::Vcopy(ncoeffs,m_vwiForcing[1],1,m_vwiForcing[3],1);
        
            // determine inverse of area normalised field. 
            m_wavePressure->GetPlane(0)->BwdTrans(m_wavePressure->GetPlane(0)->GetCoeffs(), m_wavePressure->GetPlane(0)->UpdatePhys());
            m_wavePressure->GetPlane(1)->BwdTrans(m_wavePressure->GetPlane(1)->GetCoeffs(), m_wavePressure->GetPlane(1)->UpdatePhys());
            NekDouble invnorm;

            if(m_useLinfPressureNorm)
            {
                Vmath::Vmul(npts,m_wavePressure->GetPlane(0)->UpdatePhys(),1,m_wavePressure->GetPlane(0)->UpdatePhys(),1,der1,1);
                Vmath::Vvtvp(npts,m_wavePressure->GetPlane(1)->UpdatePhys(),1,m_wavePressure->GetPlane(1)->UpdatePhys(),1,der1,1,der1,1);
                Vmath::Vsqrt(npts,der1,1,der1,1);
                
                NekDouble Linf = Vmath::Vmax(npts,der1,1);
                
                invnorm = 1.0/Linf;
            }
            else
            {
                // Determine normalisation of pressure so that |P|/A = 1
                NekDouble l2;
                l2    = m_wavePressure->GetPlane(0)->L2();
                invnorm  = l2*l2;
                l2    = m_wavePressure->GetPlane(1)->L2();
                invnorm += l2*l2;
                Vmath::Fill(2*npts,1.0,der1,1);
                NekDouble area = m_waveVelocities[0]->GetPlane(0)->PhysIntegral(der1);
                cout << "Area: " << area << endl;
                invnorm = sqrt(area/invnorm);
            }
        
            // Get hold of arrays. 
            m_waveVelocities[0]->GetPlane(0)->BwdTrans(m_waveVelocities[0]->GetPlane(0)->GetCoeffs(),m_waveVelocities[0]->GetPlane(0)->UpdatePhys());
            Array<OneD, NekDouble> u_real = m_waveVelocities[0]->GetPlane(0)->UpdatePhys();
            Vmath::Smul(npts,invnorm,u_real,1,u_real,1);
            m_waveVelocities[0]->GetPlane(1)->BwdTrans(m_waveVelocities[0]->GetPlane(1)->GetCoeffs(),m_waveVelocities[0]->GetPlane(1)->UpdatePhys());
            Array<OneD, NekDouble> u_imag = m_waveVelocities[0]->GetPlane(1)->UpdatePhys();
            Vmath::Smul(npts,invnorm,u_imag,1,u_imag,1);
            m_waveVelocities[1]->GetPlane(0)->BwdTrans(m_waveVelocities[1]->GetPlane(0)->GetCoeffs(),m_waveVelocities[1]->GetPlane(0)->UpdatePhys());
            Array<OneD, NekDouble> v_real = m_waveVelocities[1]->GetPlane(0)->UpdatePhys(); 
            Vmath::Smul(npts,invnorm,v_real,1,v_real,1);
            m_waveVelocities[1]->GetPlane(1)->BwdTrans(m_waveVelocities[1]->GetPlane(1)->GetCoeffs(),m_waveVelocities[1]->GetPlane(1)->UpdatePhys());
            Array<OneD, NekDouble> v_imag = m_waveVelocities[1]->GetPlane(1)->UpdatePhys();
            Vmath::Smul(npts,invnorm,v_imag,1,v_imag,1);

            // normalise wave for output
            Vmath::Smul(2*ncoeffs,invnorm,m_waveVelocities[0]->UpdateCoeffs(),1,m_waveVelocities[0]->UpdateCoeffs(),1);
            Vmath::Smul(2*ncoeffs,invnorm,m_waveVelocities[1]->UpdateCoeffs(),1,m_waveVelocities[1]->UpdateCoeffs(),1);
            Vmath::Smul(2*ncoeffs,invnorm,m_waveVelocities[2]->UpdateCoeffs(),1,m_waveVelocities[2]->UpdateCoeffs(),1);

            // dump field
            {
                Array<OneD, std::string> variables(3);
                Array<OneD, Array<OneD, NekDouble> > outfield(3);
                variables[0] = "u_w"; 
                variables[1] = "v_w"; 
                variables[2] = "w_w"; 
                outfield[0] = m_waveVelocities[0]->UpdateCoeffs();
                outfield[1] = m_waveVelocities[1]->UpdateCoeffs();
                outfield[2] = m_waveVelocities[2]->UpdateCoeffs();
                std::string outname = m_sessionName  + "_wave.fld";
                m_solverRoll->WriteFld(outname, m_waveVelocities[0], outfield, variables);
            }

#if 1
            int ncoeffs_p = m_wavePressure->GetPlane(0)->GetNcoeffs();
            Vmath::Smul(ncoeffs_p,invnorm,m_wavePressure->GetPlane(0)->UpdateCoeffs(),1,m_wavePressure->GetPlane(0)->UpdateCoeffs(),1);
            Vmath::Smul(ncoeffs_p,invnorm,m_wavePressure->GetPlane(1)->UpdateCoeffs(),1,m_wavePressure->GetPlane(1)->UpdateCoeffs(),1);
#else
            m_wavePressure->GetPlane(0)->BwdTrans(m_wavePressure->GetPlane(0)->GetCoeffs(),m_wavePressure->GetPlane(0)->UpdatePhys());
            Vmath::Smul(npts,invnorm,m_wavePressure->GetPlane(0)->UpdatePhys(),1,m_wavePressure->GetPlane(0)->UpdatePhys(),1);
            m_wavePressure->GetPlane(0)->FwdTrans(m_wavePressure->GetPlane(0)->UpdatePhys(),m_wavePressure->GetPlane(0)->UpdateCoeffs());
            
            m_wavePressure->GetPlane(1)->BwdTrans(m_wavePressure->GetPlane(1)->GetCoeffs(),m_wavePressure->GetPlane(1)->UpdatePhys());
            Vmath::Smul(npts,invnorm,m_wavePressure->GetPlane(1)->UpdatePhys(),1,m_wavePressure->GetPlane(1)->UpdatePhys(),1);
            m_wavePressure->GetPlane(1)->FwdTrans(m_wavePressure->GetPlane(1)->UpdatePhys(),m_wavePressure->GetPlane(1)->UpdateCoeffs());
#endif
       
            // Calculate non-linear terms for x and y directions
            // d/dx(u u* + u* u)
            Vmath::Vmul (npts,u_real,1,u_real,1,val,1);
            Vmath::Vvtvp(npts,u_imag,1,u_imag,1,val,1,val,1);
            Vmath::Smul (npts,2.0,val,1,val,1);
            m_waveVelocities[0]->GetPlane(0)->PhysDeriv(0,val,der1);
        
            // d/dy(v u* + v* u)
            Vmath::Vmul (npts,u_real,1,v_real,1,val,1);
            Vmath::Vvtvp(npts,u_imag,1,v_imag,1,val,1,val,1);
            Vmath::Smul (npts,2.0,val,1,val,1);
            m_waveVelocities[0]->GetPlane(0)->PhysDeriv(1,val,der2);
            
            Vmath::Vadd(npts,der1,1,der2,1,der1,1);
#if 1
            m_waveVelocities[projectfield]->GetPlane(0)->FwdTrans(der1,m_vwiForcing[0]);
#else
            m_waveVelocities[0]->GetPlane(0)->FwdTrans_BndConstrained(der1,m_vwiForcing[0]);
#endif
            Vmath::Smul(ncoeffs,-m_waveForceMag[0],m_vwiForcing[0],1,m_vwiForcing[0],1);

            // d/dx(u v* + u* v)
            m_waveVelocities[0]->GetPlane(0)->PhysDeriv(0,val,der1);
        
            // d/dy(v v* + v* v)
            Vmath::Vmul(npts,v_real,1,v_real,1,val,1);
            Vmath::Vvtvp(npts,v_imag,1,v_imag,1,val,1,val,1);
            Vmath::Smul (npts,2.0,val,1,val,1);
            m_waveVelocities[0]->GetPlane(0)->PhysDeriv(1,val,der2);
        
            Vmath::Vadd(npts,der1,1,der2,1,der1,1);

#if 1
            m_waveVelocities[projectfield]->GetPlane(0)->FwdTrans(der1,m_vwiForcing[1]);
#else
            m_waveVelocities[0]->GetPlane(0)->FwdTrans_BndConstrained(der1,m_vwiForcing[1]);
#endif

            Vmath::Smul(ncoeffs,-m_waveForceMag[0],m_vwiForcing[1],1,m_vwiForcing[1],1);
        
            //by default the symmetrization is on
            bool symm=true;
            m_sessionVWI->MatchSolverInfo("Symmetrization","True",symm,true);
#if 0
            if(symm== true )
            {
                
                // Symmetrise forcing
                //-> Get coordinates 
                Array<OneD, NekDouble> coord(2);
                Array<OneD, NekDouble> coord_x(npts);
                Array<OneD, NekDouble> coord_y(npts);
                
                //-> Impose symmetry (x -> -x + Lx/2, y-> -y) on coordinates
                m_waveVelocities[0]->GetPlane(0)->GetCoords(coord_x,coord_y);
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
                    Eid[i] = m_waveVelocities[0]->GetPlane(0)->GetExpIndex(coord,1e-6);
                }
                
                // Interpolate field 0 
                m_waveVelocities[0]->GetPlane(0)->BwdTrans_IterPerExp(m_vwiForcing[0],der1);
                for(i = 0; i < npts; ++i)
                {
                    physoffset = m_waveVelocities[0]->GetPlane(0)->GetPhys_Offset(Eid[i]);
                    coord[0] = coord_x[i];
                    coord[1] = coord_y[i];
                    der2 [i] = m_waveVelocities[0]->GetPlane(0)->GetExp(Eid[i])->PhysEvaluate(coord, der1 + physoffset);
                }
                //-> Average field 0 
                Vmath::Vsub(npts,der1,1,der2,1,der2,1);
                Vmath::Smul(npts,0.5,der2,1,der2,1);
#if 1
                m_waveVelocities[projectfield]->GetPlane(0)->FwdTrans(der2,m_vwiForcing[0]);
#else
                m_waveVelocities[0]->GetPlane(0)->FwdTrans_BndConstrained(der2, m_vwiForcing[0]);
#endif
                
                //-> Interpoloate field 1
                m_waveVelocities[0]->GetPlane(0)->BwdTrans_IterPerExp(m_vwiForcing[1],der1);
                for(i = 0; i < npts; ++i)
                {
                    physoffset = m_waveVelocities[0]->GetPlane(0)->GetPhys_Offset(Eid[i]);
                    coord[0] = coord_x[i];
                    coord[1] = coord_y[i];
                    der2[i]  = m_waveVelocities[0]->GetPlane(0)->GetExp(Eid[i])->PhysEvaluate(coord,                                                                         der1 + physoffset);
                }
                
                //-> Average field 1
                Vmath::Vsub(npts,der1,1,der2,1,der2,1);
                Vmath::Smul(npts,0.5,der2,1,der2,1);
#if 1
                m_waveVelocities[projectfield]->GetPlane(0)->FwdTrans(der2,m_vwiForcing[1]);
#else
                m_waveVelocities[0]->GetPlane(0)->FwdTrans_BndConstrained(der2, m_vwiForceing[1]);
#endif
            }
#else
            int i;
            if(symm== true )
            {
                cout<<"symmetrization is active"<<endl;              
                static Array<OneD, int> index = GetReflectionIndex();
                
                m_waveVelocities[0]->GetPlane(0)->BwdTrans_IterPerExp(m_vwiForcing[0],der1);
                for(i = 0; i < npts; ++i)
                {
                    if(index[i] != -1)
                    {
                        val[i] = 0.5*(der1[i] - der1[index[i]]);
                    }
                }
#if 1 
                m_waveVelocities[projectfield]->GetPlane(0)->FwdTrans(val,m_vwiForcing[0]);
#else
                m_waveVelocities[0]->GetPlane(0)->FwdTrans_BndConstrained(val, m_vwiForcing[0]);
#endif

                m_waveVelocities[0]->GetPlane(0)->BwdTrans_IterPerExp(m_vwiForcing[1],der2);
                for(i = 0; i < npts; ++i)
                {
                    if(index[i] != -1)
                    {
                        val[i] = 0.5*(der2[i] - der2[index[i]]);
                    }
                }        
#if 1 
                m_waveVelocities[projectfield]->GetPlane(0)->FwdTrans(val,m_vwiForcing[1]);
#else
                m_waveVelocities[0]->GetPlane(0)->FwdTrans_BndConstrained(val, m_vwiForcing[1]);
#endif
            }


            Vmath::Vmul(npts,der1,1,der1,1,val,1);
            Vmath::Vvtvp(npts,der2,1,der2,1,val,1,val,1);
            Vmath::Vsqrt(npts,val,1,val,1);
            cout << "F_Linf: " <<  Vmath::Vmax(npts,val,1) << endl;

#endif

            if(m_vwiRelaxation)
            {
                 Vmath::Smul(ncoeffs,1.0-m_vwiRelaxation,
                        m_vwiForcing[0],1,m_vwiForcing[0],1);
                 Vmath::Svtvp(ncoeffs,m_vwiRelaxation,m_vwiForcing[2],1,
                         m_vwiForcing[0],1,m_vwiForcing[0],1);
            
                 Vmath::Smul(ncoeffs,1.0-m_vwiRelaxation,
                        m_vwiForcing[1],1,m_vwiForcing[1],1);
                 Vmath::Svtvp(ncoeffs,m_vwiRelaxation,m_vwiForcing[3],1,
                         m_vwiForcing[1],1,m_vwiForcing[1],1);
            }

            if(m_deltaFcnApprox)
            {
                 for(int j = 0; j < 2; ++j)
                 {
                
                      m_waveVelocities[projectfield]->GetPlane(0)->BwdTrans(m_vwiForcing[j],der1);
                      for(int i = 0; i < npts; ++i)
                      {
                           der1[i] *= exp(-streak[i]*streak[i]/m_deltaFcnDecay);
                      }
                      m_waveVelocities[projectfield]->GetPlane(0)->FwdTrans(der1,m_vwiForcing[j]);
                 }
            }
            
            
            // dump output
            Array<OneD, std::string> variables(4);
            Array<OneD, Array<OneD, NekDouble> > outfield(4);
            variables[0] = "u";  variables[1] = "v"; 
            variables[2] = "pr"; variables[3] = "pi";
            outfield[0]  = m_vwiForcing[0];
            outfield[1]  = m_vwiForcing[1];
            Array<OneD,NekDouble> soln(npts,0.0);
            m_wavePressure->GetPlane(0)->BwdTrans(m_wavePressure->GetPlane(0)->GetCoeffs(),m_wavePressure->GetPlane(0)->UpdatePhys());
            outfield[2] = Array<OneD,NekDouble>(ncoeffs);
            m_waveVelocities[0]->GetPlane(0)->FwdTrans_IterPerExp(m_wavePressure->GetPlane(0)->GetPhys(),outfield[2]);
            m_wavePressure->GetPlane(1)->BwdTrans(m_wavePressure->GetPlane(1)->GetCoeffs(),m_wavePressure->GetPlane(1)->UpdatePhys());
            
            Vmath::Vmul(npts,m_wavePressure->GetPlane(0)->UpdatePhys(),1,
                        m_wavePressure->GetPlane(0)->UpdatePhys(),1,val,1);
            Vmath::Vvtvp(npts,m_wavePressure->GetPlane(1)->UpdatePhys(),1,
                         m_wavePressure->GetPlane(1)->UpdatePhys(),1,val,1,val,1);
            cout << "int P^2: " << m_wavePressure->GetPlane(0)->PhysIntegral(val) << endl;
            Vmath::Vsqrt(npts,val,1,val,1);
            cout << "PLinf: " <<  Vmath::Vmax(npts,val,1) << endl;

            outfield[3] = Array<OneD,NekDouble>(ncoeffs);
            m_waveVelocities[1]->GetPlane(0)->FwdTrans_IterPerExp(m_wavePressure->GetPlane(1)->GetPhys(),outfield[3]);
            
            std::string outname = m_sessionName  + ".vwi";
            
            m_solverRoll->WriteFld(outname, m_waveVelocities[0]->GetPlane(0), outfield, variables);

        }
    }
    
    void VortexWaveInteraction::CalcL2ToLinfPressure(void)
    {

        ExecuteWave();

        m_wavePressure->GetPlane(0)->BwdTrans(m_wavePressure->GetPlane(0)->GetCoeffs(),
                                              m_wavePressure->GetPlane(0)->UpdatePhys());
        m_wavePressure->GetPlane(1)->BwdTrans(m_wavePressure->GetPlane(1)->GetCoeffs(),
                                              m_wavePressure->GetPlane(1)->UpdatePhys());

        int npts    = m_waveVelocities[0]->GetPlane(0)->GetNpoints();
        NekDouble Linf;
        Array<OneD, NekDouble> val(2*npts,0.0);
        
        Vmath::Vmul(npts,m_wavePressure->GetPlane(0)->UpdatePhys(),1,m_wavePressure->GetPlane(0)->UpdatePhys(),1,val,1);
        Vmath::Vvtvp(npts,m_wavePressure->GetPlane(1)->UpdatePhys(),1,m_wavePressure->GetPlane(1)->UpdatePhys(),1,val,1,val,1);
        cout << "int P^2 " << m_wavePressure->GetPlane(0)->PhysIntegral(val) << endl;
        Vmath::Vsqrt(npts,val,1,val,1);
        
               
        Linf = Vmath::Vmax(npts,val,1);
        cout << "Linf: " << Linf << endl;

        NekDouble l2,norm;
        l2    = m_wavePressure->GetPlane(0)->L2();
        norm  = l2*l2;
        l2    = m_wavePressure->GetPlane(1)->L2();
        norm += l2*l2;


        Vmath::Fill(npts,1.0,val,1);
        NekDouble area = m_waveVelocities[0]->GetPlane(0)->PhysIntegral(val);

        l2 = sqrt(norm/area);

        cout << "L2:   " << l2 << endl;
        
        cout << "Ratio Linf/L2: "<< Linf/l2 << endl;
    }

    void VortexWaveInteraction::SaveFile(string file, string dir, int n)
    {
        static map<string,int> opendir;

        if(opendir.find(dir) == opendir.end())
        {
            // make directory and presume will fail if it already exists
            string mkdir = "mkdir " + dir;
            system(mkdir.c_str());

            opendir[dir] = 1;
        }
        
        string savefile = dir + "/" + file + "." + boost::lexical_cast<std::string>(n);
        string syscall  = "cp -f "  + file + " " + savefile; 

        if(system(syscall.c_str()))
        {
            ASSERTL0(false,syscall.c_str());
        }

     }


    void VortexWaveInteraction::MoveFile(string file, string dir, int n)
    {
        static map<string,int> opendir;

        if(opendir.find(dir) == opendir.end())
        {
            // make directory and presume will fail if it already exists
            string mkdir = "mkdir " + dir;
            system(mkdir.c_str());
            opendir[dir] = 1;
        }
        
        string savefile = dir + "/" + file + "." + boost::lexical_cast<std::string>(n);
        string syscall  = "mv -f "  + file + " " + savefile; 

        if(system(syscall.c_str()))
        {
            ASSERTL0(false,syscall.c_str());
        }
     }

     void VortexWaveInteraction::CopyFile(string file1end, string file2end)
     {
         string cpfile1   = m_sessionName + file1end;
         string cpfile2   = m_sessionName + file2end;
         string syscall  = "cp -f "  + cpfile1 + " " + cpfile2; 

         if(system(syscall.c_str()))
         {
             ASSERTL0(false,syscall.c_str());
         }
     }

     void VortexWaveInteraction::AppendEvlToFile(string file, int n)
     {
         FILE *fp;
         fp = fopen(file.c_str(),"a");
         fprintf(fp, "%d: %lf %16.12le  %16.12le\n",n, m_alpha[0], m_leading_real_evl[0],m_leading_imag_evl[0]);
         fclose(fp);
     }

     void VortexWaveInteraction::AppendEvlToFile(string file, NekDouble WaveForceMag)
     {
         FILE *fp;
         fp = fopen(file.c_str(),"a");
         fprintf(fp, "%lf %lf %16.12le  %16.12le\n",WaveForceMag, m_alpha[0], m_leading_real_evl[0],m_leading_imag_evl[0]);
         fclose(fp);
     }

    void VortexWaveInteraction::SaveLoopDetails(std::string SaveDir, int i)

    {
        // Save NS restart file
        SaveFile(m_sessionName + ".rst",SaveDir,i+1);
        // Save Streak Solution
        SaveFile(m_sessionName + "_streak.fld",SaveDir,i);
        // Save Wave solution output
        SaveFile(m_sessionName + ".evl",SaveDir,i);
        SaveFile(m_sessionName + "_eig_0",SaveDir,i);
        // Save new field file of eigenvalue
        SaveFile(m_sessionName + ".fld",SaveDir,i);
        if(!(m_sessionVWI->DefinesSolverInfo("INTERFACE")))
        {
            // Save new forcing file
            SaveFile(m_sessionName + ".vwi",SaveDir,i+1);
        }

    }

    void VortexWaveInteraction::ExecuteLoop(bool CalcWaveForce)
    {


	//the global info has to be written in the 
	//roll session file to use the interface loop
        if(m_sessionRoll->DefinesSolverInfo("INTERFACE"))
	{
             static int cnt=0;     
             bool skiprollstreak=false;
             if(cnt==0 && m_sessionVWI->GetParameter("rollstreakfromit")==1)
             {
                  skiprollstreak =true;
                  cout<<"skip roll-streak at the first iteration"<<endl;
             }
             
             
             if(skiprollstreak != true)
             {

                 LibUtilities::NekManager<MultiRegions::GlobalLinSysKey,MultiRegions::GlobalLinSys>::EnableManagement("GlobalLinSys");
                 ExecuteRoll();
                 LibUtilities::NekManager<MultiRegions::GlobalLinSysKey,MultiRegions::GlobalLinSys>::DisableManagement("GlobalLinSys");
#ifndef _WIN32
   	         sleep(3);
#endif
                 ExecuteStreak();
#ifndef _WIN32
	         sleep(3);
#endif
             }
              
             string syscall;
             char c[16]="";
             string movedmesh = m_sessionName + "_advPost_moved.xml";
             string movedinterpmesh = m_sessionName + "_interp_moved.xml";
             //rewrite the Rollsessionfile (we start from the waleffe forcing)
             //string meshbndjumps = m_sessionName +"_bndjumps.xml";             
             //if(cnt==0)
             //{
                 //take the conditions tag from meshbndjumps and copy into 
                 // the rolls session file
             //}


             sprintf(c,"%d",cnt);  
             //save old roll solution
             string oldroll = m_sessionName +"_roll_"+c +".fld";    
             syscall = "cp -f " + m_sessionName+"-Base.fld" + "  " + oldroll;
             cout<<syscall.c_str()<<endl;
             if(system(syscall.c_str()))
             {
                  ASSERTL0(false,syscall.c_str());
             } 
             //define file names
             string filePost   = m_sessionName + "_advPost.xml";
             string filestreak   = m_sessionName + "_streak.fld";
             string filewave    = m_sessionName + "_wave.fld";
             string filewavepressure = m_sessionName + "_wave_p_split.fld";
             string fileinterp = m_sessionName + "_interp.xml";
             string interpstreak = m_sessionName +"_interpstreak_"+ c +".fld";  
             string interwavepressure  = m_sessionName +"_wave_p_split_interp_"+ c +".fld";
             char alpchar[16]="";
cout<<"alpha = "<<m_alpha[0]<<endl;
             sprintf(alpchar, "%f", m_alpha[0]);


             if( m_sessionVWI->GetSolverInfo("INTERFACE")!="phase" )
             {
                 cout<<"zerophase"<<endl;

                 syscall  = "../../utilities/PostProcessing/Extras/MoveMesh  "
                             + filePost +"  "+ filestreak +"  "+ fileinterp + "   "+ alpchar; 

                 cout<<syscall.c_str()<<endl;
                 if(system(syscall.c_str()))
                 {
                      ASSERTL0(false,syscall.c_str());
                 }

                 //move the advPost mesh (remark update alpha!!!)
                 syscall  =  "../../utilities/PostProcessing/Extras/MoveMesh  "
                       + filePost + "  " + filestreak + "  " + filePost + "    "+ alpchar;
                 cout<<syscall.c_str()<<endl;
                 if(system(syscall.c_str()))
                 {
                     ASSERTL0(false,syscall.c_str());
                 }



                 //save oldstreak
                 string oldstreak = m_sessionName +"_streak_"+ c +".fld";            
        	 syscall = "cp -f " + filestreak + "  " + oldstreak;
                 cout<<syscall.c_str()<<endl;
                 if(system(syscall.c_str()))
                 {
                      ASSERTL0(false,syscall.c_str());
                 } 

    	         //interpolate the streak field into the new mesh
                 string movedmesh = m_sessionName + "_advPost_moved.xml";
                 string movedinterpmesh = m_sessionName + "_interp_moved.xml";

                 //create the interp streak             

                 syscall  =  "../../utilities/PostProcessing/Extras/FieldToField  "
                      + fileinterp + "  " + filestreak + "  " + movedinterpmesh + "  " 
	              + interpstreak;

                 cout<<syscall.c_str()<<endl;
                 if(system(syscall.c_str()))
                 {
                      ASSERTL0(false,syscall.c_str());
                 } 


                 //save the old mesh     
                 string meshfile = m_sessionName + ".xml";                  
                 string meshold = m_sessionName +"_"+ c +".xml";
     	         syscall = "cp -f " + meshfile + "  " + meshold;
                 cout<<syscall.c_str()<<endl;
                 if(system(syscall.c_str()))
                 {
                     ASSERTL0(false,syscall.c_str());
                 }  

                 //overwriting the meshfile with the new mesh
	         syscall = "cp -f " + movedmesh + "  " + meshfile;
                 cout<<syscall.c_str()<<endl;
                 if(system(syscall.c_str()))
                 {
                     ASSERTL0(false,syscall.c_str());
                 }

                 //overwriting the streak file!!          
                 syscall = "cp -f " + interpstreak + "  " + filestreak;
                 cout<<syscall.c_str()<<endl;
                 if(system(syscall.c_str()))
                 {
                     ASSERTL0(false,syscall.c_str());
                 } 

                 //calculate the wave
                 ExecuteWave();

                 //save the wave field:
                 string oldwave = m_sessionName +"_wave_"+c +".fld";    
	         syscall = "cp -f " + m_sessionName+".fld" + "  " + oldwave;
                 cout<<syscall.c_str()<<endl;
                 if(system(syscall.c_str()))
                 {
                      ASSERTL0(false,syscall.c_str());
                 } 

                 //save old jump conditions:
                 string ujump = m_sessionName+"_u_5.bc";
	         syscall = "cp -f " + ujump + "  " + m_sessionName+"_u_5.bc_"+c;
                 cout<<syscall.c_str()<<endl;
                 if(system(syscall.c_str()))
                 {
                      ASSERTL0(false,syscall.c_str());
                 }              

                 string vjump = m_sessionName+"_v_5.bc";
     	         syscall = "cp -f " + vjump + "  " + m_sessionName+"_v_5.bc_"+c;
                 cout<<syscall.c_str()<<endl;
                 if(system(syscall.c_str()))
                 {
                     ASSERTL0(false,syscall.c_str());
                 }    
                 cnt++;




                  //use relaxation
                  if(GetVWIIterationType()!=eFixedWaveForcingWithSubIterationOnAlpha 
                    )
                  {
                       // the critical layer should be the bnd region 3
                       int reg =3;
                       //FileRelaxation(reg);
                  }
                  char c1[16]="";
                  sprintf(c1,"%d",cnt);   
                  //calculate the jump conditions
                  string wavefile  = m_sessionName +".fld"; 
                  syscall =  "../../utilities/PostProcessing/Extras/FldCalcBCs  "
                        + movedmesh + "  " + wavefile + "  " + interpstreak + ">  data"+c1;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                       ASSERTL0(false,syscall.c_str());
                  }

                  //move the new name_interp_moved.xml into name_interp.xml
	          syscall = "cp -f " + movedinterpmesh + "  " + fileinterp;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                       ASSERTL0(false,syscall.c_str());
                  } 
                  //move the new name_advPost_moved.xml into name_advPost.xml
	          syscall = "cp -f " + movedmesh + "  " + filePost;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                       ASSERTL0(false,syscall.c_str());
                  } 
 

             }
             else if(  m_sessionVWI->GetSolverInfo("INTERFACE")=="phase" )
             {    
cout<<"phase"<<endl;
                  //determine cr:
                  NekDouble cr;
                  string cr_str;
                  stringstream st;

                  //calculate the wave
                  ExecuteWave();
   
                  //save oldstreak
                  string oldstreak = m_sessionName +"_streak_"+ c +".fld";            
	          syscall = "cp -f " + filestreak + "  " + oldstreak;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                       ASSERTL0(false,syscall.c_str());
                  } 

                  //save wave
	          syscall = "cp -f " + m_sessionName+".fld" + "  " + filewave;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                       ASSERTL0(false,syscall.c_str());
                  } 
                  //save the old mesh     
                  string meshfile = m_sessionName + ".xml";                  
                  string meshold = m_sessionName +"_"+ c +".xml";
	          syscall = "cp -f " + meshfile + "  " + meshold;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                       ASSERTL0(false,syscall.c_str());
                  }  

                  //save the oldwave field:
                  string oldwave = m_sessionName +"_wave_"+c +".fld";    
    	          syscall = "cp -f " + m_sessionName+".fld" + "  " + oldwave;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                       ASSERTL0(false,syscall.c_str());
                  }       

                  //save old jump conditions:
                  string ujump = m_sessionName+"_u_5.bc";
	          syscall = "cp -f " + ujump + "  " + m_sessionName+"_u_5.bc_"+c;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                      ASSERTL0(false,syscall.c_str());
                  }              

                  string vjump = m_sessionName+"_v_5.bc";
     	          syscall = "cp -f " + vjump + "  " + m_sessionName+"_v_5.bc_"+c;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                     ASSERTL0(false,syscall.c_str());
                  }
                  cnt++;


                  cr = m_leading_imag_evl[0]/m_alpha[0];
                  st << cr; 
                  cr_str = st.str();
cout<<"cr="<<cr_str<<endl;
                  //NB -g or NOT!!!
                  //move the mesh around the critical layer    
                  syscall  = "../../utilities/PostProcessing/Extras/MoveMesh  "
                               + filePost +"  "+ filestreak +"  "+ fileinterp + "   "+ alpchar
                               +"      "+cr_str; 
  
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                       ASSERTL0(false,syscall.c_str());
                  }
                  //NB -g or NOT!!!
                  //move the advPost mesh (remark update alpha!!!)
                  syscall  =  "../../utilities/PostProcessing/Extras/MoveMesh  "
                        + filePost + "  " + filestreak + "  " + filePost + "    "+ alpchar
                               +"      "+cr_str; 
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                      ASSERTL0(false,syscall.c_str());
                  }

                  //interp streak into the new mesh
                  syscall  =  "../../utilities/PostProcessing/Extras/FieldToField  "
                        + fileinterp + "  " + filestreak + "  " + movedinterpmesh + "  " 
	                + interpstreak;

                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                      ASSERTL0(false,syscall.c_str());
                  }

                  //split wave sol
                  syscall  =  "../../utilities/PostProcessing/Extras/SplitFld  "
                        + filePost + "  " + filewave;

                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                      ASSERTL0(false,syscall.c_str());
                  }                
                  //interp wave
                  syscall  =  "../../utilities/PostProcessing/Extras/FieldToField  "
                        + filePost + "  " + filewavepressure + "  " + movedmesh + "  " 
	                + interwavepressure;

                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                      ASSERTL0(false,syscall.c_str());
                  }




                  //use relaxation
                  if(GetVWIIterationType()!=eFixedWaveForcingWithSubIterationOnAlpha 
                    )
                  {
                      // the critical layer should be the bnd region 3
                      int reg =3;
                      //FileRelaxation(reg);
                  }
                  char c1[16]="";
                  sprintf(c1,"%d",cnt);   

                  //cp wavepressure to m_sessionName.fld(to get
                  // the right bcs names using FldCalcBCs)
                  syscall = "cp -f "+ interwavepressure +"  "+m_sessionName+".fld";
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                      ASSERTL0(false,syscall.c_str());
                  }

                  //calculate the jump conditions
                  //NB -g or NOT!!!
                  syscall =  "../../utilities/PostProcessing/Extras/FldCalcBCs  "
                        + movedmesh + "  " +m_sessionName+".fld"  + "  " 
                        + interpstreak + ">  data"+c1;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                      ASSERTL0(false,syscall.c_str());
                  }


                  //overwriting the meshfile with the new mesh
	          syscall = "cp -f " + movedmesh + "  " + meshfile;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                      ASSERTL0(false,syscall.c_str());
                  }

                  //overwriting the streak file!!    (maybe is useless)      
                  syscall = "cp -f " + interpstreak + "  " + filestreak;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                    ASSERTL0(false,syscall.c_str());
                  } 
                  //move the new name_interp_moved.xml into name_interp.xml
	          syscall = "cp -f " + movedinterpmesh + "  " + fileinterp;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                       ASSERTL0(false,syscall.c_str());
                  } 
                  //move the new name_advPost_moved.xml into name_advPost.xml
	          syscall = "cp -f " + movedmesh + "  " + filePost;
                  cout<<syscall.c_str()<<endl;
                  if(system(syscall.c_str()))
                  {
                       ASSERTL0(false,syscall.c_str());
                  }
             }
	}
	else
	{
            LibUtilities::NekManager<MultiRegions::GlobalLinSysKey,MultiRegions::GlobalLinSys>::EnableManagement("GlobalLinSys");
            ExecuteRoll();
            LibUtilities::NekManager<MultiRegions::GlobalLinSysKey,MultiRegions::GlobalLinSys>::DisableManagement("GlobalLinSys");

#ifndef _WIN32
	    sleep(3);
#endif
            ExecuteStreak();
#ifndef _WIN32
            sleep(3);
#endif
            
            if(m_moveMeshToCriticalLayer)
            {
                string syscall;
                char alpchar[16]="";
                sprintf(alpchar, "%f", m_alpha[0]);

                string filePost          = m_sessionName + "_advPost.xml";
                string filestreak        = m_sessionName + "_streak.fld";
                string filewave          = m_sessionName + "_wave.fld";
                string filewavepressure  = m_sessionName + "_wave_p_split.fld";
                string fileinterp        = m_sessionName + "_interp.xml";
                string interpstreak      = m_sessionName +"_interpstreak.fld";  
                string interwavepressure = m_sessionName +"_wave_p_split_interp.fld";
                syscall  = "../../utilities/PostProcessing/Extras/MoveMesh  "
                    + filePost +"  "+ filestreak +"  "+ fileinterp + "   "+ alpchar; 
                
                cout<<syscall.c_str()<<endl;
                if(system(syscall.c_str()))
                {
                    ASSERTL0(false,syscall.c_str());
                }
                
                //move the advPost mesh (remark update alpha!!!)
                syscall  =  "../../utilities/PostProcessing/Extras/MoveMesh  "
                    + filePost + "  " + filestreak + "  " + filePost + "    "+ alpchar;
                cout<<syscall.c_str()<<endl;
                if(system(syscall.c_str()))
                {
                    ASSERTL0(false,syscall.c_str());
                }
                
                //save oldstreak
                string oldstreak = m_sessionName +"_streak.fld";
                syscall = "cp -f " + filestreak + "  " + oldstreak;
                cout<<syscall.c_str()<<endl;
                if(system(syscall.c_str()))
                {
                    ASSERTL0(false,syscall.c_str());
                } 
                
                //interpolate the streak field into the new mesh
                string movedmesh = m_sessionName + "_advPost_moved.xml";
                string movedinterpmesh = m_sessionName + "_interp_moved.xml";
                
                //create the interp streak             
                
                syscall  =  "../../utilities/PostProcessing/Extras/FieldToField  "
                    + fileinterp + "  " + filestreak + "  " + movedinterpmesh 
                    + "  "  + interpstreak;
                
                cout<<syscall.c_str()<<endl;
                if(system(syscall.c_str()))
                {
                    ASSERTL0(false,syscall.c_str());
                } 
                
                //save the old mesh     
                string meshfile = m_sessionName + ".xml";                  
                string meshold = m_sessionName + ".xml";
                syscall = "cp -f " + meshfile + "  " + meshold;
                cout<<syscall.c_str()<<endl;
                if(system(syscall.c_str()))
                {
                    ASSERTL0(false,syscall.c_str());
                }  
                
                //overwriting the meshfile with the new mesh
                syscall = "cp -f " + movedmesh + "  " + meshfile;
                cout<<syscall.c_str()<<endl;
                if(system(syscall.c_str()))
                {
                    ASSERTL0(false,syscall.c_str());
                }
                
                //overwriting the streak file!!          
                syscall = "cp -f " + interpstreak + "  " + filestreak;
                cout<<syscall.c_str()<<endl;
                if(system(syscall.c_str()))
                {
                    ASSERTL0(false,syscall.c_str());
                } 
            }
            
            ExecuteWave();

            if(CalcWaveForce)
            {
                CalcNonLinearWaveForce();
            }
	}
    }

    bool VortexWaveInteraction::CheckEigIsStationary(bool reset)
    {
        static NekDouble previous_real_evl = -1.0; 
        static NekDouble previous_imag_evl = -1.0; 
        static int min_iter = 0;
        
        if(reset)
        {
            previous_real_evl = -1.0;
            min_iter = 0;
        }
        
        if(previous_real_evl == -1.0 || min_iter < m_minInnerIterations)
        {
            previous_real_evl = m_leading_real_evl[0];
            previous_imag_evl = m_leading_imag_evl[0];
            min_iter++;
            return false;
        }

        cout << "Growth tolerance: " << fabs((m_leading_real_evl[0] - previous_real_evl)/m_leading_real_evl[0]) << endl; 
        cout << "Phase tolerance: " << fabs((m_leading_imag_evl[0] - previous_imag_evl)/m_leading_imag_evl[0]) << endl; 

        // See if real and imaginary growth have converged to with m_eigRelTol
        if((fabs((m_leading_real_evl[0] - previous_real_evl)/m_leading_real_evl[0]) < m_eigRelTol)||(fabs(m_leading_real_evl[0]) < 1e-6))
        {
            previous_real_evl = m_leading_real_evl[0];
            previous_imag_evl = m_leading_imag_evl[0];
            if((fabs((m_leading_imag_evl[0] - previous_imag_evl)/m_leading_imag_evl[0]) < m_eigRelTol)||(fabs(m_leading_imag_evl[0]) < 1e-6))
            {
                return true;
            }
        }
        else
        {
            if(fabs(m_leading_imag_evl[0]) > 1e-2)
            {
                cout << "Warning: imaginary eigenvalue is greater than 1e-2" << endl;
            }
            previous_real_evl = m_leading_real_evl[0];
            previous_imag_evl = m_leading_imag_evl[0];
            return false;
        }
        return false;
    }

    // Check to see if leading eigenvalue is within tolerance defined
    // in m_neutralPointTol
    bool VortexWaveInteraction::CheckIfAtNeutralPoint(void)
    {
        bool returnval = false;
        if(m_sessionRoll->DefinesSolverInfo("INTERFACE"))
        {
              if(m_sessionVWI->GetSolverInfo("INTERFACE")=="phase")
              {
                  if( abs(m_leading_real_evl[0]) < 1e-4 )
                  {
                        returnval = true;
                  }
              }
              else
              {
                  if( abs(m_leading_real_evl[0]) < 1e-4  &&  abs(m_leading_imag_evl[0]) <2e-6 )
                  {
                        returnval = true;
                  }
              }
              
        }
        else
        {
               if((m_leading_real_evl[0]*m_leading_real_evl[0] + 
                  m_leading_imag_evl[0]*m_leading_imag_evl[0]) < 
                  m_neutralPointTol*m_neutralPointTol)
               {
                    returnval = true;
               }
        }
        if(fabs(m_leading_imag_evl[0]) > 1e-2)
        {
            cout << "Warning: imaginary eigenvalue is greater than 1e-2" << endl;
        }

        return returnval;
    }
    
    // Similar routine to UpdateAlpha 

    void VortexWaveInteraction::UpdateWaveForceMag(int outeriter)
    {
        NekDouble wavef_new;


        if(outeriter == 1)
        {
            m_waveForceMag[1] = m_waveForceMag[0];
            if(m_leading_real_evl[0] > 0.0)
            {
                wavef_new = m_waveForceMag[0] - m_waveForceMagStep;
            }
            else
            {
                wavef_new = m_waveForceMag[0] + m_waveForceMagStep;
            }
        }
        else
        {
            int i;
            int nstore = (m_waveForceMag.num_elements() < outeriter)? m_waveForceMag.num_elements(): outeriter;
            Array<OneD, NekDouble> WaveForce(nstore);
            Array<OneD, NekDouble> Growth(nstore);
            
            Vmath::Vcopy(nstore,m_waveForceMag,1,WaveForce,1);
            Vmath::Vcopy(nstore,m_leading_real_evl,1,Growth,1);
            
            // Sort WaveForce Growth values; 
            double store;
            int k;
            for(i = 0; i < nstore; ++i)
            {
                k = Vmath::Imin(nstore-i,&WaveForce[i],1);
                
                store    = WaveForce[i]; 
                WaveForce[i] = WaveForce[i+k];
                WaveForce[i+k] = store;

                store     = Growth[i];
                Growth[i] = Growth[i+k];
                Growth[i+k] = store; 
            }

            // See if we have any values that cross zero
            for(i = 0; i < nstore-1; ++i)
            {
                if(Growth[i]*Growth[i+1] < 0.0)
                {
                    break;
                }
            }
            
            if(i != nstore-1)
            {
                if(nstore == 2)
                {
                    wavef_new = (WaveForce[0]*Growth[1] - WaveForce[1]*Growth[0])/(Growth[1]-Growth[0]);
                }
                else
                {
                    // use a quadratic fit and step through 10000 points
                    // to find zero.
                    int     j; 
                    int     nsteps = 10000;
                    int     idx = (i == 0)?1:i;
                    double  da = WaveForce[idx+1] - WaveForce[idx-1];
                    double  gval_m1 = Growth[idx-1],a,gval;
                    double  c1 = Growth[idx-1]/(WaveForce[idx-1]-WaveForce[idx])/
                        (WaveForce[idx-1]-WaveForce[idx+1]);
                    double  c2 = Growth[idx]/(WaveForce[idx]-WaveForce[idx-1])/
                        (WaveForce[idx]-WaveForce[idx+1]);
                    double  c3 = Growth[idx+1]/(WaveForce[idx+1]-WaveForce[idx-1])/
                        (WaveForce[idx+1]-WaveForce[idx]);
                    
                    for(j = 1; j < nsteps+1; ++j)
                    {
                        a = WaveForce[i] + j*da/(double) nsteps;
                        gval = c1*(a - WaveForce[idx  ])*(a - WaveForce[idx+1]) 
                            +  c2*(a - WaveForce[idx-1])*(a - WaveForce[idx+1])
                            +  c3*(a - WaveForce[idx-1])*(a - WaveForce[idx]);
                        
                        if(gval*gval_m1 < 0.0)
                        {
                            wavef_new = ((a+da/(double)nsteps)*gval - a*gval_m1)/
                                (gval - gval_m1);
                            break;
                        }
                    }
                }
            }
            else // step backward/forward 
            {
                if(Growth[i] > 0.0)
                {
                    wavef_new = m_waveForceMag[0] - m_waveForceMagStep;
                }
                else
                {
                    wavef_new = m_waveForceMag[0] + m_waveForceMagStep;
                }
            }
        }
        
        for(int i = m_waveForceMag.num_elements()-1; i > 0; --i)
        {
            m_waveForceMag[i] = m_waveForceMag[i-1];
            m_leading_real_evl[i] = m_leading_real_evl[i-1];
            m_leading_imag_evl[i] = m_leading_imag_evl[i-1];
        }

        m_waveForceMag[0] = wavef_new;
        
    }

    void VortexWaveInteraction::UpdateDAlphaDWaveForceMag(NekDouble alpha_init)
    {
        m_dAlphaDWaveForceMag = (m_alpha[0]-alpha_init)/m_waveForceMagStep;
    }

    void VortexWaveInteraction::UpdateAlpha(int outeriter)
    {
        NekDouble alp_new;


        if(outeriter == 1)
        {
            m_alpha[1] = m_alpha[0];
            if(m_leading_real_evl[0] > 0.0)
            {
                alp_new = m_alpha[0] + m_alphaStep;
            }
            else
            {
                alp_new = m_alpha[0] - m_alphaStep;
            }
        }
        else
        {
            int i;
            int nstore = (m_alpha.num_elements() < outeriter)? m_alpha.num_elements(): outeriter;
            Array<OneD, NekDouble> Alpha(nstore);
            Array<OneD, NekDouble> Growth(nstore);
            
            Vmath::Vcopy(nstore,m_alpha,1,Alpha,1);
            Vmath::Vcopy(nstore,m_leading_real_evl,1,Growth,1);
            
            // Sort Alpha Growth values; 
            double store;
            int k;
            for(i = 0; i < nstore; ++i)
            {
                k = Vmath::Imin(nstore-i,&Alpha[i],1);
                
                store    = Alpha[i]; 
                Alpha[i] = Alpha[i+k];
                Alpha[i+k] = store;

                store     = Growth[i];
                Growth[i] = Growth[i+k];
                Growth[i+k] = store; 
            }

            // See if we have any values that cross zero
            for(i = 0; i < nstore-1; ++i)
            {
                if(Growth[i]*Growth[i+1] < 0.0)
                {
                    break;
                }
            }
            
            if(i != nstore-1)
            {
                if(nstore == 2)
                {
                    alp_new = (Alpha[0]*Growth[1] - Alpha[1]*Growth[0])/(Growth[1]-Growth[0]);
                }
                else
                {
                    // use a quadratic fit and step through 10000 points
                    // to find zero.
                    int     j; 
                    int     nsteps = 10000;
                    int     idx = (i == 0)?1:i;
                    double  da = Alpha[idx+1] - Alpha[idx-1];
                    double  gval_m1 = Growth[idx-1],a,gval;
                    double  c1 = Growth[idx-1]/(Alpha[idx-1]-Alpha[idx])/
                        (Alpha[idx-1]-Alpha[idx+1]);
                    double  c2 = Growth[idx]/(Alpha[idx]-Alpha[idx-1])/
                        (Alpha[idx]-Alpha[idx+1]);
                    double  c3 = Growth[idx+1]/(Alpha[idx+1]-Alpha[idx-1])/
                        (Alpha[idx+1]-Alpha[idx]);
                    
                    for(j = 1; j < nsteps+1; ++j)
                    {
                        a = Alpha[i] + j*da/(double) nsteps;
                        gval = c1*(a - Alpha[idx  ])*(a - Alpha[idx+1]) 
                            +  c2*(a - Alpha[idx-1])*(a - Alpha[idx+1])
                            +  c3*(a - Alpha[idx-1])*(a - Alpha[idx]);
                        
                        if(gval*gval_m1 < 0.0)
                        {
                            alp_new = ((a+da/(double)nsteps)*gval - a*gval_m1)/
                                (gval - gval_m1);
                            break;
                        }
                    }
                }
            }
            else // step backward/forward 
            {
                if(Growth[i] > 0.0)
                {
                    alp_new = m_alpha[0] + m_alphaStep;
                }
                else
                {
                    alp_new = m_alpha[0] - m_alphaStep;
                }
            }
        }
        
        for(int i = m_alpha.num_elements()-1; i > 0; --i)
        {
            m_alpha[i] = m_alpha[i-1];
            m_leading_real_evl[i] = m_leading_real_evl[i-1];
            m_leading_imag_evl[i] = m_leading_imag_evl[i-1];
        }

        m_alpha[0] = alp_new;    
    }
    
    Array<OneD, int> VortexWaveInteraction::GetReflectionIndex(void)
    {
        int i,j;
        int npts = m_waveVelocities[0]->GetPlane(0)->GetNpoints();
        int nel  = m_waveVelocities[0]->GetNumElmts();
        Array<OneD, int> index(npts);

        Array<OneD, NekDouble> coord(2);
        Array<OneD, NekDouble> coord_x(npts);
        Array<OneD, NekDouble> coord_y(npts);
        
        //-> Dermine the point which is on coordinate (x -> -x + Lx/2, y-> -y)
        m_waveVelocities[0]->GetPlane(0)->GetCoords(coord_x,coord_y);
        NekDouble xmax = Vmath::Vmax(npts,coord_x,1);
        //NekDouble tol = NekConstants::kGeomFactorsTol*NekConstants::kGeomFactorsTol;
        NekDouble tol = 1e-5;
        NekDouble xnew,ynew;

        int start  = npts-1; 
        int e_npts;

        bool useOnlyQuads = false;
        if(m_sessionVWI->DefinesSolverInfo("SymmetriseOnlyQuads"))
        {
            useOnlyQuads = true;
        }
        
        int cnt;
        for(int e = 0; e < nel; ++e)
        {
            e_npts = m_waveVelocities[0]->GetExp(e)->GetTotPoints();
            cnt = m_waveVelocities[0]->GetPhys_Offset(e);
            
            if(useOnlyQuads)
            {
                if(m_waveVelocities[0]->GetExp(e)->DetExpansionType() == StdRegions::eTriangle)
                {
                    for(i = 0; i < e_npts; ++i)
                    {
                        index[cnt+i] = -1;
                    }
                    continue;
                }
            }
            
            for(i = cnt; i < cnt+e_npts; ++i)
            {
                xnew = - coord_x[i]  + xmax;
                ynew = - coord_y[i];
                
                for(j = start; j >=0 ; --j)
                {
                    if((coord_x[j]-xnew)*(coord_x[j]-xnew) + (coord_y[j]-ynew)*(coord_y[j]-ynew) < tol)
                    {
                        index[i] = j;
                        start = j;
                        break;
                    }
                }
                
                if(j == -1)
                {
                    
                    for(j = npts-1; j > start; --j)
                    {
                        
                        if((coord_x[j]-xnew)*(coord_x[j]-xnew) + (coord_y[j]-ynew)*(coord_y[j]-ynew) < tol)
                        {
                            index[i] = j;
                            break;
                        }
                    }
                    ASSERTL0(j != start,"Failed to find matching point");
                }
            }
        }
        return index;
    }


    void VortexWaveInteraction::FileRelaxation(int reg)
    {
          cout<<"relaxation..."<<endl;
          static int cnt=0;
          Array<OneD, MultiRegions::ExpListSharedPtr> Iexp 
                                           =m_rollField[0]->GetBndCondExpansions();
          //cast to 1D explist (otherwise appenddata doesn't work)
          MultiRegions::ExpList1DSharedPtr Ilayer;  
          Ilayer = MemoryManager<MultiRegions::ExpList1D>::
                          AllocateSharedPtr(  
                          *boost::static_pointer_cast<MultiRegions::ExpList1D>(Iexp[reg]));
          int nq = Ilayer->GetTotPoints();
          if( cnt==0)
          {
                m_bcsForcing = Array<OneD, Array<OneD, NekDouble> > (4);
                m_bcsForcing[0] = Array<OneD, NekDouble> (4*nq);
                for(int i = 1; i < 4; ++i)
                {
                      m_bcsForcing[i] = m_bcsForcing[i-1] + nq;
                }           
          }

          // Read in mesh from input file

          /// ====================================================
          /// \todo Please update to use MeshGraph::Read(vSession)
          /// ====================================================
          SpatialDomains::MeshGraphSharedPtr graphShPt = 
                                     SpatialDomains::MeshGraph::Read(m_sessionName+".xml");


          std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef_u;
          std::vector<std::vector<NekDouble> > FieldData_u;
          string file = m_sessionName;


          file += "_u_5.bc"; 
          graphShPt->Import(file,FieldDef_u, FieldData_u);
          Ilayer->ExtractDataToCoeffs(FieldDef_u[0], FieldData_u[0], FieldDef_u[0]->m_fields[0],Ilayer->UpdateCoeffs());
          Ilayer->BwdTrans_IterPerExp(Ilayer->GetCoeffs(), Ilayer->UpdatePhys());
          
          if(cnt==0)
          {
               Vmath::Vcopy(nq,Ilayer->UpdatePhys(),1,m_bcsForcing[2],1);
          }
          Vmath::Vcopy(nq,Ilayer->UpdatePhys(),1,m_bcsForcing[0],1);

          //case relaxation==0 an additional array is needed
          Array<OneD, NekDouble> tmp_forcing(nq, 0.0);

          if(cnt!=0)
          {
              if(m_vwiRelaxation==1.0)
              {
                 Vmath::Vcopy(nq, m_bcsForcing[0],1, tmp_forcing,1);
              }
              Vmath::Smul(nq,1.0-m_vwiRelaxation,
                        m_bcsForcing[0],1,m_bcsForcing[0],1);
              Vmath::Svtvp(nq,m_vwiRelaxation,m_bcsForcing[2],1,
                         m_bcsForcing[0],1,Ilayer->UpdatePhys(),1);
              //generate again the bcs files:
              
    	      Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(1);   
              Ilayer->FwdTrans_IterPerExp(Ilayer->GetPhys(),Ilayer->UpdateCoeffs()); 
              fieldcoeffs[0] = Ilayer->UpdateCoeffs();		
	      std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef1  = Ilayer->GetFieldDefinitions();               
              std::vector<std::vector<NekDouble> > FieldData_1(FieldDef1.size());;
              FieldDef1[0]->m_fields.push_back("u");            	    
              Ilayer->AppendFieldData(FieldDef1[0], FieldData_1[0]);            	    
              graphShPt->Write(file,FieldDef1,FieldData_1); 
              //save the bcs for the next iteration
              if(m_vwiRelaxation!=1.0)
              {
                   Vmath::Smul(nq,1./(1.0-m_vwiRelaxation),
                        m_bcsForcing[0],1,m_bcsForcing[0],1);              
                   Vmath::Vcopy(nq,m_bcsForcing[0],1,m_bcsForcing[2],1);
              }
              else
              {
                   Vmath::Vcopy(nq, tmp_forcing,1, m_bcsForcing[2],1);                   
              }
          }
                   


          file = m_sessionName+ "_v_5.bc"; 

          std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef_v;
          std::vector<std::vector<NekDouble> > FieldData_v;
          graphShPt->Import(file,FieldDef_v, FieldData_v);
          Ilayer->ExtractDataToCoeffs(FieldDef_v[0], FieldData_v[0], FieldDef_v[0]->m_fields[0],Ilayer->UpdateCoeffs());
          Ilayer->BwdTrans_IterPerExp(Ilayer->GetCoeffs(), Ilayer->UpdatePhys());
          if(cnt==0)
          {
               Vmath::Vcopy(nq,Ilayer->UpdatePhys(),1,m_bcsForcing[3],1);
          }
          Vmath::Vcopy(nq,Ilayer->UpdatePhys(),1,m_bcsForcing[1],1);
          if(cnt!=0)
          {
              if(m_vwiRelaxation==1.0)
              {
                 Vmath::Vcopy(nq, m_bcsForcing[1],1, tmp_forcing,1);
              }
              Vmath::Smul(nq,1.0-m_vwiRelaxation,
                        m_bcsForcing[1],1,m_bcsForcing[1],1);
              Vmath::Svtvp(nq,m_vwiRelaxation,m_bcsForcing[3],1,
                         m_bcsForcing[1],1,Ilayer->UpdatePhys(),1);
              //generate again the bcs files:
    	      Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(1);   
              Ilayer->FwdTrans_IterPerExp(Ilayer->GetPhys(),Ilayer->UpdateCoeffs()); 
              fieldcoeffs[0] = Ilayer->UpdateCoeffs();		
	      std::vector<SpatialDomains::FieldDefinitionsSharedPtr>  FieldDef2  = Ilayer->GetFieldDefinitions();         
              std::vector<std::vector<NekDouble> > FieldData_2(FieldDef2.size());;      
              FieldDef2[0]->m_fields.push_back("v");            	    
              Ilayer->AppendFieldData(FieldDef2[0], FieldData_2[0]);            	             	
              graphShPt->Write(file,FieldDef2,FieldData_2); 
              //save the bcs for the next iteration
              if(m_vwiRelaxation!=1.0)
              {
                   Vmath::Smul(nq,1./(1.0-m_vwiRelaxation),
                        m_bcsForcing[1],1,m_bcsForcing[1],1);              
                   Vmath::Vcopy(nq,m_bcsForcing[1],1,m_bcsForcing[3],1);
              }
              else
              {
                   Vmath::Vcopy(nq, tmp_forcing,1, m_bcsForcing[3],1);                   
              }


          }


           cnt++;
    }
}
