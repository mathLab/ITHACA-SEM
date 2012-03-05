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
        m_sessionVWI = LibUtilities::SessionReader::CreateInstance(argc, argv, VWIFilenames);
        

        m_sessionVWI->LoadParameter("AlphaStep",m_alphaStep,0.05);
        m_sessionVWI->LoadParameter("OuterIterationStoreSize",storesize,10);
        m_sessionVWI->LoadParameter("EigenvalueRelativeTol",m_eigRelTol,1e-3);
        m_sessionVWI->LoadParameter("NeutralPointTolerance",m_neutralPointTol,1e-4);
        m_sessionVWI->LoadParameter("MaxOuterIterations",m_maxOuterIterations,100);
        
        m_sessionVWI->LoadParameter("StartIteration", m_iterStart, 0);
        m_sessionVWI->LoadParameter("EndIteration", m_iterEnd, 0);

        m_sessionVWI->LoadParameter("WaveForceMagStep",m_waveForceMagStep,0.01);
        m_sessionVWI->LoadParameter("MaxWaveForceMagIter",m_maxWaveForceMagIter,1);
        m_sessionVWI->LoadParameter("RollForceScale",m_rollForceScale,1.0);
        
        if(m_sessionVWI->DefinesSolverInfo("DeltaFcnApprox"))
        {
            m_deltaFcnApprox = true;
            m_sessionVWI->LoadParameter("DeltaFcnDecay",m_deltaFcnDecay,1.0/500);
        }
        else
        {
            m_deltaFcnApprox = false;
            m_deltaFcnDecay = 0.0;
        }

        m_alpha = Array<OneD, NekDouble> (storesize);
        m_alpha[0]        = m_sessionVWI->GetParameter("Alpha");
        m_waveForceMag      = Array<OneD, NekDouble> (storesize);
        m_waveForceMag[0]  = m_sessionVWI->GetParameter("WaveForceMag");
        
        m_leading_real_evl = Array<OneD, NekDouble> (storesize);
        m_leading_imag_evl = Array<OneD, NekDouble> (storesize);
        
        if(m_sessionVWI->DefinesParameter("Relaxation"))
        {
            m_vwiRelaxation = m_sessionVWI->GetParameter("Relaxation");
        }
        else
        {
            m_vwiRelaxation = 0.0;
        }
        
        // Initialise NS Roll solver 
        std::string IncCondFile(argv[argc-1]); 
        IncCondFile += "_IncNSCond.xml"; 
        std::vector<std::string> IncNSFilenames;
        IncNSFilenames.push_back(meshfile);   
        IncNSFilenames.push_back(IncCondFile);

        // Create Incompressible NavierStokesSolver session reader.
        m_sessionRoll = LibUtilities::SessionReader::CreateInstance(argc, argv, IncNSFilenames, m_sessionVWI->GetComm());
        std::string vEquation = m_sessionRoll->GetSolverInfo("SolverType");
        m_solverRoll = GetEquationSystemFactory().CreateInstance(vEquation, m_sessionRoll);


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
                    if(fp = fopen("OuterIter.his","r"))
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
            cout<<"eqtype="<<vEquation<<endl;
            EquationSystemSharedPtr solverRoll = GetEquationSystemFactory().CreateInstance(vEquation,m_sessionRoll);
            //the forcing terms are inserted as N bcs
              //the forcing terms are inserted as N bcs
             // Execute Roll 
             cout << "Executing Roll solver" << endl;
             solverRoll->DoInitialise();
             solverRoll->DoSolve();
             solverRoll->Output();
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
            // Read vwi file
            std::string forcefile
                = m_sessionRoll->GetFunctionFilename("BodyForce");
            
            if(forcefile != "")
            {
                m_solverRoll->ImportFld(forcefile,m_solverRoll->UpdateForces());
                
                // Scale forcing
                int npoints = m_solverRoll->UpdateForces()[0]->GetPhys().num_elements();
                for(int i = 0; i < m_solverRoll->UpdateForces().num_elements(); ++i)
                {
                    Vmath::Smul(npoints,m_rollForceScale,m_solverRoll->UpdateForces()[i]->UpdatePhys(),1,m_solverRoll->UpdateForces()[i]->UpdatePhys(),1);
                }
            }
        }

        // Execute Roll 
        cout << "Executing Roll solver" << endl;
        m_solverRoll->DoInitialise();
        m_solverRoll->DoSolve();
        m_solverRoll->Output();
        
        // Copy .fld file to .rst and base.fld
        cout << "Executing cp -f session.fld session.rst" << endl;
        CopyFile(".fld",".rst");
        cout << "Executing cp -f session.fld session-Base.fld" << endl;
        CopyFile(".fld","-Base.fld");


    }


    void VortexWaveInteraction::ExecuteStreak(void)
    {


        // Create driver
#if 0
        std::string vDriverModule;
        m_sessionStreak->LoadSolverInfo("Driver", vDriverModule, "Standard");
        
        DriverSharedPtr solverStreak = GetDriverFactory().CreateInstance(vDriverModule, m_sessionStreak); 
        solverStreak->Execute();
#else        
        // Setup and execute Advection Diffusion solver 
        string vEquation = m_sessionStreak->GetSolverInfo("EqType");
        EquationSystemSharedPtr solverStreak = GetEquationSystemFactory().CreateInstance(vEquation,m_sessionStreak);

        cout << "Executing Streak Solver" << endl;
        solverStreak->DoInitialise();
        solverStreak->DoSolve();
        solverStreak->Output();
#endif

        m_streakField = solverStreak->UpdateFields();
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
        
    }

    void VortexWaveInteraction::CalcNonLinearWaveForce(void)
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
        m_wavePressure->GetPlane(0)->BwdTrans(m_wavePressure->GetPlane(0)->GetCoeffs(),
                                              m_wavePressure->GetPlane(0)->UpdatePhys());
        m_wavePressure->GetPlane(1)->BwdTrans(m_wavePressure->GetPlane(1)->GetCoeffs(),
                                              m_wavePressure->GetPlane(1)->UpdatePhys());

        // Determine normalisation of pressure so that |P|/A = 1
        NekDouble norm = 0, l2;
        l2    = m_wavePressure->GetPlane(0)->L2();
        norm  = l2*l2;
        l2    = m_wavePressure->GetPlane(1)->L2();
        norm += l2*l2;
        Vmath::Fill(2*npts,1.0,der1,1);
        NekDouble area = m_waveVelocities[0]->GetPlane(0)->PhysIntegral(der1);
        norm = sqrt(area/norm);
        
        // Get hold of arrays. 
        m_waveVelocities[0]->GetPlane(0)->BwdTrans(m_waveVelocities[0]->GetPlane(0)->GetCoeffs(),m_waveVelocities[0]->GetPlane(0)->UpdatePhys());
        Array<OneD, NekDouble> u_real = m_waveVelocities[0]->GetPlane(0)->UpdatePhys();
        Vmath::Smul(npts,norm,u_real,1,u_real,1);
        m_waveVelocities[0]->GetPlane(1)->BwdTrans(m_waveVelocities[0]->GetPlane(1)->GetCoeffs(),m_waveVelocities[0]->GetPlane(1)->UpdatePhys());
        Array<OneD, NekDouble> u_imag = m_waveVelocities[0]->GetPlane(1)->UpdatePhys();
        Vmath::Smul(npts,norm,u_imag,1,u_imag,1);
        m_waveVelocities[1]->GetPlane(0)->BwdTrans(m_waveVelocities[1]->GetPlane(0)->GetCoeffs(),m_waveVelocities[1]->GetPlane(0)->UpdatePhys());
        Array<OneD, NekDouble> v_real = m_waveVelocities[1]->GetPlane(0)->UpdatePhys(); 
        Vmath::Smul(npts,norm,v_real,1,v_real,1);
        m_waveVelocities[1]->GetPlane(1)->BwdTrans(m_waveVelocities[1]->GetPlane(1)->GetCoeffs(),m_waveVelocities[1]->GetPlane(1)->UpdatePhys());
        Array<OneD, NekDouble> v_imag = m_waveVelocities[1]->GetPlane(1)->UpdatePhys();
        Vmath::Smul(npts,norm,v_imag,1,v_imag,1);
        
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
        
#if 0 
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
            der2 [i] = m_waveVelocities[0]->GetPlane(0)->GetExp(Eid[i])->PhysEvaluate(coord,
                                                                            der1 + physoffset);
        }
        //-> Average field 0 
        Vmath::Vsub(npts,der1,1,der2,1,der2,1);
        Vmath::Smul(npts,0.5,der2,1,der2,1);
        m_waveVelocities[0]->GetPlane(0)->FwdTrans_BndConstrained(der2, m_vwiForcing[0]);
        
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
        m_waveVelocities[0]->GetPlane(0)->FwdTrans_BndConstrained(der2, m_vwiForcing[1]);
#else
        int i;
        static Array<OneD, int> index = GetReflectionIndex();

        m_waveVelocities[0]->GetPlane(0)->BwdTrans_IterPerExp(m_vwiForcing[0],der1);
        for(i = 0; i < npts; ++i)
        {
            val[i] = 0.5*(der1[i] - der1[index[i]]);
        }

        m_waveVelocities[0]->GetPlane(0)->FwdTrans_BndConstrained(val, m_vwiForcing[0]);


        m_waveVelocities[0]->GetPlane(0)->BwdTrans_IterPerExp(m_vwiForcing[1],der1);
        for(i = 0; i < npts; ++i)
        {
            val[i] = 0.5*(der1[i] - der1[index[i]]);
        }        
            
        m_waveVelocities[0]->GetPlane(0)->FwdTrans_BndConstrained(val, m_vwiForcing[1]);
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
        Array<OneD, std::string> variables(2);
        Array<OneD, Array<OneD, NekDouble> > outfield(2);
        variables[0] = "u";   variables[1] = "v";
        outfield[0]  = m_vwiForcing[0];
        outfield[1]  = m_vwiForcing[1];
        
        std::string outname = m_sessionName  + ".vwi";
        
        m_solverRoll->WriteFld(outname, m_waveVelocities[0]->GetPlane(0), outfield, variables);
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
        ExecuteRoll();

#ifndef _WIN32
	sleep(3);
#endif
        ExecuteStreak();
#ifndef _WIN32
	sleep(3);
#endif

	//the global info has to be written in the 
	//roll session file to use the interface loop
        if(m_sessionRoll->DefinesSolverInfo("INTERFACE"))
	{
             static int cnt=0;       
             string syscall;
             //rewrite the Rollsessionfile (we start from the waleffe forcing)
             //string meshbndjumps = m_sessionName +"_bndjumps.xml";             
             //if(cnt==0)
             //{
                 //take the conditions tag from meshbndjumps and copy into 
                 // the rolls session file
             //}
             char c[16]="";
    	     sprintf(c,"%d",cnt);  
             //save old roll solution
             string oldroll = m_sessionName +"_roll_"+c +".fld";    
	     syscall = "cp -f " + m_sessionName+"-Base.fld" + "  " + oldroll;
             cout<<syscall.c_str()<<endl;
             if(system(syscall.c_str()))
             {
                  ASSERTL0(false,syscall.c_str());
             } 
	     //move the mesh around the critical layer
             string filePost   = m_sessionName + "_advPost.xml";
             string filestreak   = m_sessionName + "_streak.fld";
             string fileinterp = m_sessionName + "_interp.xml";
             syscall  = "../../../utilities/builds/PostProcessing/Extras/MoveMesh-g  "
                             + filePost +"  "+ filestreak +"  "+ fileinterp; 
/*
             syscall  = "../../../utilities/builds/PostProcessing/Extras/MoveMesh-g  "
                             + filePost +"  "+ filestreak +"  "+ filePost; 
*/
             cout<<syscall.c_str()<<endl;
             if(system(syscall.c_str()))
             {
                  ASSERTL0(false,syscall.c_str());
             }

             //move the advPost mesh
             syscall  =  "../../../utilities/builds/PostProcessing/Extras/MoveMesh-g  "
                      + filePost + "  " + filestreak + "  " + filePost;
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
             string interpstreak = m_sessionName +"_interpstreak_"+ c +".fld";  
             syscall  =  "../../../utilities/builds/PostProcessing/Extras/FieldToField-g  "
                      + fileinterp + "  " + filestreak + "  " + movedinterpmesh + "  " 
	              + interpstreak;
/*
             syscall  =  "../../../utilities/builds/PostProcessing/Extras/FieldToField-g  "
                      + filePost + "  " + filestreak + "  " + movedmesh + "  " 
	              + interpstreak;
*/
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

             //overwriting the interp_moved file
	     syscall = "cp -f " + movedinterpmesh + "  " + fileinterp;
             cout<<syscall.c_str()<<endl;
             if(system(syscall.c_str()))
             {
                  ASSERTL0(false,syscall.c_str());
             }

             //calculate the wave
             ExecuteWave();
/*
             string meshLin = m_sessionName + "_Lin.xml";
	     syscall = "../IncNavierStokesSolver/IncNavierStokesSolver-g " + meshLin;
             cout<<syscall.c_str()<<endl;
             if(system(syscall.c_str()))
             {
                ASSERTL0(false,syscall.c_str());
             }
*/
cout << "Growth =" <<m_leading_real_evl[0]<<endl; 
cout << "Phase =" <<m_leading_imag_evl[0]<<endl; 
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
             char c1[16]="";
    	     sprintf(c1,"%d",cnt);   
             //calculate the jump conditions
             string wavefile  = m_sessionName +".fld"; 
             syscall =  "../../../utilities/builds/PostProcessing/Extras/FldCalcBCs-g  "
                     + movedmesh + "  " + wavefile + "  " + interpstreak + ">  data"+c1;
             cout<<syscall.c_str()<<endl;
             if(system(syscall.c_str()))
             {
                  ASSERTL0(false,syscall.c_str());
             }

             //use relaxation
             int reg =3;
             FileRelaxation(reg);

             //interpolate the .rst into the new mesh
             string rstfile = m_sessionName + ".rst";
             syscall  =  "../../../utilities/builds/PostProcessing/Extras/FieldToField-g  "
                      + filePost + "  " + rstfile + "  " + movedmesh + "  " 
	              + rstfile;
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
	else
	{
            ExecuteWave();

            if(CalcWaveForce)
            {
                CalcNonLinearWaveForce();
            }
	}







    }

    bool VortexWaveInteraction::CheckEigIsStationary(void)
    {
        static NekDouble previous_real_evl = -1.0; 
        static NekDouble previous_imag_evl = -1.0; 
        
        if(previous_real_evl == -1.0)
        {
            previous_real_evl = m_leading_real_evl[0];
            previous_imag_evl = m_leading_imag_evl[0];
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

        if((m_leading_real_evl[0]*m_leading_real_evl[0] + 
            m_leading_imag_evl[0]*m_leading_imag_evl[0]) < 
           m_neutralPointTol*m_neutralPointTol)
        {
            returnval = true;
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
          SpatialDomains::MeshGraphSharedPtr graphShPt = 
                                     SpatialDomains::MeshGraph::Read(m_sessionName+".xml");
          std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef_u;
          std::vector<std::vector<NekDouble> > FieldData_u;
          string file = m_sessionName;


          file += "_u_5.bc"; 
          graphShPt->Import(file,FieldDef_u, FieldData_u);
          Ilayer->ExtractDataToCoeffs(FieldDef_u[0], FieldData_u[0], FieldDef_u[0]->m_fields[0]);
          Ilayer->BwdTrans_IterPerExp(Ilayer->GetCoeffs(), Ilayer->UpdatePhys());
          
          if(cnt==0)
          {
               Vmath::Vcopy(nq,Ilayer->UpdatePhys(),1,m_bcsForcing[2],1);
          }
          Vmath::Vcopy(nq,Ilayer->UpdatePhys(),1,m_bcsForcing[0],1);



          if(cnt!=0)
          {
cout<<"ucnt="<<cnt<<endl;
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
              Vmath::Smul(nq,1./(1.0-m_vwiRelaxation),
                        m_bcsForcing[0],1,m_bcsForcing[0],1);              
              Vmath::Vcopy(nq,m_bcsForcing[0],1,m_bcsForcing[2],1);
          }
                   


          file = m_sessionName+ "_v_5.bc"; 

          std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef_v;
          std::vector<std::vector<NekDouble> > FieldData_v;
          graphShPt->Import(file,FieldDef_v, FieldData_v);
          Ilayer->ExtractDataToCoeffs(FieldDef_v[0], FieldData_v[0], FieldDef_v[0]->m_fields[0]);
          Ilayer->BwdTrans_IterPerExp(Ilayer->GetCoeffs(), Ilayer->UpdatePhys());
          if(cnt==0)
          {
               Vmath::Vcopy(nq,Ilayer->UpdatePhys(),1,m_bcsForcing[3],1);
          }
          Vmath::Vcopy(nq,Ilayer->UpdatePhys(),1,m_bcsForcing[1],1);
          if(cnt!=0)
          {
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
              Vmath::Smul(nq,1./(1.0-m_vwiRelaxation),
                        m_bcsForcing[1],1,m_bcsForcing[1],1);              
              Vmath::Vcopy(nq,m_bcsForcing[1],1,m_bcsForcing[3],1);


          }


           cnt++;

    }
}
    
