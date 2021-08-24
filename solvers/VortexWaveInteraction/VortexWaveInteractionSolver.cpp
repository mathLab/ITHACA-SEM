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

#include<VortexWaveInteraction/VortexWaveInteraction.h>

using namespace std;
using namespace Nektar;

#if defined(_MSC_VER)
#include <windows.h>
#undef MoveFile
#endif

void DoFixedForcingIteration(VortexWaveInteraction &vwi);
void Mvdir(string dir, NekDouble dir_ending);

int main(int argc, char *argv[])
{
    try
    { 
        VortexWaveInteraction vwi(argc,argv);

        for(int i = 0; i < vwi.GetMaxWaveForceMagIter(); ++i)
        {
            DoFixedForcingIteration(vwi);
            
            NekDouble WaveForce = vwi.GetWaveForceMag();
            
            vwi.AppendEvlToFile("ConvergedSolns",WaveForce);            

            vwi.SetWaveForceMag(WaveForce + vwi.GetWaveForceMagStep());
            vwi.SetPrevAlpha(vwi.GetAlpha());
            //vwi.SetAlpha(vwi.GetAlpha() + vwi.GetDAlphaDWaveForceMag()*vwi.GetWaveForceMagStep());
            vwi.SetAlpha(vwi.GetAlpha() + (vwi.GetWaveForceMagStep() > 0?vwi.GetAlphaStep():(-vwi.GetAlphaStep())));
            
            // Save data directories. 
            if(vwi.GetVWIIterationType() == eFixedAlphaWaveForcing)
            {
                Mvdir("Save",WaveForce);
            }
            else
            {
                Mvdir("Save_Outer",WaveForce);
                // Execute Another loop so that not same initial conditions as last iteration
                vwi.ExecuteLoop();
            }
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
}

void Mvdir(string dir, NekDouble dir_ending)
{
    // save OuterIter.his if exists
    string saveOuterIter = "mv -f OuterIter.his "+ dir;
    if(system(saveOuterIter.c_str()))
    {
        ASSERTL0(false,saveOuterIter.c_str());
    }
    
    // Mv directory
    string newdir  = dir + boost::lexical_cast<std::string>(dir_ending);
    string syscall = "mv -f " + dir + " " + newdir;

    if(system(syscall.c_str()))
    {
        ASSERTL0(false,syscall.c_str());
    }
    
    // make new directory
    syscall = "mkdir " + dir;
    ASSERTL0(system(syscall.c_str()) == 0,
             "Failed to make directory '" + dir + "'");
}

void DoFixedForcingIteration(VortexWaveInteraction &vwi)
{

    // Reset eigenvalue checker in case used in previous iterations 
    vwi.CheckEigIsStationary(true);

    switch(vwi.GetVWIIterationType())
    {
    case eFixedAlpha:
        {
            int i;
            int nouter_iter = vwi.GetNOuterIterations();
            bool exit_iteration = false;
            
            while(exit_iteration == false)
            {
                // Reset eigenvalue checker in case used in previous iterations 
                vwi.CheckEigIsStationary(true);

                for(i = vwi.GetIterStart(); i < vwi.GetIterEnd(); ++i)
                {
                    vwi.ExecuteLoop();
                    vwi.SaveLoopDetails("Save", i);
                    vwi.AppendEvlToFile("conv.his",i);            
                        
                    if(vwi.CheckEigIsStationary())
                    {
                        vwi.SaveLoopDetails("Save_Outer", nouter_iter);
                        break;
                    }
                }
                
                // check to see if growth was converged. 
                if(i == vwi.GetIterEnd())
                {
                    cout << "Failed to converge growth rate in" << 
                        " inner iteration after " << vwi.GetIterEnd() 
                         << " loops" << endl;
                    exit(1);
                }
                
                vwi.AppendEvlToFile("OuterIter.his",nouter_iter++);            
                exit_iteration = vwi.CheckIfAtNeutralPoint();
                if(exit_iteration == false)
                {
                    vwi.UpdateWaveForceMag(nouter_iter);
                }
                
                if(nouter_iter >= vwi.GetMaxOuterIterations())
                {
                    cerr << "Failed to converge after "<< vwi.GetMaxOuterIterations() << " outer iterations" << endl;
                    exit_iteration = true;
                }
            }
        }
        break;
    case eFixedWaveForcing:
        {
            int i;
            int nouter_iter = vwi.GetNOuterIterations();
            bool exit_iteration = false;
	    NekDouble saveEigRelTol = vwi.GetEigRelTol();
	    int saveMinIters = vwi.GetMinInnerIterations();
	    int init_search = 0;
            
            
            if(init_search)
            {
                // initial set m_eigelTol to 1e-1 and inner iterations to 1 for 
                // quick search 
                vwi.SetEigRelTol(1e-1);
                vwi.SetMinInnerIterations(2);
            }

            while(exit_iteration == false)
            {
                // Reset eigenvalue checker in case used in previous iterations 
                vwi.CheckEigIsStationary(true);

                for(i = vwi.GetIterStart(); i < vwi.GetIterEnd(); ++i)
                {
                    vwi.ExecuteLoop();
                    vwi.SaveLoopDetails("Save", i);
                    vwi.AppendEvlToFile("conv.his",i);                    

		    if(vwi.CheckEigIsStationary())
                    {
                        vwi.SaveLoopDetails("Save_Outer", nouter_iter);
                        break;
                    }
                }
                
                // check to see if growth was converged. 
                if(i == vwi.GetIterEnd() )
                {
                    cout << "Failed to converge growth rate in" << 
                        " inner iteration after " << vwi.GetIterEnd() 
                         << " loops" << endl;
                    exit(1);
                }
                
                vwi.AppendEvlToFile("OuterIter.his",nouter_iter++);            
                
		exit_iteration = vwi.CheckIfAtNeutralPoint();

                cout << "m_alpha[0] is " << vwi.GetAlpha() << endl;

                if(exit_iteration == false)
                {
                    vwi.UpdateAlpha(nouter_iter);
                    if(vwi.IfIterInterface()==true)
                    {
                          vwi.CalcNonLinearWaveForce();
                    }
                }

		// Redo iteration if at first coarse search 
		if((exit_iteration == true) && (init_search == 1))
		{
		  init_search = 0;
		  vwi.SetEigRelTol(saveEigRelTol);
		  vwi.SetMinInnerIterations(saveMinIters);
		  nouter_iter = 1;
		  exit_iteration = false;
		}
                
                if(nouter_iter >= vwi.GetMaxOuterIterations())
                {
                    cerr << "Failed to converge after "<< vwi.GetMaxOuterIterations() << " outer iterations" << endl;
		    exit_iteration = true;
                }
            }
            
            vwi.UpdateDAlphaDWaveForceMag(vwi.GetPrevAlpha());

        }
        break;
    case eFixedAlphaWaveForcing:
        
        for(int i = vwi.GetIterStart(); i < vwi.GetIterEnd(); ++i)
        {
            vwi.ExecuteLoop();
            vwi.SaveLoopDetails("Save",i);
            vwi.AppendEvlToFile("conv.his",i);            
        }
        break;
    case eFixedWaveForcingWithSubIterationOnAlpha:
        {
            int i;
            int nouter_iter = vwi.GetNOuterIterations();
            bool exit_iteration = false;
            
            while(exit_iteration == false)
            {
                bool exit_alphaIter = false;
                
                vwi.ExecuteLoop(false);

                // Sub iterate Alpha
                for(i = 0; i < vwi.GetIterEnd(); ++i)
                {
                    vwi.SaveLoopDetails("Save", i);
                    
                    vwi.AppendEvlToFile("AlphaIter.his",i);            
                    
                    exit_alphaIter = vwi.CheckIfAtNeutralPoint();
                    if(exit_alphaIter == false)
                    {
                        vwi.UpdateAlpha(i+1);
                        vwi.ExecuteWave();
                    }
                    else
                    {
                        vwi.CalcNonLinearWaveForce();
                        break;
                    }
                }
                
                // check to see if growth was converged. 
                if(i == vwi.GetIterEnd())
                {
                    cout << "Failed to converge growth rate in" << 
                        " inner iteration after " << vwi.GetIterEnd() 
                         << " loops" << endl;
                    exit(1);
                }
                
                vwi.MoveFile("AlphaIter.his","Save_Outer", nouter_iter);
                vwi.SaveLoopDetails("Save_Outer", nouter_iter);
                vwi.AppendEvlToFile("OuterIter.his",nouter_iter++);            
                
                // assume that if only previous inner loop has
                // only done one iteration then we are at neutral
                // point
                if (i == 0  && vwi.IfIterInterface()==false)
                {
                    exit_iteration = true;
                }
                
                
                if(nouter_iter >= vwi.GetMaxOuterIterations())
                {
                    cerr << "Failed to converge after "<< vwi.GetMaxOuterIterations() << " outer iterations" << endl;
                    exit_iteration = true;
                }
            }
        }
        break;
    default:
        ASSERTL0(false,"Unknown iteration type");
    }
    
}


