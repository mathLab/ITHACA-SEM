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

#include<VortexWaveInteraction/VortexWaveInteraction.h>

#if defined(_MSC_VER) && BOOST_VERSION > 104700
#include <windows.h>
#undef MoveFile
#endif

void DoFixedForcingIteration(VortexWaveInteraction &vwi);
void Mvdir(string dir, NekDouble dir_ending);

int main(int argc, char *argv[])
{
    
    if(argc != 2)
    {
        cerr << "\n \t Usage: VortexWaveInteractionSolver  input \n" << endl;
        exit(1);
    }

    try
    { 
        VortexWaveInteraction vwi(argc,argv);
        
        for(int i = 0; i < vwi.GetMaxWaveForceMagIter(); ++i)
        {
            DoFixedForcingIteration(vwi);
            
            NekDouble WaveForce = vwi.GetWaveForceMag();
            
            vwi.AppendEvlToFile("ConvergedSolns",WaveForce);            

            vwi.UpdateWaveForceMag(WaveForce + vwi.GetWaveForceMagStep());
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
    system(syscall.c_str());
}

void DoFixedForcingIteration(VortexWaveInteraction &vwi)
{

    switch(vwi.GetVWIIterationType())
    {
    case eFixedAlphaWaveForcing:
        
        for(int i = vwi.GetIterStart(); i < vwi.GetIterEnd(); ++i)
        {
            vwi.ExecuteLoop();
            #ifdef _MSC_VER
            Sleep(10000);
            #else
            sleep(10);
            #endif
            vwi.SaveLoopDetails("Save",i);
            vwi.AppendEvlToFile("conv.his",i);            
        }
        break;
    case eFixedWaveForcing:
        {
            int i;
            int nouter_iter = vwi.GetNOuterIterations();
            bool exit_iteration = false;
            
            while(exit_iteration == false)
            {
                for(i = vwi.GetIterStart(); i < vwi.GetIterEnd(); ++i)
                {
                    vwi.ExecuteLoop();
                    #ifdef _MSC_VER
                    Sleep(10000);
                    #else
                    sleep(10);
                    #endif
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
                    vwi.UpdateAlpha(nouter_iter);
                }
                
                if(nouter_iter >= vwi.GetMaxOuterIterations())
                {
                    cerr << "Failed to converge after "<< vwi.GetMaxOuterIterations() << " outer iterations" << endl;
                    exit_iteration = true;
                }
            }
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
                #ifdef _MSC_VER
                Sleep(10000);
                #else
                sleep(10);
                #endif

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
                        #ifdef _MSC_VER
                        Sleep(10000);
                        #else
                        sleep(10);
                        #endif
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
                if (i == 0)
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


