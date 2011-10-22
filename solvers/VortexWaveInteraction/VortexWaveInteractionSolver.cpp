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
        
        for(int i = vwi.GetIterStart(); i < vwi.GetIterEnd(); ++i)
        {
            vwi.ExecuteRoll();
            vwi.SaveFile(".rst","Save",i);
            
            vwi.ExecuteStreak();
            vwi.SaveFile("_streak.fld","Save",i);
            
            vwi.ExecuteWave();
            vwi.SaveFile(".evl","Save",i);
            vwi.SaveFile("_eig_0","Save",i);

            cout << "Calculating Nonlinear Wave Forcing" << endl;
            vwi.CalcNonLinearWaveForce();
            vwi.SaveFile(".vwi","Save",i+1);

            vwi.AppendEvlToFile("conv.his",i);
        }
        
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

