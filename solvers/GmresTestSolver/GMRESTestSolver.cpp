///////////////////////////////////////////////////////////////////////////////
//
// File ADRSolver.cpp
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
// Description: Advection Diffusion Reaction framework solver
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/GlobalLinSysIterativeSolver.h>
// #include <MultiRegions/AssemblyMap/AssemblyMap.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr session;
    //LibUtilities::FieldIOSharedPtr       fld;
    //SpatialDomains::MeshGraphSharedPtr   graph;
    //MultiRegions::ContField2DSharedPtr   field;
    //LibUtilities::EquationSharedPtr      icond, ex_sol;
    //StdRegions::ConstFactorMap           factors;

    MultiRegions::GlobalLinSysKey pKey(StdRegions::MatrixType::eMass);
    std::weak_ptr<MultiRegions::ExpList> pExp;
    
    MultiRegions::AssemblyMapSharedPtr NoUsePtr;

    // /media/zy3416/DATADRIVE1/development/nektar++/nektar++/build/solvers/DiffusionSolver/DiffusionSolver-g

    try
    {
        // Create session reader.
        session = LibUtilities::SessionReader::CreateInstance(argc, argv);

        // Get some information about the session
        string       sessionName    = session->GetSessionName();
        string       outFile        = sessionName + ".fld";
        
        MultiRegions::GlobalLinSysIterativeSolver LinSolver(pKey,pExp,NoUsePtr);

        LinSolver.initializeLinSys(session);

        int ndim = LinSolver.getdimension();

        Array<OneD, NekDouble> soltn_cg(ndim);

        LinSolver.DoConjugateGradient(ndim, LinSolver.getrhs(), soltn_cg, NoUsePtr, 0);

        Array<OneD, NekDouble> soltn_gmres(ndim);

        LinSolver.DoConjugateGradient(ndim, LinSolver.getrhs(), soltn_gmres, NoUsePtr, 0);

        NekDouble l2error=0.0;
        NekDouble tmp;
        for(int i; i< ndim; ++i)
        {   
            tmp = soltn_cg[i]-soltn_gmres[i];
            l2error += tmp*tmp;
        }

        cout << "The L2 difference between CG&GMRES is"<<l2error<<endl;



        



        // Finalise session
        session->Finalise();
    }
    catch (const std::runtime_error& e)
    {
        return 1;
    }
    catch (const std::string& eStr)
    {
        cout << "Error: " << eStr << endl;
    }

    return 0;
}
