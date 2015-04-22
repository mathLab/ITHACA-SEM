///////////////////////////////////////////////////////////////////////////////
//
// File: CollectionTiming.cpp
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
// Description: Small demo to run timings of various operators.
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>
#include <iomanip>

#include <boost/timer/timer.hpp>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ExpList3D.h>
#include <Collections/Collection.h>
#include <SpatialDomains/MeshGraph.h>

using namespace Nektar;

using boost::timer::cpu_timer;
using boost::timer::cpu_times;
using boost::timer::nanosecond_type;
using boost::timer::format;

MultiRegions::ExpListSharedPtr SetupExpList(
    int                                  N,
    LibUtilities::SessionReaderSharedPtr session,
    SpatialDomains::MeshGraphSharedPtr   graph,
    Collections::ImplementationType      impType)
{
    graph->SetExpansionsToPolyOrder(N);

    MultiRegions::ExpListSharedPtr expList =
        MemoryManager<MultiRegions::ExpList3D>::AllocateSharedPtr(
            session, graph);

    expList->CreateCollections(impType);

    return expList;
}

void printOutput(int N, int Ntest, cpu_timer &timer, bool fmt)
{
    cpu_times times = timer.elapsed();
    nanosecond_type total = times.user + times.system;

    const double sec = 1000000000.0L;

    // Normalize timings
    double wall_sec  = times.wall / sec;
    double total_sec = total / sec;

    wall_sec /= Ntest;
    total_sec /= Ntest;

    if (fmt)
    {
        cout << setw(6)  << N-1
             << setw(18) << wall_sec
             << setw(18) << total_sec
             << endl;
    }
    else
    {
        cout << "P = " << N-1 << ": "
             << wall_sec << " (wall) "
             << total_sec << " (total), "
             << (100.0 * total_sec / wall_sec) << "% CPU"
             << endl;
    }
}

int main(int argc, char *argv[])
{
    LibUtilities::SessionReader::RegisterCmdLineFlag(
        "data", "d", "Print in data format");
    LibUtilities::SessionReaderSharedPtr session
        = LibUtilities::SessionReader::CreateInstance(argc, argv);

    bool fmt = session->DefinesCmdLineArgument("data");

    MultiRegions::ExpListSharedPtr expList;

    cpu_timer timer;

    int Ntest, maxOrder;
    session->LoadParameter("Ntest",    Ntest,    1000);
    session->LoadParameter("maxOrder", maxOrder, 10);
    maxOrder++;

    string sl = fmt ? "# " : "";

    // Read in mesh
    SpatialDomains::MeshGraphSharedPtr graph =
        SpatialDomains::MeshGraph::Read(session);

    // BwdTrans operator
    for (int imp = 1; imp < 5; ++imp)
    {
        Collections::ImplementationType impType =
            (Collections::ImplementationType)imp;
        cout << sl << "Using " << Collections::ImplementationTypeMap[imp]
             << " Collection Implementation:" << endl;

        cout << sl << "BwdTrans Op: Ntest = " << Ntest << endl;
        for (int N = 2; N < maxOrder; ++N)
        {
            expList = SetupExpList(N, session, graph, impType);
            Array<OneD, NekDouble> input (expList->GetNcoeffs());
            Array<OneD, NekDouble> output(expList->GetNpoints());

            timer.start();
            for (int i = 0; i < Ntest; ++i)
            {
                expList->BwdTrans(input, output);
            }
            timer.stop();

            printOutput(N, Ntest, timer, fmt);
        }

        cout << sl << "IProductWRTBase Op: Ntest = " << Ntest << endl;
        for (int N = 2; N < maxOrder; ++N)
        {
            expList = SetupExpList(N, session, graph, impType);
            Array<OneD, NekDouble> input (expList->GetNpoints());
            Array<OneD, NekDouble> output(expList->GetNcoeffs());

            timer.start();
            for (int i = 0; i < Ntest; ++i)
            {
                expList->IProductWRTBase(input, output);
            }
            timer.stop();

            printOutput(N, Ntest, timer, fmt);
        }

        cout << sl << "IProductWRTDerivBase Op: Ntest = " << Ntest << endl;
        for (int N = 2; N < maxOrder; ++N)
        {
            expList = SetupExpList(N, session, graph, impType);
            int nDim = expList->GetCoordim(0);
            Array<OneD, Array<OneD, NekDouble> > input (nDim);
            Array<OneD, NekDouble>               output(expList->GetNcoeffs());

            for (int i = 0; i < nDim; ++i)
            {
                input[i] = Array<OneD, NekDouble>(expList->GetNpoints());
            }

            timer.start();
            for (int i = 0; i < Ntest; ++i)
            {
                expList->IProductWRTDerivBase(input, output);
            }
            timer.stop();

            printOutput(N, Ntest, timer, fmt);
        }

        cout << sl << "PhysDeriv Op: Ntest = " << Ntest << endl;
        for (int N = 2; N < maxOrder; ++N)
        {
            expList = SetupExpList(N, session, graph, impType);
            Array<OneD, NekDouble> input  (expList->GetNpoints());
            Array<OneD, NekDouble> output0(expList->GetNpoints());
            Array<OneD, NekDouble> output1(expList->GetNpoints());
            Array<OneD, NekDouble> output2(expList->GetNpoints());

            timer.start();
            for (int i = 0; i < Ntest; ++i)
            {
                expList->PhysDeriv(input, output0, output1, output2);
            }
            timer.stop();

            printOutput(N, Ntest, timer, fmt);
        }
    }
}
