////////////////////////////////////////////////////////////////////////////////
//
//  File: XmlToTecplot.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Output Tecplot description of mesh.
//
////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;

    if(argc < 2)
    {
        fprintf(stderr,"Usage: XmlToTecplot meshfile\n");
        exit(1);
    }

    LibUtilities::SessionReader::RegisterCmdLineFlag(
        "multi-zone", "m", "Output multi-zone format (one element per zone).");

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-1]);
    SpatialDomains::MeshGraphSharedPtr graphShPt
        = SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Set up Expansion information
    SpatialDomains::ExpansionMap emap = graphShPt->GetExpansions();
    SpatialDomains::ExpansionMapIter it;

    for (it = emap.begin(); it != emap.end(); ++it)
    {
        for (int i = 0; i < it->second->m_basisKeyVector.size(); ++i)
        {
            LibUtilities::BasisKey  tmp1 = it->second->m_basisKeyVector[i];
            LibUtilities::PointsKey tmp2 = tmp1.GetPointsKey();
            it->second->m_basisKeyVector[i] = LibUtilities::BasisKey(
                tmp1.GetBasisType(), tmp1.GetNumModes(),
                LibUtilities::PointsKey(tmp1.GetNumModes(),
                                        LibUtilities::ePolyEvenlySpaced));
        }
    }
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(1);

    switch(expdim)
    {
        case 1:
        {
            MultiRegions::ExpList1DSharedPtr Exp1D;
            Exp1D = MemoryManager<MultiRegions::ExpList1D>::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] = Exp1D;
            break;
        }
        case 2:
        {
            if(vSession->DefinesSolverInfo("HOMOGENEOUS"))
            {
                std::string HomoStr = vSession->GetSolverInfo("HOMOGENEOUS");
                MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;

                ASSERTL0(
                    HomoStr == "HOMOGENEOUS1D" || HomoStr == "Homogeneous1D" ||
                    HomoStr == "1D"            || HomoStr == "Homo1D",
                    "Only 3DH1D supported for XML output currently.");

                int nplanes;
                vSession->LoadParameter("HomModesZ", nplanes);

                // choose points to be at evenly spaced points at nplanes + 1
                // points
                const LibUtilities::PointsKey Pkey(
                    nplanes + 1, LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(
                    LibUtilities::eFourier, nplanes, Pkey);
                NekDouble lz = vSession->GetParameter("LZ");

                Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>
                    ::AllocateSharedPtr(
                        vSession, Bkey, lz, false, false, graphShPt);
                Exp[0] = Exp3DH1;
            }
            else
            {
                MultiRegions::ExpList2DSharedPtr Exp2D;
                Exp2D = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(vSession,graphShPt);
                Exp[0] =  Exp2D;
            }
            break;
        }
        case 3:
        {
            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = MemoryManager<MultiRegions::ExpList3D>::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] =  Exp3D;
            break;
        }
        default:
            ASSERTL0(false,"Expansion dimension not recognised");
            break;
    }

    //-----------------------------------------------

    //----------------------------------------------
    // Write solution  depending on #define
    string   outfile(strtok(argv[argc-1],"."));
    outfile +=  ".dat";
    ofstream outstrm(outfile.c_str());

    Exp[0]->WriteTecplotHeader(outstrm);

    if (vSession->DefinesCmdLineArgument("multi-zone"))
    {
        int nExp = Exp[0]->GetExpSize();

        for (int i = 0; i < nExp; ++i)
        {
            Exp[0]->WriteTecplotZone        (outstrm, i);
            Exp[0]->WriteTecplotConnectivity(outstrm, i);
        }
    }
    else
    {
        Exp[0]->WriteTecplotZone        (outstrm);
        Exp[0]->WriteTecplotConnectivity(outstrm);
    }

    outstrm.close();
    //----------------------------------------------
    return 0;
}
