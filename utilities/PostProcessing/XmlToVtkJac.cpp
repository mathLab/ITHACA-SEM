//#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        cerr << "Usage: XmlToVtk  meshfile" << endl;
        exit(1);
    }    
    
    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv);

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-1]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = 
        SpatialDomains::MeshGraph::Read(vSession); //meshfile);
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
            Exp1D = MemoryManager<MultiRegions::ExpList1D>
                ::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] = Exp1D;
            break;
        }
        case 2:
        {
            MultiRegions::ExpList2DSharedPtr Exp2D;
            Exp2D = MemoryManager<MultiRegions::ExpList2D>
                ::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] =  Exp2D;
            break;
        }
        case 3:
        {
            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = MemoryManager<MultiRegions::ExpList3D>
                ::AllocateSharedPtr(vSession,graphShPt);
            Exp[0] =  Exp3D;
            break;
        }
        default:
        {
            ASSERTL0(false,"Expansion dimension not recognised");
            break;
        }
    }
    //----------------------------------------------
    
    Array<OneD, NekDouble> x0 (Exp[0]->GetNpoints());
    Array<OneD, NekDouble> x1 (Exp[0]->GetNpoints());
    Array<OneD, NekDouble> x2 (Exp[0]->GetNpoints());
    Exp[0]->GetCoords(x0, x1, x2);

    //----------------------------------------------
    // Write out VTK file.
    string   outname(strtok(argv[argc-1],"."));
    outname += ".vtu";
    ofstream outfile(outname.c_str());

    Array<OneD, NekDouble> tmp;
    Exp[0]->WriteVtkHeader(outfile);
    // For each field write header and footer, since there is no field data.
    for(int i = 0; i < Exp[0]->GetExpSize(); ++i)
    {
        Exp[0]->WriteVtkPieceHeader(outfile,i);
        StdRegions::StdExpansionSharedPtr e = Exp[0]->GetExp(i);
        unsigned int npts = e->GetTotPoints();
        NekDouble jac = Vmath::Vmin(e->GetMetricInfo()->GetJac().num_elements(), e->GetMetricInfo()->GetJac(), 1);
        Vmath::Fill(npts, jac, tmp = Exp[0]->UpdatePhys() + Exp[0]->GetPhys_Offset(i), 1);
        Exp[0]->WriteVtkPieceData(outfile,i, "Jac");
        Exp[0]->WriteVtkPieceFooter(outfile,i);
    }
    Exp[0]->WriteVtkFooter(outfile);
    //----------------------------------------------

    unsigned int n = Vmath::Imin(Exp[0]->GetNpoints(), Exp[0]->GetPhys(), 1);
    cout << "Minimim Jacobian: "
         << Vmath::Vmin(Exp[0]->GetNpoints(), Exp[0]->GetPhys(), 1) 
         << " at coords (" << x0[n] << ", " << x1[n] << ", " << x2[n] << ")"
         << endl;

    return 0;
}

