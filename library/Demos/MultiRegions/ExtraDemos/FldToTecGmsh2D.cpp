#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList2D.h>
#include <SpatialDomains/MeshGraph2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    MultiRegions::ExpList2DSharedPtr Exp;
    int i, j;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    if(argc != 3)
    {
#ifdef TECPLOT
        fprintf(stderr,"Usage: FldToTec2D  meshfile fieldfile\n");
#else
        fprintf(stderr,"Usage: FldToGmsh2D  meshfile fieldfile\n");
#endif
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph2D = MemoryManager<SpatialDomains::MeshGraph2D>::AllocateSharedPtr(vSession);
    //----------------------------------------------
    
    //----------------------------------------------
    // Import field file. 
    string fieldfile(argv[argc-1]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    LibUtilities::Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------

    //----------------------------------------------
    // Set up Expansion information
    vector< vector<LibUtilities::PointsType> > pointstype;
    for(i = 0; i < fielddef.size(); ++i)
    {         vector<LibUtilities::PointsType> ptype;
        for(j = 0; j < 2; ++j)
        {
            ptype.push_back(LibUtilities::ePolyEvenlySpaced);
        }
        pointstype.push_back(ptype);
    }
    graph2D->SetExpansions(fielddef,pointstype);
    //----------------------------------------------

    //----------------------------------------------        
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(vSession,graph2D);
    //----------------------------------------------  

    //----------------------------------------------
    // Copy data to file 
    for(int i = 0; i < fielddata.size(); ++i)
    {
        Exp->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_fields[0],
                                 Exp->UpdateCoeffs());
    }
    Exp->BwdTrans(Exp->GetCoeffs(),Exp->UpdatePhys());    
    //----------------------------------------------
    
    //----------------------------------------------
    // Write solution  depending on #define
#ifdef TECPLOT
    string   outfile(strtok(argv[argc-1],"."));
    string   endfile(".dat");
    outfile += endfile; 
    ofstream outstrm(outfile.c_str());
    
    Exp->WriteToFile(outstrm,eTecplot);
    outstrm.close();
#else
    string   outfile(strtok(argv[argc-1],"."));
    string   endfile(".pos");
    outfile += endfile; 
    ofstream outstrm(outfile.c_str());

    Exp->WriteToFile(outstrm,eGmsh);
    outstrm.close();
#endif
    //----------------------------------------------
    return 0;
}

