#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    MultiRegions::ExpList2DSharedPtr Exp;
    int     i, nq,  coordim;
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    if(argc != 3)
    {
        fprintf(stderr,"Usage: DatToFld meshfile file.dat\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    SpatialDomains::MeshGraph2D graph2D; 
    graph2D.ReadGeometry(meshfile);
    graph2D.ReadExpansions(meshfile);
    //----------------------------------------------


    //----------------------------------------------        
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(vSession,graph2D);
    //----------------------------------------------  

    //---------------------------------------------
    // Read Dat file
    string   in(argv[argc-1]);
    std::ifstream infile(in.c_str());
    Exp->ReadFromFile(infile,eTecplot);
    //---------------------------------------------

    //-----------------------------------------------
    // Write solution to file 
    string   out(strtok(argv[argc-1],"."));
    string   endfile(".fld");
    out += endfile; 
    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef = Exp->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size()); 
    
    for(i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("u");
        Exp->AppendFieldData(FieldDef[i], FieldData[i]);
    }
    graph2D.Write(out, FieldDef, FieldData);
    //-----------------------------------------------
    
    return 0;
}

