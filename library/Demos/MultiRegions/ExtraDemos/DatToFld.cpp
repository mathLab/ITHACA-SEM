#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList2D.h>
#include <SpatialDomains/MeshGraph2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    MultiRegions::ExpList2DSharedPtr Exp;
    int i;
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    if(argc != 3)
    {
        fprintf(stderr,"Usage: DatToFld meshfile file.dat\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    SpatialDomains::MeshGraphSharedPtr graph2D = MemoryManager<SpatialDomains::MeshGraph2D>::AllocateSharedPtr(vSession);
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
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef = Exp->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size()); 
    
    for(i = 0; i < FieldDef.size(); ++i)
    {
        FieldDef[i]->m_fields.push_back("u");
        Exp->AppendFieldData(FieldDef[i], FieldData[i]);
    }
    LibUtilities::Write(out, FieldDef, FieldData);
    //-----------------------------------------------
    
    return 0;
}

