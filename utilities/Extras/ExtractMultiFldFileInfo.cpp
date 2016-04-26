#include <cstdio>
#include <cstdlib>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>  // for ASSERTL0
#include <SpatialDomains/MeshGraph.h>   // for FieldDefinitions, etc

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        fprintf(stderr,"Usage: ExtractMultiFileInfo file.fld \n" );
        exit(1);
    }

    //default meshgraph
    SpatialDomains::MeshGraph graph; 

    //----------------------------------------------
    // Import fieldfile.
    string fieldfile(argv[argc-1]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    LibUtilities::Import(fieldfile,fielddef,fielddata);

    //----------------------------------------------

    cout << "<MultipleFldFiles  FileName=\"" << fieldfile << "\" >";
    
    cout << fielddef[0]->m_elementIDs[0]; 

    for(int f = 0; f < fielddef.size(); ++f)
    {
        for(int i = 1; i < fielddef[f]->m_elementIDs.size(); ++i)
        {
            cout << "," <<fielddef[f]->m_elementIDs[i];
        }
    }
    cout << "<\\MultipleFldFiles>" << endl;
    return 0;
}

