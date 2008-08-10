#include "Gmsh.cpp"
#include "TetGen.cpp"


using namespace std;

class TiXmlDocument;

int main(int argc, char *argv[])
{

    using namespace Utilities;
    
    if(argc != 5)
    {
        fprintf(stderr,"Usage: ./MeshConvert -g/-t -o input(mesh)file output(xml)file\n");
        exit(1);
    }

    /*----------------------------------------------
    Read in gmesh from input file
    write the xml output file from input file */
    
    char* in(argv[1]);
    char* out(argv[2]);
    char* meshfile(argv[3]);
    char* outfile(argv[4]);

    string dashg = ("-g");
    string dasht = ("-t");
    string dasho = ("-o");
    string msh = ("msh");
    string smesh = ("smesh");

    string inputArg = (in);
    string outputArg = (out);
    string inputFile = (meshfile);
    size_t found;

    found = inputFile.find_last_of(".");
    string result = inputFile.substr(found+1);

    
     if(inputArg.compare(dashg) ==0 && outputArg.compare(dasho) == 0)
      {
      
            cout << "************* Gmsh ************************" << endl;
            cout << "Start parsing.......... " << meshfile <<  endl;

            if (result.compare(msh) ==0)
            {
                Gmsh::ParseGmshFile(meshfile, outfile);

                cout << " Success parsing and writing for Gmsh ........" << outfile <<  endl;
                
            }
            else
            {
                cerr << "Parsing fails : " << endl;
                cerr << "Input file " <<inputFile << " is not a correct MSH file type !  Please use *.msh file type "<< endl;
            }
       
      }
      else if(inputArg.compare(dasht) ==0 && outputArg.compare(dasho) == 0)
      {
      
            cout << "************* TetGen ************************" << endl;
            cout << "Start parsing.......... " << meshfile <<  endl;

            if (result.compare(smesh) ==0)
            {
                TetGen::ParseTetGen(meshfile, outfile);

                cout << " Success parsing and writing for TetGen........" << outfile <<  endl;
                
            }
            else
            {
                cerr << "Parsing fails : " << endl;
                cerr << "Input file " <<inputFile << " is not a correct SMESH file type !  Please use *.smesh file type "<< endl;
            }
       
      }
       else {

                cerr << "Input argument for Gmsh must " <<  dashg << ", but we got " <<inputArg << endl;

                cerr << "Input argument for TetGen must " << dasht << ", but we got " << inputArg <<  endl;

      }




    return 0;
}
