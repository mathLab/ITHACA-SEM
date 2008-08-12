#include "Gmsh.cpp"
#include "TetGen.cpp"
#include <stdexcept>

using namespace std;

class TiXmlDocument;

int main(int argc, char *argv[])
{

    using namespace Utilities;
    
    if(argc != 5)
    {
        fprintf(stderr,"Usage: ./MeshConvert -g/-t input(mesh)file  -o output(xml)file\n");
        exit(1);
    }

    
    /*----------------------------------------------
    Read in gmesh from input file
    write the xml output file from input file */
    
    char* in(argv[1]);
    char* meshfile(argv[2]);
    char* out(argv[3]);
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


     if(inputArg.compare(dashg) ==0)
      {
      
            cout << "************* Gmsh ************************" << endl;
            cout << "Start parsing.......... " << meshfile <<  endl;

            if (result.compare(msh) ==0)
            {
                Gmsh::ParseGmshFile(meshfile, outfile);

                cout << " Success parsing and writing " << meshfile << " to " << outfile <<  endl;
                
            }
            else
            {
                stringstream error;
                error << "Parsing fails : " << endl;
                error << "Input file " <<inputFile << " is not a correct MSH file type !  Please use *.msh file type "<< endl;
                throw runtime_error(error.str());
            }
       
      }
      else if(inputArg.compare(dasht) ==0)
      {
      
            cout << "************* TetGen ************************" << endl;
            cout << "Start parsing.......... " << meshfile <<  endl;

            if (result.compare(smesh) ==0)
            {
                TetGen::ParseTetGen(meshfile, outfile);

                cout << " Success parsing and writing " << meshfile << " to " << outfile <<  endl;
                
            }
            else
            {
                stringstream error;
                error << "Parsing fails : " << endl;
                error << "Input file " <<inputFile << " is not a correct SMESH file type !  Please use *.smesh file type "<< endl;
                throw runtime_error(error.str());
            }
       
      }
       else {

                stringstream error;
                
                error << "Input argument for Gmsh must " <<  dashg << ", but we got " <<inputArg << endl;
                error << "Input argument for TetGen must " << dasht << ", but we got " << inputArg <<  endl;

                throw runtime_error(error.str());

      }


    return 0;
}
