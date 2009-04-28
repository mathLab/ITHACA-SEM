#include <tinyxml/tinyxml.h>
#include <cstring>
#include <sstream>
#include <vector>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>

using namespace std;

class TiXmlDocument;

namespace Utilities
{

 namespace TetGen
  {
    struct Node
    {
        Node(int id, double x, double y, double z, int attr, int mark):
            id(id),
            x(x), y(y), z(z),
            attr(attr), mark(mark)
            {};

        int id;
        int attr, mark;
        double x, y, z;

    };

    struct Facet
    {
        Facet( int num_corner, vector<int> corners, int mark):
            num_corner(num_corner),
            corners(corners),
            mark(mark)
            {};

        int num_corner;
        vector<int> corners;
        int mark;
    };

    struct Hole
    {
        Hole(int hole_num, double x, double y, double z):
            hole_num(hole_num),
            x(x), y(y), z(z)
            {};

        int hole_num;
        double x, y, z;

    };

    struct Region
    {
        Region(int num, double x, double y, double z, int region_num, double attr):
            num(num),
            x(x), y(y), z(z),
            region_num(region_num),
            attr(attr)
            {};

        int num;
        double x, y, z;
        int region_num;
        double attr;
    };


    void WriteTetGenToXMLFile(const char* outfile, const vector<Node> & nodes, const vector<Facet> & facets);
    void ParseTetGen(const char* infile, const char* outfile);
    static vector<string> tokenize( const std::string& str, const std::string& delimiters = ", \n\r\t");
    static string removeComment( string const& line );

   } // end of namespace TetGen

}  // end of namespace Utilities



