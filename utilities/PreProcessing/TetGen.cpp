
#include <tinyxml/tinyxml.h>
#include <cstring>
#include <sstream>
#include <vector>

#include <iostream>
#include <fstream>
#include <iomanip>

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
    static string removeComment( string const& line );
    void ParseTetGen(const char* infile, const char* outfile);
    static vector<string> tokenize( const std::string& str, const std::string& delimiters = ", \n\r\t");

    // Make string to a token so that populate the points from data fule
    vector<string> tokenize( const string& str, const string& delimiters )
    {
        // Skip delimiters at beginning.
        string::size_type lastPos = str.find_first_not_of(delimiters, 0);
        // Find first "non-delimiter".
        string::size_type pos     = str.find_first_of(delimiters, lastPos);

        vector<string> tokens;
        while (pos != string::npos  ||  lastPos != string::npos)
        {
            // Found a token, add it to the vector.
            tokens.push_back(str.substr(lastPos, pos - lastPos));

            // Skip delimiters.  Note the "not_of"
            lastPos = str.find_first_not_of(delimiters, pos);

            // Find next "non-delimiter"        
            pos = str.find_first_of(delimiters, lastPos);
        }
        return tokens;
    }

    void WriteTetGenToXMLFile(const char* outfile, const vector<Node> & nodes, const vector<Facet> & facets)
    {  
        TiXmlDocument doc;      
        TiXmlDeclaration* decl = new TiXmlDeclaration( "1.0", "utf-8", "");  
        doc.LinkEndChild( decl );  
    
        TiXmlElement * root = new TiXmlElement( "NEKTAR" );  
        doc.LinkEndChild( root );  

        TiXmlComment * comment = new TiXmlComment();
        comment->SetValue(" Embed a 3-dimensional object in a 3-dimensional space");
        root->LinkEndChild( comment );
        comment = new TiXmlComment();
        comment->SetValue( "DIM <= SPACE  "
                            "this provides a method of optimizing code for a 3-D curve embedded in 3-d space.");
        root->LinkEndChild( comment );
    
        TiXmlElement * geomTag = new TiXmlElement( "GEOMETRY" );
        geomTag->SetAttribute("DIM", 3);
        geomTag->SetAttribute("SPACE", 3);    
        root->LinkEndChild( geomTag );

        comment = new TiXmlComment();
        comment->SetValue( " Definitions that can be used below in this file. ");
        geomTag->LinkEndChild( comment );

        TiXmlElement *def = new TiXmlElement("DEF");
        
        TiXmlElement * a = new TiXmlElement( "D" );
        a->LinkEndChild( new TiXmlText(" A = 1.0 ") );
        def->LinkEndChild(a);
        
        TiXmlElement * b = new TiXmlElement( "D" );
        b->LinkEndChild( new TiXmlText(" B = 2.0 ") );
        def->LinkEndChild(b);

        TiXmlElement * c = new TiXmlElement( "D" );
        c->LinkEndChild( new TiXmlText(" C = 3.0 ") );
        def->LinkEndChild(c);

        geomTag->LinkEndChild(def);

    
        TiXmlElement* verTag = new TiXmlElement( "VERTEX" );

        comment = new TiXmlComment();
        comment->SetValue( " Always must have four values per entry. ");
        verTag->LinkEndChild( comment );

        for( int i = 0; i < nodes.size(); ++i ) {
            stringstream s;
            s << fixed << showpoint << setprecision(15) << setw(20) << nodes[i].x << " " << nodes[i].y << " " << nodes[i].z;
            TiXmlElement * v = new TiXmlElement( "V" );
            v->SetAttribute("ID",nodes[i].id);
            v->LinkEndChild( new TiXmlText(s.str()) );
            verTag->LinkEndChild(v);
        }

        geomTag->LinkEndChild( verTag );
        
        TiXmlElement* elmTag = new TiXmlElement( "ELEMENT" );
        
        comment = new TiXmlComment();
        comment->SetValue( " P - points, Q -  quads, T - triangles, S - segments, A - tetrahedron, Y - pyramid, R - prism, H- hex.");
        elmTag->LinkEndChild( comment );

        comment = new TiXmlComment();
        comment->SetValue( "Only certain element types are appropriate for the given dimension (dim on mesh)");
        elmTag->LinkEndChild( comment );

        comment = new TiXmlComment();
        comment->SetValue( " Can also use faces to define 3-D elements.  Specify with F[1] for face 1, for example. ");
        elmTag->LinkEndChild( comment );
        
        for(int i=0; i<facets.size(); ++i){
            stringstream st;
            for(int j=0; j<facets[i].corners.size(); ++j){
                st << setw(5) << facets[i].corners[j] << " ";
            }

            TiXmlElement *facet_tag;// = new TiXmlElement("T");
            switch(facets[i].corners.size()) {
            case 1: facet_tag = new TiXmlElement("P"); break; // points
            case 2: facet_tag = new TiXmlElement("S"); break; // segments
            case 3: facet_tag = new TiXmlElement("T"); break; // triangle
            case 4: facet_tag = new TiXmlElement("A"); break; // Tetrahedron
            }

            facet_tag->SetAttribute("ID", i);
            facet_tag->LinkEndChild( new TiXmlText(st.str()) );
            elmTag->LinkEndChild(facet_tag);
        }
        
        geomTag->LinkEndChild( elmTag );

        comment = new TiXmlComment();
        comment->SetValue( " V - vertex, S - segment, E - edge, F - face ");
        geomTag->LinkEndChild( comment );

        TiXmlElement * comp = new TiXmlElement ("COMPOSITE" );
        comp->LinkEndChild( new TiXmlText( "Set the composite here...." ));
        geomTag->LinkEndChild( comp );
            
        TiXmlElement * domain = new TiXmlElement ("DOMAIN" );
        domain->LinkEndChild( new TiXmlText( "Set the domain here...." ));
        geomTag->LinkEndChild( domain );

        TiXmlElement * exp = new TiXmlElement( "EXPANSIONS" );
        exp->LinkEndChild( new TiXmlText( "Set the expansion here...." ));
        root->LinkEndChild( exp );
    
        TiXmlElement * conds = new TiXmlElement( "CONDITIONS" );    
        root->LinkEndChild( conds );

        comment = new TiXmlComment();
        comment->SetValue( "Removed redundancy since we can specify any lefel of granularity in the ExpansionTypes section below");
        conds->LinkEndChild( comment );
    
        TiXmlElement* par = new TiXmlElement( "PARAMETERS" );
        par->LinkEndChild( new TiXmlText( "Set the parameters here..." ));
        conds->LinkEndChild( par );

        comment = new TiXmlComment();
        comment->SetValue( "One of these for each dimension.  These are the vector components, say, s = (u, v); comprised of two components in this example for a 3D dimension.");
        conds->LinkEndChild( comment );
    
        TiXmlElement* var = new TiXmlElement( "VARIABLES" );
        var->LinkEndChild( new TiXmlText( "Set the variables here..." ));
        conds->LinkEndChild( var );

        comment = new TiXmlComment();
        comment->SetValue( "These composites must be defined in the geometry file.");
        conds->LinkEndChild( comment );


        TiXmlElement* bregion = new TiXmlElement( "BOUNDARYREGIONS" );
        bregion->LinkEndChild( new TiXmlText( "Set the boundary region here..." ));
        conds->LinkEndChild( bregion );

        TiXmlElement* bcond = new TiXmlElement( "BOUNDARYCONDITIONS" );
        conds->LinkEndChild( bcond );

        comment = new TiXmlComment();
        comment->SetValue( "The region numbers below correspond to the regions specified in the Boundary Region definition above.");
        conds->LinkEndChild( comment );
        
        TiXmlElement* region = new TiXmlElement( "REGION" );
        region->LinkEndChild( new TiXmlText( "Set the boundary region here..." ));
        bcond->LinkEndChild(region);

        TiXmlElement* forcing = new TiXmlElement( "FORCING" );
        forcing->LinkEndChild( new TiXmlText( "Set the forcing function here..." ));
        conds->LinkEndChild( forcing );
    
        TiXmlElement* exact = new TiXmlElement( "EXACTSOLUTION" );
        exact->LinkEndChild( new TiXmlText( "Set the exact solution here..." ));
        conds->LinkEndChild( exact );

        doc.SaveFile(outfile );
    }

    // Enable to remove a comment line
    string removeComment( string const& line ) {
        string::size_type position = line.find_first_of('#');
        return line.substr( 0, position );
    }

    void ParseTetGen(const char* infile, const char* outfile)
    {
        string line, originalLine;
        vector<Node> nodes;
        vector<Facet> facets;
        vector<Hole> holes;
        vector<Region> regions;
        stringstream ss;

        int num_nodes = 0, dim = 0, attr = 0, mark =0;
        int num_facets = 0, num_holes=0, num_regions=0;
        
        fstream mshFile(infile);
        mshFile.exceptions ( fstream::badbit );

        try{
            if(mshFile.is_open())
            {
                cout <<"Start reading gmsh..." << endl;

                int lineNumber = 0, parsedLineNumber = 0;
                
                while(!mshFile.eof())
                {
                    getline(mshFile, originalLine);
                    line = removeComment(originalLine);
                    vector<string> tokens = tokenize(line);

                    if(parsedLineNumber ==0){
                        stringstream st(originalLine);
                        
                    if(tokens.size() == 2){
                        st << tokens[0] << " " << tokens[1] << endl;
                        st >> num_nodes >> dim ;
                        }
                        if(tokens.size() == 4){
                        st << tokens[0] << " " << tokens[1] << " " << tokens[2] << " " << tokens[3] << endl;
                        st >> num_nodes >> dim >> attr >> mark;
                        } else if(!tokens.empty()){
                            cerr << "Unrecognized file format on line " << lineNumber << ": " << originalLine << endl;
                            cerr << "Expected: num_nodes, dim, attr, mark" << endl;
                                        
                        } 
                    } else if((0 < parsedLineNumber) && (parsedLineNumber<= num_nodes) ){ // parse nodes
                    
                        cout << "parsedLineNumber= " << parsedLineNumber << endl;

                        int id=0, nattr=0, nmark=0;
                        double x=0, y=0, z=0;
                        // vector<Node> nodes;
                        stringstream st(originalLine);

                            if(!tokens.empty()){
                                if(tokens.size() == 4){
                                st << tokens[0] << " " << tokens[1] << " " << tokens[2] << " " << tokens[3] << endl;
                                st >> id >> x >> y >> z ;
                                }
                                if(tokens.size() == 5){
                                st << tokens[0] << " " << tokens[1] << " " << tokens[2] << " " << tokens[3] <<" "<< tokens[4] << endl;
                                st >> id >> x >> y >> z >> nattr;
                                }
                                else if(tokens.size() == 6)
                                {
                                st << tokens[0] << " " << tokens[1] << " " << tokens[2] << " " << tokens[3] <<
                                " " << tokens[4] << " " << tokens[5] << endl;
                                st >> id >> x >> y >> z >> nattr >> nmark;
                                }
                                else if( !tokens.empty() )
                                {
                                cerr << "Wrong number of coordinates provided on line " << lineNumber << ": "
                                << originalLine << endl;
                                }
                                
                            }
                        nodes.push_back( Node(id, x, y, z, nattr, nmark) );

                        cout << "id = " << id << "  x = " << x <<  "  y = " << y <<   "  z = " << z
                            << "attr = " << nattr << "  mark = " << nmark << endl;
                        cout << "/-------------------------- Done parsing nodes ------------------ " << endl;
                        
                    } else if(parsedLineNumber == num_nodes+1){ // parse facet info

                    cout << "parsedLineNumber = " << parsedLineNumber << endl; 
                        stringstream st(originalLine);
                        
                        if(tokens.size() == 2) {
                            st << tokens[0] << " " << tokens[1] << endl;
                            st >> num_facets >> mark;
                        } else if(!tokens.empty()){
                            cerr << "Unrecognized file format on line " << lineNumber << ": " << originalLine << endl;
                            cerr << "Expected: num_nodes, dim, attr, mark" << endl;                                  
                        }

                        cout << "num_facets= " << num_facets << "  mark = " << mark << endl;
                        
                    }   // parse facet
                    else if ((parsedLineNumber > num_nodes+1) && (parsedLineNumber <= num_facets+num_nodes+1)) 
                    {
                        cout << "/***************** parsing facet ------------------ " << endl;
                        cout << "parsedLineNumber = " << parsedLineNumber << endl; 

                        int id=0, nmark=0;
                        int num_corners;
                        // vector<Facet> facets;
                        vector<int> corners;
                        stringstream st(originalLine);

                        if(!tokens.empty()){

                            st << tokens[0] << " ";
                            st >> num_corners;

                            int tokenSize = size_t(num_corners+2);
                            
                            if((tokens.size() == tokenSize)||(tokens.size() == tokenSize-1) ){
                                                    
                                for(int i=0; i<num_corners; ++i){
                                
                                    int corner = 0;
                                    st << tokens[i+1] << " ";
                                    st >> corner;
                                    corners.push_back(corner);
                                    cout << "corners[" << i << "] = " << corners[i] << endl;
                                }

                            } else {
                                cerr << "Unrecognized file format on line " << lineNumber << ": " << originalLine << endl;
                                cerr << "--Expected: tokens = <list of " << num_corners+2 << "numbers>" << endl;
                            }

                        }
                        st >> nmark;

                        facets.push_back( Facet(num_corners, corners, nmark) );

                        cout  << "id = " << id << "  num_corners = " << num_corners <<"  mark = " << nmark << endl;                               
                    }
                    else if(parsedLineNumber == num_nodes+num_facets+2)
                    {
                        cout << "*********parsing holes***************" << endl;
                        cout << "parsedLineNumber = " << parsedLineNumber << endl;
                        cout << "line number = " << lineNumber << endl;
                        stringstream st(originalLine);
                        
                        if(tokens.size() == 1) {
                            st << tokens[0] << " ";
                            st >> num_holes;
                        } else if(!tokens.empty()){
                            cerr << "Unrecognized file format on line " << lineNumber << ": " << originalLine << endl;
                            cerr << "Expected: num_holes" << endl;                                  
                        }

                        cout << "num_holes= " << num_holes << endl;
                        
                    }
                    else if((parsedLineNumber > num_nodes+num_facets+2) && (parsedLineNumber <= num_nodes+num_facets+num_holes+2)) {
                                
                        cout << "parsedLineNumber= " << parsedLineNumber << endl;

                        int id=0;
                        double x=0, y=0, z=0;
                        //vector<Hole> holes;
                        stringstream st(originalLine);

                            if(!tokens.empty()){
                                if(tokens.size()==4){
                                st << tokens[0] << " " << tokens[1] << " " << tokens[2] << " " << tokens[3] << endl;
                                st >> id >> x >> y >> z ;
                                }
                                else if( !tokens.empty() )
                                {
                                cerr << "Wrong number of coordinates provided on line " << lineNumber << ": "
                                << originalLine << endl;
                                }
                                
                            }
                        holes.push_back( Hole(id, x, y, z));

                        cout << "id = " << id << "  x = " << x <<  "  y = " << y <<   "  z = " << z << endl;
                        cout << "/-------------------------- Done parsing holes ------------------ " << endl;
                        
                    }
                    else if(parsedLineNumber == num_nodes+num_facets+num_holes+3)
                    {
                        cout << "*********parsing region***************" << endl;
                        cout << "parsedLineNumber = " << parsedLineNumber << endl;               
                        stringstream st(originalLine);
                        
                        if(tokens.size() == 1) {
                            st << tokens[0] << " ";
                            st >> num_regions;
                        } else if(!tokens.empty()){
                            cerr << "Unrecognized file format on line " << lineNumber << ": " << originalLine << endl;
                            cerr << "Expected: num_holes" << endl;                                  
                        }

                        cout << "num_regions= " << num_regions << endl;
                        
                    }
                    else if((parsedLineNumber > num_nodes+num_facets+num_holes+3) &&
                            (parsedLineNumber <= num_nodes+num_facets+num_holes+num_regions+3)) {
                                
                        cout << "parsedLineNumber= " << parsedLineNumber << endl;

                        int id=0, region_num=0;
                        double rattr=0;
                        double x=0, y=0, z=0;
                        //vector<Region> regions;
                        stringstream st(originalLine);

                            if(!tokens.empty()){
                                if(tokens.size()==6){
                                st << tokens[0] << " " << tokens[1] << " " << tokens[2] << " " << tokens[3]
                                << " " << tokens[4] << " " << tokens[5] << endl;
                                st >> id >> x >> y >> z >> region_num >> rattr ;
                                }
                                else if( !tokens.empty() )
                                {
                                cerr << "Wrong number of coordinates provided on line " << lineNumber << ": "
                                << originalLine << endl;
                                }
                                
                            }
                        regions.push_back( Region(id, x, y, z, region_num, rattr));
                    
                        cout << "id = " << id << "  x = " << x <<  "  y = " << y <<   "  z = " << z
                            << "  region_num = " << region_num << "  reg-attr = " << rattr << endl;
                        cout << "/-------------------------- Done parsing Regions ------------------ " << endl;
                        
                    }
                    
                    if(!tokens.empty()){            
                        ++parsedLineNumber;
                    }
                    ++lineNumber;
                }
                    
            } else {
            
                cerr << "Unable to open file " << endl;        
            }
            }
            catch (fstream::failure e){
                cerr << "Exception opening/reading file: \n" << e.what() << endl;
            }

        mshFile.close();
        
     WriteTetGenToXMLFile(outfile, nodes, facets);    
    } // end of ParseTetGen function

   } // end of namespace TetGen

}  // end of namespace Utilities



