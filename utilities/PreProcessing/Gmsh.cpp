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
  namespace Gmsh
  {
        struct Node
        {
            Node(int id, double x, double y, double z):
                id(id),
                x(x), y(y), z(z)
            {};
            int id;
            double x, y, z;
        };

        struct Element
        {
            Element(int id, int type, vector<int> tags, vector<int> nodes ):
                    id(id),
                    type(type),
                    tags(tags),
                    nodes(nodes)
            {};
            int id;
            int type;
            vector<int> tags;
            vector<int> nodes;
        };

        void ParseGmshFile(const char* inFile, const char* outfile);
        void WriteToXMLFile(const char* outfile, const vector<Node> & nodes, const vector<Element> & elements);

        void ParseGmshFile(const char* inFile, const char* outfile)
        {
            string line;
            vector<Node> nodes;
            vector<Element> elements;
            
            fstream mshFile(inFile);

            if(mshFile.is_open())
            {
                cout <<"Start reading gmsh..." << endl;
                while(!mshFile.eof())
                {
                    getline(mshFile, line);
                    stringstream s(line);
                    string word;
                    s >> word;            

                    if(word == "$Nodes")
                    {
                        getline(mshFile, line);
                        int num_nodes = 0;
                        stringstream s(line);
                        s >> num_nodes;
                        int id=0;
                        for(int i=0; i<num_nodes; ++i){
                            getline(mshFile, line);
                            stringstream st(line);
                            double x=0, y=0, z=0;
                            st >> id >> x >> y >> z;
                            nodes.push_back( Node(id, x, y, z) );
                        }               

                    }
                    else if(word == "$Elements")
                    {
                        getline(mshFile, line);
                        int num_elements = 0;
                        stringstream s(line);
                        s >> num_elements;
                        for(int i=0; i<num_elements; ++i)
                        {
                            getline(mshFile, line);
                            stringstream st(line);
                            int id=0, elm_type=0, num_tag=0, num_nodes=0;

                            st >> id >> elm_type >> num_tag;
                            switch(elm_type) {
                                case 1:  num_nodes = 2; break;
                                case 2:  num_nodes = 3; break;
                                case 3:  num_nodes = 4; break;
                                case 4:  num_nodes = 4; break;
                                case 5:  num_nodes = 8; break;
                                case 6:  num_nodes = 6; break;
                                case 7:  num_nodes = 5; break;
                                case 15: num_nodes = 1; break;
                            }
                            vector<int> tags;
                            for(int j=0; j<num_tag; ++j){
                                int tag=0;
                                st >> tag;
                                tags.push_back(tag);
                            }
                            vector<int> nodeList;
                            for(int k=0; k<num_nodes; ++k) {
                                int node = 0;
                                st >> node;
                                nodeList.push_back(node);
                            }

                            elements.push_back( Element(id, elm_type, tags, nodeList) );
                        }

                        cout << "num_elements = " << num_elements << endl;                                  
                    }

                }
                mshFile.close();

            WriteToXMLFile(outfile, nodes, elements);
            }

            else cout << "Unable to open file";
        }

        void WriteToXMLFile(const char* outfile, const vector<Node> & nodes, const vector<Element> & elements)
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
            
            verTag = new TiXmlElement( "ELEMENT" );
            
            comment = new TiXmlComment();
            comment->SetValue( " P - points, Q -  quads, T - triangles, S - segments, A - tetrahedron, Y - pyramid, R - prism, H- hex.");
            verTag->LinkEndChild( comment );

            comment = new TiXmlComment();
            comment->SetValue( "Only certain element types are appropriate for the given dimension (dim on mesh)");
            verTag->LinkEndChild( comment );

            comment = new TiXmlComment();
            comment->SetValue( " Can also use faces to define 3-D elements.  Specify with F[1] for face 1, for example. ");
            verTag->LinkEndChild( comment );
            
            for(int i=0; i<elements.size(); ++i){
                stringstream st;
                for(int j=0; j<elements[i].nodes.size(); ++j){
                    st << setw(5) << elements[i].nodes[j] << " ";
                }

                TiXmlElement *elm_tag;
                switch(elements[i].type) {
                case 1: elm_tag = new TiXmlElement("S"); break; // Segments
                case 2: elm_tag = new TiXmlElement("T"); break; // Triangle
                case 3: elm_tag = new TiXmlElement("Q"); break; // Quads
                case 4: elm_tag = new TiXmlElement("A"); break; // Tetrahedron
                case 5: elm_tag = new TiXmlElement("H"); break; // Hexahedron
                case 6: elm_tag = new TiXmlElement("R"); break; // Prism
                case 7: elm_tag = new TiXmlElement("Y"); break; //Pyramid
                case 15:elm_tag = new TiXmlElement("P"); break; // Points
                }

                elm_tag->SetAttribute("ID", elements[i].id);
                elm_tag->LinkEndChild( new TiXmlText(st.str()) );
                verTag->LinkEndChild(elm_tag);
            }
            
            geomTag->LinkEndChild( verTag );

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
        } // end of function WriteToXMLFile
        
     } // end of namespace Gmsh
     
  } // end of namespace Utilities
     
