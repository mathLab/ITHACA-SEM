#include <PreProcessing/Gmsh.h>

#define ERROR(msg)                              \
    string strmsg(msg);                         \
    cerr << strmsg << endl;                     \
    throw runtime_error(strmsg);

namespace Utilities
{
  namespace Gmsh
  {
      
    void ParseGmshFile(const char* inFile, const char* outfile)
    {
      string line;
      int nVertices        = 0;
      int nEntities        = 0;
      int nBoundComposites = 0;
      int expDim           = 3;
     
      vector<Vertex> vertices;
      vector<Edge> edges;
      vector<Face> faces;
      vector<ZeroDElement> zeroDElements;
      vector<OneDElement> oneDElements;
      vector<TwoDElement> twoDElements; 
      vector<ThreeDElement> threeDElements; 
      vector<Composite> composites;
     

      //---------------------------------------------
      // Read the *.msh file and fill the 
      // nodes and geometric enteties structures
      
      fstream mshFile(inFile);
      
      if(!mshFile.is_open())
	{
	  ERROR("Unable to find msh file");
	}
      else
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
		  stringstream s(line);
		  s >> nVertices;
		  int id=0;
		  for(int i=0; i<nVertices; ++i)
		    {
		      getline(mshFile, line);
		      stringstream st(line);
		      double x=0, y=0, z=0;
		      st >> id >> x >> y >> z;
		      id -= 1; // counter starts at 0 
		      vertices.push_back( Vertex(id, x, y, z) );
		    }               
		}
	      else if(word == "$Elements")
		{
		  int zeroDid=0, oneDid=0, twoDid=0, threeDid=0;
		  getline(mshFile, line);
		  stringstream s(line);
		  s >> nEntities;
		  for(int i=0; i<nEntities; ++i)
		    {
		      getline(mshFile, line);
		      stringstream st(line);
		      int id=0, elm_type=0, num_tag=0, num_nodes=0;
		      
		      st >> id >> elm_type >> num_tag;
		      id -= 1; // counter starts at 0 
		      
		      num_nodes = GetNnodes(elm_type);
		      
		      vector<int> tags;
		      for(int j=0; j<num_tag; ++j)
			{
			  int tag=0;
			  st >> tag;
			  tags.push_back(tag);
			  
			  // physical entities (used to define boundary entities)
			  // are stored in tags[0] 
			  nBoundComposites = max(nBoundComposites,tags[0]);
			}
		      
		      vector<int> edgeList;
		      vector<int> nodeList;
		      for(int k=0; k<num_nodes; ++k) 
			{
			  int node = 0;
			  st >> node;
			  node -= 1; // counter starts at 0 
			  nodeList.push_back(node);
			}
		      
		      switch(elm_type)
			{
			case 15:
			  zeroDElements.push_back( ZeroDElement(zeroDid++,elm_type, tags, nodeList) );
			  break;
			case 1:
			  oneDElements.push_back( OneDElement(oneDid++,elm_type, tags, nodeList) );
			  break;
			case 2:
			case 3:
			  twoDElements.push_back( TwoDElement(twoDid++,elm_type, tags, nodeList, edgeList) );
			  break;
			default:
			  {
			    ERROR("Gmsh element type not yet supported");
			  }
			}
		    }
		  
		  // find out what expansion dimension
		  if (threeDElements.empty() && twoDElements.empty() && oneDElements.empty())
		    {
		       ERROR("Illegal expansion dimension");
		     }
		  else if (threeDElements.empty() && twoDElements.empty())
		    expDim = 1;
		  else if (threeDElements.empty())
		    expDim = 2;
		  else 
		    expDim = 3;
		  
		 		  
		  cout << "Expansion dimension is " << expDim << endl;
		  cout << "Read " << nVertices << " vertices" << endl;
		  cout << "Read " << nEntities << " geometric entities" << endl;
		  cout << "Read " << nBoundComposites << " boundary entities"  << endl;                                  
		}
	      
	    }
	  mshFile.close();
	  //---------------------------------------------
	}	  

      switch(expDim)
	{
	case 1:
	  SortZeroDElements(zeroDElements,vertices);
	  SortOneDComposites(oneDElements, zeroDElements, composites, nBoundComposites);
	  break;
	case 2:
	  SortEdgeToVertex(twoDElements,edges);
	  SortOneDElements(oneDElements,edges);
	  SortTwoDComposites(twoDElements, oneDElements, composites, nBoundComposites);
	  break;
	case 3:
	  {
	    ERROR("3D not supported yet");
	  }
	  break;
	default:
	  {
	    ERROR("illegal expansion dimension");
	  }
	}

      WriteToXMLFile(outfile, expDim, vertices, edges, faces, oneDElements, twoDElements, threeDElements, composites);

    }
    

    // note: this function is only working for 2nodes edges.
    // Any of Gmsh higher order edges are not supported 
    int GetEdge(vector<int> &vert, vector<Edge>& edges)
    {
      int i;
      int size = edges.size();
      for(i = 0; i < size; i++)
	{
	  if( ((vert[0] == edges[i].vert[0])&&
	       (vert[1] == edges[i].vert[1])) 
                    ||
	      ((vert[1] == edges[i].vert[0])&&
	       (vert[0] == edges[i].vert[1])) )
	    {
	      return i;
	    }
	}
      
      // for the last edge of an element
      vector<int> edgevert;
      edgevert.push_back(vert[0]);
      edgevert.push_back(vert[1]);
      edges.push_back( Edge(size,edgevert) );
      return size;
    }

    

    void SortEdgeToVertex(vector<TwoDElement> & elements, vector<Edge> & edges)
    {     
      int i, j;
      int size = elements.size();
      
      // fill the edge vector in the elements struct with unique edges
      for (i = 0; i < size; ++i)
	{
	  for (j = 0; j < elements[i].vert.size()-1; ++j)
	    {
	      int edgeid;
	      vector<int> vert;
	      vert.push_back(elements[i].vert[j]);
	      vert.push_back(elements[i].vert[j+1]);
	      edgeid = GetEdge(vert,edges);
	      elements[i].edge.push_back(edgeid);
	    }
	  
	  {
	    int edgeid;
	    vector<int> vert;
	    vert.push_back(elements[i].vert[elements[i].vert.size()-1]);
	    vert.push_back(elements[i].vert[0]);
	    edgeid = GetEdge(vert,edges);
	    elements[i].edge.push_back(edgeid);
	  }
	}
      cout << "...done sorting EdgeToVertex relations" << endl;
    } 



    // Get the correct edge id's from the edge struct
    void SortZeroDElements(vector<ZeroDElement> & points,const vector<Vertex> & vertices)
    {
      int i, j;
      int vertsize = vertices.size();
      int pointsize  = points.size();
      for (i = 0; i < pointsize; ++i)
	{
	  for(j = 0; j < vertsize; j++)
	    {
	      if( points[i].vert[0] == vertices[j].x )
		{
		  points[i].id = vertices[j].id;
		  break;
		}
	    }
	}
      cout << "...done sorting points ids" << endl;

    }
 
    
 
    // Get the correct edge id's from the edge struct
    void SortOneDElements(vector<OneDElement> &segments,const vector<Edge>& edges)
    {
      int i, j;
      int edgesize = edges.size();
      int segsize  = segments.size();
      for (i = 0; i < segsize; ++i)
	{
	  for(j = 0; j < edgesize; j++)
	    {
	      if( ((segments[i].vert[0] == edges[j].vert[0])&&
		   (segments[i].vert[1] == edges[j].vert[1])) 
		  ||
		  ((segments[i].vert[1] == edges[j].vert[0])&&
		   (segments[i].vert[0] == edges[j].vert[1])) )
		{
		  segments[i].id = edges[j].id;
		  break;
		}
	    }
	}
      cout << "...done sorting segment ids" << endl;

    }
 
    
 
    void SortTwoDComposites(const vector<TwoDElement> & elements, const vector<OneDElement> & edges,
			    vector<Composite> & composites, int num_composites)
    {
      int i;

      // set elements composite 
      // for these composites we do not store
      // the individual element ids (as the composites
      // still will be defined from 0 to nElm-1), we just store the 
      // total number of tris or quads 
      int nTri  = 0;
      int nQuad = 0;
      int twoDcomp = 0;
      for (int i = 0; i < elements.size(); ++i)
	{
	  if (elements[i].type == 2)
	    nTri++;
	  if (elements[i].type == 3)
	    nQuad++;
	}
      
      if (nTri != 0)
	{
	  list<int> eid;
	  eid.push_back(nTri-1); // -1 since start from 0
	  composites.push_back(Composite(twoDcomp,2,eid));
	  twoDcomp++;
	}
      if (nQuad != 0)
	{ 
	  list<int> eid;
	  eid.push_back(nQuad-1); // -1 since start from 0
	  composites.push_back(Composite(twoDcomp,3,eid));
	  twoDcomp++;
	}
      
      int comp;
      
      // Initialize edge composite
      for (i = twoDcomp; i < num_composites+twoDcomp; ++i)
	{	
	  list<int> eid;
	  composites.push_back(Composite(i,0,eid));
	}
      
      // Boundary comp are stored as physical entities and are stored in elements.tags[0] 
      // should be numbered from 1 to num_composites (zero mean not included in composite)
      for (i = 0; i < edges.size(); ++i)
	{
	  comp = edges[i].tags[0];
	  if (comp != 0)
	    {
	      // need to sort 
	      switch(edges[i].type)
		{
		case 1:
		  {
		    composites[comp-1+twoDcomp].eid.push_back(edges[i].id);
		    composites[comp-1+twoDcomp].type = edges[i].type;
		  }
		  break;
		}
	    }
	}
      
      // sort the eid list 
      // composites must be increasing numbers...
      for (i = 0; i < num_composites+twoDcomp; ++i)
	{	
	  composites[i].eid.sort();
	}
      
      cout << "...done sorting composites" << endl;

    }
    

    void SortOneDComposites(const vector<OneDElement> & elements, const vector<ZeroDElement> & points,
			  vector<Composite> & composites, int num_composites)
    {
      int i;

      // set elements composite 
      // for these composites we do not store
      // the individual element ids (as the composites
      // still will start from 0), we just store the 
      // total number of tris or quads 
      int nSeg  = elements.size();
      list<int> eid;
      eid.push_back(nSeg-1); // -1 since start from 0
      composites.push_back(Composite(0,1,eid));
           
      int comp;
      
      // Initialize edge composite
      for (i = 1; i < num_composites+1; ++i)
	{	
	  list<int> eid;
	  composites.push_back(Composite(i,0,eid));
	}
      
      // Boundary comp are stored as physical entities and are stored in elements.tag[0] 
      // should be numbered from 1 to num_composites (zero mean not included in composite)
      for (i = 0; i < points.size(); ++i)
	{
	  comp = points[i].tags[0];
	  if (comp != 0)
	    {
	      // need to sort 
	      switch(points[i].type)
		{
		case 15:
		  {
		    composites[comp].eid.push_back(points[i].id);
		    composites[comp].type = points[i].type;
		  }
		  break;
		}
	    }
	}
      
      // sort the eid list 
      // composites must be increasing numbers...
      for (i = 0; i < num_composites+1; ++i)
	{	
	  composites[i].eid.sort();
	}
      
      cout << "...done sorting composites" << endl;

    }

    
    void WriteToXMLFile(const char* outfile, int expDim, const vector<Vertex> & nodes, const vector<Edge> & edges, 
			const vector<Face> & faces, const vector<OneDElement> & oneDElements,
			const vector<TwoDElement> & twoDElements, const vector<ThreeDElement> & threeDElements,
			const vector<Composite> & composites)
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
       
      //---------------------------------------------
      // Write DIM and SPACE
      
      TiXmlElement * geomTag = new TiXmlElement( "GEOMETRY" );
      geomTag->SetAttribute("DIM", expDim);
      geomTag->SetAttribute("SPACE", expDim);    
      root->LinkEndChild( geomTag );
      
      comment = new TiXmlComment();
      comment->SetValue( " Definitions that can be used below in this file. ");
      geomTag->LinkEndChild( comment );
      //---------------------------------------------

      
      //---------------------------------------------
      // Write VERTEX 
      
      TiXmlElement* verTag = new TiXmlElement( "VERTEX" );
      
      comment = new TiXmlComment();
      comment->SetValue( " Always must have four values per entry. ");
      verTag->LinkEndChild( comment );
      
      for( int i = 0; i < nodes.size(); ++i ) {
	stringstream s;
	//  s << fixed << showpoint << setprecision(15) << setw(20) << nodes[i].x << " " << nodes[i].y << " " << nodes[i].z;
	s << scientific << setprecision(3) <<  nodes[i].x << " " << nodes[i].y << " " << nodes[i].z;
	TiXmlElement * v = new TiXmlElement( "V" );
	v->SetAttribute("ID",nodes[i].id);
	v->LinkEndChild( new TiXmlText(s.str()) );
	verTag->LinkEndChild(v);
      }
      geomTag->LinkEndChild( verTag );
      //---------------------------------------------
      
      
      //---------------------------------------------
      // Write EDGE
      
      // only if expansion dimension >= 2
      if (expDim >= 2)
	{
	  int edgecnt = 0;
	  int small, large;
	  verTag = new TiXmlElement( "EDGE" );
	  
	  comment = new TiXmlComment();
	  comment->SetValue( "Edges are vertex pairs ");
	  verTag->LinkEndChild( comment );
	  
	  for( int i = 0; i < edges.size(); ++i ) {
	    stringstream s;
	    
	    s << setw(5) << edges[i].vert[0] << "  " << edges[i].vert[1] << "   ";
	    TiXmlElement * e = new TiXmlElement( "E" );
	    e->SetAttribute("ID",edgecnt++);
	    e->LinkEndChild( new TiXmlText(s.str()) );
	    verTag->LinkEndChild(e);
	  }
	  geomTag->LinkEndChild( verTag );
	}					
      //---------------------------------------------
      
      
      //--------------------------------------------
      // Write FACES
      // only if expansion dimension = 3
      //  if (expDim == 3)
      // 
      //--------------------------------------------
      
      
      //--------------------------------------------
      // Write ELEMENT
      
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
      
      switch(expDim){
	// for 1D elements are Segments 
      case 1:
	{
	  for(int i=0; i<oneDElements.size(); ++i){
	    if (oneDElements[i].type == 1){ // 1 denotes Seg
	      
	      stringstream st;
	      for(int j=0; j<oneDElements[i].vert.size(); ++j){
		st << setw(5) << oneDElements[i].vert[j] << " ";
	      }
	      
	      TiXmlElement *elm_tag;
	      elm_tag = new TiXmlElement("S"); 
	      elm_tag->SetAttribute("ID", oneDElements[i].id);
	      elm_tag->LinkEndChild( new TiXmlText(st.str()) );
	      verTag->LinkEndChild(elm_tag);
	    }
	  }
	}
	break;
	// for 2D elements are Tri or Quads
      case 2:
	{
	  int elmcnt = 0;
	  
	  for(int i=0; i<twoDElements.size(); ++i){
	    if (twoDElements[i].type == 2 || twoDElements[i].type == 3){ // 2 denotes Tri/3 denotes Quad
	      
	      stringstream st;
	      for(int j=0; j<twoDElements[i].edge.size(); ++j){
		st << setw(5) << twoDElements[i].edge[j] << " ";
	      }
	      
	      TiXmlElement *elm_tag;
	      if (twoDElements[i].type == 2) {
		elm_tag = new TiXmlElement("T"); 
	      }
	      else if (twoDElements[i].type == 3) {
		elm_tag = new TiXmlElement("Q");
	      }
	      else{
		cout << "No valid 2D element" << endl;
		exit(1);
	      }
	      
	      elm_tag->SetAttribute("ID", elmcnt++);
	      elm_tag->LinkEndChild( new TiXmlText(st.str()) );
	      verTag->LinkEndChild(elm_tag);
	    }
	  }
	}
	break;
	// for 3D elements are Tet, Hex, Prism or Pyr
      case 3:
	{
	  int elmcnt = 0;
	  
	  for(int i=0; i<threeDElements.size(); ++i){
	    if (threeDElements[i].type == 4 || threeDElements[i].type == 5 ||
		threeDElements[i].type == 6 || threeDElements[i].type == 7    ){ 
	      
	      stringstream st;
	      for(int j=0; j< threeDElements[i].face.size(); ++j){
		st << setw(5) << threeDElements[i].face[j] << " ";
	      }
	      
	      TiXmlElement *elm_tag;
	      switch(threeDElements[i].type) {
	      case 4: elm_tag = new TiXmlElement("A"); break; // Tetrahedron
	      case 5: elm_tag = new TiXmlElement("H"); break; // Hexahedron
	      case 6: elm_tag = new TiXmlElement("R"); break; // Prism
	      case 7: elm_tag = new TiXmlElement("Y"); break; //Pyramid
	      default: cout << "No valid 3D element" << endl; exit(1);
	      }
	      
	      
	      elm_tag->SetAttribute("ID", elmcnt++);
	      elm_tag->LinkEndChild( new TiXmlText(st.str()) );
	      verTag->LinkEndChild(elm_tag);
	    }
	  }
	}
	break;
      default:
	{
	  ERROR("Dimension > 3");
	}
      }
      
      geomTag->LinkEndChild( verTag );
      //--------------------------------------------------
      
      
      //--------------------------------------------------
      // Write COMPOSITE
      
      comment = new TiXmlComment();
      comment->SetValue( " V - vertex, S - segment, E - edge, F - face ");
      geomTag->LinkEndChild( comment );
      
      verTag = new TiXmlElement("COMPOSITE");
      
      list<int>::iterator it;
      list<int> compList;
      
      for(int i=0; i<composites.size(); ++i){
	stringstream st,st_start,st_end;
	compList = composites[i].eid;
	
	for (it = compList.begin(); it != compList.end(); ++it){
	  st  << *it;
	  if ( *it == compList.back())
	    st << "]";
	  else
	    st << ",";
	}
	
	TiXmlElement *comp_tag = new TiXmlElement("C"); // Composite
	TiXmlElement *elm_tag;
	
	switch(composites[i].type) 
	  {
	  case 1:
	    { 
	      if (expDim == 1)
		st_start << " S[0-"; // Segment
	      else
		st_start << " E[";  // Segment -> Edges
	    }
	    break; 
	  case 2: 
	    {
	      if (expDim == 2)
		st_start << " T[0-";  // Triangle
	      else
		st_start << " F[";    // Triangle -> Face
	    }
	    break; 
	  case 3:
	    {
	      if (expDim == 2)
		st_start << " Q[0-";  // Quads
	      else
		st_start << " F[";    // Quads -> Face
	    }
	    break; // Quads
	  case 4: st_start << " A[0-"; break;   // Tetrahedron
	  case 5: st_start << " H[0-"; break;   // Hexahedron
	  case 6: st_start << " R[0-"; break;   // Prism
	  case 7: st_start << " Y[0-"; break;   // Pyramid
	  case 15:st_start << " V["; break;   // Points -> Vertex
	  }
	
	comp_tag->SetAttribute("ID", composites[i].id);
	comp_tag->LinkEndChild( new TiXmlText(st_start.str()) );
	comp_tag->LinkEndChild( new TiXmlText(st.str()) );
	verTag->LinkEndChild(comp_tag);
      }
      
      geomTag->LinkEndChild( verTag );
      //--------------------------------------------------
      
      
      
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
    


    
    int GetNnodes(int GmshEntity)
    {
      int nNodes;
      
      switch(GmshEntity) 
	{
	case 1:  nNodes = 2;  break;
	case 2:  nNodes = 3;  break;
	case 3:  nNodes = 4;  break;
	case 4:  nNodes = 4;  break;
	case 5:  nNodes = 8;  break;
	case 6:  nNodes = 6;  break;
	case 7:  nNodes = 5;  break;
	case 8:  nNodes = 3;  break;
	case 9:  nNodes = 6;  break;
	case 10: nNodes = 9;  break;
	case 11: nNodes = 10; break;
	case 12: nNodes = 27; break;
	case 13: nNodes = 18; break;
	case 14: nNodes = 14; break;
	case 15: nNodes = 1;  break;
	case 16: nNodes = 8;  break;
	case 17: nNodes = 20; break;
	case 18: nNodes = 15; break;
	case 19: nNodes = 13; break;
	case 20: nNodes = 9;  break;
	case 21: nNodes = 10; break;
	case 22: nNodes = 12; break;
	case 23: nNodes = 15; break;
	case 24: nNodes = 15; break;
	case 25: nNodes = 21; break;
	case 26: nNodes = 4;  break;
	case 27: nNodes = 5;  break;
	case 28: nNodes = 6;  break;
	case 29: nNodes = 20; break;
	case 30: nNodes = 35; break;
	case 31: nNodes = 56; break;
	default:
	  ERROR("unknown Gmsh element type");
	}
      
      return nNodes; 
    }
     
   
  } //end of namespace Gmsh
     
} //end of namespace Utilities
     
