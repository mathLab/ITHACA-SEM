#ifndef H_GRAPH
#define H_GRAPH

#include <list>
#include <utility>

#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>  
#include <boost/graph/adjacency_list.hpp>
#include <boost/shared_ptr.hpp>

namespace LibUtilities{
  namespace NekGraph{
    
  
    // ---------------------------------------------------------------------------
    /// GraphVertexObject

    
    class GraphVertexObject{
    public:
      GraphVertexObject(){};      
      virtual ~GraphVertexObject(){};
          
    };//end GraphVertexObject
    

    // ---------------------------------------------------------------------------
    /// GraphEdgeObject


    class GraphEdgeObject: public std::pair<GraphVertexObject*,GraphVertexObject*>{
    public:
      GraphEdgeObject(GraphVertexObject *v1, GraphVertexObject *v2){
	first = v1;
	second = v2;
      }
      
      virtual ~GraphEdgeObject(){
	first = (GraphVertexObject*)NULL;
	second = (GraphVertexObject*)NULL;
      }
      
    protected:
      
    };//end GraphEdgeObject

    
    // ---------------------------------------------------------------------------
    /// NekGraph


    typedef boost::adjacency_list<boost::listS, boost::listS, 
				  boost::bidirectionalS,boost::property<boost::vertex_index_t,int>, 
				  boost::property<boost::edge_weight_t,double> > Graph;
    
    class NekGraph{
    public:
      NekGraph(){
      
      }

      virtual ~NekGraph(){

      }


      int AddVertex(boost::shared_ptr<GraphVertexObject> ptr){

      }
    
  
      void InsertEdge(GraphEdgeObject &eobj){

      }


      
    protected:  
      Graph _graph;
      
    };//end NekGraph
    
    
  }//end namespace
}//end namespace
#endif
