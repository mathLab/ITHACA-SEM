#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <boost/lexical_cast.hpp>

#include<LibUtilities/BasicUtils/ErrorUtil.hpp>

using namespace std;

void  PrintConditions(ofstream& output);

int main(int argc, char *argv[])
{
    vector<double> xc,yc; 
    int            nx = 0;
    int            ny = 0;
    int            nummodes = 7;
    int            type = 3;
    int            i,j,k;
    string         output_file;
    ofstream       output;

    if(argc != 6)
    {
        //                     argv: 1  2   3    4        5
        cerr << "Usage "<<argv[0]<<" nx ny type nummodes output_file.xml\n";
        cerr << "where nx is the number of points in x direction,\n";
        cerr << "      ny is the number of points in y direction,\n";
        cerr << "      nummodes is the number of boundary modes.\n";
        cerr << "It generates regular mesh of triangles (type = 2) or quadrilaterals (type = 3).\n";
        cerr << "All vertices are evenly distributed within unit square\n";
        exit(1);
    }

    try{
        nx          = boost::lexical_cast<int>(argv[1]);
        ny          = boost::lexical_cast<int>(argv[2]);
        type        = boost::lexical_cast<int>(argv[3]);
        nummodes    = boost::lexical_cast<int>(argv[4]);
        output_file = boost::lexical_cast<string>(argv[5]);

        if ((type != 2) && (type != 3))
        {
            cerr << "Wrong mesh type";
            throw 1;
        }

        output.open(output_file.c_str());

        for (i = 0; i < nx; i++)
        {
            xc.push_back( (double)i * (1.0 / (nx-1)) );
        }
        for (i = 0; i < ny; i++)
        {
            yc.push_back( (double)i * (1.0 / (ny-1)) );
        }

        output<< "<?xml version=\"1.0\" encoding=\"utf-8\" ?>" << endl;
        output<< "<NEKTAR>" << endl;


        output << "<EXPANSIONS>" << endl;
        output << "<E COMPOSITE=\"C[0]\" NUMMODES=\"" << nummodes << "\" FIELDS=\"u\" TYPE=\"MODIFIED\" />" <<endl;
        output << "</EXPANSIONS>\n" << endl;

        PrintConditions(output);

        //Vertices 
        output << "<GEOMETRY DIM=\"2\" SPACE=\"2\">" << endl;
        output << "  <VERTEX>" << endl;

        nx = xc.size();
        ny = yc.size();
        for(j = 0; j < ny; ++j) 
        {
            for(k = 0; k < nx; ++k)
            {
                output << "    <V ID=\"" << j*nx+k << "\">\t";
                output << xc[k] << " "   << yc[j] << " 0.0";
                output << "  </V>" << endl; 
            }
        }

        output << "  </VERTEX>\n" << endl;

        // Edges. By default it generates edges for
        // quadrilateral mesh.
        output << "  <EDGE>" << endl;
        int cnt = 0;
        for(j = 0; j < ny-1; ++j)
        {
            for(i = 0; i < nx-1; ++i)
            {
                output << "    <E ID=\"" << cnt++ << "\">\t";
                output << j*nx+i <<"  " <<  j*nx + i+1 ;
                output << "  </E>" << endl; 
            }

            for(i = 0; i < nx; ++i)
            {
                output << "    <E ID=\"" << cnt++ << "\">\t";
                output << j*nx+i <<"  " << (j+1)*nx+i;
                output << "  </E>" << endl; 
            }
        }

        for(i = 0; i < nx-1; ++i)
        {
            output << "    <E ID=\"" << cnt++ << "\">\t";
            output << j*nx+i <<"  " <<  j*nx + i+1 ;
            output << "  </E>" << endl; 
        }

        // total number of quad edges. Useful in renumbering
        // diagonal edges later on.
        int total_quad_edges = cnt-1;

        // Triangular mesh is made by adding diagonal segments
        if (type == 2)
        {
            // generating diagonal edges

            for(j = 0; j < ny-1; ++j)
            {
                for(i = 0; i < nx-1; ++i)
                {
                    output << "    <E ID=\"" << cnt++ << "\">\t";
                    output << j*nx+i <<"  " <<  (j+1)*nx + i+1 ;
                    output << "  </E>" << endl; 
                }
            }
        }
        output << "  </EDGE>\n" << endl;

        output << "  <ELEMENT>" << endl;
        cnt = 0;
        switch(type)
        {
        // quadrilaterals
        case 3:
            for(j = 0; j < ny-1; ++j)
            {
                for(i = 0; i < nx-1; ++i)
                {
                    output << "    <Q ID=\"" << cnt++ << "\">\t";
                    output << j*((nx-1)+nx)+i <<"  " <<  j*((nx-1)+nx) + (nx-1)+i+1 << "   " ;
                    output << (j+1)*((nx-1)+nx)+i <<"  " <<  j*((nx-1)+nx) + (nx-1)+i ;
                    output << "  </Q>" << endl; 
                }
            }
            break;
        // triangles
        case 2:
            for(j = 0; j < ny-1; ++j)
            {
                for(i = 0; i < nx-1; ++i)
                {
                    output << "    <T ID=\"" << cnt++ << "\">\t";
                    output << total_quad_edges + 1 + (j*(nx-1)+i) << "   " ;
                    output << (j+1)*((nx-1)+nx)+i <<"  " <<  j*((nx-1)+nx) + (nx-1)+i ;
                    output << "  </T>" << endl; 

                    output << "    <T ID=\"" << cnt++ << "\">\t";
                    output << j*((nx-1)+nx)+i <<"  " <<  j*((nx-1)+nx) + (nx-1)+i+1 << "   " ;
                    output << total_quad_edges + 1 + (j*(nx-1)+i);
                    output << "  </T>" << endl; 
                }
            }
            break;
        default:
            cerr << "unknown element type " << type << "\n";
            throw 1;
        }
        output << "  </ELEMENT>\n" << endl;


        output << "<COMPOSITE>" << endl;

        switch(type)
        {
        // quadrilaterals
        case 3:
            output << "<C ID=\"0\"> Q[0-" << (nx-1)*(ny-1)-1 << "] </C>" << endl;
            break;
        // triangles
        case 2:
            output << "<C ID=\"0\"> T[0-" << 2*(nx-1)*(ny-1)-1 << "] </C>" << endl;
            break;
        default:
            cerr << "unknown element type\n";
            throw 1;
        }

        // boundary composites coincide for both mesh element types
        output << "<C ID=\"1\"> E[";
        for(i = 0; i < nx-1; ++i)
        {
            output << i;
            if(i !=  nx-2) 
            {
                output << ",";
            }
        }
        output << "] </C>   // south border" << endl;

        output << "<C ID=\"2\"> E[";
        for(i = 0; i < ny-1; ++i)
        {
            output << (nx-1)*(i+1) + nx*i;
            if(i != ny-2)
            {
                output << ",";
            }
        }
        output << "] </C>   // west border" << endl;

        output << "<C ID=\"3\"> E[";
        for(i = 0; i < nx-1; ++i)
        {
            output << (nx-1)*(ny-1)+ nx*(ny-1)+ i;
            if(i != nx-2)
            {
                output << ",";
            }
        }
        output << "] </C>   // north border" << endl;

        output << "<C ID=\"4\"> E[";
        for(i = 0; i < ny-1; ++i)
        {
            output << (nx-1)*(i+1) + nx*i + nx-1;
            if(i != ny-2)
            {
                output << ",";
            }
        }
        output << "] </C>   // East border" << endl;


        output << "</COMPOSITE>\n" << endl;

             
        output << "<DOMAIN> C[0] </DOMAIN>\n" << endl;
        output << "</GEOMETRY>\n" << endl;

        output << "</NEKTAR>" << endl;

    }
    catch(...)
    {
        return 1;
    }

    return 0;

}



void  PrintConditions(ofstream& output)
{
    output << "<CONDITIONS>" << endl;

    output << "<SOLVERINFO>" << endl;
    output << "<I PROPERTY=\"GlobalSysSoln\"        VALUE=\"DirectFull\"/>" << endl;
    output << "</SOLVERINFO>\n" << endl;

    output << "<PARAMETERS>" << endl;
    output << "<P> TimeStep      = 0.002  </P>" << endl;
    output << "<P> Lambda        = 1      </P>" << endl;
    output << "</PARAMETERS>\n" << endl;
    
    output << "<VARIABLES>" << endl;
    output << "  <V ID=\"0\"> u </V>" << endl; 
    output << "</VARIABLES>\n" << endl;
    
    output << "<BOUNDARYREGIONS>" << endl;
    output << "<B ID=\"0\"> C[1] </B>" << endl;
    output << "<B ID=\"1\"> C[2] </B>" << endl;
    output << "<B ID=\"2\"> C[3] </B>" << endl;
    output << "<B ID=\"3\"> C[4] </B>" << endl;
    output << "</BOUNDARYREGIONS>\n" << endl;
    
    output << "<BOUNDARYCONDITIONS>" << endl;
    output << "  <REGION REF=\"0\"> // South border " << endl;
    output << "     <D VAR=\"u\" VALUE=\"sin(PI/2*x)*sin(PI/2*y)\" />"  << endl;
    output << "  </REGION>" << endl;
            
    output << "  <REGION REF=\"1\"> // West border " << endl;
    output << "     <D VAR=\"u\" VALUE=\"sin(PI/2*x)*sin(PI/2*y)\" />"  << endl;
    output << "  </REGION>" << endl;
    
    output << "  <REGION REF=\"2\"> // North border " << endl;
    output << "     <D VAR=\"u\" VALUE=\"sin(PI/2*x)*sin(PI/2*y)\" />"  << endl;
    output << "  </REGION>" << endl;
    
    output << "  <REGION REF=\"3\"> // East border " << endl;
    output << "     <D VAR=\"u\" VALUE=\"sin(PI/2*x)*sin(PI/2*y)\" />"  << endl;
    output << "  </REGION>" << endl;
    output << "</BOUNDARYCONDITIONS>" << endl;

    output << "  <FUNCTION NAME=\"Forcing\">" << endl;
    output << "     <E VAR=\"u\" VALUE=\"-(Lambda + 2*PI*PI/4)*sin(PI/2*x)*sin(PI/2*y)\" />" << endl;
    output << "  </FUNCTION>" << endl;

    output << "  <FUNCTION NAME=\"ExactSolution\">" << endl;
    output << "     <E VAR=\"u\" VALUE=\"sin(PI/2*x)*sin(PI/2*y)\" />" << endl;
    output << "  </FUNCTION>" << endl;

    output << "</CONDITIONS>" << endl;
}

