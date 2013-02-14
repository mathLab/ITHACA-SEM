#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

#include<LibUtilities/BasicUtils/ErrorUtil.hpp>

// Use the stl version, primarily for string.
#ifndef TIXML_USE_STL
#define TIXML_USE_STL
#endif

#include <tinyxml/tinyxml.h>


void Header(FILE *, int nel);
void Middle(FILE *);
void End(FILE *);

using namespace std;

void  PrintConditions(void);

int main(int argc, char *argv[])
{
    vector<double> xc,yc; 
    int            nx = 0, ny = 0;
    int            i,j,k;
    
    if(argc != 2)
    {
        fprintf(stderr,"Usage RectanuglarMesh file\n");
        exit(1);
    }
    
    try{
        
        TiXmlDocument doc(argv[argc-1]);
        bool loadOkay = doc.LoadFile();
        
        std::stringstream errstr;
        errstr << "Unable to load file: " << argv[argc-1] << " (";
        errstr << doc.ErrorDesc() << ", line " << doc.ErrorRow()
               << ", column " << doc.ErrorCol() << ")";
        ASSERTL0(loadOkay, errstr.str());
        
        TiXmlHandle docHandle(&doc);
        TiXmlElement* master = NULL;   
        TiXmlElement* block = NULL;
        
        master = doc.FirstChildElement("NEKBLOCK");
        ASSERTL0(master, "Unable to find NEKBLOCK tag in file.");
        
        // Find the Mesh tag and same the dim and space attributes
        block = master->FirstChildElement("XBLOCK");
        
        ASSERTL0(block, "Unable to find XBLOCK tag in file.");
        TiXmlElement *val = block->FirstChildElement("X");    
        while (val)
        {
            TiXmlNode *xval = val->FirstChild();
            
            std::istringstream valDataStrm(xval->ToText()->Value());
            
            try
            {
                while(!valDataStrm.fail())
                {
                    double x_val;
                    valDataStrm >> x_val;
                    
                    if (!valDataStrm.fail())
                    {
                        xc.push_back(x_val);
                    }
                }
            }
            catch(...)
            {
                ASSERTL0(false, "Unable to read Xval data.");
            }
            
            val= val->NextSiblingElement("X");
        }
        
        block = master->FirstChildElement("YBLOCK");
        
        val = block->FirstChildElement("Y");
        while (val)
        {
            TiXmlNode *yval = val->FirstChild();
            
            std::istringstream valDataStrm(yval->ToText()->Value());
            
            try
            {
                while(!valDataStrm.fail())
                {
                    double y_val;
                    valDataStrm >> y_val;
                    
                    if (!valDataStrm.fail())
                    {
                        yc.push_back(y_val);
                    }
                }
            }
            catch(...)
            {
                ASSERTL0(false, "Unable to read Yval data.");
            }
            
            val= val->NextSiblingElement("Y");
        }
        
        cout << "<?xml version=\"1.0\" encoding=\"utf-8\" ?>" << endl;
        cout << "<NEKTAR>" << endl;
 

        cout << "<EXPANSIONS>" << endl;
        cout << "<E COMPOSITE=\"C[0]\" NUMMODES=\"7\" FIELDS=\"u\" TYPE=\"MODIFIED\" />" <<endl;
        cout << "</EXPANSIONS>\n" << endl;

        PrintConditions();

        //Vertices 
        cout << "<GEOMETRY DIM=\"2\" SPACE=\"2\">" << endl;
        cout << "  <VERTEX>" << endl;
        
        nx = xc.size();
        ny = yc.size();
        for(j = 0; j < ny; ++j) 
        {
            for(k = 0; k < nx; ++k)
            {
                cout << "    <V ID=\"" << j*nx+k << "\">\t";
                cout << std::setprecision(8)<< xc[k] << " "   << yc[j] << " 0.0";
                cout << "  </V>" << endl; 
            }
        }
    
        cout << "  </VERTEX>\n" << endl;

        cout << "  <EDGE>" << endl;
        int cnt = 0;
        for(j = 0; j < ny-1; ++j)
        {
            for(i = 0; i < nx-1; ++i)
            {
                cout << "    <E ID=\"" << cnt++ << "\">\t";
                cout << j*nx+i <<"  " <<  j*nx + i+1 ;
                cout << "  </E>" << endl; 
            }

            for(i = 0; i < nx; ++i)
            {
                cout << "    <E ID=\"" << cnt++ << "\">\t";
                cout << j*nx+i <<"  " << (j+1)*nx+i;
                cout << "  </E>" << endl; 
            }
        }
        
        for(i = 0; i < nx-1; ++i)
        {
            cout << "    <E ID=\"" << cnt++ << "\">\t";
            cout << j*nx+i <<"  " <<  j*nx + i+1 ;
            cout << "  </E>" << endl; 
        }
        cout << "  </EDGE>\n" << endl;



        cout << "  <ELEMENT>" << endl;
        cnt = 0;
        for(j = 0; j < ny-1; ++j)
        {
            for(i = 0; i < nx-1; ++i)
            {
                cout << "    <Q ID=\"" << cnt++ << "\">\t";
                cout << j*((nx-1)+nx)+i <<"  " <<  j*((nx-1)+nx) + (nx-1)+i+1 << "   " ;
                cout << (j+1)*((nx-1)+nx)+i <<"  " <<  j*((nx-1)+nx) + (nx-1)+i ;
                cout << "  </Q>" << endl; 
            }
        }
        cout << "  </ELEMENT>\n" << endl;


        cout << "<COMPOSITE>" << endl;
        cout << "<C ID=\"0\"> Q[0-" << (nx-1)*(ny-1)-1 << "] </C>" << endl;

        cout << "<C ID=\"1\"> E[";
        for(i = 0; i < nx-1; ++i)
        {
            cout << i;
            if(i !=  nx-2) 
            {
                cout << ",";
            }
        }
        cout << "] </C>   // south border" << endl;

        cout << "<C ID=\"2\"> E[";
        for(i = 0; i < ny-1; ++i)
        {
            cout << (nx-1)*(i+1) + nx*i;
            if(i != ny-2)
            {
                cout << ",";
            }
        }
        cout << "] </C>   // west border" << endl;

        cout << "<C ID=\"3\"> E[";
        for(i = 0; i < nx-1; ++i)
        {
            cout << (nx-1)*(ny-1)+ nx*(ny-1)+ i;
            if(i != nx-2)
            {
                cout << ",";
            }
        }
        cout << "] </C>   // north border" << endl;

        cout << "<C ID=\"4\"> E[";
        for(i = 0; i < ny-1; ++i)
        {
            cout << (nx-1)*(i+1) + nx*i + nx-1;
            if(i != ny-2)
            {
                cout << ",";
            }
        }
        cout << "] </C>   // East border" << endl;


        cout << "</COMPOSITE>\n" << endl;

             
        cout << "<DOMAIN> C[0] </DOMAIN>\n" << endl;
        cout << "</GEOMETRY>\n" << endl;

        cout << "</NEKTAR>" << endl;

    }
    catch(...)
    {
        return 1;
    }

    return 0;

}



void  PrintConditions(void)
{
    cout << "<CONDITIONS>" << endl;
    
    cout << "<SOLVERINFO>" << endl;
    cout << "<I PROPERTY=\"SolverType\"        VALUE=\" \"/>" << endl;
    cout << "</SOLVERINFO>\n" << endl;
            
    cout << "<PARAMETERS>" << endl;
    cout << "<P> TimeStep      = 0.002  </P>" << endl;
    cout << "</PARAMETERS>\n" << endl;
    
    cout << "<VARIABLES>" << endl;
    cout << "  <V ID=\"0\"> u </V>" << endl; 
    cout << "</VARIABLES>\n" << endl;
    
    cout << "<BOUNDARYREGIONS>" << endl;
    cout << "<B ID=\"0\"> C[1] </B>" << endl;
    cout << "<B ID=\"1\"> C[2] </B>" << endl;
    cout << "<B ID=\"2\"> C[3] </B>" << endl;
    cout << "<B ID=\"3\"> C[4] </B>" << endl;
    cout << "</BOUNDARYREGIONS>\n" << endl;
    
    cout << "<BOUNDARYCONDITIONS>" << endl;
    cout << "  <REGION REF=\"0\"> // South border " << endl;
    cout << "     <D VAR=\"u\" VALUE=\"0\" />"  << endl;
    cout << "  </REGION>" << endl;
            
    cout << "  <REGION REF=\"1\"> // West border " << endl;
    cout << "     <D VAR=\"u\" VALUE=\"0\" />"  << endl;
    cout << "  </REGION>" << endl;
    
    cout << "  <REGION REF=\"2\"> // North border " << endl;
    cout << "     <D VAR=\"u\" VALUE=\"0\" />"  << endl;
    cout << "  </REGION>" << endl;
    
    cout << "  <REGION REF=\"3\"> // East border " << endl;
    cout << "     <D VAR=\"u\" VALUE=\"0\" />"  << endl;
    cout << "  </REGION>" << endl;
    cout << "</BOUNDARYCONDITIONS>" << endl;

    cout << "</CONDITIONS>" << endl;
}

