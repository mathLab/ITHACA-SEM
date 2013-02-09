#include <cstdio>
#include <cstdlib>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>  // for ASSERTL0
#include <SpatialDomains/MeshGraph.h>   // for FieldDefinitions, etc

using namespace Nektar;

int main(int argc, char *argv[])
{
    NekDouble scal1,scal2;

    if(argc != 6)
    {
        fprintf(stderr,"Usage: FldAddFld scal1 scal2 fieldfile1 fieldfile2 outfield\n"
                "\t produces scal1*fieldfiel1 + scal2*fieldfile2 in outfield\n" );
        exit(1);
    }

    scal1  = boost::lexical_cast<double>(argv[argc-5]);
    scal2  = boost::lexical_cast<double>(argv[argc-4]);

    //default meshgraph
    SpatialDomains::MeshGraph graph; 

    //----------------------------------------------
    // Import fieldfile1.
    string fieldfile1(argv[argc-3]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef1;
    vector<vector<NekDouble> > fielddata1;
    graph.Import(fieldfile1,fielddef1,fielddata1);
    //----------------------------------------------

    //----------------------------------------------
    // Import fieldfile2.
    string fieldfile2(argv[argc-2]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef2;
    vector<vector<NekDouble> > fielddata2;
    graph.Import(fieldfile2,fielddef2,fielddata2);
    //----------------------------------------------


    ASSERTL0(fielddata1.size() == fielddata2.size(),"Inner has different size");

    //----------------------------------------------
    // Add fielddata2 to fielddata1 using m_fields definition to align data. 
    

    for(int i = 0; i < fielddata1.size(); ++i)
    {
        int j;
        int datalen1 = fielddata1[i].size()/fielddef1[i]->m_fields.size();
        int datalen2 = fielddata2[i].size()/fielddef2[i]->m_fields.size();

        ASSERTL0(datalen1 == datalen2,"Data per field is of different length");
        
        for(int k = 0; k < fielddef1[i]->m_fields.size(); ++k)
        {
            int offset = 0;
            for(j = 0; j < fielddef2[i]->m_fields.size(); ++j)
            {
                if(fielddef1[i]->m_fields[k] == fielddef2[i]->m_fields[j])
                {
                    break;
                }
                offset  += datalen1;
            }
            
            if(j == fielddef2[i]->m_fields.size())
            {
                for(j = 0; j < datalen1; ++j)
                {
                    fielddata1[i][datalen1*k+j] *= scal1;
                }
            }
            else // add fields 
            {
                for(j = 0; j < datalen1; ++j)
                {
                    fielddata1[i][datalen1*k+j] *= scal1;
                    fielddata1[i][datalen1*k+j] += scal2*fielddata2[i][offset + j];
                }
            }

        }

        // now check to see if any field in fielddef2[i]->m_fields is
        // not defined in fielddef1[i]->m_fields
        for(int k = 0; k < fielddef2[i]->m_fields.size(); ++k)
        {
            for(j = 0; j < fielddef1[i]->m_fields.size(); ++j)
            {
                if(fielddef2[i]->m_fields[k] == fielddef1[i]->m_fields[j])
                {
                    break;
                }
            }
            
            if(j == fielddef1[i]->m_fields.size())
            {
                for(j = 0; j < datalen2; ++j)
                {
                    fielddata2[i][datalen2*k+j] *= scal2;
                }

                // add this field to fielddata1
                fielddef1[i]->m_fields.push_back(fielddef2[i]->m_fields[k]);
                fielddata1[i].insert(fielddata1[i].end(),&(fielddata2[i][k*datalen2]), 
                                     &(fielddata2[i][k*datalen2])+datalen1);
            }

        }
        
    }
    //----------------------------------------------

    //-----------------------------------------------
    // Write out datafile. 
    graph.Write(argv[argc-1], fielddef1, fielddata1);
    //-----------------------------------------------

    return 0;
}

