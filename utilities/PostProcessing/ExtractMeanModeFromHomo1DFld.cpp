#include <cstdio>
#include <cstdlib>
#include <SpatialDomains/MeshGraph.h>   // for FieldDefinitions, etc
#include <StdRegions/StdTriExp.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    if(argc != 3)
    {
        fprintf(stderr,"Usage: ExtractmeanModeFromHomo1DFld fieldfile outfield\n");
        exit(1);
    }
    
    //default meshgraph
    SpatialDomains::MeshGraph graph; 

    //----------------------------------------------
    // Import fieldfile.
    string fieldfile(argv[argc-2]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    LibUtilities::Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------


    vector<vector<NekDouble> > combineddata;
    vector<LibUtilities::FieldDefinitionsSharedPtr> newfielddef;

    //----------------------------------------------
    // put mean data consecutively 
    for(int i = 0; i < fielddata.size(); ++i)
    {
        ASSERTL0(fielddef[i]->m_numHomogeneousDir == 1,"Expected fieldfile to have one homogeneous direction");

        if(fielddef[i]->m_homogeneousZIDs[0] != 0)
        {
            continue;
        }
        else
        {
            int nz = fielddef[i]->m_homogeneousZIDs.size();

            fielddef[i]->m_numHomogeneousDir = 0;
            fielddef[i]->m_basis.resize(2);
            newfielddef.push_back(fielddef[i]);
            
        
            // Determine the number of coefficients per element
            int ncoeffs;
            switch(fielddef[i]->m_shapeType)
            {
            case  LibUtilities::eTriangle:
                ncoeffs = LibUtilities::StdTriData::getNumberOfCoefficients(fielddef[i]->m_numModes[0], fielddef[i]->m_numModes[1]);
                break;
            case LibUtilities::eQuadrilateral:
                ncoeffs = fielddef[i]->m_numModes[0]*fielddef[i]->m_numModes[1];
                break;
            default:
                ASSERTL0(false,"Shape not recognised");
                break;
            }
            
            std::vector<NekDouble>::iterator vec_iter; 
            
            vector<NekDouble> newdata;
            vec_iter = fielddata[i].begin();
                        
            for(int k = 0; k < fielddef[i]->m_fields.size(); ++k)
            {
                // copy data from each field into consecutive order 
                for(int n = 0; n < fielddef[i]->m_elementIDs.size(); ++n)
                {
                    // put zero mode into newdata
                    newdata.insert(newdata.end(),vec_iter, vec_iter+ncoeffs);
                    vec_iter += nz*ncoeffs; 
                }
             }
            combineddata.push_back(newdata);
        }
    }
    //----------------------------------------------

    //-----------------------------------------------
    // Write out datafile. 
    LibUtilities::Write(argv[argc-1], newfielddef, combineddata);
    //-----------------------------------------------

    return 0;
}

