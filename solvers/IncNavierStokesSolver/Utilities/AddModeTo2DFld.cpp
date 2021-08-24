#include <cstdio>
#include <cstdlib>
#include <SpatialDomains/MeshGraph.h>   // for FieldDefinitions, etc
#include <StdRegions/StdTriExp.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    NekDouble scal1,scal2;

    if(argc != 6)
    {
        fprintf(stderr,"Usage: AddModeTo2DFld scal1 scal2 2Dfieldfile1 fieldfile2 outfield\n"
                "\t produces scal1*2Dfieldfiel1 + scal2*fieldfile2 in outfield\n" );
        exit(1);
    }

    scal1  = boost::lexical_cast<double>(argv[argc-5]);
    scal2  = boost::lexical_cast<double>(argv[argc-4]);

    //----------------------------------------------
    // Import fieldfile1.
    string fieldfile1(argv[argc-3]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef1;
    vector<vector<NekDouble> > fielddata1;
    LibUtilities::Import(fieldfile1,fielddef1,fielddata1);
    //----------------------------------------------

    //----------------------------------------------
    // Import fieldfile2.
    string fieldfile2(argv[argc-2]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef2;
    vector<vector<NekDouble> > fielddata2;
    LibUtilities::Import(fieldfile2,fielddef2,fielddata2);
    //----------------------------------------------

    vector<vector<NekDouble> > combineddata;

    ASSERTL0(fielddata1.size() == fielddata2.size(),"Inner has different size");
    //----------------------------------------------
    // Add fielddata2 to fielddata1 using m_fields definition to align data.

    int i = 0;
    int j = 0;
    int k = 0;
    int n = 0;

    for(i = 0; i < fielddata2.size(); ++i)
    {
        ASSERTL0(fielddef2[i]->m_numHomogeneousDir == 1,"Expected second fld to have one homogeneous direction");
        ASSERTL0(fielddef2[i]->m_numModes[2] == 2,"Expected Fourier field to have 2 modes");

        int datalen1 = fielddata1[i].size()/fielddef1[i]->m_fields.size();
        int datalen2 = fielddata2[i].size()/fielddef2[i]->m_fields.size();

        ASSERTL0(datalen1*2 == datalen2,"Data per fields is note compatible");

        // Determine the number of coefficients per element
        int ncoeffs = 0;
        switch(fielddef2[i]->m_shapeType)
        {
        case LibUtilities::eTriangle:
            ncoeffs = LibUtilities::StdTriData::getNumberOfCoefficients(fielddef2[i]->m_numModes[0], fielddef2[i]->m_numModes[1]);
            break;
        case LibUtilities::eQuadrilateral:
            ncoeffs = fielddef2[i]->m_numModes[0]*fielddef2[i]->m_numModes[1];
            break;
        default:
                ASSERTL0(false,"Shape not recognised");
            break;
        }

        // array for zero packing
        Array<OneD,NekDouble> Zero(ncoeffs,0.0);

        // scale first and second fields
        Vmath::Smul(fielddata1[i].size(), scal1, &fielddata1[i][0], 1,
                                                 &fielddata1[i][0], 1);
        Vmath::Smul(fielddata2[i].size(), scal2, &fielddata2[i][0], 1,
                                                 &fielddata2[i][0], 1);

        vector<NekDouble> newdata;
        auto vec_iter = fielddata2[i].begin();

        for(k = 0; k < fielddef2[i]->m_fields.size(); ++k)
        {
            // get location of 2D field information in order of field2 ordering
            int offset = 0;
            for(j = 0; j < fielddef1[i]->m_fields.size(); ++j)
            {
                if(fielddef1[i]->m_fields[j] == fielddef2[i]->m_fields[k])
                {
                    break;
                }
                offset  += datalen1;
            }

            if(j != fielddef1[i]->m_fields.size())
            {
                for(n = 0; n < fielddef2[i]->m_elementIDs.size(); ++n)
                {
                    // Real zero component
                    newdata.insert(newdata.end(),
                                   &(fielddata1[i][offset+n*ncoeffs]),
                                   &(fielddata1[i][offset+n*ncoeffs])
                                        + ncoeffs);

                    // Imaginary zero component;
                    newdata.insert(newdata.end(),&Zero[0],&Zero[0] + ncoeffs);

                    // Put orginal mode in here.
                    newdata.insert(newdata.end(),vec_iter, vec_iter+2*ncoeffs);
                    vec_iter += 2*ncoeffs;
                }
            }
            else
            {

                for(n = 0; n < fielddef2[i]->m_elementIDs.size(); ++n)
                {
                    // Real & Imag zero component
                    newdata.insert(newdata.end(),&Zero[0],&Zero[0] + ncoeffs);
                    newdata.insert(newdata.end(),&Zero[0],&Zero[0] + ncoeffs);

                    // Put orginal mode in here.
                    newdata.insert(newdata.end(),vec_iter, vec_iter+2*ncoeffs);
                    vec_iter += 2*ncoeffs;
                }
            }
        }
        combineddata.push_back(newdata);
        fielddef2[i]->m_numModes[2] += 2;
        fielddef2[i]->m_homogeneousZIDs.push_back(2);
        fielddef2[i]->m_homogeneousZIDs.push_back(3);

        // check to see if any field in fielddef1[i]->m_fields is
        // not defined in fielddef2[i]->m_fields
        for(k = 0; k < fielddef1[i]->m_fields.size(); ++k)
        {
            for(j = 0; j < fielddef2[i]->m_fields.size(); ++j)
            {
                if(fielddef1[i]->m_fields[k] == fielddef2[i]->m_fields[j])
                {
                    break;
                }
            }

            if(j == fielddef2[i]->m_fields.size())
            {
                cout << "Warning: Field \'" << fielddef1[i]->m_fields[k]
                     << "\' was not included in output " << endl;
            }

        }

    }
    //----------------------------------------------

    //-----------------------------------------------
    // Write out datafile.
    LibUtilities::Write(argv[argc-1], fielddef2, combineddata);
    //-----------------------------------------------

    return 0;
}

