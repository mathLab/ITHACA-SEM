#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int     i, j, nq,  coordim;
    Array<OneD,NekDouble>  fce;
    Array<OneD,NekDouble>  xc0,xc1,xc2;

    if(argc != 3)
    {
        fprintf(stderr,"Usage: FldToVtk  meshfile fieldfile\n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    SpatialDomains::MeshGraph graph;
    SpatialDomains::MeshGraphSharedPtr graphShPt = graph.Read(meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // Import field file.
    string fieldfile(argv[argc-1]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graphShPt->Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------

    //----------------------------------------------
    // Set up Expansion information
    vector< vector<LibUtilities::PointsType> > pointstype;
    for(i = 0; i < fielddef.size(); ++i)
    {
        vector<LibUtilities::PointsType> ptype;
        for(j = 0; j < 2; ++j)
        {
            ptype.push_back(LibUtilities::ePolyEvenlySpaced);
        }
        pointstype.push_back(ptype);
    }
    graphShPt->SetExpansions(fielddef,pointstype);
    //----------------------------------------------


    //----------------------------------------------
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_Fields.size();
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields);

    switch(expdim)
    {
    case 1:
        {
            SpatialDomains::MeshGraph1DSharedPtr mesh;

            if(!(mesh = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph1D>(graphShPt)))
            {
                ASSERTL0(false,"Dynamics cast failed");
            }

            MultiRegions::ExpList1DSharedPtr Exp1D;
            Exp1D = MemoryManager<MultiRegions::ExpList1D>
                                                    ::AllocateSharedPtr(*mesh);
            Exp[0] = Exp1D;
            for(i = 1; i < nfields; ++i)
            {
                Exp[i] = MemoryManager<MultiRegions::ExpList1D>
                                                    ::AllocateSharedPtr(*Exp1D);
            }
        }
        break;
    case 2:
        {
            SpatialDomains::MeshGraph2DSharedPtr mesh;

            if(!(mesh = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph2D>(graphShPt)))
            {
                ASSERTL0(false,"Dynamics cast failed");
            }

            MultiRegions::ExpList2DSharedPtr Exp2D;
            Exp2D = MemoryManager<MultiRegions::ExpList2D>
                                                    ::AllocateSharedPtr(*mesh);
            Exp[0] =  Exp2D;

            for(i = 1; i < nfields; ++i)
            {
                Exp[i] = MemoryManager<MultiRegions::ExpList2D>
                                                    ::AllocateSharedPtr(*Exp2D);
            }
        }
        break;
    case 3:
        ASSERTL0(false,"3D not set up");
        break;
    default:
        ASSERTL0(false,"Expansion dimension not recognised");
        break;
    }
    //----------------------------------------------

    //----------------------------------------------
    // Copy data from field file
    for(j = 0; j < nfields; ++j)
    {
        for(int i = 0; i < fielddata.size(); ++i)
        {
            Exp[j]->ExtractDataToCoeffs(fielddef[i],
                                        fielddata[i],
                                        fielddef[i]->m_Fields[j]);
        }
        Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
    }
    //----------------------------------------------

    //----------------------------------------------
    // Write solution
//    std::string var = "";

//    for(int j = 0; j < Exp.num_elements(); ++j)
//    {
//        var = var + ", " + fielddef[0]->m_Fields[j];
//    }

    string   outname(strtok(argv[argc-1],"."));
    outname += ".vtu";
    ofstream outfile(outname.c_str());
    cout << "Writing file: " << outname << " ... ";

    Exp[0]->WriteVtkHeader(outfile);
    for(int i = 0; i < Exp[0]->GetExpSize(); ++i)
    {
        Exp[0]->WriteVtkPieceHeader(outfile,i);
        for(int j = 0; j < Exp.num_elements(); ++j)
        {
            Exp[j]->WriteVtkPieceData(outfile,i, fielddef[0]->m_Fields[j]);
        }
        Exp[0]->WriteVtkPieceFooter(outfile,i);
    }
    Exp[0]->WriteVtkFooter(outfile);
    cout << "Done " << endl;
    //----------------------------------------------
    return 0;
}

