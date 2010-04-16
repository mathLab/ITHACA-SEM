#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int i,j;

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
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_Fields.size();
    int addfields = 1;
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields + addfields);

    switch(expdim)
    {
    case 1:
        {
            SpatialDomains::MeshGraph1DSharedPtr mesh;

            if(!(mesh = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph1D>(graphShPt)))
            {
                ASSERTL0(false,"Dynamic cast failed");
            }

            MultiRegions::ExpList1DSharedPtr Exp1D;
            Exp1D = MemoryManager<MultiRegions::ExpList1D>
                                                    ::AllocateSharedPtr(*mesh);
            Exp[0] = Exp1D;
            for(i = 1; i < nfields + addfields; ++i)
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
                ASSERTL0(false,"Dynamic cast failed");
            }

            MultiRegions::ExpList2DSharedPtr Exp2D;
            Exp2D = MemoryManager<MultiRegions::ExpList2D>
                                                    ::AllocateSharedPtr(*mesh);
            Exp[0] =  Exp2D;

            for(i = 1; i < nfields + addfields; ++i)
            {
                Exp[i] = MemoryManager<MultiRegions::ExpList2D>
                                                    ::AllocateSharedPtr(*Exp2D);
            }
        }
        break;
    case 3:
        {
            SpatialDomains::MeshGraph3DSharedPtr mesh;

            if(!(mesh = boost::dynamic_pointer_cast<
                                    SpatialDomains::MeshGraph3D>(graphShPt)))
            {
                ASSERTL0(false,"Dynamic cast failed");
            }

            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = MemoryManager<MultiRegions::ExpList3D>
                                                    ::AllocateSharedPtr(*mesh);
            Exp[0] =  Exp3D;

            for(i = 1; i < nfields + addfields; ++i)
            {
                Exp[i] = MemoryManager<MultiRegions::ExpList3D>
                                                    ::AllocateSharedPtr(*Exp3D);
            }
        }
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
            Exp[j]->ExtractDataToCoeffs(fielddef [i],
                                        fielddata[i],
                                        fielddef [i]->m_Fields[j]);
        }
        Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
    }
    //----------------------------------------------

    //----------------------------------------------
    // Compute gradients of fields and compute reentricity
    ASSERTL0(nfields == 2, "Need two fields (u,v) to add reentricity");
    int nq = Exp[0]->GetNpoints();
    Array<OneD, NekDouble> grad_u[2], grad_v[2];
    Array<OneD, NekDouble> mag_cross(nq);
    for (i = 0; i < 2; ++i)
    {
    	grad_u[i] = Array<OneD,NekDouble>(nq);
    	grad_v[i] = Array<OneD,NekDouble>(nq);
    }
    Exp[0]->PhysDeriv(Exp[0]->GetPhys(), grad_u[0], grad_u[1]);
    Exp[1]->PhysDeriv(Exp[1]->GetPhys(), grad_v[0], grad_v[1]);

    // Compute cross product magnitude
    for (i = 0; i < nq; ++i)
    {
    	mag_cross[i] = grad_u[0][i] * grad_v[1][i] - grad_u[1][i] * grad_v[0][i];
    }
    Exp[nfields]->FwdTrans(mag_cross, Exp[nfields]->UpdateCoeffs());


    //-----------------------------------------------
    // Write solution to file with additional computed fields
    string   fldfilename(argv[2]);
    string   out = fldfilename.substr(0, fldfilename.find_last_of("."));
    string   endfile("_add.fld");
    out += endfile;
    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                                                = Exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for(j = 0; j < nfields + addfields; ++j)
    {
		for(i = 0; i < FieldDef.size(); ++i)
		{
			if (j >= nfields)
			{
				FieldDef[i]->m_Fields.push_back("w");
			}
			else
			{
				FieldDef[i]->m_Fields.push_back(fielddef[i]->m_Fields[j]);
			}
			Exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
		}
    }
    graph.Write(out, FieldDef, FieldData);
    //-----------------------------------------------

    return 0;
}

