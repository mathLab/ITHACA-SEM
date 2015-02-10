#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <SolverUtils/Mapping/Mapping.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int i,j;

    if(argc != 3)
    {
        fprintf(stderr,"Usage: FldAddVort  meshfile infld \n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);


    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[1]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);//meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // Import field file.
    string fieldfile(argv[2]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    LibUtilities::Import(fieldfile,fielddef,fielddata);
    bool useFFT = false;
    vSession->MatchSolverInfo("USEFFT", "FFTW", useFFT, false);
    bool dealiasing = false;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_fields.size();
    int veldim = expdim;
    if (fielddef[0]->m_numHomogeneousDir > 0)
    {
	veldim = veldim + fielddef[0]->m_numHomogeneousDir;
    }
    int addfields = nfields*veldim;
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields + addfields);

    switch(expdim)
    {
    case 1:
        {
            ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 2,"Quasi-3D approach is only set up for 1 or 2 homogeneous directions");

            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                MultiRegions::ExpList2DHomogeneous1DSharedPtr Exp2DH1;

                // Define Homogeneous expansion
                //int nplanes = fielddef[0]->m_numModes[1];
                int nplanes;
                vSession->LoadParameter("HomModesZ",nplanes,fielddef[0]->m_numModes[1]);

                // choose points to be at evenly spaced points at
                // nplanes points
                const LibUtilities::PointsKey Pkey(nplanes,LibUtilities::eFourierEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[1],nplanes,Pkey);
                NekDouble ly = fielddef[0]->m_homogeneousLengths[0];

                Exp2DH1 = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,ly,useFFT,dealiasing,graphShPt);
                Exp[0] = Exp2DH1;

                for(i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(*Exp2DH1);
                }
            }
            else if(fielddef[0]->m_numHomogeneousDir == 2)
            {
                MultiRegions::ExpList3DHomogeneous2DSharedPtr Exp3DH2;

                // Define Homogeneous expansion
                //int nylines = fielddef[0]->m_numModes[1];
                //int nzlines = fielddef[0]->m_numModes[2];

				int nylines;
				int nzlines;
				vSession->LoadParameter("HomModesY",nylines,fielddef[0]->m_numModes[1]);
				vSession->LoadParameter("HomModesZ",nzlines,fielddef[0]->m_numModes[2]);

                // choose points to be at evenly spaced points at
                // nplanes points
                const LibUtilities::PointsKey PkeyY(nylines,LibUtilities::eFourierEvenlySpaced);
                const LibUtilities::BasisKey  BkeyY(fielddef[0]->m_basis[1],nylines,PkeyY);

                const LibUtilities::PointsKey PkeyZ(nylines,LibUtilities::eFourierEvenlySpaced);
                const LibUtilities::BasisKey  BkeyZ(fielddef[0]->m_basis[2],nzlines,PkeyZ);

                NekDouble ly = fielddef[0]->m_homogeneousLengths[0];
                NekDouble lz = fielddef[0]->m_homogeneousLengths[1];

                Exp3DH2 = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::AllocateSharedPtr(vSession,BkeyY,BkeyZ,ly,lz,useFFT,dealiasing,graphShPt);
                Exp[0] = Exp3DH2;

                for(i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::AllocateSharedPtr(*Exp3DH2);
                }
            }
            else
            {
                MultiRegions::ExpList1DSharedPtr Exp1D;
                Exp1D = MemoryManager<MultiRegions::ExpList1D>
                                                        ::AllocateSharedPtr(vSession,graphShPt);
                Exp[0] = Exp1D;
                for(i = 1; i < nfields + addfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList1D>
                                                        ::AllocateSharedPtr(*Exp1D);
                }
            }
        }
        break;
    case 2:
        {
            ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 1,"NumHomogeneousDir is only set up for 1");

            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;

                // Define Homogeneous expansion
                //int nplanes = fielddef[0]->m_numModes[2];

                int nplanes;
                vSession->LoadParameter("HomModesZ",nplanes,fielddef[0]->m_numModes[2]);

                // choose points to be at evenly spaced points at
                // nplanes points
                const LibUtilities::PointsKey Pkey(nplanes,LibUtilities::eFourierEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[2],nplanes,Pkey);
                NekDouble lz = fielddef[0]->m_homogeneousLengths[0];

                Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,lz,useFFT,dealiasing,graphShPt);
                Exp[0] = Exp3DH1;

                for(i = 1; i < nfields + addfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(*Exp3DH1);

                }
            }
            else
            {
                MultiRegions::ExpList2DSharedPtr Exp2D;
                Exp2D = MemoryManager<MultiRegions::ExpList2D>
                                                        ::AllocateSharedPtr(vSession,graphShPt);
                Exp[0] =  Exp2D;

                for(i = 1; i < nfields + addfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2D>
                        ::AllocateSharedPtr(*Exp2D);
                }
            }
        }
        break;
    case 3:
        {
            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = MemoryManager<MultiRegions::ExpList3D>
                                                    ::AllocateSharedPtr(vSession,graphShPt);
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
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;

    std::vector<std::vector<NekDouble> > FieldData;

    //----------------------------------------------
    // Copy data from field file
    for(j = 0; j < nfields; ++j)
    {
        for(i = 0; i < fielddata.size(); ++i)
        {
        Exp[j]->ExtractDataToCoeffs(fielddef [i],
                                        fielddata[i],
                                    fielddef [i]->m_fields[j],
                                        Exp[j]->UpdateCoeffs());
        }
        Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
    }
    //----------------------------------------------
    // Check if mapping was defined
    SolverUtils::MappingSharedPtr mapping = SolverUtils::Mapping::Load(vSession, Exp);
    if( mapping)
    {
        //Convert velocity to Cartesian system
        int nq = Exp[0]->GetNpoints();
        Array<OneD, Array<OneD, NekDouble> > tmp1 (veldim);
        Array<OneD, Array<OneD, NekDouble> > tmp2 (veldim);
        for ( i =0; i<veldim; ++i )
        {
            tmp1[i] = Array<OneD, NekDouble> (nq);
            tmp2[i] = Array<OneD, NekDouble> (nq);                    
            Vmath::Vcopy(nq, Exp[i]->GetPhys(), 1, tmp1[i], 1);
        }
        mapping->ContravarToCartesian(tmp1, tmp2);
        for ( i =0; i<veldim; ++i )
        {
            Vmath::Vcopy(nq, tmp2[i], 1, Exp[i]->UpdatePhys(), 1);
        }
    }

    //----------------------------------------------
    // Compute gradients of fields
    int nq = Exp[0]->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble> > outfield(addfields);
    for(i = 0; i < addfields; ++i)
    {
            outfield[i] = Array<OneD, NekDouble>(nq, 0.0);
    }            
    
    // Calculate Gradients
    for(i = 0; i < nfields; ++i)
    {
        for(j = 0; j < veldim; ++j)
        {
            Exp[i]->PhysDeriv(MultiRegions::DirCartesianMap[j],Exp[i]->GetPhys(),outfield[i*veldim+j]);
        }      
    }

    if( mapping)
    {
        // Calculate gradient wrt Cartesian coordinates
        Array<OneD, Array<OneD, NekDouble> > tmp1 (veldim);
        Array<OneD, Array<OneD, NekDouble> > tmp2 (veldim);
        for(i = 0; i < nfields; ++i)
        {                    
            for ( j =0; j<veldim; ++j )
            {
                tmp1[j] = Array<OneD, NekDouble> (nq);
                tmp2[j] = Array<OneD, NekDouble> (nq);
                Vmath::Vcopy(nq, outfield[i*veldim+j], 1, tmp1[j], 1);
            }
            mapping->ContravarToCartesian(tmp1, tmp2);
            for ( j =0; j<veldim; ++j )
            {
                Vmath::Vcopy(nq, tmp2[j], 1, outfield[i*veldim+j], 1);
            }                    
        }
    }

    for (i = 0; i < addfields; ++i)
    {
        Exp[nfields + i]->FwdTrans(outfield[i], Exp[nfields+i]->UpdateCoeffs());
    }
    //-----------------------------------------------
    // Write solution to file with additional computed fields
    string   fldfilename(argv[2]);
    string   out = fldfilename.substr(0, fldfilename.find_last_of("."));
    string   endfile("_with_grad.fld");
    out += endfile;
    
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = Exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    vector<string > outname;

    for (i = 0; i<nfields; ++i)
    {
	if(veldim == 1)
	{
	    outname.push_back(fielddef[0]->m_fields[i]+"_x");
	}
	else if (veldim == 2)
	{
	    outname.push_back(fielddef[0]->m_fields[i]+"_x");
	    outname.push_back(fielddef[0]->m_fields[i]+"_y");
	}
	else if (veldim == 3)
	{
	    outname.push_back(fielddef[0]->m_fields[i]+"_x");
	    outname.push_back(fielddef[0]->m_fields[i]+"_y");
	    outname.push_back(fielddef[0]->m_fields[i]+"_z");
	}
    }

    for(j = 0; j < nfields + addfields; ++j)
    {
            for(i = 0; i < FieldDef.size(); ++i)
            {
            if (j >= nfields)
            {
                    FieldDef[i]->m_fields.push_back(outname[j-nfields]);
            }
            else
            {
                    FieldDef[i]->m_fields.push_back(fielddef[i]->m_fields[j]);
            }
            Exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
            }
    }
    LibUtilities::Write(out, FieldDef, FieldData);
    //-----------------------------------------------

    return 0;
}

