#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int i,j;

    if(argc != 3)
    {
        fprintf(stderr,"Usage: AddExprToField  meshfile fieldfile\n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);


    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);//meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // Import field file.
    string fieldfile(argv[argc-1]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graphShPt->Import(fieldfile,fielddef,fielddata);
    bool useFFT = false;
    bool deal = false;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_fields.size();
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields);
	
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
                const LibUtilities::PointsKey Pkey(nplanes+1,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[1],nplanes,Pkey);
                NekDouble ly = fielddef[0]->m_homogeneousLengths[0];

                Exp2DH1 = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,ly,useFFT,deal,graphShPt);
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
                const LibUtilities::PointsKey PkeyY(nylines+1,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  BkeyY(fielddef[0]->m_basis[1],nylines,PkeyY);
				
				const LibUtilities::PointsKey PkeyZ(nzlines+1,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  BkeyZ(fielddef[0]->m_basis[2],nzlines,PkeyZ);
                
				NekDouble ly = fielddef[0]->m_homogeneousLengths[0];
				NekDouble lz = fielddef[0]->m_homogeneousLengths[1];
				
                Exp3DH2 = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>::AllocateSharedPtr(vSession,BkeyY,BkeyZ,ly,lz,useFFT,deal,graphShPt);
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
                for(i = 1; i < nfields ; ++i)
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
                // nplanes + 1 points
                const LibUtilities::PointsKey Pkey(nplanes+1,LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(fielddef[0]->m_basis[2],nplanes,Pkey);
                NekDouble lz = fielddef[0]->m_homogeneousLengths[0];

                Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>::AllocateSharedPtr(vSession,Bkey,lz,useFFT,deal,graphShPt);
                Exp[0] = Exp3DH1;

                for(i = 1; i < nfields; ++i)
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

                for(i = 1; i < nfields; ++i)
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

            for(i = 1; i < nfields; ++i)
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
                                        fielddef [i]->m_fields[j],
                                        Exp[j]->UpdateCoeffs());
        }
        Exp[j]->BwdTrans_IterPerExp(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
    }
    //----------------------------------------------

    //----------------------------------------------
    // Add expression to field
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    int nq = Exp[0]->GetNpoints();
    Array<OneD, NekDouble> x(nq);
    Array<OneD, NekDouble> y(nq);
    Exp[0]->GetCoords(x,y);
    //Array<OneD, NekDouble> tmp(nq, 1000);
    NekDouble pi = 3.14159265;
    NekDouble tmp;
    //the Exp.num_elements()==1!!!
    ASSERTL0(Exp.num_elements()==1, "the field is not a streak");
cout<<"before  Exp[0][1]="<<Exp[0]->GetPhys()[9]<<endl;
    for (int i = 0; i < nq; ++i)
    {
           //sin(pi*y)*cos(x)
           tmp = 0.01*sin(pi*y[i])*cos(x[i]);
//cout<<"add[0][1]="<<tmp<<endl;
           //Vmath::Vadd(nq, Exp[0]->GetPhys(),1,tmp,1,Exp[0]->UpdatePhys(),1);    
           Exp[0]->UpdatePhys()[i] = Exp[0]->GetPhys()[i] +tmp;
    }
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    //-----------------------------------------------
    // Write solution to file with additional computed fields
    string   fldfilename(argv[2]);
    string   out = fldfilename.substr(0, fldfilename.find_last_of("."));
    string   endfile("_add.fld");
    out += endfile;
    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                                                = Exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
    Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(Exp.num_elements());   	

    for(int j = 0; j < nfields ; ++j)
    {
                Exp[j]->FwdTrans_IterPerExp(Exp[j]->GetPhys(),Exp[j]->UpdateCoeffs());	
cout<<"  Exp[0][0]="<<Exp[0]->GetPhys()[9]<<endl;	     
		fieldcoeffs[j] = Exp[j]->UpdateCoeffs();
		for(int i = 0; i < FieldDef.size(); ++i)
		{
			FieldDef[i]->m_fields.push_back(fielddef[i]->m_fields[j]);
			Exp[j]->AppendFieldData(FieldDef[i], FieldData[i], fieldcoeffs[j]);
		}
    }
    graphShPt->Write(out, FieldDef, FieldData);
    //-----------------------------------------------

    return 0;
}

