#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int i,j;
    int surfID;

    if(argc != 5)
    {
        fprintf(stderr,"Usage: FldAddScalGrad meshfile infld outfld BoundaryID\n");
        exit(1);
    }

    surfID = boost::lexical_cast<int>(argv[argc - 1]);

    argv[argc -1] = argv[argc - 2];

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-4]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Import field file.
    string fieldfile(argv[argc-3]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    LibUtilities::Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = 1;
    int addfields = 7;
    Array<OneD, MultiRegions::ExpListSharedPtr> exp(nfields + addfields);
    MultiRegions::AssemblyMapCGSharedPtr m_locToGlobalMap;

    switch(expdim)
    {
        case 1:
        {
            ASSERTL0(false,"Expansion dimension not recognised");
        }
        break;
        case 2:
        {
            ASSERTL0(false,"Expansion dimension not recognised");
        }
        break;
        case 3:
        {
            MultiRegions::ContField3DSharedPtr originalfield =
                MemoryManager<MultiRegions::ContField3D>
                ::AllocateSharedPtr(vSession, graphShPt, 
                                    vSession->GetVariable(0));

            m_locToGlobalMap = originalfield->GetLocalToGlobalMap();
            
            exp[0] = originalfield;
            for (i=0; i<addfields; i++)
            {
                exp[i+1] = MemoryManager<MultiRegions::ContField3D>
                    ::AllocateSharedPtr(*originalfield, graphShPt, 
                                        vSession->GetVariable(0));
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
    for(j = 0; j < nfields+addfields; ++j)
    {
        for(int i = 0; i < fielddata.size(); ++i)
        {
            exp[j]->ExtractDataToCoeffs(fielddef [i],
                                       fielddata[i],
                                        fielddef [i]->m_fields[0],
                                        exp[j]->UpdateCoeffs());
        }
        
        exp[j]->BwdTrans(exp[j]->GetCoeffs(),exp[j]->UpdatePhys());
    }


    //----------------------------------------------

    //----------------------------------------------
    int n, cnt, elmtid, nq, offset, nt, boundary, nfq;
    nt = exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > grad(expdim);
    Array<OneD, Array<OneD, NekDouble> > fgrad(expdim);
    Array<OneD, Array<OneD, NekDouble> > values(addfields);
   
    // Set up mapping from Boundary condition to element details.
    StdRegions::StdExpansionSharedPtr elmt;
    StdRegions::StdExpansion2DSharedPtr bc;
    Array<OneD, int> BoundarytoElmtID;
    Array<OneD, int> BoundarytoTraceID;
    Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> > BndExp(addfields);
    Array<OneD, const NekDouble> U(nt);
    Array<OneD, NekDouble> outvalues;
    
    exp[0]->GetBoundaryToElmtMap(BoundarytoElmtID,BoundarytoTraceID); 

    //get boundary expansions for each field
    for (i = 0; i<addfields; i++)
    {
        BndExp[i] = exp[i]->GetBndCondExpansions();
    }

    // loop over the types of boundary conditions
    for(cnt = n = 0; n < BndExp[0].num_elements(); ++n)
    {   
        // identify boundary which the user wanted
        if(n == surfID)
        {   
            for(i = 0; i < BndExp[0][n]->GetExpSize(); ++i, cnt++)
            {
                // find element and face of this expansion.
                elmtid = BoundarytoElmtID[cnt];
                elmt   = exp[0]->GetExp(elmtid);
                nq     = elmt->GetTotPoints();
                offset = exp[0]->GetPhys_Offset(elmtid);
                
                // Initialise local arrays for the velocity gradients
                // size of total number of quadrature points for each element (hence local).
                for(j = 0; j < expdim; ++j)
                {
                    grad[j] = Array<OneD, NekDouble>(nq);
                }
                
                if(expdim == 2)
                { 
                }
                else
                {   
                    for (j = 0; j< addfields; j++)
                    {
                        values[j] = BndExp[j][n]->UpdateCoeffs() + BndExp[j][n]->GetCoeff_Offset(i);
                    }
                   
                    // Get face 2D expansion from element expansion
                    bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion2D> (BndExp[0][n]->GetExp(i));

                    // Number of face quadrature points
                    nfq = bc->GetTotPoints();

                    //identify boundary of element
                    boundary = BoundarytoTraceID[cnt];
                    
                    //Extract scalar field
                    U = exp[0]->GetPhys() + offset;

                    //Compute gradients
                    elmt->PhysDeriv(U,grad[0],grad[1],grad[2]);
                    
                    if(i ==0)
                    {
                        for (j = 0; j< nq; j++)
                        {
                            cout << "element grad: " << grad[0][j] << endl;
                        }
                    }

                    for(j = 0; j < expdim; ++j)
                    {
                        fgrad[j] = Array<OneD, NekDouble>(nfq);
                    }


                    // Get gradient at the quadrature points of the face
                    for(j = 0; j < expdim; ++j)
                    {
                        elmt->GetFacePhysVals(boundary,bc,grad[j],fgrad[j]); 
                        bc->FwdTrans(fgrad[j],values[j]);
                    }

                    if(i ==0)
                    {
                        for (j = 0; j< nfq; j++)
                        {
                            cout << "face grad: " << fgrad[0][j] << endl;
                        }
                    }

                    const SpatialDomains::GeomFactorsSharedPtr m_metricinfo=bc->GetMetricInfo();

                    const Array<OneD, const Array<OneD, NekDouble> > normals
                                = elmt->GetFaceNormal(boundary);

                    Array<OneD, NekDouble>  gradnorm(nfq);

                    if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                    {
                        Vmath::Vvtvvtp(nfq,normals[0],1,fgrad[0],1,
                                       normals[1],1,fgrad[1],1,gradnorm,1);
                        Vmath::Vvtvp  (nfq,normals[2],1,fgrad[2],1,gradnorm,1,gradnorm,1);
                    }
                    else
                    {
                        Vmath::Svtsvtp(nfq,normals[0][0],fgrad[0],1,
                                          normals[1][0],fgrad[1],1,gradnorm,1);
                        Vmath::Svtvp(nfq,normals[2][0],fgrad[2],1,gradnorm,1,gradnorm,1);
                    }
                    
                    for(j = 0; j<expdim; j++)
                    {
                        bc->FwdTrans(normals[j],values[j+expdim]);
                    }

                    //gradient (grad(u) n)
                    Vmath::Smul(nfq,-1.0,gradnorm,1,gradnorm,1);
                    bc->FwdTrans(gradnorm,values[expdim*2]);
                }
            }
            
        }
        else 
        {
            cnt += BndExp[0][n]->GetExpSize();
        }
    }

    for(int j = 0; j < addfields; ++j)
    {
        int ncoeffs = exp[0]->GetNcoeffs();
        Array<OneD, NekDouble> output(ncoeffs);
        
        output=exp[j+1]->UpdateCoeffs();
        
        int nGlobal=m_locToGlobalMap->GetNumGlobalCoeffs();
        Array<OneD, NekDouble> outarray(nGlobal,0.0);
        
        int bndcnt=0;
        
        const Array<OneD,const int>& map = m_locToGlobalMap->GetBndCondCoeffsToGlobalCoeffsMap();
        NekDouble sign;
        
        for(int i = 0; i < BndExp[j].num_elements(); ++i)
        {
            if(i==surfID)
            {
                const Array<OneD,const NekDouble>& coeffs = BndExp[j][i]->GetCoeffs();
                for(int k = 0; k < (BndExp[j][i])->GetNcoeffs(); ++k)
                {
                    sign = m_locToGlobalMap->GetBndCondCoeffsToGlobalCoeffsSign(bndcnt);
                    outarray[map[bndcnt++]] = sign * coeffs[k];
                }
            }
            else
            {
                bndcnt += BndExp[j][i]->GetNcoeffs();
            }
        }
        m_locToGlobalMap->GlobalToLocal(outarray,output);
    }
    

    //-----------------------------------------------
    // Write solution to file with additional computed fields
    string   out(argv[argc-2]);
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
                                                = exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
    
    vector<string > outname;
    
    outname.push_back("du/dx");
    outname.push_back("du/dy");
    outname.push_back("du/dz");
    outname.push_back("nx");
    outname.push_back("ny");
    outname.push_back("nz");
    outname.push_back("gradient");
    
    
    for(j = 0; j < nfields+addfields; ++j)
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
            exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }

    LibUtilities::Write(out, FieldDef, FieldData);
    //-----------------------------------------------

    return 0;
}
