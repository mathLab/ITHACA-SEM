#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <iostream>

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
    int nfiles, nStart, surfID;
    if(argc != 5)
    {
        fprintf(stderr,"Usage: FldAddMultiShear meshfile nfiles FirstInFile BoundaryID\n");
        exit(1);
    }

    nfiles = boost::lexical_cast<int>(argv[2]);
    surfID = boost::lexical_cast<int>(argv[4]);

    vector<string> infiles(nfiles);
    string outfile;

    // Output file name: multishear.fld
    stringstream filename2;
    filename2 <<  "multishear.fld";
    filename2 >> outfile;

    // starting checkpoint file name: name_tn_wss.fld. n can be any number.
    string basename = argv[3];
    basename = basename.substr(basename.find_last_of("t")+1, basename.find_last_of(".")-basename.find_last_of("t"));
    stringstream filename3;
    filename3 << basename;
    filename3 >> nStart;

    for (i = 0; i< nfiles; ++i)
    {
        basename = argv[3];
        string extension = ".fld";
        basename = basename.substr(0, basename.find_first_of("_"));
        stringstream filename;
        filename << basename << "_t" << i+nStart << "_wss.fld";
        filename >> infiles[i];
        cout << infiles[i]<<endl;
    }

    argv[2] = argv[3];
    argv[4] = argv[3];

    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv);


    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[1]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------


    //----------------------------------------------
    // Import first field file.
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    LibUtilities::Import(infiles[0],fielddef,fielddata);
    //----------------------------------------------


    //----------------------------------------------
    // Define Expansion
    //----------------------------------------------
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_fields.size();
    int addfields = 7;
    int sfields = nfields - expdim;
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields+addfields), shear(sfields), extraVar(addfields);
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
            i = 0;
            MultiRegions::ContField3DSharedPtr firstfield =
                MemoryManager<MultiRegions::ContField3D>
                ::AllocateSharedPtr(vSession, graphShPt,
                                    vSession->GetVariable(i));

            m_locToGlobalMap = firstfield->GetLocalToGlobalMap();

            Exp[0] = firstfield;
            for(i = 1; i < expdim; ++i)
            {
                Exp[i] = MemoryManager<MultiRegions::ContField3D>
                    ::AllocateSharedPtr(*firstfield, graphShPt,
                                        vSession->GetVariable(i));
            }

            for(i = 0; i < sfields; ++i)
            {
                Exp[i+expdim] = MemoryManager<MultiRegions::ContField3D>
                    ::AllocateSharedPtr(*firstfield, graphShPt,
                                        vSession->GetVariable(i));
            }

            for(i = 0; i < addfields; ++i)
            {
                Exp[i+nfields] = MemoryManager<MultiRegions::ContField3D>
                    ::AllocateSharedPtr(*firstfield, graphShPt,
                                        vSession->GetVariable(0));
            }

        }
        break;
        default:
            ASSERTL0(false,"Expansion dimension not recognised");
            break;
    }


    //-------------------------------------
    // FIRST FILE - INITIALISING EVERYTHING
    //-------------------------------------

    // Set up mapping from Boundary condition to element details.
    StdRegions::StdExpansionSharedPtr elmt;
    StdRegions::StdExpansion2DSharedPtr bc;
    Array<OneD, int> BoundarytoElmtID;
    Array<OneD, int> BoundarytoTraceID;
    Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >  BndExp(nfields+addfields);

    // Copy data from field file
    for(i = 0; i < fielddata.size(); ++i)
    {
        Exp[0]->ExtractDataToCoeffs(fielddef[i],
                                    fielddata[i],
                                    fielddef[i]->m_fields[0],
                                    Exp[0]->UpdateCoeffs());
    }
    Exp[0]->BwdTrans(Exp[0]->GetCoeffs(),Exp[0]->UpdatePhys());

    Exp[0]->GetBoundaryToElmtMap(BoundarytoElmtID,BoundarytoTraceID);
    BndExp[0]   = Exp[0]->GetBndCondExpansions();

    // Get face 2D expansion from element expansion
    bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion2D> (BndExp[0][surfID]->GetExp(0));

    int nfq=  bc->GetTotPoints();
    int nelem = BndExp[0][surfID]->GetExpSize();
    int nbq = nelem*nfq;

    int nt = Exp[0]->GetNpoints();
    Array<OneD, const NekDouble> testing(nt);


    // Define local arrays for wss components, and outputs (TAwss, osi, trs)
    int n, cnt, elmtid, offset, boundary, bndOffset;
    Array<OneD, NekDouble> Sx(nbq), Sy(nbq), Sz(nbq), S(nbq), Sxr(nbq), Syr(nbq), Szr(nbq);
    Array<OneD, NekDouble> Save(nbq), temp2(nbq), trs(nbq), TAwss(nbq), osi(nbq);
    Array<OneD, Array<OneD, NekDouble> > temp(sfields), values(nfields);

    for (i = 0; i < sfields; i++)
    {
        temp[i]= Array<OneD, NekDouble>(nfq);
    }

    Vmath::Zero (nbq, Sx, 1);
    Vmath::Zero (nbq, Sy, 1);
    Vmath::Zero (nbq, Sz, 1);

    // -----------------------------------------------------
    // Compute temporal average wall shear stress vector,
    // and the spatial average of the temporal average.
    // -----------------------------------------------------

    for (int fileNo = 0; fileNo < nfiles ; ++fileNo)
    {
        // Import field file.
        fielddef.clear();
        fielddata.clear();
        LibUtilities::Import(infiles[fileNo],fielddef,fielddata);

        // Copy data from field file
        for(j = 0; j < nfields; ++j)
        {
            for(int i = 0; i < fielddata.size(); ++i)
            {
                Exp[j]->ExtractDataToCoeffs(fielddef[i],
                                            fielddata[i],
                                            fielddef[i]->m_fields[j],
                                            Exp[j]->UpdateCoeffs());
            }
            Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
        }

        Exp[0]->GetBoundaryToElmtMap(BoundarytoElmtID,BoundarytoTraceID);
        BndExp[0]   = Exp[0]->GetBndCondExpansions();

        for(cnt = n = 0; n < BndExp[0].num_elements(); ++n)
        {
            if(n == surfID)
            {
                for(i = 0; i < nelem; ++i, cnt++)
                {
                    // find element and face of this expansion.
                    elmtid = BoundarytoElmtID[cnt];
                    elmt   = Exp[0]->GetExp(elmtid);
                    offset = Exp[0]->GetPhys_Offset(elmtid);
                    bndOffset = nfq*i;

                    if(expdim == 2)
                    {
                        // Not implemented in 2D.
                    }
                    else
                    {
                        // Get face 2D expansion from element expansion
                        bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion2D> (BndExp[0][n]->GetExp(i));

                        //identify boundary of element looking at.
                        boundary = BoundarytoTraceID[cnt];

                        // Get face stress values.
                        for (int t = 0; t< expdim; t++)
                        {
                            elmt->GetFacePhysVals(boundary,bc, Exp[t+expdim]->GetPhys() + offset, temp[t]);
                        }

                        for (int m = 0; m < nfq; m++)
                        {
                            Sxr[bndOffset + m] = temp[0][m];
                            Syr[bndOffset + m] = temp[1][m];
                            Szr[bndOffset + m] = temp[2][m];
                        }
                    }
                }

                // Sx = Sx + Sxr;
                Vmath::Vadd (nbq, Sxr, 1, Sx, 1, Sx, 1);
                Vmath::Vadd (nbq, Syr, 1, Sy, 1, Sy, 1);
                Vmath::Vadd (nbq, Szr, 1, Sz, 1, Sz, 1);

                Vmath::Zero (nbq, Sxr, 1);
                Vmath::Zero (nbq, Syr, 1);
                Vmath::Zero (nbq, Szr, 1);
            }
            else
            {
                cnt += BndExp[0][n]->GetExpSize();
            }
        }
    }

    // Temporal averages of each wss component: Sx = Sx / nfiles. I.e. mean temporal wss
    Vmath::Smul(nbq, 1.0/nfiles , Sx, 1, Sx, 1);
    Vmath::Smul(nbq, 1.0/nfiles , Sy, 1, Sy, 1);
    Vmath::Smul(nbq, 1.0/nfiles , Sz, 1, Sz, 1);

    // Spatial average of the temporal averaged wss vector: Save = sqrt(sx^2 + Sy^2 + Sz^2);
    // i.e magnitude of mean temporal wss.
    Vmath::Vvtvvtp(nbq, Sx, 1, Sx, 1, Sy, 1, Sy, 1, Save, 1);
    Vmath::Vvtvp(nbq, Sz, 1, Sz, 1, Save, 1, Save, 1);
    Vmath::Vsqrt(nbq, Save, 1, Save, 1);

    //temporal averaged vector / spatial average (sort of normalisation): Sxr = Sx/Save;
    //unit vector of mean temporal wss (t).
    Vmath::Vdiv(nbq, Sx, 1, Save, 1, Sxr, 1);
    Vmath::Vdiv(nbq, Sy, 1, Save, 1, Syr, 1);
    Vmath::Vdiv(nbq, Sz, 1, Save, 1, Szr, 1);

    Vmath::Zero (nbq, Sx, 1);
    Vmath::Zero (nbq, Sy, 1);
    Vmath::Zero (nbq, Sz, 1);
    Vmath::Zero (nbq, trs, 1);
    Vmath::Zero (nbq, TAwss, 1);

    // -----------------------------------------------------
    // Loop through files again, and compute transverse wss
    // -----------------------------------------------------

    for (int fileNo = 0; fileNo < nfiles ; ++fileNo)
    {
        Vmath::Zero (nbq, Sx, 1);
        Vmath::Zero (nbq, Sy, 1);
        Vmath::Zero (nbq, Sz, 1);

        // Import field file.
        fielddef.clear();
        fielddata.clear();
        LibUtilities::Import(infiles[fileNo],fielddef,fielddata);

        // Copy data from field file
        for(j = 0; j < nfields; ++j)
        {
            for(int i = 0; i < fielddata.size(); ++i)
            {
                Exp[j]->ExtractDataToCoeffs(fielddef[i],
                                            fielddata[i],
                                            fielddef[i]->m_fields[j],
                                            Exp[j]->UpdateCoeffs());
            }
            Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
        }

        Exp[0]->GetBoundaryToElmtMap(BoundarytoElmtID,BoundarytoTraceID);
        BndExp[0]   = Exp[0]->GetBndCondExpansions();

        for(cnt = n = 0; n < BndExp[0].num_elements(); ++n)
        {
            if(n == surfID)
            {
                for(i = 0; i < nelem; ++i, cnt++)
                {
                    // find element and face of this expansion.
                    elmtid = BoundarytoElmtID[cnt];
                    elmt   = Exp[0]->GetExp(elmtid);
                    offset = Exp[0]->GetPhys_Offset(elmtid);
                    bndOffset = nfq*i;

                    if(expdim == 2)
                    {
                    }
                    else
                    {
                        // Get face 2D expansion from element expansion
                        bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion2D> (BndExp[0][n]->GetExp(i));

                        //identify boundary of element looking at.
                        boundary = BoundarytoTraceID[cnt];

                        // Get face stress values.
                        for (int t = 0; t< sfields; t++)
                        {
                            elmt->GetFacePhysVals(boundary,bc,Exp[t+expdim]->GetPhys() + offset, temp[t]);
                        }

                        for (int m = 0; m < nfq; m++)
                        {
                            Sx[bndOffset + m] = temp[0][m];
                            Sy[bndOffset + m] = temp[1][m];
                            Sz[bndOffset + m] = temp[2][m];
                            //S [bndOffset + m] = temp[3][m];
                        }
                    }
                }

                Vmath::Vvtvvtp(nbq, Sx, 1, Sx, 1, Sy, 1, Sy, 1, S, 1);
                Vmath::Vvtvp(nbq, Sz, 1, Sz, 1, S, 1, S, 1);
                Vmath::Vsqrt(nbq, S, 1, S, 1);

                // transverse wall shear stress
                // trs = trs + sqrt(s^2 - (Sx*Sxr + Sy*Syr + Sz*Szr)^2)
                Vmath::Vvtvvtp(nbq, Sx, 1, Sxr, 1, Sy, 1, Syr, 1, temp2, 1);
                Vmath::Vvtvp(nbq, Sz, 1, Szr, 1, temp2, 1, temp2, 1);
                Vmath::Vmul(nbq, temp2, 1, temp2, 1, temp2, 1);
                Vmath::Vvtvm(nbq, S, 1, S, 1, temp2, 1, temp2, 1);
                for (int m = 0; m < nbq; m++)
                {
                    if (temp2[m]> 0.0)
                    {
                        trs[m] = trs[m] + sqrt(temp2[m]);
                    }  /*
                    else
                    {
                        trs[m] = 0.0;
                        }*/
                }

                Vmath::Zero (nbq, temp2, 1);

                //Time averaged wss TAwss = TAwss + S;
                Vmath::Vadd (nbq, S, 1, TAwss, 1, TAwss, 1);
            }
            else
            {
                cnt += BndExp[0][n]->GetExpSize();
            }
        }
    }


    // Final step in computing trs and TAwss: trs = trs/nfiles;
    Vmath::Smul (nbq, (1.0/nfiles), trs, 1, trs, 1);
    Vmath::Smul (nbq, (1.0/nfiles), TAwss, 1, TAwss, 1);

    //Compute osi = 0.5*(1- Save/TAwss)
    Vmath::Vdiv(nbq, Save, 1, TAwss, 1, osi, 1);
    Vmath::Smul(nbq, -0.5, osi, 1, osi, 1);
    Vmath::Sadd(nbq, 0.5, osi, 1, osi, 1);

    fielddef.clear();
    fielddata.clear();
    LibUtilities::Import(infiles[nfiles-1],fielddef,fielddata);

    // Copy u,v,w, data from last field file for the fields
    for(j = 0; j < nfields; ++j)
    {
        for(int i = 0; i < fielddata.size(); ++i)
        {
            Exp[j]->ExtractDataToCoeffs(fielddef[i],
                                        fielddata[i],
                                        fielddef[i]->m_fields[j],
                                        Exp[j]->UpdateCoeffs());
        }
        Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
    }

    for(j = 0; j < addfields; ++j)
    {
        for(int i = 0; i < fielddata.size(); ++i)
        {
            Exp[j+nfields]->ExtractDataToCoeffs(fielddef[i],
                                        fielddata[i],
                                        fielddef[i]->m_fields[j],
                                        Exp[j+nfields]->UpdateCoeffs());
        }
        Exp[j+nfields]->BwdTrans(Exp[j+nfields]->GetCoeffs(),Exp[j+nfields]->UpdatePhys());
    }

    Exp[0]->GetBoundaryToElmtMap(BoundarytoElmtID,BoundarytoTraceID);

    //get boundary expansions for each field
    for(j = 0; j < nfields+addfields; ++j)
    {
        BndExp[j]   = Exp[j]->GetBndCondExpansions();
    }

    for(cnt = n = 0; n < BndExp[0].num_elements(); ++n)
    {
        if(n == surfID)
        {
            for(i = 0; i < nelem; ++i, cnt++)
            {
                // find element and face of this expansion.
                elmtid = BoundarytoElmtID[cnt];
                offset = Exp[0]->GetPhys_Offset(elmtid);
                bndOffset = nfq*i;

                if(expdim == 2)
                {
                }
                else
                {
                    // Get face 2D expansion from element expansion
                    bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion2D> (BndExp[0][n]->GetExp(i));

                    //identify boundary of element looking at.
                    boundary = BoundarytoTraceID[cnt];

                    //Update outfield coefficients in the elemental boundary expansion
                    for (j = 0; j < nfq; j++)
                    {
                        temp[0][j] = trs[bndOffset + j];
                        temp[1][j] = TAwss[bndOffset + j];
                        temp[2][j] = osi[bndOffset + j];
                    }

                    for (j = 0; j < 3; j++)
                    {
                        values[j] = BndExp[j+nfields][n]->UpdateCoeffs() + BndExp[j+nfields][n]->GetCoeff_Offset(i);
                        bc->FwdTrans(temp[j], values[j]);

                        Vmath::Zero(nbq, temp[j],1);
                    }

                    for (j = 0; j < nfq; j++)
                    {
                        temp[0][j] = Sxr[bndOffset + j];
                        temp[1][j] = Syr[bndOffset + j];
                        temp[2][j] = Szr[bndOffset + j];
                    }

                    for (j = 0; j < 3; j++)
                    {
                        values[j] = BndExp[j+nfields+3][n]->UpdateCoeffs() + BndExp[j+nfields+3][n]->GetCoeff_Offset(i);
                        bc->FwdTrans(temp[j], values[j]);

                        Vmath::Zero(nbq, temp[j],1);
                    }


                    for (j = 0; j < nfq; j++)
                    {
                        temp[0][j] = Save[bndOffset + j];
                    }

                    values[0] = BndExp[nfields+addfields-1][n]->UpdateCoeffs() + BndExp[nfields+addfields-1][n]->GetCoeff_Offset(i);
                    bc->FwdTrans(temp[0], values[0]);

                }
            }
        }
        else
        {
            cnt += BndExp[0][n]->GetExpSize();
        }
    }

    for(j = 0; j < nfields + addfields; ++j)
    {
        int ncoeffs = Exp[j]->GetNcoeffs();
        Array<OneD, NekDouble> output(ncoeffs);

        output=Exp[j]->UpdateCoeffs();

        int nGlobal=m_locToGlobalMap->GetNumGlobalCoeffs();
        Array<OneD, NekDouble> outarray(nGlobal,0.0);

        int bndcnt=0;

        const Array<OneD,const int>& map = m_locToGlobalMap->GetBndCondCoeffsToGlobalCoeffsMap();
        NekDouble sign;

        for(i = 0; i < BndExp[j].num_elements(); ++i)
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
    // Write solution to file
    // ----------------------------------------------
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = Exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    vector<string > outname;

    outname.push_back("TransWSS");
    outname.push_back("TAWSS");
    outname.push_back("OSI");
    outname.push_back("norm_mean_x");
    outname.push_back("norm_mean_y");
    outname.push_back("norm_mean_z");
    outname.push_back("mean_mag");

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
            Exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }

    LibUtilities::Write(outfile, FieldDef, FieldData);

    return 0;
}

