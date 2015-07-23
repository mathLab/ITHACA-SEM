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
    int i,j,k;
    int surfID, nfiles, nStart;

    if(argc != 5)
    {
        fprintf(stderr,"Usage: FldAddWSS meshfile nfiles infld BoundaryID\n");
        exit(1);
    }

    surfID = boost::lexical_cast<int>(argv[argc - 1]);
    nfiles = boost::lexical_cast<int>(argv[2]);

    vector<string> infiles(nfiles), outfiles(nfiles);

    /* format of files if 1 file: name.fld. Output: name_wss.fld
       format of files if n checkpoint files: name_tn.fld. Output: name_tn_wss.fld, n=0,1,2...
     */


    if (nfiles == 1)
    {
        infiles[0] = argv[3];
        string basename = argv[3];
        basename = basename.substr(0, basename.find_last_of("."));
        stringstream filename2;
        filename2 << basename << "_wss.fld";
        filename2 >> outfiles[0];
    }
    else {
        
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
            stringstream filename, filename2;
            filename << basename << "_t" << i + nStart << extension;
            filename >> infiles[i];
            filename2 << basename << "_t" << i +nStart << "_wss.fld";
            filename2 >> outfiles[i];
        }
    }

    argv[argc - 1] = argv[argc - 2];
    argv[argc - 3] = argv[argc - 2]; 

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv);

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-4]);
    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    //Get viscosity from file
    NekDouble m_kinvis;
    m_kinvis = vSession->GetParameter("Kinvis");
    //----------------------------------------------
    
    //----------------------------------------------
    // Import first field file.
    string fieldfile(argv[3]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    LibUtilities::Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_fields.size()-1;
    int addfields = (nfields == 3)? 4:3;
    int nstress = (nfields == 3)? 6:3;
    Array<OneD, MultiRegions::ExpListSharedPtr> vel(nfields), shear(addfields);
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
            
            ASSERTL0(false,"Not implemented in 2D");

            i=0;
            MultiRegions::ContField2DSharedPtr cont2D;
            cont2D = MemoryManager<MultiRegions::ContField2D>
                ::AllocateSharedPtr(vSession,graphShPt,vSession->GetVariable(i));
            vel[0] =  cont2D;
            
            for(i = 1; i < nfields; ++i)
            {
                vel[i] = MemoryManager<MultiRegions::ContField2D>
                    ::AllocateSharedPtr(*cont2D);
            }
        }
        break;
        case 3:
        {
            MultiRegions::ContField3DSharedPtr firstfield =
                MemoryManager<MultiRegions::ContField3D>
                ::AllocateSharedPtr(vSession, graphShPt, 
                                    vSession->GetVariable(0));

            m_locToGlobalMap = firstfield->GetLocalToGlobalMap();

            vel[0] = firstfield;
            for(i = 1; i < nfields; ++i)
            {
                vel[i] = MemoryManager<MultiRegions::ContField3D>
                    ::AllocateSharedPtr(*firstfield, graphShPt, 
                                        vSession->GetVariable(i));
            }

            for(i = 0; i < addfields; ++i)
            {
                shear[i] = MemoryManager<MultiRegions::ContField3D>
                    ::AllocateSharedPtr(*firstfield, graphShPt, 
                                        vSession->GetVariable(0));
            }
        }
        break;
        default:
            ASSERTL0(false,"Expansion dimension not recognised");
            break;
    }
    //----------------------------------------------


    // Define arrays for WSS (global array), stress components and gradients (local arrays) 
    int n, cnt, elmtid, nq, offset, nt, boundary, nfq;
    Array<OneD, Array<OneD, NekDouble> > grad(nfields*nfields);
    Array<OneD, Array<OneD, NekDouble> > stress(nstress), fstress(nstress);   
    Array<OneD, Array<OneD, NekDouble> > values(addfields);
    
    // Set up mapping from Boundary condition to element details.
    StdRegions::StdExpansionSharedPtr elmt;
    StdRegions::StdExpansion2DSharedPtr bc;
    Array<OneD, int> BoundarytoElmtID, BoundarytoTraceID;
    Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >  BndExp(addfields);
   

    for (int fileNo = 0; fileNo < nfiles ; ++fileNo)
    {
        //----------------------------------------------
        // Import field file. 
        fielddef.clear();
        fielddata.clear();
        LibUtilities::Import(infiles[fileNo],fielddef,fielddata);
        //----------------------------------------------

        //----------------------------------------------
        // Copy data from field file
        for(j = 0; j < nfields; ++j)
        {
            for(i = 0; i < fielddata.size(); ++i)
            {
                vel[j]->ExtractDataToCoeffs(fielddef[i],
                                            fielddata[i],
                                            fielddef[i]->m_fields[j],
                                            vel[j]->UpdateCoeffs());
            }
            
            vel[j]->BwdTrans(vel[j]->GetCoeffs(),vel[j]->UpdatePhys());
        }
        
        for(j = 0; j < addfields; ++j)
        {
            for(i = 0; i < fielddata.size(); ++i)
            {
                shear[j]->ExtractDataToCoeffs(fielddef[i],
                                              fielddata[i],
                                              fielddef[i]->m_fields[0],
                                              shear[j]->UpdateCoeffs());

            }
            
            shear[j]->BwdTrans(shear[j]->GetCoeffs(),shear[j]->UpdatePhys());
        }
        
        //----------------------------------------------
        
        //----------------------------------------------
        // Compute WSS for each element on boundary
        // Define total quadrature points, and velocity fields. 
        nt = shear[0]->GetNpoints();
        Array<OneD, const NekDouble> U(nt),V(nt),W(nt);     

        for(j = 0; j < addfields; ++j)
        {
            values[j] = Array<OneD, NekDouble>(nt);
        }

        
        shear[0]->GetBoundaryToElmtMap(BoundarytoElmtID,BoundarytoTraceID);        

        //get boundary expansions for each field
        for(j = 0; j < addfields; ++j)
        {
            BndExp[j]   = shear[j]->GetBndCondExpansions();
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
                    elmt   = shear[0]->GetExp(elmtid);
                    nq     = elmt->GetTotPoints();
                    offset = shear[0]->GetPhys_Offset(elmtid);

                    // Initialise local arrays for the velocity gradients, and stress components
                    // size of total number of quadrature points for each element (hence local).
                    for(j = 0; j < nfields*nfields; ++j)
                    {
                        grad[j] = Array<OneD, NekDouble>(nq);
                    }
                                        
                    for(j = 0; j < nstress; ++j)
                    {
                        stress[j] = Array<OneD, NekDouble>(nq);
                    }

                   
                    if(nfields == 2)
                    { 
                        //Not implemented in 2D.   
                    }
                    else
                    {   
                        // Get face 2D expansion from element expansion
                        bc =  boost::dynamic_pointer_cast<StdRegions::StdExpansion2D> (BndExp[0][n]->GetExp(i));
                        nfq=  bc->GetTotPoints();
                        
                        //identify boundary of element looking at.
                        boundary = BoundarytoTraceID[cnt];
                        
                        //Get face normals
                        const SpatialDomains::GeomFactorsSharedPtr m_metricinfo=bc->GetMetricInfo();
                        
                        const Array<OneD, const Array<OneD, NekDouble> > normals
                            = elmt->GetFaceNormal(boundary);

                        // initialise arrays
                        for(j = 0; j < nstress; ++j)
                        {
                            fstress[j] = Array<OneD, NekDouble>(nfq);
                        }
                        Array<OneD, NekDouble> values2(nfq), S(nfq), Sx(nfq), Sy(nfq), Sz(nfq);
                 
                        //Extract Velocities
                        U = vel[0]->GetPhys() + offset;
                        V = vel[1]->GetPhys() + offset;
                        W = vel[2]->GetPhys() + offset;
                        
                        //Compute gradients (velocity correction scheme method)
                        elmt->PhysDeriv(U,grad[0],grad[1],grad[2]);
                        elmt->PhysDeriv(V,grad[3],grad[4],grad[5]);
                        elmt->PhysDeriv(W,grad[6],grad[7],grad[8]);
                        
                        //Compute stress component terms
                        // t_xx = 2.mu.Ux
                        Vmath::Smul (nq,(2*m_kinvis),grad[0],1,stress[0],1);
                        // tyy = 2.mu.Vy
                        Vmath::Smul (nq,(2*m_kinvis),grad[4],1,stress[1],1);
                        // tzz = 2.mu.Wz
                        Vmath::Smul (nq,(2*m_kinvis),grad[8],1,stress[2],1);
                        // txy = mu.(Uy+Vx)
                        Vmath::Vadd (nq,grad[1],1,grad[3],1,stress[3],1);
                        Vmath::Smul (nq,m_kinvis,stress[3],1,stress[3],1);
                        // txz = mu.(Uz+Wx)
                        Vmath::Vadd (nq,grad[2],1,grad[6],1,stress[4],1);
                        Vmath::Smul (nq,m_kinvis,stress[4],1,stress[4],1);
                        // tyz = mu.(Vz+Wy)
                        Vmath::Vadd (nq,grad[5],1,grad[7],1,stress[5],1);
                        Vmath::Smul (nq,m_kinvis,stress[5],1,stress[5],1);

                        // Get face stress values.
                        for(j = 0; j < nstress; ++j)
                        {
                            elmt->GetFacePhysVals(boundary,bc,stress[j],fstress[j]); 
                        }
                     
                        //calcuate wss, and update velocity coefficients in the elemental boundary expansion
                        for (j = 0; j< addfields; j++)
                        {
                            values[j] = BndExp[j][n]->UpdateCoeffs() + BndExp[j][n]->GetCoeff_Offset(i);
                        }

                        if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                        {
                            // Sx
                            Vmath::Vvtvvtp(nfq,normals[0],1,fstress[0],1,
                                           normals[1],1,fstress[3],1,Sx,1);
                            Vmath::Vvtvp  (nfq,normals[2],1,fstress[4],1,Sx,1,Sx,1);

                            // Sy
                            Vmath::Vvtvvtp(nfq,normals[0],1,fstress[3],1,
                                           normals[1],1,fstress[1],1,Sy,1);
                            Vmath::Vvtvp  (nfq,normals[2],1,fstress[5],1,Sy,1,Sy,1);
                 
                            // Sz
                            Vmath::Vvtvvtp(nfq,normals[0],1,fstress[4],1,
                                           normals[1],1,fstress[5],1,Sz,1);
                            Vmath::Vvtvp  (nfq,normals[2],1,fstress[2],1,Sz,1,Sz,1);
                        }
                        else
                        {
                            // Sx
                            Vmath::Svtsvtp(nfq,normals[0][0],fstress[0],1,
                                           normals[1][0],fstress[3],1,Sx,1);
                            Vmath::Svtvp(nfq,normals[2][0],fstress[4],1,Sx,1,Sx,1);
                       
                            // Sy
                            Vmath::Svtsvtp(nfq,normals[0][0],fstress[3],1,
                                           normals[1][0],fstress[1],1,Sy,1);
                            Vmath::Svtvp(nfq,normals[2][0],fstress[5],1,Sy,1,Sy,1);
                      
                            // Sz
                            Vmath::Svtsvtp(nfq,normals[0][0],fstress[4],1,
                                           normals[1][0],fstress[5],1,Sz,1);
                            Vmath::Svtvp(nfq,normals[2][0],fstress[2],1,Sz,1,Sz,1);
                        }
                        
                        // T = T - (T.n)n
                        if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                        {
                            Vmath::Vvtvvtp(nfq,normals[0],1,Sx,1,
                                           normals[1],1, Sy,1,values2,1);
                            Vmath::Vvtvp  (nfq,normals[2],1, Sz,1,values2,1,values2,1);

                            // values2 = (T.n)n
                            Vmath::Smul(nfq, -1.0, values2, 1, values2, 1);

                            //Sx
                            Vmath::Vvtvp(nfq,normals[0],1, values2,1,Sx,1,Sx,1);
                            bc->FwdTrans(Sx, values[0]);
                           
                            //Sy
                            Vmath::Vvtvp(nfq,normals[1],1, values2,1,Sy,1,Sy,1);
                            bc->FwdTrans(Sy, values[1]);
                           
                            //Sz
                            Vmath::Vvtvp(nfq,normals[2],1, values2,1,Sz,1,Sz,1);
                            bc->FwdTrans(Sz, values[2]);
                                                        
                        }
                        else
                        {
                            Vmath::Svtsvtp(nfq,normals[0][0],Sx,1,
                                           normals[1][0],Sy,1,values2,1);
                            Vmath::Svtvp(nfq,normals[2][0],Sz,1,values2,1,values2,1);

                            Vmath::Smul(nfq, -1.0, values2, 1, values2, 1);

                            //Sx
                            Vmath::Svtvp(nfq,normals[0][0],values2,1,Sx,1,Sx,1);
                            bc->FwdTrans(Sx, values[0]);
                           
                            //Sy
                            Vmath::Svtvp(nfq,normals[1][0],values2,1,Sy,1,Sy,1);
                            bc->FwdTrans(Sy, values[1]);
                           
                            //Sz
                            Vmath::Svtvp(nfq,normals[2][0],values2,1,Sz,1,Sz,1);
                            bc->FwdTrans(Sz, values[2]);
                        }

                       
                        //Tw 
                        Vmath::Vvtvvtp(nfq, Sx, 1, Sx, 1, Sy, 1, Sy, 1, S, 1);
                        Vmath::Vvtvp(nfq, Sz, 1, Sz, 1, S, 1, S, 1);
                        Vmath::Vsqrt(nfq, S, 1, S, 1);
                        bc->FwdTrans(S, values[3]); 
                    }
                }
            }
            else 
            {
                cnt += BndExp[0][n]->GetExpSize();
            }
        }
        

        for(j = 0; j < addfields; ++j)
        {
            int ncoeffs = shear[j]->GetNcoeffs();
            Array<OneD, NekDouble> output(ncoeffs);
            
            output=shear[j]->UpdateCoeffs();

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
                    for(k = 0; k < (BndExp[j][i])->GetNcoeffs(); ++k)
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
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
            = shear[0]->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
        
        vector<string > outname;
        
        if(addfields == 3)
        {
            // Not implemented in 2D.
        }
        else
        {
            outname.push_back("Tx");
            outname.push_back("Ty");
            outname.push_back("Tz");
            outname.push_back("Tw");
        }

        for(j = 0; j < nfields + addfields; ++j)
        {
            for(i = 0; i < FieldDef.size(); ++i)
            {
                if (j < nfields)
                {
                    FieldDef[i]->m_fields.push_back(fielddef[i]->m_fields[j]);
                    vel[j]->AppendFieldData(FieldDef[i], FieldData[i]);
                }
                else
                {
                    FieldDef[i]->m_fields.push_back(outname[j-nfields]);
                    shear[j-nfields]->AppendFieldData(FieldDef[i], FieldData[i]);
                }
            }
        }
        LibUtilities::Write(outfiles[fileNo], FieldDef, FieldData);
        
        FieldDef.clear();
        FieldData.clear();
        //-----------------------------------------------
    }
    
    return 0;
}
