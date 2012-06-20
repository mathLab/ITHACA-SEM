#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <iomanip>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <SpatialDomains/SpatialData.h>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList0D.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpListHomogeneous1D.h>
#include <MultiRegions/ExpListHomogeneous2D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <LocalRegions/MatrixKey.h>

#include <tinyxml/tinyxml.h>

#include <boost/math/special_functions/fpclassify.hpp>

using namespace Nektar;

void SetFields(SpatialDomains::MeshGraphSharedPtr &mesh,
               vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef,
               LibUtilities::SessionReaderSharedPtr &session,
               Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,
               const vector<std::string>  &variables,
               bool homogeneous);

void GenerateField(Array<OneD,MultiRegions::ExpListSharedPtr> &field0,
                   Array<OneD, NekDouble> x1,
                   Array<OneD, NekDouble> y1,
                   Array<OneD,MultiRegions::ExpListSharedPtr> &field1);

void GenerateFieldHomo(MultiRegions::ExpListSharedPtr field0,
                       Array<OneD, NekDouble> x1,
                       Array<OneD, NekDouble> y1,
                       MultiRegions::ExpListSharedPtr field1);

bool Checkbndmeshes(Array<OneD, NekDouble> x0,
                    Array<OneD, NekDouble> y0,       	       
                    Array<OneD, NekDouble> x1,
                    Array<OneD, NekDouble> y1); 

void Writefield(LibUtilities::SessionReaderSharedPtr vSession,
                const vector<std::string> &variables,
                string fieldfile, SpatialDomains::MeshGraphSharedPtr &graph,
                Array<OneD, MultiRegions::ExpListSharedPtr> &outfield);    


int main(int argc, char *argv[])
{


    if(argc != 5)
    {
        fprintf(stderr,"Usage: ./FieldToField  meshfile0 fieldfile0  meshfile1  fieldfile1\n");
        exit(1);
    }
    
    
    //----------------------------------------------
    string meshfile0  (argv[argc-4]); 
    string fieldfile0 (argv[argc-3]);
    //argc=nfiles0 =2 
    //in this way only the mesh0 is taken to create vSession
    int nfiles0 = 2;
    //nfiles0=5;
    std::vector<std::string> filenames0;
    filenames0.push_back(meshfile0);
    filenames0.push_back(meshfile0);
    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv, filenames0);
            //= LibUtilities::SessionReader::CreateInstance(2, argv);
    // Read in mesh from input file0

    /// ===================================================
    /// \todo Please update using MeshGraph::Read(vSession)
    /// (it's now possible having multiple sessions)
    /// ===================================================

    SpatialDomains::MeshGraphSharedPtr graphShPt = SpatialDomains::MeshGraph::Read(meshfile0);
    //----------------------------------------------          
    // Import fieldfile0.
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graphShPt->Import(fieldfile0,fielddef,fielddata);
    //----------------------------------------------    
    //read info from fldfile
    const vector<std::string> variables = fielddef[0]->m_fields;
    

    bool homo=(fielddef[0]->m_numHomogeneousDir)? true: false;

    // Define Expansion    
    int nfields; 
    nfields = variables.size();
    Array<OneD, MultiRegions::ExpListSharedPtr> fields; 
    fields= Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);  
    SetFields(graphShPt,fielddef, vSession,fields,variables,homo);
    int nq;//pointsper plane
    //-----------------------------------------------
    // Copy data from file:fill fields with the fielddata
    if(fielddef[0]->m_numHomogeneousDir == 1)
    {
        nq = fields[0]->GetPlane(1)->GetTotPoints();
        //THE IM PHYS VALUES ARE WRONG USING bwdTrans !!!
        for(int j = 0; j < nfields; ++j)
        {  
            for(int i = 0; i < fielddata.size(); ++i)
            {
                fields[j]->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_fields[j]);
            } 
            //bwd plane 0
            fields[j]->GetPlane(0)->BwdTrans_IterPerExp(fields[j]->GetPlane(0)->GetCoeffs(),  
                          fields[j]->GetPlane(0)->UpdatePhys() );

            //bwd plane 1
            fields[j]->GetPlane(1)->BwdTrans_IterPerExp(fields[j]->GetPlane(1)->GetCoeffs(), 
                          fields[j]->GetPlane(1)->UpdatePhys() );
            

        }    
    }
    else
    {
        nq = fields[0]->GetTotPoints();
        for(int j = 0; j < nfields; ++j)
        {  	    
            for(int i = 0; i < fielddata.size(); ++i)
            {
                fields[j]->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_fields[j]);
            }             
            fields[j]->BwdTrans_IterPerExp(fields[j]->GetCoeffs(),fields[j]->UpdatePhys());      
        }
    }
    //----------------------------------------------    
/*    
         for(int g=0; g<fields[0]->GetPlane(1)->GetTotPoints(); g++)
         {
cout<<"g="<<g<<"  phys f0="<<fields[0]->GetPlane(0)->GetPhys()[g]<<" f1="<<fields[0]->GetPlane(1)->GetPhys()[g]<<endl;
         }
*/


    // store mesh0 quadrature points    
    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> y0(nq);  
    Array<OneD, NekDouble> z0(nq);

    if(
        //vSession->DefinesSolverInfo("HOMOGENEOUS")
        fielddef[0]->m_numHomogeneousDir == 1

      )
    {
       fields[0]->GetPlane(1)->GetCoords(x0,y0,z0);    
    }
    else
    {
       fields[0]->GetCoords(x0,y0);  
    }
  
    //----------------------------------------------    
  
    //Read in the mesh1
    string meshfile1(argv[argc-2]);    
    //Name of the output fieldfile
    string fieldfile1(argv[argc-1]);           
    //-----------------------------------------------    
  
    //define the output field:
    Array<OneD, MultiRegions::ExpListSharedPtr> outfield; 

    outfield = Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);   
    //set output fields over the graphShPt1 graph
    SpatialDomains::MeshGraphSharedPtr graphShPt1;     
    //remark: homo cases malloc() error or  segmentation fault..
    //remember: there is static cnt in stfields function to define 
    // the homo quantities only for the mesh0 case

    /// ===================================================
    /// \todo Please update using MeshGraph::Read(vSession)
    /// (it's now possible having multiple sessions)
    /// ===================================================

    graphShPt1 = SpatialDomains::MeshGraph::Read(meshfile1);  


    //define second session1..	
    std::vector<std::string> filenames1;       
    filenames1.push_back(meshfile1);   
    filenames1.push_back(meshfile1);     
    //cout<<"mesh1 read"<<endl;
    argc=2;
    argv[1]=argv[3];
    argv[2]=argv[4];
    //cout<<argv[1]<<"  a="<<argv[2]<<endl;
    LibUtilities::SessionReaderSharedPtr vSession1
        = LibUtilities::SessionReader::CreateInstance(argc, argv, filenames1, vSession->GetComm());
    //--------------------------------------------------

    //set explist 1
    SetFields(graphShPt1,fielddef, vSession1, outfield,variables,homo);    	    
    //------------------------------------------------ 

    // store the new points:
    int nq1;
    if( fielddef[0]->m_numHomogeneousDir == 1)
    {
        nq1 = outfield[0]->GetPlane(1)->GetTotPoints();
    }
    else
    {
        nq1 = outfield[0]->GetTotPoints();
    }
    Array<OneD, NekDouble> x1(nq1);
    Array<OneD, NekDouble> y1(nq1);
    Array<OneD, NekDouble> z1(nq1);
    
    if(fielddef[0]->m_numHomogeneousDir == 1)
    {
       outfield[0]->GetPlane(1)->GetCoords(x1,y1,z1);    
    }
    else
    {
       outfield[0]->GetCoords(x1,y1);
    }


    //------------------------------------------------
    
    //check 2Dmeshes compatibilities
    // same max min values x,y...
    bool check;    
    check = Checkbndmeshes(x0, y0, x1, y1);    
    ASSERTL0(check, "meshes not compatible (different borders)");   
    //-----------------------------------------------
    
    //generate the new fields
    if(fielddef[0]->m_numHomogeneousDir == 1)
    {
        for(int t=0; t< nfields; t++)
        {        
    	    GenerateFieldHomo(fields[t], x1, y1, outfield[t]);
        }
    }
    else
    {
        GenerateField(fields, x1, y1, outfield);
    }

    //------------------------------------------------
    
    //write fieldfile
    Writefield(vSession, variables, fieldfile1, graphShPt1,outfield);        
    //------------------------------------------------
}


// Define Expansion       		
void SetFields(SpatialDomains::MeshGraphSharedPtr &mesh,
               vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef,
               LibUtilities::SessionReaderSharedPtr &session,
               Array<OneD,MultiRegions::ExpListSharedPtr> &Exp,
               const vector<std::string> &variables,
               bool homogeneous)
{
    //session reader stuff has to be evaluated only from the
    // first session which refers to mesh0
    static int cnt=0;		
    // Setting parameteres for homogenous problems
    MultiRegions::GlobalSysSolnType solnType;
    NekDouble static LhomX;   ///< physical length in X direction (if homogeneous) 
    NekDouble static LhomY;   ///< physical length in Y direction (if homogeneous)
    NekDouble static LhomZ;   ///< physical length in Z direction (if homogeneous)
    
    bool static DeclareCoeffPhysArrays = true;		
    int static npointsX;   ///< number of points in X direction (if homogeneous)
    int static npointsY;   ///< number of points in Y direction (if homogeneous)
    int static npointsZ;   ///< number of points in Z direction (if homogeneous)	
    int static HomoDirec       = 0;
    bool static useFFT = false;
    bool static deal = false;

    ///Parameter for homogeneous expansions		
    enum HomogeneousType
    {
        eHomogeneous1D,
        eHomogeneous2D,
        eHomogeneous3D,
        eNotHomogeneous
    };
    
    enum HomogeneousType HomogeneousType = eNotHomogeneous;
    
    if(cnt==0)
    { 	
        if(fielddef[0]->m_numHomogeneousDir == 1)
        {	
            //only homo1D is working
            HomogeneousType = eHomogeneous1D;
            npointsZ = fielddef[0]->m_numModes[2];
            LhomZ = fielddef[0]->m_homogeneousLengths[0];
            HomoDirec       = 1;
            
        }
    }		
    cnt++;
    int i;		
    int expdim   = mesh->GetMeshDimension();
    int nvariables = variables.size();
    Exp= Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);  
    // Continuous Galerkin projection
    
    switch(expdim)
    {
    case 1:
        {
            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                const LibUtilities::PointsKey PkeyY(npointsY,LibUtilities::eFourierEvenlySpaced);
                const LibUtilities::BasisKey  BkeyY(LibUtilities::eFourier,npointsY,PkeyY);
                
                for(i = 0 ; i < nvariables; i++)
                {
                    Exp[i] = MemoryManager<MultiRegions::ContField3DHomogeneous1D>
                        ::AllocateSharedPtr(session,BkeyY,LhomY,useFFT,deal,mesh,variables[i]);
                }
            }
            else
            {
                for(i = 0 ; i < nvariables; i++)
                {
                    Exp[i] = MemoryManager<MultiRegions::ContField1D>
                        ::AllocateSharedPtr(session,mesh,variables[i]);
                }
            }
            
            break;
        }
    case 2:
        {   
            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                cout<<"homo"<<endl;
                const LibUtilities::PointsKey PkeyZ(npointsZ,LibUtilities::eFourierEvenlySpaced);
                const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourier,npointsZ,PkeyZ);
                for(i = 0 ; i < nvariables; i++)
                {                        	
                    Exp[i] = MemoryManager<MultiRegions::ContField3DHomogeneous1D>
                        ::AllocateSharedPtr(session,BkeyZ,LhomZ,useFFT,deal,mesh,variables[i]);                                    
                }
            }
            else
            {    
                //cout<<" norm field"<<endl;               	    
                i = 0;
                MultiRegions::ExpList2DSharedPtr firstfield;
                firstfield = MemoryManager<MultiRegions::ExpList2D>
                    ::AllocateSharedPtr(session,mesh,DeclareCoeffPhysArrays,variables[i]);
                
                Exp[0] = firstfield;
                for(i = 1 ; i < nvariables; i++)
                {                        	
                    Exp[i] = MemoryManager<MultiRegions::ExpList2D>
                        ::AllocateSharedPtr(*firstfield,DeclareCoeffPhysArrays);
                }
                    }
            
            break;
        }
    case 3:
        {
            if(fielddef[0]->m_numHomogeneousDir == 1)
            {
                ASSERTL0(false,"3D fully periodic problems not implemented yet");
            }
            else
            {
                i = 0;
                MultiRegions::ContField3DSharedPtr firstfield =
                    MemoryManager<MultiRegions::ContField3D>
                    ::AllocateSharedPtr(session,mesh,variables[i]);
                
                Exp[0] = firstfield;
                for(i = 1 ; i < nvariables; i++)
                {
                    Exp[i] = MemoryManager<MultiRegions::ContField3D>
                        ::AllocateSharedPtr(*firstfield,mesh,variables[i]);
                }
            }
            break;
        }
    default:
        ASSERTL0(false,"Expansion dimension not recognised");
        break;
    }   
    
}

void GenerateField(Array<OneD,MultiRegions::ExpListSharedPtr> &field0,
                   Array<OneD, NekDouble> x1,
                   Array<OneD, NekDouble> y1,
                   Array<OneD,MultiRegions::ExpListSharedPtr> &field1)
{
    Array<OneD, NekDouble> coords(2);
    int nq1 = field1[0]->GetTotPoints();
    int elmtid, offset;
    cerr << "Interpolating " << field0.num_elements() << " fields ["; 

    for(int r=0; r< nq1; r++)
    {
        if(!(r%100))
        {
            cerr << ".";
        }
        coords[0] = x1[r];
        coords[1] = y1[r];
        
        elmtid = field0[0]->GetExpIndex(coords, 0.00001);
        offset = field0[0]->GetPhys_Offset(elmtid);
 
        for(int s = 0; s < field0.num_elements(); ++s)
        {
            field1[s]->UpdatePhys()[r] = field0[s]->GetExp(elmtid)->
                PhysEvaluate(coords, field0[s]->GetPhys() +offset);    
            if( boost::math::isnan(field1[s]->UpdatePhys()[r]) )
            {            
                cout<<"x="<<x1[r]<<"   y="<<y1[r]<<"    offset="<<offset<<"  elmtid="<<elmtid<<endl;                  
                cout<<"new val="<<field1[s]->UpdatePhys()[r]<<endl;
                //ASSERTL0( abs(field1->UpdatePhys()[r])<10000000000, "interp failed");
            }
        }
    }
    cerr << "]"<< endl;
}	

void GenerateFieldHomo(MultiRegions::ExpListSharedPtr field0,
    	               Array<OneD, NekDouble> x1,
    	               Array<OneD, NekDouble> y1,
    	               MultiRegions::ExpListSharedPtr field1)
{
    Array<OneD, NekDouble> coords(2);
    int nq1 = field1->GetPlane(1)->GetTotPoints();
    ASSERTL0(nq1 == field1->GetPlane(1)->GetTotPoints(), "problem");
    int elmtid, offset;
    
    //plane 0
    for(int r=0; r< nq1; r++)
    {
        coords[0] = x1[r];
        coords[1] = y1[r];
        
        elmtid = field0->GetPlane(0)->GetExpIndex(coords, 0.00001);
        offset = field0->GetPlane(0)->GetPhys_Offset(elmtid);
        field1->GetPlane(0)->UpdatePhys()[r] = field0->GetPlane(0)->GetExp(elmtid)->
            PhysEvaluate(coords, field0->GetPlane(0)->GetPhys() +offset);    
        if( boost::math::isnan(field1->GetPlane(0)->UpdatePhys()[r]) )
        {            
            cout<<"x="<<x1[r]<<"   y="<<y1[r]<<"    offset="<<offset<<"  elmtid="<<elmtid<<endl;                  
            cout<<"new val="<<field1->GetPlane(0)->UpdatePhys()[r]<<endl;
            //ASSERTL0( abs(field1->UpdatePhys()[r])<10000000000, "interp failed");
        }
        
    } 
    
    
    //plane1
    for(int r=0; r< nq1; r++)
    {
        coords[0] = x1[r];
        coords[1] = y1[r];
        
        elmtid = field0->GetPlane(1)->GetExpIndex(coords, 0.00001);
        offset = field0->GetPlane(1)->GetPhys_Offset(elmtid);
        field1->GetPlane(1)->UpdatePhys()[r] = field0->GetPlane(1)->GetExp(elmtid)->
            PhysEvaluate(coords, field0->GetPlane(1)->GetPhys() +offset);    
        if( boost::math::isnan(field1->GetPlane(1)->UpdatePhys()[r]) )
        {            
            cout<<"x="<<x1[r]<<"   y="<<y1[r]<<"    offset="<<offset<<"  elmtid="<<elmtid<<endl;                  
            cout<<"new val="<<field1->GetPlane(1)->UpdatePhys()[r]<<endl;
            //ASSERTL0( abs(field1->UpdatePhys()[r])<10000000000, "interp failed");
        }
        
    }              
}


bool Checkbndmeshes(  Array<OneD, NekDouble> x0,
                      Array<OneD, NekDouble> y0,       	       
                      Array<OneD, NekDouble> x1,
                      Array<OneD, NekDouble> y1)
{
    NekDouble x0min,x0max,y0min, y0max;
    NekDouble x1min,x1max,y1min, y1max;       	       
    NekDouble tol = 0.0000001;
    NekDouble tol1 = 0.00001;
    x0min = Vmath::Vmin(x0.num_elements(),x0,1);
    x0max = Vmath::Vmax(x0.num_elements(),x0,1);
    y0min = Vmath::Vmin(y0.num_elements(),y0,1);
    y0max = Vmath::Vmax(y0.num_elements(),y0,1);
    
    x1min = Vmath::Vmin(x1.num_elements(),x1,1);
    x1max = Vmath::Vmax(x1.num_elements(),x1,1);
    y1min = Vmath::Vmin(y1.num_elements(),y1,1);
    y1max = Vmath::Vmax(y1.num_elements(),y1,1);       	       
    
    //cout<<std::setprecision(8)<<"x0max="<<x0max<<endl;
    //cout<<std::setprecision(8)<<"x1max="<<x1max<<endl;
    //cout<<std::setprecision(8)<<"y0max="<<y0max<<endl;
    //cout<<std::setprecision(8)<<"y1max="<<y1max<<endl;
    //cout<<"abs(x0min-x1min )="<<abs(x0min-x1min )<<endl;
    //cout<<"abs(x0max-x1max)="<<abs(x0max-x1max)<<endl;
    //cout<<"abs(y0min-y1min)="<<abs(y0min-y1min)<<endl;
    //cout<<"abs(y0max-y1max)="<<abs(y0max-y1max)<<endl;
    
    if(  abs(x0min-x1min )< tol1
         && abs(x0max-x1max)< tol1
         && abs(y0min-y1min)< tol1
         && abs(y0max-y1max)< tol1
         )
    {
        if(abs(x0min-x1min )> tol
           || abs(x0max-x1max)> tol
           || abs(y0min-y1min)> tol
           || abs(y0max-y1max)> tol  )
        {
            cout<<"Warning: mesh boundary points differ more than 10^-7"<<endl;
        }

        return true;                  
    }
    else
    { 
               	   return false;
    }
    

}       	       


void Writefield(LibUtilities::SessionReaderSharedPtr vSession,
                const vector<std::string> &variables,
                string fieldfile, 
                SpatialDomains::MeshGraphSharedPtr &graph,  	    
                Array<OneD, MultiRegions::ExpListSharedPtr> &outfield)
{
    string var;
    std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
        = outfield[0]->GetFieldDefinitions();  			
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());    		
    Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(outfield.num_elements());   	
    
    for(int j=0; j< fieldcoeffs.num_elements(); ++j)
    {  
        if(FieldDef[0]->m_numHomogeneousDir)
        {
            //plane 0
            outfield[j]->GetPlane(0)->FwdTrans_IterPerExp(outfield[j]->GetPlane(0)->GetPhys(),outfield[j]->GetPlane(0)->UpdateCoeffs());
            
            //plane 1
            outfield[j]->GetPlane(1)->FwdTrans_IterPerExp(outfield[j]->GetPlane(1)->GetPhys(),outfield[j]->GetPlane(1)->UpdateCoeffs());
            
        }
        else
        {  
            outfield[j]->FwdTrans_IterPerExp(outfield[j]->GetPhys(),outfield[j]->UpdateCoeffs());
                     }
 	
        fieldcoeffs[j] = outfield[j]->UpdateCoeffs();	
        
        for(int i=0; i< FieldDef.size(); i++)
        {		     	     
            FieldDef[i]->m_fields.push_back(variables[j]);   
            outfield[0]->AppendFieldData(FieldDef[i], FieldData[i], fieldcoeffs[j]);  
        }
    }
    graph->Write(fieldfile,FieldDef,FieldData);		
}    		
