#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>


#include <MultiRegions/ExpList.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>

#include <tinyxml.h>

#include <boost/math/special_functions/fpclassify.hpp>

using namespace Nektar;

int main(int argc, char *argv[])
{
    void SetFields(
        SpatialDomains::MeshGraphSharedPtr              &mesh,
        vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef,
        LibUtilities::SessionReaderSharedPtr            &session,
        Array<OneD,MultiRegions::ExpListSharedPtr>      &Exp, 
        int                                             nvariables,
        const vector<std::string>                       &variables, 
        bool                                            homogeneous);
    
    void InterpolateField(
        Array<OneD, MultiRegions::ExpListSharedPtr> &field0,
        Array<OneD, MultiRegions::ExpListSharedPtr> &field1,
        Array<OneD, NekDouble>                      x1,
        Array<OneD, NekDouble>                      y1,
        Array<OneD, NekDouble>                      z1 = NullNekDouble1DArray,
        NekDouble                                   clamp_low = -10000000,
        NekDouble                                   clamp_up = 10000000);
    
    void InterpolateFieldHomo(
        MultiRegions::ExpListSharedPtr              field0,
        Array<OneD, NekDouble>                      x1,
        Array<OneD, NekDouble>                      y1,
        MultiRegions::ExpListSharedPtr              field1);
    
    bool Checkbndmeshes2D(    
        Array<OneD, NekDouble>                      x0,
        Array<OneD, NekDouble>                      y0,       	       
        Array<OneD, NekDouble>                      x1,
        Array<OneD, NekDouble>                      y1);
    
    bool Checkbndmeshes3D(    
        Array<OneD, NekDouble>                      x0,
        Array<OneD, NekDouble>                      y0,       	       
        Array<OneD, NekDouble>                      z0,       	       
        Array<OneD, NekDouble>                      x1,
        Array<OneD, NekDouble>                      y1,       	       
        Array<OneD, NekDouble>                      z1);
    
    void Writefield(
        LibUtilities::SessionReaderSharedPtr        vSession,
        const vector<std::string>                   &variables,
        string                                      fieldfile, 
        SpatialDomains::MeshGraphSharedPtr          &graph,  	    
        Array<OneD, MultiRegions::ExpListSharedPtr> &outfield);    

    if (argc != 5)
    {
        fprintf(stderr, "Usage: ./FieldToField  meshfile0 fieldfile0  "
                        "meshfile1  fieldfile1\n");
        exit(1);
    }
    
    //----------------------------------------------
    string meshfile0(argv[argc-4]); 
    string fieldfile0(argv[argc-3]);
    
    std::vector<std::string> filenames0;
    filenames0.push_back(meshfile0);
    filenames0.push_back(meshfile0);
    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv, filenames0);
  
    // Read in mesh from input file0
    SpatialDomains::MeshGraphSharedPtr graphShPt 
        = SpatialDomains::MeshGraph::Read(meshfile0);
    int expdim = graphShPt->GetMeshDimension();
    
    //----------------------------------------------          
    // Import fieldfile0.
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    LibUtilities::Import(fieldfile0, fielddef, fielddata);
    //----------------------------------------------    

    //read info from fldfile
   // const std::vector<std::string> variables = fielddef[0]->m_fields; 
    const std::vector<std::string> variables = vSession->GetVariables(); 
    bool homo = (fielddef[0]->m_numHomogeneousDir > 0)? true: false;
    
    // Define Expansion    
    int nfields; 
    nfields = variables.size();
    Array<OneD, MultiRegions::ExpListSharedPtr> fields; 
    fields= Array<OneD, MultiRegions::ExpListSharedPtr>(nfields);  
    SetFields(graphShPt,fielddef, vSession,fields,nfields,variables,homo);
    int nq; //pointsper plane
    
    //-----------------------------------------------
    // Copy data from file:fill fields with the fielddata
    if (fielddef[0]->m_numHomogeneousDir == 1)
    {
        nq = fields[0]->GetPlane(1)->GetTotPoints();
        
        //THE IM PHYS VALUES ARE WRONG USING bwdTrans !!!
        for (int j = 0; j < nfields; ++j)
        {  
            for (int i = 0; i < fielddata.size(); ++i)
            {
                fields[j]->ExtractDataToCoeffs(fielddef[i], 
                                               fielddata[i], 
                                               fielddef[i]->m_fields[j], 
                                               fields[j]->UpdateCoeffs());
            }
            
            //bwd plane 0
            fields[j]->GetPlane(0)->BwdTrans_IterPerExp(
                                        fields[j]->GetPlane(0)->GetCoeffs(),  
                                        fields[j]->GetPlane(0)->UpdatePhys());
            
            //bwd plane 1
            fields[j]->GetPlane(1)->BwdTrans_IterPerExp(
                                        fields[j]->GetPlane(1)->GetCoeffs(), 
                                        fields[j]->GetPlane(1)->UpdatePhys());
        } 
    }
    else
    {
        nq = fields[0]->GetTotPoints();
        for (int j = 0; j < nfields; ++j)
        {  	    
            for (int i = 0; i < fielddata.size(); ++i)
            {
                fields[j]->ExtractDataToCoeffs(fielddef[i], 
                                               fielddata[i], 
                                               fielddef[i]->m_fields[j], 
                                               fields[j]->UpdateCoeffs());
            }             
            fields[j]->BwdTrans_IterPerExp(fields[j]->GetCoeffs(), 
                                           fields[j]->UpdatePhys());      
        }
    }
    

    // store mesh0 quadrature points    
    Array<OneD, NekDouble> x0(nq);
    Array<OneD, NekDouble> y0(nq);  
    Array<OneD, NekDouble> z0(nq);
    
    if (fielddef[0]->m_numHomogeneousDir == 1)
    {
        fields[0]->GetPlane(1)->GetCoords(x0, y0, z0);    
    }
    else
    {
        if (expdim == 2)
        {
            fields[0]->GetCoords(x0, y0);  
        }
        else if (expdim == 3)
        {
            fields[0]->GetCoords(x0, y0, z0);
        }
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
    //the homo quantities only for the mesh0 case
    
    graphShPt1 = SpatialDomains::MeshGraph::Read(meshfile1);  
    
    //define second session1..	
    std::vector<std::string> filenames1;       
    filenames1.push_back(meshfile1);   
    filenames1.push_back(meshfile1);     

    argc = 2;
    argv[1] = argv[3];
    argv[2] = argv[4];

    LibUtilities::SessionReaderSharedPtr vSession1
        = LibUtilities::SessionReader::CreateInstance(argc, argv, 
                                                      filenames1, 
                                                      vSession->GetComm());
    //--------------------------------------------------
    
    //set explist 1
    SetFields(graphShPt1, fielddef, vSession1, outfield, nfields, 
              variables, homo);
    //------------------------------------------------ 
    // store the new points:
    int nq1;
    if (fielddef[0]->m_numHomogeneousDir == 1)
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
    
    if (fielddef[0]->m_numHomogeneousDir == 1)
    {
        outfield[0]->GetPlane(1)->GetCoords(x1, y1, z1);    
    }
    else
    {
        if (expdim == 2)
        {
            outfield[0]->GetCoords(x1, y1);
        }
        else if (expdim == 3)
        {
            outfield[0]->GetCoords(x1, y1, z1);
        }
    }

    //------------------------------------------------
    //check 2Dmeshes compatibilities
    // same max min values x,y...
    //bool check = false;
    bool check = true;
    if (expdim == 2)
    {
        //check = Checkbndmeshes2D(x0, y0, x1, y1);
    }
    else if (expdim == 3)
    {
        check = Checkbndmeshes3D(x0, y0, z0, x1, y1, z1);
    }
    ASSERTL0(check, "meshes not compatible (different borders)");   
    //-----------------------------------------------
    
    if (fielddef[0]->m_numHomogeneousDir == 1)
    {
        for (int t = 0; t < nfields; t++)
        {        
    	    InterpolateFieldHomo(fields[t], x1, y1, outfield[t]);
        }
    }
    else
    {
        cout << "Interpolating [" << flush;
                                
        if (expdim == 2)
        {
            InterpolateField(fields, outfield, x1, y1);
        }
        else if (expdim == 3)
        {
            NekDouble clamp_up, clamp_low;
            vSession->LoadParameter("ClampToUpperValue", clamp_up,   10000000);
            vSession->LoadParameter("ClampToLowerValue", clamp_low, -10000000);
            InterpolateField(fields, outfield, x1, y1, z1, clamp_low, clamp_up);
        }
        cout << "]" << endl;
    }


    //------------------------------------------------
    //write fieldfile
    Writefield(vSession, variables, fieldfile1, graphShPt1, outfield);        
    //------------------------------------------------

}



// Define Expansion       		
void SetFields(
    SpatialDomains::MeshGraphSharedPtr              &mesh,
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef,
    LibUtilities::SessionReaderSharedPtr            &session,
    Array<OneD,MultiRegions::ExpListSharedPtr>      &Exp,
    int                                             nvariables,
    const vector<std::string>                       &variables, 
    bool                                            homogeneous)
{
    // Session reader stuff has to be evaluated only from the
    // first session which refers to mesh0
    static int cnt = 0;
    
    // Setting parameteres for homogenous problems
    NekDouble static LhomY;///< physical length in Y direction (if homogeneous)
    NekDouble static LhomZ;///< physical length in Z direction (if homogeneous)
		
    bool static DeclareCoeffPhysArrays = true;
	
    int static npointsY;///< number of points in Y direction (if homogeneous)
    int static npointsZ;///< number of points in Z direction (if homogeneous)	
    
    bool static useFFT = false;
    bool static deal   = false;
    
    cnt++;
    int i;		
    int expdim = mesh->GetMeshDimension();
    
    switch (expdim)
    {
        case 1:
        {
            if (fielddef[0]->m_numHomogeneousDir == 1)
            {
                const LibUtilities::PointsKey PkeyY(
                                        npointsY, 
                                        LibUtilities::eFourierEvenlySpaced);
                const LibUtilities::BasisKey  BkeyY(
                                        LibUtilities::eFourier, 
                                        npointsY, PkeyY);
                
                for (i = 0 ; i < nvariables; i++)
                {
                    Exp[i] = MemoryManager<MultiRegions::
                        ContField3DHomogeneous1D>::AllocateSharedPtr(
                                                    session, BkeyY, LhomY, 
                                                    useFFT, deal, mesh, 
                                                    variables[i]);
                }
            }
            else
            {
                for (i = 0 ; i < nvariables; i++)
                {
                    Exp[i] = MemoryManager<MultiRegions::ContField1D>
                        ::AllocateSharedPtr(session, mesh, variables[i]);
                }
            }
            break;
        }
        case 2:
        {   
            if (fielddef[0]->m_numHomogeneousDir == 1)
            {
                const LibUtilities::PointsKey PkeyZ(
                                        npointsZ, 
                                        LibUtilities::eFourierEvenlySpaced);
                const LibUtilities::BasisKey  BkeyZ(
                                        LibUtilities::eFourier,
                                        npointsZ, PkeyZ);
                
                MultiRegions::ContField3DHomogeneous1DSharedPtr firstfield = 
                    MemoryManager<MultiRegions::ContField3DHomogeneous1D>::
                        AllocateSharedPtr(session, BkeyZ, LhomZ, useFFT, 
                                          deal, mesh, variables[0]);
                
                Exp[0] = firstfield;
                for (i = 1 ; i < nvariables; i++)
                {                        	
                    Exp[i] = MemoryManager<MultiRegions::
                        ContField3DHomogeneous1D>::
                            AllocateSharedPtr(*firstfield);                                   
                }
            }
            else
            {    
                i = 0;
                MultiRegions::DisContField2DSharedPtr firstfield;
                firstfield = MemoryManager<MultiRegions::DisContField2D>::
                    AllocateSharedPtr(session, mesh, variables[i],
                                      DeclareCoeffPhysArrays);
                
                Exp[0] = firstfield;
                for (i = 1 ; i < nvariables; i++)
                {                        	
                    Exp[i] = MemoryManager<MultiRegions::DisContField2D>
                        ::AllocateSharedPtr(*firstfield, mesh, variables[i],
                                            DeclareCoeffPhysArrays);
                }
            }
            break;
        }
        case 3:
        {
            if (fielddef[0]->m_numHomogeneousDir == 1)
            {
                ASSERTL0(false, 
                         "3D fully periodic problems not implemented yet");
            }
            else
            {
                i = 0;
                MultiRegions::ExpList3DSharedPtr firstfield =
                    MemoryManager<MultiRegions::ExpList3D>::
                        AllocateSharedPtr(session, mesh);
                
                Exp[0] = firstfield;
                for (i = 1 ; i < nvariables; i++)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3D>
                        ::AllocateSharedPtr(*firstfield);
                }
            }
            break;
        }
        default:
            ASSERTL0(false, "Expansion dimension not recognised");
            break;
    }   
    
}



void InterpolateField(
    Array<OneD, MultiRegions::ExpListSharedPtr> &field0,
    Array<OneD, MultiRegions::ExpListSharedPtr> &field1,
    Array<OneD, NekDouble>                      x,
    Array<OneD, NekDouble>                      y,
    Array<OneD, NekDouble>                      z,
    NekDouble                                   clamp_low,
    NekDouble                                   clamp_up)
{
    int expdim = (z == NullNekDouble1DArray)? 2: 3;
    
    Array<OneD, NekDouble> coords(expdim), Lcoords(expdim);
    int nq1 = field1[0]->GetTotPoints();
    int elmtid, offset;
    int r, f;
    static int intpts = 0;
    
    ASSERTL0(field0.num_elements() == field1.num_elements(), 
             "Input field dimension must be same as output dimension");

    for (r = 0; r < nq1; r++)
    {
        coords[0] = x[r];
        coords[1] = y[r];
        if (expdim == 3)
        {
            coords[2] = z[r];
        }

        // Obtain Element and LocalCoordinate to interpolate
        elmtid = field0[0]->GetExpIndex(coords, Lcoords, 1e-3);

        offset = field0[0]->GetPhys_Offset(field0[0]->
                                           GetOffset_Elmt_Id(elmtid));

        for (f = 0; f < field1.num_elements(); ++f)
        {
            NekDouble value;
            value = field0[f]->GetExp(elmtid)->
                StdPhysEvaluate(Lcoords, field0[f]->GetPhys() +offset);    

            if ((boost::math::isnan)(value))
            {            
                ASSERTL0(false, "new value is not a number");
            }
            else
            {
                value = (value > clamp_up)? clamp_up : 
                        ((value < clamp_low)? clamp_low :
                        value);
                
                field1[f]->UpdatePhys()[r] = value;
            }
        }
        
        if (intpts%1000 == 0)
        {
            cout <<"." << flush;
        }
        intpts ++;
    }    
}



void InterpolateFieldHomo(
    MultiRegions::ExpListSharedPtr field0,
    Array<OneD, NekDouble>         x1,
    Array<OneD, NekDouble>         y1,
    MultiRegions::ExpListSharedPtr field1)
{
    Array<OneD, NekDouble> coords(2);
    int nq1 = field1->GetPlane(1)->GetTotPoints();
    int elmtid, offset;
    
    //plane 0
    for (int r = 0; r < nq1; r++)
    {
        coords[0] = x1[r];
        coords[1] = y1[r];
        
        elmtid = field0->GetPlane(0)->GetExpIndex(coords, 1e-3);
        offset = field0->GetPlane(0)->GetPhys_Offset(elmtid);
        field1->GetPlane(0)->UpdatePhys()[r] = field0->GetPlane(0)->
            GetExp(elmtid)->PhysEvaluate(
                                    coords, 
                                    field0->GetPlane(0)->GetPhys() +offset);
        
        if ((boost::math::isnan)(field1->GetPlane(0)->UpdatePhys()[r]))
        {            
            //ASSERTL0(abs(field1->UpdatePhys()[r])<10000000000, 
            //         "interp failed");
        }
    } 
    
    
    // Plane1
    for (int r = 0; r < nq1; r++)
    {
        coords[0] = x1[r];
        coords[1] = y1[r];
        
        elmtid = field0->GetPlane(1)->GetExpIndex(coords, 1e-3);
        offset = field0->GetPlane(1)->GetPhys_Offset(elmtid);
        field1->GetPlane(1)->UpdatePhys()[r] = field0->GetPlane(1)->
            GetExp(elmtid)->PhysEvaluate(
                                    coords, 
                                    field0->GetPlane(1)->GetPhys() +offset);
        
        if((boost::math::isnan)(field1->GetPlane(1)->UpdatePhys()[r]))
        {            
            //ASSERTL0(abs(field1->UpdatePhys()[r])<10000000000, 
            //         "interp failed");
        }
        
    }              
}



bool Checkbndmeshes2D(  
    Array<OneD, NekDouble> x0,
    Array<OneD, NekDouble> y0,       	       
    Array<OneD, NekDouble> x1,
    Array<OneD, NekDouble> y1)
{
    NekDouble x0min, x0max, y0min, y0max;
    NekDouble x1min, x1max, y1min, y1max;       	       
    NekDouble tol  = 0.0000001;
    NekDouble tol1 = 0.00001;
    x0min = Vmath::Vmin(x0.num_elements(), x0, 1);
    x0max = Vmath::Vmax(x0.num_elements(), x0, 1);
    y0min = Vmath::Vmin(y0.num_elements(), y0, 1);
    y0max = Vmath::Vmax(y0.num_elements(), y0, 1);
    
    x1min = Vmath::Vmin(x1.num_elements(), x1, 1);
    x1max = Vmath::Vmax(x1.num_elements(), x1, 1);
    y1min = Vmath::Vmin(y1.num_elements(), y1, 1);
    y1max = Vmath::Vmax(y1.num_elements(), y1, 1);       	       
    
    if (abs(x0min-x1min) < tol1 &&
        abs(x0max-x1max) < tol1 &&
        abs(y0min-y1min) < tol1 &&
        abs(y0max-y1max) < tol1)
    {
        if (abs(x0min-x1min) > tol ||
            abs(x0max-x1max) > tol ||
            abs(y0min-y1min) > tol ||
            abs(y0max-y1max) > tol)
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



bool Checkbndmeshes3D(
    Array<OneD, NekDouble> x0,
    Array<OneD, NekDouble> y0,
    Array<OneD, NekDouble> z0,       	       
    Array<OneD, NekDouble> x1,
    Array<OneD, NekDouble> y1,
    Array<OneD, NekDouble> z1)
{
    NekDouble x0min, x0max, y0min, y0max, z0min, z0max;
    NekDouble x1min, x1max, y1min, y1max, z1min, z1max;       	       
    NekDouble tol  = 0.0000001;
    NekDouble tol1 = 0.00001;
    x0min = Vmath::Vmin(x0.num_elements(), x0, 1);
    x0max = Vmath::Vmax(x0.num_elements(), x0, 1);
    y0min = Vmath::Vmin(y0.num_elements(), y0, 1);
    y0max = Vmath::Vmax(y0.num_elements(), y0, 1);
    z0min = Vmath::Vmin(z0.num_elements(), z0, 1);
    z0max = Vmath::Vmax(z0.num_elements(), z0, 1);

    x1min = Vmath::Vmin(x1.num_elements(), x1, 1);
    x1max = Vmath::Vmax(x1.num_elements(), x1, 1);
    y1min = Vmath::Vmin(y1.num_elements(), y1, 1);
    y1max = Vmath::Vmax(y1.num_elements(), y1, 1);
    z1min = Vmath::Vmin(z1.num_elements(), z1, 1);
    z1max = Vmath::Vmax(z1.num_elements(), z1, 1);       	       
    
    if (abs(x0min-x1min) < tol1 &&
        abs(x0max-x1max) < tol1 &&
        abs(y0min-y1min) < tol1 &&
        abs(y0max-y1max) < tol1 &&
        abs(z0min-z1min) < tol1 &&
        abs(z0max-z1max) < tol1)
    {
        if (abs(x0min-x1min) > tol ||
            abs(x0max-x1max) > tol ||
            abs(y0min-y1min) > tol ||
            abs(y0max-y1max) > tol ||
            abs(z0min-z1min) > tol ||
            abs(z0max-z1max) > tol)
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



void Writefield(
    LibUtilities::SessionReaderSharedPtr        vSession,
    const std::vector<std::string>              &variables,
    string                                      fieldfile, 
    SpatialDomains::MeshGraphSharedPtr          &graph,  	    
    Array<OneD, MultiRegions::ExpListSharedPtr> &outfield)
{
    string var;
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = outfield[0]->GetFieldDefinitions();  			
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());    		
    Array<OneD, Array<OneD, NekDouble> > fieldcoeffs(outfield.num_elements());   	
    
    for (int j = 0; j < fieldcoeffs.num_elements(); ++j)
    {  
        if (FieldDef[0]->m_numHomogeneousDir == 1)
        {
            // plane 0
            outfield[j]->GetPlane(0)->FwdTrans_IterPerExp(
                                    outfield[j]->GetPlane(0)->GetPhys(),
                                    outfield[j]->GetPlane(0)->UpdateCoeffs());
            
            // plane 1
            outfield[j]->GetPlane(1)->FwdTrans_IterPerExp(
                                    outfield[j]->GetPlane(1)->GetPhys(),
                                    outfield[j]->GetPlane(1)->UpdateCoeffs());
        }
        else
        {  
            outfield[j]->FwdTrans_IterPerExp(
                                    outfield[j]->GetPhys(), 
                                    outfield[j]->UpdateCoeffs());
        }
 	
        fieldcoeffs[j] = outfield[j]->UpdateCoeffs();	
        
        for (int i = 0; i < FieldDef.size(); i++)
        {		     	     
            FieldDef[i]->m_fields.push_back(variables[j]);   
            outfield[0]->AppendFieldData(FieldDef[i], FieldData[i], 
                                         fieldcoeffs[j]);  
        }
    }
    LibUtilities::Write(fieldfile, FieldDef, FieldData);		
}    		
