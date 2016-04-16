/**
 * This function calculate the vorticity vector starting from an .fld file.
 * It is meant to be used with solutions produced by the incompressible Navier-Stokes solver.
 * To use it with solutions coming form another solver further generalisations are required.
 */
#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>


#include <SolverUtils/EquationSystem.h>

#include <LocalRegions/MatrixKey.h>
#include <LibUtilities/BasicUtils/Equation.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>

using namespace std;
using namespace Nektar;

int main(int argc, char *argv[])
{
    int i,j;
    
    if(argc != 3)
    {
        fprintf(stderr,"Usage: ./StreamFunction2D file.xml file.fld\n");
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
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    LibUtilities::Import(fieldfile,fielddef,fielddata);
    
    //----------------------------------------------
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_fields.size();
    int vorticitydim;

    vorticitydim = 2;
    
    
    
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields + vorticitydim);
    
    Array<OneD, MultiRegions::ExpListSharedPtr> FieldVar(1);
    
    switch(expdim)
    {

        case 2:
        {

            {
                MultiRegions::ExpList2DSharedPtr Exp2D;
                Exp2D = MemoryManager<MultiRegions::ExpList2D>
                ::AllocateSharedPtr(vSession,graphShPt);
                Exp[0] =  Exp2D;
                
                FieldVar[0] = MemoryManager<MultiRegions::ContField2D>::AllocateSharedPtr(vSession, graphShPt, "v");
                                
                for(i = 1; i < nfields + vorticitydim; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2D>
                    ::AllocateSharedPtr(*Exp2D);
                }
            }
        }
        break;

        default:
            ASSERTL0(false,"The input file must be two-dimensional");
            break;
    }
    
    //----------------------------------------------
    // Copy data from field file
    for(j = 0; j < nfields; ++j)
    {
        for(unsigned int i = 0; i < fielddata.size(); ++i)
        {
            Exp[j]->ExtractDataToCoeffs(fielddef [i],
                                        fielddata[i],
                                        fielddef [i]->m_fields[j],
                                        Exp[j]->UpdateCoeffs());
        }
        Exp[j]->BwdTrans(Exp[j]->GetCoeffs(),Exp[j]->UpdatePhys());
    }
    //----------------------------------------------
    
    int nq = Exp[0]->GetNpoints();
    
    Array<OneD, NekDouble> Uy(nq);
    Array<OneD, NekDouble> Vx(nq);

    Array<OneD, NekDouble> StreamFunc(nq);
    Array<OneD, NekDouble> Qz(nq);
    
    switch(expdim)
    {

        case 2:
        {

            {
                Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[0],Exp[1]->GetPhys(),Vx);
                Exp[0]->PhysDeriv(MultiRegions::DirCartesianMap[1],Exp[0]->GetPhys(),Uy);
                
                Vmath::Vsub(nq,Vx,1,Uy,1,Qz,1);
                
                //The Vorticity is stored
                Exp[3]->FwdTrans(Qz, Exp[3]->UpdateCoeffs());
                
                //We now calculate the Stream Function as the solution of the
                //Poisson equation: Vorticity = - \nabla^2 StreamFunction
                StdRegions::ConstFactorMap factor;
                factor[StdRegions::eFactorLambda] = 0.0;
                
                Vmath::Smul(nq,-1.0,Qz,1,Qz,1);
                Exp[4]->SetPhys(Qz);
                
                FieldVar[0]->HelmSolve(Qz, Exp[4]->UpdateCoeffs(), NullFlagList, factor);
            }
        }
        break;

        default:
        {
            ASSERTL0(false,"The input file must be two-dimensional");
        }
        break;
    }
    
    //-----------------------------------------------
    // Write solution to file with additional computed fields
    string   fldfilename(argv[2]);
    string   out = fldfilename.substr(0, fldfilename.find_last_of("."));
    string   endfile("_with_2DstremFunction.fld");
    out += endfile;
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
    = Exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());
    
    for(j = 0; j < nfields + vorticitydim; ++j)
    {
        for(i = 0; i < FieldDef.size(); ++i)
        {
            if (j >= nfields)
            {
                if(j == 4)
                {
                    FieldDef[i]->m_fields.push_back("StreamFunc");
                }
                else
                {
                    FieldDef[i]->m_fields.push_back("Qz");
                }
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

