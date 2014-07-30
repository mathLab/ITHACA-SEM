#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField2D.h>
#include <SpatialDomains/MeshGraph2D.h>

#include <Collections/Collection.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ContField2DSharedPtr Exp;

    int Ntest = 100;
    if(argc < 2)
    {
        fprintf(stderr,"Usage: Collection2D ntest meshfile  or Collection2D meshfile \n");
        exit(1);
    }
    
    // Read in mesh and set up ExpList
    SpatialDomains::MeshGraphSharedPtr graph2D = 
        SpatialDomains::MeshGraph::Read(vSession);

    Exp = MemoryManager<MultiRegions::ContField2D>::
        AllocateSharedPtr(vSession,graph2D,vSession->GetVariable(0));

    int nelmt = Exp->GetNumElmts();

    for(int imp = 0; imp < 3; ++imp)
    {

        // set up different collection implementations:
        switch(imp){
        case 2:
            {
                Exp->CreateCollections(Collections::eSumFac);
                cout << endl << "Using SumFac Collection Implementation" << endl;
            }
            break;
        case 1:
            {
                Exp->CreateCollections(Collections::eStdMat);
                cout << endl << "Using StdMat Collection Implementation" << endl;
            }
            break;
        default:
            {
                cout <<"Using IterPerExp Collection Implementation" << endl;
            }
            break;
        }
        
        //BwdTrans comparison
        {
            
            cout << "BwdTrans Op: Ntest = " << Ntest << endl;
            
            Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(Exp->GetNpoints());
            Array<OneD, NekDouble> phys2(Exp->GetNpoints());
            
            Timer t;
            t.Start();
            
            // Do test by calling every element in loop 
            for (int i = 0; i < Ntest; ++i)
            {
                for(int j = 0; j < nelmt; ++j)
                {
                    Exp->GetExp(j)->BwdTrans(coeffs+Exp->GetCoeff_Offset(j), tmp = phys1+Exp->GetPhys_Offset(j));
                }
            }
            t.Stop();
            NekDouble orig = t.TimePerTest(Ntest);
            cout << "\t ExpList   : " << orig << endl; 
            
            t.Start();
            // call collection implementation in thorugh ExpList. 
            for (int i = 0; i < Ntest; ++i)
            {
                Exp->BwdTrans(coeffs, phys2);
            }
            t.Stop();
            NekDouble col = t.TimePerTest(Ntest);
            cout << "\t Collection: " << t.TimePerTest(Ntest) << endl;
            Vmath::Vsub(phys1.num_elements(), phys1, 1, phys2, 1, phys1, 1);
            cout << "\t Difference: "<< Vmath::Vmax(phys1.num_elements(), phys1, 1) << endl;
            
            cout << "\t Ratio: " << (orig/col) << endl;
            
        }
        
        // IProductWRTBase Comparison 
        {
            cout << "IProductWRTBase Op: Ntest = " << Ntest << endl;
            
            const int nq = Exp->GetNpoints();
            const int nc = Exp->GetNcoeffs();
            Array<OneD, NekDouble> xc(nq), yc(nq), tmp,tmp1;
            Array<OneD, NekDouble> input(nq), output1(nc), output2(nc);
            
            Exp->GetCoords(xc, yc);
            for (int i = 0; i < nq; ++i)
            {
                input[i] = sin(xc[i])*cos(yc[i]);
            }
            
            Timer t;
            t.Start();
            for (int i = 0; i < Ntest; ++i)
            {
                // Do test by calling every element in loop 
                for (int i = 0; i < Ntest; ++i)
                {
                    for(int j = 0; j < nelmt; ++j)
                    {
                        Exp->GetExp(j)->IProductWRTBase(input + Exp->GetPhys_Offset(j), 
                                                        tmp  = output1 + Exp->GetCoeff_Offset(j));
                    }
                }
            }
            t.Stop();
            NekDouble orig = t.TimePerTest(Ntest);
            cout << "\t ExpList: " << orig << endl; 
            t.Start();
            Exp->IProductWRTBase(input, output2);
            t.Stop();
            NekDouble col = t.TimePerTest(Ntest);
            cout << "\t Collection: " << t.TimePerTest(Ntest) << endl;
            Vmath::Vsub(nc, output1, 1, output2, 1, output1, 1);
            cout << "\t Difference: " << Vmath::Vmax(nc, output1, 1) << endl;
            cout << "\t Ratio: " << (orig/col) << endl;
            
        }

        // PhysDeriv Comparison 
        {
            cout << "PhysDeriv Op: Ntest = " << Ntest << endl;
            
            const int nq = Exp->GetNpoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), tmp,tmp1;
            Array<OneD, NekDouble> input(nq), diff1_0(nq), diff1_1(nq);
            Array<OneD, NekDouble> diff2_0(nq), diff2_1(nq);
            
            Exp->GetCoords(xc, yc);
            
            for (int i = 0; i < nq; ++i)
            {
                input[i] = sin(xc[i])*cos(yc[i]);
            }
            
            Timer t;
            t.Start();
            for (int i = 0; i < Ntest; ++i)
            {
                // Do test by calling every element in loop 
                for (int i = 0; i < Ntest; ++i)
                {
                    for(int j = 0; j < nelmt; ++j)
                    {
                        Exp->GetExp(j)->PhysDeriv(input +Exp->GetPhys_Offset(j), 
                                                  tmp  = diff1_0 +Exp->GetPhys_Offset(j), 
                                                  tmp1 = diff1_1 +Exp->GetPhys_Offset(j));
                    }
                }
            }
            t.Stop();
            NekDouble orig = t.TimePerTest(Ntest);
            cout << "\t ExpList: " << orig << endl; 
            t.Start();
            Exp->PhysDeriv(input, diff2_0, diff2_1);
            t.Stop();
            NekDouble col = t.TimePerTest(Ntest);
            cout << "\t Collection: " << t.TimePerTest(Ntest) << endl;
            Vmath::Vsub(nq, diff1_0, 1, diff2_0, 1, diff1_0, 1);
            Vmath::Vsub(nq, diff1_1, 1, diff2_1, 1, diff1_1, 1);
            cout << "\t Difference: " << Vmath::Vmax(nq, diff1_0, 1) << "," << Vmath::Vmax(nq, diff1_0, 1) << endl;
            cout << "\t Ratio: " << (orig/col) << endl;
            
        }
    }
    vSession->Finalise();

    return 0;
}

