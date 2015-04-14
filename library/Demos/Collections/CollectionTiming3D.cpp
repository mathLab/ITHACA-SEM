#include <cstdio>
#include <cstdlib>
#include <iomanip>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ExpList3D.h>
#include <SpatialDomains/MeshGraph3D.h>

#include <Collections/Collection.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ExpList3DSharedPtr Exp;

    if(argc < 2)
    {
        fprintf(stderr,"Usage: Collection3D meshfile \n");
        exit(1);
    }

    int Ntest;
    vSession->LoadParameter("Ntest", Ntest, 1000);

    // Read in mesh and set up ExpList
    SpatialDomains::MeshGraphSharedPtr graph3D =
        SpatialDomains::MeshGraph::Read(vSession);

    Collections::ImplementationType ImpType;
    for(int imp = 0; imp < 4; ++imp)
    {
        // set up different collection implementations:
        switch(imp){
        case 3:
            {
                ImpType = Collections::eSumFac;
                cout << endl << "Using SumFac Collection Implementation:" << endl;
            }
            break;
        case 2:
            {
                ImpType = Collections::eStdMat;
                cout << endl << "Using StdMat Collection Implementation:" << endl;
            }
            break;
        case 1:
            {
                ImpType = Collections::eIterPerExp;
                cout <<"Using IterPerExp Collection Implementation:" << endl;
            }
            break;
        default:
            {
                ImpType = Collections::eNoCollection;
                cout <<"Using No Collection:" << endl;
            }
            break;
        }

        //BwdTrans comparison
        cout << "BwdTrans Op: Ntest = " << Ntest << endl;

        for(int N = 2; N < 11; ++N)
        {
            graph3D->SetExpansionsToPolyOrder(N);

            Exp = MemoryManager<MultiRegions::ExpList3D>::
                AllocateSharedPtr(vSession,graph3D);

            Exp->CreateCollections(ImpType);

            int nelmt = Exp->GetNumElmts();

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
            cout << "N: "<< std::setw(2) << N << " ExpList: " <<std::setw(11) <<  orig;
            t.Start();
            // call collection implementation in thorugh ExpList.
            for (int i = 0; i < Ntest; ++i)
            {
                Exp->BwdTrans(coeffs, phys2);
            }
            t.Stop();
            NekDouble col = t.TimePerTest(Ntest);
            cout << "   Collection: " << std::setw(11) << t.TimePerTest(Ntest);
            cout << "   Ratio: " << (orig/col) << endl;
        }


        // IProductWRTBase Comparison
        cout << "IProductWRTBase Op: Ntest = " << Ntest << endl;

        for(int N = 2; N < 11; ++N)
        {
            graph3D->SetExpansionsToPolyOrder(N);

            Exp = MemoryManager<MultiRegions::ExpList3D>::
                AllocateSharedPtr(vSession,graph3D);

            Exp->CreateCollections(ImpType);

            int nelmt = Exp->GetNumElmts();

            const int nq = Exp->GetNpoints();
            const int nc = Exp->GetNcoeffs();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq), tmp,tmp1;
            Array<OneD, NekDouble> input(nq), output1(nc), output2(nc);

            Exp->GetCoords(xc, yc, zc);
            for (int i = 0; i < nq; ++i)
            {
                input[i] = sin(xc[i])*cos(yc[i])*cos(zc[i]);
            }

            Timer t;
            t.Start();
            for (int i = 0; i < Ntest; ++i)
            {
                // Do test by calling every element in loop
                for(int j = 0; j < nelmt; ++j)
                {
                    Exp->GetExp(j)->IProductWRTBase(input + Exp->GetPhys_Offset(j),
                                                    tmp  = output1 + Exp->GetCoeff_Offset(j));
                }
            }
            t.Stop();
            NekDouble orig = t.TimePerTest(Ntest);
            cout << "N: "<< std::setw(2) << N << " ExpList: " <<std::setw(11) <<  orig;
            t.Start();
            for (int i = 0; i < Ntest; ++i)
            {
                Exp->IProductWRTBase(input, output2);
            }
            t.Stop();
            NekDouble col = t.TimePerTest(Ntest);
            cout << "   Collection: " << std::setw(11) << t.TimePerTest(Ntest);
            cout << "   Ratio: " << (orig/col) << endl;


        }

        // PhysDeriv Comparison
        cout << "PhysDeriv Op: Ntest = " << Ntest << endl;

        for(int N = 2; N < 11; ++N)
        {
            graph3D->SetExpansionsToPolyOrder(N);

            Exp = MemoryManager<MultiRegions::ExpList3D>::
                AllocateSharedPtr(vSession,graph3D);

            Exp->CreateCollections(ImpType);

            int nelmt = Exp->GetNumElmts();

            const int nq = Exp->GetNpoints();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq), tmp,tmp1,tmp2;
            Array<OneD, NekDouble> input(nq), diff1_0(nq), diff1_1(nq), diff1_2(nq);
            Array<OneD, NekDouble> diff2_0(nq), diff2_1(nq),  diff2_2(nq);

            Exp->GetCoords(xc, yc, zc);

            for (int i = 0; i < nq; ++i)
            {
                input[i] = sin(xc[i])*cos(yc[i])*cos(zc[i]);
            }

            Timer t;
            t.Start();
            for (int i = 0; i < Ntest; ++i)
            {
                // Do test by calling every element in loop
                for(int j = 0; j < nelmt; ++j)
                {
                    Exp->GetExp(j)->PhysDeriv(input +Exp->GetPhys_Offset(j),
                                              tmp  = diff1_0 +Exp->GetPhys_Offset(j),
                                              tmp1 = diff1_1 +Exp->GetPhys_Offset(j),
                                              tmp2 = diff1_2 +Exp->GetPhys_Offset(j));
                }
            }
            t.Stop();
            NekDouble orig = t.TimePerTest(Ntest);
            cout << "N: "<< std::setw(2) << N << " ExpList: " <<std::setw(11) <<  orig;            t.Start();
            for (int i = 0; i < Ntest; ++i)
            {
                Exp->PhysDeriv(input, diff2_0, diff2_1,diff2_2);
            }
            t.Stop();
            NekDouble col = t.TimePerTest(Ntest);
            cout << "   Collection: " << std::setw(11) << t.TimePerTest(Ntest);
            cout << "   Ratio: " << (orig/col) << endl;
        }

        // IProductWRTDerivBase Comparison
        cout << "IProductWRTDerivBase Op: Ntest = " << Ntest << endl;
        for(int N = 2; N < 11; ++N)
        {
            graph3D->SetExpansionsToPolyOrder(N);

            Exp = MemoryManager<MultiRegions::ExpList3D>::
                AllocateSharedPtr(vSession,graph3D);

            Exp->CreateCollections(ImpType);

            int nelmt = Exp->GetNumElmts();

            const int nq = Exp->GetNpoints();
            const int nc = Exp->GetNcoeffs();
            Array<OneD, NekDouble> xc(nq), yc(nq), zc(nq),tmp,tmp1;
            Array<OneD, NekDouble> input1(nq), input2(nq), input3(nq), output1(nc), output2(nc);
            Array<OneD, Array<OneD, NekDouble> > input(3);

            input[0] = input1; input[1] = input2; input[2] = input3;

            Exp->GetCoords(xc, yc, zc);
            for (int i = 0; i < nq; ++i)
            {
                input1[i] = sin(xc[i])*cos(yc[i])*sin(zc[i]);
                input2[i] = cos(xc[i])*sin(yc[i])*cos(zc[i]);
                input3[i] = cos(xc[i])*sin(yc[i])*sin(zc[i]);
            }

            Timer t;
            t.Start();
            for (int i = 0; i < Ntest; ++i)
            {
                // Do test by calling every element in loop
                for(int j = 0; j < nelmt; ++j)
                {
                    Exp->GetExp(j)->IProductWRTDerivBase(0,input1 + Exp->GetPhys_Offset(j),
                                                    tmp  = output1 + Exp->GetCoeff_Offset(j));
                    Exp->GetExp(j)->IProductWRTDerivBase(1,input2 + Exp->GetPhys_Offset(j),
                                                         tmp  = output2 + Exp->GetCoeff_Offset(j));
                }
                Vmath::Vadd(nc,output1,1,output2,1,output1,1);

                for(int j = 0; j < nelmt; ++j)
                {
                    Exp->GetExp(j)->IProductWRTDerivBase(2,input3 + Exp->GetPhys_Offset(j),
                                                         tmp  = output2 + Exp->GetCoeff_Offset(j));
                }
                Vmath::Vadd(nc,output1,1,output2,1,output1,1);
            }
            t.Stop();
            NekDouble orig = t.TimePerTest(Ntest);
            cout << "N: "<< std::setw(2) << N << " ExpList: " <<std::setw(11) <<  orig;
            t.Start();
            for (int i = 0; i < Ntest; ++i)
            {
                Exp->IProductWRTDerivBase(input,output2);
            }
            t.Stop();
            NekDouble col = t.TimePerTest(Ntest);
            cout << "   Collection: " << std::setw(11) << t.TimePerTest(Ntest);
            cout << "   Ratio: " << (orig/col) << endl;

        }
    }
    vSession->Finalise();

    return 0;
}

