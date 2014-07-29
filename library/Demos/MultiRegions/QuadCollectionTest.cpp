#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField2D.h>
#include <SpatialDomains/MeshGraph2D.h>

#include <Collections/Collection.h>

#define NBWD 10000

using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ContField2DSharedPtr Exp;

    if(argc < 2)
    {
        fprintf(stderr,"Usage: QuadCollection meshfile [SysSolnType]   or   \n");
        exit(1);
    }

    SpatialDomains::MeshGraphSharedPtr graph2D = 
        SpatialDomains::MeshGraph::Read(vSession);

    Exp = MemoryManager<MultiRegions::ContField2D>::
        AllocateSharedPtr(vSession,graph2D,vSession->GetVariable(0));

    vector<SpatialDomains::GeometrySharedPtr> geom(Exp->GetExpSize());

    for (int i = 0; i < Exp->GetExpSize(); ++i)
    {
        geom[i] = Exp->GetExp(i)->GetGeom();
    }
    
    Collections::Collection c(Exp->GetExp(0), geom,Collections::eSumFac);

#if 1
    {
        Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
        Array<OneD, NekDouble> phys1(Exp->GetNpoints());
        Array<OneD, NekDouble> phys2(Exp->GetNpoints());

        Timer t;
        t.Start();
        for (int i = 0; i < NBWD; ++i)
        {
#if 0
            for(int j = 0; j < Exp->GetNumElmts(); ++j)
            {
                Exp->GetExp(j)->BwdTrans(coeffs+Exp->GetCoeff_Offset(j), tmp = phys1+Exp->GetPhys_Offset(j));
            }
#else
            Exp->BwdTrans(coeffs, phys1);
#endif
        }
        t.Stop();
        NekDouble orig = t.TimePerTest(NBWD);
        cout << "ExpList: " << orig << endl; 

        t.Start();
        for (int i = 0; i < NBWD; ++i)
        {
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);
        }
        t.Stop();
        NekDouble col = t.TimePerTest(NBWD);
        cout << "Collection: " << t.TimePerTest(NBWD) << endl;

        cout << "Ratio: " << (orig/col) << endl;

        Vmath::Vsub(phys1.num_elements(), phys1, 1, phys2, 1, phys1, 1);

        cout << Vmath::Vmax(phys1.num_elements(), phys1, 1) << endl;
    }
#endif

#if 0 
    {
        const int nq = Exp->GetNpoints();
        Array<OneD, NekDouble> xc(nq), yc(nq);
        Array<OneD, NekDouble> input(nq), outexp0(nq), outexp1(nq), outcol(2*nq);

        Exp->GetCoords(xc, yc);
        
        for (int i = 0; i < nq; ++i)
        {
            input[i] = sin(xc[i])*cos(yc[i]);
        }

        Timer t;
        t.Start();
        for (int i = 0; i < NBWD; ++i)
        {
            Exp->PhysDeriv(input, outexp0, NullNekDouble1DArray);
        }
        t.Stop();
        NekDouble orig = t.TimePerTest(NBWD);
        cout << "ExpList: " << orig << endl; 

        t.Start();
        for (int i = 0; i < NBWD; ++i)
        {
            c.ApplyOperator(Collections::ePhysDeriv, input, outcol);
        }
        t.Stop();
        NekDouble col = t.TimePerTest(NBWD);
        cout << "Collection: " << t.TimePerTest(NBWD) << endl;
        cout << "Ratio: " << (orig/col) << endl;

        Vmath::Vsub(nq, outexp0, 1, outcol, 1, outexp0, 1);
        cout << "Error: " << Vmath::Vmax(nq, outexp0, 1) << endl;
    }
#endif
    vSession->Finalise();

    return 0;
}

