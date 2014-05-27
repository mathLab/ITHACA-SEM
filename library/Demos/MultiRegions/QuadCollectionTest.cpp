#include <cstdio>
#include <cstdlib>

#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField2D.h>
#include <SpatialDomains/MeshGraph2D.h>

#include <Collections/Collection.h>

#define NBWD 1000000

using namespace Nektar;

int main(int argc, char *argv[])
{
    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv);

    MultiRegions::ContField2DSharedPtr Exp;

    if(argc < 2)
    {
        fprintf(stderr,"Usage: Helmholtz2D meshfile [SysSolnType]   or   \n");
        exit(1);
    }

    SpatialDomains::MeshGraphSharedPtr graph2D = 
        SpatialDomains::MeshGraph::Read(vSession);

    Exp = MemoryManager<MultiRegions::ContField2D>::
        AllocateSharedPtr(vSession,graph2D,vSession->GetVariable(0));

    Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0);
    Array<OneD, NekDouble> phys1(Exp->GetNpoints());
    Array<OneD, NekDouble> phys2(Exp->GetNpoints());

    Timer t;
    t.Start();
    for (int i = 0; i < NBWD; ++i)
    {
        Exp->BwdTrans(coeffs, phys1);
    }
    t.Stop();
    NekDouble orig = t.TimePerTest(NBWD);
    cout << "ExpList: " << orig << endl; 

    vector<SpatialDomains::GeometrySharedPtr> geom(Exp->GetExpSize());

    for (int i = 0; i < Exp->GetExpSize(); ++i)
    {
        geom[i] = Exp->GetExp(i)->GetGeom();
    }
    
    Collections::Collection c(Exp->GetExp(0), geom);

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

    vSession->Finalise();

    return 0;
}

