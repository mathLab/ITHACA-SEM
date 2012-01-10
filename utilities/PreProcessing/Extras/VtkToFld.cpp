#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ExpList2D.h>

#include <vtk-5.4/vtkPolyDataReader.h>
#include <vtk-5.4/vtkPolyData.h>
#include <vtk-5.4/vtkPointData.h>
#include <vtk-5.4/vtkPoints.h>
#include <vtk-5.4/vtkCellArray.h>
#include <vtk-5.4/vtkCellDataToPointData.h>

// Usage: VtkToFld session.xml input.vtk output.fld

int main(int argc, char* argv[])
{
    MultiRegions::ExpList2DSharedPtr Exp;

    std::vector<std::string> vFilenames;
    vFilenames.push_back(std::string(argv[1]));

    LibUtilities::SessionReaderSharedPtr vSession
            = LibUtilities::SessionReader::CreateInstance(argc, argv, vFilenames);
cout << "Loaded session" << endl;
    try
    {
        //----------------------------------------------
        // Read in mesh from input file
        SpatialDomains::MeshGraphSharedPtr graph2D = MemoryManager<SpatialDomains::MeshGraph2D>::AllocateSharedPtr(vSession);
        //----------------------------------------------
cout << "Created graph" << endl;
        //----------------------------------------------
        // Define Expansion
        Exp = MemoryManager<MultiRegions::ExpList2D>::
            AllocateSharedPtr(vSession,graph2D);
        //----------------------------------------------
cout << "Created ExpList" << endl;
        //----------------------------------------------
        // Set up coordinates of mesh for Forcing function evaluation
        int coordim = Exp->GetCoordim(0);
        int nq      = Exp->GetTotPoints();

        Array<OneD, NekDouble> xc0(nq,0.0);
        Array<OneD, NekDouble> xc1(nq,0.0);
        Array<OneD, NekDouble> xc2(nq,0.0);

        switch(coordim)
        {
        case 2:
            Exp->GetCoords(xc0,xc1);
            break;
        case 3:
            Exp->GetCoords(xc0,xc1,xc2);
            break;
        default:
            ASSERTL0(false,"Coordim not valid");
            break;
        }
        //----------------------------------------------

        vtkPolyDataReader *vtkMeshReader = vtkPolyDataReader::New();
        vtkMeshReader->SetFileName(argv[2]);
        vtkMeshReader->Update();
        vtkPolyData *vtkMesh = vtkMeshReader->GetOutput();

        vtkPoints *vtkPoints = vtkMesh->GetPoints();
        vtkCellArray *vtkPolys = vtkMesh->GetPolys();
        vtkCellDataToPointData* c2p = vtkCellDataToPointData::New();
        c2p->SetInput(vtkMesh);
        vtkPolyData *vtkDataAtPoints = c2p->GetPolyDataOutput();
        vtkPointData *vtkPointData = vtkDataAtPoints->GetPointData();

        for(int i = 0; i < nq; ++i)
        {
            cout << "Point: " << i << endl;
            for (int j = 0; j < vtkDataAtPoints->GetNumberOfPoints(); ++i)
            {
                double p[3];
                vtkDataAtPoints->GetPoint(j, p);
                if (fabs(p[0]-xc0[i]) + fabs(p[1]-xc1[i]) + fabs(p[2]-xc2[i]) < 1e-06)
                {
                    cout << "Found value" << endl;
                    Exp->UpdatePhys()[i] = vtkPointData->GetScalars("Image_Intensity")->GetTuple1(j);
                }
            }
        }

        //-----------------------------------------------
        // Write solution to file
        string   out(argv[3]);
        if (vSession->GetComm()->GetSize() > 1)
        {
            out += "." + boost::lexical_cast<string>(vSession->GetComm()->GetRank());
        }
        std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef
                                                    = Exp->GetFieldDefinitions();
        std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

        for(int i = 0; i < FieldDef.size(); ++i)
        {
            FieldDef[i]->m_fields.push_back("u");
            Exp->AppendFieldData(FieldDef[i], FieldData[i]);
        }
        graph2D->Write(out, FieldDef, FieldData);
        //-----------------------------------------------
    }
    catch (...) {
        cout << "ERROR" << endl;
    }
}
