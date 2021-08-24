#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/VtkUtil.hpp>

#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>

int main(int argc, char* argv[])
{
    vtkIdType npts;
    vtkIdType* pts;

    // Read mesh
    vtkPolyDataReader *vtkMeshReader = vtkPolyDataReader::New();
    vtkMeshReader->SetFileName(argv[1]);
    vtkMeshReader->Update();
    vtkPolyData *vtkMesh = vtkMeshReader->GetOutput();
    vtkPoints *vtkPoints = vtkMesh->GetPoints();
    vtkCellArray *vtkStrips = vtkMesh->GetStrips();

    // Check we found points and strips in the file.
    ASSERTL0(vtkPoints, "ERROR: cannot get points from mesh.");
    ASSERTL0(vtkStrips, "ERROR: cannot get triangle strips from mesh.");

    // Create new cell array for polygons
    vtkCellArray *vtkPolys = vtkCellArray::New();

    // Generate the polygons from the triangle strips
    vtkStrips->InitTraversal();
    for (int i = 0; vtkStrips->GetNextCell(npts, pts); ++i)
    {
        for (int j = 0; j < npts - 2; ++j)
        {
            vtkPolys->InsertNextCell(3, &pts[j]);
        }
    }

    // Create the new poly data
    vtkPolyData *vtkNewMesh = vtkPolyData::New();
    vtkNewMesh->SetPoints(vtkPoints);
    vtkNewMesh->SetPolys(vtkPolys);

    // Copy data across
    for (int i = 0; i < vtkMesh->GetPointData()->GetNumberOfArrays(); ++i)
    {
        vtkNewMesh->GetPointData()->SetScalars(
                                        vtkMesh->GetPointData()->GetArray(i));
    }

    // Write out the new mesh
    vtkPolyDataWriter *vtkMeshWriter = vtkPolyDataWriter::New();
    vtkMeshWriter->SetFileName(argv[2]);
#if VTK_MAJOR_VERSION <= 5
    vtkMeshWriter->SetInput(vtkNewMesh);
#else
    vtkMeshWriter->SetInputData(vtkNewMesh);
#endif
    vtkMeshWriter->Write();
}
