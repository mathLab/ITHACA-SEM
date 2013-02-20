/*
 * VtkToPng.cpp
 *
 *  Created on: 20 Feb 2013
 *      Author: cc
 */

//#include <vtkVersion.h>
//#if VTK_MAJOR_VERSION == 6
//int main(int, char *argv[])
//{
//  std::cout << argv[0] << " requires VTK 5.10 or earlier. This VTK version is " << vtkVersion::GetVTKVersion() << std::endl;
//  return EXIT_SUCCESS;
//}
//#else
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#include <vtkGraphicsFactory.h>
#include <vtkImagingFactory.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkCamera.h>

#include <iostream>
using namespace std;

int main(int argc, char * argv[])
{
  // Setup offscreen rendering
  vtkSmartPointer<vtkGraphicsFactory> graphics_factory =
    vtkSmartPointer<vtkGraphicsFactory>::New();
  graphics_factory->SetOffScreenOnlyMode( 1);
  graphics_factory->SetUseMesaClasses( 1 );

  vtkSmartPointer<vtkImagingFactory> imaging_factory =
    vtkSmartPointer<vtkImagingFactory>::New();
  imaging_factory->SetUseMesaClasses( 1 );

  // Create a poly data reader
  cout << "Reading file: " << argv[1] << endl;
  vtkXMLUnstructuredGridReader* reader = vtkXMLUnstructuredGridReader::New();
  reader->SetFileName (argv[1]);
  reader->Update();
  vtkDataSet* data = reader->GetOutputAsDataSet();

  data->GetPointData()->SetActiveScalars("u");

  // Lookup table
  vtkSmartPointer<vtkLookupTable> lookup = vtkSmartPointer<vtkLookupTable>::New();
  lookup->SetHueRange(0.0,1.0);
  lookup->SetSaturationRange(1,1);
  lookup->SetTableRange(data->GetScalarRange());
  lookup->SetValueRange(1,1);
  lookup->Build();

  // Create a mapper and actor
  vtkSmartPointer<vtkDataSetMapper> mapper =
    vtkSmartPointer<vtkDataSetMapper>::New();
  mapper->SetInput(data);
  mapper->ImmediateModeRenderingOn();
  mapper->ScalarVisibilityOn();
  mapper->SetScalarModeToUsePointData();
  mapper->SetScalarRange(data->GetScalarRange());
  mapper->SetLookupTable(lookup);


  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  // The usual rendering stuff.
  vtkCamera *camera = vtkCamera::New();
      camera->SetPosition(0.0,-1.0,1.0);
      camera->SetFocalPoint(0,0,0);


  // A renderer and render window
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetOffScreenRendering( 1 );
  renderWindow->AddRenderer(renderer);

  // Add the actors to the scene
  renderer->AddActor(actor);
  renderer->SetBackground(0,0,0); // Background color white
  renderer->SetActiveCamera(camera);
  renderer->ResetCamera();

  renderWindow->Render();

  vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
    vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowToImageFilter->SetInput(renderWindow);
  windowToImageFilter->SetMagnification(4);
  windowToImageFilter->Update();

  cout << "Writing file: " << argv[2] << endl;
  vtkSmartPointer<vtkPNGWriter> writer =
    vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetFileName(argv[2]);
  writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  writer->Write();

  return EXIT_SUCCESS;
}
//#endif
