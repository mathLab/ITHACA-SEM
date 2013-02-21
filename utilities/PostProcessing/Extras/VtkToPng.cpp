///////////////////////////////////////////////////////////////////////////////
//
// File VtkToPng.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Render a VTK Unstructured Grid file as a PNG
//
///////////////////////////////////////////////////////////////////////////////

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
    if (argc != 2 && argc != 4)
    {
        cout << "Usage: VtkToPng vtk-file [lower upper]" << endl;
        exit(-1);
    }

    string vInput = argv[1];
    string vOutput = vInput.substr(0, vInput.find_last_of('.')) + ".png";

    // Setup offscreen rendering
    vtkSmartPointer<vtkGraphicsFactory> graphics_factory =
            vtkSmartPointer<vtkGraphicsFactory>::New();
    graphics_factory->SetOffScreenOnlyMode( 1);
    graphics_factory->SetUseMesaClasses( 1 );

    vtkSmartPointer<vtkImagingFactory> imaging_factory =
            vtkSmartPointer<vtkImagingFactory>::New();
    imaging_factory->SetUseMesaClasses( 1 );

    // Create a poly data reader and retrieve dataset from file
    vtkXMLUnstructuredGridReader* reader = vtkXMLUnstructuredGridReader::New();
    reader->SetFileName(vInput.c_str());
    reader->Update();

    vtkDataSet* data = reader->GetOutputAsDataSet();
    data->GetPointData()->SetActiveScalars("u");

    double scalar_range[2];
    data->GetScalarRange(scalar_range);
    if (argc == 4)
    {
        scalar_range[0] = atof(argv[2]);
        scalar_range[1] = atof(argv[3]);
    }

    // Lookup table
    vtkSmartPointer<vtkLookupTable> lookup = vtkSmartPointer<vtkLookupTable>::New();
    lookup->SetHueRange(0.0,1.0);
    lookup->SetSaturationRange(1,1);
    lookup->SetTableRange(scalar_range);
    lookup->SetValueRange(1,1);
    lookup->Build();

    // Create a mapper and actor
    vtkSmartPointer<vtkDataSetMapper> mapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInput(data);
    mapper->ImmediateModeRenderingOn();
    mapper->ScalarVisibilityOn();
    mapper->SetScalarModeToUsePointData();
    mapper->UseLookupTableScalarRangeOn();
    //mapper->SetScalarRange(data->GetScalarRange());
    mapper->SetLookupTable(lookup);


    vtkSmartPointer<vtkActor> actor =
            vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);

    // Configure camera position and direction
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

    // Create an image of scene
    vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
            vtkSmartPointer<vtkWindowToImageFilter>::New();
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetMagnification(4);
    windowToImageFilter->Update();

    // Write image to PNG
    vtkSmartPointer<vtkPNGWriter> writer =
            vtkSmartPointer<vtkPNGWriter>::New();
    writer->SetFileName(vOutput.c_str());
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();

    return EXIT_SUCCESS;
}
