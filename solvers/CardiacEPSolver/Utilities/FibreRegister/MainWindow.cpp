#include <QtGui/QtGui>

#include <vtkCamera.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkPointLocator.h>
#include <vtkMath.h>
//#include <vtkTextProperty.h>
//#include <vtkTable.h>
#include "vtkEventQtSlotConnect.h"
#include "vtkPolyDataWriter.h"

#include "MainWindow.h"

using namespace std;

MainWindow::MainWindow( QWidget* parent, Qt::WindowFlags fl )
    : QMainWindow( parent, fl ) {

    Draw();

    // Set up VTK stuff
    mLeftData = vtkPolyData::New();
    mRightData = vtkPolyData::New();

    mLeftRenderer = vtkRenderer::New();
    mLeftRenderer->SetBackground(0.0,0.0,0.0);
    mRightRenderer = vtkRenderer::New();
    mRightRenderer->SetBackground(0.0,0.0,0.0);

    mLeftFilterSmooth = vtkSmoothPolyDataFilter::New();
    mLeftFilterSmooth->SetNumberOfIterations(200);
    mLeftFilterSmooth->SetInput(mLeftData);
    mRightFilterSmooth = vtkSmoothPolyDataFilter::New();
    mRightFilterSmooth->SetNumberOfIterations(200);
    mRightFilterSmooth->SetInput(mRightData);

    mLeftFilterDepthSort = vtkDepthSortPolyData::New();
    mLeftFilterDepthSort->SetInputConnection(mLeftFilterSmooth->GetOutputPort());
    mLeftFilterDepthSort->SetDirectionToBackToFront();
    mLeftFilterDepthSort->SetVector(1, 1, 1);
    mLeftFilterDepthSort->SetCamera(mLeftRenderer->GetActiveCamera());
    mLeftFilterDepthSort->SortScalarsOn();
    mRightFilterDepthSort = vtkDepthSortPolyData::New();
    mRightFilterDepthSort->SetInputConnection(mRightFilterSmooth->GetOutputPort());
    mRightFilterDepthSort->SetDirectionToBackToFront();
    mRightFilterDepthSort->SetVector(1, 1, 1);
    mRightFilterDepthSort->SetCamera(mRightRenderer->GetActiveCamera());
    mRightFilterDepthSort->SortScalarsOn();

    mLeftMapper = vtkPolyDataMapper::New();
    mLeftMapper->SetInputConnection(mLeftFilterDepthSort->GetOutputPort());
    mLeftMapper->SetScalarModeToUsePointData();
    mLeftMapper->SetColorModeToMapScalars();
    //mLeftMapper->SetLookupTable(mColourTransfer);
    mLeftMapper->SetScalarRange(-10,500);
    mRightMapper = vtkPolyDataMapper::New();
    mRightMapper->SetInputConnection(mLeftFilterDepthSort->GetOutputPort());
    mRightMapper->SetScalarModeToUsePointData();
    mRightMapper->SetColorModeToMapScalars();
    //mRightMapper->SetLookupTable(mColourTransfer);
    mRightMapper->SetScalarRange(-10,500);

    // VTK Actor
    mLeftActor = vtkActor::New();
    mLeftActor->SetMapper(mLeftMapper);
    mRightActor = vtkActor::New();
    mRightActor->SetMapper(mLeftMapper);

    mLeftPointsData = vtkPolyData::New();
    mLeftSphere = vtkSphereSource::New();
    mLeftFilterGlyph = vtkGlyph3D::New();
    mLeftFilterGlyph->SetInput(mLeftPointsData);
    mLeftFilterGlyph->SetScaleModeToDataScalingOff();
    mLeftFilterGlyph->SetSourceConnection(mLeftSphere->GetOutputPort());
    mLeftPointsMapper = vtkPolyDataMapper::New();
    mLeftPointsMapper->SetInputConnection(mLeftFilterGlyph->GetOutputPort());
    mLeftPointsActor = vtkActor::New();
    mLeftPointsActor->SetMapper(mLeftPointsMapper);
    mLeftPointsActor->GetProperty()->SetColor(1.0,0.0,0.0);
    mRightPointsData = vtkPolyData::New();
    mRightSphere = vtkSphereSource::New();
    mRightFilterGlyph = vtkGlyph3D::New();
    mRightFilterGlyph->SetInput(mRightPointsData);
    mRightFilterGlyph->SetScaleModeToDataScalingOff();
    mRightFilterGlyph->SetSourceConnection(mRightSphere->GetOutputPort());
    mRightPointsMapper = vtkPolyDataMapper::New();
    mRightPointsMapper->SetInputConnection(mRightFilterGlyph->GetOutputPort());
    mRightPointsActor = vtkActor::New();
    mRightPointsActor->SetMapper(mRightPointsMapper);
    mRightPointsActor->GetProperty()->SetColor(1.0,0.0,0.0);

    mLeftRenderer->AddActor(mLeftActor);
    mLeftRenderer->AddActor(mLeftPointsActor);
    mLeftRenderer->ResetCamera();
    mLeftVtk->GetRenderWindow()->AddRenderer(mLeftRenderer);

    mRightRenderer->AddActor(mRightActor);
    mRightRenderer->AddActor(mRightPointsActor);
    mRightRenderer->ResetCamera();
    mRightVtk->GetRenderWindow()->AddRenderer(mRightRenderer);

    Connections = vtkEventQtSlotConnect::New();
    Connections->Connect(mLeftVtk->GetRenderWindow()->GetInteractor(),
                         vtkCommand::RightButtonPressEvent,
                         this,
                         SLOT(CreateLeftPoint( vtkObject*, unsigned long, void*, void*, vtkCommand*)));

}

MainWindow::~MainWindow() {

}


void MainWindow::Draw() {
    resize(800, 600);

    mLeftVtk = new QVTKWidget(this, QFlag(0) );
    mRightVtk = new QVTKWidget(this, QFlag(0) );

    QLabel* vFileLeftLabel = new QLabel(tr("Left:"));
    QLabel* vFilePointsLabel = new QLabel(tr("Landmarks:"));
    QLabel* vFileRightLabel = new QLabel(tr("Right:"));

    mFileLeftEditBox = new QLineEdit(tr(""));
    mFileLeftBrowse = new QPushButton(tr("Browse..."));
    mFileRightEditBox = new QLineEdit(tr(""));
    mFileRightBrowse = new QPushButton(tr("Browse..."));
    mFilePointsEditBox = new QLineEdit(tr(""));
    mFilePointsBrowse = new QPushButton(tr("Browse..."));
    connect(mFileLeftBrowse, SIGNAL(clicked()), this, SLOT(BrowseLeft()));
    connect(mFileRightBrowse, SIGNAL(clicked()), this, SLOT(BrowseRight()));
    connect(mFilePointsBrowse, SIGNAL(clicked()), this, SLOT(BrowsePoints()));

    mFileLoadButton = new QPushButton(tr("Load"));
    connect(mFileLoadButton, SIGNAL(clicked()), this, SLOT(Load()));
    mFileExportPointsButton = new QPushButton(tr("Export points..."));
    connect(mFileExportPointsButton, SIGNAL(clicked()), this, SLOT(ExportRightPoints()));

    mFileGrid = new QGridLayout;
    mFileGrid->addWidget(vFileLeftLabel, 0, 0);
    mFileGrid->addWidget(vFilePointsLabel, 1, 0);
    mFileGrid->addWidget(vFileRightLabel, 2, 0);
    mFileGrid->addWidget(mFileLeftEditBox, 0, 1);
    mFileGrid->addWidget(mFilePointsEditBox, 1, 1);
    mFileGrid->addWidget(mFileRightEditBox, 2, 1);
    mFileGrid->addWidget(mFileLeftBrowse, 0, 2);
    mFileGrid->addWidget(mFilePointsBrowse, 1, 2);
    mFileGrid->addWidget(mFileRightBrowse, 2, 2);
    mFileGrid->addWidget(mFileLoadButton, 3, 1);
    mFileGrid->addWidget(mFileExportPointsButton, 4, 1);

    mFileBox = new QGroupBox(tr("Files"));
    mFileBox->setLayout(mFileGrid);
    mFileBox->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);

    mSettingsGrid = new QVBoxLayout;
    mSettingsGrid->addWidget(mFileBox);
    mSettingsGrid->addStretch();

    mSettingsBox = new QGroupBox(tr("Settings"));
    mSettingsBox->setLayout(mSettingsGrid);

    mRootGrid = new QGridLayout;
    mRootGrid->setColumnMinimumWidth(0,300);
    mRootGrid->setColumnStretch(1,1);
    mRootGrid->setColumnStretch(2,1);
    mRootGrid->addWidget(mLeftVtk, 0, 1);
    mRootGrid->addWidget(mRightVtk, 0, 2);
    mRootGrid->addWidget(mSettingsBox, 0, 0);

    mRootWidget = new QWidget;
    mRootWidget->setLayout(mRootGrid);

    setCentralWidget(mRootWidget);

}


void MainWindow::BrowseLeft() {
    QString geoFile = QFileDialog::getOpenFileName(this,
            tr("Load Fibre Orientation"), "", tr("Source Fibre Orientation (*.vtk)"));
    mFileLeftEditBox->setText(geoFile);
}

void MainWindow::BrowseRight() {
    QString dataFile = QFileDialog::getOpenFileName(this,
            tr("Load Target Geometry"), "", tr("Target Geometry (*.vtk)"));
    mFileRightEditBox->setText(dataFile);
}

void MainWindow::BrowsePoints() {
    QString pointsFile = QFileDialog::getOpenFileName(this,
            tr("Load Source Landmark Points"), "", tr("Landmark Points (*.vtk)"));
    mFilePointsEditBox->setText(pointsFile);
}

void MainWindow::Load() {
    vtkPolyDataReader* vLeftReader = vtkPolyDataReader::New();
    vLeftReader->SetFileName(mFileLeftEditBox->text().toStdString().c_str());
    vLeftReader->Update();
    mLeftData = vLeftReader->GetOutput();
    mLeftFilterSmooth->SetInput(mLeftData);

    vtkPolyDataReader* vRightReader = vtkPolyDataReader::New();
    vRightReader->SetFileName(mFileRightEditBox->text().toStdString().c_str());
    vRightReader->Update();
    mRightData = vRightReader->GetOutput();
    mRightFilterSmooth->SetInput(mRightData);

    vtkPoints* vPointList = vtkPoints::New();
    mLeftPointsData->SetPoints(vPointList);

    mLeftRenderer->ResetCamera(mLeftActor->GetBounds());
    mRightRenderer->ResetCamera(mRightActor->GetBounds());

    Update();
}

void MainWindow::Update() {
    mLeftVtk->update();
    mRightVtk->update();
}

void MainWindow::CreateLeftPoint(vtkObject* caller, unsigned long vtk_event, void* client_data, void* call_data, vtkCommand* command) {
    double p[3];
    double d;

    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(caller);

    // consume event so the interactor style doesn't get it
    command->AbortFlagOn();
    int* position = iren->GetEventPosition();

    // Use World point picker to find a 3-space point on atrial surface
    // Point will be on camera bound if not on surface
    vtkWorldPointPicker *picker = vtkWorldPointPicker::New();
    picker->Pick(position[0], position[1], position[2], mLeftRenderer);
    picker->GetPickPosition(p);

    // Use Point locator to find nearest vertex id on atrial surface within a
    // given radius
    vtkPointLocator* vLocator = vtkPointLocator::New();
    vLocator->SetDataSet(mLeftData);
    vLocator->SetTolerance(0.1);
    int nearestPointId = vLocator->FindClosestPointWithinRadius(5.0, p, d);

    // If we are on the surface (i.e. id != -1) then add this vertex to the list
    // of landmark points.
    if (nearestPointId >= 0) {
        mLeftData->GetPoints()->GetPoint(nearestPointId, p);
        mLeftPointsData->GetPoints()->InsertNextPoint(p);
        mLeftPointsData->Modified();
    }

    // Update display
    Update();
}

void MainWindow::ExportRightPoints() {
    QString pointsFile = QFileDialog::getSaveFileName(this,
            tr("Export Landmark Points"), "", tr("Landmark Points (*.vtk)"));

    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(mLeftPointsData);
    writer->SetFileName(pointsFile.toStdString().c_str());
    writer->Write();
}
