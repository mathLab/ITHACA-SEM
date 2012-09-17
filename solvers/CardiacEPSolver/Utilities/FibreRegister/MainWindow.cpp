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
    mSourceData = vtkPolyData::New();
    mTargetData = vtkPolyData::New();

    mSourceRenderer = vtkRenderer::New();
    mSourceRenderer->SetBackground(0.0,0.0,0.0);
    mTargetRenderer = vtkRenderer::New();
    mTargetRenderer->SetBackground(0.0,0.0,0.0);

    mSourceFilterSmooth = vtkSmoothPolyDataFilter::New();
    mSourceFilterSmooth->SetNumberOfIterations(200);
    mSourceFilterSmooth->SetInput(mSourceData);
    mTargetFilterSmooth = vtkSmoothPolyDataFilter::New();
    mTargetFilterSmooth->SetNumberOfIterations(200);
    mTargetFilterSmooth->SetInput(mTargetData);

    mSourceFilterDepthSort = vtkDepthSortPolyData::New();
    mSourceFilterDepthSort->SetInputConnection(mSourceFilterSmooth->GetOutputPort());
    mSourceFilterDepthSort->SetDirectionToBackToFront();
    mSourceFilterDepthSort->SetVector(1, 1, 1);
    mSourceFilterDepthSort->SetCamera(mSourceRenderer->GetActiveCamera());
    mSourceFilterDepthSort->SortScalarsOn();
    mTargetFilterDepthSort = vtkDepthSortPolyData::New();
    mTargetFilterDepthSort->SetInputConnection(mTargetFilterSmooth->GetOutputPort());
    mTargetFilterDepthSort->SetDirectionToBackToFront();
    mTargetFilterDepthSort->SetVector(1, 1, 1);
    mTargetFilterDepthSort->SetCamera(mTargetRenderer->GetActiveCamera());
    mTargetFilterDepthSort->SortScalarsOn();

    mSourceMapper = vtkPolyDataMapper::New();
    mSourceMapper->SetInputConnection(mSourceFilterDepthSort->GetOutputPort());
    mSourceMapper->ScalarVisibilityOff();
    mTargetMapper = vtkPolyDataMapper::New();
    mTargetMapper->SetInputConnection(mTargetFilterDepthSort->GetOutputPort());
    mTargetMapper->ScalarVisibilityOff();

    // VTK Actor
    mSourceActor = vtkActor::New();
    mSourceActor->SetMapper(mSourceMapper);
    mSourceActor->GetProperty()->SetColor(0.9,0.9,0.9);
    mTargetActor = vtkActor::New();
    mTargetActor->SetMapper(mTargetMapper);

    mSourcePointsData = vtkPolyData::New();
    mSourceSphere = vtkSphereSource::New();
    mSourceFilterGlyph = vtkGlyph3D::New();
    mSourceFilterGlyph->SetInput(mSourcePointsData);
    mSourceFilterGlyph->SetScaleModeToDataScalingOff();
    mSourceFilterGlyph->SetSourceConnection(mSourceSphere->GetOutputPort());
    mSourcePointsMapper = vtkPolyDataMapper::New();
    mSourcePointsMapper->SetInputConnection(mSourceFilterGlyph->GetOutputPort());
    mSourcePointsActor = vtkActor::New();
    mSourcePointsActor->SetMapper(mSourcePointsMapper);
    mSourcePointsActor->GetProperty()->SetColor(1.0,0.0,0.0);
    mTargetPointsData = vtkPolyData::New();
    mTargetSphere = vtkSphereSource::New();
    mTargetFilterGlyph = vtkGlyph3D::New();
    mTargetFilterGlyph->SetInput(mTargetPointsData);
    mTargetFilterGlyph->SetScaleModeToDataScalingOff();
    mTargetFilterGlyph->SetSourceConnection(mTargetSphere->GetOutputPort());
    mTargetPointsMapper = vtkPolyDataMapper::New();
    mTargetPointsMapper->SetInputConnection(mTargetFilterGlyph->GetOutputPort());
    mTargetPointsActor = vtkActor::New();
    mTargetPointsActor->SetMapper(mTargetPointsMapper);
    mTargetPointsActor->GetProperty()->SetColor(1.0,0.0,0.0);

    mSourceRenderer->AddActor(mSourceActor);
    mSourceRenderer->AddActor(mSourcePointsActor);
    mSourceRenderer->ResetCamera();
    mSourceVtk->GetRenderWindow()->AddRenderer(mSourceRenderer);

    mTargetRenderer->AddActor(mTargetActor);
    mTargetRenderer->AddActor(mTargetPointsActor);
    mTargetRenderer->ResetCamera();
    mTargetVtk->GetRenderWindow()->AddRenderer(mTargetRenderer);

    Connections = vtkEventQtSlotConnect::New();
    Connections->Connect(mTargetVtk->GetRenderWindow()->GetInteractor(),
                         vtkCommand::RightButtonPressEvent,
                         this,
                         SLOT(CreateTargetPoint( vtkObject*, unsigned long, void*, void*, vtkCommand*)));

}

MainWindow::~MainWindow() {

}


void MainWindow::Draw() {
    resize(800, 600);

    mSourceVtk = new QVTKWidget(this, QFlag(0) );
    mTargetVtk = new QVTKWidget(this, QFlag(0) );

    QLabel* vFileSourceLabel = new QLabel(tr("Source:"));
    QLabel* vFilePointsLabel = new QLabel(tr("Landmarks:"));
    QLabel* vFileTargetLabel = new QLabel(tr("Target:"));

    mFileSourceEditBox = new QLineEdit(tr(""));
    mFileSourceBrowse = new QPushButton(tr("Browse..."));
    mFileTargetEditBox = new QLineEdit(tr(""));
    mFileTargetBrowse = new QPushButton(tr("Browse..."));
    mFileLandmarksEditBox = new QLineEdit(tr(""));
    mFileLandmarksBrowse = new QPushButton(tr("Browse..."));
    connect(mFileSourceBrowse, SIGNAL(clicked()), this, SLOT(BrowseSource()));
    connect(mFileTargetBrowse, SIGNAL(clicked()), this, SLOT(BrowseTarget()));
    connect(mFileLandmarksBrowse, SIGNAL(clicked()), this, SLOT(BrowseLandmarks()));

    mFileLoadButton = new QPushButton(tr("Load"));
    connect(mFileLoadButton, SIGNAL(clicked()), this, SLOT(Load()));
    mFileExportLandmarksButton = new QPushButton(tr("Export landmarks..."));
    connect(mFileExportLandmarksButton, SIGNAL(clicked()), this, SLOT(ExportTargetPoints()));

    mFileGrid = new QGridLayout;
    mFileGrid->addWidget(vFileSourceLabel, 0, 0);
    mFileGrid->addWidget(vFilePointsLabel, 1, 0);
    mFileGrid->addWidget(vFileTargetLabel, 2, 0);
    mFileGrid->addWidget(mFileSourceEditBox, 0, 1);
    mFileGrid->addWidget(mFileLandmarksEditBox, 1, 1);
    mFileGrid->addWidget(mFileTargetEditBox, 2, 1);
    mFileGrid->addWidget(mFileSourceBrowse, 0, 2);
    mFileGrid->addWidget(mFileLandmarksBrowse, 1, 2);
    mFileGrid->addWidget(mFileTargetBrowse, 2, 2);
    mFileGrid->addWidget(mFileLoadButton, 3, 1);
    mFileGrid->addWidget(mFileExportLandmarksButton, 4, 1);

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
    mRootGrid->addWidget(mSourceVtk, 0, 1);
    mRootGrid->addWidget(mTargetVtk, 0, 2);
    mRootGrid->addWidget(mSettingsBox, 0, 0);

    mRootWidget = new QWidget;
    mRootWidget->setLayout(mRootGrid);

    setCentralWidget(mRootWidget);

}


void MainWindow::BrowseSource() {
    QString geoFile = QFileDialog::getOpenFileName(this,
            tr("Load Source Geometry"), "", tr("Source Geometry (*.vtk)"));
    mFileSourceEditBox->setText(geoFile);
}

void MainWindow::BrowseTarget() {
    QString dataFile = QFileDialog::getOpenFileName(this,
            tr("Load Target Geometry"), "", tr("Target Geometry (*.vtk)"));
    mFileTargetEditBox->setText(dataFile);
}

void MainWindow::BrowseLandmarks() {
    QString pointsFile = QFileDialog::getOpenFileName(this,
            tr("Load Source Landmark Points"), "", tr("Landmark Points (*.vtk)"));
    mFileLandmarksEditBox->setText(pointsFile);
}

void MainWindow::Load() {
    vtkPolyDataReader* vSourceReader = vtkPolyDataReader::New();
    vSourceReader->SetFileName(mFileSourceEditBox->text().toStdString().c_str());
    vSourceReader->Update();
    mSourceData = vSourceReader->GetOutput();
    mSourceFilterSmooth->SetInput(mSourceData);

    vtkPolyDataReader* vTargetReader = vtkPolyDataReader::New();
    vTargetReader->SetFileName(mFileTargetEditBox->text().toStdString().c_str());
    vTargetReader->Update();
    mTargetData = vTargetReader->GetOutput();
    mTargetFilterSmooth->SetInput(mTargetData);

    vtkPolyDataReader* vSourceLandmarks = vtkPolyDataReader::New();
    vSourceLandmarks->SetFileName(mFileLandmarksEditBox->text().toStdString().c_str());
    vSourceLandmarks->Update();
    mSourcePointsData = vSourceLandmarks->GetOutput();
    mSourceFilterGlyph->SetInput(mSourcePointsData);

    vtkPoints* vTargetPoints = vtkPoints::New();
    mTargetPointsData->SetPoints(vTargetPoints);

    mSourceRenderer->ResetCamera(mSourceActor->GetBounds());
    mTargetRenderer->ResetCamera(mTargetActor->GetBounds());

    Update();
}

void MainWindow::Update() {
    mSourceVtk->update();
    mTargetVtk->update();
}

void MainWindow::CreateTargetPoint(vtkObject* caller, unsigned long vtk_event, void* client_data, void* call_data, vtkCommand* command) {
    double p[3];
    double d;

    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(caller);

    // consume event so the interactor style doesn't get it
    command->AbortFlagOn();
    int* position = iren->GetEventPosition();

    // Use World point picker to find a 3-space point on atrial surface
    // Point will be on camera bound if not on surface
    vtkWorldPointPicker *picker = vtkWorldPointPicker::New();
    picker->Pick(position[0], position[1], position[2], mTargetRenderer);
    picker->GetPickPosition(p);

    // Use Point locator to find nearest vertex id on atrial surface within a
    // given radius
    vtkPointLocator* vLocator = vtkPointLocator::New();
    vLocator->SetDataSet(mTargetData);
    vLocator->SetTolerance(0.1);
    int nearestPointId = vLocator->FindClosestPointWithinRadius(5.0, p, d);

    // If we are on the surface (i.e. id != -1) then add this vertex to the list
    // of landmark points.
    if (nearestPointId >= 0) {
        mTargetData->GetPoints()->GetPoint(nearestPointId, p);
        mTargetPointsData->GetPoints()->InsertNextPoint(p);
        mTargetPointsData->Modified();
    }

    // Update display
    Update();
}

void MainWindow::ExportTargetPoints() {
    QString pointsFile = QFileDialog::getSaveFileName(this,
            tr("Export Landmark Points"), "", tr("Landmark Points (*.vtk)"));

    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(mTargetPointsData);
    writer->SetFileName(pointsFile.toStdString().c_str());
    writer->Write();
}
