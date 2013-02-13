#include <QtGui/QtGui>

#include <vtkCamera.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkPointLocator.h>
#include <vtkMath.h>
#include "vtkEventQtSlotConnect.h"
#include "vtkPolyDataWriter.h"
#include <vtkLookupTable.h>
#include <vtkFloatArray.h>

#include "MainWindow.h"

using namespace std;

// Maximum value to use for potential function
#define HEIGHT_MAX 99


/**
 * Constructs the graphical interface and VTK pipeline
 */
MainWindow::MainWindow( QWidget* parent, Qt::WindowFlags fl )
: QMainWindow( parent, fl ) {

    // Draw UI components
    Draw();
    
    // -----------------------------------------------
    // SOURCE
    // -----------------------------------------------
    mSourceData = vtkPolyData::New();
    
    mSourceRenderer = vtkRenderer::New();
    mSourceRenderer->SetBackground(0.0,0.0,0.0);
    
    mSourceFilterSmooth = vtkSmoothPolyDataFilter::New();
    mSourceFilterSmooth->SetNumberOfIterations(200);
    mSourceFilterSmooth->SetInput(mSourceData);
    
    mSourceNormals = vtkPolyDataNormals::New();
    mSourceNormals->SetInputConnection(mSourceFilterSmooth->GetOutputPort());
    mSourceNormals->SetFeatureAngle(60.0);
    mSourceNormals->ComputeCellNormalsOff();
    mSourceNormals->ComputePointNormalsOn();

    mSourceFilterDepthSort = vtkDepthSortPolyData::New();
    mSourceFilterDepthSort->SetInputConnection(mSourceNormals->GetOutputPort());
    mSourceFilterDepthSort->SetDirectionToBackToFront();
    mSourceFilterDepthSort->SetVector(1, 1, 1);
    mSourceFilterDepthSort->SetCamera(mSourceRenderer->GetActiveCamera());
    mSourceFilterDepthSort->SortScalarsOn();
    
    mSourceLookupTable = vtkLookupTable::New();
    mSourceLookupTable->SetTableRange(0,HEIGHT_MAX);
    mSourceLookupTable->SetHueRange(0.6,0.6);
    mSourceLookupTable->SetSaturationRange(0,1);
    mSourceLookupTable->SetValueRange(1,1);
    mSourceLookupTable->Build();
    
    mSourceMapper = vtkPolyDataMapper::New();
    mSourceMapper->SetInputConnection(mSourceFilterDepthSort->GetOutputPort());
    mSourceMapper->ImmediateModeRenderingOn();
    mSourceMapper->ScalarVisibilityOn();
    mSourceMapper->SetScalarModeToUsePointData();
    //mSourceMapper->SetColorModeToMapScalars();
    mSourceMapper->SetScalarRange(0, HEIGHT_MAX);
    mSourceMapper->SetLookupTable(mSourceLookupTable);
    
    mSourceActor = vtkActor::New();
    mSourceActor->SetMapper(mSourceMapper);
    
    // --- LANDMARK POINTS -----
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
    
    // ---- LANDMARK LABELS ----
    mSourcePointsIds = vtkIdFilter::New();
    mSourcePointsIds->SetInput(mSourcePointsData);
    mSourcePointsIds->PointIdsOn();
    mSourcePointsVisible = vtkSelectVisiblePoints::New();
    mSourcePointsVisible->SetInputConnection(mSourcePointsIds->GetOutputPort());
    mSourcePointsVisible->SetRenderer(mSourceRenderer);
    mSourcePointsLabelMapper = vtkLabeledDataMapper::New();
    mSourcePointsLabelMapper->SetInputConnection(mSourcePointsVisible->GetOutputPort());
    mSourcePointsLabelMapper->SetLabelModeToLabelFieldData();
    mSourcePointsLabelActor = vtkActor2D::New();
    mSourcePointsLabelActor->SetMapper(mSourcePointsLabelMapper);
    
    // ---- HEIGHT POINTS -----
    mSourceHeightPointData = vtkPolyData::New();
    mSourceHeightPointsVisible = vtkSelectVisiblePoints::New();
    mSourceHeightPointsVisible->SetInput(mSourceHeightPointData);
    mSourceHeightPointsVisible->SetRenderer(mSourceRenderer);
    mSourceHeightPointsLabelMapper = vtkLabeledDataMapper::New();
    mSourceHeightPointsLabelMapper->SetInputConnection(mSourceHeightPointsVisible->GetOutputPort());
    mSourceHeightPointsLabelMapper->SetLabelModeToLabelFieldData();
    mSourceHeightPointsLabelMapper->SetLabelFormat("%2.0f");
    mSourceHeightPointsLabelActor = vtkActor2D::New();
    mSourceHeightPointsLabelActor->SetMapper(mSourceHeightPointsLabelMapper);

    mSourceHeightGlyph = vtkGlyph3D::New();
    mSourceHeightGlyph->SetInput(mSourceHeightPointData);
    mSourceHeightGlyph->SetScaleModeToDataScalingOff();
    mSourceHeightGlyph->SetSourceConnection(mSourceSphere->GetOutputPort());
    mSourceHeightPointsMapper = vtkPolyDataMapper::New();
    mSourceHeightPointsMapper->SetInputConnection(mSourceHeightGlyph->GetOutputPort());
    mSourceHeightPointsMapper->SetScalarRange(0.0, HEIGHT_MAX);
    mSourceHeightPointsMapper->SetLookupTable(mSourceLookupTable);
    mSourceHeightPointsActor = vtkActor::New();
    mSourceHeightPointsActor->SetMapper(mSourceHeightPointsMapper);

    mSourceHeightContours = vtkContourFilter::New();
    mSourceHeightContours->SetInput(mSourceData);
    mSourceHeightContours->GenerateValues(21, 0.0, HEIGHT_MAX);
    mSourceHeightContourMapper = vtkPolyDataMapper::New();
    mSourceHeightContourMapper->SetInput(mSourceHeightContours->GetOutput());
    mSourceHeightContourMapper->SetScalarRange(0.0, HEIGHT_MAX);
    mSourceHeightContourActor = vtkActor::New();
    mSourceHeightContourActor->SetMapper(mSourceHeightContourMapper);
    mSourceHeightContourActor->GetProperty()->SetLineWidth(2.0);

    mSourceRenderer->AddActor(mSourceActor);
    mSourceRenderer->AddActor(mSourcePointsActor);
    mSourceRenderer->AddActor(mSourcePointsLabelActor);
    mSourceRenderer->AddActor(mSourceHeightPointsLabelActor);
    mSourceRenderer->AddActor(mSourceHeightPointsActor);
    mSourceRenderer->AddActor(mSourceHeightContourActor);
    mSourceRenderer->ResetCamera();
    mSourceVtk->GetRenderWindow()->AddRenderer(mSourceRenderer);

    Connections_s = vtkEventQtSlotConnect::New();
    Connections_s->Connect(mSourceVtk->GetRenderWindow()->GetInteractor(),
                           vtkCommand::RightButtonPressEvent,
                           this,
                           SLOT(CreateSourcePoint( vtkObject*, unsigned long,
                                                   void*, void*, vtkCommand*)));

    // -------------------------------------------


    // -------------------------------------------
    // TARGET
    // -------------------------------------------
    mTargetData = vtkPolyData::New();

    mTargetRenderer = vtkRenderer::New();
    mTargetRenderer->SetBackground(0.0,0.0,0.0);

    mTargetFilterSmooth = vtkSmoothPolyDataFilter::New();
    mTargetFilterSmooth->SetNumberOfIterations(200);
    mTargetFilterSmooth->SetInput(mTargetData);

    mTargetFilterDepthSort = vtkDepthSortPolyData::New();
    mTargetFilterDepthSort->SetInputConnection(mTargetFilterSmooth->GetOutputPort());
    mTargetFilterDepthSort->SetDirectionToBackToFront();
    mTargetFilterDepthSort->SetVector(1, 1, 1);
    mTargetFilterDepthSort->SetCamera(mTargetRenderer->GetActiveCamera());
    mTargetFilterDepthSort->SortScalarsOn();

    mTargetMapper = vtkPolyDataMapper::New();
    mTargetMapper->SetInputConnection(mTargetFilterDepthSort->GetOutputPort());
    mTargetMapper->ScalarVisibilityOff();

    mTargetActor = vtkActor::New();
    mTargetActor->SetMapper(mTargetMapper);

    // --- LANDMARK POINTS -----
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

    // ---- LANDMARK LABELS ----
    mTargetPointsIds = vtkIdFilter::New();
    mTargetPointsIds->SetInput(mTargetPointsData);
    mTargetPointsIds->PointIdsOn();
    mTargetPointsVisible = vtkSelectVisiblePoints::New();
    mTargetPointsVisible->SetInputConnection(mTargetPointsIds->GetOutputPort());
    mTargetPointsVisible->SetRenderer(mTargetRenderer);
    mTargetPointsLabelMapper = vtkLabeledDataMapper::New();
    mTargetPointsLabelMapper->SetInputConnection(mTargetPointsVisible->GetOutputPort());
    mTargetPointsLabelMapper->SetLabelModeToLabelFieldData();
    mTargetPointsLabelActor = vtkActor2D::New();
    mTargetPointsLabelActor->SetMapper(mTargetPointsLabelMapper);
    
    mTargetRenderer->AddActor(mTargetActor);
    mTargetRenderer->AddActor(mTargetPointsActor);
    mTargetRenderer->AddActor(mTargetPointsLabelActor);
    mTargetRenderer->ResetCamera();
    mTargetVtk->GetRenderWindow()->AddRenderer(mTargetRenderer);
    
    Connections = vtkEventQtSlotConnect::New();
    Connections->Connect(mTargetVtk->GetRenderWindow()->GetInteractor(),
                         vtkCommand::RightButtonPressEvent,
                         this,
                         SLOT(CreateTargetPoint( vtkObject*, unsigned long, void*, void*, vtkCommand*)));

    newHeightValue = mHeightSlider->minimum();
}


/**
 * Main window destructor
 */
MainWindow::~MainWindow() {
    
}


/**
 * Draw the GUI components.
 */
void MainWindow::Draw() {
    resize(800, 600);
    
    mSourceVtk = new QVTKWidget(this, QFlag(0) );
    mTargetVtk = new QVTKWidget(this, QFlag(0) );
    
    // File box
    QLabel* vFileSourceLabel = new QLabel(tr("Source:"));
    QLabel* vFilePointsLabel = new QLabel(tr("Landmarks:"));
    QLabel* vFileTargetLabel = new QLabel(tr("Target:"));
    
    mFileSourceEditBox = new QLineEdit(tr(""));
    mFileSourceBrowse = new QPushButton(tr("Browse..."));
    mFileTargetEditBox = new QLineEdit(tr(""));
    
    mFileTargetBrowse = new QPushButton(tr("Browse..."));
    connect(mFileSourceBrowse, SIGNAL(clicked()), this, SLOT(BrowseSource()));
    connect(mFileTargetBrowse, SIGNAL(clicked()), this, SLOT(BrowseTarget()));
    
    mFileLoadButton = new QPushButton(tr("Load"));
    connect(mFileLoadButton, SIGNAL(clicked()), this, SLOT(Load()));
    
    mFileExportSourceButton = new QPushButton(tr("Export Source..."));
    connect(mFileExportSourceButton, SIGNAL(clicked()), this, SLOT(ExportSource()));

    mFileGrid = new QGridLayout;
    mFileGrid->addWidget(vFileSourceLabel, 0, 0);
    mFileGrid->addWidget(vFileTargetLabel, 1, 0);
    mFileGrid->addWidget(mFileSourceEditBox, 0, 1);
    mFileGrid->addWidget(mFileTargetEditBox, 1, 1);
    mFileGrid->addWidget(mFileSourceBrowse, 0, 2);
    mFileGrid->addWidget(mFileTargetBrowse, 1, 2);
    mFileGrid->addWidget(mFileLoadButton, 2, 1); //Shifted downwards
    mFileGrid->addWidget(mFileExportSourceButton, 3, 1);
    
    mFileBox = new QGroupBox(tr("Files"));
    mFileBox->setLayout(mFileGrid);
    mFileBox->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);

    // Landmark Box
    mLandmarksLoadButton = new QPushButton(tr("Load landmark points."));
    connect(mLandmarksLoadButton, SIGNAL(clicked()), this, SLOT(LoadSourceLandmarks()));
    mUndoLandmarkButton = new QPushButton(tr("Undo last landmark point"));
    connect(mUndoLandmarkButton, SIGNAL(clicked()), this, SLOT(UndoLastLandmarkPoint()));
    mLandmarksExportButton = new QPushButton(tr("Export landmarks..."));
    connect(mLandmarksExportButton, SIGNAL(clicked()), this, SLOT(ExportTargetPoints()));

    mLandmarkGrid = new QGridLayout;
    mLandmarkGrid->addWidget(mLandmarksLoadButton, 0, 0);
    mLandmarkGrid->addWidget(mUndoLandmarkButton, 1, 0);
    mLandmarkGrid->addWidget(mLandmarksExportButton, 2, 0);//Shifted downward

    mLandmarkBox = new QGroupBox(tr("Landmarks"));
    mLandmarkBox->setLayout(mLandmarkGrid);
    mLandmarkBox->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);

    // Height box
    mLoadHeightButton = new QPushButton(tr("Load height points."));
    connect(mLoadHeightButton, SIGNAL(clicked()), this, SLOT(LoadHeightPoints()));
    mHeightSliderLabel = new QLabel(tr("Height (0):")); //Height label
    mHeightSlider = new QSlider;
    mHeightSlider->setRange(0, HEIGHT_MAX);
    mHeightSlider->setOrientation(Qt::Horizontal);
    connect(mHeightSlider, SIGNAL(sliderMoved(int)), this, SLOT(HeightValueChanged(int)));
    QLabel* vHeightInterpLabel = new QLabel(tr("Interp Range"));
    mHeightInterpRange = new QSlider;
    mHeightInterpRange->setRange(0, 30);
    mHeightInterpRange->setOrientation(Qt::Horizontal);
    connect(mHeightInterpRange, SIGNAL(sliderMoved(int)), this, SLOT(HeightInterpChanged(int)));
    mUndoHeightButton = new QPushButton(tr("Undo last point"));
    connect(mUndoHeightButton, SIGNAL(clicked()), this, SLOT(UndoLastHeightPoint()));
    mHeightExport = new QPushButton(tr("Export height points."));
    connect(mHeightExport, SIGNAL(clicked()), this, SLOT(ExportHeightPoints()));
    mHeightShowLabels = new QCheckBox(tr("Show height points"));
    mHeightShowLabels->setChecked(true);
    connect(mHeightShowLabels, SIGNAL(clicked(bool)), this, SLOT(HeightLabelsSelect(bool)));

    mHeightGrid = new QGridLayout;
    mHeightGrid->addWidget(mLoadHeightButton, 0, 1);
    mHeightGrid->addWidget(mHeightSliderLabel, 1, 0);
    mHeightGrid->addWidget(mHeightSlider, 1, 1);//Add widget of height edit box
    mHeightGrid->addWidget(mUndoHeightButton, 2, 1);
    mHeightGrid->addWidget(vHeightInterpLabel, 3, 0);
    mHeightGrid->addWidget(mHeightInterpRange, 3, 1);
    mHeightGrid->addWidget(mHeightExport, 4, 1);
    mHeightGrid->addWidget(mHeightShowLabels, 5, 1);

    mHeightBox = new QGroupBox(tr("Height function"));
    mHeightBox->setLayout(mHeightGrid);
    mHeightBox->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);


    mSettingsGrid = new QVBoxLayout;
    mSettingsGrid->addWidget(mFileBox);
    mSettingsGrid->addWidget(mLandmarkBox);
    mSettingsGrid->addWidget(mHeightBox);
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


/**
 * Prompt for a VTK file to load as the source geometry and store the filename
 * in the source filename textbox.
 */
void MainWindow::BrowseSource() {
    QString geoFile = QFileDialog::getOpenFileName(this,
            tr("Load Source Geometry"), "", tr("Source Geometry (*.vtk)"));
    mFileSourceEditBox->setText(geoFile);
}


/**
 * Prompt for a VTK file to load as the target geometry and store the filename
 * in the target filename textbox.
 */
void MainWindow::BrowseTarget() {
    QString dataFile = QFileDialog::getOpenFileName(this,
            tr("Load Target Geometry"), "", tr("Target Geometry (*.vtk)"));
    mFileTargetEditBox->setText(dataFile);
}


/**
 * Load the source and target geometries from the VTK files specified in the
 * source and target filename textboxes. The source and target data is reset
 * before load and the displays updated after loading the files.
 */
void MainWindow::Load() {
    ResetData();

    LoadSource();
    LoadTarget();

    Update();
}


/**
 * Reset all VTK polydata objects to the original state (i.e. empty)
 */
void MainWindow::ResetData() {
    mSourceData->Reset();
    mTargetData->Reset();
    mTargetPointsData->Reset();
    mSourceHeightPointData->Reset();

    vtkPoints* vTargetPoints = vtkPoints::New();
    mTargetPointsData->SetPoints(vTargetPoints);
    mTargetPointsIds->SetInput(mTargetPointsData);

    vtkPoints* vSourcePoints = vtkPoints::New();
    mSourceHeightPointData->SetPoints(vSourcePoints);

    vtkDoubleArray* Heights= vtkDoubleArray::New();
    Heights->SetName("height"); //setting name of the array
    mSourceHeightPointData->GetPointData()->AddArray(Heights);
    mSourceHeightPointsVisible->SetInput(mSourceHeightPointData);
}


/**
 * Load a source geometry from the file specified by the source filename textbox
 * into the source data VTK polydata structure. The file is loaded into a
 * temporary dataset and the points, cells and data copied into the main dataset
 * to overcome the issue of being unable to update the scalar data after
 * loading.
 */
void MainWindow::LoadSource() {
    // Load file into temporary poly data
    vtkPolyData* vTmpData = vtkPolyData::New();
    vtkPolyDataReader* vSourceReader = vtkPolyDataReader::New();
    vSourceReader->SetFileName(mFileSourceEditBox->text().toStdString().c_str());
    vSourceReader->Update();
    vTmpData = vSourceReader->GetOutput();

    // Populate mSourceData
    vtkPoints* vPts = vtkPoints::New();
    vtkCellArray* vCells = vtkCellArray::New();
    vPts->DeepCopy(vTmpData->GetPoints());
    vCells->DeepCopy(vTmpData->GetPolys());
    mSourceData->Reset();
    mSourceData->SetPoints(vPts);
    mSourceData->SetPolys(vCells);

    if (vTmpData->GetPointData()->HasArray("height"))
    {
        vtkDoubleArray* vData = vtkDoubleArray::New();
        vData->DeepCopy(vTmpData->GetPointData()->GetScalars("height"));
        mSourceData->GetPointData()->SetScalars(vData);
    }

    mSourceRenderer->ResetCamera(mSourceActor->GetBounds());
}


/**
 * Load the target geometry from the file specified by the target filename
 * textbox into the target data VTK polydata structure.
 */
void MainWindow::LoadTarget() {
    vtkPolyDataReader* vTargetReader = vtkPolyDataReader::New();
    vTargetReader->SetFileName(mFileTargetEditBox->text().toStdString().c_str());
    vTargetReader->Update();
    mTargetData = vTargetReader->GetOutput();
    mTargetFilterSmooth->SetInput(mTargetData);

    mTargetRenderer->ResetCamera(mTargetActor->GetBounds());
}


/**
 * Load the reference landmark points used for surface registration.
 */
void MainWindow::LoadSourceLandmarks() {
    QString pointsFile = QFileDialog::getOpenFileName(this,
            tr("Load Source Landmark Points"), "",
            tr("Landmark Points (*.vtk)"));

    vtkPolyDataReader* vSourceLandmarks = vtkPolyDataReader::New();
    vSourceLandmarks->SetFileName(pointsFile.toStdString().c_str());
    vSourceLandmarks->Update();
    mSourcePointsData = vSourceLandmarks->GetOutput();
    mSourceFilterGlyph->SetInput(mSourcePointsData);
    mSourcePointsIds->SetInput(mSourcePointsData);

    Update();
}


/**
 * Load the height potential function points used for prescribing the
 * interpolated height potential function and thereafter the fibre vector field.
 */
void MainWindow::LoadHeightPoints() {
    QString pointsFile = QFileDialog::getOpenFileName(this,
            tr("Load Height Points"), "", tr("Landmark Points (*.vtk)"));

    vtkPolyDataReader* vPointReader = vtkPolyDataReader::New();
    vPointReader->SetFileName(pointsFile.toStdString().c_str());
    vPointReader->Update();
    vtkPolyData* vTmp = vPointReader->GetOutput();
    unsigned int nPts = vTmp->GetNumberOfPoints();

    vtkPoints* pointList = vtkPoints::New();
    pointList->DeepCopy(vTmp->GetPoints());
    vtkDoubleArray* pointScalars = vtkDoubleArray::New();
    pointScalars->SetNumberOfTuples(nPts);
    pointScalars->SetName("height");
    pointScalars->DeepCopy(vTmp->GetPointData()->GetScalars("height"));

    mSourceHeightPointData->Reset();
    mSourceHeightPointData->SetPoints(pointList);
    mSourceHeightPointData->GetPointData()->SetScalars(pointScalars);
    mSourceHeightPointData->Modified();

    Update();
}


/**
 * Update the value of the slider label when the value changes.
 */
void MainWindow::HeightValueChanged(int value)
{
    mHeightSliderLabel->setText("Height (" + QString::number(value) + ") :");
}


/**
 * Update the display when the interpolation distance changes.
 */
void MainWindow::HeightInterpChanged(int value)
{
    newHeightValue = value;
    Update();
}


/**
 * Toggle the display of height potential point labels. Hiding these speeds up
 * the VTK rendering and display.
 */
void MainWindow::HeightLabelsSelect(bool value)
{
    if (value) {
        mSourceHeightPointsLabelActor->VisibilityOn();
        mSourceHeightPointsActor->VisibilityOn();
    }
    else {
        mSourceHeightPointsLabelActor->VisibilityOff();
        mSourceHeightPointsActor->VisibilityOff();
    }
    Update();
}


/**
 * Update all VTK displays. This includes interpolating the height potential
 * function and computing fibre directions.
 */
void MainWindow::Update() {
    mSourceData->GetPointData()->RemoveArray("HeightSurface");
    mSourceData->GetPointData()->RemoveArray("Gradient");

    vtkDoubleArray* vSurfaceData=InterpSurface( mSourceHeightPointData,
                                                mSourceData,
                                                newHeightValue);
    mSourceData->GetPointData()->SetScalars(vSurfaceData);

    ComputeFibreDirection();

    mSourceData->Modified();

    mSourceLookupTable->SetTableRange(0, HEIGHT_MAX);
    mSourceLookupTable->Build();
    mSourceMapper->SetScalarRange(0, HEIGHT_MAX);
    mSourceMapper->SetLookupTable(mSourceLookupTable);


    mSourceVtk->update();
    mTargetVtk->update();
}


/**
 * Adds a new landmark point to the target geometry.
 */
void MainWindow::CreateTargetPoint(vtkObject* caller, unsigned long vtk_event,
        void* client_data, void* call_data, vtkCommand* command) {

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


/**
 * Given the prescription of a function at a set of points, interpolates the
 * value of the function at all points in the given mesh.
 * @param   pPointData      Set of points with prescribed scalar values.
 * @param   pSurface        Surface points on which to set interpolated scalars.
 * @param   pInterpDistance Maximum distance around points to interpolate.
 * @returns                 Interpolated scalar values at mesh points.
 */
vtkDoubleArray* MainWindow::InterpSurface(vtkPolyData* pPointData,
        vtkPolyData* pSurface, int pInterpDistance) {

    const unsigned int nPoints = 20;
    double p[3];
    double q[3];
    vtkIdList* list = vtkIdList::New();

    // Number of points in surface
    int nSurfacePoints = pSurface->GetNumberOfPoints();
    
    vtkDoubleArray* vPointDataValues = vtkDoubleArray::SafeDownCast(
            pPointData->GetPointData()->GetScalars("height")
        );

    // Create output data array
    vtkDoubleArray* vOutput = vtkDoubleArray::New();
    vOutput->SetNumberOfComponents(1);
    vOutput->SetName("HeightSurface");
    vOutput->SetNumberOfValues(nSurfacePoints);
    
    // Return empty array if less than 4 points selected
    if (pPointData->GetNumberOfPoints() < nPoints)
    {
        return vOutput;
    }
    
    // Create a point locator to find the nearest points in the point list
    // to a given point in the surface.
    vtkPointLocator* vLocator = vtkPointLocator::New();
    vLocator->SetDataSet(pPointData);
    
    // Loop over every point in the mesh and
    for (int i = 0; i < nSurfacePoints; ++i) {

        // Retrieve the x,y,z coordinates of this surface point
        pSurface->GetPoint(i, p);

        // Seek the IDs of the nearest 4 points in the point list to this
        // coordinate.
        vLocator->FindClosestNPoints(nPoints, p, list);
        
        // Perform a weighted average of values at the 4 closest points
        double val = 0;
        double w = 0;
        int c = 0;

        for (int j = 0; j < nPoints; ++j) {
            // Get the j-th ID from the list.
            int id = list->GetId(j);
            
            // Get the coordiate of this point in the point list.
            pPointData->GetPoint(id, q);
            
            // Find the distance between the points.
            double d_pq = sqrt(vtkMath::Distance2BetweenPoints(p, q));
            
            // Use the point, only if it is within the interpolation distance.
            if (d_pq <= pInterpDistance) {
                // Add 1e10: So we dont divide by 0.
                val += vPointDataValues->GetValue(id) * (1.0 - (d_pq / pInterpDistance));

                w   += 1.0 - (d_pq / pInterpDistance);
                c   += 1;
            }
        }

        // If there was at least one point within the interpolation distance
        // we finish by dividing by the sum of the weights. If no points within
        // distance, we set the value to zero.
        if (c >= 1) {
            val /= w;
        }
        else {
            val = 0.0;
        }

        // Set the value in the output array.
        vOutput->SetValue(i, val);
    }
    return vOutput;
}


/**
 * Creates a new height potential function point with a value as chosen using
 * the height potential value slider.
 */
void MainWindow::CreateSourcePoint(vtkObject* caller, unsigned long vtk_event,
        void* client_data, void* call_data, vtkCommand* command) {

    double p_s[3];
    double d_s;
    
    vtkRenderWindowInteractor* iren
            = vtkRenderWindowInteractor::SafeDownCast(caller);
    
    // consume event so the interactor style doesn't get it
    command->AbortFlagOn();
    int* pos = iren->GetEventPosition();
    
    // Use World point picker to find a 3-space point on atrial surface
    // Point will be on camera bound if not on surface
    vtkWorldPointPicker *picker_source = vtkWorldPointPicker::New();
    picker_source->Pick(pos[0], pos[1], pos[2], mSourceRenderer);
    picker_source->GetPickPosition(p_s);
    
    // Use Point locator to find nearest vertex id on atrial surface within a
    // given radius
    vtkPointLocator* vLocator_source = vtkPointLocator::New();
    vLocator_source->SetDataSet(mSourceData);
    vLocator_source->SetTolerance(0.1);
    int nearestPointId_source
            = vLocator_source->FindClosestPointWithinRadius(5.0, p_s, d_s);
    
    double h = double(mHeightSlider->value());
    
    // If we are on the surface (i.e. id != -1) then add this vertex to the list
    // of height potential function points.
    if (nearestPointId_source >= 0) {
        mSourceData->GetPoints()->GetPoint(nearestPointId_source, p_s);

        vtkPoints* oldPointList = mSourceHeightPointData->GetPoints();
        vtkDataArray* oldPointScalars
                = mSourceHeightPointData->GetPointData()->GetScalars("height");

        vtkPoints* pointList = AddPoint(oldPointList, p_s);
        vtkDataArray* pointScalars = AddScalar(oldPointScalars, h);

        mSourceHeightPointData->GetPointData()->RemoveArray("height");
        mSourceHeightPointData->SetPoints(pointList);
        mSourceHeightPointData->GetPointData()->SetScalars(pointScalars);
        mSourceHeightPointData->Modified();
    }
    
    // Update display
    Update();
}


/**
 * Given a scalar function \phi, computes the fibre direction as the
 * cross-product of the \nabla\phi and the surface normal vector. Therefore, the
 * fibre direction lies along the contours of constant \phi.
 */
void MainWindow::ComputeFibreDirection()
{
    int nPts = mSourceData->GetNumberOfPoints();
    double b[3];
    double a[3], p[3];

    vtkDoubleArray* vGradientData = vtkDoubleArray::New();
    vGradientData->SetName("Gradient");
    vGradientData->SetNumberOfComponents(3);
    vGradientData->SetNumberOfTuples(nPts);

    vtkDoubleArray* vFibreData = vtkDoubleArray::New();
    vFibreData->SetName("FibreDirection");
    vFibreData->SetNumberOfComponents(3);
    vFibreData->SetNumberOfTuples(nPts);

    if (mSourceHeightPointData->GetNumberOfPoints() < 4)
    {
        cout << "Not enough points to compute fibre direction." << endl;
    }
    else
    {
        vtkGradientFilter* vGradient = vtkGradientFilter::New();
        vGradient->SetInput(mSourceData);
        vGradient->SetInputScalars(vtkDataObject::FIELD_ASSOCIATION_POINTS,
                                   "HeightSurface");
        vGradient->SetResultArrayName("Gradient");
        vGradient->Update();

        vtkDoubleArray* vGradVector = vtkDoubleArray::SafeDownCast(
                vGradient->GetOutput()->GetPointData()->GetVectors("Gradient")
            );
        vtkDataArray* vNormVector = vtkFloatArray::SafeDownCast(
                mSourceNormals->GetOutput()->GetPointData()->GetNormals()
            );

        vGradientData->DeepCopy(vGradVector);

        if (!vGradVector || !vNormVector)
        {
            cout << "One of the gradient or normals are not defined" << endl;
        }
        else
        {
            for (unsigned int i = 0; i < nPts; ++i)
            {
                vGradVector->GetTupleValue(i, a);
                vNormVector->GetTuple(i, b);
                vtkMath::Cross(a, b, p);
                vFibreData->SetTuple(i, p);
            }
        }
    }

    mSourceData->GetPointData()->AddArray(vGradientData);
    mSourceData->GetPointData()->AddArray(vFibreData);
}


/**
 * Generic routine to add a new point p to a dataset array.
 */
vtkPoints* MainWindow::AddPoint(vtkPoints* array, double* p)
{
    double q[3];
    vtkPoints* vPts = vtkPoints::New();
    int nPts = array->GetNumberOfPoints();
    vPts->SetNumberOfPoints(nPts+1);
    for (int i = 0; i < nPts; ++i)
    {
        array->GetPoint(i, q);
        vPts->SetPoint(i, q);
    }
    vPts->SetPoint(nPts, p);
    return vPts;
}


/**
 * Generic routine to add a new scalar value v to a double array.
 */
vtkDoubleArray* MainWindow::AddScalar(vtkDataArray* array, double v)
{
    vtkDoubleArray* vVals = vtkDoubleArray::New();
    vVals->SetName("height");
    if (!array)
    {
        cout << "WARNING: Null array given when adding scalar." << endl;
        return vVals;
    }

    int nPts = array->GetNumberOfTuples();
    vVals->SetNumberOfTuples(nPts+1);
    for (int i = 0; i < nPts; ++i)
    {
        vVals->SetTuple1(i, array->GetTuple1(i));
    }
    vVals->SetTuple1(nPts, v);
    return vVals;
}


/**
 * Prompt for a filename and write the source geometry and scalar/vector data
 * to the specified file.
 */
void MainWindow::ExportSource() {
    QString file = QFileDialog::getSaveFileName(this,
            tr("Export Source"), "", tr("Geometry (*.vtk)"));

    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(mSourceData);
    writer->SetFileName(file.toStdString().c_str());
    writer->Write();
}


/**
 * Prompt for a filename and write the target landmark points to the specified
 * file.
 */
void MainWindow::ExportTargetPoints() {
    QString pointsFile = QFileDialog::getSaveFileName(this,
            tr("Export Landmark Points"), "", tr("Landmark Points (*.vtk)"));
    
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(mTargetPointsData);
    writer->SetFileName(pointsFile.toStdString().c_str());
    writer->Write();
}


/**
 * Prompt for a filename and write the height potential function points
 * prescribed by the user to the specified file.
 */
void MainWindow::ExportHeightPoints() {
    QString pointsFile = QFileDialog::getSaveFileName(this,
            tr("Export Height Points"), "", tr("Height Points (*.vtk)"));

    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(mSourceHeightPointData);
    writer->SetFileName(pointsFile.toStdString().c_str());
    writer->Write();
}


/**
 * Removes the previously added landmark point from the target landmark point
 * dataset.
 */
void MainWindow::UndoLastLandmarkPoint() {
    if (mTargetPointsData->GetPoints()->GetNumberOfPoints() > 0)
    {
        double * p = new double[3];
        int i = 0;
        unsigned int n = 0;

        vtkPoints* vTargetPoints = vtkPoints::New();
        n = mTargetPointsData->GetPoints()->GetNumberOfPoints();
        for (i = 0; i < n - 1; ++i)
        {
            p = mTargetPointsData->GetPoint(i);
            vTargetPoints->InsertNextPoint(p);
        }
        mTargetPointsData->SetPoints(vTargetPoints);

        Update();

        delete[] p;
    }
}


/**
 * Removes the previously added height potential function point and its
 * corresponding scalar value.
 */
void MainWindow::UndoLastHeightPoint() {
    if (mSourceHeightPointData->GetPoints()->GetNumberOfPoints() > 0)
    {
        double * p = new double[3];
        int i = 0;
        unsigned int n = 0;

        vtkPoints* vPts = vtkPoints::New();
        vtkDataArray* vData = vtkDoubleArray::New();
        vData->SetName("height");
        vData->SetNumberOfTuples(n);
        n = mSourceHeightPointData->GetPoints()->GetNumberOfPoints();
        for (i = 0; i < n - 1; ++i)
        {
            p = mSourceHeightPointData->GetPoint(i);
            vPts->InsertNextPoint(p);
            vData->SetTuple1(i, mSourceHeightPointData->GetPointData()->
                                        GetScalars("height")->GetTuple1(i));
        }
        mSourceHeightPointData->GetPointData()->RemoveArray("height");
        mSourceHeightPointData->SetPoints(vPts);
        mSourceHeightPointData->GetPointData()->SetScalars(vData);

        Update();

        delete[] p;
    }
}
