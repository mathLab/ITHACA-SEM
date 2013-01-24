#include <QtGui/QtGui>

#include <vtkCamera.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkPointLocator.h>
#include <vtkMath.h>
#include "vtkEventQtSlotConnect.h"
#include "vtkPolyDataWriter.h"
#include <vtkLookupTable.h>

#include "MainWindow.h"

using namespace std;

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
    
    mSourceFilterDepthSort = vtkDepthSortPolyData::New();
    mSourceFilterDepthSort->SetInputConnection(mSourceFilterSmooth->GetOutputPort());
    mSourceFilterDepthSort->SetDirectionToBackToFront();
    mSourceFilterDepthSort->SetVector(1, 1, 1);
    mSourceFilterDepthSort->SetCamera(mSourceRenderer->GetActiveCamera());
    mSourceFilterDepthSort->SortScalarsOn();
    
    mSourceLookupTable = vtkLookupTable::New();
    mSourceLookupTable->SetTableRange(0,10);
    mSourceLookupTable->SetNumberOfColors(16);
    mSourceLookupTable->SetHueRange(0.0,0.6667);
    mSourceLookupTable->Build();
    
    mSourceMapper = vtkPolyDataMapper::New();
    mSourceMapper->SetInputConnection(mSourceFilterDepthSort->GetOutputPort());
    mSourceMapper->ImmediateModeRenderingOn();
    mSourceMapper->ScalarVisibilityOn();
    mSourceMapper->SetScalarModeToUsePointData();
    //mSourceMapper->SetColorModeToMapScalars();
    mSourceMapper->SetScalarRange(mSourceData->GetScalarRange());
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
    mSourceHeightPointsLabelActor = vtkActor2D::New();
    mSourceHeightPointsLabelActor->SetMapper(mSourceHeightPointsLabelMapper);

    mSourceRenderer->AddActor(mSourceActor);
    mSourceRenderer->AddActor(mSourcePointsActor);
    mSourceRenderer->AddActor(mSourcePointsLabelActor);
    mSourceRenderer->AddActor(mSourceHeightPointsLabelActor);
    mSourceRenderer->ResetCamera();
    mSourceVtk->GetRenderWindow()->AddRenderer(mSourceRenderer);

    Connections_s = vtkEventQtSlotConnect::New();
    Connections_s->Connect(mSourceVtk->GetRenderWindow()->GetInteractor(),
                           vtkCommand::RightButtonPressEvent,
                           this,
                           SLOT(CreateSourcePoint( vtkObject*, unsigned long, void*, void*, vtkCommand*)));

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
    QLabel* vFileHeightLabel = new QLabel(tr("Height:")); //Height label
    
    mFileSourceEditBox = new QLineEdit(tr(""));
    mFileSourceBrowse = new QPushButton(tr("Browse..."));
    mFileTargetEditBox = new QLineEdit(tr(""));
    
    mFileHeightEditBox = new QLineEdit(tr(""));//Edit box
    
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
    mFileGrid->addWidget(vFileHeightLabel, 3, 0);//Add widget of height label
    mFileGrid->addWidget(mFileSourceEditBox, 0, 1);
    mFileGrid->addWidget(mFileLandmarksEditBox, 1, 1);
    mFileGrid->addWidget(mFileTargetEditBox, 2, 1);
    mFileGrid->addWidget(mFileHeightEditBox, 3, 1);//Add widget of height edit box
    mFileGrid->addWidget(mFileSourceBrowse, 0, 2);
    mFileGrid->addWidget(mFileLandmarksBrowse, 1, 2);
    mFileGrid->addWidget(mFileTargetBrowse, 2, 2);
    mFileGrid->addWidget(mFileLoadButton, 4, 1); //Shifted downwards
    mFileGrid->addWidget(mFileExportLandmarksButton, 5, 1);//Shifted downward
    
    mFileBox = new QGroupBox(tr("Files"));
    mFileBox->setLayout(mFileGrid);
    mFileBox->setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    
    mUndoButton = new QPushButton(tr("Undo last point"));
    connect(mUndoButton, SIGNAL(clicked()), this, SLOT(UndoLastPoint()));
    
    mSettingsGrid = new QVBoxLayout;
    mSettingsGrid->addWidget(mFileBox);
    mSettingsGrid->addWidget(mUndoButton);
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
    ResetData();

    LoadSource();
    LoadSourceLandmarks();
    LoadTarget();

    Update();
}

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
    mSourceHeightPointData->GetPointData()->AddArray(Heights); //get the points dataa of the mSourcePoints Data, and we are then adding a new scalar array to it.
    mSourceHeightPointsVisible->SetInput(mSourceHeightPointData);
}


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

    mSourceRenderer->ResetCamera(mSourceActor->GetBounds());
}

void MainWindow::LoadTarget() {
    vtkPolyDataReader* vTargetReader = vtkPolyDataReader::New();
    vTargetReader->SetFileName(mFileTargetEditBox->text().toStdString().c_str());
    vTargetReader->Update();
    mTargetData = vTargetReader->GetOutput();
    mTargetFilterSmooth->SetInput(mTargetData);

    mTargetRenderer->ResetCamera(mTargetActor->GetBounds());
}

void MainWindow::LoadSourceLandmarks() {
    vtkPolyDataReader* vSourceLandmarks = vtkPolyDataReader::New();
    vSourceLandmarks->SetFileName(mFileLandmarksEditBox->text().toStdString().c_str());
    vSourceLandmarks->Update();
    mSourcePointsData = vSourceLandmarks->GetOutput();
    mSourceFilterGlyph->SetInput(mSourcePointsData);
    mSourcePointsIds->SetInput(mSourcePointsData);
}


void MainWindow::Update() {
    // Height Interpolation
    int vInterpDistance=10;
    vtkDoubleArray* vSurfaceData=InterpSurface(mSourceHeightPointData, mSourceData, vInterpDistance);
    mSourceData->GetPointData()->RemoveArray("HeightSurface");
    mSourceData->GetPointData()->SetScalars(vSurfaceData);
    mSourceData->Modified();
    mSourceLookupTable->SetTableRange(mSourceData->GetScalarRange());
    mSourceLookupTable->Build();
    mSourceMapper->SetLookupTable(mSourceLookupTable);

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

//Interp Surface

vtkDoubleArray* MainWindow::InterpSurface(vtkPolyData* pPointData, vtkPolyData* pSurface, int pInterpDistance) {
    //p point data:points which have been selected.
    //psurface: the points on the mesh
    //Below we set the interp distance as 
    int nSurfacePoints = pSurface->GetNumberOfPoints();//Return number of points in array
    vtkDoubleArray* vPointDataValues = dynamic_cast<vtkDoubleArray*>(pPointData->GetPointData()->GetArray("height"));
    
    
    //for vtkDoubleArray: you are able to set  number of components, values, name, insert next value,set value etc
    
    vtkDoubleArray* vOutput = vtkDoubleArray::New();
    vOutput->SetNumberOfComponents(1);
    vOutput->SetName("HeightSurface");
    vOutput->SetNumberOfValues(nSurfacePoints);
    vOutput->FillComponent(0, 0.0);
    
    if (pPointData->GetNumberOfPoints() < 4) return vOutput;//We check to see if the number of points which we have clicked on is more than 4- since we need the four cloest neighbours.
    
    vtkPointLocator* vLocator = vtkPointLocator::New();
    vLocator->SetDataSet(pPointData);//We set the data set as the poitns which we have selected on the mesh
    
    
   
    for (int i = 0; i < nSurfacePoints; ++i) { //For all the points on the mesh
        double p[3]; //third array of the double p
        vtkIdList* list = vtkIdList::New(); //we create a list, lists are used to pass id's between objects
        
        pSurface->GetPoint(i, p);//Copy point components (eg co-ordinates) into the array p[3] for each point in the mesh
        //This basically gives you the co-ordinate of each point from the mesh
        
        vLocator->FindClosestNPoints(4, p, list); //Find 4 closest points to p and store in list,
        //points are sorted from closest to furthest
        
        //int numPoints= min(4, int(list->GetNumberOfIds()));
        
        
        //cout << list->GetNumberOfIds() << endl;
        
        //Weighted average of 4 closest points
        double val = 0;
        double w = 0;
        int c = 0;
        for (int j = 0; j < 4; ++j) { //We set j to have the values 0,1,2,3
            int id = list->GetId(j);  //Re: the 4 closest points are stored in the array list- so we get the id of these 4 positions
            
            //This follows the procedure which we did above, here we get the co-ordinate of the clicked points
            double q[3];
            pPointData->GetPoint(id, q); //co-ordinate of the height points
            
            double d_pq = sqrt(vtkMath::Distance2BetweenPoints(p, q));
            //Gets the squared distance between the two points, i.e. (y2-y1)^2 + (x2-x1)^2: this result is then square-rooted. We basically find length of the line.
            
            if (d_pq <= pInterpDistance) {//If the distabce between p and q is less than the chosen interp distance
                val += vPointDataValues->GetValue(id) / (d_pq + 1E-10);//Getvalue gets the "data" at that particular index
                //1e10: So we dont divide by 10
                //vPointDataValues represent the height which you input on the screet
                w   += 1.0 / (d_pq + 1E-10);
                c   += 1;
            }
        }
        
        if (c >= 1) {
            val /= w;
        }
        else {
            val = 0.0;
        }

        vOutput->SetValue(i, val);//Sets value of point in mesh as weighted average
    }
    return vOutput;
}

//End Interp surface

//Create Source Point

void MainWindow::CreateSourcePoint(vtkObject* caller, unsigned long vtk_event, void* client_data, void* call_data, vtkCommand* command) {
    double p_s[3];
    double d_s;
    
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::SafeDownCast(caller);
    
    // consume event so the interactor style doesn't get it
    command->AbortFlagOn();
    int* position_s = iren->GetEventPosition();
    
    // Use World point picker to find a 3-space point on atrial surface
    // Point will be on camera bound if not on surface
    vtkWorldPointPicker *picker_source = vtkWorldPointPicker::New();
    picker_source->Pick(position_s[0], position_s[1], position_s[2], mSourceRenderer);
    picker_source->GetPickPosition(p_s);
    
    // Use Point locator to find nearest vertex id on atrial surface within a
    // given radius
    vtkPointLocator* vLocator_source = vtkPointLocator::New();
    vLocator_source->SetDataSet(mSourceData);
    vLocator_source->SetTolerance(0.1);
    int nearestPointId_source = vLocator_source->FindClosestPointWithinRadius(5.0, p_s, d_s);
    
    double Heights = mFileHeightEditBox->text().toDouble();
    
    // If we are on the surface (i.e. id != -1) then add this vertex to the list
    // of landmark points.
    if (nearestPointId_source >= 0) {
        mSourceData->GetPoints()->GetPoint(nearestPointId_source, p_s);//get all vertices in the mesh, get nearest point using the id., put in ps
        vtkPoints* pointList = AddPoint(mSourceHeightPointData->GetPoints(), p_s);
        vtkDataArray* pointScalars = AddScalar(mSourceHeightPointData->GetPointData()->GetScalars("height"), Heights);
        //mSourceHeightPointData->Reset();
        mSourceHeightPointData->GetPointData()->RemoveArray("height");
        mSourceHeightPointData->SetPoints(pointList);
        mSourceHeightPointData->GetPointData()->SetScalars(pointScalars);
//        mSourceHeightPointData->GetPoints()->InsertNextPoint(p_s);//list of point which were already clicked on, add point which you have just clicked on. We now want to add the scalar value of the hight to this list.
//        mSourceHeightPointData->GetPointData()->GetScalars("height")-> InsertNextTuple1(Heights);//Tuple is size 1- scalar.
        mSourceHeightPointData->Modified();
        cout << "Updated points" << endl;
    }
    
    // Update display
    Update();
}

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

vtkDoubleArray* MainWindow::AddScalar(vtkDataArray* array, double v)
{
    vtkDoubleArray* vVals = vtkDoubleArray::New();
    vVals->SetName("height");
    if (!array)
    {
        cout << "Null array." << endl;
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

void MainWindow::ExportTargetPoints() {
    QString pointsFile = QFileDialog::getSaveFileName(this,
                                                      tr("Export Landmark Points"), "", tr("Landmark Points (*.vtk)"));
    
    vtkPolyDataWriter *writer = vtkPolyDataWriter::New();
    writer->SetInput(mTargetPointsData);
    writer->SetFileName(pointsFile.toStdString().c_str());
    writer->Write();
}

void MainWindow::UndoLastPoint() {
    if (mTargetPointsData->GetPoints()->GetNumberOfPoints() > 0)
    {
        double * p = new double[3];
        vtkPoints* vTargetPoints = vtkPoints::New();
        for (unsigned int i = 0; i < mTargetPointsData->GetPoints()->GetNumberOfPoints() - 1; ++i)
        {
            p = mTargetPointsData->GetPoint(i);
            vTargetPoints->InsertNextPoint(p);
        }
        mTargetPointsData->SetPoints(vTargetPoints);
        delete[] p;
        Update();
    }
}
