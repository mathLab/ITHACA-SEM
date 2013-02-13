#ifndef CLASS_MAINWINDOW_H
#define CLASS_MAINWINDOW_H

#include <vector>
#include <string>
#include <sstream>

#include <QtGui>
#include <QVTKWidget.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkPointData.h>
#include <vtkActor.h>
#include <vtkActor2D.h>
#include <vtkProperty.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkDepthSortPolyData.h>
#include <vtkGlyph3D.h>
#include <vtkSphereSource.h>
#include <vtkWorldPointPicker.h>
#include <vtkIdFilter.h>
#include <vtkSelectVisiblePoints.h>
#include <vtkLabeledDataMapper.h>
#include <vtkContourFilter.h>
#include <vtkGradientFilter.h>
#include <vtkPolyDataNormals.h>

class vtkEventQtSlotConnect;


class MainWindow : public QMainWindow
{
    Q_OBJECT

    public:
        MainWindow( QWidget* parent = 0, Qt::WindowFlags fl = 0 );
        ~MainWindow();

    private slots:
        void BrowseSource();
        void BrowseTarget();

        void ResetData();
        void Load();
        void LoadSource();
        void LoadSourceLandmarks();
        void LoadHeightPoints();
        void LoadTarget();

        void HeightValueChanged(int value);
        void HeightInterpChanged(int value);
        void HeightLabelsSelect(bool value);

        void Update();
        
        void CreateTargetPoint(vtkObject*, unsigned long, void*, void*, vtkCommand*);
    
        void CreateSourcePoint(vtkObject*, unsigned long, void*, void*, vtkCommand*);

        void ExportSource();
        void ExportTargetPoints();
        void ExportHeightPoints();
        void UndoLastLandmarkPoint();
        void UndoLastHeightPoint();

    private:
        // GUI widgets
        QWidget* mRootWidget;
        QGridLayout* mRootGrid;
        QGroupBox* mSettingsBox;
        QVBoxLayout* mSettingsGrid;
        QGroupBox* mFileBox;
        QGridLayout* mFileGrid;
        QLineEdit* mFileSourceEditBox;
        QPushButton* mFileSourceBrowse;
        QLineEdit* mFileTargetEditBox;
        QPushButton* mFileTargetBrowse;
        QPushButton* mFileLoadButton;
        QPushButton* mFileExportSourceButton;

        QGroupBox* mLandmarkBox;
        QGridLayout* mLandmarkGrid;
        QPushButton* mLandmarksLoadButton;
        QPushButton* mUndoLandmarkButton;
        QPushButton* mLandmarksExportButton;

        QGroupBox* mHeightBox;
        QGridLayout* mHeightGrid;
        QPushButton* mLoadHeightButton;
        QPushButton* mUndoHeightButton;
        QLabel* mHeightSliderLabel;
        QSlider* mHeightSlider;
        QSlider* mHeightInterpRange;
        QCheckBox* mHeightShowLabels;
        QPushButton* mHeightExport;

        QVTKWidget* mSourceVtk;
        QVTKWidget* mTargetVtk;
    

        // VTK
        // Left Surface pipeline
        vtkPolyData* mSourceData;
        vtkLookupTable* mSourceLookupTable;
        vtkSmoothPolyDataFilter* mSourceFilterSmooth;
        vtkPolyDataNormals* mSourceNormals;
        vtkDepthSortPolyData* mSourceFilterDepthSort;
        vtkPolyDataMapper* mSourceMapper;
        vtkActor* mSourceActor;
        vtkRenderer* mSourceRenderer;

        // Right Surface pipeline
        vtkPolyData* mTargetData;
        vtkSmoothPolyDataFilter* mTargetFilterSmooth;
        vtkDepthSortPolyData* mTargetFilterDepthSort;
        vtkPolyDataMapper* mTargetMapper;
        vtkActor* mTargetActor;
        vtkRenderer* mTargetRenderer;

        // Source Points pipeline
        vtkPolyData* mSourcePointsData;
        vtkPolyData* mSourceHeightPointData;
        vtkSphereSource* mSourceSphere;
        vtkGlyph3D* mSourceFilterGlyph;
        vtkPolyDataMapper* mSourcePointsMapper;
        vtkActor* mSourcePointsActor;

        // Source Points Labels pipeline
        vtkIdFilter* mSourcePointsIds;
        vtkSelectVisiblePoints* mSourcePointsVisible;
        vtkSelectVisiblePoints* mSourceHeightPointsVisible;
        vtkLabeledDataMapper* mSourcePointsLabelMapper;
        vtkLabeledDataMapper* mSourceHeightPointsLabelMapper;
        vtkActor2D* mSourcePointsLabelActor;
        vtkActor2D* mSourceHeightPointsLabelActor;
    
        vtkGlyph3D* mSourceHeightGlyph;
        vtkPolyDataMapper* mSourceHeightPointsMapper;
        vtkActor* mSourceHeightPointsActor;

        // Source height contours
        vtkContourFilter* mSourceHeightContours;
        vtkPolyDataMapper* mSourceHeightContourMapper;
        vtkActor* mSourceHeightContourActor;

        // Source height gradient
        vtkPolyDataMapper* mSourceHeightVectorMapper;
        vtkActor* mSourceHeightVectorActor;

        // Target Points pipeline
        vtkPolyData* mTargetPointsData;
        vtkSphereSource* mTargetSphere;
        vtkGlyph3D* mTargetFilterGlyph;
        vtkPolyDataMapper* mTargetPointsMapper;
        vtkActor* mTargetPointsActor;

        // Target Points Labels pipeline
        vtkIdFilter* mTargetPointsIds;
        vtkSelectVisiblePoints* mTargetPointsVisible;
        vtkLabeledDataMapper* mTargetPointsLabelMapper;
        vtkActor2D* mTargetPointsLabelActor;

        vtkEventQtSlotConnect* Connections;
        vtkEventQtSlotConnect* Connections_s;

        unsigned int newHeightValue;

        void Draw();

        vtkDoubleArray* InterpSurface(vtkPolyData* pPointData, vtkPolyData* pSurface, int pInterpDistance);
        void ComputeFibreDirection();

        vtkPoints* AddPoint(vtkPoints* array, double* p);
        vtkDoubleArray* AddScalar(vtkDataArray* array, double v);


};

#endif
