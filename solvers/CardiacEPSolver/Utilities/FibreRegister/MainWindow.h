#ifndef CLASS_MAINWINDOW_H
#define CLASS_MAINWINDOW_H

#include <vector>
#include <string>
#include <sstream>

//#include <QtGui/QMainWindow>
//#include <QtGui/QLineEdit>
//#include <QtGui/QGridLayout>
#include <QtGui>
#include <QVTKWidget.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataMapper.h>
#include <vtkPointData.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkDepthSortPolyData.h>
#include <vtkGlyph3D.h>
#include <vtkSphereSource.h>
#include <vtkWorldPointPicker.h>
//#include <vtkLookupTable.h>
//#include <vtkColorTransferFunction.h>
//#include <vtkScalarBarActor.h>
//#include <vtkTextMapper.h>
//#include <vtkActor2D.h>

class vtkEventQtSlotConnect;


class MainWindow : public QMainWindow
{
    Q_OBJECT

    public:
        MainWindow( QWidget* parent = 0, Qt::WindowFlags fl = 0 );
        ~MainWindow();

    private slots:
        void BrowseLeft();
        void BrowseRight();
        void Load();
        void Update();
        void CreateLeftPoint(vtkObject*, unsigned long, void*, void*, vtkCommand*);

    private:
        // GUI widgets
        QWidget* mRootWidget;
        QGridLayout* mRootGrid;
        QGroupBox* mSettingsBox;
        QGridLayout* mSettingsGrid;
        QGroupBox* mFileBox;
        QGridLayout* mFileGrid;
        QLineEdit* mFileLeftEditBox;
        QPushButton* mFileLeftBrowse;
        QLineEdit* mFileRightEditBox;
        QPushButton* mFileRightBrowse;
        QPushButton* mFileLoadButton;
        QVTKWidget* mLeftVtk;
        QVTKWidget* mRightVtk;

        // VTK
        // Left Surface pipeline
        vtkPolyData* mLeftData;
        vtkSmoothPolyDataFilter* mLeftFilterSmooth;
        vtkDepthSortPolyData* mLeftFilterDepthSort;
        vtkPolyDataMapper* mLeftMapper;
        vtkActor* mLeftActor;
        vtkRenderer* mLeftRenderer;

        // Right Surface pipeline
        vtkPolyData* mRightData;
        vtkSmoothPolyDataFilter* mRightFilterSmooth;
        vtkDepthSortPolyData* mRightFilterDepthSort;
        vtkPolyDataMapper* mRightMapper;
        vtkActor* mRightActor;
        vtkRenderer* mRightRenderer;

        // Left Points pipeline
        vtkPolyData* mLeftPointsData;
        vtkSphereSource* mLeftSphere;
        vtkGlyph3D* mLeftFilterGlyph;
        vtkPolyDataMapper* mLeftPointsMapper;
        vtkActor* mLeftPointsActor;

        // Right Points pipeline
        vtkPolyData* mRightPointsData;
        vtkSphereSource* mRightSphere;
        vtkGlyph3D* mRightFilterGlyph;
        vtkPolyDataMapper* mRightPointsMapper;
        vtkActor* mRightPointsActor;

        vtkEventQtSlotConnect* Connections;

        void Draw();


};

#endif
