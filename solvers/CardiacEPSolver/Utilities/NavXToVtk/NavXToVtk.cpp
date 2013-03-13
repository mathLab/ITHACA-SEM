/*
 * NavXToVtk.cpp
 *
 *  Created on: 19 Jan 2013
 *      Author: cc
 */

#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>

#include "tinyxml/tinyxml.h"

#include <vector>
#include <map>
#include <string>
#include <sstream>
using namespace std;

vtkPolyData* LoadGeoFile(string nxgeo);
vtkPolyData* LoadDataFile(string nxdxl);
std::vector<std::string> Split(std::string& input);

int main(int argc, char* argv[])
{
    if (argc != 5)
    {
        cout << "Usage: NavXToVtk DxLandmarkGeo.xml DxL_1.csv geometry.vtk points.vtk" << endl;
        exit(-1);
    }
    string nxgeo = argv[1];
    string nxdxl = argv[2];
    string outgeo = argv[3];
    string outdxl = argv[4];

    cout << "NavX Geometry: " << nxgeo << endl;
    cout << "NavX DxL Data: " << nxdxl << endl;
    cout << "VTK Output:    " << outgeo << endl;
    cout << "VTK Output:    " << outdxl << endl;

    vtkPolyData* mesh = LoadGeoFile(nxgeo);
    vtkPolyData* points = LoadDataFile(nxdxl);

    if (!mesh) {
        cout << "Error loading geometry file." << endl;
        exit(-1);
    }
    if (!points) {
        cout << "Error loading DxL file." << endl;
        exit(-1);
    }

    cout << "Writing VTK geometry file " << outgeo << endl;
    vtkPolyDataWriter *writergeo = vtkPolyDataWriter::New();
    writergeo->SetInput(mesh);
    writergeo->SetFileName(outgeo.c_str());
    writergeo->Write();

    cout << "Writing VTK electrogram points file " << outdxl << endl;
    vtkPolyDataWriter *writerdxl = vtkPolyDataWriter::New();
    writerdxl->SetInput(points);
    writerdxl->SetFileName(outdxl.c_str());
    writerdxl->Write();

}


vtkPolyData* LoadGeoFile(string nxgeo) {
    double x, y, z;
    TiXmlDocument doc;
    bool ok = doc.LoadFile(nxgeo.c_str());
    if (!ok) {
        cout << "Failed to load Geo file." << endl;
        return 0;
    }
    vtkPolyData* mGeometryData = vtkPolyData::New();
    mGeometryData->Reset();

    vtkPoints* vPts = vtkPoints::New();
    vtkCellArray* vCells = vtkCellArray::New();

    TiXmlElement* dif = doc.FirstChildElement("DIF");
    if (!dif) {cout << "Error." << endl; return 0;}
    TiXmlElement* difbody = dif->FirstChildElement("DIFBody");
    if (!difbody) {cout << "Error." << endl; return 0;}
    TiXmlElement* volumes = difbody->FirstChildElement("Volumes");
    if (!volumes) {cout << "Error." << endl; return 0;}
    TiXmlElement* volume = volumes->FirstChildElement("Volume");
    while (volume) {
        TiXmlElement* vertices = volume->FirstChildElement("Vertices");
        if (!vertices) {cout << "Error." << endl; return 0;};
        stringstream s;
        s.str(vertices->GetText());
        while (s) {
            s >> x;
            s >> y;
            s >> z;
            if (!s) break;
            vPts->InsertNextPoint(x,y,z);
        }

        TiXmlElement* polygons = volume->FirstChildElement("Polygons");
        vtkIdType pts[3];
        if (!polygons) {cout << "Error." << endl; return 0;}
        s.clear();
        s.str(polygons->GetText());
        while (s) {
            s >> pts[0]; pts[0]--;
            s >> pts[1]; pts[1]--;
            s >> pts[2]; pts[2]--;
            if (!s) break;
            vCells->InsertNextCell(3, pts);
        }
        mGeometryData->SetPoints(vPts);
        mGeometryData->SetPolys(vCells);
        volume = volume->NextSiblingElement("Volume");
    }
    return mGeometryData;
}

vtkPolyData* LoadDataFile(string nxdxl) {
    vtkPoints* vPts = vtkPoints::New();
    ifstream f(nxdxl.c_str());
    string input;
    if (!f.good()) {
        cout << "Data file cannot be opened." << endl;
        return 0;
    }
    vtkPolyData* mPointsData = vtkPolyData::New();
    mPointsData->Reset();

    std::map<std::string, std::vector<std::string> > values;
    std::vector<std::string> tmp;
    while (!f.eof()) {
        getline(f,input);
        if (input.length() >= 10 && input.substr(0,10) == "Begin data") {
            break;
        }
    }

    while (!f.eof()) {
        getline(f,input);
        if (input.length() > 16 && input.substr(0,16) == "Exported seconds") {
            break;
        }
        tmp = Split(input);
        string name = tmp[0].substr(0,tmp[0].find_first_of(':'));
        tmp.erase(tmp.begin());
        values[name] = tmp;
    }

    vector<string> labels;
    while (!f.eof()) {
        getline(f, input);
        if (input.length() > 9 && input.substr(0,9) == "rov trace") {
            labels = Split(input);
            for (int i = 1; i < labels.size(); ++i) {
                vPts->InsertNextPoint(
                                        atof(values["roving x"][i-1].c_str()),
                                        atof(values["roving y"][i-1].c_str()),
                                        atof(values["roving z"][i-1].c_str()));
            }
            break;
        }
    }
    mPointsData->SetPoints(vPts);

    return mPointsData;
}

std::vector<std::string> Split(std::string& input) {
    std::istringstream s(input);
    std::string tmp;
    std::vector<std::string> output;
    while (s) {
        if (!getline(s, tmp, ',')) {
            break;
        }
        output.push_back(tmp);
    }
    return output;
}
