////////////////////////////////////////////////////////////////////////////////
//
//  File: Opencascade.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: occ headers.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_CADSYSTEM_OCC
#define NEKMESHUTILS_CADSYSTEM_OCC

/// This is a list of OpenCascade headers required for use with nektar

/// IO classes
#include <STEPCAFControl_Reader.hxx>
#include <StepRepr_RepresentationItem.hxx>
#include <TDocStd_Document.hxx>
#include <XSControl_WorkSession.hxx>
#include <XSControl_TransferReader.hxx>
#include <XCAFDoc_DocumentTool.hxx>
#include <Storage.hxx>

/// STL classes
#include <BRepMesh_IncrementalMesh.hxx>


/// Shape Analysis / exploration classes
#include <BRepTopAdaptor_FClass2d.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <TopExp.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <TopExp_Explorer.hxx>
#include <GeomAdaptor_HSurface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <GeomLProp_CLProps.hxx>
#include <GeomLProp_SLProps.hxx>
#include <ShapeAnalysis_Curve.hxx>
#include <BRepBndLib.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <GCPnts_AbscissaPoint.hxx>

/// Shape fixing classes
#include <ShapeFix_Face.hxx>

/// Shape Building classes
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <Geom_BSplineCurve.hxx>
#include <GeomAPI_PointsToBSpline.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <Geom_TrimmedCurve.hxx>

/// Data structure classes
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_DataMapOfShapeShape.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <Transfer_Binder.hxx>
#include <TransferBRep.hxx>
#include <Transfer_TransientProcess.hxx>
#include <Interface_InterfaceModel.hxx>
#include <TCollection_HAsciiString.hxx>

/// CORE SHAPE classes
#include <TopoDS_Shape.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS.hxx>

/// GP clasases
#include <gp_Pnt.hxx>
#include <gp_Ax1.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Trsf.hxx>

#include <Geom_BSplineCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Ellipse.hxx>
#include <gp_Circ.hxx>
#include <gp_Elips.hxx>

#endif
