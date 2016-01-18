////////////////////////////////////////////////////////////////////////////////
//
//  File: InputNekpp.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
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
//  Description: GMSH converter.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/algorithm/string.hpp>

#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#include "InputPYFR.h"

using namespace std;

namespace Nektar
{
namespace Utilities
{

ModuleKey InputPYFR::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eInputModule, "pyfrm"), InputPYFR::create,
        "PYFR mesh files.");

/**
 * @brief Set up InputNekpp object.
 */
InputPYFR::InputPYFR(MeshSharedPtr m) : InputModule(m)
{

}

InputPYFR::~InputPYFR()
{

}

void InputPYFR::Process()
{
    string filename = m_config["infile"].as<string>();

    H5std_string FILE_NAME( filename );
    H5File file( FILE_NAME, H5F_ACC_RDONLY );

    vector<string> elements;
    vector<string> bcs;
    string con, uuid;

    Group g = file.openGroup("/");
    for(int i = 0; i < g.getNumObjs(); i++)
    {
        string obj = g.getObjnameByIdx(i);
        string pre = obj.substr(0, obj.find("_"));

        if(pre == "bcon")
        {
            bcs.push_back(obj);
        }
        else if(pre == "spt")
        {
            elements.push_back(obj);
        }
        else if(pre == "mesh")
        {
            uuid = obj;
        }
        else if(pre == "con")
        {
            con = obj;
        }
        else
        {
            ASSERTL0(false, "hdf5 object not recognised");
        }
    }

    m_mesh->m_expDim = 2;
    m_mesh->m_spaceDim = 2;
    m_mesh->m_nummode = 2;

    map<string, vector<ElementSharedPtr> > elmap;

    int nct = 0;
    for(int i = 0; i < elements.size(); i++)
    {
        // all this is to just get the data out of hdf5, why so complicated!!!
        DataSet d = g.openDataSet(elements[i]);
        DataSpace s = d.getSpace();
        int rank = s.getSimpleExtentNdims();
        hsize_t* dims_out = new hsize_t[rank];
        int ndims = s.getSimpleExtentDims( dims_out, NULL);
        hsize_t* offset = new hsize_t[ndims];	// hyperslab offset in the file
        memset(offset, 0, rank * sizeof(hsize_t)) ;
        s.selectHyperslab( H5S_SELECT_SET, dims_out, offset);
        DataSpace memspace( ndims, dims_out); //describe hyperslab in memory space
        memspace.selectHyperslab( H5S_SELECT_SET, dims_out, offset);
        hsize_t totSize = 1 ;
        for(int j = 0 ; j < ndims ; j++)
        {
            totSize *= dims_out[j];
        }
        float* data = new float[totSize] ;
        d.read( data, PredType::NATIVE_FLOAT, memspace, s );

        Array<ThreeD, NekDouble> dataArray(dims_out[0],dims_out[1],dims_out[2]);
        for(int j = 0; j < dims_out[0]; j++)
        {
            for(int k = 0; k < dims_out[1]; k++)
            {
                for(int l = 0; l < dims_out[2]; l++)
                {
                    dataArray[j][k][l] = data[j*dims_out[1]*dims_out[2] + k*dims_out[2] + l];
                }
            }
        }

        vector<string> tmp;
        boost::split(tmp, elements[i], boost::is_any_of("_"));

        LibUtilities::ShapeType st;
        if(boost::iequals(tmp[1], "quad"))
        {
            st = LibUtilities::eQuadrilateral;
        }
        else if(boost::iequals(tmp[1], "tri"))
        {
            st = LibUtilities::eTriangle;
        }
        else
        {
            ASSERTL0(false,"element type not recognised");
        }

        map<int, int> nodemap; //nektar order to pyfr
        int nv;
        if(st == LibUtilities::eQuadrilateral)
        {
            nodemap[0] = 0;
            nodemap[1] = 1;
            nodemap[2] = 3;
            nodemap[3] = 2;
            nv = 4;
        }
        else if(st == LibUtilities::eTriangle)
        {
            nodemap[0] = 0;
            nodemap[1] = 1;
            nodemap[2] = 2;
            nv = 3;
        }

        for(int j = 0; j < dims_out[1]; j++) //loop over all elements
        {
            vector<NodeSharedPtr> ns(nv);
            for(int k = 0; k < nv; k++)
            {
                ns[k] = boost::shared_ptr<Node>(new Node(nct++, dataArray[nodemap[k]][j][0],
                                                dataArray[nodemap[k]][j][1], 0.0));
            }
            vector<int> tags;
            tags.push_back(i);
            ElmtConfig conf(st, 1, false, false);
            ElementSharedPtr E = GetElementFactory().CreateInstance(
                                    st, conf, ns, tags);
            m_mesh->m_element[m_mesh->m_expDim].push_back(E);
            elmap[tmp[1]].push_back(E);
        }
    }

    typedef struct conec
    {
        char      el[5];
        int       id;
        int       fc;
        int       bl;
    } conec;
    StrType strdatatype(PredType::C_S1, 5);

    //con
    {
        // all this is to just get the data out of hdf5, why so complicated!!!
        DataSet d = g.openDataSet(con);
        DataSpace s = d.getSpace();
        int rank = s.getSimpleExtentNdims();
        hsize_t* dims_out = new hsize_t[rank];
        int ndims = s.getSimpleExtentDims( dims_out, NULL);
        hsize_t* offset = new hsize_t[ndims];	// hyperslab offset in the file
        memset(offset, 0, rank * sizeof(hsize_t)) ;
        s.selectHyperslab( H5S_SELECT_SET, dims_out, offset);
        DataSpace memspace( ndims, dims_out); //describe hyperslab in memory space
        memspace.selectHyperslab( H5S_SELECT_SET, dims_out, offset);

        CompType cn( sizeof(conec) );
        cn.insertMember( "f0", HOFFSET(conec, el), strdatatype);
        cn.insertMember( "f1", HOFFSET(conec, id), PredType::NATIVE_INT);
        cn.insertMember( "f2", HOFFSET(conec, fc), PredType::NATIVE_INT);
        cn.insertMember( "f3", HOFFSET(conec, bl), PredType::NATIVE_INT);

        hsize_t totSize = 1 ;
        for(int j = 0 ; j < ndims ; j++)
        {
            totSize *= dims_out[j];
        }
        conec* data = new conec[totSize];
        d.read( data, cn );

        for(int k = 0; k < dims_out[1]; k++)
        {
            ElementSharedPtr e1 = elmap[data[0*dims_out[1] + k].el][data[0*dims_out[1] + k].id];
            ElementSharedPtr e2 = elmap[data[1*dims_out[1] + k].el][data[1*dims_out[1] + k].id];

            e2->SetEdge(data[1*dims_out[1] + k].fc, e1->GetEdge(data[0*dims_out[1] + k].fc));
        }

    }
    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();

    for(int i = 0; i < bcs.size(); i++)
    {
        // all this is to just get the data out of hdf5, why so complicated!!!
        DataSet d = g.openDataSet(bcs[i]);
        DataSpace s = d.getSpace();
        int rank = s.getSimpleExtentNdims();
        hsize_t* dims_out = new hsize_t[rank];
        int ndims = s.getSimpleExtentDims( dims_out, NULL);
        hsize_t* offset = new hsize_t[ndims];	// hyperslab offset in the file
        memset(offset, 0, rank * sizeof(hsize_t)) ;
        s.selectHyperslab( H5S_SELECT_SET, dims_out, offset);
        DataSpace memspace( ndims, dims_out); //describe hyperslab in memory space
        memspace.selectHyperslab( H5S_SELECT_SET, dims_out, offset);

        CompType cn( sizeof(conec) );
        cn.insertMember( "f0", HOFFSET(conec, el), strdatatype);
        cn.insertMember( "f1", HOFFSET(conec, id), PredType::NATIVE_INT);
        cn.insertMember( "f2", HOFFSET(conec, fc), PredType::NATIVE_INT);
        cn.insertMember( "f3", HOFFSET(conec, bl), PredType::NATIVE_INT);

        hsize_t totSize = 1 ;
        for(int j = 0 ; j < ndims ; j++)
        {
            totSize *= dims_out[j];
        }

        conec* data = new conec[totSize];
        d.read( data, cn );

        for(int k = 0; k < dims_out[0]; k++)
        {
            ElementSharedPtr e1 = elmap[data[k].el][data[k].id];
            EdgeSharedPtr e = e1->GetEdge(data[k].fc);

            vector<NodeSharedPtr> ns(2);
            ns[0] = e->m_n1;
            ns[1] = e->m_n2;

            vector<int> tags;
            tags.push_back(elements.size()+i);
            ElmtConfig conf(LibUtilities::eSegment, 1, false, false);
            ElementSharedPtr E = GetElementFactory().CreateInstance(
                                    LibUtilities::eSegment, conf, ns, tags);
            m_mesh->m_element[m_mesh->m_expDim-1].push_back(E);
        }
    }

    ProcessVertices();
    ProcessEdges();
    ProcessFaces();
    ProcessElements();
    ProcessComposites();
}

}
}
