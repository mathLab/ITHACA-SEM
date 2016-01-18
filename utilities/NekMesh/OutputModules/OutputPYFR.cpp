////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputGmsh.cpp
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
//  Description: Gmsh file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/MeshElements/MeshElements.h>

#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#include "OutputPYFR.h"

using namespace std;
using namespace Nektar::NekMeshUtils;

namespace Nektar
{
namespace Utilities
{

ModuleKey OutputPYFR::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eOutputModule, "pyfrm"), OutputPYFR::create,
        "Writes PYFR pyfrm file.");

OutputPYFR::OutputPYFR(MeshSharedPtr m) : OutputModule(m)
{

}

OutputPYFR::~OutputPYFR()
{

}

/**
 * @brief Process a mesh to output to Gmsh MSH format.
 *
 * Gmsh output is fairly straightforward. The file first contains a
 * list of nodes, followed by a list of elements. Since
 * Mesh::vertexSet only contains vertices of the linear elements, we
 * first loop over the elements so that any high-order vertices can be
 * enumerated and then added to the node list. We then print out the
 * list of nodes and finally print the element list.
 */
void OutputPYFR::Process()
{
    if (m_mesh->m_verbose)
    {
        cout << "OutputPYFR: Writing file..." << endl;
    }

    string filename = m_config["outfile"].as<string>();

    H5std_string FILE_NAME( filename );
    H5File* file = new H5File( filename, H5F_ACC_TRUNC );

    CompositeMap cm = m_mesh->m_composite;
    CompositeMap::iterator it;
    for(it = cm.begin(); it != cm.end(); it++)
    {
        cout << it->first << " " << it->second->m_tag << endl;

        int nv=0;
        map<int, int> nodemap;
        string dsname;
        if(it->second->m_tag == "Q")
        {
            nv = 4;
            dsname = "spt_quad_p0";
            nodemap[0] = 0;
            nodemap[1] = 1;
            nodemap[2] = 3;
            nodemap[3] = 2;
        }
        else if(it->second->m_tag == "T")
        {
            nv = 3;
            dsname = "spt_tri_p0";
            nodemap[0] = 0;
            nodemap[1] = 1;
            nodemap[2] = 2;
        }
        else
        {
            continue;
        }

        hsize_t dimsf[] = {nv, it->second->m_items.size(), m_mesh->m_expDim};

        DataSpace dataspace( 3, dimsf );

        DataSet* dataset = new DataSet(file->createDataSet(dsname, PredType::NATIVE_DOUBLE, dataspace));

        double* data = new double[nv*it->second->m_items.size()*m_mesh->m_expDim];

        for(int i = 0; i < it->second->m_items.size(); i++)
        {
            vector<NodeSharedPtr> ns = it->second->m_items[i]->GetVertexList();
            for(int j = 0; j < nv; j++)
            {
                data[j*it->second->m_items.size()*m_mesh->m_expDim + i*m_mesh->m_expDim + 0] = ns[nodemap[j]]->m_x;
                data[j*it->second->m_items.size()*m_mesh->m_expDim + i*m_mesh->m_expDim + 1] = ns[nodemap[j]]->m_y;
            }
        }

        dataset->write( data, PredType::NATIVE_DOUBLE );
    }

    typedef struct conec
    {
        char      el[5];
        int       id;
        int       fc;
        int       bl;
    } conec;
    StrType strdatatype(PredType::C_S1, 5);

    { //con
        EdgeSet::iterator eit;
        EdgeSet interiorcons;
        for(eit = m_mesh->m_edgeSet.begin(); eit != m_mesh->m_edgeSet.end(); eit++)
        {
            ASSERTL0((*eit)->m_elLink.size() == 2 || (*eit)->m_elLink.size() == 1, "not enough element links");

            if((*eit)->m_elLink.size() == 2)
            {
                interiorcons.insert(*eit);
            }
        }

        hsize_t dimsf[] = {2, interiorcons.size()};

        DataSpace dataspace( 2, dimsf );

        CompType cn( sizeof(conec) );
        cn.insertMember( "f0", HOFFSET(conec, el), strdatatype);
        cn.insertMember( "f1", HOFFSET(conec, id), PredType::NATIVE_INT32);
        cn.insertMember( "f2", HOFFSET(conec, fc), PredType::NATIVE_INT8);
        cn.insertMember( "f3", HOFFSET(conec, bl), PredType::NATIVE_INT8);

        DataSet* dataset = new DataSet(file->createDataSet("con_p0", cn, dataspace));

        conec* cons = new conec[2*interiorcons.size()];

        for(eit = interiorcons.begin(); eit != interiorcons.end(); eit++)
        {
            
        }
    }


}

}
}
