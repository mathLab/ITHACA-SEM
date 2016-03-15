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

#include <NekMeshUtils/MeshElements/Element.h>

#include <LibUtilities/Foundations/ManagerAccess.h>

#include <H5Cpp.h>
#ifndef H5_NO_NAMESPACE
    using namespace H5;
#endif

#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators
#include <boost/uuid/uuid_io.hpp>         // streaming operators etc.

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

vector<Array<OneD, NekDouble> > GetPyFrStandardTet(int n)
{
    vector<Array<OneD, NekDouble> > ret;

    NekDouble du = 2.0 / (n-1);

    for(int z = 0; z < n; z++)
    {
        for(int y = 0; y < n - z; y++)
        {
            for(int x = 0; x < n - y - z; x++)
            {
                Array<OneD, NekDouble> uvw(3);
                uvw[0] = -1.0 + x * du;
                uvw[1] = -1.0 + y * du;
                uvw[2] = -1.0 + z * du;
                ret.push_back(uvw);
            }
        }
    }

    return ret;
}

vector<Array<OneD, NekDouble> > GetPyFrStandardPri(int n)
{
    vector<Array<OneD, NekDouble> > ret;

    NekDouble du = 2.0 / (n-1);

    for(int z = 0; z < n; z++)
    {
        int ct = 0;
        for(int y = 0; y < n; y++)
        {
            for(int x = 0; x < n - y; x++)
            {
                Array<OneD, NekDouble> uvw(3);
                uvw[0] = -1.0 + x * du;
                uvw[1] = -1.0 + y * du;
                uvw[2] = -1.0 + z * du;
                ret.push_back(uvw);
            }
        }
    }

    return ret;
}

vector<Array<OneD, NekDouble> > GetPyFRNodes(ElementSharedPtr el, int nm)
{
    vector<Array<OneD, NekDouble> > ret;

    if(el->GetConf().m_e == LibUtilities::eTetrahedron)
    {
        vector<Array<OneD, NekDouble> > Tet = GetPyFrStandardTet(nm);

        SpatialDomains::GeometrySharedPtr geom = el->GetGeom(3);
        geom->FillGeom();
        StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();
        Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
        Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
        Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

        Array<OneD, NekDouble> xc(xmap->GetTotPoints());
        Array<OneD, NekDouble> yc(xmap->GetTotPoints());
        Array<OneD, NekDouble> zc(xmap->GetTotPoints());
        xmap->BwdTrans(coeffs0,xc);
        xmap->BwdTrans(coeffs1,yc);
        xmap->BwdTrans(coeffs2,zc);

        for(int i = 0; i < Tet.size(); i++)
        {
            Array<OneD, NekDouble> xp = Tet[i];

            Array<OneD, NekDouble> xyz(3);
            xyz[0] = xmap->PhysEvaluate(xp, xc);
            xyz[1] = xmap->PhysEvaluate(xp, yc);
            xyz[2] = xmap->PhysEvaluate(xp, zc);
            ret.push_back(xyz);
        }

    }
    else if(el->GetConf().m_e == LibUtilities::ePrism)
    {
        vector<Array<OneD, NekDouble> > Pri = GetPyFrStandardPri(nm);

        SpatialDomains::GeometrySharedPtr geom = el->GetGeom(3);
        geom->FillGeom();
        StdRegions::StdExpansionSharedPtr xmap = geom->GetXmap();
        Array<OneD, NekDouble> coeffs0 = geom->GetCoeffs(0);
        Array<OneD, NekDouble> coeffs1 = geom->GetCoeffs(1);
        Array<OneD, NekDouble> coeffs2 = geom->GetCoeffs(2);

        Array<OneD, NekDouble> xc(xmap->GetTotPoints());
        Array<OneD, NekDouble> yc(xmap->GetTotPoints());
        Array<OneD, NekDouble> zc(xmap->GetTotPoints());
        xmap->BwdTrans(coeffs0,xc);
        xmap->BwdTrans(coeffs1,yc);
        xmap->BwdTrans(coeffs2,zc);

        for(int i = 0; i < Pri.size(); i++)
        {
            Array<OneD, NekDouble> xp = Pri[i];

            Array<OneD, NekDouble> xyz(3);
            xyz[0] = xmap->PhysEvaluate(xp, xc);
            xyz[1] = xmap->PhysEvaluate(xp, yc);
            xyz[2] = xmap->PhysEvaluate(xp, zc);
            ret.push_back(xyz);
        }
    }

    return ret;
}

string GetDSName(string n)
{
    string ret = "Null";
    if(n == "A")
    {
        ret = "spt_tet_p0";
    }
    else if(n == "R")
    {
        ret = "spt_pri_p0";
    }
    else
    {
        cout << "unknown type" << endl;
        abort();
    }

    return ret;
}

int GetNump(string n, int nm)
{
    int ret = 0;
    if(n == "A")
    {
        LibUtilities::PointsKey pkey(nm, LibUtilities::eNodalTetEvenlySpaced);

        ret = LibUtilities::PointsManager()[pkey]->GetTotNumPoints();
    }
    else if(n == "R")
    {
        LibUtilities::PointsKey pkey(nm, LibUtilities::eNodalPrismEvenlySpaced);

        ret = LibUtilities::PointsManager()[pkey]->GetTotNumPoints();
    }
    return ret;
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

    int numtri = 0;
    int numquad = 0;

    int numhex = 0;
    int numpri = 0;
    int numpry = 0;
    int numtet = 0;

    map<int, int> nekidtopyid;

    CompositeMap cm = m_mesh->m_composite;
    CompositeMap::iterator it;
    for(it = cm.begin(); it != cm.end(); it++)
    {
        cout << it->first << " " << it->second->m_tag;

        if(it->second->m_tag == "F" || it->second->m_tag == "E")
        {
            //not interested in boundary faces or edges here
            continue;
        }

        string dsname = GetDSName(it->second->m_tag);

        int np = GetNump(it->second->m_tag, m_mesh->m_nummode);
        cout << " points: " << np << endl;

        hsize_t dimsf[] = {np, it->second->m_items.size(), m_mesh->m_expDim};

        DataSpace dataspace( 3, dimsf );

        DataSet* dataset = new DataSet(file->createDataSet(dsname, PredType::NATIVE_DOUBLE, dataspace));

        double* data = new double[np * it->second->m_items.size() * m_mesh->m_expDim];

        for(int i = 0; i < it->second->m_items.size(); i++)
        {
            if(it->second->m_tag == "T")
            {
                nekidtopyid[it->second->m_items[i]->GetId()] = numtri++;
            }
            else if(it->second->m_tag == "Q")
            {
                nekidtopyid[it->second->m_items[i]->GetId()] = numquad++;
            }
            else if(it->second->m_tag == "H")
            {
                nekidtopyid[it->second->m_items[i]->GetId()] = numhex++;
            }
            else if(it->second->m_tag == "R")
            {
                nekidtopyid[it->second->m_items[i]->GetId()] = numpri++;
            }
            else if(it->second->m_tag == "P")
            {
                nekidtopyid[it->second->m_items[i]->GetId()] = numpry++;
            }
            else if(it->second->m_tag == "A")
            {
                nekidtopyid[it->second->m_items[i]->GetId()] = numtet++;
            }

            vector<Array<OneD, NekDouble> > ns = GetPyFRNodes(it->second->m_items[i], m_mesh->m_nummode);

            for(int j = 0; j < ns.size(); j++)
            {
                int offset = j*it->second->m_items.size()*m_mesh->m_expDim + i*m_mesh->m_expDim;
                data[offset + 0] = ns[j][0];
                data[offset + 1] = ns[j][1];
                if(m_mesh->m_expDim == 3)
                {
                    data[offset + 2] = ns[j][2];
                }
            }
        }

        dataset->write( data, PredType::NATIVE_DOUBLE );
    }

    typedef struct conec
    {
        char      el[4];
        int       id;
        int       fc;
        int       bl;
    } conec;
    StrType strdatatype(PredType::C_S1, 4);

    //map from nektar face to pyfr face
    map<string,map<int,int> > linkmap;
    map<int,int> lm;
    lm[0] = 0;
    lm[1] = 1;
    lm[2] = 3;
    lm[3] = 2;
    linkmap["tet"] = lm;
    lm.clear();
    lm[0] = 2;
    lm[1] = 1;
    lm[2] = 3;
    lm[3] = 0;
    lm[4] = 4;
    linkmap["pri"] = lm;

    { //con
        if(m_mesh->m_expDim == 2)
        {
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

            int ct = 0;
            for(eit = interiorcons.begin(); eit != interiorcons.end(); eit++)
            {
                conec c1,c2;
                ElementSharedPtr e1 = (*eit)->m_elLink[0].first;
                ElementSharedPtr e2 = (*eit)->m_elLink[1].first;

                string str1, str2;

                if(e1->GetConf().m_e == LibUtilities::eTriangle)
                {
                    str1 = "tri";
                }
                else if(e1->GetConf().m_e == LibUtilities::eQuadrilateral)
                {
                    str1 = "quad";
                }

                if(e2->GetConf().m_e == LibUtilities::eTriangle)
                {
                    str2 = "tri";
                }
                else if(e2->GetConf().m_e == LibUtilities::eQuadrilateral)
                {
                    str2 = "quad";
                }

                strcpy(c1.el, str1.c_str());
                strcpy(c2.el, str2.c_str());

                c1.id = nekidtopyid[e1->GetId()];
                c2.id = nekidtopyid[e2->GetId()];

                c1.fc = (*eit)->m_elLink[0].second;
                c2.fc = (*eit)->m_elLink[1].second;

                c1.bl = 0; c2.bl = 0;

                cons[ct] = c1;
                cons[ct + interiorcons.size()] = c2;
                ct++;
            }

            dataset->write( cons, cn );
        }
        else if(m_mesh->m_expDim == 3)
        {
            FaceSet::iterator fit;
            FaceSet interiorcons;
            for(fit = m_mesh->m_faceSet.begin(); fit != m_mesh->m_faceSet.end(); fit++)
            {
                ASSERTL0((*fit)->m_elLink.size() == 2 || (*fit)->m_elLink.size() == 1, "not enough element links");

                if((*fit)->m_elLink.size() == 2)
                {
                    interiorcons.insert(*fit);
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

            int ct = 0;
            for(fit = interiorcons.begin(); fit != interiorcons.end(); fit++)
            {
                conec c1,c2;
                ElementSharedPtr e1 = (*fit)->m_elLink[0].first;
                ElementSharedPtr e2 = (*fit)->m_elLink[1].first;

                string str1, str2;

                if(e1->GetConf().m_e == LibUtilities::eTetrahedron)
                {
                    str1 = "tet";
                }
                else if(e1->GetConf().m_e == LibUtilities::ePrism)
                {
                    str1 = "pri";
                }

                if(e2->GetConf().m_e == LibUtilities::eTetrahedron)
                {
                    str2 = "tet";
                }
                else if(e2->GetConf().m_e == LibUtilities::ePrism)
                {
                    str2 = "pri";
                }

                map<string,map<int, int> >::iterator lmap1 = linkmap.find(str1);
                map<string,map<int, int> >::iterator lmap2 = linkmap.find(str2);
                ASSERTL0(lmap1 != linkmap.end() && lmap2 != linkmap.end(),"unsupported element");

                strcpy(c1.el, str1.c_str());
                strcpy(c2.el, str2.c_str());

                c1.id = nekidtopyid[e1->GetId()];
                c2.id = nekidtopyid[e2->GetId()];

                c1.fc = lmap1->second[(*fit)->m_elLink[0].second];
                c2.fc = lmap2->second[(*fit)->m_elLink[1].second];

                c1.bl = 0; c2.bl = 0;

                cons[ct] = c1;
                cons[ct + interiorcons.size()] = c2;
                ct++;
            }

            dataset->write( cons, cn );
        }

    }

    StrType strdatatypevar(PredType::C_S1, H5T_VARIABLE);
    { //uuid
        hsize_t dimsf[] = {1};

        DataSpace dataspace( 1, dimsf );

        DataSet* dataset = new DataSet(file->createDataSet("mesh_uuid", strdatatypevar, dataspace));

        boost::uuids::uuid uuid = boost::uuids::random_generator()();
        string uuidstr =  boost::uuids::to_string(uuid);

        dataset->write( H5std_string(uuidstr), strdatatypevar );
    }

    //bcons
    for(it = cm.begin(); it != cm.end(); it++)
    {
        string dsname;
        if(m_mesh->m_expDim == 2)
        {
            if(it->second->m_tag == "E")
            {
                stringstream ss;
                ss << "bcon_C" << it->first << "_p0";
                dsname = ss.str();
            }
            else
            {
                continue;
            }
        }
        else if(m_mesh->m_expDim == 3)
        {
            if(it->second->m_tag == "F")
            {
                stringstream ss;
                ss << "bcon_C" << it->first << "_p0";
                dsname = ss.str();
            }
            else
            {
                continue;
            }
        }

        hsize_t dimsf[] = {it->second->m_items.size()};

        DataSpace dataspace( 1, dimsf );

        CompType cn( sizeof(conec) );
        cn.insertMember( "f0", HOFFSET(conec, el), strdatatype);
        cn.insertMember( "f1", HOFFSET(conec, id), PredType::NATIVE_INT32);
        cn.insertMember( "f2", HOFFSET(conec, fc), PredType::NATIVE_INT8);
        cn.insertMember( "f3", HOFFSET(conec, bl), PredType::NATIVE_INT8);

        DataSet* dataset = new DataSet(file->createDataSet(dsname, cn, dataspace));

        conec* cons = new conec[it->second->m_items.size()];

        for(int i = 0; i < it->second->m_items.size(); i++)
        {
            pair<ElementSharedPtr,int> el;
            if(m_mesh->m_expDim ==2 )
            {
                ASSERTL0(it->second->m_items[i]->GetEdgeLink()->m_elLink.size()==1,"el link makes no sense");
                el = it->second->m_items[i]->GetEdgeLink()->m_elLink[0];
            }
            else if(m_mesh->m_expDim == 3)
            {
                ASSERTL0(it->second->m_items[i]->GetFaceLink()->m_elLink.size()==1,"el link makes no sense");
                el = it->second->m_items[i]->GetFaceLink()->m_elLink[0];
            }

            conec c;

            string str;

            if(el.first->GetConf().m_e == LibUtilities::eTriangle)
            {
                str = "tri";
            }
            else if(el.first->GetConf().m_e == LibUtilities::eQuadrilateral)
            {
                str = "quad";
            }
            else if(el.first->GetConf().m_e == LibUtilities::eTetrahedron)
            {
                str = "tet";
            }
            else if(el.first->GetConf().m_e == LibUtilities::ePrism)
            {
                str = "pri";
            }

            strcpy(c.el, str.c_str());

            map<string,map<int, int> >::iterator lmap = linkmap.find(str);
            ASSERTL0(lmap != linkmap.end() ,"unsupported element");

            c.id = nekidtopyid[el.first->GetId()];
            c.fc = lmap->second[el.second];
            c.bl = 0;

            cons[i] = c;
        }

        dataset->write( cons, cn );
    }




}

}
}
