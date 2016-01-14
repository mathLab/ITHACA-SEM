////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessIsoContour.cpp
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
//  Description: Generate isocontours from field data.
//
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <iostream>
using namespace std;

#include "ProcessIsoContour.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessIsoContour::className =
    GetModuleFactory().RegisterCreatorFunction(
                        ModuleKey(eProcessModule, "isocontour"),
                        ProcessIsoContour::create,
                        "Extract an isocontour of fieldid variable and at "
                        "value fieldvalue, Optionally fieldstr can be "
                        "specified for a string defiition or smooth for "
                        "smoothing");

ProcessIsoContour::ProcessIsoContour(FieldSharedPtr f) :
    ProcessEquiSpacedOutput(f)
{

    m_config["fieldstr"]           = ConfigOption(false, "NotSet",
                                        "string of isocontour to be extracted");
    m_config["fieldname"]          = ConfigOption(false, "isocon",
                                        "name for isocontour if fieldstr "
                                        "specified, default is isocon");

    m_config["fieldid"]            = ConfigOption(false, "NotSet",
                                        "field id to extract");
    m_config["fieldvalue"]         = ConfigOption(false, "NotSet",
                                        "field value to extract");

    m_config["globalcondense"]     = ConfigOption(true, "NotSet",
                                        "Globally condense contour to unique "
                                        "values");

    m_config["smooth"]             = ConfigOption(true, "NotSet",
                                        "Smooth isocontour (implies global "
                                        "condense)");
    m_config["smoothiter"]         = ConfigOption(false, "100",
                                        "Number of smoothing cycle, default = "
                                        "100");
    m_config["smoothposdiffusion"] = ConfigOption(false,"0.5",
                                        "Postive diffusion coefficient "
                                        "(0 < lambda < 1), default = 0.5");
    m_config["smoothnegdiffusion"] = ConfigOption(false,"0.495",
                                        "Negative diffusion coefficient "
                                        "(0 < mu < 1), default = 0.495");

    m_config["removesmallcontour"] = ConfigOption(false,"0",
                                        "Remove contours with less than specified number of triangles."
                                                    "Only valid with GlobalCondense or Smooth options.");


}

ProcessIsoContour::~ProcessIsoContour(void)
{
}

void ProcessIsoContour::Process(po::variables_map &vm)
{
    vector<IsoSharedPtr> iso;

    // extract all fields to equi-spaced
    SetupEquiSpacedField();

    int     fieldid;
    NekDouble value;

    if(m_config["fieldstr"].m_beenSet) //generate field of interest
    {
        fieldid = m_f->m_fieldPts->GetNFields();

        Array<OneD, NekDouble> pts(m_f->m_fieldPts->GetNpoints());

        // evaluate new function
        LibUtilities::AnalyticExpressionEvaluator strEval;
        string varstr = "x y z";
        vector<Array<OneD, const NekDouble> > interpfields;

        for(int i = 0; i < m_f->m_fieldPts->GetDim(); ++i)
        {
            interpfields.push_back(m_f->m_fieldPts->GetPts(i));
        }
        for(int i = 0; i < m_f->m_fieldPts->GetNFields(); ++i)
        {
            varstr += " " + m_f->m_fieldPts->GetFieldName(i);
            interpfields.push_back(m_f->m_fieldPts->GetPts(i+3));
        }

        int ExprId  = -1;
        std::string str = m_config["fieldstr"].as<string>();
        ExprId = strEval.DefineFunction(varstr.c_str(), str);

        strEval.Evaluate(ExprId, interpfields, pts);

        // set up field name if provided otherwise called "isocon" from default
        string fieldName = m_config["fieldname"].as<string>();

        m_f->m_fieldPts->AddField(pts, fieldName);
    }
    else
    {
        ASSERTL0(m_config["fieldid"].as<string>() != "NotSet", "fieldid must be specified");
        fieldid = m_config["fieldid"].as<int>();
    }

    ASSERTL0(m_config["fieldvalue"].as<string>() != "NotSet", "fieldvalue must be specified");
    value   = m_config["fieldvalue"].as<NekDouble>();

    iso = ExtractContour(fieldid,value);
    bool smoothing      = m_config["smooth"].m_beenSet;
    bool globalcondense = m_config["globalcondense"].m_beenSet;
    if(smoothing||globalcondense)
    {
        vector<IsoSharedPtr> glob_iso;
        int nfields = m_f->m_fieldPts->GetNFields() + m_f->m_fieldPts->GetDim();
        IsoSharedPtr g_iso = MemoryManager<Iso>::AllocateSharedPtr(nfields-3);
        int mincontour = 0; 

        g_iso->globalcondense(iso);

        if(smoothing)
        {
            int  niter = m_config["smoothiter"].as<int>();
            NekDouble lambda = m_config["smoothposdiffusion"].as<NekDouble>();
            NekDouble mu     = m_config["smoothnegdiffusion"].as<NekDouble>();
            g_iso->smooth(niter,lambda,-mu);
        }


        if((mincontour = m_config["removesmallcontour"].as<int>()))
        {
            g_iso->separate_regions(glob_iso,mincontour);
        }
        else
        {
            glob_iso.push_back(g_iso);
        }
        ResetFieldPts(glob_iso);
    }
    else
    {
        ResetFieldPts(iso);
    }
}


/*********************/
/* Function TwoPairs */
/*********************/

void TwoPairs (Array<OneD, NekDouble> &cx,
               Array<OneD, NekDouble> &cy,
               Array<OneD, NekDouble> &cz,
               int &pr)
/*  The logic of the following SWITCH statement may only be
    understood with a "peusdo truth-table" */
{
    if (((cx[0]-cx[1])==0.0)&&
        ((cy[0]-cy[1])==0.0)&&
        ((cz[0]-cz[1])==0.0))
    {
        if (((cx[2]-cx[3])==0.0)&&
            ((cy[2]-cy[3])==0.0)&&
            ((cz[2]-cz[3])==0.0))
        {
            pr=4;
        }
        else
        {
            pr=3;
        }
    }
    else
    {
        pr=1;
    }
}


/*************************/
/* Function ThreeSimilar */
/*************************/
void ThreeSimilar (const int i, const int j, const int k,
                   int &pr, int &ps)
/*  The logic of the following SWITCH statement may be
    understood with a "peusdo truth-table" */
{
    switch (i + j + k)
    {
        case (3):
            pr = 3;
            ps = 4;
            break;
        case (4):
            pr = 2;
            ps = 4;
            break;
        case (5):
            if (j == 1)
            {
                pr = 2;
                ps = 3;
            }
            else
            {
                pr = 1;
                ps = 4;
            }
            break;
        case (6):
            if (i == 0)
            {
                pr = 1;
                ps = 3;
            }
            else
            {
                pr = 0;
                ps = 4;
            }
            break;
        case (7):
            if (i == 0)
            {
                pr = 1;
                ps = 2;
            }
            else
            {
                pr = 0;
                ps = 3;
            }
            break;
        case (8):
            pr = 0;
            ps = 2;
            break;
        case (9):
            pr = 0;
            ps = 1;
            break;
        default:
            ASSERTL0(false,"Error in 5-point triangulation in ThreeSimilar");
            break;
    }
}

vector<IsoSharedPtr> ProcessIsoContour::ExtractContour(
        const int fieldid,
        const NekDouble val)
{
    vector<IsoSharedPtr> returnval;

    int coordim = m_f->m_fieldPts->GetDim();
    int nfields = m_f->m_fieldPts->GetNFields() + coordim;

    ASSERTL0(m_f->m_fieldPts->GetPtsType() == LibUtilities::ePtsTetBlock,
             "This methods is currently only set up for 3D fields");
    ASSERTL1(coordim + fieldid < nfields,
             "field id is larger than number contained in FieldPts");
    Array<OneD, Array<OneD, NekDouble> > fields;
    m_f->m_fieldPts->GetPts(fields);

    Array<OneD, NekDouble> c = fields[coordim + fieldid];

    int i, j, k, ii, jj, kk, r, s, n, counter, boolean;
    Array<OneD, Array<OneD, NekDouble> > intfields(nfields);
    intfields[0] = Array<OneD, NekDouble>(5*nfields);
    for(i = 1; i < nfields; ++i)
    {
        intfields[i] = intfields[i-1] + 5;
    }
    Array<OneD, NekDouble> cx = intfields[0];
    Array<OneD, NekDouble> cy = intfields[1];
    Array<OneD, NekDouble> cz = intfields[2];

    vector<Array<OneD, int> > ptsConn;
    m_f->m_fieldPts->GetConnectivity(ptsConn);
    for(int zone = 0; zone < ptsConn.size(); ++zone)
    {
        IsoSharedPtr iso;

        iso = MemoryManager<Iso>::AllocateSharedPtr(nfields-3);

        int nelmt = ptsConn[zone].num_elements()
            /(coordim+1);

        Array<OneD, int> conn = ptsConn[zone];

        for (n = 0, i = 0; i < nelmt; ++i)
        {
            // check to see if val is between vertex values
            if(!(((c[conn[i*4]]  >val)&&(c[conn[i*4+1]]>val)&&
                  (c[conn[i*4+2]]>val)&&(c[conn[i*4+3]]>val))||
                 ((c[conn[i*4  ]]<val)&&(c[conn[i*4+1]]<val)&&
                  (c[conn[i*4+2]]<val)&&(c[conn[i*4+3]]<val))))
            {

                // loop over all edges and interpolate if
                // contour is between vertex values
                for (counter = 0, j=0; j<=2; j++)
                {
                    for (k=j+1; k<=3; k++)
                    {
                        if (((c[conn[i*4+j]]>=val)&&
                             (val>=c[conn[i*4+k]]))||
                            ((c[conn[i*4+j]]<=val)&&
                             (val<=c[conn[i*4+k]])))
                        {
                            // linear interpolation of fields
                            // (and coords).
                            NekDouble cj = c[conn[i*4+j]];
                            NekDouble ck = c[conn[i*4+k]];
                            NekDouble factor =  (val-cj)/(ck-cj);

                            if(fabs(cj-ck) > 1e-12)
                            {
                                // interpolate coordinates and fields
                                for(int f = 0; f < nfields; ++f)
                                {
                                    if(counter == 5)
                                    {
                                        ASSERTL0(false,"Counter is 5");
                                    }
                                    intfields[f][counter] =
                                        fields[f][conn[4*i+j]] +
                                        factor*(fields[f][conn[4*i+k]] -
                                                fields[f][conn[4*i+j]]);
                                }
                                ++counter;
                            }
                        }
                    }
                }

                switch(counter)
                {
                case 3:
                    n+=1;
                    iso->resize_fields(3*n);

                    for(j = 0; j < 3; ++j)
                    {
                        iso->set_fields(3*(n-1)+j,intfields,j);
                    }
                    break;
                case 4:
                    n+=2;
                    iso->resize_fields(3*n);

                    for(j = 0; j < 3; ++j)
                    {
                        iso->set_fields(3*(n-2)+j,intfields,j);
                        iso->set_fields(3*(n-1)+j,intfields,j+1);
                    }
                    break;
                case 5:
                    n+=1;
                    iso->resize_fields(3*n);

                    boolean=0;
                    for (ii=0;ii<=2;ii++)
                    {
                        for (jj=ii+1;jj<=3;jj++)
                        {
                            for (kk=jj+1;kk<=4;kk++)
                            {
                                if((((cx[ii]-cx[jj])==0.0)&&
                                    ((cy[ii]-cy[jj])==0.0)&&
                                    ((cz[ii]-cz[jj])==0.0))&&
                                   (((cx[ii]-cx[kk])==0.0)&&
                                    ((cy[ii]-cy[kk])==0.0)&&
                                    ((cz[ii]-cz[kk])==0.0)))
                                {
                                    boolean+=1;
                                    ThreeSimilar (ii,jj,kk,r,s);

                                    iso->set_fields(3*(n-1)  ,intfields,ii);
                                    iso->set_fields(3*(n-1)+1,intfields,r);
                                    iso->set_fields(3*(n-1)+2,intfields,s);
                                }
                                else
                                {
                                    boolean+=0;
                                }
                            }
                        }
                    }

                    if (boolean==0)
                    {
                        TwoPairs (cx,cy,cz,r);

                        iso->set_fields(3*(n-1)  ,intfields,0);
                        iso->set_fields(3*(n-1)+1,intfields,2);
                        iso->set_fields(3*(n-1)+2,intfields,r);
                    }
                    break;
                }
            }
        }
        iso->set_ntris(n);

        // condense the information in this elemental extraction.
        iso->condense();
        returnval.push_back(iso);
    }

    return returnval;
}

// reset m_fieldPts with values from iso;
void ProcessIsoContour::ResetFieldPts(vector<IsoSharedPtr> &iso)
{
    int nfields = m_f->m_fieldPts->GetNFields() + m_f->m_fieldPts->GetDim();

    // set output to triangle block.
    m_f->m_fieldPts->SetPtsType(LibUtilities::ePtsTriBlock);

    Array<OneD, Array<OneD, NekDouble> > newfields(nfields);

    // count up number of points
    int npts = 0;
    for(int i =0; i < iso.size(); ++i)
    {
        npts += iso[i]->get_nvert();
    }

    // set up coordinate in new field
    newfields[0] = Array<OneD, NekDouble>(npts);
    newfields[1] = Array<OneD, NekDouble>(npts);
    newfields[2] = Array<OneD, NekDouble>(npts);

    int cnt = 0;
    for(int i =0; i < iso.size(); ++i)
    {
        for(int j = 0; j < iso[i]->get_nvert(); ++j,++cnt)
        {
            newfields[0][cnt] = iso[i]->get_x(j);
            newfields[1][cnt] = iso[i]->get_y(j);
            newfields[2][cnt] = iso[i]->get_z(j);
        }
    }


    // set up fields
    for(int f = 0; f < nfields-3; ++f)
    {
        newfields[f+3] = Array<OneD, NekDouble>(npts);

        cnt = 0;
        for(int i =0; i < iso.size(); ++i)
        {
            for(int j = 0; j < iso[i]->get_nvert(); ++j,++cnt)
            {
                newfields[f+3][cnt] = iso[i]->get_fields(f,j);
            }
        }
    }

    m_f->m_fieldPts->SetPts(newfields);

    // set up connectivity data.
    vector<Array<OneD, int> > ptsConn;
    m_f->m_fieldPts->GetConnectivity(ptsConn);
    cnt = 0;
    ptsConn.clear();
    for(int i =0; i < iso.size(); ++i)
    {
        int ntris = iso[i]->get_ntris();
        Array<OneD, int> conn(ntris*3);

        for(int j = 0; j < 3*ntris; ++j)
        {
            conn[j] = cnt + iso[i]->get_vid(j);
        }
        ptsConn.push_back(conn);
        cnt += iso[i]->get_nvert();
    }
    m_f->m_fieldPts->SetConnectivity(ptsConn);
}

void Iso::condense(void)
{
    register int i,j,cnt;
    IsoVertex v;
    vector<IsoVertex> vert;
    vector<IsoVertex>::iterator pt;

    if(!m_ntris) return;

    if(m_condensed) return;
    m_condensed = true;

    vert.reserve(m_ntris/6);

    m_vid = Array<OneD, int>(3*m_ntris);

    // fill first 3 points and initialise fields
    v.m_fields.resize(m_fields.size());
    for(cnt =0, i = 0; i < 3; ++i)
    {
        v.m_x = m_x[i];
        v.m_y = m_y[i];
        v.m_z = m_z[i];
        for(int f = 0; f < m_fields.size(); ++f)
        {
            v.m_fields[f] = m_fields[f][i];
        }
        v.m_id = cnt;
        vert.push_back(v);
        m_vid[i] = v.m_id;
        ++cnt;
    }

    for(i = 1; i < m_ntris; ++i)
    {
        for(j = 0; j < 3; ++j)
        {
            v.m_x = m_x[3*i+j];
            v.m_y = m_y[3*i+j];
            v.m_z = m_z[3*i+j];

            pt = find(vert.begin(),vert.end(),v);
            if(pt != vert.end())
            {
                m_vid[3*i+j] = pt[0].m_id;
            }
            else
            {
                v.m_id = cnt;

                for(int f = 0; f < m_fields.size(); ++f)
                {
                    v.m_fields[f] = m_fields[f][3*i+j];
                }

                vert.push_back(v);

                m_vid[3*i+j] = v.m_id;
                ++cnt;
            }
        }
    }

    // remove elements with multiple vertices
    for(i = 0,cnt=0; i < m_ntris;)
    {
        if((m_vid[3*i]  ==m_vid[3*i+1])||
           (m_vid[3*i]  ==m_vid[3*i+2])||
           (m_vid[3*i+1]==m_vid[3*i+2]))
        {
            cnt++;
            for(j = 3*i; j < 3*(m_ntris-1); ++j)
            {
                m_vid[j] = m_vid[j+3];
            }
            m_ntris--;
        }
        else
        {
            ++i;
        }
    }

    m_nvert  = vert.size();

    m_x.resize(m_nvert);
    m_y.resize(m_nvert);
    m_z.resize(m_nvert);

    for(int f = 0; f < m_fields.size(); ++f)
    {
        m_fields[f].resize(m_nvert);
    }

    for(i = 0; i < m_nvert; ++i)
    {
        m_x[i] = vert[i].m_x;
        m_y[i] = vert[i].m_y;
        m_z[i] = vert[i].m_z;
        for(int f = 0; f < m_fields.size(); ++f)
        {
            m_fields[f][i] = vert[i].m_fields[f];
        }
    }
}


NekDouble SQ_PNT_TOL=1e-16;

// define == if point is within 1e-4
bool operator == (const IsoVertex& x, const IsoVertex& y)
{
    return ((x.m_x-y.m_x)*(x.m_x-y.m_x) + (x.m_y-y.m_y)*(x.m_y-y.m_y) +
            (x.m_z-y.m_z)*(x.m_z-y.m_z) < SQ_PNT_TOL)? true:false;
}

// define != if point is outside 1e-4
bool operator != (const IsoVertex& x, const IsoVertex& y)
{
    return ((x.m_x-y.m_x)*(x.m_x-y.m_x) + (x.m_y-y.m_y)*(x.m_y-y.m_y) +
            (x.m_z-y.m_z)*(x.m_z-y.m_z) < SQ_PNT_TOL)? 0:1;
}


bool same(NekDouble x1, NekDouble y1, NekDouble z1,
          NekDouble x2, NekDouble y2, NekDouble z2)
{
    if((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) < SQ_PNT_TOL)
    {
        return true;
    }

    return false;
}

void Iso::globalcondense(vector<IsoSharedPtr> &iso)
{
    int    i,j,n;
    int    nvert,nelmt;
    int    niso=iso.size();
    int    id1,id2;
    Array<OneD, Array<OneD, int> > vidmap;

    if(m_condensed) return;
    m_condensed = true;

    vidmap = Array<OneD, Array<OneD, int> > (niso);

    m_ntris = 0;
    for(i = 0; i < niso; ++i)
    {
        if(iso[i]->m_ntris)
        {
            m_ntris += iso[i]->m_ntris;
        }
    }

    m_vid = Array<OneD, int>(3*m_ntris);

    m_nvert = 0;
    for(i = 0; i < niso; ++i)
    {
        if(iso[i]->m_ntris)
        {
            m_nvert += iso[i]->m_nvert;
        }
    }

    vector< vector<int> > isocon;
    isocon.resize(niso);

    // identify which iso are connected by at least one point;
    // find min x,y,z and max x,y,z and see if overlap to select
    // which zones should be connected
    {
        vector<Array<OneD, NekDouble> > sph(niso);
        Array<OneD, NekDouble> rng(6);
        for(i = 0; i < niso; ++i)
        {
            sph[i] = Array<OneD, NekDouble>(4);

            // find max and min of isocontour
            rng[0] = rng[3] = iso[i]->m_x[0];
            rng[1] = rng[4] = iso[i]->m_x[1];
            rng[2] = rng[5] = iso[i]->m_x[2];

            for(id1 = 1; id1 < iso[i]->m_nvert;++id1)
            {
                rng[0] = min(rng[0],iso[i]->m_x[i]);
                rng[1] = min(rng[1],iso[i]->m_y[i]);
                rng[2] = min(rng[2],iso[i]->m_z[i]);

                rng[3] = max(rng[3],iso[i]->m_x[i]);
                rng[4] = max(rng[4],iso[i]->m_y[i]);
                rng[5] = max(rng[5],iso[i]->m_z[i]);
            }

            // centroid
            sph[i][0] = (rng[3]+rng[0])/2.0;
            sph[i][1] = (rng[4]+rng[1])/2.0;
            sph[i][2] = (rng[5]+rng[2])/2.0;

            // radius;
            sph[i][3] = sqrt((rng[3]-rng[0])*(rng[3]-rng[0]) +
                             (rng[4]-rng[1])*(rng[4]-rng[1]) +
                             (rng[5]-rng[2])*(rng[5]-rng[2]));
        }

        for(i = 0; i < niso; ++i)
        {
            for(j = i+1; j < niso; ++j)
            {
                NekDouble diff=sqrt((sph[i][0]-sph[j][0])*(sph[i][0]-sph[j][0])+
                          (sph[i][1]-sph[j][1])*(sph[i][1]-sph[j][1])+
                          (sph[i][2]-sph[j][2])*(sph[i][2]-sph[j][2]));

                // if centroid is closer than added radii
                if(diff < sph[i][3] + sph[j][3])
                {
                    isocon[i].push_back(j);
                }
            }
        }
    }


    for(i = 0; i < niso; ++i)
    {
        vidmap[i] = Array<OneD, int>(iso[i]->m_nvert,-1);
    }
    nvert = 0;
    // identify which vertices are connected to tolerance
    cout << "GlobalCondense: Matching Vertices [" << endl << flush;
    int cnt = 0;
    // count up amount of checking to be done
    NekDouble totpts = 0; 
    for(i = 0; i < niso; ++i)
    {
        for(n = 0; n < isocon[i].size(); ++n)
        {
            int con = isocon[i][n];
            totpts += iso[i]->m_nvert*iso[con]->m_nvert;
        }
    }
    int totchk  = totpts/50;
    int cnt_out = 0; 
    for(i = 0; i < niso; ++i)
    {
        for(n = 0; n < isocon[i].size(); ++n)
        {
            int con = isocon[i][n];
            for(id1 = 0; id1 < iso[i]->m_nvert; ++id1)
            {
                for(id2 = 0; id2 < iso[con]->m_nvert; ++id2, ++cnt)
                {
                    if(cnt%totchk == 0)
                    {
                        cout <<cnt_out << "%" << '\r' << flush;
                        cnt_out +=2; 
                    }

                    if((vidmap[con][id2] == -1)||(vidmap[i][id1] == -1))
                    {
                        if(same(iso[i]->m_x[id1],  iso[i]->m_y[id1],
                                iso[i]->m_z[id1],  iso[con]->m_x[id2],
                                iso[con]->m_y[id2],iso[con]->m_z[id2]))
                        {
                            if((vidmap[i][id1] == -1) &&
                               (vidmap[con][id2] != -1))
                            {
                                vidmap[i][id1] = vidmap[con][id2];
                            }
                            else if((vidmap[con][id2] == -1) &&
                                    (vidmap[i][id1] != -1))
                            {
                                vidmap[con][id2] = vidmap[i][id1];
                            }
                            else if((vidmap[con][id2] == -1) &&
                                    (vidmap[i][id1] == -1))
                            {
                                vidmap[i][id1] = vidmap[con][id2] = nvert++;
                            }
                        }
                    }
                }
            }
        }

        for(id1 = 0; id1 < iso[i]->m_nvert;++id1)
        {
            if(vidmap[i][id1] == -1)
            {
                vidmap[i][id1] = nvert++;
            }
        }
    }
    cout <<endl << "]"<<endl;
    m_nvert = nvert;

    nelmt = 0;
    // reset m_vid;
    for(n = 0; n < niso; ++n)
    {
        for(i = 0; i < iso[n]->m_ntris; ++i,nelmt++)
        {
            for(j=0; j < 3;++j)
            {
                m_vid[3*nelmt+j] = vidmap[n][iso[n]->m_vid[3*i+j]];
            }
        }
    }

    m_ntris = nelmt;

    m_x.resize(m_nvert);
    m_y.resize(m_nvert);
    m_z.resize(m_nvert);

    m_fields.resize(iso[0]->m_fields.size());
    for(i = 0; i < iso[0]->m_fields.size(); ++i)
    {
        m_fields[i].resize(m_nvert);
    }


    // reset coordinate and fields.
    for(n = 0; n < niso; ++n)
    {
        for(i = 0; i < iso[n]->m_nvert; ++i)
        {
            m_x[vidmap[n][i]] = iso[n]->m_x[i];
            m_y[vidmap[n][i]] = iso[n]->m_y[i];
            m_z[vidmap[n][i]] = iso[n]->m_z[i];

            for(j = 0; j < m_fields.size(); ++j)
            {
                m_fields[j][vidmap[n][i]] = iso[n]->m_fields[j][i];
            }
        }
    }
}

void Iso::smooth(int n_iter, NekDouble lambda, NekDouble mu)
{
    int   iter,i,j;
    NekDouble del_v[3];
    NekDouble w;
    Array<OneD, NekDouble>  xtemp, ytemp, ztemp;
    vector< vector<int> > adj,vertcon;
    vector<int>::iterator iad;
    vector<int>::iterator ipt;

    // determine elements around each vertex
    vertcon.resize(m_nvert);
    for(i = 0; i < m_ntris; ++i)
    {
        for(j = 0; j < 3; ++j)
        {
            vertcon[m_vid[3*i+j]].push_back(i);
        }
    }

    // determine vertices around each vertex
    adj.resize(m_nvert);

    for(i =0; i < m_nvert; ++i)
    {
        for(ipt = vertcon[i].begin(); ipt != vertcon[i].end(); ++ipt)
        {
            for(j = 0; j < 3; ++j)
            {
                // make sure not adding own vertex
                if(m_vid[3*(*ipt)+j] != i)
                {
                    // check to see if vertex has already been added
                    for(iad = adj[i].begin(); iad != adj[i].end();++iad)
                    {
                        if(*iad == (m_vid[3*(*ipt)+j])) break;
                    }

                    if(iad == adj[i].end())
                    {
                        adj[i].push_back(m_vid[3*(*ipt)+j]);
                    }
                }
            }
        }
    }

    xtemp = Array<OneD, NekDouble>(m_nvert);
    ytemp = Array<OneD, NekDouble>(m_nvert);
    ztemp = Array<OneD, NekDouble>(m_nvert);

    // smooth each point
    for (iter=0;iter<n_iter;iter++)
    {
        // compute first weighted average
        for(i=0;i< m_nvert;++i)
        {
            w = 1.0/(NekDouble)(adj[i].size());

            del_v[0] = del_v[1] = del_v[2] = 0.0;

            for(iad = adj[i].begin(); iad != adj[i].end(); ++iad)
            {
                del_v[0] =  del_v[0] + (m_x[*iad]-m_x[i])*w;
                del_v[1] =  del_v[1] + (m_y[*iad]-m_y[i])*w;
                del_v[2] =  del_v[2] + (m_z[*iad]-m_z[i])*w;
            }

            xtemp[i] = m_x[i] + del_v[0] * lambda;
            ytemp[i] = m_y[i] + del_v[1] * lambda;
            ztemp[i] = m_z[i] + del_v[2] * lambda;
        }

        // compute second weighted average
        for(i=0;i< m_nvert;++i)
        {

            w = 1.0/(NekDouble)(adj[i].size());
            del_v[0] = del_v[1] = del_v[2] = 0;

            for(iad = adj[i].begin(); iad != adj[i].end(); ++iad)
            {
                del_v[0] =  del_v[0] + (m_x[*iad]-m_x[i])*w;
                del_v[1] =  del_v[1] + (m_y[*iad]-m_y[i])*w;
                del_v[2] =  del_v[2] + (m_z[*iad]-m_z[i])*w;
            }

            m_x[i] = xtemp[i] + del_v[0] * mu;
            m_y[i] = ytemp[i] + del_v[1] * mu;
            m_z[i] = ztemp[i] + del_v[2] * mu;
        }
    }
}


    void Iso::separate_regions(vector<IsoSharedPtr> &sep_iso, int minsize)
    {
        int i,j,k,id;
        Array<OneD, vector<int> >vertcon(m_nvert);
        vector<int>::iterator ipt;
        list<int> tocheck;  
        list<int>::iterator cid;

        Array<OneD, bool> viddone(m_nvert,false);
        
        // make list of connecting tris around each vertex
        for(i = 0; i < m_ntris; ++i)
        {
            for(j = 0; j < 3; ++j)
            {
                vertcon[m_vid[3*i+j]].push_back(i);
            }
        }
        
        Array<OneD, int> vidregion(m_nvert,-1);
  
        int nregions = -1;

        cout << "Identifying separate regions [." << flush ;
        // check all points are assigned;
        for(k = 0; k < m_nvert; ++k)
        {
            if(vidregion[k] == -1)
            {
                vidregion[k] = ++nregions;
      
                // find all elmts around this.. vertex  that need to be checked
                for(ipt = vertcon[k].begin(); ipt != vertcon[k].end(); ++ipt)
                {
                    for(i = 0; i < 3; ++i)
                    {
                        if(vidregion[id = m_vid[3*(ipt[0])+i]] == -1)
                        {
                            tocheck.push_back(id);
                            vidregion[id] = nregions;
                        }
                    }
                }
                viddone[k] = 1;
                
                // check all other neighbouring vertices
                while(tocheck.size())
                {
                    cid = tocheck.begin();
                    while(cid != tocheck.end())
                    {
                        if(!viddone[*cid])
                        {
                            for(ipt = vertcon[*cid].begin(); ipt != vertcon[*cid].end(); ++ipt)
                            {
                                for(i = 0; i < 3; ++i)
                                {
                                    if(vidregion[id = m_vid[3*(ipt[0])+i]] == -1)
                                    {
                                        tocheck.push_back(id);
                                        vidregion[id] = nregions;
                                    }
                                }
                            }
                            viddone[*cid] = 1;
                            ++cid;
                            tocheck.pop_front();
                        }
                    }
                }
            }
        }
        nregions++;

        cout << "]" << endl <<  flush ;
        
        Array<OneD, int> nvert(nregions,0);
        Array<OneD, int> nelmt(nregions,0);

        // count nverts
        for(i = 0; i < m_nvert; ++i)
        {
            nvert[vidregion[i]] +=1;
        }

        // count nelmts
        for(i = 0; i < m_ntris; ++i)
        {
            nelmt[vidregion[m_vid[3*i]]] +=1;
        }

        Array<OneD, int> vidmap(m_nvert);
        // generate new list of isocontour
        for(int n = 0; n < nregions; ++n)
        {
            if(nelmt[n] > minsize)
            {
                int nfields = m_fields.size();
                IsoSharedPtr iso = MemoryManager<Iso>::AllocateSharedPtr(nfields);
                
                iso->m_ntris = nelmt[n];
                iso->m_vid = Array<OneD, int>(3*nelmt[n]);
                
                iso->m_nvert = nvert[n];
                iso->m_x.resize(nvert[n]);
                iso->m_y.resize(nvert[n]);
                iso->m_z.resize(nvert[n]);
                
                iso->m_fields.resize(nfields);
                for(i = 0; i < nfields; ++i)
                {
                    iso->m_fields[i].resize(nvert[n]);
                }


                int cnt = 0; 
                // generate vid map; 
                Vmath::Zero(m_nvert,vidmap,1);
                for(i = 0; i < m_nvert; ++i)
                {
                    if(vidregion[i] == n)
                    {
                        vidmap[i] = cnt++;
                    }
                }

                cnt = 0; 
                for(i = 0; i < m_ntris; ++i)
                {
                    if(vidregion[m_vid[3*i]] == n)
                    {
                        for(j = 0; j < 3; ++j)
                        {
                            iso->m_vid[3*cnt+j] = vidmap[m_vid[3*i+j]];
                            
                            iso->m_x[vidmap[m_vid[3*i+j]]] = m_x[m_vid[3*i+j]];
                            iso->m_y[vidmap[m_vid[3*i+j]]] = m_y[m_vid[3*i+j]];
                            iso->m_z[vidmap[m_vid[3*i+j]]] = m_z[m_vid[3*i+j]];

                            for(k = 0; k < nfields; ++k)
                            {
                                iso->m_fields[k][vidmap[m_vid[3*i+j]]] = m_fields[k][m_vid[3*i+j]];
                            }
                        }
                        cnt++;
                    }
                }

                WARNINGL0(cnt == nelmt[n],"Number of elements do not match");
                sep_iso.push_back(iso);
            }
        }
    }
}
}


