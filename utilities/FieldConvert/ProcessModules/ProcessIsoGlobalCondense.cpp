///////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessIsoGlobalCondense.cpp
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
//  Description: Globally condense isocontour - default versoin 
//
///////////////////////////////////////////////////////////////////////////////
#include <string>
#include <iostream>
using namespace std;

#include "ProcessIsoContour.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <LibUtilities/BasicUtils/Progressbar.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace Nektar
{
namespace Utilities
{

bool same(NekDouble x1, NekDouble y1, NekDouble z1,
          NekDouble x2, NekDouble y2, NekDouble z2)
{
    if((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) < SQ_PNT_TOL)
    {
        return true;
    }

    return false;
}
    
    
void Iso::globalcondense(vector<IsoSharedPtr> &iso, bool verbose)
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
            for(j = i; j < niso; ++j)
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
    int cnt = 0;
    // count up amount of checking to be done
    NekDouble totiso = 0; 
    for(i = 0; i < niso; ++i)
    {
        totiso += isocon[i].size();
    }


    if(verbose)
    {
        cout << "Progress Bar totiso: " << totiso << endl;
    }

    int progcnt = -1; 
    for(i = 0; i < niso; ++i)
    {
        for(n = 0; n < isocon[i].size(); ++n, ++cnt)
        {
            
            if(verbose && totiso >= 40)
            {
                progcnt = LibUtilities::PrintProgressbar(cnt,totiso,"Condensing verts",progcnt);
            }

            int con = isocon[i][n];
            for(id1 = 0; id1 < iso[i]->m_nvert; ++id1)
            {

                if(verbose && totiso < 40)
                {
                     if(cnt % (int)(totiso/200) == 0)
                     {
                          progcnt =  LibUtilities::PrintProgressbar(id1,iso[i]->m_nvert,"isocon",progcnt);
                     }
                }

                int start  = 0; 
                if(con == i)
                {
                    start = id1+1;
                }
                for(id2 = start; id2 < iso[con]->m_nvert; ++id2)
                {
                    
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
    cout << endl;
}
    
}
}


