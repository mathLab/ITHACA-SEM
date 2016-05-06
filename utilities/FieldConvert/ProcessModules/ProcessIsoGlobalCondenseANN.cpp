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

#include "ANN/ANN.h"

namespace Nektar
{
namespace Utilities
{
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
    vector<int> global_to_unique_map;
    global_to_unique_map.resize(m_nvert);

    vector< pair<int,int> > global_to_iso_map;
    global_to_iso_map.resize(m_nvert);

    //Create kdtree
    int n_neighbs = 5;
    int neighbs_max = 100; 
    ANNpointArray dataPts = annAllocPts(m_nvert, 3);
    ANNidxArray   nnIdx = new ANNidx [neighbs_max];
    ANNdistArray  dists = new ANNdist[neighbs_max];
    
    //Fill vertex array into libAnn format
    id2 = 0;
    for(i = 0; i < niso; ++i)
    {
        for(id1 = 0; id1 < iso[i]->m_nvert; ++id1)
        {
            dataPts[id2][0] = iso[i]->m_x[id1];
            dataPts[id2][1] = iso[i]->m_y[id1];
            dataPts[id2][2] = iso[i]->m_z[id1];
            global_to_unique_map[id2]=id2;
            global_to_iso_map[id2] = make_pair(i,id1);
            id2++;
        }
    }
    
    //Build tree
    ANNkd_tree* kdTree;
    // build search structure
    kdTree = new ANNkd_tree( dataPts,	  // the data points
                             m_nvert,	  // number of points
                             3);	  // dimension of space
    
    //Find neipghbours
    ANNpoint queryPt = annAllocPt(3);
    int      unique_index = 0;
    bool     unique_index_found = false;
    int      prog; 
    for(i = 0; i < m_nvert; ++i)
    {
        if(verbose)
        {
            prog = LibUtilities::PrintProgressbar(i,m_nvert,"Nearest verts",prog);
        }

        n_neighbs  = 5; 
        queryPt[0] = dataPts[i][0];
        queryPt[1] = dataPts[i][1];
        queryPt[2] = dataPts[i][2];
        kdTree->annkSearch(queryPt, n_neighbs, nnIdx, dists, 0); //eps set to zero


        while((dists[n_neighbs-1]<SQ_PNT_TOL) && (n_neighbs*2 < neighbs_max))
        {
            n_neighbs *=2; 
            kdTree->annkSearch(queryPt, n_neighbs, nnIdx, dists, 0); //eps set to zero
        }

        WARNINGL0(n_neighbs*2 < neighbs_max,"Failed to find less than 100 neighbouring points");

        id1 = 0;
        unique_index_found = false;

        int nptsfound = 0; 
        for(id1 = 0; id1 < n_neighbs; ++id1)
        {
            if(dists[id1]<SQ_PNT_TOL)
            {
                id2 = nnIdx[id1];
                nptsfound ++;
                if(global_to_unique_map[id2] <unique_index) 
                {
                    // point has already been defined
                    continue; 
                }
                else
                {
                    global_to_unique_map[id2] = unique_index;
                    unique_index_found = true;
                }
            }
        }

        
        if(unique_index_found)
        {
            unique_index++;
        }
    }
    if(verbose)
    {
        cout << endl;
    }

    m_nvert = unique_index;

    nelmt = 0;
    // reset m_vid;
    int cnt = 0; 
    for(n = 0; n < niso; ++n)
    {
        for(i = 0; i < iso[n]->m_ntris; ++i,nelmt++)
        {
            for(j=0; j < 3;++j)
            {
                m_vid[3*nelmt+j] = global_to_unique_map[iso[n]->m_vid[3*i+j]+cnt]; 
            }
        }
        cnt += iso[n]->m_nvert; 
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

    for(n = 0; n < global_to_unique_map.size(); ++n)
    {
        m_x[global_to_unique_map[n]] = dataPts[n][0];
        m_y[global_to_unique_map[n]] = dataPts[n][1];
        m_z[global_to_unique_map[n]] = dataPts[n][2];
        
        int isoid = global_to_iso_map[n].first;
        int ptid  = global_to_iso_map[n].second;
        for(j = 0; j < m_fields.size(); ++j)
        {
            m_fields[j][global_to_unique_map[n]] = iso[isoid]->
                m_fields[j][ptid];
        }
    }
}    
}
}


