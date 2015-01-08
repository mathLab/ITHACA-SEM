////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessIsoContourUtils.cpp
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
//  Description: Set up fields as interpolation to equispaced output
//
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <iostream>
using namespace std;

#include "ProcessIsoContourUtils.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessIsoContourUtils::className =
    GetModuleFactory().RegisterCreatorFunction(
                        ModuleKey(eProcessModule, "isocontourutil"),
                        ProcessIsoContourUtils::create,
                        "process isocontour with different utilities");


ProcessIsoContourUtils::ProcessIsoContourUtils(FieldSharedPtr f)
{
    m_f = f;
    m_config["globalcondense"]     = ConfigOption(true, "NotSet",
                                        "Globally condense contour to unique "
                                        "values");

    m_config["smooth"]             = ConfigOption(true, "NotSet",
                                        "Smooth isocontour (implies global "
                                        "condense)");
    m_config["smoothiter"]         = ConfigOption(false, "100",
                                        "Number of smoothing cycle, "
                                        "default = 100");
    m_config["smoothposdiffusion"] = ConfigOption(false, "0.5",
                                        "Postive diffusion coefficient "
                                        "(0 < lambda < 1), default = 0.5");
    m_config["smoothnegdiffusion"] = ConfigOption(false, "0.495",
                                        "Negative diffusion coefficient "
                                        "(0 < mu < 1), default = 0.495");
}

ProcessIsoContourUtils::~ProcessIsoContourUtils(void)
{
}

void ProcessIsoContourUtils::Process(po::variables_map &vm)
{
    vector<IsoSharedPtr> iso;

    ASSERTL0(m_f->m_fieldPts != NullFieldPts,
             "Need to specify an input .dat file with isocontourutil module");

    LoadFieldPtsToIso(iso);

    bool smoothing      = m_config["smooth"].m_beenSet;
    bool globalcondense = m_config["globalcondense"].m_beenSet;
    if(smoothing||globalcondense)
    {
        vector<IsoSharedPtr> glob_iso;
        int nfields = m_f->m_fieldPts->m_pts.num_elements();
        IsoSharedPtr g_iso = MemoryManager<Iso>::AllocateSharedPtr(nfields-3);

        g_iso->globalcondense(iso);

        if(smoothing)
        {
            int  niter = m_config["smoothiter"].as<int>();
            NekDouble lambda = m_config["smoothposdiffusion"].as<NekDouble>();
            NekDouble mu     = m_config["smoothnegdiffusion"].as<NekDouble>();
            g_iso->smooth(niter,lambda,-mu);
        }

        glob_iso.push_back(g_iso);

        ResetFieldPts(glob_iso);
    }
    else
    {
        WARNINGL0(false,"No actions set");
    }
}

void ProcessIsoContourUtils::LoadFieldPtsToIso(vector<IsoSharedPtr> &iso)
{
    int nfields = m_f->m_fieldPts->m_pts.num_elements();

    // set output to triangle block.
    int nvertoffset = 0;
    for(int n = 0; n < m_f->m_fieldPts->m_ptsConn.size(); ++n)
    {
        IsoSharedPtr newiso = MemoryManager<Iso>::AllocateSharedPtr(nfields-3);

        int nconn = m_f->m_fieldPts->m_ptsConn[n].num_elements();
        newiso->set_ntris(nconn/3);

        int nvert = Vmath::Vmax(nconn,m_f->m_fieldPts->m_ptsConn[n],1)+1;
        nvert -= nvertoffset;

        newiso->resize_fields(nvert);

        newiso->set_nvert(nvert);

        for(int i = 0; i < nvert; ++i)
        {
            newiso->set_fields(i,m_f->m_fieldPts->m_pts,nvertoffset+i);
        }

        newiso->resize_vid(nconn);

        for(int i = 0; i < nconn; ++i)
        {
            newiso->set_vid(i,m_f->m_fieldPts->m_ptsConn[n][i]-nvertoffset);
        }

        nvertoffset += nvert;

        iso.push_back(newiso);
    }
}

}
}


