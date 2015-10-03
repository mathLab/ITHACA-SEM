///////////////////////////////////////////////////////////////////////
//
//  File: OutputMTCProp.cpp
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
//  Description: mtc file in propriety format. 
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
#include <iomanip>
using namespace std;

#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/PtsField.h>

#include "OutputMTCProp.h"

namespace Nektar
{
namespace Utilities
{

ModuleKey OutputMTCProp::m_className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eOutputModule, "mtc"), OutputMTCProp::create,
        "Writes Juliens special format.");

OutputMTCProp::OutputMTCProp(FieldSharedPtr f) : OutputModule(f)
{

    m_config["ascii"] = ConfigOption(true,"0","write output in ascii format");
}

OutputMTCProp::~OutputMTCProp()
{
}

void OutputMTCProp::Process(po::variables_map &vm)
{
    LibUtilities::PtsFieldSharedPtr fPts = m_f->m_fieldPts;


    if (m_f->m_verbose)
    {
        cout << "OutputMTCProp: Writing file..." << endl;
    }

    // Do nothing if no expansion defined
    if (fPts == LibUtilities::NullPtsField)
    {
        cout << "No field points specified: Have you used -m interppoints?" << endl;
        return;
    }

    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();

    int nprocs = m_f->m_comm->GetSize();
    int rank  = m_f->m_comm->GetRank();
    // Amend for parallel output if required
    if(nprocs != 1)
    {
        int    dot = filename.find_last_of('.');
        string ext = filename.substr(dot,filename.length()-dot);
        string procId = "_P" + boost::lexical_cast<std::string>(rank);
        string start = filename.substr(0,dot);
        filename = start + procId + ext;
    }

    int i   = 0;
    int j   = 0;
    
    if(fPts->GetNpoints() == 0)
    {
        cout << "No points found" << endl;
        return;
    }
    // Write solution.
    ofstream outfile(filename.c_str());
        
    // points type
    LibUtilities::PtsType pType = fPts->GetPtsType();
    
    ASSERTL0(pType == LibUtilities::ePtsBox,"MTC Propriety format only suitable for box version of interppoints");

    if(m_config["ascii"].as<bool>())
    {
        //ascii output for debugging. 
        // write header information 
        if(rank == 0)
        {
            int nvar = fPts->GetNFields();
            int ni = fPts->GetPointsPerEdge(0);
            int nj = fPts->GetPointsPerEdge(1);
            int nk = fPts->GetPointsPerEdge(0);
            
            outfile << nvar << " "; 
            outfile << ni  << " "; 
            outfile << nj  << " "; 
            outfile << nk  << " "; 
            
            float x0,y0,z0,dx,dy,dz;
            
            vector<NekDouble> boxsize = fPts->GetBoxSize();
            
            ASSERTL0(boxsize.size() == 6,"Insufficient points specified in box size");
            
            x0 = boxsize[0];
            y0 = boxsize[2];
            z0 = boxsize[4];
            
            outfile << x0 << " ";
            outfile << y0 << " ";
            outfile << z0 << " ";

            dx = (boxsize[1] - boxsize[0])/(double)ni;
            dy = (boxsize[3] - boxsize[2])/(double)nj;
            dz = (boxsize[5] - boxsize[4])/(double)nk;
            
            outfile << dx << " ";
            outfile << dy << " ";
            outfile << dz << " ";
            
            // write out number of components of each variable
            int one = 1;
            for(int i = 0; i < nvar; ++i)
            {
                outfile << one << " "; 
            }
            
            for(i = 0; i < fPts->GetNFields(); ++i)
            {
                string var("0",80);
                var = fPts->GetFieldName(i);
                //var.resize(80);
                outfile << var << " "; 
            }
        }
        
        for(i = 0; i < fPts->GetNpoints(); ++i)
        {
            for(j = 0; j < fPts->GetNFields(); ++j)
            {
                float data = m_f->m_data[j][i];
                outfile << data << " "; 
            }
        }
    }
    else
    {
        // write header information 
        if(rank == 0)
        {
            int nvar = fPts->GetNFields();
            int ni = fPts->GetPointsPerEdge(0);
            int nj = fPts->GetPointsPerEdge(1);
            int nk = fPts->GetPointsPerEdge(0);
            
            
            outfile.write((char*)& nvar,sizeof(int));
            outfile.write((char*)& ni,  sizeof(int));
            outfile.write((char*)& nj,  sizeof(int));
            outfile.write((char*)& nk,  sizeof(int));
            
            float x0,y0,z0,dx,dy,dz;
            
            vector<NekDouble> boxsize = fPts->GetBoxSize();
            
            ASSERTL0(boxsize.size() == 6,"Insufficient points specified in box size");
            
            x0 = boxsize[0];
            y0 = boxsize[2];
            z0 = boxsize[4];
            
            outfile.write((char*)& x0,  sizeof(float));
            outfile.write((char*)& y0,  sizeof(float));
            outfile.write((char*)& z0,  sizeof(float));
            
            dx = (boxsize[1] - boxsize[0])/(double)ni;
            dy = (boxsize[3] - boxsize[2])/(double)nj;
            dz = (boxsize[5] - boxsize[4])/(double)nk;
            
            outfile.write((char*)& dx,  sizeof(float));
            outfile.write((char*)& dy,  sizeof(float));
            outfile.write((char*)& dz,  sizeof(float));
            
            // write out number of components of each variable
            int one = 1;
            for(int i = 0; i < nvar; ++i)
            {
                outfile.write((char*)& one,  sizeof(int));
            }
            
            for(i = 0; i < fPts->GetNFields(); ++i)
            {
                string var = fPts->GetFieldName(i);
                var.resize(80);
                outfile.write((char*)&var[0],80*sizeof(char));
            }
        }    
        
        for(i = 0; i < fPts->GetNpoints(); ++i)
        {
            for(j = 0; j < fPts->GetNFields(); ++j)
            {
                float data = m_f->m_data[j][i];
                outfile.write((char *) &data,sizeof(float));
            }
        }
    }

    cout << "Written file: " << filename << endl;
}

}
}




