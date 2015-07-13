////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessMultiShear.cpp
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
//  Description: Computes tawss, osi, transwss, afi, cfi fields.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
#include <sstream>
using namespace std;

#include "ProcessMultiShear.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessMultiShear::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "shear"),
        ProcessMultiShear::create, "Computes shear stress metrics.");

ProcessMultiShear::ProcessMultiShear(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["N"] = ConfigOption(false,"1","Number of chk or fld files");
    m_config["fromfld"] = ConfigOption(false, "NotSet",
                                       "First fld file. First underscore flags position of id in name.");

    ASSERTL0(m_config["fromfld"].as<string>().compare("NotSet") != 0,
             "Need to specify fromfld=file.fld ");

    m_f->m_fldToBnd = false;
}

ProcessMultiShear::~ProcessMultiShear()
{
}

void ProcessMultiShear::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessMultiShear: Calculating shear stress metrics..." << endl;
    }

    int nstart, i, j, nfields;
    bool wssg = false;
    NekDouble nfld = m_config["N"].as<NekDouble>();
    string fromfld, basename, endname, nstartStr;
    stringstream filename;
    vector<string> infiles(nfld);
    vector< boost::shared_ptr<Field> > m_fromField(nfld);

    // Set up list of input fld files.
    fromfld = m_config["fromfld"].as<string>();
    basename = fromfld.substr(0, fromfld.find_first_of("_")+1);
    filename << fromfld.substr(fromfld.find_first_of("_")+1, fromfld.size());
    filename >> nstart;
    filename.str("");
    filename << nstart;
    filename >> nstartStr;
    filename.str("");
    endname = fromfld.substr(fromfld.find(nstartStr)+nstartStr.size(), fromfld.size());

    for (i=0; i<nfld; ++i)
    {
        stringstream filename;
        filename << basename << i+nstart << endname;
        filename >> infiles[i];
        cout << infiles[i]<<endl;
    }

    for ( i = 0; i<nfld; ++i)
    {
        m_fromField[i] = boost::shared_ptr<Field>(new Field());
        m_fromField[i]->m_session = m_f->m_session;
        m_fromField[i]->m_graph = m_f->m_graph;
        m_fromField[i]->m_fld = MemoryManager<LibUtilities::FieldIO>
            ::AllocateSharedPtr(m_fromField[0]->m_session->GetComm());
    }

    //Import all fld files.
    for (i = 0; i < nfld; ++i)
    {
        if(m_f->m_exp.size())
        {
            // Set up ElementGIDs in case of parallel processing
            Array<OneD,int> ElementGIDs(m_f->m_exp[0]->GetExpSize());
            for (j = 0; j < m_f->m_exp[0]->GetExpSize(); ++j)
            {
                ElementGIDs[j] = m_f->m_exp[0]->GetExp(j)->GetGeom()->GetGlobalID();
            }
            m_fromField[i]->m_fld->Import(infiles[i],m_fromField[i]->m_fielddef,
                                          m_fromField[i]->m_data,
                                          LibUtilities::NullFieldMetaDataMap,
                                          ElementGIDs);
        }
        else
        {
            m_fromField[i]->m_fld->Import(infiles[i],m_fromField[i]->m_fielddef,
                                          m_fromField[i]->m_data,
                                          LibUtilities::NullFieldMetaDataMap);
        }

        nfields    = m_fromField[i]->m_fielddef[0]->m_fields.size();
        int NumHomogeneousDir = m_fromField[i]->m_fielddef[0]->m_numHomogeneousDir;

        if (nfields == 5 || nfields == 7)
        {
            wssg = true;
        }

        // Set up Expansion information to use mode order from field
        m_fromField[i]->m_graph->SetExpansions(m_fromField[i]->m_fielddef);

        //Set up expansions, and extract data.
        m_fromField[i]->m_exp.resize(nfields);
        m_fromField[i]->m_exp[0] = m_fromField[i]->SetUpFirstExpList(NumHomogeneousDir,true);

        for(j = 1; j < nfields; ++j)
        {
            m_fromField[i]->m_exp[j] = m_f->AppendExpList(NumHomogeneousDir);
        }

        for (j = 0; j < nfields; ++j)
        {
            for (int k = 0; k < m_fromField[i]->m_data.size(); ++k)
            {
                m_fromField[i]->m_exp[j]->ExtractDataToCoeffs(
                    m_fromField[i]->m_fielddef[k],
                    m_fromField[i]->m_data[k],
                    m_fromField[i]->m_fielddef[k]->m_fields[j],
                    m_fromField[i]->m_exp[j]->UpdateCoeffs());
            }
            m_fromField[i]->m_exp[j]->BwdTrans(m_fromField[i]->m_exp[j]->GetCoeffs(),
                                               m_fromField[i]->m_exp[j]->UpdatePhys());
        }
    }


    int spacedim   = m_f->m_graph->GetSpaceDimension();
    if ((m_fromField[0]->m_fielddef[0]->m_numHomogeneousDir) == 1 ||
        (m_fromField[0]->m_fielddef[0]->m_numHomogeneousDir) == 2)
    {
        spacedim = 3;
    }


    int nout = 5; // TAWSS, OSI, transWSS, AFI, CFI (optional WSSG)
    if (wssg) { nout = 6; }

    int npoints = m_fromField[0]->m_exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > normTemporalMeanVec(spacedim),normCrossDir(spacedim),  outfield(nout), dtemp(spacedim);
    Array<OneD, NekDouble> TemporalMeanMag(npoints,0.0), DotProduct(npoints,0.0), temp(npoints,0.0);


    for (i = 0; i < spacedim; ++i)
    {
        normTemporalMeanVec[i] = Array<OneD, NekDouble>(npoints);
        normCrossDir[i] = Array<OneD, NekDouble>(npoints);
        dtemp[i] = Array<OneD, NekDouble>(npoints);
        Vmath::Zero(npoints, normTemporalMeanVec[i],1);
        Vmath::Zero(npoints, normCrossDir[i],1);
    }
    
    for (i = 0; i < nout; ++i)
    {
        outfield[i] = Array<OneD, NekDouble>(npoints, 0.0);
    }

    // -----------------------------------------------------
    // Compute temporal average wall shear stress vector,
    // it's spatial average, and normalise it.
    for (i = 0; i < nfld; ++i)
    {
        for (j = 0; j < spacedim; ++j)
        {
            Vmath::Vadd(npoints,m_fromField[i]->m_exp[j]->GetPhys(), 1, normTemporalMeanVec[j], 1, normTemporalMeanVec[j], 1); 
        }
    }

    for (i = 0; i < spacedim; ++i)
    {
        Vmath::Smul(npoints, 1.0/nfld, normTemporalMeanVec[i], 1, normTemporalMeanVec[i], 1);
        Vmath::Vvtvp(npoints, normTemporalMeanVec[i], 1, normTemporalMeanVec[i], 1,
                     TemporalMeanMag, 1, TemporalMeanMag, 1);
    }

    Vmath::Vsqrt(npoints, TemporalMeanMag, 1, TemporalMeanMag, 1);

    for (i = 0; i < spacedim; ++i)
    {
        Vmath::Vdiv(npoints, normTemporalMeanVec[i], 1, TemporalMeanMag, 1, normTemporalMeanVec[i], 1);
    }
    //---------------------------------------------------

    if (wssg) //cross product with normals to obtain direction parallel to temporal mean vector. 
    {
        Vmath::Vmul(npoints,m_fromField[0]->m_exp[nfields-1]->GetPhys(),1,normTemporalMeanVec[1],1,normCrossDir[0],1);
        Vmath::Vvtvm(npoints,m_fromField[0]->m_exp[nfields-2]->GetPhys(),1,normTemporalMeanVec[2],1,normCrossDir[0],1,
                     normCrossDir[0],1);

        Vmath::Vmul(npoints,m_fromField[0]->m_exp[nfields-3]->GetPhys(),1,normTemporalMeanVec[2],1,normCrossDir[1],1);
        Vmath::Vvtvm(npoints,m_fromField[0]->m_exp[nfields-1]->GetPhys(),1,normTemporalMeanVec[0],1,normCrossDir[1],1,
                     normCrossDir[1],1);
        
        Vmath::Vmul(npoints,m_fromField[0]->m_exp[nfields-2]->GetPhys(),1,normTemporalMeanVec[0],1,normCrossDir[2],1);
        Vmath::Vvtvm(npoints,m_fromField[0]->m_exp[nfields-3]->GetPhys(),1,normTemporalMeanVec[1],1,normCrossDir[2],1,
                     normCrossDir[2],1);
    }


    // Compute tawss, trs,  osi, taafi, tacfi, WSSG.
    for (i = 0; i < nfld; ++i)
    {
        for (j = 0; j < spacedim; ++j)
        {
            Vmath::Vvtvp(npoints, m_fromField[i]->m_exp[j]->GetPhys(), 1, normTemporalMeanVec[j], 1,
                         DotProduct, 1, DotProduct, 1);
        }
        
        //TAWSS
        Vmath::Vadd(npoints, m_fromField[i]->m_exp[spacedim]->GetPhys(), 1, outfield[0], 1, outfield[0], 1);
        
        //transWSS
        Vmath::Vmul(npoints, DotProduct, 1, DotProduct, 1, temp,1);
        Vmath::Vvtvm(npoints, m_fromField[i]->m_exp[spacedim]->GetPhys(), 1, 
                     m_fromField[i]->m_exp[spacedim]->GetPhys(), 1, temp, 1, temp, 1);
        
        for (j = 0; j < npoints; ++j)
        {
            if(temp[j] > 0.0)
            {
                outfield[1][j] = outfield[1][j] + sqrt(temp[j]);
            }
        }
        
        //TAAFI
        Vmath::Vdiv(npoints, DotProduct, 1, m_fromField[i]->m_exp[spacedim]->GetPhys(), 1, temp, 1);
        Vmath::Vadd(npoints, temp, 1, outfield[3], 1, outfield[3], 1);
        
        //TACFI
        for (j = 0; j < npoints; ++j)
        {
            temp[j] = 1 - temp[j]*temp[j];
            if(temp[j] > 0.0)
            {
                outfield[4][j] = outfield[4][j] + sqrt(temp[j]);
            }
        }
        
        // WSSG
        if (wssg)
        {
            //parallel component:
            Vmath::Zero(npoints, temp,1);

            m_fromField[i]->m_exp[0]->PhysDeriv(DotProduct,dtemp[0],dtemp[1],dtemp[2]);
            for(j=0; j<spacedim; j++)
            {
                Vmath::Vvtvp(npoints,dtemp[j],1,normTemporalMeanVec[j],1,temp,1,temp,1); 
            }
            Vmath::Vmul(npoints,temp,1,temp,1,temp,1);


            //perpendicular component:
            Vmath::Zero(npoints, DotProduct,1);
           
            for(j = 0; j < spacedim; ++j)
            {
                Vmath::Vvtvp(npoints,m_fromField[i]->m_exp[j]->GetPhys(),1,normCrossDir[j],1,
                             DotProduct, 1, DotProduct, 1);
                
            }
            m_fromField[i]->m_exp[0]->PhysDeriv(DotProduct,dtemp[0],dtemp[1],dtemp[2]);
            Vmath::Zero(npoints, DotProduct,1);

            for(j=0; j<spacedim; j++)
            {
                Vmath::Vvtvp(npoints,dtemp[j],1,normCrossDir[j],1,DotProduct,1,DotProduct,1); 
            }
            Vmath::Vvtvp(npoints,DotProduct,1,DotProduct,1,temp,1, temp,1);

            //Final wssg computation
            Vmath::Vsqrt(npoints,temp,1,temp,1);
            Vmath::Vadd(npoints, temp, 1, outfield[5], 1, outfield[5], 1);
        }
        
        
        Vmath::Zero(npoints, DotProduct,1);
    }

    //Divide by nfld
    Vmath::Smul(npoints, 1.0/nfld, outfield[0], 1, outfield[0], 1);
    Vmath::Smul(npoints, 1.0/nfld, outfield[1], 1, outfield[1], 1);
    Vmath::Smul(npoints, 1.0/nfld, outfield[3], 1, outfield[3], 1);
    Vmath::Smul(npoints, 1.0/nfld, outfield[4], 1, outfield[4], 1);

    //OSI
    for (i = 0; i < npoints; ++i)
    {
        outfield[2][i] = 0.5 * (1 - TemporalMeanMag[i]/outfield[0][i]);
    }

    /* TAWSS = sum(wss)/nfld
     * transWSS = sum( sqrt( wss^2 - (wss . normTempMean)^2) )/nfld.
     * OSI = 0.5*(1-TemporalMeanMag/TAWSS)
     * TAAFI = sum(cos)/nfld
     * TACFI = sum(sin)/nfld = sum( sqrt(1-cos^2) )/nfld.
     */

    m_f->m_exp.resize(nout);
    m_f->m_fielddef = m_fromField[0]->m_fielddef;
    m_f->m_exp[0] = m_f->SetUpFirstExpList(m_f->m_fielddef[0]->m_numHomogeneousDir,true);

    for(i = 1; i < nout; ++i)
    {
        m_f->m_exp[i] = m_f->AppendExpList(m_f->m_fielddef[0]->m_numHomogeneousDir);
    }

    m_f->m_fielddef[0]->m_fields.resize(nout);
    m_f->m_fielddef[0]->m_fields[0] = "TAWSS";
    m_f->m_fielddef[0]->m_fields[1] = "transWSS";
    m_f->m_fielddef[0]->m_fields[2] = "OSI";
    m_f->m_fielddef[0]->m_fields[3] = "TAAFI";
    m_f->m_fielddef[0]->m_fields[4] = "TACFI";

    if (wssg)
    { 
        m_f->m_fielddef[0]->m_fields[5] = "|WSSG|"; 
        Vmath::Smul(npoints, 1.0/nfld, outfield[5], 1, outfield[5], 1);
    }


    for(i = 0; i < nout; ++i)
    {
        m_f->m_exp[i]->FwdTrans(outfield[i],
                                m_f->m_exp[i]->UpdateCoeffs());
        m_f->m_exp[i]->BwdTrans(m_f->m_exp[i]->GetCoeffs(),
                                m_f->m_exp[i]->UpdatePhys());
    }


    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = m_fromField[0]->m_exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    for( i = 0; i < nout; ++i)
    {
        for ( j = 0; j < FieldDef.size(); ++j)
        {
            FieldDef[j]->m_fields.push_back(m_f->m_fielddef[0]->m_fields[i]);
            m_f->m_exp[i]->AppendFieldData(FieldDef[j], FieldData[j]);
        }
    }

    m_f->m_fielddef = FieldDef;
    m_f->m_data     = FieldData;
}

}
}
