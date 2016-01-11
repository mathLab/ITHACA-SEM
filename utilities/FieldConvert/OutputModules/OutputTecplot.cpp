////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputTecplot.cpp
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
//  Description: Dat file format output.
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
#include <iomanip>
using namespace std;

#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/PtsField.h>

#include "OutputTecplot.h"

namespace Nektar
{
namespace Utilities
{

ModuleKey OutputTecplot::m_className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eOutputModule, "dat"), OutputTecplot::create,
        "Writes a Tecplot file.");

OutputTecplot::OutputTecplot(FieldSharedPtr f) : OutputModule(f)
{
    if(f->m_setUpEquiSpacedFields)
    {
        m_outputType = eFullBlockZoneEquiSpaced;
    }
    else
    {
        m_requireEquiSpaced = true;
        m_outputType = eFullBlockZone;
    }
}

OutputTecplot::~OutputTecplot()
{
}

void OutputTecplot::Process(po::variables_map &vm)
{
    LibUtilities::PtsFieldSharedPtr fPts = m_f->m_fieldPts;

    m_doError = (vm.count("error") == 1)?  true: false;

    if (m_f->m_verbose)
    {
        cout << "OutputTecplot: Writing file..." << endl;
    }
    // Do nothing if no expansion defined
    if (fPts == LibUtilities::NullPtsField && !m_f->m_exp.size())
    {
        return;
    }


    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();

    // Amend for parallel output if required
    if(m_f->m_comm->GetSize() != 1)
    {
        int    dot = filename.find_last_of('.');
        string ext = filename.substr(dot,filename.length()-dot);
        string procId = "_P" + boost::lexical_cast<std::string>(
                                 m_f->m_comm->GetRank());
        string start = filename.substr(0,dot);
        filename = start + procId + ext;
    }

    if(fPts != LibUtilities::NullPtsField)
    {
        int i   = 0;
        int j   = 0;
        int dim = fPts->GetDim();

        if(fPts->GetNpoints() == 0)
        {
            return;
        }

        // Write solution.
        ofstream outfile(filename.c_str());

        switch(dim)
        {
        case 1:
            outfile << "VARIABLES = x";
            break;
        case 2:
            outfile << "VARIABLES = x,y";
            break;
        case 3:
            outfile << "VARIABLES = x,y,z";
            break;
        }

        vector<Array<OneD, int> > ptsConn;
        fPts->GetConnectivity(ptsConn);

        for(i = 0; i < fPts->GetNFields(); ++i)
        {
            outfile << "," << fPts->GetFieldName(i);
        }
        outfile << endl;
        bool DumpAsFEPoint = true;
        switch(fPts->GetPtsType())
        {
            case LibUtilities::ePtsFile:
            case LibUtilities::ePtsLine:
            {
                outfile << " ZONE I="
                        << fPts->GetNpoints()
                        << " F=POINT" << endl;
                break;
            }
            case LibUtilities::ePtsPlane:
            {
                outfile << " ZONE I=" << fPts->GetPointsPerEdge(0)
                        <<      " J=" << fPts->GetPointsPerEdge(1)
                        << " F=POINT" << endl;
                break;
            }
            case LibUtilities::ePtsTriBlock:
            {
                int numBlocks = 0;
                for(i = 0; i < ptsConn.size(); ++i)
                {
                    numBlocks +=
                        ptsConn[i].num_elements()/3;
                }
                outfile << "Zone, N="
                        << fPts->GetNpoints()
                        << ", E=" << numBlocks
                        << ", F=FEBlock" << ", ET=TRIANGLE"
                        << std::endl;
                DumpAsFEPoint = false;
                break;
            }
            case LibUtilities::ePtsTetBlock:
            {
                int numBlocks = 0;
                for(i = 0; i < ptsConn.size(); ++i)
                {
                    numBlocks +=
                        ptsConn[i].num_elements()/4;
                }
                outfile << "Zone, N="
                        << fPts->GetNpoints()
                        << ", E=" << numBlocks
                        << ", F=FEBlock" << ", ET=TETRAHEDRON"
                        << std::endl;
                DumpAsFEPoint = false;
                break;
            }
            default:
                ASSERTL0(false, "ptsType not supported yet.");
        }

        if(DumpAsFEPoint) // dump in point format
        {
            for(i = 0; i < fPts->GetNpoints(); ++i)
            {
                for(j = 0; j < dim; ++j)
                {
                    outfile << std::setw(12)
                            << fPts->GetPointVal(j, i) << " ";
                }

                for(j = 0; j < fPts->GetNFields(); ++j)
                {
                    outfile << std::setw(12)
                            << m_f->m_data[j][i] << " ";
                }
                outfile << endl;
            }
        }
        else // dump in block format
        {
            for(j = 0; j < dim + fPts->GetNFields(); ++j)
            {
                for(i = 0; i < fPts->GetNpoints(); ++i)
                {
                    outfile <<  fPts->GetPointVal(j, i) << " ";
                    if((!(i % 1000))&&i)
                    {
                        outfile << std::endl;
                    }
                }
                outfile << endl;
            }

            // dump connectivity data if it exists
            for(i = 0; i < ptsConn.size();++i)
            {
                for(j = 0; j < ptsConn[i].num_elements(); ++j)
                {
                    outfile << ptsConn[i][j] +1 << " ";
                    if( ( !(j % 10 * dim) ) && j )
                    {
                        outfile << std::endl;
                    }
                }
            }
        }
    }
    else
    {

        // Write solution.
        ofstream outfile(filename.c_str());
        std::string var;
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> fDef =
                                                        m_f->m_fielddef;
        if (fDef.size())
        {
            var = fDef[0]->m_fields[0];

            for (int j = 1; j < fDef[0]->m_fields.size(); ++j)
            {
                var = var + ", " + fDef[0]->m_fields[j];
            }
        }

        WriteTecplotHeader(outfile,var);
        WriteTecplotZone(outfile);
        if(var.length()) // see if any variables are defined
        {
            for(int j = 0; j < m_f->m_exp.size(); ++j)
            {
                WriteTecplotField(j,outfile);
            }
        }

        WriteTecplotConnectivity(outfile);
    }

    cout << "Written file: " << filename << endl;
}

/**
 * Write Tecplot Files Header
 * @param   outfile Output file name.
 * @param   var                 variables names
 */
void OutputTecplot::WriteTecplotHeader(std::ofstream &outfile,
                                       std::string var)
{
    int coordim  = m_f->m_exp[0]->GetExp(0)->GetCoordim();
    MultiRegions::ExpansionType HomoExpType = m_f->m_exp[0]->GetExpType();

    if(HomoExpType == MultiRegions::e3DH1D)
    {
        if(m_f->m_session->DefinesSolverInfo("ModeType")&&
           boost::iequals(m_f->m_session->GetSolverInfo("ModeType"),"HalfMode"))
        { // turn off for half mode case
        }
        else
        {
            coordim +=1;
        }
    }
    else if (HomoExpType == MultiRegions::e3DH2D)
    {
        coordim += 2;
    }

    outfile << "Variables = x";

    if(coordim == 2)
    {
        outfile << ", y";
    }
    else if (coordim == 3)
    {
        outfile << ", y, z";
    }

    if(var.length())
    {
        outfile << ", "<< var << std::endl << std::endl;
    }
    else
    {
        outfile << std::endl << std::endl;
    }
}


/**
 * Write Tecplot Files Zone
 * @param   outfile    Output file name.
 * @param   expansion  Expansion that is considered
 */
void OutputTecplot::WriteTecplotZone(std::ofstream &outfile)
{
    switch(m_outputType)
    {
        case eFullBlockZone: //write as full block zone
        {
            int i,j;
            int coordim   = m_f->m_exp[0]->GetCoordim(0);
            int totpoints = m_f->m_exp[0]->GetTotPoints();
            MultiRegions::ExpansionType HomoExpType = m_f->m_exp[0]->GetExpType();

            Array<OneD,NekDouble> coords[3];

            coords[0] = Array<OneD,NekDouble>(totpoints);
            coords[1] = Array<OneD,NekDouble>(totpoints);
            coords[2] = Array<OneD,NekDouble>(totpoints);

            m_f->m_exp[0]->GetCoords(coords[0],coords[1],coords[2]);

            if (m_doError)
            {
                NekDouble l2err;
                std::string coordval[] = {"x","y","z"};
                int rank = m_f->m_comm->GetRank();

                for(int i = 0; i < coordim; ++i)
                {
                    l2err = m_f->m_exp[0]->L2(coords[i]);
                    if(rank == 0)
                    {
                        cout << "L 2 error (variable "
                             << coordval[i]  << ") : " << l2err  << endl;
                    }
                }
            }

            int numBlocks = GetNumTecplotBlocks();
            int nBases = m_f->m_exp[0]->GetExp(0)->GetNumBases();

            if (HomoExpType == MultiRegions::e3DH1D)
            {
                int nPlanes = m_f->m_exp[0]->GetZIDs().num_elements();
                if(nPlanes == 1) // halfMode case
                {
                    // do nothing
                }
                else
                {
                    nBases  += 1;
                    coordim += 1;
                    NekDouble tmp = numBlocks * (nPlanes-1);
                    numBlocks = (int)tmp;
                }
            }
            else if (HomoExpType == MultiRegions::e3DH2D)
            {
                nBases  += 2;
                coordim += 1;
            }

            outfile << "Zone, N=" << totpoints << ", E="<<
                numBlocks << ", F=FEBlock" ;

            switch(nBases)
            {
            case 1:
                outfile << ", ET=LINESEG" << std::endl;
                break;
            case 2:
                outfile << ", ET=QUADRILATERAL" << std::endl;
                break;
            case 3:
                outfile << ", ET=BRICK" << std::endl;
                break;
            }

            // write out coordinates in block format
            for(j = 0; j < coordim; ++j)
            {
                for(i = 0; i < totpoints; ++i)
                {
                    outfile << coords[j][i] << " ";
                    if((!(i % 1000))&&i)
                    {
                        outfile << std::endl;
                    }
                }
                outfile << std::endl;
            }
            break;
        }
        case eSeperateZones:
        {
            for(int i = 0; i < m_f->m_exp[0]->GetExpSize(); ++i)
            {
                m_f->m_exp[0]->WriteTecplotZone(outfile,i);
            }
            break;
        }
        case eFullBlockZoneEquiSpaced:
            ASSERTL0(false,
                     "Should not have this option in this method");
            break;
    }
}

int OutputTecplot::GetNumTecplotBlocks(void)
{
    int returnval = 0;

    if(m_f->m_exp[0]->GetExp(0)->GetNumBases() == 1)
    {
        for(int i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
        {
            returnval += (m_f->m_exp[0]->GetExp(i)->GetNumPoints(0)-1);
        }
    }
    else  if(m_f->m_exp[0]->GetExp(0)->GetNumBases() == 2)
    {
        for(int i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
        {
            returnval += (m_f->m_exp[0]->GetExp(i)->GetNumPoints(0)-1)*
                (m_f->m_exp[0]->GetExp(i)->GetNumPoints(1)-1);
        }
    }
    else
    {
        for(int i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
        {
            returnval += (m_f->m_exp[0]->GetExp(i)->GetNumPoints(0)-1)*
                (m_f->m_exp[0]->GetExp(i)->GetNumPoints(1)-1)*
                (m_f->m_exp[0]->GetExp(i)->GetNumPoints(2)-1);
        }
    }

    return returnval;
}


/**
 * Write Tecplot Files Field
 * @param   outfile    Output file name.
 * @param   expansion  Expansion that is considered
 */
void OutputTecplot::WriteTecplotField(const int field,
                                      std::ofstream &outfile)
{

    if(m_outputType == eFullBlockZone) //write as full block zone
    {
        int totpoints = m_f->m_exp[0]->GetTotPoints();

        if(m_f->m_exp[field]->GetPhysState() == false)
        {
            m_f->m_exp[field]->BwdTrans(m_f->m_exp[field]->GetCoeffs(),
                                        m_f->m_exp[field]->UpdatePhys());
        }

        if (m_doError)
        {
            NekDouble l2err = m_f->m_exp[0]->L2(m_f->m_exp[field]->UpdatePhys());

            if(m_f->m_comm->GetRank() == 0)
            {
                cout << "L 2 error (variable "
                     << m_f->m_fielddef[0]->m_fields[field]  << ") : "
                     << l2err  << endl;
            }
        }
        else
        {
            for(int i = 0; i < totpoints; ++i)
            {
                outfile << m_f->m_exp[field]->GetPhys()[i] << " ";
                if((!(i % 1000))&&i)
                {
                    outfile << std::endl;
                }
            }
        }
        outfile << std::endl;
    }
    else
    {
        for(int e = 0; e < m_f->m_exp[field]->GetExpSize(); ++e)
        {
            m_f->m_exp[field]->WriteTecplotField(outfile,e);
        }
    }
}

void  OutputTecplot::WriteTecplotConnectivity(std::ofstream &outfile)
{
    int i,j,k,l;
    int nbase = m_f->m_exp[0]->GetExp(0)->GetNumBases();
    int cnt = 0;

    for(i = 0; i < m_f->m_exp[0]->GetNumElmts(); ++i)
    {
        if(nbase == 1)
        {
            int np0 = m_f->m_exp[0]->GetExp(i)->GetNumPoints(0);

            for(k = 1; k < np0; ++k)
            {
                outfile << cnt + k +1 << " ";
                outfile << cnt + k    << endl;
            }

            cnt += np0;
        }
        else if(nbase == 2)
        {
            int np0 = m_f->m_exp[0]->GetExp(i)->GetNumPoints(0);
            int np1 = m_f->m_exp[0]->GetExp(i)->GetNumPoints(1);
            int totPoints = m_f->m_exp[0]->GetTotPoints();
            int nPlanes = 1;

            if (m_f->m_exp[0]->GetExpType() == MultiRegions::e3DH1D)
            {
                nPlanes = m_f->m_exp[0]->GetZIDs().num_elements();

                if(nPlanes > 1) // default to 2D case for HalfMode when nPlanes = 1
                {
                    totPoints = m_f->m_exp[0]->GetPlane(0)->GetTotPoints();


                    for(int n = 1; n < nPlanes; ++n)
                    {
                        for(j = 1; j < np1; ++j)
                        {
                            for(k = 1; k < np0; ++k)
                            {
                                outfile << cnt + (n-1)*totPoints + (j-1)*np0 + k
                                        << " ";
                                outfile << cnt + (n-1)*totPoints + (j-1)*np0 + k + 1
                                        << " ";
                                outfile << cnt + (n-1)*totPoints + j*np0 + k + 1
                                        << " ";
                                outfile << cnt + (n-1)*totPoints + j*np0 + k
                                    << " ";
                                
                                outfile << cnt + n*totPoints + (j-1)*np0 + k
                                        << " ";
                                outfile << cnt + n*totPoints + (j-1)*np0 + k + 1
                                        << " ";
                                outfile << cnt + n*totPoints + j*np0 + k + 1
                                        << " ";
                                outfile << cnt + n*totPoints + j*np0 + k    << endl;
                            }
                        }
                    }
                    cnt += np0*np1;
                }
            }

            if(nPlanes == 1)
            {
                for(j = 1; j < np1; ++j)
                {
                    for(k = 1; k < np0; ++k)
                    {
                        outfile << cnt + (j-1)*np0 + k  << " ";
                        outfile << cnt + (j-1)*np0 + k +1 << " ";
                        outfile << cnt + j*np0 + k +1 << " ";
                        outfile << cnt + j*np0 + k    << endl;
                    }
                }
                cnt += np0*np1;
            }
        }
        else if(nbase == 3)
        {
            int np0 = m_f->m_exp[0]->GetExp(i)->GetNumPoints(0);
            int np1 = m_f->m_exp[0]->GetExp(i)->GetNumPoints(1);
            int np2 = m_f->m_exp[0]->GetExp(i)->GetNumPoints(2);

            for(j = 1; j < np2; ++j)
            {
                for(k = 1; k < np1; ++k)
                {
                    for(l = 1; l < np0; ++l)
                    {
                        outfile << cnt + (j-1)*np0*np1 + (k-1)*np0 + l  << " ";
                        outfile << cnt + (j-1)*np0*np1 + (k-1)*np0 + l +1 << " ";
                        outfile << cnt + (j-1)*np0*np1 +  k*np0 + l +1 << " ";
                        outfile << cnt + (j-1)*np0*np1 +  k*np0 + l  << " ";

                        outfile << cnt + j*np0*np1 + (k-1)*np0 + l  << " ";
                        outfile << cnt + j*np0*np1 + (k-1)*np0 + l +1 << " ";
                        outfile << cnt + j*np0*np1 +  k*np0 + l +1 << " ";
                        outfile << cnt + j*np0*np1 +  k*np0 + l  << endl;
                    }
                }
            }
            cnt += np0*np1*np2;
        }
        else
        {
            ASSERTL0(false,"Not set up for this dimension");
        }

    }
}

}
}




