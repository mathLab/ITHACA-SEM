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
    m_requireEquiSpaced = true;
}

OutputTecplot::~OutputTecplot()
{
}

void OutputTecplot::Process(po::variables_map &vm)
{
    m_doError = (vm.count("error") == 1)?  true: false;

    if (m_f->m_verbose)
    {
        cout << "OutputTecplot: Writing file..." << endl;
    }

    // Do nothing if no expansion defined
    if(m_f->m_fieldPts == NullFieldPts &&!m_f->m_exp.size())
    {
        return;
    }

    // Extract the output filename and extension
    string filename = m_config["outfile"].as<string>();

    if(m_f->m_fieldPts != NullFieldPts)
    {
        int dim = m_f->m_fieldPts->m_ptsDim;
        // Write solution.
        ofstream outfile(filename.c_str());

        switch(dim)
        {
        case 1:
            outfile << "VARIABLES = x";
            break;
        case 2:
            outfile << "VARIABLES = x y";
            break;
        case 3:
            outfile << "VARIABLES = x y z";
            break;
        }

        for(int i = 0; i < m_f->m_fieldPts->m_fields.size(); ++i)
        {
            outfile << " " << m_f->m_fieldPts->m_fields[i];
        }
        outfile << endl;

        switch(m_f->m_fieldPts->m_ptype)
        {
        case Utilities::ePtsFile:
        case Utilities::ePtsLine:
            outfile << " ZONE I="
                    << m_f->m_fieldPts->m_pts[0].num_elements()
                    << " F=POINT" << endl;
            break;
        case Utilities::ePtsPlane:
            outfile << " ZONE I=" << m_f->m_fieldPts->m_npts[0]
                    <<      " J=" << m_f->m_fieldPts->m_npts[1]
                    << " F=POINT" << endl;
            break;
        default:
            ASSERTL0(false, "Points type not supported yet.");
        }

        for(int i = 0; i < m_f->m_fieldPts->m_pts[0].num_elements(); ++i)
        {
            for(int j = 0; j < dim; ++j)
            {
                outfile << std::setw(12)
                        << m_f->m_fieldPts->m_pts[j][i] << " ";
            }

            for(int j = 0; j < m_f->m_fieldPts->m_fields.size(); ++j)
            {
                outfile << std::setw(12) << m_f->m_data[j][i] << " ";
            }
            outfile << endl;
        }
    }
    else
    {
        // Amend for parallel output if required
        if(m_f->m_session->GetComm()->GetSize() != 1)
        {
            int    dot = filename.find_last_of('.');
            string ext = filename.substr(dot,filename.length()-dot);
            string procId = "_P" + boost::lexical_cast<std::string>(
                                        m_f->m_session->GetComm()->GetRank());
            string start = filename.substr(0,dot);
            filename = start + procId + ext;
        }

        // Write solution.
        ofstream outfile(filename.c_str());
        std::string var;
        if(m_f->m_fielddef.size())
        {
            var = m_f->m_fielddef[0]->m_fields[0];

            for (int j = 1; j < m_f->m_fielddef[0]->m_fields.size(); ++j)
            {
                var = var + ", " + m_f->m_fielddef[0]->m_fields[j];
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
        coordim +=1;
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
void OutputTecplot::WriteTecplotZone(std::ofstream &outfile,
                                     int expansion)
{
    if(expansion == -1) //write as full block zone
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
            int rank = m_f->m_session->GetComm()->GetRank();

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
            nBases  += 1;
            coordim += 1;
            int nPlanes = m_f->m_exp[0]->GetZIDs().num_elements();
            NekDouble tmp = numBlocks * (nPlanes-1);
            numBlocks = (int)tmp;
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

    }
    else
    {
        m_f->m_exp[0]->WriteTecplotZone(outfile,expansion);
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
                                      std::ofstream &outfile,
                                      int expansion)
{

    if(expansion == -1) //write as full block zone
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

            if(m_f->m_session->GetComm()->GetRank() == 0)
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
        m_f->m_exp[field]->WriteTecplotField(outfile,expansion);
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
            else
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




