////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessScalGrad.cpp
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
//  Description: Computes scalar gradient field.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessScalGrad.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessScalGrad::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "scalargrad"),
        ProcessScalGrad::create, "Computes scalar gradient field.");

ProcessScalGrad::ProcessScalGrad(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["bnd"] = ConfigOption(false,"All","Boundary to be extracted");
    f->m_writeBndFld = true;
    f->m_declareExpansionAsContField = true;
    m_f->m_fldToBnd = false;
}

ProcessScalGrad::~ProcessScalGrad()
{
}

void ProcessScalGrad::Process(po::variables_map &vm)
{
    int i, j, k;
    if (m_f->m_verbose)
    {
        cout << "ProcessScalGrad: Calculating scalar gradient..." << endl;
    }

    // Set up Field options to output boundary fld
    string bvalues =  m_config["bnd"].as<string>();

    if(bvalues.compare("All") == 0)
    {
        Array<OneD, const MultiRegions::ExpListSharedPtr>
            BndExp = m_f->m_exp[0]->GetBndCondExpansions();

        for(i = 0; i < BndExp.num_elements(); ++i)
        {
            m_f->m_bndRegionsToWrite.push_back(i);
        }
    }
    else
    {
        ASSERTL0(ParseUtils::GenerateOrderedVector(bvalues.c_str(),
                                                   m_f->m_bndRegionsToWrite),"Failed to interpret range string");
    }

    int spacedim  = m_f->m_graph->GetSpaceDimension();
    if ((m_f->m_fielddef[0]->m_numHomogeneousDir) == 1 ||
        (m_f->m_fielddef[0]->m_numHomogeneousDir) == 2)
    {
        spacedim = 3;
    }

    int nfields = m_f->m_fielddef[0]->m_fields.size();
    //ASSERTL0(nfields == 1,"Implicit assumption that input is in ADR format of (u)");

    if (spacedim == 1)
    {
        ASSERTL0(false, "Error: scalar gradient for a 1D problem cannot "
                 "be computed");
    }


    int ngrad     = spacedim;
    int n, cnt, elmtid, nq, offset, boundary, nfq;
    int npoints = m_f->m_exp[0]->GetNpoints();
    string var;
    Array<OneD, NekDouble> scalar;
    Array<OneD, Array<OneD, NekDouble> > grad(ngrad), fgrad(ngrad), outfield(nfields);

    StdRegions::StdExpansionSharedPtr elmt;
    StdRegions::StdExpansion2DSharedPtr bc;
    Array<OneD, int> BoundarytoElmtID, BoundarytoTraceID;
    Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> > BndExp(nfields);

    m_f->m_exp[0]->GetBoundaryToElmtMap(BoundarytoElmtID, BoundarytoTraceID);

    for (i = 0; i < nfields; i++)
    {
        var = m_f->m_fielddef[0]->m_fields[i];
        stringstream filename;
        filename << var << "_scalar_gradient";
        filename >> var;
        m_f->m_fielddef[0]->m_fields[i] = var;
        
        BndExp[i] = m_f->m_exp[i]->GetBndCondExpansions();
        outfield[i] = Array<OneD, NekDouble>(npoints);
    }
    
    // loop over the types of boundary conditions
    for(cnt = n = 0; n < BndExp[0].num_elements(); ++n)
    {
        bool doneBnd = false;
        // identify if boundary has been defined
        for(int b = 0; b < m_f->m_bndRegionsToWrite.size(); ++b)
        {
            if(n == m_f->m_bndRegionsToWrite[b])
            {
                doneBnd = true;
                for(i = 0; i < BndExp[0][n]->GetExpSize(); ++i, cnt++)
                {
                    // find element and face of this expansion.
                    elmtid = BoundarytoElmtID[cnt];
                    elmt   = m_f->m_exp[0]->GetExp(elmtid);
                    nq     = elmt->GetTotPoints();
                    offset = m_f->m_exp[0]->GetPhys_Offset(elmtid);

                    // Initialise local arrays for the velocity gradients, and stress components
                    // size of total number of quadrature points for each element (hence local).
                    for(j = 0; j < ngrad; ++j)
                    {
                        grad[j] = Array<OneD, NekDouble>(nq);
                    }

                    if(spacedim == 2)
                    {
                        //Not implemented in 2D.
                    }
                    else
                    {
                        for(j = 0; j < nfields; j++)
                        {
                            outfield[j] = BndExp[j][n]->UpdateCoeffs() + BndExp[j][n]->GetCoeff_Offset(i);
                        }

                        // Get face 2D expansion from element expansion
                        bc  =  boost::dynamic_pointer_cast<StdRegions::StdExpansion2D> (BndExp[0][n]->GetExp(i));
                        nfq =  bc->GetTotPoints();

                        //identify boundary of element looking at.
                        boundary = BoundarytoTraceID[cnt];

                        //Get face normals
                        const SpatialDomains::GeomFactorsSharedPtr m_metricinfo = bc->GetMetricInfo();

                        const Array<OneD, const Array<OneD, NekDouble> > normals
                            = elmt->GetFaceNormal(boundary);

                        // initialise arrays
                        for(j = 0; j < ngrad; ++j)
                        {
                            fgrad[j] = Array<OneD, NekDouble>(nfq);
                        }
                        Array<OneD, NekDouble>  gradnorm(nfq);

                        for(k = 0; k < nfields; k++)
                        {
                            Vmath::Zero(nfq, gradnorm, 1);

                            scalar = m_f->m_exp[k]->GetPhys() + offset;
                            elmt->PhysDeriv(scalar, grad[0],grad[1],grad[2]);
                            
                            for(j = 0; j < ngrad; ++j)
                            {
                                elmt->GetFacePhysVals(boundary,bc,grad[j],fgrad[j]);
                            }
                            
                            //surface curved
                            if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                            {
                                for (j=0; j<ngrad; j++)
                                {
                                    Vmath::Vvtvp(nfq, normals[j], 1, fgrad[j], 1, gradnorm, 1, gradnorm, 1);
                                }
                            }
                            else
                            {
                                for (j=0; j<ngrad; j++)
                                {
                                    Vmath::Svtvp(nfq, normals[j][0], fgrad[j], 1, gradnorm, 1, gradnorm, 1);
                                }
                            }
                            bc->FwdTrans(gradnorm, outfield[k]);
                        }

                    }
                }
            }
        }
        if(doneBnd == false)
        {
            cnt += BndExp[0][n]->GetExpSize();
        }
    }

    for(j = 0; j < nfields; ++j)
    {
        for(int b = 0; b < m_f->m_bndRegionsToWrite.size(); ++b)
        {
            m_f->m_exp[j]->UpdateBndCondExpansion(m_f->m_bndRegionsToWrite[b]) = BndExp[j][m_f->m_bndRegionsToWrite[b]];
        }

    }
}

}
}
