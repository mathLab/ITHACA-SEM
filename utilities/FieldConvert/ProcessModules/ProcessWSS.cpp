////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessWSS.cpp
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
//  Description: Computes wss field.
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include "ProcessWSS.h"

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.hpp>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
namespace Utilities
{

ModuleKey ProcessWSS::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "wss"),
        ProcessWSS::create, "Computes wall shear stress field.");

ProcessWSS::ProcessWSS(FieldSharedPtr f) : ProcessModule(f)
{
    m_config["bnd"] = ConfigOption(false,"All","Boundary to be extracted");
    m_config["addnormals"] = ConfigOption(true,"NotSet","Add normals to output");
    f->m_writeBndFld = true;
    f->m_declareExpansionAsContField = true;
    m_f->m_fldToBnd = false;
}

ProcessWSS::~ProcessWSS()
{
}

void ProcessWSS::Process(po::variables_map &vm)
{
    if (m_f->m_verbose)
    {
        cout << "ProcessWSS: Calculating wall shear stress..." << endl;
    }

    m_f->m_addNormals = m_config["addnormals"].m_beenSet;

    // Set up Field options to output boundary fld
    string bvalues =  m_config["bnd"].as<string>();

    if(bvalues.compare("All") == 0)
    {
        Array<OneD, const MultiRegions::ExpListSharedPtr>
            BndExp = m_f->m_exp[0]->GetBndCondExpansions();

        for(int i = 0; i < BndExp.num_elements(); ++i)
        {
            m_f->m_bndRegionsToWrite.push_back(i);
        }
    }
    else
    {
        ASSERTL0(ParseUtils::GenerateOrderedVector(bvalues.c_str(),
                                                   m_f->m_bndRegionsToWrite),"Failed to interpret range string");
    }


    NekDouble m_kinvis;
    m_kinvis = m_f->m_session->GetParameter("Kinvis");

    int i, j;
    int spacedim  = m_f->m_graph->GetSpaceDimension();
    if ((m_f->m_fielddef[0]->m_numHomogeneousDir) == 1 ||
        (m_f->m_fielddef[0]->m_numHomogeneousDir) == 2)
    {
        spacedim = 3;
    }

    int nfields = m_f->m_fielddef[0]->m_fields.size();
    ASSERTL0(nfields == spacedim +1,"Implicit assumption that input is in incompressible format of (u,v,p) or (u,v,w,p)");

    nfields = nfields - 1;
    if (spacedim == 1)
    {
        ASSERTL0(false, "Error: wss for a 1D problem cannot "
                 "be computed");
    }

    int newfields = (spacedim == 2)? 3:4;
    int nshear    = (spacedim == 2)? 3:4;
    int nstress   = (spacedim == 2)? 3:6;
    int ngrad     = nfields*nfields;

    int n, cnt, elmtid, nq, offset, boundary, nfq;
    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > velocity(nfields), grad(ngrad), fgrad(ngrad);
    Array<OneD, Array<OneD, NekDouble> > stress(nstress), fstress(nstress);
    Array<OneD, Array<OneD, NekDouble> > outfield(newfields), fshear(nshear);

    StdRegions::StdExpansionSharedPtr elmt;
    StdRegions::StdExpansion2DSharedPtr bc;
    Array<OneD, int> BoundarytoElmtID, BoundarytoTraceID;
    Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr> > BndExp(newfields);

    m_f->m_exp.resize(newfields);
    string var = "u";
    for(i = nfields+1; i < newfields; ++i)
    {
        m_f->m_exp[i] = m_f->AppendExpList(m_f->m_fielddef[0]->m_numHomogeneousDir, var);
    }

    m_f->m_fielddef[0]->m_fields.resize(newfields);
    if(spacedim == 2)
    {
        m_f->m_fielddef[0]->m_fields[0] = "Shear_x";
        m_f->m_fielddef[0]->m_fields[1] = "Shear_y";
        m_f->m_fielddef[0]->m_fields[2] = "Shear_mag";
    }
    else
    {
        m_f->m_fielddef[0]->m_fields[0] = "Shear_x";
        m_f->m_fielddef[0]->m_fields[1] = "Shear_y";
        m_f->m_fielddef[0]->m_fields[2] = "Shear_z";
        m_f->m_fielddef[0]->m_fields[3] = "Shear_mag";
    }


    for (i = 0; i < newfields; ++i)
    {
        outfield[i] = Array<OneD, NekDouble>(npoints);
    }
    for (i = 0; i < nfields; ++i)
    {
        velocity[i] = Array<OneD, NekDouble>(npoints);
    }

    m_f->m_exp[0]->GetBoundaryToElmtMap(BoundarytoElmtID, BoundarytoTraceID);
    //get boundary expansions for each field
    for(int j = 0; j < newfields; ++j)
    {
        BndExp[j]   = m_f->m_exp[j]->GetBndCondExpansions();
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
                for(int i = 0; i < BndExp[0][n]->GetExpSize(); ++i, cnt++)
                {
                    // find element and face of this expansion.
                    elmtid = BoundarytoElmtID[cnt];
                    elmt   = m_f->m_exp[0]->GetExp(elmtid);
                    nq     = elmt->GetTotPoints();
                    offset = m_f->m_exp[0]->GetPhys_Offset(elmtid);

                    // Initialise local arrays for the velocity gradients, and stress components
                    // size of total number of quadrature points for each element (hence local).
                    for(int j = 0; j < ngrad; ++j)
                    {
                        grad[j] = Array<OneD, NekDouble>(nq);
                    }

                    for(int j = 0; j < nstress; ++j)
                    {
                        stress[j] = Array<OneD, NekDouble>(nq);
                    }

                    if(nfields == 2)
                    {
                        ASSERTL0(false, "Error: not implemented in 2D.");
                    }
                    else
                    {
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
                        for(int j = 0; j < nstress; ++j)
                        {
                            fstress[j] = Array<OneD, NekDouble>(nfq);
                        }

                        for(int j = 0; j < nfields*nfields; ++j)
                        {
                            fgrad[j] = Array<OneD, NekDouble>(nfq);
                        }

                        for(int j = 0; j < nshear; ++j)
                        {
                            fshear[j] = Array<OneD, NekDouble>(nfq);
                        }


                        //Extract Velocities
                        for(int j = 0; j < nfields; ++j)
                        {
                            velocity[j] = m_f->m_exp[j]->GetPhys() + offset;
                        }

                        //Compute gradients (velocity correction scheme method)
                        elmt->PhysDeriv(velocity[0],grad[0],grad[1],grad[2]);
                        elmt->PhysDeriv(velocity[1],grad[3],grad[4],grad[5]);
                        elmt->PhysDeriv(velocity[2],grad[6],grad[7],grad[8]);

                         //Compute stress component terms
                        // t_xx = 2.mu.Ux
                        Vmath::Smul (nq,(2*m_kinvis),grad[0],1,stress[0],1);
                        // tyy = 2.mu.Vy
                        Vmath::Smul (nq,(2*m_kinvis),grad[4],1,stress[1],1);
                        // tzz = 2.mu.Wz
                        Vmath::Smul (nq,(2*m_kinvis),grad[8],1,stress[2],1);
                        // txy = mu.(Uy+Vx)
                        Vmath::Vadd (nq,grad[1],1,grad[3],1,stress[3],1);
                        Vmath::Smul (nq,m_kinvis,stress[3],1,stress[3],1);
                        // txz = mu.(Uz+Wx)
                        Vmath::Vadd (nq,grad[2],1,grad[6],1,stress[4],1);
                        Vmath::Smul (nq,m_kinvis,stress[4],1,stress[4],1);
                        // tyz = mu.(Vz+Wy)
                        Vmath::Vadd (nq,grad[5],1,grad[7],1,stress[5],1);
                        Vmath::Smul (nq,m_kinvis,stress[5],1,stress[5],1);


                        // Get face stress values.
                        for(j = 0; j < nstress; ++j)
                        {
                            elmt->GetFacePhysVals(boundary,bc,stress[j],fstress[j]);
                        }

                        //calcuate wss, and update velocity coeffs in the elemental boundary expansion
                        for (j = 0; j< newfields; j++)
                        {
                            outfield[j] = BndExp[j][n]->UpdateCoeffs() + BndExp[j][n]->GetCoeff_Offset(i);
                        }

                        //surface curved
                        if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                        {
                            // Sx
                            Vmath::Vvtvvtp(nfq,normals[0],1,fstress[0],1,
                                           normals[1],1,fstress[3],1,fshear[0],1);
                            Vmath::Vvtvp  (nfq,normals[2],1,fstress[4],1,fshear[0],1,fshear[0],1);

                            // Sy
                            Vmath::Vvtvvtp(nfq,normals[0],1,fstress[3],1,
                                           normals[1],1,fstress[1],1,fshear[1],1);
                            Vmath::Vvtvp  (nfq,normals[2],1,fstress[5],1,fshear[1],1,fshear[1],1);

                            // Sz
                            Vmath::Vvtvvtp(nfq,normals[0],1,fstress[4],1,
                                           normals[1],1,fstress[5],1,fshear[2],1);
                            Vmath::Vvtvp  (nfq,normals[2],1,fstress[2],1,fshear[2],1,fshear[2],1);
                        }
                        else
                        {
                            // Sx
                            Vmath::Svtsvtp(nfq,normals[0][0],fstress[0],1,
                                           normals[1][0],fstress[3],1,fshear[0],1);
                            Vmath::Svtvp(nfq,normals[2][0],fstress[4],1,fshear[0],1,fshear[0],1);

                            // Sy
                            Vmath::Svtsvtp(nfq,normals[0][0],fstress[3],1,
                                           normals[1][0],fstress[1],1,fshear[1],1);
                            Vmath::Svtvp(nfq,normals[2][0],fstress[5],1,fshear[1],1,fshear[1],1);

                            // Sz
                            Vmath::Svtsvtp(nfq,normals[0][0],fstress[4],1,
                                           normals[1][0],fstress[5],1,fshear[2],1);
                            Vmath::Svtvp(nfq,normals[2][0],fstress[2],1,fshear[2],1,fshear[2],1);
                        }

                        // T = T - (T.n)n
                        if (m_metricinfo->GetGtype() == SpatialDomains::eDeformed)
                        {
                            Vmath::Vvtvvtp(nfq,normals[0],1,fshear[0],1,
                                           normals[1],1, fshear[1],1,fshear[3],1);
                            Vmath::Vvtvp  (nfq,normals[2],1, fshear[2],1,fshear[3],1,fshear[3],1);
                            Vmath::Smul(nfq, -1.0, fshear[3], 1, fshear[3], 1);

                            for (j = 0; j < nfields; j++)
                            {
                                Vmath::Vvtvp(nfq,normals[j], 1, fshear[3], 1, fshear[j], 1, fshear[j], 1);
                                bc->FwdTrans(fshear[j], outfield[j]);
                            }
                        }
                        else
                        {
                            Vmath::Svtsvtp(nfq,normals[0][0],fshear[0],1,
                                           normals[1][0],fshear[1],1,fshear[3],1);
                            Vmath::Svtvp(nfq,normals[2][0],fshear[2],1,fshear[3],1,fshear[3],1);
                            Vmath::Smul(nfq, -1.0, fshear[3], 1,fshear[3], 1);

                            for (j = 0; j < nfields; j++)
                            {
                                Vmath::Svtvp(nfq,normals[j][0],fshear[3],1,fshear[j],1,fshear[j],1);
                                bc->FwdTrans(fshear[j], outfield[j]);
                            }
                        }

                        // Tw
                        Vmath::Vvtvvtp(nfq, fshear[0], 1, fshear[0], 1, fshear[1], 1, fshear[1], 1, fshear[3], 1);
                        Vmath::Vvtvp(nfq, fshear[2], 1, fshear[2], 1, fshear[3], 1, fshear[3], 1);
                        Vmath::Vsqrt(nfq, fshear[3], 1, fshear[3], 1);
                        bc->FwdTrans(fshear[3], outfield[3]);

                    }
                }
            }
        }
        if(doneBnd == false)
        {
            cnt += BndExp[0][n]->GetExpSize();
        }
    }


    for(int j = 0; j < newfields; ++j)
    {
        for(int b = 0; b < m_f->m_bndRegionsToWrite.size(); ++b)
        {
            m_f->m_exp[j]->UpdateBndCondExpansion(m_f->m_bndRegionsToWrite[b]) = BndExp[j][m_f->m_bndRegionsToWrite[b]];
        }
    }
}

}
}
