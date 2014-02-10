////////////////////////////////////////////////////////////////////////////////
//
//  File: FldQCriterion.cpp
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
//  Description: Add Q criterion to field.
//
////////////////////////////////////////////////////////////////////////////////

#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ExpList1D.h>
#include <MultiRegions/ExpList2D.h>
#include <MultiRegions/ExpList3D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous1D.h>
#include <MultiRegions/ExpList3DHomogeneous2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    int i, j;

    if (argc != 4)
    {
        fprintf(stderr, "Usage: FldQCriterion  meshfile infld outfld\n");
        exit(1);
    }

    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv);


    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc - 3]);
    SpatialDomains::MeshGraphSharedPtr graphShPt =
        SpatialDomains::MeshGraph::Read(vSession);
    //----------------------------------------------

    //----------------------------------------------
    // Import field file.
    string fieldfile(argv[argc - 2]);
    vector<LibUtilities::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    LibUtilities::Import(fieldfile, fielddef, fielddata);
    bool useFFT = false;
    bool dealiasing = false;
    //----------------------------------------------

    //----------------------------------------------
    // Define Expansion
    int expdim  = graphShPt->GetMeshDimension();
    int nfields = fielddef[0]->m_fields.size();
    int addfields = 1;
    Array<OneD, MultiRegions::ExpListSharedPtr> Exp(nfields + addfields);

    switch (expdim)
    {
        case 1:
        {
            ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 2,
                     "Quasi-3D approach is only set up for 1 or 2 homogeneous "
                     "directions");

            if (fielddef[0]->m_numHomogeneousDir == 1)
            {
                MultiRegions::ExpList2DHomogeneous1DSharedPtr Exp2DH1;

                // Define Homogeneous expansion
                int nplanes;
                vSession->LoadParameter(
                    "HomModesZ", nplanes, fielddef[0]->m_numModes[1]);

                // choose points to be at evenly spaced points at
                // nplanes points
                const LibUtilities::PointsKey Pkey(
                    nplanes, LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(
                    fielddef[0]->m_basis[1], nplanes, Pkey);
                NekDouble ly = fielddef[0]->m_homogeneousLengths[0];

                Exp2DH1 = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>
                    ::AllocateSharedPtr(
                        vSession, Bkey, ly, useFFT, dealiasing, graphShPt);
                Exp[0] = Exp2DH1;

                for (i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>
                        ::AllocateSharedPtr(*Exp2DH1);
                }
            }
            else if (fielddef[0]->m_numHomogeneousDir == 2)
            {
                MultiRegions::ExpList3DHomogeneous2DSharedPtr Exp3DH2;

                // Define Homogeneous expansion
                //int nylines = fielddef[0]->m_numModes[1];
                //int nzlines = fielddef[0]->m_numModes[2];

                int nylines;
                int nzlines;
                vSession->LoadParameter(
                    "HomModesY", nylines, fielddef[0]->m_numModes[1]);
                vSession->LoadParameter(
                    "HomModesZ", nzlines, fielddef[0]->m_numModes[2]);

                // choose points to be at evenly spaced points at
                // nplanes points
                const LibUtilities::PointsKey PkeyY(
                    nylines, LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  BkeyY(
                    fielddef[0]->m_basis[1], nylines, PkeyY);

                const LibUtilities::PointsKey PkeyZ(
                    nzlines, LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  BkeyZ(
                    fielddef[0]->m_basis[2], nzlines, PkeyZ);

                NekDouble ly = fielddef[0]->m_homogeneousLengths[0];
                NekDouble lz = fielddef[0]->m_homogeneousLengths[1];

                Exp3DH2 = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>
                    ::AllocateSharedPtr(vSession, BkeyY, BkeyZ, ly, lz,
                                        useFFT, dealiasing, graphShPt);
                Exp[0] = Exp3DH2;
                for (i = 1; i < nfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous2D>
                        ::AllocateSharedPtr(*Exp3DH2);
                }
            }
            else
            {
                MultiRegions::ExpList1DSharedPtr Exp1D;
                Exp1D = MemoryManager<MultiRegions::ExpList1D>
                    ::AllocateSharedPtr(vSession, graphShPt);
                Exp[0] = Exp1D;
                for (i = 1; i < nfields + addfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList1D>
                        ::AllocateSharedPtr(*Exp1D);
                }
            }
        }
        break;
        case 2:
        {

            ASSERTL0(fielddef[0]->m_numHomogeneousDir <= 1,
                     "NumHomogeneousDir is only set up for 1");

            if (fielddef[0]->m_numHomogeneousDir == 1)
            {
                MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;

                // Define Homogeneous expansion
                //int nplanes = fielddef[0]->m_numModes[2];

                int nplanes;
                vSession->LoadParameter(
                    "HomModesZ", nplanes, fielddef[0]->m_numModes[2]);

                // choose points to be at evenly spaced points at
                // nplanes points
                const LibUtilities::PointsKey Pkey(
                    nplanes, LibUtilities::ePolyEvenlySpaced);
                const LibUtilities::BasisKey  Bkey(
                    fielddef[0]->m_basis[2], nplanes, Pkey);
                NekDouble lz = fielddef[0]->m_homogeneousLengths[0];

                Exp3DH1 = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>
                    ::AllocateSharedPtr(vSession, Bkey, lz, useFFT,
                                        dealiasing, graphShPt);
                Exp[0] = Exp3DH1;

                for (i = 1; i < nfields + addfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList3DHomogeneous1D>
                        ::AllocateSharedPtr(*Exp3DH1);
                }
            }
            else
            {
                MultiRegions::ExpList2DSharedPtr Exp2D;
                Exp2D = MemoryManager<MultiRegions::ExpList2D>
                    ::AllocateSharedPtr(vSession, graphShPt);
                Exp[0] =  Exp2D;

                for (i = 1; i < nfields + addfields; ++i)
                {
                    Exp[i] = MemoryManager<MultiRegions::ExpList2D>
                        ::AllocateSharedPtr(*Exp2D);
                }
            }
        }
        break;
        case 3:
        {

            MultiRegions::ExpList3DSharedPtr Exp3D;
            Exp3D = MemoryManager<MultiRegions::ExpList3D>
                ::AllocateSharedPtr(vSession, graphShPt);
            Exp[0] =  Exp3D;

            for (i = 1; i < nfields + addfields; ++i)
            {
                Exp[i] = MemoryManager<MultiRegions::ExpList3D>
                    ::AllocateSharedPtr(*Exp3D);
            }
        }
        break;
        default:
            ASSERTL0(false, "Expansion dimension not recognised");
            break;
    }
    //----------------------------------------------

    //----------------------------------------------
    // Copy data from field file
    for (j = 0; j < nfields; ++j)
    {
        for (int i = 0; i < fielddata.size(); ++i)
        {
            Exp[j]->ExtractDataToCoeffs(fielddef [i],
                                        fielddata[i],
                                        fielddef [i]->m_fields[j],
                                        Exp[j]->UpdateCoeffs());
        }
        Exp[j]->BwdTrans(Exp[j]->GetCoeffs(), Exp[j]->UpdatePhys());
    }
    //----------------------------------------------

    //----------------------------------------------
    // Compute gradients of fields
    int nq = Exp[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble> > grad(nfields * nfields);

    Array<OneD, Array<OneD, NekDouble> > omega(nfields * nfields);
    Array<OneD, Array<OneD, NekDouble> > S    (nfields * nfields);
    Array<OneD, Array<OneD, NekDouble> > Q    (nfields * nfields);

    Array<OneD, Array<OneD, NekDouble> > outfield (addfields);
    Array<OneD, Array<OneD, NekDouble> > outfield1(addfields);
    Array<OneD, Array<OneD, NekDouble> > outfield2(addfields);
    Array<OneD, Array<OneD, NekDouble> > outfield3(addfields);

    for (i = 0; i < nfields * nfields; ++i)
    {
        grad[i] = Array<OneD, NekDouble>(nq);
    }

    for (i = 0; i < addfields; ++i)
    {
        outfield [i] = Array<OneD, NekDouble>(nq);
        outfield1[i] = Array<OneD, NekDouble>(nq);
        outfield2[i] = Array<OneD, NekDouble>(nq);
        outfield3[i] = Array<OneD, NekDouble>(nq);

        omega[i] = Array<OneD, NekDouble>(nq);
        S[i] = Array<OneD, NekDouble>(nq);
        Q[i] = Array<OneD, NekDouble>(nq);
    }

    // Calculate Gradient & Vorticity
    for (i = 0; i < nfields; ++i)
    {
        Exp[i]->PhysDeriv(
            Exp[i]->GetPhys(), grad[i * nfields], grad[i * nfields + 1],
            grad[i * nfields + 2]);
    }

    // W_x = Wy - Vz
    Vmath::Vsub(nq, grad[2 * nfields + 1], 1, grad[1 * nfields + 2], 1,
                outfield1[0], 1);
    // W_x^2
    Vmath::Vmul(nq, outfield1[0], 1, outfield1[0], 1, outfield1[0], 1);

    // W_y = Uz - Wx
    Vmath::Vsub(nq, grad[0 * nfields + 2], 1, grad[2 * nfields + 0], 1,
                outfield2[0], 1);
    // W_y^2
    Vmath::Vmul(nq, outfield2[0], 1, outfield2[0], 1, outfield2[0], 1);

    // W_z = Vx - Uy
    Vmath::Vsub(nq, grad[1 * nfields + 0], 1, grad[0 * nfields + 1], 1,
                outfield3[0], 1);
    // W_z^2
    Vmath::Vmul(nq, outfield3[0], 1, outfield3[0], 1, outfield3[0], 1);

    // add fields omega = 0.5*(W_x^2 + W_y^2 + W_z^2)

    NekDouble fac = 0.5;
    Vmath::Vadd(nq, &outfield1[0][0], 1, &outfield2[0][0], 1, &omega[0][0], 1);
    Vmath::Vadd(nq, &omega[0][0], 1, &outfield3[0][0], 1, &omega[0][0], 1);

    for (int k = 0; k < addfields; ++k)
    {
        Vmath::Smul(nq, fac, &omega[k][0], 1, &omega[k][0], 1);
    }

    Vmath::Zero(nq, &outfield1[0][0], 1);
    Vmath::Zero(nq, &outfield2[0][0], 1);
    Vmath::Zero(nq, &outfield3[0][0], 1);

    Vmath::Vmul(nq, grad[0 * nfields + 0], 1, grad[0 * nfields + 0], 1,
                outfield1[0], 1);
    Vmath::Vmul(nq, grad[1 * nfields + 1], 1, grad[1 * nfields + 1], 1,
                outfield2[0], 1);
    Vmath::Vmul(nq, grad[2 * nfields + 2], 1, grad[2 * nfields + 2], 1,
                outfield3[0], 1);

    Vmath::Vadd(nq, &outfield1[0][0], 1, &outfield2[0][0], 1, &S[0][0], 1);
    Vmath::Vadd(nq, &S[0][0], 1, &outfield3[0][0], 1, &S[0][0], 1);

    // W_y + V_z
    Vmath::Vadd(nq, grad[2 * nfields + 1], 1, grad[1 * nfields + 2], 1,
                outfield1[0], 1);
    Vmath::Vmul(nq, &outfield1[0][0], 1, &outfield1[0][0], 1,
                &outfield1[0][0], 1);

    // U_z + W_x
    Vmath::Vadd(nq, grad[0 * nfields + 2], 1, grad[2 * nfields + 0], 1,
                outfield2[0], 1);
    Vmath::Vmul(nq, &outfield2[0][0], 1, &outfield2[0][0], 1,
                &outfield2[0][0], 1);

    // V_x + U_y
    Vmath::Vadd(nq, grad[1 * nfields + 0], 1, grad[0 * nfields + 1], 1,
                outfield3[0], 1);
    Vmath::Vmul(nq, &outfield3[0][0], 1, &outfield3[0][0], 1,
                &outfield3[0][0], 1);

    Vmath::Vadd(nq, &outfield1[0][0], 1, &outfield2[0][0], 1,
                &outfield2[0][0], 1);
    Vmath::Vadd(nq, &outfield2[0][0], 1, &outfield3[0][0], 1,
                &outfield3[0][0], 1);

    for (int k = 0; k < addfields; ++k)
    {
        Vmath::Smul(nq, fac, &outfield3[k][0], 1, &outfield3[k][0], 1);
    }

    Vmath::Vadd(nq, &outfield3[0][0], 1, &S[0][0], 1, &S[0][0], 1);
    Vmath::Vsub(nq, omega[0], 1, S[0], 1, Q[0], 1);

    for (int k = 0; k < addfields; ++k)
    {
        Vmath::Smul(nq, fac, &Q[k][0], 1, &Q[k][0], 1);
    }

    for (i = 0; i < addfields; ++i)
    {
        Exp[nfields + i]->FwdTrans(Q[i], Exp[nfields + i]->UpdateCoeffs());
    }

    //-----------------------------------------------
    // Write solution to file with additional computed fields
    string out(argv[argc - 1]);
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef
        = Exp[0]->GetFieldDefinitions();
    std::vector<std::vector<NekDouble> > FieldData(FieldDef.size());

    vector<string> outname;
    outname.push_back("Q");

    for (j = 0; j < nfields + addfields; ++j)
    {
        for (i = 0; i < FieldDef.size(); ++i)
        {
            if (j >= nfields)
            {
                FieldDef[i]->m_fields.push_back(outname[j - nfields]);
            }
            else
            {
                FieldDef[i]->m_fields.push_back(fielddef[i]->m_fields[j]);
            }
            Exp[j]->AppendFieldData(FieldDef[i], FieldData[i]);
        }
    }
    LibUtilities::Write(out, FieldDef, FieldData);
    //-----------------------------------------------

    return 0;
}

