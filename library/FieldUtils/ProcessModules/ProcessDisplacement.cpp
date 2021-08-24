////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessDisplacement.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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
//  Description: Deforms a mesh given input field defining displacement.
//
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include "ProcessDisplacement.h"

#include <LibUtilities/BasicUtils/HashUtils.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LocalRegions/SegExp.h>
#include <LocalRegions/TriExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdTriExp.h>

namespace Nektar
{
namespace FieldUtils
{
struct TriFaceIDs
{
    TriFaceIDs(int a, int b, int c) : a(a), b(b), c(c)
    {
    }
    int a;
    int b;
    int c;
};

struct TriFaceHash : std::unary_function<TriFaceIDs, std::size_t>
{
    std::size_t operator()(TriFaceIDs const &p) const
    {
        std::vector<int> ids(3);

        ids[0] = p.a;
        ids[1] = p.b;
        ids[2] = p.c;

        std::sort(ids.begin(), ids.end());
        return hash_combine(ids[0], ids[1], ids[2]);
    }
};

bool operator==(TriFaceIDs const &p1, TriFaceIDs const &p2)
{
    std::vector<int> ids1(3), ids2(3);

    ids1[0] = p1.a;
    ids1[1] = p1.b;
    ids1[2] = p1.c;
    ids2[0] = p2.a;
    ids2[1] = p2.b;
    ids2[2] = p2.c;

    std::sort(ids1.begin(), ids1.end());
    std::sort(ids2.begin(), ids2.end());

    return ids1[0] == ids2[0] && ids1[1] == ids2[1] && ids1[2] == ids2[2];
}

typedef std::unordered_map<TriFaceIDs, int, TriFaceHash> TriFaceMap;

ModuleKey ProcessDisplacement::className =
    GetModuleFactory().RegisterCreatorFunction(
        ModuleKey(eProcessModule, "displacement"),
        ProcessDisplacement::create,
        "Deform a mesh given an input field defining displacement");

ProcessDisplacement::ProcessDisplacement(FieldSharedPtr f)
    : ProcessBoundaryExtract(f)
{
    m_config["to"] =
        ConfigOption(false, "", "Name of file containing high order boundary");

    m_config["usevertexids"] = ConfigOption(
        true, "0", "Use vertex IDs instead of face IDs for matching");
}

ProcessDisplacement::~ProcessDisplacement()
{
}

void ProcessDisplacement::Process(po::variables_map &vm)
{
    ProcessBoundaryExtract::Process(vm);
    ASSERTL0( !boost::iequals(m_config["bnd"].as<string>(), "All"),
        "ProcessDisplacement needs bnd parameter with a single id.");

    string toFile = m_config["to"].as<string>();

    if (toFile == "")
    {
        cout << "ProcessDisplacement: you must provide a file" << endl;
        return;
    }

    bool useVertexIds = m_config["usevertexids"].as<bool>();

    vector<string> files;
    files.push_back(toFile);
    LibUtilities::SessionReaderSharedPtr bndSession =
        LibUtilities::SessionReader::CreateInstance(0, NULL, files);
    SpatialDomains::MeshGraphSharedPtr bndGraph =
        SpatialDomains::MeshGraph::Read(bndSession);

    // Try to find boundary condition expansion.
    int bndCondId = m_config["bnd"].as<int>();

    // FIXME: We should be storing boundary condition IDs
    // somewhere...
    if (bndGraph->GetMeshDimension() == 1)
    {
        m_f->m_exp.push_back(m_f->AppendExpList(0, "v"));
        m_f->m_variables.push_back("v");

        MultiRegions::ExpListSharedPtr bndCondExpU =
            m_f->m_exp[0]->GetBndCondExpansions()[bndCondId];
        MultiRegions::ExpListSharedPtr bndCondExpV =
            m_f->m_exp[1]->GetBndCondExpansions()[bndCondId];

        map<int, int> bndCondIds;
        for (int i = 0; i < bndCondExpU->GetExpSize(); ++i)
        {
            bndCondIds[bndCondExpU->GetExp(i)->GetGeom()->GetGlobalID()] = i;
        }

        const SpatialDomains::SegGeomMap &tmp = bndGraph->GetAllSegGeoms();

        for (auto &sIt : tmp)
        {
            auto mIt = bndCondIds.find(sIt.first);

            if (mIt == bndCondIds.end())
            {
                cout << "Warning: couldn't find element " << sIt.first << endl;
                continue;
            }

            int e = mIt->second;

            SpatialDomains::SegGeomSharedPtr from =
                std::dynamic_pointer_cast<SpatialDomains::SegGeom>(
                    bndCondExpU->GetExp(e)->GetGeom());

            SpatialDomains::SegGeomSharedPtr to = sIt.second;

            // Create temporary SegExp
            LocalRegions::SegExpSharedPtr toSeg =
                MemoryManager<LocalRegions::SegExp>::AllocateSharedPtr(
                    bndCondExpU->GetExp(e)->GetBasis(0)->GetBasisKey(), to);

            const int offset = bndCondExpU->GetPhys_Offset(e);
            const int nq     = toSeg->GetTotPoints();

            Array<OneD, NekDouble> xL(nq), xC(nq), yL(nq), yC(nq), tmp;

            bndCondExpU->GetExp(e)->GetCoords(xC, yC);
            toSeg->GetCoords(xL, yL);

            Vmath::Vsub(nq, xL, 1, xC, 1,
                        tmp = bndCondExpU->UpdatePhys() + offset, 1);
            Vmath::Vsub(nq, yL, 1, yC, 1,
                        tmp = bndCondExpV->UpdatePhys() + offset, 1);
        }

        // bndconstrained?
        bndCondExpU->FwdTrans_BndConstrained(bndCondExpU->GetPhys(),
                                             bndCondExpU->UpdateCoeffs());
        bndCondExpV->FwdTrans_BndConstrained(bndCondExpV->GetPhys(),
                                             bndCondExpV->UpdateCoeffs());
    }
    else if (bndGraph->GetMeshDimension() == 2)
    {
        m_f->m_exp.push_back(m_f->AppendExpList(0, "v"));
        m_f->m_exp.push_back(m_f->AppendExpList(0, "w"));
        m_f->m_variables.push_back("v");
        m_f->m_variables.push_back("w");

        MultiRegions::ExpListSharedPtr bndCondExpU =
            m_f->m_exp[0]->GetBndCondExpansions()[bndCondId];
        MultiRegions::ExpListSharedPtr bndCondExpV =
            m_f->m_exp[1]->GetBndCondExpansions()[bndCondId];
        MultiRegions::ExpListSharedPtr bndCondExpW =
            m_f->m_exp[2]->GetBndCondExpansions()[bndCondId];

        map<int, int> bndCondIds;
        for (int i = 0; i < bndCondExpU->GetExpSize(); ++i)
        {
            bndCondIds[bndCondExpU->GetExp(i)->GetGeom()->GetGlobalID()] = i;
        }

        TriFaceMap vertexFaceMap;

        if (useVertexIds)
        {
            for (int i = 0; i < bndCondExpU->GetExpSize(); ++i)
            {
                SpatialDomains::TriGeomSharedPtr from =
                    std::dynamic_pointer_cast<SpatialDomains::TriGeom>(
                        bndCondExpU->GetExp(i)->GetGeom());

                TriFaceIDs t(from->GetVid(0), from->GetVid(1), from->GetVid(2));
                vertexFaceMap[t] = i;
            }
        }

        const SpatialDomains::TriGeomMap &tmp = bndGraph->GetAllTriGeoms();

        for (auto &sIt : tmp)
        {
            int e;

            if (useVertexIds)
            {
                TriFaceIDs t(sIt.second->GetVid(0), sIt.second->GetVid(1),
                             sIt.second->GetVid(2));

                auto tIt = vertexFaceMap.find(t);
                e = tIt == vertexFaceMap.end() ? -1 : tIt->second;
            }
            else
            {
                auto mIt = bndCondIds.find(sIt.first);
                e = mIt == bndCondIds.end() ? -1 : mIt->second;
            }

            if (e == -1)
            {
                cout << "Warning: couldn't find element " << sIt.first << endl;
                continue;
            }

            SpatialDomains::TriGeomSharedPtr from =
                std::dynamic_pointer_cast<SpatialDomains::TriGeom>(
                    bndCondExpU->GetExp(e)->GetGeom());

            SpatialDomains::TriGeomSharedPtr to = sIt.second;

            // Create temporary SegExp
            LocalRegions::TriExpSharedPtr toSeg =
                MemoryManager<LocalRegions::TriExp>::AllocateSharedPtr(
                    bndCondExpU->GetExp(e)->GetBasis(0)->GetBasisKey(),
                    bndCondExpV->GetExp(e)->GetBasis(1)->GetBasisKey(), to);

            const int offset = bndCondExpU->GetPhys_Offset(e);
            const int nq     = toSeg->GetTotPoints();

            Array<OneD, NekDouble> xL(nq), xC(nq), yL(nq), yC(nq), tmp;
            Array<OneD, NekDouble> zL(nq), zC(nq);

            bndCondExpU->GetExp(e)->GetCoords(xC, yC, zC);
            toSeg->GetCoords(xL, yL, zL);

            Vmath::Vsub(nq, xL, 1, xC, 1,
                        tmp = bndCondExpU->UpdatePhys() + offset, 1);
            Vmath::Vsub(nq, yL, 1, yC, 1,
                        tmp = bndCondExpV->UpdatePhys() + offset, 1);
            Vmath::Vsub(nq, zL, 1, zC, 1,
                        tmp = bndCondExpW->UpdatePhys() + offset, 1);
        }

        // bndconstrained?
        bndCondExpU->FwdTrans_BndConstrained(bndCondExpU->GetPhys(),
                                             bndCondExpU->UpdateCoeffs());
        bndCondExpV->FwdTrans_BndConstrained(bndCondExpV->GetPhys(),
                                             bndCondExpV->UpdateCoeffs());
        bndCondExpW->FwdTrans_BndConstrained(bndCondExpW->GetPhys(),
                                             bndCondExpW->UpdateCoeffs());
    }
}
}
}
