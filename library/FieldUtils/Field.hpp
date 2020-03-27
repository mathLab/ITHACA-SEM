////////////////////////////////////////////////////////////////////////////////
//
//  File: Field.hpp
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
//  Description: Field converter module base classes.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_FIELD
#define FIELDUTILS_FIELD

#include <memory>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/PtsField.h>
#include <LibUtilities/BasicUtils/PtsIO.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <SpatialDomains/MeshGraph.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>

#include "FieldUtilsDeclspec.h"

namespace Nektar
{
namespace FieldUtils
{

struct Field
{
    FIELD_UTILS_EXPORT Field()
        : m_verbose(false), m_declareExpansionAsContField(false),
          m_declareExpansionAsDisContField(false),
          m_requireBoundaryExpansion(false), m_writeBndFld(false),
          m_addNormals(false), m_fieldPts(LibUtilities::NullPtsField)
    {
    }

    FIELD_UTILS_EXPORT ~Field()
    {
        if (m_comm)
        {
            m_comm->Finalise();
        }
    }
    bool m_verbose;
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> m_fielddef;
    std::vector<std::vector<double> > m_data;
    std::vector<MultiRegions::ExpListSharedPtr> m_exp;
    std::vector<std::string> m_variables;

    int m_numHomogeneousDir;

    bool m_declareExpansionAsContField;
    bool m_declareExpansionAsDisContField;

    bool m_requireBoundaryExpansion;

    bool m_useFFT;

    LibUtilities::CommSharedPtr m_comm;
    LibUtilities::SessionReaderSharedPtr m_session;
    SpatialDomains::MeshGraphSharedPtr m_graph;
    std::map<std::string, std::vector<std::string> > m_inputfiles;

    bool m_writeBndFld;
    std::vector<unsigned int> m_bndRegionsToWrite;
    bool m_addNormals;

    LibUtilities::PtsFieldSharedPtr m_fieldPts;

    LibUtilities::FieldMetaDataMap m_fieldMetaDataMap;

    FIELD_UTILS_EXPORT MultiRegions::ExpListSharedPtr SetUpFirstExpList(
        int NumHomogeneousDir, bool fldfilegiven = false)
    {

        MultiRegions::ExpListSharedPtr exp;

        // Set up expansion list
        int expdim      = m_graph->GetMeshDimension();
        bool dealiasing = false;

        m_session->MatchSolverInfo("USEFFT", "FFTW", m_useFFT, false);

        switch (expdim)
        {
            case 1:
            {
                ASSERTL0(NumHomogeneousDir <= 2,
                         "Quasi-3D approach is only set up for 1 or 2 "
                         "homogeneous directions");

                if (NumHomogeneousDir == 1)
                {
                    MultiRegions::ExpList2DHomogeneous1DSharedPtr Exp2DH1;

                    // Define Homogeneous expansion
                    int nplanes;
                    NekDouble ly;
                    LibUtilities::BasisType btype;

                    if (fldfilegiven)
                    {
                        nplanes = m_fielddef[0]->m_numModes[1];
                        ly      = m_fielddef[0]->m_homogeneousLengths[0];
                        btype   = m_fielddef[0]->m_basis[1];
                    }
                    else
                    {
                        m_session->LoadParameter("HomModesZ", nplanes);
                        m_session->LoadParameter("LY", ly);
                        btype = LibUtilities::eFourier;
                    }

                    // Choose points to be at evenly spaced points at
                    // nplanes points
                    const LibUtilities::PointsKey Pkey(
                        nplanes, LibUtilities::eFourierEvenlySpaced);

                    const LibUtilities::BasisKey Bkey(btype, nplanes, Pkey);

                    if (m_declareExpansionAsContField ||
                        m_declareExpansionAsDisContField)
                    {
                        ASSERTL0(false, "ContField2DHomogeneous1D or "
                                        "DisContField2DHomogenenous1D has "
                                        "not been implemented");
                    }

                    Exp2DH1 =
                        MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::
                            AllocateSharedPtr(m_session, Bkey, ly, m_useFFT,
                                              dealiasing, m_graph,
                                              Collections::eNoCollection);
                    exp = Exp2DH1;
                }
                else if (NumHomogeneousDir == 2)
                {
                    MultiRegions::ExpList3DHomogeneous2DSharedPtr Exp3DH2;

                    int nylines, nzlines;
                    NekDouble ly, lz;
                    LibUtilities::BasisType btype1, btype2;

                    if (fldfilegiven)
                    {
                        nylines = m_fielddef[0]->m_numModes[1];
                        nzlines = m_fielddef[0]->m_numModes[2];
                        ly      = m_fielddef[0]->m_homogeneousLengths[0];
                        lz      = m_fielddef[0]->m_homogeneousLengths[1];
                        btype1  = m_fielddef[0]->m_basis[1];
                        btype2  = m_fielddef[0]->m_basis[2];
                    }
                    else
                    {
                        m_session->LoadParameter("HomModesY", nylines);
                        m_session->LoadParameter("HomModesZ", nzlines);
                        m_session->LoadParameter("LY", ly);
                        m_session->LoadParameter("LZ", lz);
                        btype1 = LibUtilities::eFourier;
                        btype2 = LibUtilities::eFourier;
                    }

                    // Choose points to be at evenly spaced points at
                    // nplanes points
                    const LibUtilities::PointsKey PkeyY(
                        nylines, LibUtilities::eFourierEvenlySpaced);
                    const LibUtilities::BasisKey BkeyY(btype1, nylines, PkeyY);

                    const LibUtilities::PointsKey PkeyZ(
                        nzlines, LibUtilities::eFourierEvenlySpaced);
                    const LibUtilities::BasisKey BkeyZ(btype2, nzlines, PkeyZ);

                    if (m_declareExpansionAsContField)
                    {
                        Exp3DH2 = MemoryManager<
                            MultiRegions::ContField3DHomogeneous2D>::
                            AllocateSharedPtr(m_session, BkeyY, BkeyZ, ly, lz,
                                              m_useFFT, dealiasing, m_graph,
                                              m_session->GetVariable(0),
                                              Collections::eNoCollection);
                    }
                    else if (m_declareExpansionAsDisContField)
                    {
                        Exp3DH2 = MemoryManager<
                            MultiRegions::DisContField3DHomogeneous2D>::
                            AllocateSharedPtr(m_session, BkeyY, BkeyZ, ly, lz,
                                              m_useFFT, dealiasing, m_graph,
                                              m_session->GetVariable(0),
                                              Collections::eNoCollection);
                    }
                    else
                    {
                        Exp3DH2 = MemoryManager<
                            MultiRegions::ExpList3DHomogeneous2D>::
                            AllocateSharedPtr(m_session, BkeyY, BkeyZ, ly, lz,
                                              m_useFFT, dealiasing, m_graph,
                                              Collections::eNoCollection);
                    }

                    exp = Exp3DH2;
                }
                else
                {
                    MultiRegions::ExpList1DSharedPtr Exp1D;

                    if (m_declareExpansionAsContField)
                    {
                        Exp1D = MemoryManager<MultiRegions::ContField1D>::
                            AllocateSharedPtr(m_session, m_graph,
                                              m_session->GetVariable(0),
                                              Collections::eNoCollection);
                    }
                    else if (m_declareExpansionAsDisContField)
                    {
                        Exp1D = MemoryManager<MultiRegions::DisContField1D>::
                            AllocateSharedPtr(m_session, m_graph,
                                              m_session->GetVariable(0),
                                              true,
                                              Collections::eNoCollection);
                    }
                    else
                    {
                        Exp1D = MemoryManager<MultiRegions::ExpList1D>::
                            AllocateSharedPtr(m_session, m_graph,
                                              true,
                                              Collections::eNoCollection);
                    }

                    exp = Exp1D;
                }
            }
            break;
            case 2:
            {
                ASSERTL0(NumHomogeneousDir <= 1,
                         "NumHomogeneousDir is only set up for 1");

                if (NumHomogeneousDir == 1)
                {
                    MultiRegions::ExpList3DHomogeneous1DSharedPtr Exp3DH1;

                    // Define Homogeneous expansion
                    int nplanes;
                    NekDouble lz;
                    LibUtilities::BasisType  btype;
                    LibUtilities::PointsType ptype =
                            LibUtilities::eFourierEvenlySpaced;

                    if (fldfilegiven)
                    {
                        nplanes = m_fielddef[0]->m_numModes[2];
                        lz      = m_fielddef[0]->m_homogeneousLengths[0];
                        btype   = m_fielddef[0]->m_basis[2];

                        if (btype == LibUtilities::eFourierSingleMode)
                        {
                            btype = LibUtilities::eFourier;
                            m_fielddef[0]->m_basis[2] =
                                LibUtilities::eFourierSingleMode;
                            if (nplanes <= 2)
                            {
                                nplanes = 4;
                            }
                        }
                        else if (btype == LibUtilities::eFourierHalfModeRe &&
                                 nplanes == 1)
                        {
                            ptype = LibUtilities::ePolyEvenlySpaced;
                        }
                    }
                    else
                    {
                        m_session->LoadParameter("HomModesZ", nplanes);
                        m_session->LoadParameter("LZ", lz);
                        btype = LibUtilities::eFourier;
                    }

                    // Choose points to be at evenly spaced points at
                    // nplanes points
                    const LibUtilities::PointsKey Pkey(
                        nplanes, ptype);

                    const LibUtilities::BasisKey Bkey(btype, nplanes, Pkey);

                    if (m_declareExpansionAsContField)
                    {
                        Exp3DH1 = MemoryManager<
                            MultiRegions::ContField3DHomogeneous1D>::
                            AllocateSharedPtr(m_session, Bkey, lz, m_useFFT,
                                              dealiasing, m_graph,
                                              m_session->GetVariable(0),
                                              Collections::eNoCollection);
                    }
                    else if (m_declareExpansionAsDisContField)
                    {
                        Exp3DH1 = MemoryManager<
                            MultiRegions::DisContField3DHomogeneous1D>::
                            AllocateSharedPtr(m_session, Bkey, lz, m_useFFT,
                                              dealiasing, m_graph,
                                              m_session->GetVariable(0),
                                              Collections::eNoCollection);
                    }
                    else
                    {
                        Exp3DH1 = MemoryManager<
                            MultiRegions::ExpList3DHomogeneous1D>::
                            AllocateSharedPtr(m_session, Bkey, lz, m_useFFT,
                                              dealiasing, m_graph,
                                              "DefaultVar",
                                              Collections::eNoCollection);
                    }
                    exp = Exp3DH1;
                }
                else
                {
                    MultiRegions::ExpList2DSharedPtr Exp2D;

                    if (m_declareExpansionAsContField)
                    {
                        Exp2D = MemoryManager<MultiRegions::ContField2D>::
                            AllocateSharedPtr(m_session, m_graph,
                                              m_session->GetVariable(0),
                                              true,false,
                                              Collections::eNoCollection);
                    }
                    else if (m_declareExpansionAsDisContField)
                    {
                        Exp2D = MemoryManager<MultiRegions::DisContField2D>::
                            AllocateSharedPtr(m_session, m_graph,
                                              m_session->GetVariable(0),
                                              true,true,
                                              Collections::eNoCollection);
                    }
                    else
                    {
                        Exp2D = MemoryManager<MultiRegions::ExpList2D>::
                            AllocateSharedPtr(m_session, m_graph,
                                              true, 
                                              "DefaultVar",
                                              Collections::eNoCollection);
                    }

                    exp = Exp2D;
                }
            }
            break;
            case 3:
            {
                MultiRegions::ExpList3DSharedPtr Exp3D;

                if (m_declareExpansionAsContField)
                {
                    Exp3D = MemoryManager<MultiRegions::ContField3D>::
                        AllocateSharedPtr(m_session, m_graph,
                                          m_session->GetVariable(0),
                                          false,
                                          Collections::eNoCollection);
                }
                else if (m_declareExpansionAsDisContField)
                {
                    Exp3D = MemoryManager<MultiRegions::DisContField3D>::
                        AllocateSharedPtr(m_session, m_graph,
                                          m_session->GetVariable(0),
                                          true,
                                          Collections::eNoCollection);
                }
                else
                {
                    Exp3D = MemoryManager<
                        MultiRegions::ExpList3D>::AllocateSharedPtr(
                                               m_session,
                                               m_graph,
                                               "DefaultVar",
                                               Collections::eNoCollection);
                }

                exp = Exp3D;
            }
            break;
            default:
                ASSERTL0(false, "Expansion dimension not recognised");
                break;
        }

        return exp;
    };

    /**
     * @brief Construct a FieldIO object for the file @p filename.
     *
     * This routine constructs an appropriate FieldIO object for a filename
     * through the LibUtilities::FieldIO::GetFileType function to detect the
     * file format. The result is then cached in Field::m_fld to avoid needing
     * to repeatedly construct the object.
     *
     * @param filename  Filename to open.
     * @return Reader for @p filename.
     */
    FIELD_UTILS_EXPORT LibUtilities::FieldIOSharedPtr FieldIOForFile(
        std::string filename)
    {
        LibUtilities::CommSharedPtr c = m_comm;
        std::string fmt = LibUtilities::FieldIO::GetFileType(filename, c);
        auto it = m_fld.find(fmt);

        if (it == m_fld.end())
        {
            LibUtilities::FieldIOSharedPtr fld =
                LibUtilities::GetFieldIOFactory().CreateInstance(fmt, c, true);
            m_fld[fmt] = fld;
            return fld;
        }
        else
        {
            return it->second;
        }
    }

    FIELD_UTILS_EXPORT MultiRegions::ExpListSharedPtr AppendExpList(
        int NumHomogeneousDir,
        std::string var = "DefaultVar",
        bool NewField = false)
    {
        if (var.compare("DefaultVar") == 0 && m_requireBoundaryExpansion)
        {
            if (m_session->GetVariables().size())
            {
                var = m_session->GetVariables()[0];
            }
        }
        MultiRegions::ExpListSharedPtr tmp;
        switch (m_graph->GetMeshDimension())
        {
            case 1:
            {
                if (NumHomogeneousDir == 1)
                {
                    ASSERTL0(!(m_declareExpansionAsContField ||
                               m_declareExpansionAsDisContField),
                             "ContField2DHomogeneous1D or "
                             "DisContField2DHomogenenous1D has not been "
                             "implemented");

                    MultiRegions::ExpList2DHomogeneous1DSharedPtr tmp2 =
                        std::dynamic_pointer_cast<
                            MultiRegions::ExpList2DHomogeneous1D>(m_exp[0]);

                    tmp = MemoryManager<MultiRegions::ExpList2DHomogeneous1D>::
                        AllocateSharedPtr(*tmp2);
                }
                else if (NumHomogeneousDir == 2)
                {
                    if (m_declareExpansionAsContField)
                    {
                        MultiRegions::ContField3DHomogeneous2DSharedPtr tmp2 =
                            std::dynamic_pointer_cast<
                                MultiRegions::ContField3DHomogeneous2D>(
                                m_exp[0]);

                        tmp = MemoryManager<
                            MultiRegions::ContField3DHomogeneous2D>::
                            AllocateSharedPtr(*tmp2);
                    }
                    else if (m_declareExpansionAsDisContField)
                    {
                        MultiRegions::DisContField3DHomogeneous2DSharedPtr
                            tmp2 = std::dynamic_pointer_cast<
                                MultiRegions::DisContField3DHomogeneous2D>(
                                m_exp[0]);

                        tmp = MemoryManager<
                            MultiRegions::DisContField3DHomogeneous2D>::
                            AllocateSharedPtr(*tmp2);
                    }
                    else
                    {
                        MultiRegions::ExpList3DHomogeneous2DSharedPtr tmp2 =
                            std::dynamic_pointer_cast<
                                MultiRegions::ExpList3DHomogeneous2D>(m_exp[0]);

                        tmp = MemoryManager<
                            MultiRegions::ExpList3DHomogeneous2D>::
                            AllocateSharedPtr(*tmp2);
                    }
                }
                else
                {
                    if (m_declareExpansionAsContField)
                    {
                        MultiRegions::ContField1DSharedPtr tmp2 =
                            std::dynamic_pointer_cast<
                                MultiRegions::ContField1D>(m_exp[0]);

                        tmp = MemoryManager<MultiRegions::ContField1D>::
                            AllocateSharedPtr(m_session, m_graph, var);
                    }
                    else if (m_declareExpansionAsDisContField)
                    {
                        MultiRegions::DisContField1DSharedPtr tmp2 =
                            std::dynamic_pointer_cast<
                                MultiRegions::DisContField1D>(m_exp[0]);

                        tmp = MemoryManager<MultiRegions::DisContField1D>::
                            AllocateSharedPtr(m_session, m_graph, var);
                    }
                    else
                    {
                        MultiRegions::ExpList1DSharedPtr tmp2 =
                            std::dynamic_pointer_cast<
                                MultiRegions::ExpList1D>(m_exp[0]);

                        tmp = MemoryManager<
                            MultiRegions::ExpList1D>::AllocateSharedPtr(*tmp2);
                    }
                }
            }
            break;
            case 2:
            {
                if (NumHomogeneousDir == 1)
                {
                    if (m_declareExpansionAsContField)
                    {
                        if (NewField)
                        {
                            bool dealiasing = false;

                            tmp = MemoryManager<
                                MultiRegions::ContField3DHomogeneous1D>::
                                AllocateSharedPtr(
                                    m_session, m_exp[0]
                                                   ->GetHomogeneousBasis()
                                                   ->GetBasisKey(),
                                    m_exp[0]->GetHomoLen(), m_useFFT,
                                    dealiasing, m_graph, var);
                        }
                        else
                        {
                            MultiRegions::ContField3DHomogeneous1DSharedPtr
                                tmp2 = std::dynamic_pointer_cast<
                                    MultiRegions::ContField3DHomogeneous1D>(
                                    m_exp[0]);

                            ASSERTL0(tmp2, "Failed to type cast m_exp[0]");
                            tmp = MemoryManager<
                                MultiRegions::ContField3DHomogeneous1D>::
                                AllocateSharedPtr(*tmp2, m_graph, var);
                        }
                    }
                    else if (m_declareExpansionAsDisContField)
                    {
                        if (NewField)
                        {
                            bool dealiasing = false;

                            tmp = MemoryManager<
                                MultiRegions::DisContField3DHomogeneous1D>::
                                AllocateSharedPtr(
                                    m_session, m_exp[0]
                                                   ->GetHomogeneousBasis()
                                                   ->GetBasisKey(),
                                    m_exp[0]->GetHomoLen(), m_useFFT,
                                    dealiasing, m_graph, var);
                        }
                        else
                        {
                            MultiRegions::DisContField3DHomogeneous1DSharedPtr
                                tmp2 = std::dynamic_pointer_cast<
                                    MultiRegions::DisContField3DHomogeneous1D>(
                                    m_exp[0]);
                            ASSERTL0(tmp2, "Failed to type cast m_exp[0]");

                            tmp = MemoryManager<
                                MultiRegions::DisContField3DHomogeneous1D>::
                                AllocateSharedPtr(*tmp2);
                        }
                    }
                    else
                    {
                        if (NewField)
                        {
                            bool dealiasing = false;

                            tmp = MemoryManager<
                                MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(
                                    m_session, m_exp[0]
                                                   ->GetHomogeneousBasis()
                                                   ->GetBasisKey(),
                                    m_exp[0]->GetHomoLen(), m_useFFT,
                                    dealiasing, m_graph);
                        }
                        else
                        {
                            MultiRegions::ExpList3DHomogeneous1DSharedPtr tmp2 =
                                std::dynamic_pointer_cast<
                                    MultiRegions::ExpList3DHomogeneous1D>(
                                    m_exp[0]);
                            ASSERTL0(tmp2, "Failed to type cast m_exp[0]");

                            tmp = MemoryManager<
                                MultiRegions::ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*tmp2);
                        }
                    }
                }
                else
                {
                    if (m_declareExpansionAsContField)
                    {
                        if (NewField)
                        {
                            tmp = MemoryManager<MultiRegions::ContField2D>::
                                AllocateSharedPtr(m_session, m_graph, var);
                        }
                        else // call copy constructor
                        {

                            MultiRegions::ContField2DSharedPtr tmp2 =
                                std::dynamic_pointer_cast<
                                    MultiRegions::ContField2D>(m_exp[0]);

                            tmp = MemoryManager<MultiRegions::ContField2D>::
                                AllocateSharedPtr(*tmp2, m_graph, var);
                        }
                    }
                    else if (m_declareExpansionAsDisContField)
                    {
                        if (NewField)
                        {
                            tmp = MemoryManager<MultiRegions::DisContField2D>::
                                AllocateSharedPtr(m_session, m_graph, var);
                        }
                        else // call copy constructor
                        {
                            MultiRegions::DisContField2DSharedPtr tmp2 =
                                std::dynamic_pointer_cast<
                                    MultiRegions::DisContField2D>(m_exp[0]);

                            tmp = MemoryManager<MultiRegions::DisContField2D>::
                                AllocateSharedPtr(*tmp2, m_graph, var);
                        }
                    }
                    else
                    {
                        MultiRegions::ExpList2DSharedPtr tmp2 =
                            std::dynamic_pointer_cast<
                                MultiRegions::ExpList2D>(m_exp[0]);

                        tmp = MemoryManager<
                            MultiRegions::ExpList2D>::AllocateSharedPtr(*tmp2);
                    }
                }
            }
            break;
            case 3:
            {
                if (m_declareExpansionAsContField)
                {
                    if (NewField)
                    {
                        tmp = MemoryManager<MultiRegions::ContField3D>::
                            AllocateSharedPtr(m_session, m_graph, var);
                    }
                    else
                    {
                        MultiRegions::ContField3DSharedPtr tmp2 =
                            std::dynamic_pointer_cast<
                                MultiRegions::ContField3D>(m_exp[0]);

                        tmp = MemoryManager<MultiRegions::ContField3D>::
                            AllocateSharedPtr(*tmp2, m_graph, var);
                    }
                }
                else if (m_declareExpansionAsDisContField)
                {
                    if (NewField)
                    {
                        tmp = MemoryManager<MultiRegions::DisContField3D>::
                            AllocateSharedPtr(m_session, m_graph, var);
                    }
                    else
                    {
                        MultiRegions::DisContField3DSharedPtr tmp2 =
                            std::dynamic_pointer_cast<
                                MultiRegions::DisContField3D>(m_exp[0]);

                        tmp = MemoryManager<MultiRegions::DisContField3D>::
                            AllocateSharedPtr(*tmp2, m_graph, var);
                    }
                }
                else
                {
                    MultiRegions::ExpList3DSharedPtr tmp2 =
                        std::dynamic_pointer_cast<MultiRegions::ExpList3D>(
                            m_exp[0]);

                    tmp = MemoryManager<
                        MultiRegions::ExpList3D>::AllocateSharedPtr(*tmp2);
                }
            }
            break;
            default:
                ASSERTL0(false, "Expansion dimension not recognised");
                break;
        }

        return tmp;
    }

    FIELD_UTILS_EXPORT void ClearField()
    {
        m_session  = LibUtilities::SessionReaderSharedPtr();
        m_graph    = SpatialDomains::MeshGraphSharedPtr();
        m_fieldPts = LibUtilities::NullPtsField;
        m_exp.clear();
        m_fielddef = std::vector<LibUtilities::FieldDefinitionsSharedPtr>();
        m_data     = std::vector<std::vector<NekDouble> > ();
        m_variables.clear();
    }

private:
    /// Map to store FieldIO instances. Key is the reader type, value is the
    /// FieldIO object.
    std::map<std::string, LibUtilities::FieldIOSharedPtr> m_fld;
};

typedef std::shared_ptr<Field> FieldSharedPtr;
}
}

#endif
