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
//  Description: Field converter module base classes.
//
////////////////////////////////////////////////////////////////////////////////

#include <boost/shared_ptr.hpp>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/MeshGraph.h>

#include <MultiRegions/ExpList.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ExpList2DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <MultiRegions/DisContField3DHomogeneous2D.h>
#include <MultiRegions/DisContField2D.h>
#include <MultiRegions/DisContField3D.h>


using namespace std;

namespace Nektar
{
namespace Utilities
{

enum PtsType{
    ePtsFile,
    ePtsLine,
    ePtsPlane,
    ePtsBlock
};

struct FieldPts
{
    FieldPts(void): m_ptsDim(0),
                    m_ptype(ePtsFile) {}

    int                                     m_ptsDim;
    int                                     m_nFields;
    Array<OneD, Array<OneD, NekDouble> >    m_pts;
    PtsType                                 m_ptype;
    vector<int>                             m_npts;
    vector<std::string>                     m_fields;


    // Interpolate field_id (which is the id after the coordinates)
    void Interp1DPts(
            const NekDouble          coord,
            Array<OneD, NekDouble > &intfields)
    {
        // currently assume first field is coordinate
        WARNINGL1(m_ptsDim == 1,
                 "Assumed only one coordinate given taking first coordinate "
                 "for interpolation");
        int npts = m_pts[0].num_elements();
        int i;

        for(i = 0; i < npts-1; ++i)
        {
            if((m_pts[0][i] <= coord) && (coord <= m_pts[0][i+1]))
            {
                NekDouble pdiff = m_pts[0][i+1]-m_pts[0][i];

                if(npts <= 2)
                {
                    //  linear interpolation
                    for(int j = 0; j < m_nFields; ++j)
                    {
                        intfields[j] = m_pts[m_ptsDim+j][i]
                                        * (m_pts[0][i+1] - coord) / pdiff
                                     + m_pts[m_ptsDim+j][i+1]
                                        * (coord - m_pts[0][i]) / pdiff;
                    }
                }
                else // quadratic interpolation
                {
                    if(i < npts-2)
                    { // forwards stencil
                        NekDouble pdiff2 = m_pts[0][i+2] - m_pts[0][i+1];

                        NekDouble h1 = (m_pts[0][i+1]-coord)
                                        * (m_pts[0][i+2] - coord)
                                        / (pdiff * (pdiff+pdiff2));
                        NekDouble h2 = (coord-m_pts[0][i])
                                        * (m_pts[0][i+2] - coord)
                                        / (pdiff * pdiff2);
                        NekDouble h3 = (coord-m_pts[0][i])
                                        * (coord - m_pts[0][i+1])
                                        / ((pdiff + pdiff2) * pdiff2);
                        for(int j = 0; j < m_nFields; ++j)
                        {
                            intfields[j] = m_pts[m_ptsDim+j][i] * h1
                                        +  m_pts[m_ptsDim+j][i+1] * h2
                                        +  m_pts[m_ptsDim+j][i+2] * h3;
                        }
                    }
                    else
                    { // backwards stencil
                        NekDouble pdiff2 = m_pts[0][i] - m_pts[0][i-1];

                        NekDouble h1 = (m_pts[0][i+1]-coord)
                                        * (coord - m_pts[0][i-1])
                                        / (pdiff * pdiff2);
                        NekDouble h2 = (coord - m_pts[0][i])
                                        * (coord - m_pts[0][i-1])
                                        / (pdiff * (pdiff + pdiff2));
                        NekDouble h3 = (m_pts[0][i]-coord)
                                        * (m_pts[0][i+1] - coord)
                                        / ((pdiff + pdiff2) * pdiff);
                        for(int j = 0; j < m_nFields; ++j)
                        {
                            intfields[j] = m_pts[m_ptsDim+j][i] * h1
                                        +  m_pts[m_ptsDim+j][i+1] * h2
                                        +  m_pts[m_ptsDim+j][i-1] * h3;
                        }
                    }
                }
                break;
            }
        }
        ASSERTL0(i != npts-1, "Failed to find coordinate " +
                              boost::lexical_cast<string>(coord) +
                              " within provided input points");
    };

};

typedef boost::shared_ptr<FieldPts> FieldPtsSharedPtr;
static FieldPtsSharedPtr NullFieldPts;

struct Field {
    Field() : m_verbose(false),
              m_declareExpansionAsContField(false),
              m_declareExpansionAsDisContField(false),
              m_writeBndFld(false),
              m_fieldPts(NullFieldPts){}

    ~Field()
    {
        if (m_session)
        {
            m_session->Finalise();
        }
    }
    bool                                    m_verbose;
    vector<LibUtilities::FieldDefinitionsSharedPtr> m_fielddef;
    vector<vector<double> >                 m_data;
    vector<MultiRegions::ExpListSharedPtr>  m_exp;

    bool                                    m_declareExpansionAsContField;
    bool                                    m_declareExpansionAsDisContField;

    LibUtilities::CommSharedPtr             m_comm;
    LibUtilities::SessionReaderSharedPtr    m_session;
    SpatialDomains::MeshGraphSharedPtr      m_graph;
    LibUtilities::FieldIOSharedPtr          m_fld;
    map<string, vector<string> >            m_inputfiles;

    bool                                    m_writeBndFld;
    vector<unsigned int>                    m_bndRegionsToWrite;
    bool                                    m_fldToBnd;

    FieldPtsSharedPtr                       m_fieldPts;


    MultiRegions::ExpListSharedPtr SetUpFirstExpList(int NumHomogeneousDir,
                                                     bool fldfilegiven = false)
    {

        MultiRegions::ExpListSharedPtr exp;

        // Set up expansion list
        int expdim  = m_graph->GetMeshDimension();

        bool useFFT     = false;
        bool dealiasing = false;

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

                    if(fldfilegiven)
                    {
                        nplanes = m_fielddef[0]->m_numModes[1];
                        ly      = m_fielddef[0]->m_homogeneousLengths[0];
                        btype   = m_fielddef[0]->m_basis[1];
                    }
                    else
                    {
                        m_session->LoadParameter("HomModesZ", nplanes);
                        m_session->LoadParameter("LY",ly);
                        btype = LibUtilities::eFourier;
                    }

                    // Choose points to be at evenly spaced points at
                    // nplanes points
                    const LibUtilities::PointsKey
                        Pkey(nplanes, LibUtilities::ePolyEvenlySpaced);

                    const LibUtilities::BasisKey Bkey(btype, nplanes, Pkey);



                    if(m_declareExpansionAsContField||
                       m_declareExpansionAsDisContField)
                    {
                        ASSERTL0(false,"ContField2DHomogeneous1D or "
                                 "DisContField2DHomogenenous1D has "
                                 "not been implemented");
                    }

                    Exp2DH1 = MemoryManager<MultiRegions::
                        ExpList2DHomogeneous1D>::
                        AllocateSharedPtr(m_session, Bkey, ly,
                                          useFFT,  dealiasing,
                                          m_graph);
                    exp = Exp2DH1;
                }
                else if (NumHomogeneousDir == 2)
                {
                    MultiRegions::ExpList3DHomogeneous2DSharedPtr Exp3DH2;

                    int nylines,nzlines;
                    NekDouble ly,lz;
                    LibUtilities::BasisType btype1,btype2;

                    if(fldfilegiven)
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
                        m_session->LoadParameter("LY",ly);
                        m_session->LoadParameter("LZ",lz);
                        btype1 = LibUtilities::eFourier;
                        btype2 = LibUtilities::eFourier;
                    }

                    // Choose points to be at evenly spaced points at
                    // nplanes points
                    const LibUtilities::PointsKey
                        PkeyY(nylines, LibUtilities::ePolyEvenlySpaced);
                    const LibUtilities::BasisKey BkeyY(btype1, nylines, PkeyY);

                    const LibUtilities::PointsKey
                        PkeyZ(nzlines, LibUtilities::ePolyEvenlySpaced);
                    const LibUtilities::BasisKey BkeyZ(btype2, nzlines, PkeyZ);

                    if(m_declareExpansionAsContField)
                    {
                        Exp3DH2 = MemoryManager<MultiRegions::
                            ContField3DHomogeneous2D>::
                            AllocateSharedPtr(m_session, BkeyY, BkeyZ,
                                              ly, lz, useFFT, dealiasing,
                                              m_graph,
                                              m_session->GetVariable(0));
                    }
                    else if(m_declareExpansionAsDisContField)
                    {
                        Exp3DH2 = MemoryManager<MultiRegions::
                            DisContField3DHomogeneous2D>::
                            AllocateSharedPtr(m_session, BkeyY, BkeyZ,
                                              ly, lz, useFFT, dealiasing,
                                              m_graph,
                                              m_session->GetVariable(0));
                    }
                    else
                    {
                        Exp3DH2 = MemoryManager<MultiRegions::
                            ExpList3DHomogeneous2D>::
                            AllocateSharedPtr(m_session, BkeyY, BkeyZ,
                                              ly, lz, useFFT, dealiasing,
                                              m_graph);
                    }

                    exp = Exp3DH2;
                }
                else
                {
                    MultiRegions::ExpList1DSharedPtr Exp1D;

                    if(m_declareExpansionAsContField)
                    {
                        Exp1D = MemoryManager<MultiRegions::ContField1D>
                            ::AllocateSharedPtr(m_session, m_graph,
                                                m_session->GetVariable(0));
                    }
                    else if(m_declareExpansionAsDisContField)
                    {
                        Exp1D = MemoryManager<MultiRegions::DisContField1D>
                            ::AllocateSharedPtr(m_session, m_graph,
                                                m_session->GetVariable(0));
                    }
                    else
                    {
                        Exp1D = MemoryManager<MultiRegions::ExpList1D>
                            ::AllocateSharedPtr(m_session, m_graph);
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
                    LibUtilities::BasisType btype;

                    if(fldfilegiven)
                    {
                        nplanes =  m_fielddef[0]->m_numModes[2];
                        lz      = m_fielddef[0]->m_homogeneousLengths[0];
                        btype   = m_fielddef[0]->m_basis[2];
                    }
                    else
                    {
                        m_session->LoadParameter("HomModesZ", nplanes);
                        m_session->LoadParameter("LZ",lz);
                        btype = LibUtilities::eFourier;
                    }

                    // Choose points to be at evenly spaced points at
                    // nplanes points
                    const LibUtilities::PointsKey
                        Pkey(nplanes, LibUtilities::ePolyEvenlySpaced);

                    const LibUtilities::BasisKey  Bkey(btype, nplanes, Pkey);

                    if(m_declareExpansionAsContField)
                    {
                        Exp3DH1 = MemoryManager<MultiRegions::
                            ContField3DHomogeneous1D>::
                            AllocateSharedPtr(m_session, Bkey, lz, useFFT,
                                              dealiasing, m_graph,
                                              m_session->GetVariable(0));
                    }
                    else if (m_declareExpansionAsDisContField)
                    {
                        Exp3DH1 = MemoryManager<MultiRegions::
                            DisContField3DHomogeneous1D>::
                            AllocateSharedPtr(m_session,
                                              Bkey, lz, useFFT,
                                              dealiasing, m_graph,
                                              m_session->GetVariable(0));
                    }
                    else
                    {
                        Exp3DH1 = MemoryManager<MultiRegions::
                            ExpList3DHomogeneous1D>::
                            AllocateSharedPtr(m_session, Bkey, lz, useFFT,
                                              dealiasing, m_graph);
                    }
                    exp = Exp3DH1;
                }
                else
                {
                    MultiRegions::ExpList2DSharedPtr Exp2D;

                    if(m_declareExpansionAsContField)
                    {
                        Exp2D = MemoryManager<MultiRegions::ContField2D>
                            ::AllocateSharedPtr(m_session,m_graph,
                                                m_session->GetVariable(0));
                    }
                    else if(m_declareExpansionAsDisContField)
                    {
                        Exp2D = MemoryManager<MultiRegions::DisContField2D>
                            ::AllocateSharedPtr(m_session,m_graph,
                                                m_session->GetVariable(0));
                    }
                    else
                    {
                        Exp2D = MemoryManager<MultiRegions::ExpList2D>
                            ::AllocateSharedPtr(m_session,m_graph);
                    }

                    exp = Exp2D;
                }
            }
            break;
        case 3:
            {
                MultiRegions::ExpList3DSharedPtr Exp3D;
                
                if(m_declareExpansionAsContField)
                {
                    Exp3D = MemoryManager<MultiRegions::ContField3D>
                        ::AllocateSharedPtr(m_session,m_graph,
                                            m_session->GetVariable(0));
                }
                else if(m_declareExpansionAsDisContField)
                {
                    Exp3D = MemoryManager<MultiRegions::DisContField3D>
                    ::AllocateSharedPtr(m_session,m_graph,
                                        m_session->GetVariable(0));
                }
                else
                {
                    Exp3D = MemoryManager<MultiRegions::ExpList3D>
                        ::AllocateSharedPtr(m_session, m_graph);
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

    MultiRegions::ExpListSharedPtr AppendExpList(int NumHomogeneousDir,
                                                 string var = "DefaultVar",
                                                 bool NewField = false)
    {
        MultiRegions::ExpListSharedPtr tmp;
        switch (m_graph->GetMeshDimension())
        {
        case 1:
            {
                if (NumHomogeneousDir == 1)
                {
                    ASSERTL0(m_declareExpansionAsContField ||
                             m_declareExpansionAsDisContField,
                             "ContField2DHomogeneous1D or "
                             "DisContField2DHomogenenous1D has not been "
                             "implemented");

                    MultiRegions::ExpList2DHomogeneous1DSharedPtr tmp2 =
                        boost::dynamic_pointer_cast<MultiRegions::
                        ExpList2DHomogeneous1D>(m_exp[0]);

                    tmp = MemoryManager<MultiRegions::
                        ExpList2DHomogeneous1D>::
                        AllocateSharedPtr(*tmp2);

                }
                else if (NumHomogeneousDir == 2)
                {
                    if(m_declareExpansionAsContField)
                    {
                        MultiRegions::ContField3DHomogeneous2DSharedPtr tmp2 =
                            boost::dynamic_pointer_cast<MultiRegions::
                            ContField3DHomogeneous2D>(m_exp[0]);

                        tmp = MemoryManager<MultiRegions::
                            ContField3DHomogeneous2D>::
                            AllocateSharedPtr(*tmp2);
                    }
                    else  if(m_declareExpansionAsDisContField)
                    {
                        MultiRegions::DisContField3DHomogeneous2DSharedPtr tmp2 =
                            boost::dynamic_pointer_cast<MultiRegions::
                            DisContField3DHomogeneous2D>(m_exp[0]);

                        tmp = MemoryManager<MultiRegions::
                            DisContField3DHomogeneous2D>::
                            AllocateSharedPtr(*tmp2);
                    }
                    else
                    {
                        MultiRegions::ExpList3DHomogeneous2DSharedPtr tmp2 =
                            boost::dynamic_pointer_cast<MultiRegions::
                            ExpList3DHomogeneous2D>(m_exp[0]);

                        tmp = MemoryManager<MultiRegions::
                            ExpList3DHomogeneous2D>::
                            AllocateSharedPtr(*tmp2);
                    }


                }
                else
                {
                    if(m_declareExpansionAsContField)
                    {
                        MultiRegions::ContField1DSharedPtr tmp2 =
                            boost::dynamic_pointer_cast<MultiRegions::
                            ContField1D>(m_exp[0]);

                        tmp = MemoryManager<MultiRegions::ContField1D>::
                            AllocateSharedPtr(m_session,m_graph,var);
                    }
                    else if(m_declareExpansionAsDisContField)
                    {
                        MultiRegions::DisContField1DSharedPtr tmp2 =
                            boost::dynamic_pointer_cast<MultiRegions::
                            DisContField1D>(m_exp[0]);
                        
                        tmp = MemoryManager<MultiRegions::DisContField1D>::
                            AllocateSharedPtr(m_session,m_graph,var);
                    }
                    else
                    {
                        MultiRegions::ExpList1DSharedPtr tmp2 =
                            boost::dynamic_pointer_cast<MultiRegions::
                            ExpList1D>(m_exp[0]);
                        
                        tmp = MemoryManager<MultiRegions::ExpList1D>::
                            AllocateSharedPtr(*tmp2);
                    }

                }
            }
            break;
        case 2:
            {
                if (NumHomogeneousDir == 1)
                {
                    if(m_declareExpansionAsContField)
                    {
                        if(NewField)
                        {
                            bool useFFT     = false;
                            bool dealiasing = false;

                            tmp  = MemoryManager<MultiRegions::
                                ContField3DHomogeneous1D>::AllocateSharedPtr(
                                        m_session,
                                        m_exp[0]->GetHomogeneousBasis()
                                                    ->GetBasisKey(),
                                        m_exp[0]->GetHomoLen(),
                                        useFFT, dealiasing, m_graph, var);
                        }
                        else
                        {
                            MultiRegions::ContField3DHomogeneous1DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                ContField3DHomogeneous1D>(m_exp[0]);
                            
                            ASSERTL0(tmp2,"Failed to type cast m_exp[0]");
                            tmp = MemoryManager<MultiRegions::
                                ContField3DHomogeneous1D>::
                                AllocateSharedPtr(*tmp2);
                        }
                    }
                    else  if(m_declareExpansionAsDisContField)
                    {
                        if(NewField)
                        {
                            bool useFFT     = false;
                            bool dealiasing = false;

                            tmp  = MemoryManager<MultiRegions::
                                DisContField3DHomogeneous1D>::AllocateSharedPtr(
                                        m_session,
                                        m_exp[0]->GetHomogeneousBasis()
                                                    ->GetBasisKey(),
                                        m_exp[0]->GetHomoLen(),
                                        useFFT, dealiasing, m_graph,var);
                        }
                        else
                        {
                            MultiRegions::DisContField3DHomogeneous1DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                DisContField3DHomogeneous1D>(m_exp[0]);
                            ASSERTL0(tmp2,"Failed to type cast m_exp[0]");
                            
                            tmp = MemoryManager<MultiRegions::
                                DisContField3DHomogeneous1D>::
                                AllocateSharedPtr(*tmp2);
                        }
                    }
                    else
                    {
                        if(NewField)
                        {
                            bool useFFT     = false;
                            bool dealiasing = false;
                            
                            tmp  = MemoryManager<MultiRegions::
                                ExpList3DHomogeneous1D>::AllocateSharedPtr(
                                        m_session,
                                        m_exp[0]->GetHomogeneousBasis()
                                                    ->GetBasisKey(),
                                        m_exp[0]->GetHomoLen(),
                                        useFFT, dealiasing, m_graph);
                        }
                        else
                        {
                            MultiRegions::ExpList3DHomogeneous1DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                ExpList3DHomogeneous1D>(m_exp[0]);
                            ASSERTL0(tmp2,"Failed to type cast m_exp[0]");

                            tmp = MemoryManager<MultiRegions::
                                ExpList3DHomogeneous1D>::
                                AllocateSharedPtr(*tmp2);
                        }
                    }

                }
                else
                {
                    if(m_declareExpansionAsContField)
                    {
                        if(NewField)
                        {
                            tmp = MemoryManager<MultiRegions::ContField2D>::
                                AllocateSharedPtr(m_session,m_graph,var);
                        }
                        else // call copy constructor
                        {

                            MultiRegions::ContField2DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                ContField2D>(m_exp[0]);

                            tmp = MemoryManager<MultiRegions::ContField2D>::
                                AllocateSharedPtr(*tmp2,m_graph,var);
                        }
                    }
                    else if(m_declareExpansionAsDisContField)
                    {
                        if(NewField)
                        {
                            tmp = MemoryManager<MultiRegions::DisContField2D>::
                                AllocateSharedPtr(m_session,m_graph,var);
                        }
                        else // call copy constructor
                        {
                            MultiRegions::DisContField2DSharedPtr tmp2 =
                                boost::dynamic_pointer_cast<MultiRegions::
                                DisContField2D>(m_exp[0]);
                            
                            tmp = MemoryManager<MultiRegions::DisContField2D>::
                                AllocateSharedPtr(*tmp2,m_graph,var);
                        }
                    }
                    else
                    {
                        MultiRegions::ExpList2DSharedPtr tmp2 =
                            boost::dynamic_pointer_cast<MultiRegions::
                            ExpList2D>(m_exp[0]);

                        tmp = MemoryManager<MultiRegions::ExpList2D>::
                            AllocateSharedPtr(*tmp2);
                    }
                }
            }
            break;
        case 3:
            {
                if(m_declareExpansionAsContField)
                {
                    if(NewField)
                    {
                        tmp = MemoryManager<MultiRegions::ContField3D>::
                            AllocateSharedPtr(m_session,m_graph,var);
                    }
                    else
                    {
                        MultiRegions::ContField3DSharedPtr tmp2 =
                            boost::dynamic_pointer_cast<MultiRegions::
                            ContField3D>(m_exp[0]);

                        tmp = MemoryManager<MultiRegions::ContField3D>::
                            AllocateSharedPtr(*tmp2,m_graph,var);
                    }
                }
                else if(m_declareExpansionAsDisContField)
                {
                    if(NewField)
                    {
                        tmp = MemoryManager<MultiRegions::DisContField3D>::
                            AllocateSharedPtr(m_session,m_graph,var);
                    }
                    else
                    {
                        MultiRegions::DisContField3DSharedPtr tmp2 =
                            boost::dynamic_pointer_cast<MultiRegions::
                            DisContField3D>(m_exp[0]);

                        tmp = MemoryManager<MultiRegions::DisContField3D>::
                            AllocateSharedPtr(*tmp2,m_graph,var);
                    }
                }
                else
                {
                    MultiRegions::ExpList3DSharedPtr tmp2 =
                        boost::dynamic_pointer_cast<MultiRegions::
                        ExpList3D>(m_exp[0]);

                    tmp = MemoryManager<MultiRegions::ExpList3D>::
                        AllocateSharedPtr(*tmp2);
                }
            }
            break;
        default:
            ASSERTL0(false, "Expansion dimension not recognised");
            break;
        }

        return tmp;
    }

};

typedef boost::shared_ptr<Field> FieldSharedPtr;

}
}

