///////////////////////////////////////////////////////////////////////////////
//
// File LinearisedAdvection.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Evaluation of the linearised advective term
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/AdvectionTerms/LinearisedAdvection.h>
#include <StdRegions/StdSegExp.h>

using namespace std;

namespace Nektar
{
string LinearisedAdvection::className
        = SolverUtils::GetAdvectionFactory().RegisterCreatorFunction(
                "Linearised",
                LinearisedAdvection::create);

/**
 * Constructor. Creates ...
 *
 * \param
 * \param
 */

LinearisedAdvection::LinearisedAdvection():
    Advection()
{
}


void LinearisedAdvection::v_InitObject(
        LibUtilities::SessionReaderSharedPtr        pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
{
    Advection::v_InitObject(pSession, pFields);

    m_session            = pSession;
    m_boundaryConditions = MemoryManager<SpatialDomains::BoundaryConditions>
        ::AllocateSharedPtr(m_session, pFields[0]->GetGraph());
    m_spacedim           = pFields[0]->GetGraph()->GetSpaceDimension();
    m_expdim             = pFields[0]->GetGraph()->GetMeshDimension();

    //Setting parameters for homogeneous problems
    m_HomoDirec          = 0;
    m_useFFT             = false;
    m_HomogeneousType    = eNotHomogeneous;
    m_singleMode         = false;
    m_halfMode           = false;
    m_multipleModes      = false;

    if(m_session->DefinesSolverInfo("HOMOGENEOUS"))
    {
        std::string HomoStr = m_session->GetSolverInfo("HOMOGENEOUS");
        m_spacedim          = 3;

        if((HomoStr == "HOMOGENEOUS1D")||(HomoStr == "Homogeneous1D")||
           (HomoStr == "1D")||(HomoStr == "Homo1D"))
        {
            m_HomogeneousType = eHomogeneous1D;
            m_LhomZ           = m_session->GetParameter("LZ");
            m_HomoDirec       = 1;

            ASSERTL0(m_session->DefinesSolverInfo("ModeType"),
                     "Need to specify ModeType as HalfMode,SingleMode or "
                     "MultipleModes");

            m_session->MatchSolverInfo("ModeType",      "SingleMode",
                                       m_singleMode,    false);
            m_session->MatchSolverInfo("ModeType",      "HalfMode",
                                       m_halfMode,      false);
            m_session->MatchSolverInfo("ModeType",      "MultipleModes",
                                       m_multipleModes, false);

            if(m_singleMode)
            {
                m_npointsZ = 2;
            }
            else if(m_halfMode)
            {
                m_npointsZ = 1;
            }
            else if(m_multipleModes)
            {
                m_npointsZ = m_session->GetParameter("HomModesZ");
            }
        }

        if((HomoStr == "HOMOGENEOUS2D")||(HomoStr == "Homogeneous2D")||
           (HomoStr == "2D")||(HomoStr == "Homo2D"))
        {
            m_HomogeneousType = eHomogeneous2D;
            m_session->LoadParameter("HomModesY", m_npointsY);
            m_session->LoadParameter("LY",        m_LhomY);
            m_session->LoadParameter("HomModesZ", m_npointsZ);
            m_session->LoadParameter("LZ",        m_LhomZ);
            m_HomoDirec       = 2;
        }

        if((HomoStr == "HOMOGENEOUS3D")||(HomoStr == "Homogeneous3D")||
           (HomoStr == "3D")||(HomoStr == "Homo3D"))
        {
            m_HomogeneousType = eHomogeneous3D;
            m_session->LoadParameter("HomModesX",m_npointsX);
            m_session->LoadParameter("LX",       m_LhomX   );
            m_session->LoadParameter("HomModesY",m_npointsY);
            m_session->LoadParameter("LY",       m_LhomY   );
            m_session->LoadParameter("HomModesZ",m_npointsZ);
            m_session->LoadParameter("LZ",       m_LhomZ   );
            m_HomoDirec       = 3;
        }

        if(m_session->DefinesSolverInfo("USEFFT"))
        {
            m_useFFT = true;
        }
    }
    else
    {
        m_npointsZ = 1; // set to default value so can use to identify 2d or 3D (homogeneous) expansions
    }

    int nvar = m_session->GetVariables().size();
    m_baseflow = Array<OneD, Array<OneD, NekDouble> >(nvar);
    for (int i = 0; i < nvar; ++i)
    {
        m_baseflow[i] = Array<OneD, NekDouble>(pFields[i]->GetTotPoints(), 0.0);
    }

    int nBaseDerivs = (m_halfMode || m_singleMode) ? 2 : m_spacedim;
    m_gradBase = Array<OneD, Array<OneD, NekDouble> >(nvar*nBaseDerivs);
    for (int i = 0; i < nvar; ++i)
    {
        for (int j = 0; j < nBaseDerivs; ++j)
        {
            m_gradBase[i*nBaseDerivs + j ] = Array<OneD, NekDouble>
                                            (pFields[i]->GetTotPoints(), 0.0);
        }
    }

    ASSERTL0(m_session->DefinesFunction("BaseFlow"),
             "Base flow must be defined for linearised forms.");
    string file = m_session->GetFunctionFilename("BaseFlow", 0);


    //Periodic base flows
    if (m_session->DefinesParameter("N_slices"))
    {
        m_session->LoadParameter("N_slices",m_slices);
        if (m_slices > 1)
        {
            ASSERTL0(m_session->GetFunctionType("BaseFlow", 0)
                == LibUtilities::eFunctionTypeFile,
                "Base flow should be a sequence of files.");
            m_session->LoadParameter("BaseFlow_interporder", m_interporder, 0);
            ASSERTL0(m_interporder <= m_slices,
                "BaseFlow_interporder should be smaller than or equal to N_slices.");
            m_isperiodic = m_interporder < 2;
            m_session->LoadParameter("N_start", m_start, 0);
            m_session->LoadParameter("N_skip", m_skip, 1);
            DFT(file,pFields,m_slices);
        }
        else
        {
            ASSERTL0(false, "Number of slices must be a positive number "
                            "greater than 1");
        }
    }
    //Steady base-flow
    else
    {
        m_slices=1;

        //BaseFlow from file
        if (m_session->GetFunctionType("BaseFlow", m_session->GetVariable(0))
            == LibUtilities::eFunctionTypeFile)
        {
            ImportFldBase(file,pFields,1);

        }
        //analytic base flow
        else
        {
            int nq = pFields[0]->GetNpoints();
            Array<OneD,NekDouble> x0(nq);
            Array<OneD,NekDouble> x1(nq);
            Array<OneD,NekDouble> x2(nq);

            // get the coordinates (assuming all fields have the same
            // discretisation)
            pFields[0]->GetCoords(x0,x1,x2);

            for(unsigned int i = 0 ; i < pFields.size(); i++)
            {
                LibUtilities::EquationSharedPtr ifunc
                    = m_session->GetFunction("BaseFlow", i);

                ifunc->Evaluate(x0,x1,x2,m_baseflow[i]);
            }
        }
    }

    for (int i = 0; i < nvar; ++i)
    {
        UpdateGradBase(i, pFields[i]);
    }

    if (m_session->DefinesParameter("period"))
    {
        m_period=m_session->GetParameter("period");
    }
    else
    {
        m_period=(m_session->GetParameter("TimeStep")*m_slices)/(m_slices-1.);
    }
    if (m_session->GetComm()->GetRank() == 0)
    {
        cout << "baseflow info : interpolation order " << m_interporder
             << ", period " << m_period << ", periodicity ";
        if (m_isperiodic)
        {
            cout << "yes\n";
        }
        else
        {
            cout << "no\n";
        }
        cout << "baseflow info : files from " << m_start << " to "
             << (m_start + (m_slices-1)*m_skip)
             << " (skip " << m_skip  << ") with " << (m_slices -  (m_interporder > 1))
             << " time intervals" << endl;
    }
}

LinearisedAdvection::~LinearisedAdvection()
{
}


//Advection function

void LinearisedAdvection::v_Advect(
        const int nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &advVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        Array<OneD, Array<OneD, NekDouble> >              &outarray,
        const NekDouble                                   &time,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
{
    ASSERTL1(nConvectiveFields == inarray.size(),
             "Number of convective fields and Inarray are not compatible");

    int nPointsTot  = fields[0]->GetNpoints();
    int ndim        = advVel.size();
    int nBaseDerivs = (m_halfMode || m_singleMode) ? 2 : m_spacedim;
    int nDerivs     = (m_halfMode) ? 2 : m_spacedim;

    Array<OneD, Array<OneD, NekDouble> > velocity(ndim);
    for(int i = 0; i < ndim; ++i)
    {
        if(fields[i]->GetWaveSpace() && !m_singleMode && !m_halfMode)
        {
            velocity[i] = Array<OneD, NekDouble>(nPointsTot,0.0);
            fields[i]->HomogeneousBwdTrans(advVel[i],velocity[i]);
        }
        else
        {
            velocity[i] = advVel[i];
        }
    }

    Array<OneD, Array<OneD, NekDouble> > grad (nDerivs);
    for( int i = 0; i < nDerivs; ++i)
    {
        grad[i] = Array<OneD, NekDouble> (nPointsTot);
    }

    // Evaluation of the base flow for periodic cases
    if (m_slices > 1)
    {
        for (int i = 0; i < ndim; ++i)
        {
            UpdateBase(m_slices, m_interp[i], m_baseflow[i],
                       time, m_period);
            UpdateGradBase(i, fields[i]);
        }
    }

    //Evaluate the linearised advection term
    for( int i = 0; i < nConvectiveFields; ++i)
    {
        // Calculate gradient
        switch(nDerivs)
        {
            case 1:
            {
                fields[i]->PhysDeriv(inarray[i], grad[0]);
            }
            break;
            case 2:
            {
                fields[i]->PhysDeriv(inarray[i], grad[0], grad[1]);
            }
            break;
            case 3:
            {
                fields[i]->PhysDeriv(inarray[i], grad[0], grad[1], grad[2]);
                if(m_multipleModes)
                {
                    // transform gradients into physical Fourier space
                    fields[i]->HomogeneousBwdTrans(grad[0], grad[0]);
                    fields[i]->HomogeneousBwdTrans(grad[1], grad[1]);
                    fields[i]->HomogeneousBwdTrans(grad[2], grad[2]);
                }
            }
            break;
        }

        // Calculate U_j du'_i/dx_j
        Vmath::Vmul(nPointsTot,grad[0], 1, m_baseflow[0], 1, outarray[i], 1);
        for( int j = 1; j < nDerivs; ++j)
        {
            Vmath::Vvtvp(nPointsTot,grad[j], 1,
                                    m_baseflow[j], 1,
                                    outarray[i], 1,
                                    outarray[i], 1);
        }

        // Add u'_j dU_i/dx_j
        int lim = (m_halfMode || m_singleMode) ? 2 : ndim;
        if (m_halfMode && i==2)
        {
            lim = 0;
        }
        for( int j = 0; j < lim; ++j)
        {
            Vmath::Vvtvp(nPointsTot,m_gradBase[i*nBaseDerivs + j], 1,
                                    velocity[j], 1,
                                    outarray[i], 1,
                                    outarray[i], 1);
        }

        if(m_multipleModes)
        {
            fields[i]->HomogeneousFwdTrans(outarray[i],outarray[i]);
        }
        Vmath::Neg(nPointsTot,outarray[i],1);
    }
}

void LinearisedAdvection::v_SetBaseFlow(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
{
    if (m_session->GetSolverInfo("EqType") == "UnsteadyNavierStokes")
    {
        // The SFD method is only applied to the velocity variables in
        // incompressible
        ASSERTL1(inarray.size() == (m_baseflow.size() - 1),
                 "Number of base flow variables does not match what is "
                 "expected.");
    }
    else
    {
        ASSERTL1(inarray.size() == (m_baseflow.size()),
             "Number of base flow variables does not match what is expected.");
    }

    int npts = inarray[0].size();

    for (int i = 0; i < inarray.size(); ++i)
    {
        ASSERTL1(npts == m_baseflow[i].size(),
             "Size of base flow array does not match expected.");
        Vmath::Vcopy(npts, inarray[i], 1, m_baseflow[i], 1);
        UpdateGradBase(i, fields[i]);
    }
}


/**
 * Import field from infile and load into \a m_fields. This routine will
 * also perform a \a BwdTrans to ensure data is in both the physical and
 * coefficient storage.
 * @param   pInFile          Filename to read.
 * @param   pFields          Array of expansion lists
 */
void LinearisedAdvection::ImportFldBase(
        std::string                                  pInfile,
        Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        int                                          pSlice)
{
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
    std::vector<std::vector<NekDouble> >                 FieldData;

    int nqtot = m_baseflow[0].size();
    Array<OneD, NekDouble> tmp_coeff(pFields[0]->GetNcoeffs(), 0.0);

    int numexp = pFields[0]->GetExpSize();
    Array<OneD,int> ElementGIDs(numexp);

    // Define list of global element ids
    for(int i = 0; i < numexp; ++i)
    {
        ElementGIDs[i] = pFields[0]->GetExp(i)->GetGeom()->GetGlobalID();
    }

    //Get Homogeneous
    LibUtilities::FieldIOSharedPtr fld = LibUtilities::FieldIO::CreateForFile(
        m_session, pInfile);
    fld->Import(pInfile, FieldDef, FieldData,
                LibUtilities::NullFieldMetaDataMap,
                ElementGIDs);

    int nSessionVar     = m_session->GetVariables().size();
    int nSessionConvVar = nSessionVar - 1;
    int nFileVar        = FieldDef[0]->m_fields.size();
    int nFileConvVar    = nFileVar - 1; // Ignore pressure
    if (m_halfMode)
    {
        ASSERTL0(nFileVar == 3, "For half mode, expect 2D2C base flow.");
        nFileConvVar = 2;
    }

    for(int j = 0; j < nFileConvVar; ++j)
    {
        for(int i = 0; i < FieldDef.size(); ++i)
        {
            bool flag = FieldDef[i]->m_fields[j] ==
                m_session->GetVariable(j);

            ASSERTL0(flag, (std::string("Order of ") + pInfile
                            + std::string(" data and that defined in "
                            "the session file differs")).c_str());

            pFields[j]->ExtractDataToCoeffs(
                                FieldDef[i],
                                FieldData[i],
                                FieldDef[i]->m_fields[j],
                                tmp_coeff);
        }

        if(m_singleMode || m_halfMode)
        {
            pFields[j]->GetPlane(0)->BwdTrans(tmp_coeff, m_baseflow[j]);

            if(m_singleMode)
            {
                //copy the bwd trans into the second plane for single
                //Mode Analysis
                int ncplane=(pFields[0]->GetNpoints())/m_npointsZ;
                Vmath::Vcopy(ncplane,&m_baseflow[j][0],1,&m_baseflow[j][ncplane],1);
            }
        }
        else // fully 3D base flow - put in physical space.
        {
            bool oldwavespace = pFields[j]->GetWaveSpace();
            pFields[j]->SetWaveSpace(false);
            pFields[j]->BwdTrans(tmp_coeff, m_baseflow[j]);
            pFields[j]->SetWaveSpace(oldwavespace);
        }
    }

    // Zero unused fields (e.g. w in a 2D2C base flow).
    for (int j = nFileConvVar; j < nSessionConvVar; ++j)
    {
        Vmath::Fill(nqtot, 0.0, m_baseflow[j], 1);
    }

    // If time-periodic, put loaded data into the slice storage.
    if (m_slices > 1)
    {
        for(int i = 0; i < nSessionConvVar; ++i)
        {
            Vmath::Vcopy(nqtot, &m_baseflow[i][0], 1, &m_interp[i][pSlice*nqtot], 1);
        }
    }
}


void LinearisedAdvection::UpdateBase(
        const NekDouble                     m_slices,
        const Array<OneD, const NekDouble> &inarray,
        Array<OneD, NekDouble>             &outarray,
        const NekDouble                     m_time,
        const NekDouble                     m_period)
{
    int npoints     = m_baseflow[0].size();
    if (m_isperiodic)
    {
        NekDouble BetaT = 2*M_PI*fmod (m_time, m_period) / m_period;
        NekDouble phase;
        Array<OneD, NekDouble> auxiliary(npoints);

        Vmath::Vcopy(npoints,&inarray[0],1,&outarray[0],1);
        Vmath::Svtvp(npoints, cos(0.5*m_slices*BetaT),&inarray[npoints],1,&outarray[0],1,&outarray[0],1);

        for (int i = 2; i < m_slices; i += 2)
        {
            phase = (i>>1) * BetaT;

            Vmath::Svtvp(npoints, cos(phase),&inarray[i*npoints],1,&outarray[0],1,&outarray[0],1);
            Vmath::Svtvp(npoints, -sin(phase), &inarray[(i+1)*npoints], 1, &outarray[0], 1,&outarray[0],1);
        }
    }
    else
    {
        NekDouble x = m_time;
        x = x/m_period*(m_slices-1);
        int ix = x;
        if (ix < 0)
        {
            ix = 0;
        }
        if (ix > m_slices-2)
        {
            ix = m_slices-2;
        }
        int padleft = m_interporder/2 - 1;
        if (padleft > ix)
        {
            padleft = ix;
        }
        int padright = m_interporder - 1 - padleft;
        if (padright > m_slices-1-ix)
        {
            padright = m_slices-1-ix;
        }
        padleft = m_interporder - 1 - padright;
        Array<OneD, NekDouble> coeff(m_interporder, 1.);
        for (int i=0; i<m_interporder; ++i)
        {
            for (int j=0; j<m_interporder; ++j)
            {
                if (i!=j)
                {
                    coeff[i] *= (x - ix + padleft - (NekDouble)j) /
                                ((NekDouble)i - (NekDouble)j);
                }
            }
        }
        Vmath::Zero(npoints, &outarray[0], 1);
        for (int i = ix-padleft; i < ix+padright+1; ++i)
        {
            Vmath::Svtvp(npoints, coeff[i-ix+padleft], &inarray[i*npoints], 1,
                        &outarray[0], 1, &outarray[0], 1);
        }
    }
}

void LinearisedAdvection::UpdateGradBase(
        const int                                var,
        const MultiRegions::ExpListSharedPtr     &field)
{
    int nBaseDerivs = (m_halfMode || m_singleMode) ? 2 : m_spacedim;
    switch(m_spacedim)
    {
        case 1:         // 1D
        {
            field->PhysDeriv(m_baseflow[var], m_gradBase[var*nBaseDerivs + 0]);
        }
        break;
        case 2:  //2D
        {
            field->PhysDeriv(m_baseflow[var], m_gradBase[var*nBaseDerivs + 0],
                                              m_gradBase[var*nBaseDerivs + 1]);
        }
        break;
        case 3:
        {
            if(m_halfMode) // can assume W = 0 and  d/dz == 0
            {
                if( var < 2)
                {
                    field->PhysDeriv(m_baseflow[var],
                                              m_gradBase[var*nBaseDerivs + 0],
                                              m_gradBase[var*nBaseDerivs + 1]);
                }
            }
            else if(m_singleMode) // single mode where d/dz = 0
            {
                field->PhysDeriv(m_baseflow[var], m_gradBase[var*nBaseDerivs + 0],
                                                  m_gradBase[var*nBaseDerivs + 1]);
            }
            else
            {
                // Differentiate base flow in physical space
                bool oldwavespace = field->GetWaveSpace();
                field->SetWaveSpace(false);
                field->PhysDeriv(m_baseflow[var], m_gradBase[var*nBaseDerivs + 0],
                                                  m_gradBase[var*nBaseDerivs + 1],
                                                  m_gradBase[var*nBaseDerivs + 2]);
                field->SetWaveSpace(oldwavespace);
            }
        }
        break;
    }
}

DNekBlkMatSharedPtr LinearisedAdvection::GetFloquetBlockMatrix(FloquetMatType mattype, bool UseContCoeffs) const
{
    DNekMatSharedPtr    loc_mat;
    DNekBlkMatSharedPtr BlkMatrix;
    int n_exp = 0;

    n_exp = m_baseflow[0].size(); // will operatore on m_phys

    Array<OneD,unsigned int> nrows(n_exp);
    Array<OneD,unsigned int> ncols(n_exp);

    nrows = Array<OneD, unsigned int>(n_exp,m_slices);
    ncols = Array<OneD, unsigned int>(n_exp,m_slices);

    MatrixStorage blkmatStorage = eDIAGONAL;
    BlkMatrix = MemoryManager<DNekBlkMat>
        ::AllocateSharedPtr(nrows,ncols,blkmatStorage);


    const LibUtilities::PointsKey Pkey(m_slices,LibUtilities::eFourierEvenlySpaced);
    const LibUtilities::BasisKey  BK(LibUtilities::eFourier,m_slices,Pkey);
    StdRegions::StdSegExp StdSeg(BK);

    StdRegions::StdMatrixKey matkey(StdRegions::eFwdTrans,
                                    StdSeg.DetShapeType(),
                                    StdSeg);

    loc_mat = StdSeg.GetStdMatrix(matkey);

    // set up array of block matrices.
    for(int i = 0; i < n_exp; ++i)
    {
        BlkMatrix->SetBlock(i,i,loc_mat);
    }

    return BlkMatrix;
}

//Discrete Fourier Transform for Floquet analysis
void LinearisedAdvection::DFT(const string file,
        Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
        const NekDouble m_slices)
{
    int ConvectedFields = m_baseflow.size()-1;
    int npoints         = m_baseflow[0].size();
    m_interp            = Array<OneD, Array<OneD, NekDouble> > (ConvectedFields);

    for (int i = 0; i < ConvectedFields; ++i)
    {
        m_interp[i] = Array<OneD,NekDouble>(npoints*m_slices, 0.0);
    }

    // Import the slides into auxiliary vector
    // The base flow should be stored in the form "filename_%d.ext"
    // A subdirectory can also be included, such as "dir/filename_%d.ext"
    size_t found = file.find("%d");
    ASSERTL0(found != string::npos && file.find("%d", found+1) == string::npos,
             "Since N_slices is specified, the filename provided for function "
             "'BaseFlow' must include exactly one instance of the format "
             "specifier '%d', to index the time-slices.");
    char* buffer = new char[file.length() + 8];
    int nstart = m_start;
    for (int i = nstart; i < nstart+m_slices*m_skip; i+=m_skip)
    {
        sprintf(buffer, file.c_str(), i);
        ImportFldBase(buffer,pFields,(i-nstart)/m_skip);
        if(m_session->GetComm()->GetRank() == 0)
        {
            cout << "read base flow file " << buffer << endl;
        }
    }
    delete[] buffer;
    if(!m_isperiodic)
    {
        return;
    }

    // Discrete Fourier Transform of the fields
    for(int k=0; k< ConvectedFields;++k)
    {
#ifdef NEKTAR_USING_FFTW

        //Discrete Fourier Transform using FFTW
        Array<OneD, NekDouble> fft_in(npoints*m_slices);
        Array<OneD, NekDouble> fft_out(npoints*m_slices);

        Array<OneD, NekDouble> m_tmpIN(m_slices);
        Array<OneD, NekDouble> m_tmpOUT(m_slices);

        //Shuffle the data
        for(int j= 0; j < m_slices; ++j)
        {
            Vmath::Vcopy(npoints,&m_interp[k][j*npoints],1,&(fft_in[j]),m_slices);
        }

        m_FFT = LibUtilities::GetNektarFFTFactory().CreateInstance("NekFFTW", m_slices);

        //FFT Transform
        for(int i=0; i<npoints; i++)
        {
            m_FFT->FFTFwdTrans(m_tmpIN =fft_in + i*m_slices, m_tmpOUT =fft_out + i*m_slices);

        }

        //Reshuffle data
        for(int s = 0; s < m_slices; ++s)
        {
            Vmath::Vcopy(npoints,&fft_out[s],m_slices,&m_interp[k][s*npoints],1);

        }

        Vmath::Zero(fft_in.size(),&fft_in[0],1);
        Vmath::Zero(fft_out.size(),&fft_out[0],1);
#else
        //Discrete Fourier Transform using MVM
        DNekBlkMatSharedPtr blkmat;
        blkmat = GetFloquetBlockMatrix(eForwardsPhys);

        int nrows = blkmat->GetRows();
        int ncols = blkmat->GetColumns();

        Array<OneD, NekDouble> sortedinarray(ncols);
        Array<OneD, NekDouble> sortedoutarray(nrows);

        //Shuffle the data
        for(int j= 0; j < m_slices; ++j)
        {
            Vmath::Vcopy(npoints,&m_interp[k][j*npoints],1,&(sortedinarray[j]),m_slices);
        }

        // Create NekVectors from the given data arrays
        NekVector<NekDouble> in (ncols,sortedinarray,eWrapper);
        NekVector<NekDouble> out(nrows,sortedoutarray,eWrapper);

        // Perform matrix-vector multiply.
        out = (*blkmat)*in;

        //Reshuffle data
        for(int s = 0; s < m_slices; ++s)
        {
            Vmath::Vcopy(npoints,&sortedoutarray[s],m_slices,&m_interp[k][s*npoints],1);
        }

        for(int r=0; r<sortedinarray.size();++r)
        {
            sortedinarray[0]=0;
            sortedoutarray[0]=0;
        }

#endif
    }

}

} //end of namespace

