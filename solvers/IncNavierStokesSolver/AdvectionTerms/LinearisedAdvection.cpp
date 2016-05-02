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
// License for the specific language governing rights and limitations under
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

#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/DisContField2D.h>
#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <MultiRegions/DisContField3DHomogeneous2D.h>

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

    ASSERTL0(m_session->DefinesFunction("BaseFlow"),
             "Base flow must be defined for linearised forms.");
    string file = m_session->GetFunctionFilename("BaseFlow", 0);


    //Periodic base flows
    if(m_session->DefinesParameter("N_slices"))
    {
        m_session->LoadParameter("N_slices",m_slices);
        if(m_slices>1)
        {
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

            for(unsigned int i = 0 ; i < pFields.num_elements(); i++)
            {
                LibUtilities::EquationSharedPtr ifunc
                    = m_session->GetFunction("BaseFlow", i);

                ifunc->Evaluate(x0,x1,x2,m_baseflow[i]);
            }
        }
    }

    if(m_session->DefinesParameter("period"))
    {
        m_period=m_session->GetParameter("period");
    }
    else
    {
        m_period=(m_session->GetParameter("TimeStep")*m_slices)/(m_slices-1.);
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
        const NekDouble                                   &time)
{
    ASSERTL1(nConvectiveFields == inarray.num_elements(),
             "Number of convective fields and Inarray are not compatible");

    int nPointsTot = fields[0]->GetNpoints();
    int ndim       = advVel.num_elements();

    Array<OneD, NekDouble > Deriv
                = Array<OneD, NekDouble> (nPointsTot*nConvectiveFields);

    Array<OneD, NekDouble> grad0,grad1,grad2;

    // Evaluation of the gradient of each component of the base flow
    // \nabla U
    Array<OneD, NekDouble> grad_base_u0,grad_base_u1,grad_base_u2;

    // \nabla V
    Array<OneD, NekDouble> grad_base_v0,grad_base_v1,grad_base_v2;

    // \nabla W
    Array<OneD, NekDouble> grad_base_w0,grad_base_w1,grad_base_w2;

    grad0        = Array<OneD, NekDouble> (nPointsTot);
    grad_base_u0 = Array<OneD, NekDouble> (nPointsTot);
    grad_base_v0 = Array<OneD, NekDouble> (nPointsTot);
    grad_base_w0 = Array<OneD, NekDouble> (nPointsTot);

    // Evaluation of the base flow for periodic cases
    if (m_slices > 1)
    {
        ASSERTL0(m_session->GetFunctionType("BaseFlow", 0)
                 == LibUtilities::eFunctionTypeFile,
                 "Base flow should be a sequence of files.");

        for (int i = 0; i < ndim; ++i)
        {
            UpdateBase(m_slices, m_interp[i], m_baseflow[i],
                       time, m_period);
        }
    }


    //Evaluate the linearised advection term
    switch(ndim)
    {
        case 1:         // 1D
            ASSERTL0(false,"Not set up for 1D");
            break;
        case 2:  //2D
        {
            grad1 = Array<OneD, NekDouble> (nPointsTot);
            grad_base_u1 = Array<OneD, NekDouble> (nPointsTot);
            grad_base_v1 = Array<OneD, NekDouble> (nPointsTot);

            //Derivates of the base flow
            fields[0]-> PhysDeriv(m_baseflow[0], grad_base_u0, grad_base_u1);
            fields[0]-> PhysDeriv(m_baseflow[1], grad_base_v0, grad_base_v1);

            //x-equation
            fields[0]->PhysDeriv(inarray[0],grad0,grad1);
            // Evaluate U du'/dx
            Vmath::Vmul (nPointsTot,grad0,1,m_baseflow[0],1,outarray[0],1);
            //Evaluate U du'/dx+ V du'/dy
            Vmath::Vvtvp(nPointsTot,grad1,1,m_baseflow[1],1,outarray[0],1,
                         outarray[0],1);
            //Evaluate (U du'/dx+ V du'/dy)+u' dU/dx
            Vmath::Vvtvp(nPointsTot,grad_base_u0,1,advVel[0],1,outarray[0],1,
                         outarray[0],1);
            //Evaluate (U du'/dx+ V du'/dy +u' dU/dx)+v' dU/dy
            Vmath::Vvtvp(nPointsTot,grad_base_u1,1,advVel[1],1,outarray[0],1,
                         outarray[0],1);
            Vmath::Neg(nPointsTot,outarray[0],1);

            //y-equation
            fields[0]->PhysDeriv(inarray[1],grad0,grad1);
            // Evaluate U dv'/dx
            Vmath::Vmul (nPointsTot,grad0,1,m_baseflow[0],1,outarray[1],1);
            //Evaluate U dv'/dx+ V dv'/dy
            Vmath::Vvtvp(nPointsTot,grad1,1,m_baseflow[1],1,outarray[1],1,
                         outarray[1],1);
            //Evaluate (U dv'/dx+ V dv'/dy)+u' dV/dx
            Vmath::Vvtvp(nPointsTot,grad_base_v0,1,advVel[0],1,outarray[1],1,
                         outarray[1],1);
            //Evaluate (U dv'/dx+ V dv'/dy +u' dv/dx)+v' dV/dy
            Vmath::Vvtvp(nPointsTot,grad_base_v1,1,advVel[1],1,outarray[1],1,
                         outarray[1],1);
            Vmath::Neg(nPointsTot,outarray[1],1);
        }
        break;
        case 3:  //3D
        {
            grad1 = Array<OneD, NekDouble> (nPointsTot);
            grad2 = Array<OneD, NekDouble> (nPointsTot);

            grad_base_u1 = Array<OneD, NekDouble> (nPointsTot);
            grad_base_v1 = Array<OneD, NekDouble> (nPointsTot);
            grad_base_w1 = Array<OneD, NekDouble> (nPointsTot,0.0);

            grad_base_u2 = Array<OneD, NekDouble> (nPointsTot,0.0);
            grad_base_v2 = Array<OneD, NekDouble> (nPointsTot,0.0);
            grad_base_w2 = Array<OneD, NekDouble> (nPointsTot,0.0);

            // note base flow should be moved to GetBaseFlow method
            if(m_halfMode) // can assume W = 0 and  d/dz == 0
            {
                // note base flow should be moved to GetBaseFlow method
                fields[0]->PhysDeriv(m_baseflow[0],
                                     grad_base_u0, grad_base_u1);
                fields[0]->PhysDeriv(m_baseflow[1],
                                     grad_base_v0, grad_base_v1);
            }
            else if(m_singleMode) // single mode where d/dz = 0
            {
                fields[0]->PhysDeriv(m_baseflow[0], grad_base_u0,
                                     grad_base_u1);
                fields[0]->PhysDeriv(m_baseflow[1], grad_base_v0,
                                     grad_base_v1);
                fields[0]->PhysDeriv(m_baseflow[2], grad_base_w0,
                                     grad_base_w1);
            }
            else if(m_multipleModes)
            {
                // Differentiate base flow in physical space
                bool oldwavespace = fields[0]->GetWaveSpace();
                fields[0]->SetWaveSpace(false);
                fields[0]->PhysDeriv(m_baseflow[0], grad_base_u0,
                                     grad_base_u1,  grad_base_u2);
                fields[0]->PhysDeriv(m_baseflow[1], grad_base_v0,
                                     grad_base_v1,  grad_base_v2);
                fields[0]->PhysDeriv(m_baseflow[2], grad_base_w0,
                                     grad_base_w1,  grad_base_w2);
                fields[0]->SetWaveSpace(oldwavespace);
            }
            else
            {
                ASSERTL0(false, "ERROR: Must be one of half, single or multiple modes");
            }

            //x-equation
            //
            // Could probably clean up independent looping if clean up
            // base flow derivative evaluation
            if(m_multipleModes)
            {
                fields[0]->PhysDeriv(inarray[0], grad0, grad1, grad2);
                // transform gradients into physical Fourier space
                fields[0]->HomogeneousBwdTrans(grad0, grad0);
                fields[0]->HomogeneousBwdTrans(grad1, grad1);
                fields[0]->HomogeneousBwdTrans(grad2, grad2);
            }
            else
            {
                if(m_halfMode) //W = 0 so no need for d/dz
                {
                    fields[0]->PhysDeriv(inarray[0], grad0, grad1);
                }
                else // Single mode
                {
                    fields[0]->PhysDeriv(inarray[0], grad0, grad1, grad2);
                }
            }

            //Evaluate:  U du'/dx
            Vmath::Vmul (nPointsTot, grad0, 1,  m_baseflow[0], 1,
                         outarray[0], 1);
            //Evaluate and add: V du'/dy
            Vmath::Vvtvp(nPointsTot, grad1,  1,  m_baseflow[1], 1,
                         outarray[0], 1,  outarray[0],   1);
            //Evaluate and add: u' dU/dx
            Vmath::Vvtvp(nPointsTot,grad_base_u0,1,advVel[0],1,
                         outarray[0],1,outarray[0],1);
            //Evaluate and add: v' dU/dy
            Vmath::Vvtvp(nPointsTot,grad_base_u1,1,advVel[1],1,
                         outarray[0],1,outarray[0],1);

            if(!m_halfMode)
            {
                //Evaluate an add W du'/dz
                Vmath::Vvtvp(nPointsTot, grad2, 1, m_baseflow[2],
                             1,outarray[0], 1, outarray[0], 1);
            }

            if(m_multipleModes)
            {
                //Evaluate and add w' dU/dz
                Vmath::Vvtvp(nPointsTot,grad_base_u2,1,
                             advVel[2],1,outarray[0],1,outarray[0],1);
                fields[0]->HomogeneousFwdTrans(outarray[0],outarray[0]);
            }
            Vmath::Neg(nPointsTot,outarray[0],1);

            //y-equation------------
            if(m_multipleModes)
            {
                fields[0]->PhysDeriv(inarray[1], grad0, grad1, grad2);
                // transform gradients into physical fouier space
                fields[0]->HomogeneousBwdTrans(grad0, grad0);
                fields[0]->HomogeneousBwdTrans(grad1, grad1);
                fields[0]->HomogeneousBwdTrans(grad2, grad2);
            }
            else
            {
                if(m_halfMode) //W = 0 so no need for d/dz
                {
                    fields[0]->PhysDeriv(inarray[1], grad0, grad1);
                }
                else // Single mode
                {
                    fields[0]->PhysDeriv(inarray[1], grad0, grad1, grad2);
                }
            }

            //Evaluate U dv'/dx
            Vmath::Vmul (nPointsTot, grad0,     1,  m_baseflow[0], 1,
                         outarray[1], 1);
            //Evaluate V dv'/dy
            Vmath::Vvtvp(nPointsTot, grad1,     1,  m_baseflow[1], 1,
                         outarray[1], 1,  outarray[1],   1);
            //Evaluate u' dV/dx
            Vmath::Vvtvp(nPointsTot,grad_base_v0,1,advVel[0],1
                         ,outarray[1],1,outarray[1],1);
            //Evaluate v' dV/dy
            Vmath::Vvtvp(nPointsTot,grad_base_v1,1,advVel[1],1,
                         outarray[1],1,outarray[1],1);

            if(!m_halfMode)
            {
                //Evaluate  W du'/dz
                Vmath::Vvtvp(nPointsTot, grad2,       1,  m_baseflow[2], 1,
                             outarray[1], 1,  outarray[1],   1);
            }

            if(m_multipleModes)
            {
                //Evaluate  w' dV/dz
                Vmath::Vvtvp(nPointsTot,grad_base_v2,1,advVel[2],1,
                             outarray[1],1,outarray[1],1);
                fields[0]->HomogeneousFwdTrans(outarray[1],outarray[1]);
            }

            Vmath::Neg(nPointsTot,outarray[1],1);

            //z-equation ------------
            if(m_multipleModes)
            {
                fields[0]->PhysDeriv(inarray[2], grad0, grad1, grad2);
                // transform gradients into physical fouier space
                fields[0]->HomogeneousBwdTrans(grad0, grad0);
                fields[0]->HomogeneousBwdTrans(grad1, grad1);
                fields[0]->HomogeneousBwdTrans(grad2, grad2);
            }
            else
            {
                if(m_halfMode) //W = 0 so no need for d/dz
                {
                    fields[0]->PhysDeriv(inarray[2], grad0, grad1);
                }
                else // Single mode
                {
                    fields[0]->PhysDeriv(inarray[2], grad0, grad1, grad2);
                }
            }

            //Evaluate U dw'/dx
            Vmath::Vmul (nPointsTot, grad0,      1, m_baseflow[0], 1,
                         outarray[2],  1);
            //Evaluate V dw'/dx
            Vmath::Vvtvp(nPointsTot, grad1,      1, m_baseflow[1], 1,
                         outarray[2],  1, outarray[2],   1);

            if(!m_halfMode) // since if halfmode W = 0
            {
                //Evaluate u' dW/dx
                Vmath::Vvtvp(nPointsTot,grad_base_w0,1,advVel[0],1,
                             outarray[2],1,outarray[2],1);
                //Evaluate v' dW/dy
                Vmath::Vvtvp(nPointsTot,grad_base_w1,1,advVel[1],1,
                             outarray[2],1,outarray[2],1);

                //Evaluate W dw'/dz
                Vmath::Vvtvp(nPointsTot, grad2,       1, m_baseflow[2], 1,
                             outarray[2],  1, outarray[2],   1);
            }

            if(m_multipleModes)
            {
                //Evaluate w' dW/dz
                Vmath::Vvtvp(nPointsTot,grad_base_w2,1,advVel[2],1,
                             outarray[2],1,outarray[2],1);
                fields[0]->HomogeneousFwdTrans(outarray[2],outarray[2]);
            }
            Vmath::Neg(nPointsTot,outarray[2],1);
        }
    }
}

void LinearisedAdvection::v_SetBaseFlow(
        const Array<OneD, Array<OneD, NekDouble> >    &inarray)
{
    if (m_session->GetSolverInfo("EqType") == "UnsteadyNavierStokes")
    {
        // The SFD method is only applied to the velocity variables in
        // incompressible
        ASSERTL1(inarray.num_elements() == (m_baseflow.num_elements() - 1),
                 "Number of base flow variables does not match what is "
                 "expected.");
    }
    else
    {
        ASSERTL1(inarray.num_elements() == (m_baseflow.num_elements()),
             "Number of base flow variables does not match what is expected.");
    }

    int npts = inarray[0].num_elements();

    for (int i = 0; i < inarray.num_elements(); ++i)
    {
        ASSERTL1(npts == m_baseflow[i].num_elements(),
             "Size of base flow array does not match expected.");
        Vmath::Vcopy(npts, inarray[i], 1, m_baseflow[i], 1);
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

    int nqtot = m_baseflow[0].num_elements();
    Array<OneD, NekDouble> tmp_coeff(pFields[0]->GetNcoeffs(), 0.0);

    int numexp = pFields[0]->GetExpSize();
    Array<OneD,int> ElementGIDs(numexp);

    // Define list of global element ids
    for(int i = 0; i < numexp; ++i)
    {
        ElementGIDs[i] = pFields[0]->GetExp(i)->GetGeom()->GetGlobalID();
    }

    LibUtilities::FieldIOSharedPtr fld =
    MemoryManager<LibUtilities::FieldIO>::AllocateSharedPtr(
                                                    m_session->GetComm());
    fld->Import(pInfile, FieldDef, FieldData,
                LibUtilities::NullFieldMetaDataMap,
                ElementGIDs);

    int nSessionVar     = m_session->GetVariables().size();
    int nSessionConvVar = nSessionVar - 1;
    int nFileVar        = FieldDef[0]->m_fields.size();
    int nFileConvVar    = nFileVar - 1; // Ignore pressure
    if (m_halfMode)
    {
        ASSERTL0(nFileVar == 3, "For half/single mode, expect 2D2C base flow.");
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
    for (int j = nFileConvVar; j < nSessionConvVar; ++j) {
        Vmath::Fill(nqtot, 0.0, m_baseflow[j], 1);
    }

    // If time-periodic, put loaded data into the slice storage.
    if(m_session->DefinesParameter("N_slices"))
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
    int npoints     = m_baseflow[0].num_elements();
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

DNekBlkMatSharedPtr LinearisedAdvection::GetFloquetBlockMatrix(FloquetMatType mattype, bool UseContCoeffs) const
{
    DNekMatSharedPtr    loc_mat;
    DNekBlkMatSharedPtr BlkMatrix;
    int n_exp = 0;

    n_exp = m_baseflow[0].num_elements(); // will operatore on m_phys

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
    int ConvectedFields = m_baseflow.num_elements()-1;
    int npoints         = m_baseflow[0].num_elements();
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
    for (int i = 0; i < m_slices; ++i)
    {
        sprintf(buffer, file.c_str(), i);
        ImportFldBase(buffer,pFields,i);
    }
    delete[] buffer;


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

        Vmath::Zero(fft_in.num_elements(),&fft_in[0],1);
        Vmath::Zero(fft_out.num_elements(),&fft_out[0],1);
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

        for(int r=0; r<sortedinarray.num_elements();++r)
        {
            sortedinarray[0]=0;
            sortedoutarray[0]=0;
        }

#endif

        //scaling of the Fourier coefficients
        NekDouble j=-1;
        for (int i = 2; i < m_slices; i += 2)
        {
            Vmath::Smul(2*npoints,j,&m_interp[k][i*npoints],1,&m_interp[k][i*npoints],1);
            j=-j;

        }

    }

}

} //end of namespace

