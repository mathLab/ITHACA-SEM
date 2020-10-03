///////////////////////////////////////////////////////////////////////////////
//
// File: Advection.cpp
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
// Description: Abstract base class for advection.
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <SolverUtils/Advection/Advection.h>

namespace Nektar
{
namespace SolverUtils
{

/**
 * @returns The advection factory.
 */
AdvectionFactory& GetAdvectionFactory()
{
    static AdvectionFactory instance;
    return instance;
}


/**
 * @param   pSession            Session configuration data.
 * @param   pFields             Array of ExpList objects.
 */
void Advection::InitObject(
    const LibUtilities::SessionReaderSharedPtr        pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
{
    v_InitObject(pSession, pFields);
}


/**
 * @param   nConvectiveFields   Number of velocity components.
 * @param   pFields             Expansion lists for scalar fields.
 * @param   pAdvVel             Advection velocity.
 * @param   pInarray            Scalar data to advect.
 * @param   pOutarray           Advected scalar data.
 * @param   pTime               Simulation time.
 */
void Advection::Advect(
    const int                                          nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble> >        &pAdvVel,
    const Array<OneD, Array<OneD, NekDouble> >        &pInarray,
    Array<OneD, Array<OneD, NekDouble> >              &pOutarray,
    const NekDouble                                   &pTime,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
{
    v_Advect(nConvectiveFields, pFields, pAdvVel, pInarray,
            pOutarray, pTime, pFwd, pBwd);
}

template<typename DataType, typename TypeNekBlkMatSharedPtr>
void Advection::AddTraceJacToMat(
        const int                                           nConvectiveFields,
        const int                                           nSpaceDim,
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields,
        const Array<OneD, TypeNekBlkMatSharedPtr>           &TracePntJacCons,
        Array<OneD, Array<OneD, TypeNekBlkMatSharedPtr> >   &gmtxarray,
        const Array<OneD, TypeNekBlkMatSharedPtr>           &TracePntJacGrad        ,
        const Array<OneD, Array<OneD, DataType> >           &TracePntJacGradSign    )
{
    MultiRegions::ExpListSharedPtr tracelist = pFields[0]->GetTrace();
    std::shared_ptr<LocalRegions::ExpansionVector> traceExp= tracelist->GetExp();
    int ntotTrac            = (*traceExp).size();
    int nTracePnts          = tracelist->GetTotPoints();
    int nTracPnt,nTracCoef;

    DNekMatSharedPtr tmp2Add;
    Array<OneD, DataType>    MatData0;
    Array<OneD, NekDouble>    MatData1;
    Array<OneD, NekDouble>    MatData2;

    Array<OneD, DNekMatSharedPtr>  TraceJacFwd(ntotTrac);
    Array<OneD, DNekMatSharedPtr>  TraceJacBwd(ntotTrac);

    for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
    {
        nTracCoef           = (*traceExp)[nelmt]->GetNcoeffs();
        nTracPnt            = (*traceExp)[nelmt]->GetTotPoints();
        TraceJacFwd[nelmt]    =MemoryManager<DNekMat>
                            ::AllocateSharedPtr(nTracCoef, nTracPnt,0.0);
        TraceJacBwd[nelmt]    =MemoryManager<DNekMat>
                            ::AllocateSharedPtr(nTracCoef, nTracPnt,0.0);
    }

    std::shared_ptr<LocalRegions::ExpansionVector> fieldExp= pFields[0]->GetExp();
    int ntotElmt            = (*fieldExp).size();
    int nElmtPnt,nElmtCoef;

    Array<OneD, DNekMatSharedPtr>  ElmtJacQuad(ntotElmt);
    Array<OneD, DNekMatSharedPtr>  ElmtJacCoef(ntotElmt);

    Array<OneD, NekDouble> SymmMatData;
    
    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
    {
        nElmtCoef          = (*fieldExp)[nelmt]->GetNcoeffs();
        nElmtPnt           = (*fieldExp)[nelmt]->GetTotPoints();
        
        ElmtJacQuad[nelmt]    =MemoryManager<DNekMat>
                            ::AllocateSharedPtr(nElmtCoef, nElmtPnt,0.0);
        ElmtJacCoef[nelmt]    =MemoryManager<DNekMat>
                            ::AllocateSharedPtr(nElmtCoef, nElmtCoef,0.0);
    }

    bool TracePntJacGradflag = true;

    Array<OneD, Array<OneD, DataType> > TraceJacConsSign(2);
    for(int i = 0; i < 2; i++)
    {
        TraceJacConsSign[i] =   Array<OneD, DataType>(nTracePnts,1.0);
    }

    if(0==TracePntJacGrad.size())
    {
        TracePntJacGradflag = false;
    }

    for(int m = 0; m < nConvectiveFields; m++)
    {
        for(int n = 0; n < nConvectiveFields; n++)
        {
            // ElmtJacCons to set 0
            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                (*ElmtJacCoef[nelmt]) =  0.0;
                (*ElmtJacQuad[nelmt]) =  0.0;
            }

            if(TracePntJacGradflag)
            {
                // TODO: only 1 BJac is stored here because BJac = - FJac
                for(int ndir = 0; ndir < nSpaceDim; ndir++)
                {
                    // ElmtJacGrad to set 0
                    for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                    {
                        (*ElmtJacQuad[nelmt]) =  0.0;
                    }
                    int ngrad = n*nSpaceDim+ndir;
                    
                    CalcJacobTraceInteg(pFields,m,ngrad,TracePntJacGrad,TracePntJacGradSign,TraceJacFwd,TraceJacBwd);
                    pFields[0]->AddTraceJacToElmtJac(TraceJacFwd,TraceJacBwd,ElmtJacQuad);
                    //TODO:: 3 directions together to lower down cost
                    pFields[0]->AddRightIPTPhysDerivBase(ndir,ElmtJacQuad,ElmtJacCoef);
                }
                for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
                {
                    nElmtCoef          = (*fieldExp)[nelmt]->GetNcoeffs();
                    if(SymmMatData.size()<nElmtCoef)
                    {
                        SymmMatData = Array<OneD, NekDouble> (nElmtCoef*nElmtCoef);
                    }
                    MatData1 =  ElmtJacCoef[nelmt]->GetPtr();
                    for(int i=0;i<nElmtCoef;i++)
                    {
                        Vmath::Vcopy(nElmtCoef,&MatData1[i*nElmtCoef],1,&SymmMatData[i],nElmtCoef);
                    }
                    Vmath::Vadd(nElmtCoef*nElmtCoef,SymmMatData,1,MatData1,1,MatData1,1);
                    // Vmath::Smul(nElmtCoef*nElmtCoef,0.5,MatData1,1,MatData1,1);
                }
            }

            CalcJacobTraceInteg(pFields,m,n,TracePntJacCons, TraceJacConsSign,TraceJacFwd,TraceJacBwd);
            //TODO: To modify CalcJacobTraceInteg&AddTraceJacToElmtJac to Bwd after Fwd not at the same time
            pFields[0]->AddTraceJacToElmtJac(TraceJacFwd, TraceJacBwd,ElmtJacQuad);
            //TODO: To check whether it is ok to reuse ElmtJacQuad as output
            pFields[0]->AddRightIPTBaseMatrix(ElmtJacQuad,ElmtJacCoef);

            // add ElmtJacCons to gmtxarray[m][n]
            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                tmp2Add         = ElmtJacCoef[nelmt];
                MatData0        = gmtxarray[m][n]->GetBlock(nelmt,nelmt)->GetPtr();
                MatData1        = tmp2Add->GetPtr();
                for(int i=0;i<MatData0.size();i++)
                {
                    MatData0[i] += DataType(MatData1[i]);
                }

                // Vmath::Vadd(MatData0.size(),MatData0,1,MatData1,1,MatData0,1);
                // (*tmpGmtx)      = (*tmpGmtx)  +   (*tmp2Add);
            }
        }
    }
}

template SOLVER_UTILS_EXPORT void Advection::AddTraceJacToMat( 
        const int                                           nConvectiveFields,
        const int                                           nSpaceDim,
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields,
        const Array<OneD, DNekBlkMatSharedPtr>              &TracePntJacCons,
        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >      &gmtxarray,
        const Array<OneD, DNekBlkMatSharedPtr>              &TracePntJacGrad        ,
        const Array<OneD, Array<OneD, NekDouble> >          &TracePntJacGradSign    );
template SOLVER_UTILS_EXPORT void Advection::AddTraceJacToMat( 
        const int                                           nConvectiveFields,
        const int                                           nSpaceDim,
        const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields,
        const Array<OneD, SNekBlkMatSharedPtr>              &TracePntJacCons,
        Array<OneD, Array<OneD, SNekBlkMatSharedPtr> >      &gmtxarray,
        const Array<OneD, SNekBlkMatSharedPtr>              &TracePntJacGrad        ,
        const Array<OneD, Array<OneD, NekSingle> >          &TracePntJacGradSign    );


// TODO: To change it into one by one
template<typename DataType, typename TypeNekBlkMatSharedPtr>
void Advection::CalcJacobTraceInteg(
    const Array<OneD, MultiRegions::ExpListSharedPtr>   &pFields,
    const int                                             m,
    const int                                             n,
    const Array<OneD, const TypeNekBlkMatSharedPtr>       & PntJac,
    const Array<OneD, const Array<OneD, DataType> >       & PntJacSign,
    Array<OneD, DNekMatSharedPtr>                         & TraceJacFwd,
    Array<OneD, DNekMatSharedPtr>                         & TraceJacBwd)
{
    MultiRegions::ExpListSharedPtr tracelist = pFields[0]->GetTrace();
    std::shared_ptr<LocalRegions::ExpansionVector> traceExp= tracelist->GetExp();
    int ntotTrac            = (*traceExp).size();
    int nTracPnt,noffset,pntoffset;

    Array<OneD, int > tracepnts(ntotTrac);
    Array<OneD, Array<OneD, NekDouble> > JacFwd(ntotTrac);
    Array<OneD, Array<OneD, NekDouble> > JacBwd(ntotTrac);

    for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
    {
        nTracPnt            = (*traceExp)[nelmt]->GetTotPoints();
        tracepnts[nelmt]     =   nTracPnt;

        JacFwd[nelmt]     =Array<OneD, NekDouble>(nTracPnt,0.0);
        JacBwd[nelmt]     =Array<OneD, NekDouble>(nTracPnt,0.0);
    }

    // DNekMatSharedPtr FtmpMat,BtmpMat;
    DataType ftmp,btmp;
    for(int  nelmt = 0; nelmt < ntotTrac; nelmt++)
    {
        nTracPnt            =   tracepnts[nelmt];
        noffset             =   tracelist->GetPhys_Offset(nelmt);
        for(int npnt = 0; npnt < nTracPnt; npnt++)
        {
            pntoffset = noffset+npnt;
            ftmp                    =   (*PntJac[0]->GetBlock(pntoffset,pntoffset))(m,n);
            JacFwd[nelmt][npnt]     =   NekDouble(PntJacSign[0][pntoffset]*ftmp);

            btmp                    =   (*PntJac[1]->GetBlock(pntoffset,pntoffset))(m,n);
            JacBwd[nelmt][npnt]     =   NekDouble(PntJacSign[1][pntoffset]*btmp);
        }
    }
    tracelist->GetDiagMatIpwrtBase(JacFwd,TraceJacFwd);
    tracelist->GetDiagMatIpwrtBase(JacBwd,TraceJacBwd);
}


void Advection::v_AdvectVolumeFlux(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble>>         &pAdvVel,
    const Array<OneD, Array<OneD, NekDouble>>         &pInarray,
    TensorOfArray3D<NekDouble>                        &pVolumeFlux,
    const NekDouble                                   &pTime)
{
    boost::ignore_unused(nConvectiveFields, pFields, pAdvVel, pInarray,
                        pVolumeFlux, pTime);
    ASSERTL0(false, "Not defined for AdvectVolumeFlux.");
}

/**
 * @brief calculate the advection flux in the cell the trace  integration
 */
void Advection::v_AdvectTraceFlux(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble>>         &pAdvVel,
    const Array<OneD, Array<OneD, NekDouble>>         &pInarray,
    Array<OneD, Array<OneD, NekDouble>>               &pTraceFlux,
    const NekDouble                                   &pTime,
    const Array<OneD, Array<OneD, NekDouble>>         &pFwd,
    const Array<OneD, Array<OneD, NekDouble>>         &pBwd)
{
    boost::ignore_unused(nConvectiveFields, pFields, pAdvVel, pInarray,
                        pTraceFlux, pTime, pFwd, pBwd);
    ASSERTL0(false, "Not defined for AdvectTraceFlux.");
}

/**
 * @brief Similar with Advection::Advect(): calculate the advection flux
 * The difference is in the outarray:
 *  it is the coefficients of basis for AdvectCoeffs()
 *  it is the physics on quadrature points for Advect()
 *
 * @param   nConvectiveFields   Number of velocity components.
 * @param   pFields             Expansion lists for scalar fields.
 * @param   pAdvVel             Advection velocity.
 * @param   pInarray            Scalar data to advect.
 * @param   pOutarray           Advected scalar data.
 * @param   pTime               Simulation time.
 */
void Advection::AdvectCoeffs(
    const int                                          nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble> >        &pAdvVel,
    const Array<OneD, Array<OneD, NekDouble> >        &pInarray,
    Array<OneD, Array<OneD, NekDouble> >              &pOutarray,
    const NekDouble                                   &pTime,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
{
    v_AdvectCoeffs(nConvectiveFields, pFields, pAdvVel, pInarray,
                   pOutarray, pTime, pFwd, pBwd);
}

/**
 * This function should be overridden in derived classes to initialise the
 * specific advection data members. However, this base class function should
 * be called as the first statement of the overridden function to ensure the
 * base class is correctly initialised in order.
 *
 * @param   pSession            Session information.
 * @param   pFields             Expansion lists for scalar fields.
 */
void Advection::v_InitObject(
    const LibUtilities::SessionReaderSharedPtr        pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr>       pFields)
{
    m_spaceDim = pFields[0]->GetCoordim(0);

    if (pSession->DefinesSolverInfo("HOMOGENEOUS"))
    {
        std::string HomoStr = pSession->GetSolverInfo("HOMOGENEOUS");
        if (HomoStr == "HOMOGENEOUS1D" || HomoStr == "Homogeneous1D" ||
            HomoStr == "1D"            || HomoStr == "Homo1D" ||
            HomoStr == "HOMOGENEOUS2D" || HomoStr == "Homogeneous2D" ||
            HomoStr == "2D"            || HomoStr == "Homo2D")
        {
            m_spaceDim++;
        }
        else
        {
            NEKERROR(ErrorUtil::efatal,
                     "Only 1D homogeneous dimension supported.");
        }
    }
}


/**
 *
 */
void Advection::v_SetBaseFlow(
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
{
    boost::ignore_unused(inarray, fields);
    NEKERROR(ErrorUtil::efatal,
            "A baseflow is not appropriate for this advection type.");
}

void Advection::v_AdvectCoeffs(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble> >        &advVel,
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
    Array<OneD, Array<OneD, NekDouble> >              &outarray,
    const NekDouble                                   &time,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
{
    boost::ignore_unused(nConvectiveFields, fields, advVel, inarray, outarray,
                        time, pFwd, pBwd);
    ASSERTL0(false, "v_AdvectCoeffs not defined");
}

void Advection::v_NumCalRiemFluxJac( 
        const int                                          nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &AdvVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
        DNekBlkMatSharedPtr &FJac,
        DNekBlkMatSharedPtr &BJac)
{
    boost::ignore_unused(nConvectiveFields, fields, AdvVel, inarray,  pFwd, pBwd,
                         FJac,BJac);
    ASSERTL0(false, "v_NumCalRiemFluxJac no defined");
}

void Advection::v_AddVolumJacToMat( 
        const Array<OneD, MultiRegions::ExpListSharedPtr>     &pFields,
        const int                                             &nConvectiveFields,
        const Array<OneD, const Array<OneD,  Array<OneD,
        Array<OneD,  Array<OneD,  NekDouble> > > > >          &ElmtJacArray,
        Array<OneD, Array<OneD, DNekBlkMatSharedPtr> >        &gmtxarray)
{
    boost::ignore_unused(pFields,nConvectiveFields, ElmtJacArray,gmtxarray);
    ASSERTL0(false,"v_AddVolumJacToMat NOT SPECIFIED");
    return;
}

void Advection::v_AddVolumJacToMat( 
        const Array<OneD, MultiRegions::ExpListSharedPtr>     &pFields,
        const int                                             &nConvectiveFields,
        const Array<OneD, const Array<OneD,  Array<OneD,
        Array<OneD,  Array<OneD,  NekDouble> > > > >          &ElmtJacArray,
        Array<OneD, Array<OneD, SNekBlkMatSharedPtr> >        &gmtxarray)
    {
        boost::ignore_unused(pFields,nConvectiveFields,ElmtJacArray,gmtxarray);
        ASSERTL0(false,"v_AddVolumJacToMat NOT SPECIFIED");
        return;
    }

    void Advection::CalcTraceJac(
        const int                                          nConvectiveFields,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble> >        &AdvVel,
        const Array<OneD, Array<OneD, NekDouble> >        &inarray,
        const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
        const Array<OneD, Array<OneD, NekDouble> >        &pBwd,
        DNekBlkMatSharedPtr &FJac,
        DNekBlkMatSharedPtr &BJac)
    {
        boost::ignore_unused(AdvVel);
        // int nPointsTot      = fields[0]->GetTotPoints();
        // int nCoeffs         = fields[0]->GetNcoeffs();
        int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
        
        // Store forwards/backwards space along trace space
        Array<OneD, Array<OneD, NekDouble> > Fwd    (nConvectiveFields);
        Array<OneD, Array<OneD, NekDouble> > Bwd    (nConvectiveFields);

        if (pFwd == NullNekDoubleArrayofArray ||
            pBwd == NullNekDoubleArrayofArray)
        {
            for(int i = 0; i < nConvectiveFields; ++i)
            {
                Fwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                Bwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
            }
        }
        else
        {
            for(int i = 0; i < nConvectiveFields; ++i)
            {
                Fwd[i]     = pFwd[i];
                Bwd[i]     = pBwd[i];
            }
        }

        ASSERTL1(m_riemann,
                    "Riemann solver must be provided for AdvectionWeakDG.");
    
        m_riemann->CalcFluxJacobian(m_spaceDim, Fwd, Bwd, FJac,BJac);
    }

}
}
