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
// Description: Abstract base class for advection.
//
///////////////////////////////////////////////////////////////////////////////

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


void Advection::v_AddTraceJac2Mat(
    const int                                          nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, DNekBlkMatSharedPtr>            &TraceJac,
    DNekBlkMatSharedPtr &gmtx)
{
    ASSERTL0(false,"v_AddTraceJac2Mat NOT SPECIFIED");
    return;
}


void Advection::v_AddVolumJac2Mat( const int nConvectiveFields,
                                const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                                const   Array<OneD, const  Array<OneD, DNekMatSharedPtr> >&ElmtJac,
                                const int nDirctn, DNekBlkMatSharedPtr &gmtx)
{
    ASSERTL0(false,"v_AddVolumJac2Mat NOT SPECIFIED");
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

     /*      

    if(FALSE)
    {
        Array<OneD, NekDouble> PointFwd(nConvectiveFields,0.0),PointBwd(nConvectiveFields,0.0);
        Array<OneD, NekDouble> PntFluxFwd(nConvectiveFields,0.0);
        Array<OneD, NekDouble> PntFluxBwd(nConvectiveFields,0.0);
        Array<OneD, NekDouble> PointFlux(nConvectiveFields,0.0);
        int cnt=0;
        for(int i=0; i<nTracePointsTot;  i++)
        {
            for(int j=0; j<nConvectiveFields;j++)
            {
                PointFwd[j] = Fwd[j][i];
                PointBwd[j] = Bwd[j][i];
            }
            NekVector<NekDouble> VectFluxFwd(nConvectiveFields,PntFluxFwd,eWrapper);
            NekVector<NekDouble> VectFluxBwd(nConvectiveFields,PntFluxBwd,eWrapper);

            DNekMat &MF = (*FJac->GetBlock(i,i));
            NekVector<NekDouble> VectFwd(nConvectiveFields,PointFwd,eWrapper);
            VectFluxFwd = MF * VectFwd;


            DNekMat &MB = (*BJac->GetBlock(i,i));
            NekVector<NekDouble> VectBwd(nConvectiveFields,PointBwd,eWrapper);
            VectFluxBwd = MB * VectBwd;

            NekDouble error=0.0;
            for(int j=0;j<nConvectiveFields;j++)
            {
                PointFlux[j] = PntFluxFwd[j]+PntFluxBwd[j];
                error += abs(PointFlux[j]-numflux[j][i]);
            }

            if(error>1.0E-7)
            {
                cnt++;
                std::cout   <<std::scientific<<std::setw(12)<<std::setprecision(5)
                        <<"abs(PointFlux[0]-numflux[0][i])   =   "<<abs(PointFlux[0]-numflux[0][i])<<"    "<<PointFlux[0]<<"    "<<numflux[0][i]<<std::endl
                        <<"abs(PointFlux[1]-numflux[1][i])   =   "<<abs(PointFlux[1]-numflux[1][i])<<"    "<<PointFlux[1]<<"    "<<numflux[1][i]<<std::endl
                        <<"abs(PointFlux[2]-numflux[2][i])   =   "<<abs(PointFlux[2]-numflux[2][i])<<"    "<<PointFlux[2]<<"    "<<numflux[2][i]<<std::endl
                        <<"abs(PointFlux[3]-numflux[3][i])   =   "<<abs(PointFlux[3]-numflux[3][i])<<"    "<<PointFlux[3]<<"    "<<numflux[3][i]<<std::endl;
            }
        }
        
        std::cout   <<"cnt= "<<cnt<<std::endl;
        
        int j = 0;
    }

 */
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
            ASSERTL0(false, "Only 1D homogeneous dimension supported.");
        }
    }
}


/**
 *
 */
void Advection::v_SetBaseFlow(
        const Array<OneD, Array<OneD, NekDouble> >    &inarray,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields)
{
    ASSERTL0(false,
            "A baseflow is not appropriate for this advection type.");
}

}
}
