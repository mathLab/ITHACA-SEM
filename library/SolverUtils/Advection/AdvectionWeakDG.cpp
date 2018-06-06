///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionWeakDG.cpp
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
// Description: Weak DG advection class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Advection/AdvectionWeakDG.h>
#include <LocalRegions/MatrixKey.h>
#include <iostream>
#include <iomanip>
#include <math.h>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string AdvectionWeakDG::type = GetAdvectionFactory().
            RegisterCreatorFunction("WeakDG", AdvectionWeakDG::create);

        AdvectionWeakDG::AdvectionWeakDG()
        {
        }

        /**
         * @brief Initialise AdvectionWeakDG objects and store them before
         * starting the time-stepping.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void AdvectionWeakDG::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            Advection::v_InitObject(pSession, pFields);
        }

        /**
         * @brief Compute the advection term at each time-step using the
         * Discontinuous Galerkin approach (DG).
         *
         * @param nConvectiveFields   Number of fields.
         * @param fields              Pointer to fields.
         * @param advVel              Advection velocities.
         * @param inarray             Solution at the previous time-step.
         * @param outarray            Advection term to be passed at the
         *                            time integration class.
         */
        void AdvectionWeakDG::v_Advect(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &advVel,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray,
            const NekDouble                                   &time,
            const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
            const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
        {
            int nPointsTot      = fields[0]->GetTotPoints();
            int nCoeffs         = fields[0]->GetNcoeffs();
            int nTracePointsTot = fields[0]->GetTrace()->GetTotPoints();
            int i, j;
            
            Array<OneD, Array<OneD, NekDouble> > tmp(nConvectiveFields);
            Array<OneD, Array<OneD, Array<OneD, NekDouble> > > fluxvector(
                nConvectiveFields);

            // Allocate storage for flux vector F(u).
            for (i = 0; i < nConvectiveFields; ++i)
            {
                fluxvector[i] =
                    Array<OneD, Array<OneD, NekDouble> >(m_spaceDim);
                for (j = 0; j < m_spaceDim; ++j)
                {
                    fluxvector[i][j] = Array<OneD, NekDouble>(nPointsTot);
                }
            }

            ASSERTL1(m_riemann,
                     "Riemann solver must be provided for AdvectionWeakDG.");

            m_fluxVector(inarray, fluxvector);

            // Get the advection part (without numerical flux)
            for(i = 0; i < nConvectiveFields; ++i)
            {
                tmp[i] = Array<OneD, NekDouble>(nCoeffs, 0.0);

                fields[i]->IProductWRTDerivBase(fluxvector[i],tmp[i]);
            }

            // Store forwards/backwards space along trace space
            Array<OneD, Array<OneD, NekDouble> > Fwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > Bwd    (nConvectiveFields);
            Array<OneD, Array<OneD, NekDouble> > numflux(nConvectiveFields);

            if (pFwd == NullNekDoubleArrayofArray ||
                pBwd == NullNekDoubleArrayofArray)
            {
                for(i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    Bwd[i]     = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    numflux[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                    fields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
                }
            }
            else
            {
                for(i = 0; i < nConvectiveFields; ++i)
                {
                    Fwd[i]     = pFwd[i];
                    Bwd[i]     = pBwd[i];
                    numflux[i] = Array<OneD, NekDouble>(nTracePointsTot, 0.0);
                }
            }

            m_riemann->Solve(m_spaceDim, Fwd, Bwd, numflux);

  /*           

            if(FALSE)
            {

                Array<OneD, unsigned int> n_blks(nTracePointsTot);
                for(int i=0;i<nTracePointsTot;i++)
                {
                    n_blks[i]    = nConvectiveFields;
                }
                DNekBlkMatSharedPtr FJac = MemoryManager<DNekBlkMat>
                    ::AllocateSharedPtr(n_blks, n_blks, eDIAGONAL);
                DNekBlkMatSharedPtr BJac = MemoryManager<DNekBlkMat>
                    ::AllocateSharedPtr(n_blks, n_blks, eDIAGONAL);
                m_riemann->CalcFluxJacobian(m_spaceDim, Fwd, Bwd, FJac,BJac);




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

            // Evaulate <\phi, \hat{F}\cdot n> - OutField[i]
            for(i = 0; i < nConvectiveFields; ++i)
            {
                Vmath::Neg                      (nCoeffs, tmp[i], 1);
                fields[i]->AddTraceIntegral     (numflux[i], tmp[i]);
                fields[i]->MultiplyByElmtInvMass(tmp[i], tmp[i]);
                fields[i]->BwdTrans             (tmp[i], outarray[i]);
            }
        }

        
        void AdvectionWeakDG::v_AddVolumJac2Mat( const int nConvectiveFields,
                                        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                                        const   Array<OneD, const  Array<OneD, DNekMatSharedPtr> >&ElmtJac,
                                        const int nDirctn, DNekBlkMatSharedPtr &gmtx)
        {

            std::shared_ptr<LocalRegions::ExpansionVector> pexp;
            pexp        = pFields[0]->GetExp();
            int ntotElmt            = (*pexp).size();
            int nElmtPnt            = (*pexp)[0]->GetTotPoints();
            int nElmtCoef           = (*pexp)[0]->GetNcoeffs();
            // int nCoefVar            =   nElmtCoef*nConvectiveFields;

            StdRegions::ConstFactorMap  factors;
            StdRegions::VarCoeffMap     VarCoeff;

            NekDouble tmp;

            DNekMatSharedPtr    ElmtMat = MemoryManager<DNekMat>
                ::AllocateSharedPtr(nElmtCoef, nElmtCoef);

            Array<OneD, NekDouble>  tmpPhys(nElmtPnt,0.0);
            Array<OneD, NekDouble>  tmpCoef(nElmtCoef,0.0);

            // DNekScalMatSharedPtr    loc_ScalmatNvar;
            // DNekMatSharedPtr        loc_matNvar;
            DNekMatSharedPtr        tmpElmtMat ;

            

            for(int  nelmt = 0; nelmt < ntotElmt; nelmt++)
            {
                nElmtCoef           = (*pexp)[nelmt]->GetNcoeffs();
                nElmtPnt            = (*pexp)[nelmt]->GetTotPoints();

                tmpElmtMat       =   gmtx->GetBlock(nelmt,nelmt);
                // tmpElmtMat          =   tmpElmtSclMat->GetOwnedMatrix();

                // int nrowsVars   = nElmtCoef*nConvectiveFields;
                // int ncolsVars   = nElmtCoef*nConvectiveFields;
                // loc_matNvar = MemoryManager<DNekMat>::AllocateSharedPtr(nrowsVars,ncolsVars,0.0);

                if (tmpPhys.num_elements()!=nElmtPnt||tmpCoef.num_elements()!=nElmtCoef) 
                {
                    ElmtMat = MemoryManager<DNekMat>
                            ::AllocateSharedPtr(nElmtCoef, nElmtCoef);
                    
                    tmpPhys     =   Array<OneD, NekDouble>(nElmtPnt,0.0);
                    tmpCoef     =   Array<OneD, NekDouble>(nElmtCoef,0.0);

                    // nCoefVar    =   nElmtCoef*nConvectiveFields;
                }

                // LocalRegions::MatrixKey matkey(StdRegions::eBwdTrans,
                //                             (*pexp)[nelmt]->DetShapeType(),
                //                             *((*pexp)[nelmt]), factors, VarCoeff);

                StdRegions::StdMatrixKey matkey(StdRegions::eFwdTrans,
                                            (*pexp)[nelmt]->.DetShapeType(),
                                             *((*pexp)[nelmt]));

                DNekScalMatSharedPtr BwdTransMat =  (*pexp)[nelmt]->GetStdMatrix(matkey);

                for(int m = 0; m < nConvectiveFields; m++)
                {
                    for(int n = 0; n < nConvectiveFields; n++)
                    {
                        for(int ncl = 0; ncl < nElmtCoef; ncl++)
                        {
                            for(int npnt = 0; npnt < nElmtPnt; npnt++)
                            {
                                tmpPhys[npnt]   =   (*(ElmtJac[nelmt][npnt]))(m,n)*(*BwdTransMat)(npnt,ncl);
                            }

                            (*pexp)[nelmt]->IProductWRTDerivBase(nDirctn,tmpPhys,tmpCoef);

                            for(int nrw = 0; nrw < nElmtCoef; nrw++)
                            {
                                (*ElmtMat)(nrw,ncl)   =   tmpCoef[nrw];
                            }
                        }

                        for(int ncl = 0; ncl < nElmtCoef; ncl++)
                        {
                            int nclVar = ncl*nConvectiveFields+m;

                            for(int nrw = 0; nrw < nElmtCoef; nrw++)
                            {
                                int nrwVar = nrw*nConvectiveFields+n;
                                tmp   =   (*tmpElmtMat)(nrwVar,nclVar) + (*ElmtMat)(nrw,ncl);
                                tmpElmtMat->SetValue(nrwVar,nclVar,tmp);
                            }
                        }

                    }
                    
                }

                // loc_ScalmatNvar = MemoryManager<DNekScalMat>::
                //             AllocateSharedPtr(1.0,loc_matNvar);
                // gmtx->SetBlock(nelmt,nelmt,loc_ScalmatNvar);
                
            }
      
        }

        void AdvectionWeakDG::v_AddTraceJac2Mat(
            const int                                          nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
            const Array<OneD, DNekBlkMatSharedPtr>            &TraceJac,
            DNekBlkMatSharedPtr &gmtx)
        {
            std::shared_ptr<LocalRegions::ExpansionVector> pexp;
            pexp                = pFields[0]->GetTrace();
            int ntotTrace       = (*pexp).size();
            int nTracePointsTot = pexp->GetTotPoints();
            int nTraceCoeffsTot = pexp->GetNcoeffs();
            
            Array<TwoD, int      > BlkId(nTraceCoeffsTot,0  );
            Array<TwoD, int      > IdArray(nTraceCoeffsTot,0  );

            int MaxElmtNcoef   =   (*pexp)[0]->GetNcoeffs();
            int MaxElmtNphys   =   (*pexp)[0]->GetTotPoints();
            
            for(int  nelmt = 0; nelmt < ntotTrace; nelmt++)
            {
                nElmtCoef           = (*pexp)[nelmt]->GetNcoeffs();
                nElmtPnt            = (*pexp)[nelmt]->GetTotPoints();

                MaxElmtNcoef = std::max(MaxElmtNcoef,nElmtCoef);
                MaxElmtNphys = std::max(MaxElmtNphys,nElmtPnt);
            }
            WARNINGL0(  (*pexp)[0]->GetTotPoints()==MaxElmtNphys   ,    "traces have nonuniform phys points, needs test");
            WARNINGL0(  (*pexp)[0]->GetNcoeffs()==MaxElmtNcoef   ,    "traces have nonuniform modes coeff, needs test");
           
            Array<OneD, Array<OneD, NekDouble> >    tmpPhys(MaxElmtNcoef);
            for(int i = 0; i < MaxElmtNcoef; i++)
            {
                tmpPhys[i]  =   Array<OneD, NekDouble> tmpPhys(nTracePointsTot,0.0);
            }
            
            int ncoef_offset,nphys_offset;
            Array<OneD, DNekBlkMatSharedPtr> tmpelmtJac;
            
            for(int m = 0; m < nConvectiveFields; ++m)
            {
                for(int n = 0; n < nConvectiveFields; ++n)
                {
                    for(int  nelmt = 0; nelmt < ntotTrace; nelmt++)
                    {
                        nElmtCoef           = (*pexp)[nelmt]->GetNcoeffs();
                        nElmtPnt            = (*pexp)[nelmt]->GetTotPoints();

                        ncoef_offset        =  pexp->GetCoeff_Offset(nelmt);
                        nphys_offset        =  pexp->GetPhys_Offset(nelmt);
                        tmpelmtJac          =  TraceJac[nphys_offset];

                        LocalRegions::MatrixKey matkey(StdRegions::eBwdTrans,
                                                    (*pexp)[nelmt]->DetShapeType(),
                                                    *((*pexp)[nelmt]), factors, VarCoeff);

                        DNekScalMatSharedPtr BwdTransMat =  (*pexp)[nelmt]->GetLocMatrix(matkey);

                        for(int ncolms = 0; ncolms < nElmtCoef; ncolms++)
                        {
                            for(int nrows = 0; nrows < nElmtPnt; nrows++)
                            {
                                tmpPhys[ncolms][nrows] = (*tmpelmtJac[nElmtPnt])(m,n)*(*BwdTransMat)(nElmtPnt,nElmtCoef);
                            }
                        }

                        /// in case nonuniform trace coeff 
                        for(int ncolms = nElmtCoef; ncolms < MaxElmtNcoef; ncolms++)
                        {
                            for(int nrows = 0; nrows < nElmtPnt; nrows++)
                            {
                                tmpPhys[ncolms][nrows] = 0.0;
                            }
                        }
                    }

                    fields[m]->AddTraceIntegral     (tmpPhys[i], tmp[m]);





                    fields[m]->MultiplyByElmtInvMass(tmp[m], tmp[m]);
                    fields[m]->BwdTrans             (tmp[m], outarray[m]);
                }
               
            }
        }
    }//end of namespace SolverUtils
}//end of namespace Nektar
