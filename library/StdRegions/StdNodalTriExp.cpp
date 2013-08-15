///////////////////////////////////////////////////////////////////////////////
//
// File StdNodalTriExp.cpp
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
// Description: Nodal triangle routines built upon StdExpansion2D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdNodalTriExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>  // for PointsManager, etc

namespace Nektar
{
    namespace StdRegions
    {
        StdNodalTriExp::StdNodalTriExp():
            StdTriExp(),
            m_nodalPointsKey()
        {
        }

        StdNodalTriExp::StdNodalTriExp(
            const LibUtilities::BasisKey &Ba, 
            const LibUtilities::BasisKey &Bb, 
            LibUtilities::PointsType Ntype):
            StdExpansion  (LibUtilities::StdTriData::getNumberOfCoefficients(
                               Ba.GetNumModes(),
                               Bb.GetNumModes()),
                           2,Ba,Bb),
            StdExpansion2D(LibUtilities::StdTriData::getNumberOfCoefficients(
                               Ba.GetNumModes(),
                               Bb.GetNumModes()),
                           Ba,Bb),
            StdTriExp     (Ba,Bb),
            m_nodalPointsKey()
        {
            ASSERTL0(m_base[0]->GetNumModes() == m_base[1]->GetNumModes(),
                     "Nodal basis initiated with different orders in the a "
                     "and b directions");   
            int nummodes =  Ba.GetNumModes();
            m_nodalPointsKey = MemoryManager<LibUtilities::PointsKey>::
                AllocateSharedPtr(nummodes,Ntype);
        }

        StdNodalTriExp::StdNodalTriExp(const StdNodalTriExp &T):
            StdExpansion(T),
            StdExpansion2D(T),
            StdTriExp(T),
            m_nodalPointsKey(T.m_nodalPointsKey)
        {
        }

        StdNodalTriExp::~StdNodalTriExp()
        { 
        }
        
        
        //-------------------------------
        // Nodal basis specific routines
        //-------------------------------
        
        void StdNodalTriExp::NodalToModal()
        {
            NodalToModal(m_coeffs,m_coeffs); 
        }

        void StdNodalTriExp::NodalToModal(
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            StdMatrixKey   Nkey(eInvNBasisTrans, DetShapeType(), *this,
                                NullConstFactorMap, NullVarCoeffMap,
                                m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  inv_vdm = GetStdMatrix(Nkey);

            NekVector<NekDouble> nodal(m_ncoeffs,inarray,eWrapper);
            NekVector<NekDouble> modal(m_ncoeffs,outarray,eWrapper);
            modal = (*inv_vdm) * nodal;
        }

        void StdNodalTriExp::NodalToModalTranspose()
        {
            NodalToModalTranspose(m_coeffs,m_coeffs); 
        }

        // Operate with transpose of NodalToModal transformation
        void StdNodalTriExp::NodalToModalTranspose(
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            StdMatrixKey   Nkey(eInvNBasisTrans, DetShapeType(), *this,
                                NullConstFactorMap, NullVarCoeffMap,
                                m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  inv_vdm = GetStdMatrix(Nkey);

            NekVector<NekDouble> nodal(m_ncoeffs,inarray,eCopy);
            NekVector<NekDouble> modal(m_ncoeffs,outarray,eWrapper);
            modal = Transpose(*inv_vdm) * nodal;
        }

        void StdNodalTriExp::ModalToNodal()
        {
            ModalToNodal(m_coeffs,m_coeffs);
        }

        void StdNodalTriExp::ModalToNodal(
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            StdMatrixKey      Nkey(eNBasisTrans, DetShapeType(), *this,
                                    NullConstFactorMap, NullVarCoeffMap,
                                    m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  vdm = GetStdMatrix(Nkey);

            // Multiply out matrix
            NekVector<NekDouble> modal(m_ncoeffs,inarray,eWrapper);
            NekVector<NekDouble> nodal(m_ncoeffs,outarray,eWrapper);
            nodal = (*vdm)*modal;
        }

        void StdNodalTriExp::GetNodalPoints(
            Array<OneD, const NekDouble> &x, 
            Array<OneD, const NekDouble> &y)
        {
            LibUtilities::PointsManager()[*m_nodalPointsKey]->GetPoints(x,y);
        }

        DNekMatSharedPtr StdNodalTriExp::GenNBasisTransMatrix()
        {
            int             i,j;
            Array<OneD, const NekDouble>  r, s; 
            Array<OneD, NekDouble> c(2);
            DNekMatSharedPtr Mat;

            Mat = MemoryManager<DNekMat>::AllocateSharedPtr(
                m_ncoeffs, m_ncoeffs);
            GetNodalPoints(r,s);

            //Store the values of m_phys in a temporary array
            int nqtot = GetTotPoints();
            Array<OneD,NekDouble> tmp_phys(nqtot);
            Vmath::Vcopy(nqtot,m_phys,1,tmp_phys,1);

            for(i = 0; i < m_ncoeffs; ++i)
            {
                // fill physical space with mode i
                StdTriExp::v_FillMode(i,m_phys);

                // interpolate mode i to the Nodal points 'j' and
                // store in outarray
                for(j = 0; j < m_ncoeffs; ++j)
                {
                    c[0] = r[j];
                    c[1] = s[j];
                    (*Mat)(j,i) = StdTriExp::v_PhysEvaluate(c);
                }
            }
            // Restore the original values of m_phys
            Vmath::Vcopy(nqtot,tmp_phys,1,m_phys,1);
            return Mat;
        }


        //---------------------------------------
        // Transforms
        //---------------------------------------
        
        void StdNodalTriExp::v_BwdTrans(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            v_BwdTrans_SumFac(inarray,outarray);
        }

        void StdNodalTriExp::v_BwdTrans_SumFac(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            Array<OneD, NekDouble> tmp(m_ncoeffs);
            NodalToModal(inarray,tmp);
            StdTriExp::v_BwdTrans_SumFac(tmp,outarray);
        }

        void StdNodalTriExp::v_FwdTrans(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            v_IProductWRTBase(inarray,outarray);
            
            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass, DetShapeType(), *this,
                                      NullConstFactorMap, NullVarCoeffMap,
                                      m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);

            // copy inarray in case inarray == outarray
            NekVector<NekDouble> in(m_ncoeffs,outarray,eCopy);
            NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);
            
            out = (*matsys)*in;
        }


        //---------------------------------------
        // Inner product functions
        //---------------------------------------
        
        void StdNodalTriExp::v_IProductWRTBase(
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            v_IProductWRTBase_SumFac(inarray,outarray);
        }
        
        void StdNodalTriExp::v_IProductWRTBase_SumFac(
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            StdTriExp::v_IProductWRTBase_SumFac(inarray,outarray);
            NodalToModalTranspose(outarray,outarray);    
        }

        void StdNodalTriExp::v_IProductWRTDerivBase(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            v_IProductWRTDerivBase_SumFac(dir,inarray,outarray);
        }
        
        void StdNodalTriExp::v_IProductWRTDerivBase_SumFac(
            const int                           dir, 
            const Array<OneD, const NekDouble>& inarray, 
                  Array<OneD,       NekDouble>& outarray)
        {
            StdTriExp::v_IProductWRTDerivBase_SumFac(dir,inarray,outarray);
            NodalToModalTranspose(outarray,outarray);
        }
        
        //---------------------------------------
        // Evaluation functions
        //---------------------------------------
        
        void StdNodalTriExp::v_FillMode(
            const int               mode, 
            Array<OneD, NekDouble> &outarray)
        {
            ASSERTL2(mode >= m_ncoeffs, 
                "calling argument mode is larger than total expansion order");

            Vmath::Zero(m_ncoeffs, outarray, 1);
            outarray[mode] = 1.0;
            v_BwdTrans(outarray,outarray);
        }


        //---------------------------
        // Helper functions
        //---------------------------

        int StdNodalTriExp::v_NumBndryCoeffs() const
        {
            return 3 + (GetBasisNumModes(0)-2) + 2*(GetBasisNumModes(1)-2);
        } 

        void StdNodalTriExp::v_WriteToFile(
            std::ofstream &outfile, 
            OutputFormat   format, 
            const bool     dumpVar, 
            std::string    var)
        { 
            if(format==eTecplot)
            {
                int i,j;
                int nquad0 = m_base[0]->GetNumPoints();
                int nquad1 = m_base[1]->GetNumPoints();
                Array<OneD, const NekDouble> z0 = m_base[0]->GetZ();
                Array<OneD, const NekDouble> z1 = m_base[1]->GetZ();
                                                
                if(dumpVar)
                {
                    outfile << "Variables = z1,  z2"; 
                    outfile << ", "<< var << std::endl << std::endl;
                }
                
                outfile << "Zone, I=" << nquad0 << ", J=" << 
                    nquad1 <<", F=Point" << std::endl;

                for(j = 0; j < nquad1; ++j)
                {
                    for(i = 0; i < nquad0; ++i)
                    {
                        outfile << 0.5*(1+z0[i])*(1.0-z1[j])-1 <<  " " << 
                            z1[j] << " " << m_phys[j*nquad0+i] << std::endl;
                    }
                }
            }
            else if(format==eGmsh)
            {   
                if(dumpVar)
                {
                    outfile<<"View.MaxRecursionLevel = 4;"<<endl;
                    outfile<<"View.TargetError = 0.00;"<<endl;
                    outfile<<"View.AdaptVisualizationGrid = 1;"<<endl;
                    outfile<<"View \" \" {"<<endl;
                }

                outfile<<"ST("<<endl;                
                // write the coordinates of the vertices of the triangle
                outfile<<"-1.0, -1.0, 0.0,"<<endl;
                outfile<<" 1.0, -1.0, 0.0,"<<endl;
                outfile<<"-1.0,  1.0, 0.0" <<endl;
                outfile<<")"<<endl;

                // calculate the coefficients (monomial format)
                int i,j;

                Array<OneD,NekDouble> xi1(GetNcoeffs());
                Array<OneD,NekDouble> xi2(GetNcoeffs());
                GetNodalPoints(xi1,xi2);
                
                Array<OneD,NekDouble> x(GetNcoeffs());
                Array<OneD,NekDouble> y(GetNcoeffs());
                
                for(i=0;i<GetNcoeffs();i++)
                {
                    x[i] = 0.5*(1.0+xi1[i]);
                    y[i] = 0.5*(1.0+xi2[i]);
                }

                int cnt  = 0;
                int cnt2 = 0;
                int maxnummodes = max(m_base[0]->GetNumModes(),
                                      m_base[1]->GetNumModes());
                int nDumpCoeffs = maxnummodes*maxnummodes;
                Array<TwoD, int> dumpExponentMap(nDumpCoeffs,3,0);
                Array<OneD, int> indexMap(GetNcoeffs(),0);
                Array<TwoD, int> exponentMap(GetNcoeffs(),3,0);
                for(i = 0; i < maxnummodes; i++)
                {
                    for(j = 0; j < maxnummodes; j++)
                    {
                        if(j<maxnummodes-i)
                        {
                            exponentMap[cnt][0] = j;
                            exponentMap[cnt][1] = i;
                            indexMap[cnt++]  = cnt2;
                        }

                        dumpExponentMap[cnt2][0]   = j;
                        dumpExponentMap[cnt2++][1] = i;
                    }            
                }

                NekMatrix<NekDouble> vdm(GetNcoeffs(),GetNcoeffs());
                for(i = 0 ; i < GetNcoeffs(); i++)
                {
                    for(j = 0 ; j < GetNcoeffs(); j++)
                    {
                        vdm(i,j) = pow(x[i],exponentMap[j][0])*
                            pow(y[i],exponentMap[j][1]);
                    }
                } 

                vdm.Invert();  

                NekVector<NekDouble> in(GetNcoeffs(),m_coeffs,eWrapper);
                NekVector<NekDouble> out(GetNcoeffs());
                out = vdm*in;

                Array<OneD,NekDouble> dumpOut(nDumpCoeffs,0.0);
                for(i = 0 ; i < GetNcoeffs(); i++)
                {
                    dumpOut[ indexMap[i]  ] = out[i];
                }

                //write the coefficients
                outfile<<"{";
                for(i = 0; i < nDumpCoeffs; i++)
                {
                    outfile<<dumpOut[i];
                    if(i < nDumpCoeffs - 1)
                    {
                        outfile<<", ";
                    }
                }
                outfile<<"};"<<endl;
              
                if(dumpVar)
                {   
                    outfile<<"INTERPOLATION_SCHEME"<<endl;
                    outfile<<"{"<<endl;
                    for(i=0; i < nDumpCoeffs; i++)
                    {
                        outfile<<"{";
                        for(j = 0; j < nDumpCoeffs; j++)
                        {
                            if(i==j)
                            {
                                outfile<<"1.00";
                            }
                            else
                            {
                                outfile<<"0.00";
                            }
                            if(j < nDumpCoeffs - 1)
                            {
                                outfile<<", ";
                            }
                        }
                        if(i < nDumpCoeffs - 1)
                        {
                            outfile<<"},"<<endl;
                        }
                        else
                        {
                            outfile<<"}"<<endl<<"}"<<endl;
                        }
                    }
                    
                    outfile<<"{"<<endl;
                    for(i=0; i < nDumpCoeffs; i++)
                    {
                        outfile<<"{";
                        for(j = 0; j < 3; j++)
                        {
                            outfile<<dumpExponentMap[i][j];
                            if(j < 2)
                            {
                                outfile<<", ";
                            }
                        }
                        if(i < nDumpCoeffs  - 1)
                        {
                            outfile<<"},"<<endl;
                        }
                        else
                        {
                            outfile<<"}"<<endl<<"};"<<endl;
                        }
                    }
                    outfile<<"};"<<endl;
                }    
            }
            else
            {
                ASSERTL0(false, "Output routine not implemented for "
                         "requested type of output");
            }
        }      
        
        
        //--------------------------
        // Mappings
        //--------------------------

        void StdNodalTriExp::v_GetEdgeToElementMap(
            const int                  eid,
            const Orientation      edgeOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD,          int> &signarray)
        {
            ASSERTL0(eid >= 0 && eid <= 2,
                     "Local Edge ID must be between 0 and 2"); 
            
            const int nEdgeCoeffs = GetEdgeNcoeffs(eid);
            
            if (maparray.num_elements() != nEdgeCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeCoeffs);
            }

            if (signarray.num_elements() != nEdgeCoeffs)
            {
                signarray = Array<OneD, int>(nEdgeCoeffs,1);
            }
            else
            {
                fill(signarray.get(), signarray.get()+nEdgeCoeffs, 1);
            }
            
            maparray[0] = eid;
            maparray[1] = eid == 2 ? 0 : eid+1;
            for (int i = 2; i < nEdgeCoeffs; i++)
            {
                maparray[i] = eid*(nEdgeCoeffs-2)+1+i; 
            }  

            if (edgeOrient == eBackwards)
            {
                reverse(maparray.get(), maparray.get()+nEdgeCoeffs);
            }
        }

        int StdNodalTriExp::v_GetVertexMap(const int localVertexId)
        {
            ASSERTL0(localVertexId >= 0 && localVertexId <= 2,
                     "Local Vertex ID must be between 0 and 2");                
            return localVertexId;
        }

        void StdNodalTriExp::v_GetEdgeInteriorMap(
            const int                  eid,
            const Orientation      edgeOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD,          int> &signarray)
        {
            ASSERTL0(eid >= 0 && eid <= 2,
                     "Local Edge ID must be between 0 and 2"); 
            
            const int nEdgeIntCoeffs = GetEdgeNcoeffs(eid)-2;
            
            if (maparray.num_elements() != nEdgeIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeIntCoeffs);
            }

            if (signarray.num_elements() != nEdgeIntCoeffs)
            {
                signarray = Array<OneD, int>(nEdgeIntCoeffs,1);
            }
            else
            {
                fill(signarray.get(), signarray.get()+nEdgeIntCoeffs, 1);
            }
            
            for (int i = 0; i < nEdgeIntCoeffs; i++)
            {
                maparray[i] = eid*nEdgeIntCoeffs+3+i; 
            }  

            if (edgeOrient == eBackwards)
            {
                reverse(maparray.get(), maparray.get()+nEdgeIntCoeffs);
            }
        }

        void StdNodalTriExp::v_GetInteriorMap(
            Array<OneD, unsigned int>& outarray)
        {
            unsigned int i;
            if (outarray.num_elements() != GetNcoeffs()-NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(
                    GetNcoeffs()-NumBndryCoeffs());
            }

            for (i = NumBndryCoeffs(); i < GetNcoeffs(); i++)
            {
                outarray[i-NumBndryCoeffs()] = i;
            }
        }
        
        void StdNodalTriExp::v_GetBoundaryMap(
            Array<OneD, unsigned int>& outarray)
        {
            unsigned int i;
            if (outarray.num_elements()!=NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(NumBndryCoeffs());
            }
            
            for (i = 0; i < NumBndryCoeffs(); i++)
            {
                outarray[i] = i;
            }
        }


        //---------------------------------------
        // Wrapper functions
        //---------------------------------------
        
        DNekMatSharedPtr StdNodalTriExp::v_GenMatrix(const StdMatrixKey &mkey)
        {
            DNekMatSharedPtr Mat;
            
            switch(mkey.GetMatrixType())
            {
                case eNBasisTrans:
                    Mat = GenNBasisTransMatrix();
                    break;
                default:
                    Mat = StdExpansion::CreateGeneralMatrix(mkey);
                    break;
            }
            
            return Mat;
        }
        
        DNekMatSharedPtr StdNodalTriExp::v_CreateStdMatrix(
            const StdMatrixKey &mkey)
        {
            return StdNodalTriExp::v_GenMatrix(mkey);
        }


        //---------------------------------------
        // Operator evaluation functions
        //---------------------------------------

        void StdNodalTriExp::v_MassMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {              
            StdExpansion::MassMatrixOp_MatFree(inarray,outarray,mkey);
        }

        void StdNodalTriExp::v_LaplacianMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {                
            StdExpansion::LaplacianMatrixOp_MatFree_GenericImpl(
                inarray,outarray,mkey);
        }

        void StdNodalTriExp::v_LaplacianMatrixOp(
            const int                           k1,
            const int                           k2,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
            
        {           
            StdExpansion::LaplacianMatrixOp_MatFree(
                k1,k2,inarray,outarray,mkey);
        }

        void StdNodalTriExp::v_WeakDerivMatrixOp(
            const int                           i,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            StdExpansion::WeakDerivMatrixOp_MatFree(i,inarray,outarray,mkey);
        }
        
        void StdNodalTriExp::v_HelmholtzMatrixOp(
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
            const StdMatrixKey                 &mkey)
        {
            StdExpansion::HelmholtzMatrixOp_MatFree_GenericImpl(
                inarray,outarray,mkey);
        }       
        

        //---------------------------------------
        // Private helper functions
        //---------------------------------------

    } // end of namespace
} // end of namespace

