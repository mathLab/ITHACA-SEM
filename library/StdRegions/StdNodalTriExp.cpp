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

namespace Nektar
{
    namespace StdRegions
    {
        StdNodalTriExp::StdNodalTriExp(void):
            StdTriExp(),
            m_nodalPointsKey()
        {
        }

        StdNodalTriExp::StdNodalTriExp(const LibUtilities::BasisKey &Ba, 
                                       const LibUtilities::BasisKey &Bb, 
                                       LibUtilities::PointsType Ntype):
            StdExpansion  (StdTriData::getNumberOfCoefficients(Ba.GetNumModes(),(Bb.GetNumModes())),2,Ba,Bb),
            StdExpansion2D(StdTriData::getNumberOfCoefficients(Ba.GetNumModes(),(Bb.GetNumModes())),Ba,Bb),
            StdTriExp(Ba,Bb),
            m_nodalPointsKey()
        {
            ASSERTL0(m_base[0]->GetNumModes() == m_base[1]->GetNumModes(),
                     "Nodal basis initiated with different orders in the a "
                     "and b directions");   
            int nummodes =  Ba.GetNumModes();
            m_nodalPointsKey = MemoryManager<LibUtilities::PointsKey>::AllocateSharedPtr(nummodes,Ntype);
        }

        StdNodalTriExp::StdNodalTriExp(const StdNodalTriExp &T):
            StdExpansion(T),
            StdExpansion2D(T),
            StdTriExp(T),
            m_nodalPointsKey(T.m_nodalPointsKey)
        {
        }

        // Destructor
        StdNodalTriExp::~StdNodalTriExp()
        { 
        }
        
//        DNekMatSharedPtr StdNodalTriExp::GenMatrix(const StdMatrixKey &mkey)
//        {
//            DNekMatSharedPtr Mat;
//
//            switch(mkey.GetMatrixType())
//            {
//            case eNBasisTrans:
//                Mat =  GenNBasisTransMatrix();
//                break;
//            default:
//                Mat = StdExpansion::CreateGeneralMatrix(mkey);
//                break;
//            }
//
//            return Mat;
//        }

        DNekMatSharedPtr StdNodalTriExp::GenNBasisTransMatrix()
        {
            int             i,j;
            Array<OneD, const NekDouble>  r, s; 
            Array<OneD, NekDouble> c(2);
            DNekMatSharedPtr Mat;

            Mat = MemoryManager<DNekMat>::AllocateSharedPtr(m_ncoeffs,m_ncoeffs);
            GetNodalPoints(r,s);

            //Store the values of m_phys in a temporary array
            int nqtot = GetTotPoints();
            Array<OneD,NekDouble> tmp_phys(nqtot);
            Vmath::Vcopy(nqtot,m_phys,1,tmp_phys,1);

            for(i = 0; i < m_ncoeffs; ++i)
            {
                // fill physical space with mode i
                StdTriExp::FillMode(i,m_phys);

                // interpolate mode i to the Nodal points 'j' and
                // store in outarray
                for(j = 0; j < m_ncoeffs; ++j)
                {
                    c[0] = r[j];
                    c[1] = s[j];
                    (*Mat)(j,i) = StdTriExp::PhysEvaluate(c);
                }
            }
            // Restore the original values of m_phys
            Vmath::Vcopy(nqtot,tmp_phys,1,m_phys,1);
            return Mat;
        }

        void StdNodalTriExp::NodalToModal()
        {
            NodalToModal(m_coeffs,m_coeffs); 
        }

        void StdNodalTriExp::NodalToModal(const Array<OneD, const NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray)
        {
            StdMatrixKey   Nkey(eInvNBasisTrans,DetExpansionType(),*this,m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  inv_vdm = GetStdMatrix(Nkey);

            NekVector<const NekDouble> nodal(m_ncoeffs,inarray,eWrapper);
            NekVector<NekDouble> modal(m_ncoeffs,outarray,eWrapper);
            modal = (*inv_vdm) * nodal;
        }


        void StdNodalTriExp::NodalToModalTranspose()
        {
            NodalToModalTranspose(m_coeffs,m_coeffs); 
        }

        // Operate with transpose of NodalToModal transformation
        void StdNodalTriExp::NodalToModalTranspose(const Array<OneD, const NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray)
        {            
            StdMatrixKey   Nkey(eInvNBasisTrans,DetExpansionType(),*this,m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  inv_vdm = GetStdMatrix(Nkey);

            NekVector<const NekDouble> nodal(m_ncoeffs,inarray,eCopy);
            NekVector<NekDouble> modal(m_ncoeffs,outarray,eWrapper);
            modal = Transpose(*inv_vdm) * nodal;
        }


        void StdNodalTriExp::ModalToNodal()
        {
            ModalToNodal(m_coeffs,m_coeffs);
        }

        void StdNodalTriExp::ModalToNodal(const Array<OneD, const NekDouble>& inarray, 
                Array<OneD, NekDouble> &outarray)
        {
            StdMatrixKey      Nkey(eNBasisTrans,DetExpansionType(),*this,m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  vdm = GetStdMatrix(Nkey);

            // Multiply out matrix
            NekVector<const NekDouble> modal(m_ncoeffs,inarray,eWrapper);
            NekVector<NekDouble> nodal(m_ncoeffs,outarray,eWrapper);
            nodal = (*vdm)*modal;
        }

        void StdNodalTriExp::FillMode(const int mode, 
            Array<OneD, NekDouble> &outarray)
        {

            ASSERTL2(mode >= m_ncoeffs, 
                "calling argument mode is larger than total expansion order");

            Vmath::Zero(m_ncoeffs, outarray, 1);
            outarray[mode] = 1.0;
            BwdTrans(outarray,outarray);
        }

        ///////////////////////////////
        /// Evaluation Methods
        ///////////////////////////////

        void StdNodalTriExp::FwdTrans(const Array<OneD, const NekDouble>& inarray,
            Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(inarray,outarray);
            
            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetExpansionType(),*this, 
                                    m_nodalPointsKey->GetPointsType());
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);

            // copy inarray in case inarray == outarray
            NekVector<const NekDouble> in(m_ncoeffs,outarray,eCopy);
            NekVector<NekDouble> out(m_ncoeffs,outarray,eWrapper);
            
            out = (*matsys)*in;
        }

        void StdNodalTriExp::GetBoundaryMap(Array<OneD, unsigned int>& outarray)
        {
            unsigned int i;
            if(outarray.num_elements()!=NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(NumBndryCoeffs());
            }
            
            for(i = 0; i < NumBndryCoeffs(); i++)
            {
                outarray[i] = i;
            }
        }

        void StdNodalTriExp::GetInteriorMap(Array<OneD, unsigned int>& outarray)
        {
            unsigned int i;
            if(outarray.num_elements()!=GetNcoeffs()-NumBndryCoeffs())
            {
                outarray = Array<OneD, unsigned int>(GetNcoeffs()-NumBndryCoeffs());
            }

            for(i = NumBndryCoeffs(); i < GetNcoeffs(); i++)
            {
                outarray[i-NumBndryCoeffs()] = i;
            }
        }
 
        void StdNodalTriExp::GetEdgeInteriorMap(const int eid, const EdgeOrientation edgeOrient,
                                                Array<OneD, unsigned int> &maparray,
                                                Array<OneD, int> &signarray)
        {
            ASSERTL0((eid>=0)&&(eid<=2),"Local Edge ID must be between 0 and 2"); 
            const int nEdgeIntCoeffs = GetEdgeNcoeffs(eid)-2;
            
            if(maparray.num_elements() != nEdgeIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeIntCoeffs);
            }

            if(signarray.num_elements() != nEdgeIntCoeffs)
            {
                signarray = Array<OneD, int>(nEdgeIntCoeffs,1);
            }
            else
            {
                fill( signarray.get() , signarray.get()+nEdgeIntCoeffs, 1 );
            }
            
            for(int i = 0; i < nEdgeIntCoeffs; i++)
            {
                maparray[i] = eid*nEdgeIntCoeffs+3+i; 
            }  

            if(edgeOrient == eBackwards)
            {
                reverse( maparray.get() , maparray.get()+nEdgeIntCoeffs );
            }
                          
        }

        void StdNodalTriExp::GetEdgeToElementMap(const int eid, const EdgeOrientation edgeOrient,
                                             Array<OneD, unsigned int> &maparray,
                                             Array<OneD, int> &signarray)
        {
            ASSERTL0((eid>=0)&&(eid<=2),"Local Edge ID must be between 0 and 2"); 
            const int nEdgeCoeffs = GetEdgeNcoeffs(eid);
            
            if(maparray.num_elements() != nEdgeCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeCoeffs);
            }

            if(signarray.num_elements() != nEdgeCoeffs)
            {
                signarray = Array<OneD, int>(nEdgeCoeffs,1);
            }
            else
            {
                fill( signarray.get() , signarray.get()+nEdgeCoeffs, 1 );
            }
            
            maparray[0] = eid;
            maparray[1] = (eid==2) ? 0 : eid+1;
            for(int i = 2; i < nEdgeCoeffs; i++)
            {
                maparray[i] = eid*(nEdgeCoeffs-2)+1+i; 
            }  

            if(edgeOrient == eBackwards)
            {
                reverse( maparray.get() , maparray.get()+nEdgeCoeffs );
            }
        }

        void StdNodalTriExp::WriteToFile(std::ofstream &outfile, OutputFormat format, const bool dumpVar, std::string var)
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
                int maxnummodes = max(m_base[0]->GetNumModes(),m_base[1]->GetNumModes());
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
                        vdm(i,j) = pow(x[i],exponentMap[j][0])*pow(y[i],exponentMap[j][1]);
                    }
                } 

                vdm.Invert();  

                NekVector<const NekDouble> in(GetNcoeffs(),m_coeffs,eWrapper);
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
                ASSERTL0(false, "Output routine not implemented for requested type of output");
            }
        }      


    } // end of namespace
} // end of namespace

/** 
* $Log: StdNodalTriExp.cpp,v $
* Revision 1.32  2009/12/18 00:26:54  bnelson
* Fixed visual studio compiler warning.
*
* Revision 1.31  2009/04/27 21:32:45  sherwin
* Updated WriteToField method
*
* Revision 1.30  2008/11/05 16:08:15  pvos
* Added elemental optimisation functionality
*
* Revision 1.29  2008/09/17 13:46:06  pvos
* Added LocalToGlobalC0ContMap for 3D expansions
*
* Revision 1.28  2008/07/19 21:12:54  sherwin
* Removed MapTo function and made orientation convention anticlockwise in UDG routines
*
* Revision 1.27  2008/07/04 10:18:40  pvos
* Some updates
*
* Revision 1.26  2008/06/05 15:06:06  pvos
* Added documentation
*
* Revision 1.25  2008/05/30 00:33:49  delisi
* Renamed StdRegions::ShapeType to StdRegions::ExpansionType.
*
* Revision 1.24  2008/05/29 21:36:25  pvos
* Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
*
* Revision 1.23  2008/05/07 16:04:57  pvos
* Mapping + Manager updates
*
* Revision 1.22  2008/04/06 06:04:15  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.21  2008/04/02 22:18:10  pvos
* Update for 2D local to global mapping
*
* Revision 1.20  2008/03/18 14:15:45  pvos
* Update for nodal triangular helmholtz solver
*
* Revision 1.19  2007/12/17 13:03:51  sherwin
* Modified StdMatrixKey to contain a list of constants and GenMatrix to take a StdMatrixKey
*
* Revision 1.18  2007/08/11 23:42:26  sherwin
* A few changes
*
* Revision 1.17  2007/07/27 00:22:54  bnelson
* Memory manager now accepts non-const parameters to the allocate methods.
*
* Revision 1.16  2007/07/20 02:16:54  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.15  2007/07/12 12:55:16  sherwin
* Simplified Matrix Generation
*
* Revision 1.14  2007/07/11 13:35:17  kirby
* *** empty log message ***
*
* Revision 1.13  2007/07/11 06:35:24  sherwin
* Update after manager reshuffle
*
* Revision 1.12  2007/07/10 21:05:17  kirby
* even more fixes
*
* Revision 1.11  2007/07/09 15:19:15  sherwin
* Introduced an InvMassMatrix and replaced the StdLinSysManager call with a StdMatrixManager call to the inverse matrix
*
* Revision 1.10  2007/05/15 05:18:23  bnelson
* Updated to use the new Array object.
*
* Revision 1.9  2007/04/26 15:00:17  sherwin
* SJS compiling working version using SHaredArrays
*
* Revision 1.8  2007/04/18 09:44:01  sherwin
* Working version for StdNodalTri. Removed lapack.cpp from compile.
*
* Revision 1.7  2007/04/10 14:00:45  sherwin
* Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions
*
* Revision 1.6  2007/01/17 16:05:40  pvos
* updated doxygen documentation
*
* Revision 1.5  2007/01/15 21:13:46  sherwin
* Nodal stuff correction and added Field Classes
*
* Revision 1.4  2006/12/10 19:00:54  sherwin
* Modifications to handle nodal expansions
*
* Revision 1.3  2006/06/06 15:25:21  jfrazier
* Removed unreferenced variables and replaced ASSERTL0(false, ....) with
* NEKERROR.
*
* Revision 1.2  2006/06/01 14:46:16  kirby
* *** empty log message ***
*
* Revision 1.1  2006/05/04 18:58:32  kirby
* *** empty log message ***
*
* Revision 1.19  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.18  2006/04/01 21:59:27  sherwin
* Sorted new definition of ASSERT
*
* Revision 1.17  2006/03/21 09:21:32  sherwin
* Introduced NekMemoryManager
*
* Revision 1.16  2006/03/01 08:25:04  sherwin
*
* First compiling version of StdRegions
*
**/ 

