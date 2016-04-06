///////////////////////////////////////////////////////////////////////////////
//
// File StdPyrExp.cpp
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
// Description: pyramadic routines built upon StdExpansion3D
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/StdPyrExp.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <iomanip>

using namespace std;

namespace Nektar
{
    namespace StdRegions
    {
        StdPyrExp::StdPyrExp() // Deafult construct of standard expansion directly called. 
        {
        }
        
        StdPyrExp::StdPyrExp(const LibUtilities::BasisKey &Ba,
                             const LibUtilities::BasisKey &Bb,
                             const LibUtilities::BasisKey &Bc) 
            : StdExpansion  (LibUtilities::StdPyrData::getNumberOfCoefficients(
                                 Ba.GetNumModes(),
                                 Bb.GetNumModes(),
                                 Bc.GetNumModes()),
                             3, Ba, Bb, Bc),
              StdExpansion3D(LibUtilities::StdPyrData::getNumberOfCoefficients(
                                 Ba.GetNumModes(),
                                 Bb.GetNumModes(),
                                 Bc.GetNumModes()),
                             Ba, Bb, Bc)
        {
            if (Ba.GetNumModes() > Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'a' direction is higher "
                         "than order in 'c' direction");
            }
            if (Bb.GetNumModes() > Bc.GetNumModes())
            {
                ASSERTL0(false, "order in 'b' direction is higher "
                         "than order in 'c' direction");
            }

            // Set up mode mapping which takes 0\leq i\leq N_coeffs -> (p,q,r)
            // of the 3D tensor product
            const int P = Ba.GetNumModes() - 1;
            const int Q = Bb.GetNumModes() - 1;
            const int R = Bc.GetNumModes() - 1;
            int cnt = 0;

            // Vertices
            m_map[Mode(0, 0, 0, 0)] = cnt++;
            m_map[Mode(1, 0, 0, 0)] = cnt++;
            m_map[Mode(1, 1, 0, 0)] = cnt++;
            m_map[Mode(0, 1, 0, 0)] = cnt++;
            m_map[Mode(0, 0, 1, 1)] = cnt++;

            // Edge 0
            for (int i = 2; i <= P; ++i)
            {
                m_map[Mode(i, 0, 0, GetTetMode(i, 0, 0))] = cnt++;
            }

            // Edge 1
            for (int i = 2; i <= Q; ++i)
            {
                m_map[Mode(1, i, 0, GetTetMode(0, i, 0))] = cnt++;
            }

            // Edge 2
            for (int i = 2; i <= P; ++i)
            {
                m_map[Mode(i, 1, 0, GetTetMode(i, 0, 0))] = cnt++;
            }

            // Edge 3
            for (int i = 2; i <= Q; ++i)
            {
                m_map[Mode(0, i, 0, GetTetMode(0, i, 0))] = cnt++;
            }

            // Edge 4
            for (int i = 2; i <= R; ++i)
            {
                m_map[Mode(0, 0, i, i)] = cnt++;
            }

            // Edge 5
            for (int i = 2; i <= R; ++i)
            {
                m_map[Mode(1, 0, i, i)] = cnt++;
            }

            // Edge 6
            for (int i = 2; i <= R; ++i)
            {
                m_map[Mode(1, 1, i, i)] = cnt++;
            }

            // Edge 7
            for (int i = 2; i <= R; ++i)
            {
                m_map[Mode(0, 1, i, i)] = cnt++;
            }

            // Face 0 - TODO check this
            for (int j = 2; j <= Q; ++j)
            {
                for (int i = 2; i <= P; ++i)
                {
                    m_map[Mode(i, j, 0, GetTetMode((i-2+j-2) % (Q-1) + 2, 0, 0))] = cnt++;
                }
            }

            // Face 1
            for (int i = 2; i <= P; ++i)
            {
                for (int j = 1; j <= R-i; ++j)
                {
                    m_map[Mode(i, 0, j, GetTetMode(i, 0, j))] = cnt++;
                }
            }

            // Face 2
            for (int i = 2; i <= Q; ++i)
            {
                for (int j = 1; j <= R-i; ++j)
                {
                    m_map[Mode(1, i, j, GetTetMode(0, i, j))] = cnt++;
                }
            }

            // Face 3
            for (int i = 2; i <= P; ++i)
            {
                for (int j = 1; j <= R-i; ++j)
                {
                    m_map[Mode(i, 1, j, GetTetMode(i, 0, j))] = cnt++;
                }
            }

            // Face 4
            for (int i = 2; i <= Q; ++i)
            {
                for (int j = 1; j <= R-i; ++j)
                {
                    m_map[Mode(0, i, j, GetTetMode(0, i, j))] = cnt++;
                }
            }

            // Interior (tetrahedral modes)
            for (int i = 2; i <= P+1; ++i)
            {
                for (int j = 1; j <= Q-i+1; ++j)
                {
                    for (int k = 1; k <= R-i-j+1; ++k)
                    {
                        // need to go to j+1-th mode in the 'b' direction to
                        // select correct modified_a mode
                        m_map[Mode(i, j+1, k, GetTetMode(i-1, j, k))] = cnt++;
                    }
                }
            }

            ASSERTL0(m_map.size() == m_ncoeffs,
                     "Duplicate coefficient entries in map");

            map<Mode, unsigned int, cmpop>::iterator it;
            for (it = m_map.begin(); it != m_map.end(); ++it)
            {
                const int p  = it->first.get<0>();
                const int q  = it->first.get<1>();
                const int r  = it->first.get<2>();
                const int rp = it->first.get<3>();
                if (m_idxMap.find(p) == m_idxMap.end())
                {
                    m_idxMap[p] = map<int, map<int, pair<int, int> > >();
                }

                if (m_idxMap[p].find(q) == m_idxMap[p].end())
                {
                    m_idxMap[p][q] = map<int, pair<int, int> >();
                }

                if (m_idxMap[p][q].find(r) == m_idxMap[p][q].end())
                {
                    m_idxMap[p][q][r] = pair<int, int>(it->second, rp);
                }
            }
        }

        StdPyrExp::StdPyrExp(const StdPyrExp &T)
            : StdExpansion  (T),
              StdExpansion3D(T)
        {
        }


        // Destructor
        StdPyrExp::~StdPyrExp()
        {   
        } 


        //---------------------------------------
        // Differentiation/integration Methods
        //---------------------------------------
        
        /**
         * \brief Calculate the derivative of the physical points 
         *  
         * The derivative is evaluated at the nodal physical points.
         * Derivatives with respect to the local Cartesian coordinates.
         *  
         * \f$\begin{Bmatrix} \frac {\partial} {\partial \xi_1} \\ \frac
         * {\partial} {\partial \xi_2} \\ \frac {\partial} {\partial \xi_3}
         * \end{Bmatrix} = \begin{Bmatrix} \frac 2 {(1-\eta_3)} \frac \partial
         * {\partial \bar \eta_1} \\ \frac {\partial} {\partial \xi_2} \ \
         * \frac {(1 + \bar \eta_1)} {(1 - \eta_3)} \frac \partial {\partial
         * \bar \eta_1} + \frac {\partial} {\partial \eta_3} \end{Bmatrix}\f$
         */
        void StdPyrExp::v_PhysDeriv(
            const Array<OneD, const NekDouble> &u_physical,
                  Array<OneD,       NekDouble> &out_dxi1,
                  Array<OneD,       NekDouble> &out_dxi2,
                  Array<OneD,       NekDouble> &out_dxi3)
        {
            // PhysDerivative implementation based on Spen's book page 152.
            int    Qx = m_base[0]->GetNumPoints();
            int    Qy = m_base[1]->GetNumPoints();
            int    Qz = m_base[2]->GetNumPoints();

            Array<OneD, NekDouble> dEta_bar1(Qx*Qy*Qz,0.0);
            Array<OneD, NekDouble> dXi2     (Qx*Qy*Qz,0.0);
            Array<OneD, NekDouble> dEta3    (Qx*Qy*Qz,0.0);
            PhysTensorDeriv(u_physical, dEta_bar1, dXi2, dEta3);

            Array<OneD, const NekDouble> eta_x, eta_y, eta_z;
            eta_x = m_base[0]->GetZ();
            eta_y = m_base[1]->GetZ();
            eta_z = m_base[2]->GetZ();

            int i, j, k, n;

            for (k = 0, n = 0; k < Qz; ++k)
            {
                for (j = 0; j < Qy; ++j)
                {
                    for (i = 0; i < Qx; ++i, ++n)
                    {
                        if (out_dxi1.num_elements() > 0)
                            out_dxi1[n] = 2.0/(1.0 - eta_z[k]) * dEta_bar1[n];
                        if (out_dxi2.num_elements() > 0)
                            out_dxi2[n] = 2.0/(1.0 - eta_z[k]) * dXi2[n];
                        if (out_dxi3.num_elements() > 0)
                            out_dxi3[n] = (1.0+eta_x[i])/(1.0-eta_z[k])*dEta_bar1[n] +
                                (1.0+eta_y[j])/(1.0-eta_z[k])*dXi2[n] + dEta3[n];
                    } 
                }
            }
        }

        void StdPyrExp::v_PhysDeriv(const int dir,
                                    const Array<OneD, const NekDouble>& inarray,
                                          Array<OneD,       NekDouble>& outarray)
        {
            switch(dir)
            {
                case 0:
                {
                    v_PhysDeriv(inarray, outarray, NullNekDouble1DArray,
                                NullNekDouble1DArray);
                    break;
                }
                
                case 1:
                {
                    v_PhysDeriv(inarray, NullNekDouble1DArray, outarray,
                                NullNekDouble1DArray);
                    break;
                }
                
                case 2:
                {
                    v_PhysDeriv(inarray, NullNekDouble1DArray,
                                NullNekDouble1DArray, outarray);
                    break;
                }
                
                default:
                {
                    ASSERTL1(false,"input dir is out of range");
                }
                break;
            }
        }

        void StdPyrExp::v_StdPhysDeriv(const Array<OneD, const NekDouble> &inarray, 
                                             Array<OneD,       NekDouble> &out_d0,
                                             Array<OneD,       NekDouble> &out_d1,
                                             Array<OneD,       NekDouble> &out_d2)
        {
            StdPyrExp::v_PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }

        void StdPyrExp::v_StdPhysDeriv(const int dir,
                                       const Array<OneD, const NekDouble>& inarray,
                                             Array<OneD,       NekDouble>& outarray)
        {
            StdPyrExp::v_PhysDeriv(dir, inarray, outarray);
        }
        
        //---------------------------------------
        // Transforms
        //---------------------------------------
        
	/** 
         * \brief Backward transformation is evaluated at the quadrature
         * points.
         *
         * \f$ u^{\delta} (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{m(pqr)} \hat
         * u_{pqr} \phi_{pqr} (\xi_{1i}, \xi_{2j}, \xi_{3k})\f$
         * 
         * Backward transformation is three dimensional tensorial expansion
         *
         * \f$ u (\xi_{1i}, \xi_{2j}, \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_p^a
	 *  (\xi_{1i}) \lbrace { \sum_{q=0}^{Q_y} \psi_{q}^a (\xi_{2j})
	 *  \lbrace { \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c (\xi_{3k})
	 *  \rbrace} \rbrace}. \f$
	 *
         * And sumfactorizing step of the form is as:\ \ \f$ f_{pqr}
         * (\xi_{3k}) = \sum_{r=0}^{Q_z} \hat u_{pqr} \psi_{pqr}^c
         * (\xi_{3k}),\\ g_{p} (\xi_{2j}, \xi_{3k}) = \sum_{r=0}^{Q_y}
         * \psi_{p}^a (\xi_{2j}) f_{pqr} (\xi_{3k}),\\ u(\xi_{1i}, \xi_{2j},
         * \xi_{3k}) = \sum_{p=0}^{Q_x} \psi_{p}^a (\xi_{1i}) g_{p}
         * (\xi_{2j}, \xi_{3k}).  \f$
         **/
        void StdPyrExp::v_BwdTrans(const Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,       NekDouble> &outarray)
        {
            if (m_base[0]->Collocation() && 
                m_base[1]->Collocation() &&
                m_base[2]->Collocation())
            {
                Vmath::Vcopy(m_base[0]->GetNumPoints()*
                             m_base[1]->GetNumPoints()*
                             m_base[2]->GetNumPoints(),
                             inarray, 1, outarray, 1);
            }
            else
            {
                StdPyrExp::v_BwdTrans_SumFac(inarray,outarray);
            }
        }

        /**
         * Sum-factorisation implementation of the BwdTrans operation.
         */
        void StdPyrExp::v_BwdTrans_SumFac(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            Array<OneD, NekDouble> wsp;
            v_BwdTrans_SumFacKernel(m_base[0]->GetBdata(),
                                    m_base[1]->GetBdata(),
                                    m_base[2]->GetBdata(),
                                    inarray,outarray,wsp,
                                    true,true,true);
        }

        void StdPyrExp::v_BwdTrans_SumFacKernel(
            const Array<OneD, const NekDouble> &base0,
            const Array<OneD, const NekDouble> &base1,
            const Array<OneD, const NekDouble> &base2,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
                  Array<OneD,       NekDouble> &wsp,
            bool                                doCheckCollDir0,
            bool                                doCheckCollDir1,
            bool                                doCheckCollDir2)
        {
            const int Qx = m_base[0]->GetNumPoints();
            const int Qy = m_base[1]->GetNumPoints();
            const int Qz = m_base[2]->GetNumPoints();

            const NekDouble *bx = base0.get();
            const NekDouble *by = base1.get();
            const NekDouble *bz = base2.get();

            // Need to count coeffs for storage...
            map<int, map<int, map<int, pair<int, int> > > >::iterator it_p;
            map<int, map<int,          pair<int, int> > >  ::iterator it_q;
            map<int,                   pair<int, int> >    ::iterator it_r;

            int pqcnt = 0;
            for (it_p = m_idxMap.begin(); it_p != m_idxMap.end(); ++it_p)
            {
                for (it_q = it_p->second.begin(); it_q != it_p->second.end(); ++it_q)
                {
                    pqcnt++;
                }
            }

            Array<OneD, NekDouble> fpq(pqcnt);
            Array<OneD, NekDouble> fp (m_base[0]->GetNumModes());
            int i ,j, k, s = 0, cnt = 0, cnt2 = 0;

            for (k = 0; k < Qz; ++k)
            {
                NekDouble bz1 = bz[k+Qz];

                cnt = 0;
                for (it_p = m_idxMap.begin(); it_p != m_idxMap.end(); ++it_p)
                {
                    for (it_q = it_p->second.begin(); it_q != it_p->second.end(); ++it_q)
                    {
                        NekDouble sum = 0.0;
                        for (it_r = it_q->second.begin(); it_r != it_q->second.end(); ++it_r)
                        {
                            sum += inarray[it_r->second.first] * bz[k + Qz*it_r->second.second];
                        }
                        fpq[cnt++] = sum;
                    }
                }

                for (j = 0; j < Qy; ++j)
                {
                    NekDouble by0 = bz1*by[j];
                    NekDouble by1 = bz1*by[j+Qy];

                    cnt = cnt2 = 0;
                    for (it_p = m_idxMap.begin(); it_p != m_idxMap.end(); ++it_p)
                    {
                        NekDouble sum = 0.0;
                        for (it_q = it_p->second.begin(); it_q != it_p->second.end(); ++it_q)
                        {
                            sum += by[j + Qy*it_q->first] * fpq[cnt++];
                        }
                        fp[cnt2++] = sum;
                    }

                    for (i = 0; i < Qx; ++i, ++s)
                    {
                        cnt2 = 0;
                        NekDouble sum = 0.0;
                        for (it_p = m_idxMap.begin(); it_p != m_idxMap.end(); ++it_p)
                        {
                            sum += bx[i + Qx*it_p->first] * fp[cnt2++];
                        }
                        sum += inarray[4]*(by1*(bx[i] + bx[i+Qx]) + by0*bx[i+Qx]);
                        outarray[s] = sum;
                    }
                }
            }
        }

	/** \brief Forward transform from physical quadrature space
            stored in \a inarray and evaluate the expansion coefficients and
            store in \a outarray
            
            Inputs:\n
            
            - \a inarray: array of physical quadrature points to be transformed
            
            Outputs:\n
            
            - \a outarray: updated array of expansion coefficients. 
            
        */    
        void StdPyrExp::v_FwdTrans(const Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,       NekDouble> &outarray)
        {
            v_IProductWRTBase(inarray,outarray);

            // get Mass matrix inverse
            StdMatrixKey      masskey(eInvMass,DetShapeType(),*this);
            DNekMatSharedPtr  matsys = GetStdMatrix(masskey);

            // copy inarray in case inarray == outarray
            DNekVec in (m_ncoeffs, outarray);
            DNekVec out(m_ncoeffs, outarray, eWrapper);

            out = (*matsys)*in;
        }
        
        
        //---------------------------------------
        // Inner product functions
        //---------------------------------------

        /** \brief  Inner product of \a inarray over region with respect to the 
            expansion basis m_base[0]->GetBdata(),m_base[1]->GetBdata(), m_base[2]->GetBdata() and return in \a outarray 
            
            Wrapper call to StdPyrExp::IProductWRTBase
            
            Input:\n
            
            - \a inarray: array of function evaluated at the physical collocation points
            
            Output:\n
            
            - \a outarray: array of inner product with respect to each basis over region
            
        */
        void StdPyrExp::v_IProductWRTBase(
            const Array<OneD, const NekDouble> &inarray, 
                  Array<OneD,       NekDouble> &outarray)
        {
            if (m_base[0]->Collocation() &&
                m_base[1]->Collocation() &&
                m_base[2]->Collocation())
            {
                MultiplyByStdQuadratureMetric(inarray, outarray);
            }
            else
            {
                StdPyrExp::v_IProductWRTBase_SumFac(inarray,outarray);
            }
        }

        void StdPyrExp::v_IProductWRTBase_SumFac(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray,
            bool                                multiplybyweights)
        {
            Array<OneD, NekDouble> wsp;

            if(multiplybyweights)
            {
                Array<OneD, NekDouble> tmp(inarray.num_elements());

                v_MultiplyByStdQuadratureMetric(inarray, tmp);

                v_IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                               m_base[1]->GetBdata(),
                                               m_base[2]->GetBdata(),
                                               tmp,outarray,wsp,
                                               true,true,true);
            }
            else
            {
                v_IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                               m_base[1]->GetBdata(),
                                               m_base[2]->GetBdata(),
                                               inarray,outarray,wsp,
                                               true,true,true);
            }
        }

        void StdPyrExp::v_IProductWRTBase_SumFacKernel(
            const Array<OneD, const NekDouble> &base0,
            const Array<OneD, const NekDouble> &base1,
            const Array<OneD, const NekDouble> &base2,
            const Array<OneD, const NekDouble> &inarray,
                  Array<OneD,       NekDouble> &outarray,
                  Array<OneD,       NekDouble> &wsp,
            bool                                doCheckCollDir0,
            bool                                doCheckCollDir1,
            bool                                doCheckCollDir2)
        {
            int i, j, k, s;
            int Qx = m_base[0]->GetNumPoints();
            int Qy = m_base[1]->GetNumPoints();
            int Qz = m_base[2]->GetNumPoints();

            const NekDouble *bx = base0.get();
            const NekDouble *by = base1.get();
            const NekDouble *bz = base2.get();
            
            map<int, map<int, map<int, pair<int, int> > > >::iterator it_p;
            map<int, map<int,          pair<int, int> > >  ::iterator it_q;
            map<int,                   pair<int, int> >    ::iterator it_r;

            Array<OneD, NekDouble> f (Qy*Qz);
            Array<OneD, NekDouble> fb(Qz);

            for (it_p = m_idxMap.begin(); it_p != m_idxMap.end(); ++it_p)
            {
                const int p = it_p->first;
                s = 0;
                for (k = 0; k < Qz; ++k)
                {
                    for (j = 0; j < Qy; ++j)
                    {
                        NekDouble sum = 0.0;
                        for (i = 0; i < Qx; ++i, ++s)
                        {
                            sum += bx[i + Qx*p]*inarray[s];
                        }
                        f[j+Qy*k] = sum;
                    }
                }

                for (it_q = it_p->second.begin(); it_q != it_p->second.end(); ++it_q)
                {
                    const int q = it_q->first;

                    for (k = 0; k < Qz; ++k)
                    {
                        NekDouble sum = 0.0;
                        for (j = 0; j < Qy; ++j)
                        {
                            sum += by[j + Qy*q]*f[j+Qy*k];
                        }
                        fb[k] = sum;
                    }

                    for (it_r = it_q->second.begin(); it_r != it_q->second.end(); ++it_r)
                    {
                        const int rpqr = it_r->second.second;
                        NekDouble sum = 0.0;
                        for (k = 0; k < Qz; ++k)
                        {
                            sum += bz[k + Qz*rpqr]*fb[k];
                        }

                        outarray[it_r->second.first] = sum;
                    }
                }
            }

            // Correct for top mode
            s = 0;
            for (k = 0; k < Qz; ++k)
            {
                for (j = 0; j < Qy; ++j)
                {
                    for (i = 0; i < Qx; ++i, ++s)
                    {
                        outarray[4] += inarray[s] * bz[k+Qz]*(
                            bx[i+Qx]*by[j+Qy] + 
                            bx[i+Qx]*by[j   ] + 
                            bx[i   ]*by[j+Qy]);
                    }
                }
            }
        }

        void StdPyrExp::v_IProductWRTDerivBase(
            const int                           dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            StdPyrExp::v_IProductWRTDerivBase_SumFac(dir,inarray,outarray);
        }

        /**
         * @param   inarray     Function evaluated at physical collocation
         *                      points.
         * @param   outarray    Inner product with respect to each basis
         *                      function over the element.
         */
        void StdPyrExp::v_IProductWRTDerivBase_SumFac(
            const int dir,
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int i;
            int nquad0  = m_base[0]->GetNumPoints();
            int nquad1  = m_base[1]->GetNumPoints();
            int nquad2  = m_base[2]->GetNumPoints();
            int nqtot   = nquad0*nquad1*nquad2;

            Array<OneD, NekDouble> gfac0(nquad0);
            Array<OneD, NekDouble> gfac1(nquad1);
            Array<OneD, NekDouble> gfac2(nquad2);
            Array<OneD, NekDouble> tmp0 (nqtot);
            Array<OneD, NekDouble> wsp;

            const Array<OneD, const NekDouble>& z0 = m_base[0]->GetZ();
            const Array<OneD, const NekDouble>& z1 = m_base[1]->GetZ();
            const Array<OneD, const NekDouble>& z2 = m_base[2]->GetZ();

            // set up geometric factor: (1+z0)/2
            for(i = 0; i < nquad0; ++i)
            {
                gfac0[i] = 0.5*(1+z0[i]);
            }

            // set up geometric factor: (1+z1)/2
            for(i = 0; i < nquad1; ++i)
            {
                gfac1[i] = 0.5*(1+z1[i]);
            }

            // Set up geometric factor: 2/(1-z2)
            for(i = 0; i < nquad2; ++i)
            {
            	gfac2[i] = 2.0/(1-z2[i]);
            }

            // Derivative in first/second direction is always scaled as follows
            const int nq01 = nquad0*nquad1;
            for(i = 0; i < nquad2; ++i)
            {
                Vmath::Smul(nq01, gfac2[i],
                            &inarray[0] + i*nq01, 1,
                            &tmp0   [0] + i*nq01, 1);
            }

            MultiplyByStdQuadratureMetric(tmp0, tmp0);

            switch(dir)
            {
                case 0:
                {
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata (),
                                                 m_base[2]->GetBdata (),
                                                 tmp0, outarray, wsp,
                                                 false, true, true);
                    break;
                }
                case 1:
                {
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata (),
                                                 m_base[1]->GetDbdata(),
                                                 m_base[2]->GetBdata (),
                                                 tmp0, outarray, wsp,
                                                 true, false, true);
                    break;
                }
                case 2:
                {
                    Array<OneD, NekDouble> tmp3(m_ncoeffs);
                    Array<OneD, NekDouble> tmp4(m_ncoeffs);

                    // Scale eta_1 derivative by gfac0
                    for(i = 0; i < nquad1*nquad2; ++i)
                    {
                        Vmath::Vmul(nquad0, tmp0 .get() + i*nquad0, 1,
                                            gfac0.get(),            1,
                                            tmp0 .get() + i*nquad0, 1);
                    }
                    IProductWRTBase_SumFacKernel(m_base[0]->GetDbdata(),
                                                 m_base[1]->GetBdata(),
                                                 m_base[2]->GetBdata(),
                                                 tmp0,  tmp3, wsp,
                                                 false, true, true);

                    // Scale eta_2 derivative by gfac1*gfac2
                    for(i = 0; i < nquad2; ++i)
                    {
                        Vmath::Smul(nq01, gfac2[i],
                                    &inarray[0] + i*nq01, 1,
                                    &tmp0   [0] + i*nq01, 1);
                    }
                    for(i = 0; i < nquad1*nquad2; ++i)
                    {
                        Vmath::Smul(nquad0, gfac1[i%nquad1],
                                            &tmp0[0] + i*nquad0, 1,
                                            &tmp0[0] + i*nquad0, 1);
                    }

                    MultiplyByStdQuadratureMetric(tmp0, tmp0);
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                 m_base[1]->GetDbdata(),
                                                 m_base[2]->GetBdata(),
                                                 tmp0, tmp4,  wsp,
                                                 true, false, true);

                    MultiplyByStdQuadratureMetric(inarray,tmp0);
                    IProductWRTBase_SumFacKernel(m_base[0]->GetBdata(),
                                                 m_base[1]->GetBdata(),
                                                 m_base[2]->GetDbdata(),
                                                 tmp0,outarray,wsp,
                                                 true, true, false);

                    Vmath::Vadd(m_ncoeffs,&tmp3[0],1,&outarray[0],1,&outarray[0],1);
                    Vmath::Vadd(m_ncoeffs,&tmp4[0],1,&outarray[0],1,&outarray[0],1);
                    break;
                }
                default:
                {
                    ASSERTL1(false, "input dir is out of range");
                    break;
                }
            }
        }

        //---------------------------------------
        // Evaluation functions
        //---------------------------------------
        
        void StdPyrExp::v_LocCoordToLocCollapsed(
            const Array<OneD, const NekDouble>& xi,
                  Array<OneD,       NekDouble>& eta)
        {
            if (fabs(xi[2]-1.0) < NekConstants::kNekZeroTol)
            {
                // Very top point of the pyramid
                eta[0] = -1.0;
                eta[1] = -1.0;
                eta[2] = xi[2];
            }
            else  
            {
                // Below the line-singularity -- Common case
                eta[2] = xi[2]; // eta_z = xi_z
                eta[1] = 2.0*(1.0 + xi[1])/(1.0 - xi[2]) - 1.0; 
                eta[0] = 2.0*(1.0 + xi[0])/(1.0 - xi[2]) - 1.0;
            } 
        }

        void StdPyrExp::v_GetCoords(Array<OneD, NekDouble> &xi_x, 
                                    Array<OneD, NekDouble> &xi_y,
                                    Array<OneD, NekDouble> &xi_z)
        {
            Array<OneD, const NekDouble> etaBar_x = m_base[0]->GetZ();
            Array<OneD, const NekDouble> eta_y    = m_base[1]->GetZ();
            Array<OneD, const NekDouble> eta_z    = m_base[2]->GetZ();
            int Qx = GetNumPoints(0);
            int Qy = GetNumPoints(1);
            int Qz = GetNumPoints(2);

            // Convert collapsed coordinates into cartesian coordinates: eta --> xi
            for (int k = 0; k < Qz; ++k )
            {
                for (int j = 0; j < Qy; ++j) 
                {
                    for (int i = 0; i < Qx; ++i) 
                    {
                        int s = i + Qx*(j + Qy*k);

                        xi_z[s] = eta_z[k];
                        xi_y[s] = (1.0 + eta_y[j]) * (1.0 - eta_z[k]) / 2.0  -  1.0;
                        xi_x[s] = (1.0 + etaBar_x[i]) * (1.0 - eta_z[k]) / 2.0  -  1.0;
                    }
                }
            }
        }

        void StdPyrExp::v_FillMode(const int mode, Array<OneD, NekDouble> &outarray)
        {
            Array<OneD, NekDouble> tmp(m_ncoeffs, 0.0);
            tmp[mode] = 1.0;
            v_BwdTrans(tmp, outarray);
        }
        
        
        //---------------------------------------
        // Helper functions
        //---------------------------------------

        int StdPyrExp::v_GetNverts() const
        {
            return 5;
        }

        int StdPyrExp::v_GetNedges() const
        {
            return 8;
        }

        int StdPyrExp::v_GetNfaces() const
        {
            return 5;
        }
        
        LibUtilities::ShapeType StdPyrExp::v_DetShapeType() const
        {
            return LibUtilities::ePyramid;
        }

        int StdPyrExp::v_NumBndryCoeffs() const
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_C ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            
            int P = m_base[0]->GetNumModes();
            int Q = m_base[1]->GetNumModes();
            int R = m_base[2]->GetNumModes();
            
            return LibUtilities::StdPyrData::
                                    getNumberOfBndCoefficients(P, Q, R);
        }

        int StdPyrExp::v_GetEdgeNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 7, "edge id is out of range");
            
            if (i == 0 || i == 2)
            {
                return GetBasisNumModes(0);
            }
            else if (i == 1 || i == 3)
            {
                return GetBasisNumModes(1);
            }
            else
            {
                return GetBasisNumModes(2);
            }
        }

        int StdPyrExp::v_GetFaceNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
            
            if (i == 0)
            {
                return GetBasisNumModes(0)*GetBasisNumModes(1);
            }
            else if (i == 1 || i == 3)
            {
                int P = GetBasisNumModes(0)-1, Q = GetBasisNumModes(2)-1;
                return Q+1 + (P*(1 + 2*Q - P))/2;
            }
            else
            {
                int P = GetBasisNumModes(1)-1, Q = GetBasisNumModes(2)-1;
                return Q+1 + (P*(1 + 2*Q - P))/2;
            }
        }
        
        int StdPyrExp::v_GetFaceIntNcoeffs(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");

            int P = m_base[0]->GetNumModes()-1;
            int Q = m_base[1]->GetNumModes()-1;
            int R = m_base[2]->GetNumModes()-1;

            if (i == 0)
            {
                return (P-1)*(Q-1);
            }
            else if (i == 1 || i == 3)
            {
                return (P-1) * (2*(R-1) - (P-1) - 1) / 2;
            }
            else
            {
                return (Q-1) * (2*(R-1) - (Q-1) - 1) / 2;
            }
        }

        int StdPyrExp::v_GetFaceNumPoints(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
            
            if (i == 0)
            {
                return m_base[0]->GetNumPoints()*
                       m_base[1]->GetNumPoints();
            }
            else if (i == 1 || i == 3)
            {
                return m_base[0]->GetNumPoints()*
                       m_base[2]->GetNumPoints();
            }
            else
            {
                return m_base[1]->GetNumPoints()*
                       m_base[2]->GetNumPoints();
            }
        }


        const LibUtilities::BasisKey StdPyrExp::v_DetFaceBasisKey(
            const int i, const int k) const
        {
            ASSERTL2(i >= 0 && i <= 4, "face id is out of range");
            ASSERTL2(k >= 0 && k <= 1, "basis key id is out of range");

            switch(i)
            {
                case 0:
                {
                    return EvaluateQuadFaceBasisKey(k,
                                                    m_base[k]->GetBasisType(),
                                                    m_base[k]->GetNumPoints(),
                                                    m_base[k]->GetNumModes());
                    
                }
                case 1:
                case 3:
                {
                    return EvaluateTriFaceBasisKey(k,
                                                   m_base[2*k]->GetBasisType(),
                                                   m_base[2*k]->GetNumPoints(),
                                                   m_base[2*k]->GetNumModes());
                }
                case 2:
                case 4:
                {
                    return EvaluateTriFaceBasisKey(k,
                                                   m_base[k+1]->GetBasisType(),
                                                   m_base[k+1]->GetNumPoints(),
                                                   m_base[k+1]->GetNumModes());
                }
            }

            // Should never get here.
            return LibUtilities::NullBasisKey;
        }

        int StdPyrExp::v_CalcNumberOfCoefficients(
            const std::vector<unsigned int> &nummodes, 
            int &modes_offset)
        {
            int nmodes = LibUtilities::StdPyrData::getNumberOfCoefficients(
                nummodes[modes_offset],
                nummodes[modes_offset+1],
                nummodes[modes_offset+2]);
            
            modes_offset += 3;
            return nmodes;
        }
        
        LibUtilities::BasisType StdPyrExp::v_GetEdgeBasisType(const int i) const
        {
            ASSERTL2(i >= 0 && i <= 7, "edge id is out of range");
            if (i == 0 || i == 2)
            {
                return GetBasisType(0);
            }
            else if (i == 1 || i == 3)
            {
                return GetBasisType(1);
            }
            else
            {
                return GetBasisType(2);
            }
        }


        //---------------------------------------
        // Mappings
        //---------------------------------------
        
        void StdPyrExp::v_GetFaceToElementMap(
            const int                  fid, 
            const Orientation          faceOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD,          int> &signarray,
            int                        P, 
            int                        Q)
        {
            ASSERTL1(GetEdgeBasisType(0) == GetEdgeBasisType(1),
                     "Method only implemented if BasisType is identical"
                     "in x and y directions");
            ASSERTL1(GetEdgeBasisType(0) == LibUtilities::eModified_A && 
                     GetEdgeBasisType(4) == LibUtilities::eModified_C,
                     "Method only implemented for Modified_A BasisType"
                     "(x and y direction) and Modified_C BasisType (z "
                     "direction)");

            int i, j, k, p, q, r, nFaceCoeffs;
            int nummodesA, nummodesB;
            
            int order0 = m_base[0]->GetNumModes();
            int order1 = m_base[1]->GetNumModes();
            int order2 = m_base[2]->GetNumModes();

            switch (fid)
            {
            case 0:
                nummodesA = order0;
                nummodesB = order1;
                break;
            case 1:
            case 3:
                nummodesA = order0;
                nummodesB = order2;
                break;
            case 2:
            case 4:
                nummodesA = order1;
                nummodesB = order2;
                break;
            }
            
            bool CheckForZeroedModes = false;

            if (P == -1)
            {
                P = nummodesA;
                Q = nummodesB;
                nFaceCoeffs = GetFaceNcoeffs(fid);
            }
            else if (fid > 0)
            {
                nFaceCoeffs = P*(2*Q-P+1)/2;
                CheckForZeroedModes = true;
            }
            else
            {
                nFaceCoeffs = P*Q;
                CheckForZeroedModes = true;
            }

            // Allocate the map array and sign array; set sign array to ones (+)
            if (maparray.num_elements() != nFaceCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nFaceCoeffs);
            }
            
            if (signarray.num_elements() != nFaceCoeffs)
            {
                signarray = Array<OneD, int>(nFaceCoeffs,1);
            }
            else
            {
                fill(signarray.get(), signarray.get() + nFaceCoeffs, 1);
            }

            // Set up an array indexing for quads, since the ordering may need
            // to be transposed.
            Array<OneD, int> arrayindx(nFaceCoeffs,-1);

            if (fid == 0)
            {
                for (i = 0; i < Q; i++)
                {
                    for (j = 0; j < P; j++)
                    {
                        if (faceOrient < 9)
                        {
                            arrayindx[i*P+j] = i*P+j;
                        }
                        else
                        {
                            arrayindx[i*P+j] = j*Q+i;
                        }
                    }
                }
            }

            // Set up ordering inside each 2D face. Also for triangular faces,
            // populate signarray.
            int cnt = 0, cnt2;
            switch (fid) 
            {
                case 0: // Bottom quad

                    // Fill in vertices
                    maparray[arrayindx[0]]           = 0;
                    maparray[arrayindx[1]]           = 1;
                    maparray[arrayindx[P+1]] = 2;
                    maparray[arrayindx[P]]   = 3;

                    // Edge 0
                    cnt = 5;
                    for (p = 2; p < P; ++p)
                    {
                        maparray[arrayindx[p]] = p-2 + cnt;
                    }

                    // Edge 1
                    cnt += P-2;
                    for (q = 2; q < Q; ++q)
                    {
                        maparray[arrayindx[q*P+1]] = q-2 + cnt;
                    }

                    // Edge 2
                    cnt += Q-2;
                    for (p = 2; p < P; ++p)
                    {
                        maparray[arrayindx[P+p]] = p-2 + cnt;
                    }

                    // Edge 3
                    cnt += P-2;
                    for (q = 2; q < Q; ++q)
                    {
                        maparray[arrayindx[q*P]] = q-2 + cnt;
                    }

                    // Interior
                    cnt += Q-2 + 4*(P-2);
                    for (q = 2; q < Q; ++q)
                    {
                        for (p = 2; p < P; ++p)
                        {
                            maparray[arrayindx[q*P+p]] = cnt + (q-2)*P+(p-2);
                        }
                    }


                    break;

                case 1: // Left triangle
                    // Vertices
                    maparray[0]         = 0;
                    maparray[1]         = 4;
                    maparray[Q] = 1;

                    // Edge 0 (pyramid edge 0)
                    cnt = 5;
                    q   = 2*Q-1;
                    for (p = 2; p < P; q += Q-p, ++p)
                    {
                        maparray[q] = cnt++;
                        if ((int)faceOrient == 7)
                        {
                            signarray[q] = p % 2 ? -1 : 1;
                        }
                    }

                    // Edge 1 (pyramid edge 5)
                    cnt = 5 + 2*(order0-2) + 2*(order1-2) + (order2-2);
                    for (q = 2; q < Q; ++q)
                    {
                        maparray[q] = cnt++;
                        if ((int)faceOrient == 7)
                        {
                            signarray[q] = q % 2 ? -1 : 1;
                        }
                    }
                    
                    // Edge 2 (pyramid edge 4)
                    cnt = 5 + 2*(order0-2) + 2*(order1-2);
                    for (q = 2; q < Q; ++q)
                    {
                        maparray[Q+q-1] = cnt++;
                        if ((int)faceOrient == 7)
                        {
                            signarray[Q+q-1] = q % 2 ? -1 : 1;
                        }
                    }

                    // Interior
                    cnt  = 5 + 2*(order0-2) + 2*(order1-2) + 4*(order2-2)
                        + v_GetFaceIntNcoeffs(0);
                    cnt2 = 2*Q + 1;
                    for (p = 2; p < P; ++p)
                    {
                        for (r = 2; r < Q-p; ++r)
                        {
                            maparray[cnt2] = cnt++;
                            if ((int)faceOrient == 7 && p > 1)
                            {
                                signarray[cnt2++] = p % 2 ? -1 : 1;
                            }
                        }
                        cnt2++;
                    }
                    break;

                case 2:
                    // Vertices
                    maparray[0]         = 1;
                    maparray[1]         = 4;
                    maparray[Q] = 2;

                    // Edge 0 (pyramid edge 1)
                    cnt = 5 + (order0-2);
                    q   = 2*Q-1;
                    for (p = 2; p < P; q += Q-p, ++p)
                    {
                        maparray[q] = cnt++;
                        if ((int)faceOrient == 7)
                        {
                            signarray[q] = p % 2 ? -1 : 1;
                        }
                    }

                    // Edge 1 (pyramid edge 6)
                    cnt = 5 + 2*(order0-2) + 2*(order1-2) + 2*(order2-2);
                    for (q = 2; q < Q; ++q)
                    {
                        maparray[q] = cnt++;
                        if ((int)faceOrient == 7)
                        {
                            signarray[q] = q % 2 ? -1 : 1;
                        }
                    }
                    
                    // Edge 2 (pyramid edge 5)
                    cnt = 5 + 2*(order0-2) + 2*(order1-2) + (order2-2);
                    for (q = 2; q < Q; ++q)
                    {
                        maparray[Q+q-1] = cnt++;
                        if ((int)faceOrient == 7)
                        {
                            signarray[Q+q-1] = q % 2 ? -1 : 1;
                        }
                    }

                    // Interior
                    cnt  = 5 + 2*(order0-2) + 2*(order1-2) + 4*(order2-2)
                        + v_GetFaceIntNcoeffs(0) + v_GetFaceIntNcoeffs(1);
                    cnt2 = 2*Q + 1;
                    for (p = 2; p < P; ++p)
                    {
                        for (r = 2; r < Q-p; ++r)
                        {
                            maparray[cnt2] = cnt++;
                            if ((int)faceOrient == 7 && p > 1)
                            {
                                signarray[cnt2++] = p % 2 ? -1 : 1;
                            }
                        }
                        cnt2++;
                    }
                    break;

                case 3: // Right triangle
                    // Vertices
                    maparray[0]         = 3;
                    maparray[1]         = 4;
                    maparray[Q] = 2;

                    // Edge 0 (pyramid edge 2)
                    cnt = 5 + (order0-2) + (order1-2);
                    q   = 2*Q-1;
                    for (p = 2; p < P; q += Q-p, ++p)
                    {
                        maparray[q] = cnt++;
                        if ((int)faceOrient == 7)
                        {
                            signarray[q] = p % 2 ? -1 : 1;
                        }
                    }

                    // Edge 1 (pyramid edge 6)
                    cnt = 5 + 2*(order0-2) + 2*(order1-2) + 2*(order2-2);
                    for (q = 2; q < Q; ++q)
                    {
                        maparray[q] = cnt++;
                        if ((int)faceOrient == 7)
                        {
                            signarray[q] = q % 2 ? -1 : 1;
                        }
                    }
                    
                    // Edge 2 (pyramid edge 7)
                    cnt = 5 + 2*(order0-2) + 2*(order1-2) + 3*(order2-2);
                    for (q = 2; q < Q; ++q)
                    {
                        maparray[Q+q-1] = cnt++;
                        if ((int)faceOrient == 7)
                        {
                            signarray[Q+q-1] = q % 2 ? -1 : 1;
                        }
                    }

                    // Interior
                    cnt  = 5 + 2*(order0-2) + 2*(order1-2) + 4*(order2-2)
                        + v_GetFaceIntNcoeffs(0) + v_GetFaceIntNcoeffs(1)
                        + v_GetFaceIntNcoeffs(2);
                    cnt2 = 2*Q + 1;
                    for (p = 2; p < P; ++p)
                    {
                        for (r = 2; r < Q-p; ++r)
                        {
                            maparray[cnt2] = cnt++;
                            if ((int)faceOrient == 7 && p > 1)
                            {
                                signarray[cnt2++] = p % 2 ? -1 : 1;
                            }
                        }
                        cnt2++;
                    }
                    break;

                case 4: // Rear tri
                    // Vertices
                    maparray[0]         = 0;
                    maparray[1]         = 4;
                    maparray[Q] = 3;

                    // Edge 0 (pyramid edge 3)
                    cnt = 5 + 2*(order0-2) + (order1-2);
                    q   = 2*Q-1;
                    for (p = 2; p < P; q += Q-p, ++p)
                    {
                        maparray[q] = cnt++;
                        if ((int)faceOrient == 7)
                        {
                            signarray[q] = p % 2 ? -1 : 1;
                        }
                    }

                    // Edge 1 (pyramid edge 7)
                    cnt = 5 + 2*(order0-2) + 2*(order1-2) + 3*(order2-2);
                    for (q = 2; q < Q; ++q)
                    {
                        maparray[q] = cnt++;
                        if ((int)faceOrient == 7)
                        {
                            signarray[q] = q % 2 ? -1 : 1;
                        }
                    }
                    
                    // Edge 2 (pyramid edge 4)
                    cnt = 5 + 2*(order0-2) + 2*(order1-2);
                    for (q = 2; q < Q; ++q)
                    {
                        maparray[Q+q-1] = cnt++;
                        if ((int)faceOrient == 7)
                        {
                            signarray[Q+q-1] = q % 2 ? -1 : 1;
                        }
                    }

                    // Interior
                    cnt  = 5 + 2*(order0-2) + 2*(order1-2) + 4*(order2-2)
                        + v_GetFaceIntNcoeffs(0) + v_GetFaceIntNcoeffs(1)
                        + v_GetFaceIntNcoeffs(2) + v_GetFaceIntNcoeffs(3);
                    cnt2 = 2*Q + 1;
                    for (p = 2; p < P; ++p)
                    {
                        for (r = 2; r < Q-p; ++r)
                        {
                            maparray[cnt2] = cnt++;
                            if ((int)faceOrient == 7 && p > 1)
                            {
                                signarray[cnt2++] = p % 2 ? -1 : 1;
                            }
                        }
                        cnt2++;
                    }
                    break;
                    
                default:
                    ASSERTL0(false, "Face to element map unavailable.");
            }

            if (fid > 0)
            {

               if(CheckForZeroedModes)
                {
                    // zero signmap and set maparray to zero if elemental
                    // modes are not as large as face modesl
                    int idx = 0; 
                    for (j = 0; j < P; ++j)
                    {
                        idx += Q-j;
                        for (k = Q-j; k < Q-j; ++k)
                        {
                            signarray[idx]  = 0.0;
                            maparray[idx++] = maparray[0];
                        }
                    }
                    
                    for (j = P; j < P; ++j)
                    {
                        for (k = 0; k < Q-j; ++k)
                        {
                            signarray[idx]  = 0.0;
                            maparray[idx++] = maparray[0];
                        }
                    }
                }
                
                // Triangles only have one possible orientation (base
                // direction reversed); swap edge modes.
                if ((int)faceOrient == 7)
                {
                    swap(maparray[0], maparray[P]);
                    for (i = 1; i < P-1; ++i)
                    {
                        swap(maparray[i+1], maparray[P+i]);
                    }
                }
            }
            else
            {
                if(CheckForZeroedModes)
                {
                    // zero signmap and set maparray to zero if elemental
                    // modes are not as large as face modesl
                    for (j = 0; j < P; ++j)
                    {
                        for (k = Q; k < Q; ++k)
                        {
                            signarray[arrayindx[j+k*P]] = 0.0;
                            maparray[arrayindx[j+k*P]]  = maparray[0];
                        }
                    }

                    for (j = P; j < P; ++j)
                    {
                        for (k = 0; k < Q; ++k)
                        {
                            signarray[arrayindx[j+k*P]] = 0.0;
                            maparray[arrayindx[j+k*P]]  = maparray[0];
                        }
                    }                    
                }

                // The code below is exactly the same as that taken from
                // StdHexExp and reverses the 'b' and 'a' directions as
                // appropriate (1st and 2nd if statements respectively) in
                // quadrilateral faces.
                if (faceOrient == 6 || faceOrient == 8 ||
                    faceOrient == 11 || faceOrient == 12)
                {
                    if (faceOrient < 9)
                    {
                        for (i = 3; i < Q; i += 2)
                        {
                            for (j = 0; j < P; j++)
                            {
                                signarray[arrayindx[i*P+j]] *= -1;
                            }
                        }
                        
                        for (i = 0; i < P; i++)
                        {
                            swap(maparray [i], maparray [i+P]);
                            swap(signarray[i], signarray[i+P]);
                        }
                    }
                    else
                    {
                        for (i = 0; i < Q; i++)
                        {
                            for (j = 3; j < P; j += 2)
                            {
                                signarray[arrayindx[i*P+j]] *= -1;
                            }
                        }

                        for (i = 0; i < Q; i++)
                        {
                            swap (maparray [i], maparray [i+Q]);
                            swap (signarray[i], signarray[i+Q]);
                        }
                    }
                }

                if (faceOrient == 7 || faceOrient == 8 ||
                    faceOrient == 10 || faceOrient == 12)
                {
                    if (faceOrient < 9)
                    {
                        for (i = 0; i < Q; i++)
                        {
                            for (j = 3; j < P; j += 2)
                            {
                                signarray[arrayindx[i*P+j]] *= -1;
                            }
                        }

                        for(i = 0; i < Q; i++)
                        {
                            swap(maparray [i*P], maparray [i*P+1]);
                            swap(signarray[i*P], signarray[i*P+1]);
                        }
                    }
                    else
                    {
                        for (i = 3; i < Q; i += 2)
                        {
                            for (j = 0; j < P; j++)
                            {
                                signarray[arrayindx[i*P+j]] *= -1;
                            }
                        }

                        for (i = 0; i < P; i++)
                        {
                            swap(maparray [i*Q], maparray [i*Q+1]);
                            swap(signarray[i*Q], signarray[i*Q+1]);
                        }
                    }
                }
            }
        }

        int StdPyrExp::v_GetVertexMap(int vId, bool useCoeffPacking)
        {
            ASSERTL0(GetEdgeBasisType(vId) == LibUtilities::eModified_A ||
                     GetEdgeBasisType(vId) == LibUtilities::eModified_A ||
                     GetEdgeBasisType(vId) == LibUtilities::eModified_C,
                     "Mapping not defined for this type of basis");
            return vId;
        }

        void StdPyrExp::v_GetEdgeInteriorMap(
            const int                  eid,
            const Orientation          edgeOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD, int>          &signarray)
        {
            int       i;
            bool      signChange;
            const int P              = m_base[0]->GetNumModes() - 2;
            const int Q              = m_base[1]->GetNumModes() - 2;
            const int R              = m_base[2]->GetNumModes() - 2;
            const int nEdgeIntCoeffs = v_GetEdgeNcoeffs(eid) - 2;
            
            if (maparray.num_elements() != nEdgeIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nEdgeIntCoeffs);
            }
            
            if(signarray.num_elements() != nEdgeIntCoeffs)
            {
                signarray = Array<OneD, int>(nEdgeIntCoeffs,1);
            }
            else
            {
                fill(signarray.get(), signarray.get()+nEdgeIntCoeffs, 1);
            }
            
            // If edge is oriented backwards, change sign of modes which have
            // degree 2n+1, n >= 1.
            signChange = edgeOrient == eBackwards;

            int offset = 5;
            
            switch (eid)
            {
                case 0:
                    break;
                case 1:
                    offset += P;
                    break;
                case 2:
                    offset += P+Q;
                    break;
                case 3:
                    offset += 2*P+Q;
                    break;
                case 4:
                    offset += 2*(P+Q);
                    break;
                case 5:
                    offset += 2*(P+Q)+R;
                    break;
                case 6:
                    offset += 2*(P+Q+R);
                    break;
                case 7:
                    offset += 2*(P+Q)+3*R;
                    break;
                default:
                    ASSERTL0(false, "Edge not defined.");
                    break;
            }

            for (i = 0; i < nEdgeIntCoeffs; ++i)
            {
                maparray[i] = offset + i;
            }

            if (signChange)
            {
                for (i = 1; i < nEdgeIntCoeffs; i += 2)
                {
                    signarray[i] = -1;
                }
            }
        }
        
        void StdPyrExp::v_GetFaceInteriorMap(
            const int                  fid,
            const Orientation          faceOrient,
            Array<OneD, unsigned int> &maparray,
            Array<OneD, int>          &signarray)
        {
            const int P              = m_base[0]->GetNumModes() - 1;
            const int Q              = m_base[1]->GetNumModes() - 1;
            const int R              = m_base[2]->GetNumModes() - 1;
            const int nFaceIntCoeffs = v_GetFaceIntNcoeffs(fid);
            int       p, q, r, idx   = 0;
            int       nummodesA      = 0;
            int       nummodesB      = 0;
            int       i, j;

            if (maparray.num_elements() != nFaceIntCoeffs)
            {
                maparray = Array<OneD, unsigned int>(nFaceIntCoeffs);
            }
            
            if (signarray.num_elements() != nFaceIntCoeffs)
            {
                signarray = Array<OneD, int>(nFaceIntCoeffs, 1);
            }
            else
            {
                fill(signarray.get(), signarray.get()+nFaceIntCoeffs, 1);
            }

            // Set up an array indexing for quad faces, since the ordering may
            // need to be transposed depending on orientation.
            Array<OneD, int> arrayindx(nFaceIntCoeffs);
            if (fid == 0) 
            {
                nummodesA = P-1;
                nummodesB = Q-1;

                for (i = 0; i < nummodesB; i++)
                {
                    for (j = 0; j < nummodesA; j++)
                    {
                        if (faceOrient < 9)
                        {
                            arrayindx[i*nummodesA+j] = i*nummodesA+j;
                        }
                        else
                        {
                            arrayindx[i*nummodesA+j] = j*nummodesB+i;
                        }
                    }
                }
            }

            int offset = 5 + 2*(P-1) + 2*(Q-1) + 4*(R-1);

            for (i = 0; i < fid; ++i)
            {
                offset += v_GetFaceIntNcoeffs(i);
            }

            switch (fid)
            {
                case 0:
                    for (q = 2; q <= Q; ++q)
                    {
                        for (p = 2; p <= P; ++p)
                        {
                            maparray[arrayindx[(q-2)*nummodesA+(p-2)]]
                                    = offset + (q-2)*nummodesA+(p-2);
                        }
                    }
                    break;

                case 1:
                case 3:
                    for (p = 2; p <= P; ++p)
                    {
                        for (r = 1; r <= R-p; ++r, ++idx)
                        {
                            if ((int)faceOrient == 7)
                            {
                                signarray[idx] = p % 2 ? -1 : 1;
                            }
                            maparray[idx] = offset + idx;
                        }
                    }
                    break;

                case 2:
                case 4:
                    for (q = 2; q <= Q; ++q)
                    {
                        for (r = 1; r <= R-q; ++r, ++idx)
                        {
                            if ((int)faceOrient == 7)
                            {
                                signarray[idx] = q % 2 ? -1 : 1;
                            }
                            maparray[idx] = offset + idx;
                        }
                    }
                    break;

                default:
                    ASSERTL0(false, "Face interior map not available.");
            }

            // Triangular faces are processed in the above switch loop; for
            // remaining quad faces, set up orientation if necessary.
            if (fid > 0)
            {
                return;
            }

            if (faceOrient == 6 || faceOrient == 8 ||
                faceOrient == 11 || faceOrient == 12)
            {
                if (faceOrient < 9)
                {
                    for (i = 1; i < nummodesB; i += 2)
                    {
                        for (j = 0; j < nummodesA; j++)
                        {
                            signarray[arrayindx[i*nummodesA+j]] *= -1;
                        }
                    }
                }
                else
                {
                    for (i = 0; i < nummodesB; i++)
                    {
                        for (j = 1; j < nummodesA; j += 2)
                        {
                            signarray[arrayindx[i*nummodesA+j]] *= -1;
                        }
                    }
                }
            }

            if (faceOrient == 7 || faceOrient == 8 ||
                faceOrient == 10 || faceOrient == 12)
            {
                if (faceOrient < 9)
                {
                    for (i = 0; i < nummodesB; i++)
                    {
                        for (j = 1; j < nummodesA; j += 2)
                        {
                            signarray[arrayindx[i*nummodesA+j]] *= -1;
                        }
                    }
                }
                else
                {
                    for (i = 1; i < nummodesB; i += 2)
                    {
                        for (j = 0; j < nummodesA; j++)
                        {
                            signarray[arrayindx[i*nummodesA+j]] *= -1;
                        }
                    }
                }
            }
        }

        void StdPyrExp::v_GetInteriorMap(Array<OneD, unsigned int> &outarray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_C ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            const int nBndCoeffs = v_NumBndryCoeffs();
            const int nIntCoeffs = m_ncoeffs - NumBndryCoeffs();

            if (outarray.num_elements() != nIntCoeffs)
            {
                outarray = Array<OneD, unsigned int>(nIntCoeffs);
            }

            // Loop over all interior modes.
            int p, idx = 0;
            for (p = nBndCoeffs; p < m_ncoeffs; ++p)
            {
                outarray[idx++] = p;
            }
        }

        void StdPyrExp::v_GetBoundaryMap(Array<OneD, unsigned int> &maparray)
        {
            ASSERTL1(GetBasisType(0) == LibUtilities::eModified_A ||
                     GetBasisType(0) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(1) == LibUtilities::eModified_A ||
                     GetBasisType(1) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");
            ASSERTL1(GetBasisType(2) == LibUtilities::eModified_C ||
                     GetBasisType(2) == LibUtilities::eGLL_Lagrange,
                     "BasisType is not a boundary interior form");

            int idx = 0, nBndry = v_NumBndryCoeffs();

            for (idx = 0; idx < nBndry; ++idx)
            {
                maparray[idx] = idx;
            }
        }

        //---------------------------------------
        // Wrapper functions
        //---------------------------------------
        
        DNekMatSharedPtr StdPyrExp::v_GenMatrix(const StdMatrixKey &mkey)
        {
            return CreateGeneralMatrix(mkey);
        }
        
        DNekMatSharedPtr StdPyrExp::v_CreateStdMatrix(const StdMatrixKey &mkey)
        {
            return v_GenMatrix(mkey);
        }

        /**
         * @brief Number tetrahedral modes in r-direction. Much the same as
         * StdTetExp::GetTetMode but slightly simplified since we know that the
         * polynomial order is the same in each direction.
         */
        int StdPyrExp::GetTetMode(const int I, const int J, const int K)
        {
            const int R = m_base[2]->GetNumModes();
            int i, j, cnt = 0;
            for (i = 0; i < I; ++i)
            {
                cnt += (R-i)*(R-i+1)/2;
            }

            i = R-I;
            for (j = 0; j < J; ++j)
            {
                cnt += i;
                i--;
            }

            return cnt + K;
        }

        void StdPyrExp::v_MultiplyByStdQuadratureMetric(
            const Array<OneD, const NekDouble>& inarray,
                  Array<OneD,       NekDouble>& outarray)
        {
            int i, j;
            
            int  nquad0 = m_base[0]->GetNumPoints();
            int  nquad1 = m_base[1]->GetNumPoints();
            int  nquad2 = m_base[2]->GetNumPoints();
            
            const Array<OneD, const NekDouble>& w0 = m_base[0]->GetW();
            const Array<OneD, const NekDouble>& w1 = m_base[1]->GetW();
            const Array<OneD, const NekDouble>& w2 = m_base[2]->GetW();
            
            const Array<OneD, const NekDouble>& z2 = m_base[2]->GetZ();
            
            // Multiply by integration constants in x-direction
            for(i = 0; i < nquad1*nquad2; ++i)
            {
                Vmath::Vmul(nquad0, inarray.get()+i*nquad0, 1,
                            w0.get(), 1, outarray.get()+i*nquad0,1);
            }
            
            // Multiply by integration constants in y-direction
            for(j = 0; j < nquad2; ++j)
            {
                for(i = 0; i < nquad1; ++i)
                {
                    Blas::Dscal(nquad0,w1[i], &outarray[0]+i*nquad0 +
                                j*nquad0*nquad1,1);
                }
            }
            
            // Multiply by integration constants in z-direction; need to
            // incorporate factor [(1-eta_3)/2]^2 into weights, but only if
            // using GLL quadrature points.
            switch(m_base[2]->GetPointsType())
            {
                // (2,0) Jacobi inner product.
                case LibUtilities::eGaussRadauMAlpha2Beta0:
                    for(i = 0; i < nquad2; ++i)
                    {
                        Blas::Dscal(nquad0*nquad1, 0.25*w2[i],
                                    &outarray[0]+i*nquad0*nquad1, 1);
                    }
                    break;
                
                default:
                    for(i = 0; i < nquad2; ++i)
                    {
                        Blas::Dscal(nquad0*nquad1,0.125*(1-z2[i])*(1-z2[i])*w2[i],
                                    &outarray[0]+i*nquad0*nquad1,1);
                    }
                    break;
            }
        }
    }
}
