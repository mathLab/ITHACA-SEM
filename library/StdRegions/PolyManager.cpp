///////////////////////////////////////////////////////////////////////////////
//
// File PolyManager.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: Polynomial routines definition and manager 
//
///////////////////////////////////////////////////////////////////////////////

#include <StdRegions/PolyManager.h>
#include <StdRegions/StdBasis.h>


namespace Nektar
{
    namespace StdRegions 
    {

        /*====================================================================
        Routines to do with Zero and Weights management
        ====================================================================*/

        void PolyManager::GetZW(const PointsType ptype, const int  np, 
            const double * &z,     const double * &w, 
            const double alpha,    const double beta)
        {
            ZWLookup(ptype,np,z,w,alpha,beta,m_zwvals[(int)ptype]);
        }


        void PolyManager::GetZW(const BasisKey *B, const double * &z,  
            const double * &w)
        {
            ZWLookup(B->GetPointsType(),B->GetPointsOrder(),z,w,
                B->GetAlpha(),B->GetBeta(),m_zwvals[(int)B->GetPointsType()]);
        }

        void PolyManager:: ResetZW(void)
        {
            std::vector<ZeroWeights*>::iterator defzw;

            for(int i = 0; i < SIZE_PointsType; ++i)
            {
                for(defzw = m_zwvals[i].begin(); defzw != m_zwvals[i].end(); ++defzw)
                {
                    delete defzw[0];
                }
                m_zwvals[i].erase(m_zwvals[i].begin(),m_zwvals[i].end());
                m_zwvalscur[i] = m_zwvals[i].begin();
            }
        }

        void PolyManager::ZWLookup(const PointsType ptype, const int np, 
            const double* &z, const double* &w, 
            const double alpha, const double beta,
            std::vector<ZeroWeights*> &val)
        {
            std::vector<ZeroWeights*>::iterator def;
            ZeroWeightsKey key(ptype,np,alpha,beta);

            def = find(val.begin(),val.end(),m_zwvalscur[(int)ptype],key);

            if(def != val.end()){ // return found values
                z = def[0]->GetZ();
                w = def[0]->GetW();
                m_zwvalscur[(int)ptype] = def;
            }
            else{ // add new definition to list 
                ZeroWeights * tmp = new ZeroWeights(key);
                val.push_back(tmp);
                z = tmp->GetZ();
                w = tmp->GetW();
                m_zwvalscur[(int)ptype] = (--val.end());
            }
        }

        void PolyManager::GetInterpVec(const double zi, const PointsType ptype, 
            const double *z, const int np, 
            const double alpha, const double beta, 
            double *I)
        {
            switch(ptype)
            {
                case eGauss:
                    Polylib::Imgj(I,z,&zi,np,1,alpha,beta);
                    break;

                case eLobatto:
                    Polylib::Imglj(I,z,&zi,np,1,alpha,beta);
                    break;

                case eRadauM:
                case eRadauP:
                    Polylib::Imgrjm(I,z,&zi,np,1,alpha,beta);
                    break;

                case ePolyEvenSp:
                    ASSERTL0(false, "Need to define routine for equispaced polynomial weights");
                    break;

                case eFourierEvenSp:
                    ASSERTL0(!np%2, "Fourier points need to be of even order");
                    break;

                default:
                    ASSERTL0(false, "Unknown point type distribution");
            }
        }

        void ZeroWeights::SetZW(const PointsType ptype, const int np, 
            const double alpha, const double beta)
        {

            m_pointstype = ptype;
            m_pointsorder    = np;
            m_alpha = alpha;
            m_beta  = beta;

            if(m_pointsorder){
                m_z = new double[m_pointsorder];
                m_w = new double[m_pointsorder];

                switch(ptype){
                case eGauss:
                    Polylib::zwgj(m_z,m_w,m_pointsorder,m_alpha,m_beta);
                    break;
                case eLobatto:
                    Polylib::zwglj(m_z,m_w,m_pointsorder,m_alpha,m_beta);
                    break;
                case eRadauM:
                    Polylib::zwgrjm(m_z,m_w,m_pointsorder,m_alpha,m_beta);
                    break;
                case eRadauP:
                    Polylib::zwgrjp(m_z,m_w,m_pointsorder,m_alpha,m_beta);
                    break;
                case ePolyEvenSp:
                    // define points including end points 
                    for(int i = 0; i < m_pointsorder; ++i)
                    {
                        m_z[i] = -1.0 + i*2.0/(double)(m_pointsorder-1);
                    }

                    ASSERTL0(false, "Need to define routine for equispaced polynomial weights");
                    break;

                case eFourierEvenSp:
                    ASSERTL0(!m_pointsorder%2, "Fourier points need to be of even order");

                    // define points in the region [-1,1]
                    for(int i = 0; i < m_pointsorder; ++i)
                    {
                        m_z[i] = -1.0 + i*2.0/(double)m_pointsorder;
                        m_w[i] =  2.0/(double)m_pointsorder;
                    }
                    break;
                default:
                    ASSERTL0(false, "Unknown point type distribution");
                }
            }
            else
            {
                m_z = (double *) NULL; 
                m_w = (double *) NULL; 
            }
        }


        bool operator == (const ZeroWeights &x, const ZeroWeightsKey &y)
        {
            if( (x.m_pointsorder == y.m_pointsorder) &&
                (x.m_pointstype == y.m_pointstype) &&
                (x.m_alpha == y.m_alpha) &&
                (x.m_beta == y.m_beta) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }


        bool operator == (const ZeroWeightsKey &x, const ZeroWeights &y)
        {
            if( (x.m_pointsorder == y.m_pointsorder) &&
                (x.m_pointstype == y.m_pointstype) &&
                (x.m_alpha == y.m_alpha) &&
                (x.m_beta == y.m_beta) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        bool operator == (const ZeroWeightsKey &x, const ZeroWeights *y)
        {
            if( (x.m_pointsorder    == (*y).m_pointsorder) &&
                (x.m_pointstype == (*y).m_pointstype) &&
                (x.m_alpha == (*y).m_alpha) &&
                (x.m_beta  == (*y).m_beta) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        bool operator == (const ZeroWeights *x, const ZeroWeightsKey &y)
        {
            if( ((*x).m_pointsorder == y.m_pointsorder) &&
                ((*x).m_pointstype == y.m_pointstype) &&
                ((*x).m_alpha == y.m_alpha) &&
                ((*x).m_beta == y.m_beta))
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        bool operator != (const ZeroWeights &x, const ZeroWeightsKey &y)
        {
            if( (x.m_pointsorder == y.m_pointsorder) &&
                (x.m_pointstype == y.m_pointstype) &&
                (x.m_alpha == y.m_alpha) &&
                (x.m_beta == y.m_beta) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        bool operator != (const ZeroWeightsKey &x, const ZeroWeights &y)
        {
            if( (x.m_pointsorder == y.m_pointsorder) &&
                (x.m_pointstype == y.m_pointstype) &&
                (x.m_alpha == y.m_alpha) &&
                (x.m_beta == y.m_beta) )
            {   
                return false;
            }
            else
            {
                return true;
            }
        }

        bool operator != (const ZeroWeightsKey &x, const ZeroWeights *y)
        {
            if( (x.m_pointsorder    == (*y).m_pointsorder) &&
                (x.m_pointstype == (*y).m_pointstype) &&
                (x.m_alpha == (*y).m_alpha) &&
                (x.m_beta  == (*y).m_beta) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        bool operator != (const ZeroWeights *x, const ZeroWeightsKey &y)
        {
            if( ((*x).m_pointsorder == y.m_pointsorder) &&
                ((*x).m_pointstype == y.m_pointstype) &&
                ((*x).m_alpha == y.m_alpha) &&
                ((*x).m_beta == y.m_beta) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }


        /*====================================================================
        Routines to derivative matrix  management
        ===================================================================*/

        void PolyManager::GetD(const PointsType ptype, const int  np, 
            const double * &D, const double alpha, 
            const double beta)
        {
            DLookup(ptype,np,D,alpha,beta,m_dmvals[(int)ptype]);
        }


        void PolyManager::GetD(const BasisKey *B, const double * &D)
        {
            DLookup(B->GetPointsType(),B->GetPointsOrder(),D,B->GetAlpha(),
                B->GetBeta(),m_dmvals[(int)B->GetPointsType()]);
        }

        void PolyManager:: ResetD(void)
        {
            std::vector<DerivMatrix*>::iterator defder;

            for(int i = 0; i < SIZE_PointsType; ++i)
            {
                for(defder = m_dmvals[i].begin(); defder != m_dmvals[i].end(); ++defder)
                {
                    delete defder[0];
                }
                m_dmvals[i].erase(m_dmvals[i].begin(),m_dmvals[i].end());
                m_dmvalscur[i] = m_dmvals[i].begin();
            }
        }

        void PolyManager::DLookup(const PointsType ptype, const int np, 
            const double* &D, const double alpha, 
            const double beta, std::vector<DerivMatrix*> &val)
        {

            std::vector<DerivMatrix*>::iterator def;
            DerivMatrixKey key(ptype,np,alpha,beta);

            def = find(val.begin(),val.end(),m_dmvalscur[(int)ptype],key);

            if(def != val.end()){
                D = def[0]->GetD();
                m_dmvalscur[(int)ptype] = def;
            }
            else{
                DerivMatrix * tmp = new DerivMatrix(this,key);
                val.push_back(tmp);
                D = tmp->GetD();
                m_dmvalscur[(int)ptype] = (--val.end());
            }
        }

        void DerivMatrix::SetD(PolyManager *P, const PointsType ptype, const int np, 
            const double alpha, const double beta)
        {
            m_pointstype = ptype;
            m_pointsorder    = np;
            m_alpha = alpha;
            m_beta  = beta;

            const double *z,*w;

            if(m_pointsorder){
                m_dm = new double[m_pointsorder*m_pointsorder];

                // get zeros for calculation of derivative matrices
                P->GetZW(ptype,m_pointsorder,z,w,alpha,beta);

                switch(ptype){
                    case eGauss:
                        Polylib::Dgj(m_dm,z,m_pointsorder,m_alpha,m_beta);
                        break;
                    case eLobatto:
                        Polylib::Dglj(m_dm,z,m_pointsorder,m_alpha,m_beta);
                        break;
                    case eRadauM:
                        Polylib::Dgrjm(m_dm,z,m_pointsorder,m_alpha,m_beta);
                        break;
                    case eRadauP:
                        Polylib::Dgrjp(m_dm,z,m_pointsorder,m_alpha,m_beta);
                        break;
                    case ePolyEvenSp:
                        ASSERTL0(false, "Need to define routine for equispaced polynomial Derivative");
                        break;
                    case eFourierEvenSp:
                        ASSERTL0(false, "Fourier points Derivative matrix need defining for Fourier spacing");
                        break;
                    default:
                        ASSERTL0(false, "Unknown point type distribution");
                }
            }
        }


        bool operator == (const DerivMatrix &x, const DerivMatrixKey &y)
        {
            if( (x.m_pointsorder == y.m_pointsorder) &&
                (x.m_pointstype == y.m_pointstype) &&
                (x.m_alpha == y.m_alpha) &&
                (x.m_beta == y.m_beta) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }


        bool operator == (const DerivMatrixKey &x, const DerivMatrix &y)
        {
            if( (x.m_pointsorder == y.m_pointsorder) &&
                (x.m_pointstype == y.m_pointstype) &&
                (x.m_alpha == y.m_alpha) &&
                (x.m_beta == y.m_beta) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        bool operator == (const DerivMatrixKey &x, const DerivMatrix *y)
        {
            if( (x.m_pointsorder == (*y).m_pointsorder) &&
                (x.m_pointstype == (*y).m_pointstype) &&
                (x.m_alpha == (*y).m_alpha) &&
                (x.m_beta == (*y).m_beta) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        bool operator == (const DerivMatrix *x, const DerivMatrixKey &y)
        {
            if( ((*x).m_pointsorder == y.m_pointsorder) &&
                ((*x).m_pointstype == y.m_pointstype) &&
                ((*x).m_alpha == y.m_alpha) &&
                ((*x).m_beta == y.m_beta) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        bool operator != (const DerivMatrix &x, const DerivMatrixKey &y)
        {
            if( (x.m_pointsorder == y.m_pointsorder) &&
                (x.m_pointstype == y.m_pointstype) &&
                (x.m_alpha == y.m_alpha) &&
                (x.m_beta == y.m_beta) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        bool operator != (const DerivMatrixKey &x, const DerivMatrix &y)
        {
            if( (x.m_pointsorder == y.m_pointsorder) &&
                (x.m_pointstype == y.m_pointstype) &&
                (x.m_alpha == y.m_alpha) &&
                (x.m_beta == y.m_beta) )   
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        bool operator != (const DerivMatrixKey &x, const DerivMatrix *y)
        {
            if( (x.m_pointsorder == (*y).m_pointsorder) &&
                (x.m_pointstype == (*y).m_pointstype) &&
                (x.m_alpha == (*y).m_alpha) &&
                (x.m_beta == (*y).m_beta) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        bool operator != (const DerivMatrix *x, const DerivMatrixKey &y)
        {
            if( ((*x).m_pointsorder == y.m_pointsorder) &&
                ((*x).m_pointstype == y.m_pointstype) &&
                ((*x).m_alpha == y.m_alpha) &&
                ((*x).m_beta == y.m_beta) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }


        /*====================================================================
        Routines to interpolation matrix management
        ===================================================================*/

        void PolyManager::GetI(const PointsType orig_ptype, const int orig_np,  
            const double orig_alpha, const double orig_beta,
            const PointsType dest_ptype, const int dest_np,  
            const double dest_alpha, const double dest_beta, 
            const double * &I)
        {
            ILookup(orig_ptype, orig_np, orig_alpha, orig_beta, dest_ptype,
                dest_np,dest_alpha,dest_beta,I, m_imvals[(int)orig_ptype]);
        }

        void PolyManager::GetI(const BasisKey *orig_B, const BasisKey *dest_B, 
            const double * &I)
        {
            ILookup(orig_B->GetPointsType(), orig_B->GetPointsOrder(),
                orig_B->GetAlpha(), orig_B->GetBeta(), dest_B->GetPointsType(),
                dest_B->GetPointsOrder(), dest_B->GetAlpha(), dest_B->GetBeta(), 
                I, m_imvals[(int)orig_B->GetPointsType()]);
        }

        void PolyManager:: ResetI(void)
        {
            std::vector<InterpMatrix*>::iterator definterp;

            for(int i = 0; i < SIZE_PointsType; ++i)
            {
                for(definterp = m_imvals[i].begin(); definterp != m_imvals[i].end(); 
                    ++definterp)
                {
                    delete definterp[0];
                }
                m_imvals[i].erase(m_imvals[i].begin(),m_imvals[i].end());
                m_imvalscur[i] = m_imvals[i].begin();
            }
        }


        void PolyManager::ILookup( const PointsType orig_ptype, const int orig_np, 
            const double orig_alpha, const double orig_beta,
            const PointsType dest_ptype, const int dest_np, 
            const double dest_alpha,  const double dest_beta, 
            const double* &I, std::vector<InterpMatrix*> &val)
        {

            std::vector<InterpMatrix*>::iterator def;

            InterpMatrixKey key(orig_ptype, orig_np, orig_alpha, orig_beta, 
                dest_ptype, dest_np, dest_alpha, dest_beta);

            def = find(val.begin(),val.end(),m_imvalscur[(int)orig_ptype],key);

            if(def != val.end()){
                I = def[0]->GetInterp();
                m_imvalscur[(int)orig_ptype] = def;
            }
            else{
                InterpMatrix * tmp = new InterpMatrix(this,key);
                val.push_back(tmp);
                I = tmp->GetInterp();
                m_imvalscur[(int)orig_ptype] = (--val.end());
            }    
        }

        void InterpMatrix::SetI(PolyManager *P, const PointsType orig_ptype, 
            const int orig_np, const double orig_alpha, 
            const double orig_beta, const PointsType dest_ptype, 
            const int dest_np, const double dest_alpha,
            const double dest_beta)
        {
            m_orig_ptype = orig_ptype;
            m_orig_np    = orig_np;
            m_orig_alpha = orig_alpha;
            m_orig_beta  = orig_beta;
            m_dest_ptype = dest_ptype;
            m_dest_np    = dest_np;
            m_dest_alpha = dest_alpha;
            m_dest_beta  = dest_beta;

            const double *orig_z,*orig_w,*dest_z,*dest_w;

            if(m_orig_np){
                m_im = new double[m_orig_np*m_dest_np];

                // get zeros for calculation of interpolation matrices
                P->GetZW(orig_ptype,orig_np,orig_z,orig_w,orig_alpha,orig_beta);
                P->GetZW(dest_ptype,dest_np,dest_z,dest_w,dest_alpha,dest_beta);

                switch(orig_ptype){
      case eGauss:
          Polylib::Imgj(m_im,orig_z,dest_z,orig_np,dest_np,orig_alpha,orig_beta);
          break;
      case eLobatto:
          Polylib::Imglj(m_im,orig_z,dest_z,orig_np,dest_np,orig_alpha,dest_beta);
          break;
      case eRadauM:
          Polylib::Imgrjm(m_im,orig_z,dest_z,orig_np,dest_np,orig_alpha,dest_beta);
          break;
      case eRadauP:
          Polylib::Imgrjp(m_im,orig_z,dest_z,orig_np,dest_np,orig_alpha,dest_beta);
          break;
      case ePolyEvenSp:
          ASSERTL0(false, "Need to define routine for equispaced polynomial Interpolation");
          break;
      case eFourierEvenSp:
          ASSERTL0(false, "Fourier points Interpolation matrix for Fourier spacing");
          break;
defaults:
          ASSERTL0(false, "Unknown point type distribution");
                }
            }
        }


        bool operator == (const InterpMatrix &x, const InterpMatrixKey &y)
        {
            if( (x.m_orig_np == y.m_orig_np) &&
                (x.m_dest_np == y.m_dest_np) &&
                (x.m_orig_ptype == y.m_orig_ptype) &&
                (x.m_dest_ptype == y.m_dest_ptype) &&
                (x.m_orig_alpha == y.m_orig_alpha) &&
                (x.m_dest_alpha == y.m_dest_alpha) &&
                (x.m_orig_beta == y.m_orig_beta) &&
                (x.m_dest_beta == y.m_dest_beta) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }


        bool operator == (const InterpMatrixKey &x, const InterpMatrix &y)
        {
            if( (x.m_orig_np == y.m_orig_np) &&
                (x.m_dest_np == y.m_dest_np) &&
                (x.m_orig_ptype == y.m_orig_ptype) &&
                (x.m_dest_ptype == y.m_dest_ptype) &&
                (x.m_orig_alpha == y.m_orig_alpha) &&
                (x.m_dest_alpha == y.m_dest_alpha) &&
                (x.m_orig_beta == y.m_orig_beta) &&
                (x.m_dest_beta == y.m_dest_beta) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        bool operator == (const InterpMatrixKey &x, const InterpMatrix *y)
        {
            if( (x.m_orig_np == (*y).m_orig_np) &&
                (x.m_dest_np == (*y).m_dest_np) &&
                (x.m_orig_ptype == (*y).m_orig_ptype) &&
                (x.m_dest_ptype == (*y).m_dest_ptype) &&
                (x.m_orig_alpha == (*y).m_orig_alpha) &&
                (x.m_dest_alpha == (*y).m_dest_alpha) &&
                (x.m_orig_beta == (*y).m_orig_beta) &&
                (x.m_dest_beta == (*y).m_dest_beta) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        bool operator == (const InterpMatrix *x, const InterpMatrixKey &y)
        {
            if( ((*x).m_orig_np    == y.m_orig_np) &&
                ((*x).m_dest_np    == y.m_dest_np) &&
                ((*x).m_orig_ptype == y.m_orig_ptype) &&
                ((*x).m_dest_ptype == y.m_dest_ptype) &&
                ((*x).m_orig_alpha == y.m_orig_alpha) &&
                ((*x).m_dest_alpha == y.m_dest_alpha) &&
                ((*x).m_orig_beta  == y.m_orig_beta) &&
                ((*x).m_dest_beta  == y.m_dest_beta) )
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        bool operator != (const InterpMatrix &x, const InterpMatrixKey &y)
        {
            if( (x.m_orig_np == y.m_orig_np) &&
                (x.m_dest_np == y.m_dest_np) &&
                (x.m_orig_ptype == y.m_orig_ptype) &&
                (x.m_dest_ptype == y.m_dest_ptype) &&
                (x.m_orig_alpha == y.m_orig_alpha) &&
                (x.m_dest_alpha == y.m_dest_alpha) &&
                (x.m_orig_beta == y.m_orig_beta) &&
                (x.m_dest_beta == y.m_dest_beta) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        bool operator != (const InterpMatrixKey &x, const InterpMatrix &y)
        {
            if( (x.m_orig_np == y.m_orig_np) &&
                (x.m_dest_np == y.m_dest_np) &&
                (x.m_orig_ptype == y.m_orig_ptype) &&
                (x.m_dest_ptype == y.m_dest_ptype) &&
                (x.m_orig_alpha == y.m_orig_alpha) &&
                (x.m_dest_alpha == y.m_dest_alpha) &&
                (x.m_orig_beta == y.m_orig_beta) &&
                (x.m_dest_beta == y.m_dest_beta) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        bool operator != (const InterpMatrixKey &x, const InterpMatrix *y)
        {
            if( (x.m_orig_np == (*y).m_orig_np) &&
                (x.m_dest_np == (*y).m_dest_np) &&
                (x.m_orig_ptype == (*y).m_orig_ptype) &&
                (x.m_dest_ptype == (*y).m_dest_ptype) &&
                (x.m_orig_alpha == (*y).m_orig_alpha) &&
                (x.m_dest_alpha == (*y).m_dest_alpha) &&
                (x.m_orig_beta == (*y).m_orig_beta) &&
                (x.m_dest_beta == (*y).m_dest_beta) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }

        bool operator != (const InterpMatrix *x, const InterpMatrixKey &y)
        {
            if( ((*x).m_orig_np == y.m_orig_np) &&
                ((*x).m_dest_np == y.m_dest_np) &&
                ((*x).m_orig_ptype == y.m_orig_ptype) &&
                ((*x).m_dest_ptype == y.m_dest_ptype) &&
                ((*x).m_orig_alpha == y.m_orig_alpha) &&
                ((*x).m_dest_alpha == y.m_dest_alpha) &&
                ((*x).m_orig_beta == y.m_orig_beta) &&
                ((*x).m_dest_beta == y.m_dest_beta) )
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    } // end of namespace 
} // end of namespace 

/** 
* $Log: PolyManager.cpp,v $
* Revision 1.30  2006/04/25 20:23:33  jfrazier
* Various fixes to correct bugs, calls to ASSERT, etc.
*
* Revision 1.29  2006/03/12 14:20:44  sherwin
*
* First compiling version of SpatialDomains and associated modifications
*
* Revision 1.28  2006/03/05 22:11:02  sherwin
*
* Sorted out Project1D, Project2D and Project_Diff2D as well as some test scripts
*
* Revision 1.27  2006/02/27 23:47:23  sherwin
*
* Standard coding update upto compilation of StdHexExp.cpp
*
* Revision 1.26  2006/02/26 23:37:29  sherwin
*
* Updates and compiling checks upto StdExpansions1D
*
* Revision 1.25  2006/02/19 13:26:12  sherwin
*
* Coding standard revisions so that libraries compile
*
* Revision 1.24  2006/02/15 08:06:35  sherwin
*
* Put files into coding standard (although they do not compile)
*
*
**/ 

