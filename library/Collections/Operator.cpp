///////////////////////////////////////////////////////////////////////////////
//
// File: Operator.h
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
// Description: Operator top class definition
//
///////////////////////////////////////////////////////////////////////////////

#include <loki/Singleton.h>
#include <Collections/Operator.h>
#include <Collections/Collection.h>

namespace Nektar {
    namespace Collections {

        bool operator< (OperatorKey const &p1, OperatorKey const &p2)
        {
            if (boost::get<0>(p1) < boost::get<0>(p2))
            {
                return true;
            }
            if (boost::get<0>(p1) > boost::get<0>(p2))
            {
                return false;
            }
            if (boost::get<1>(p1) < boost::get<1>(p2))
            {
                return true;
            }
            if (boost::get<1>(p1) > boost::get<1>(p2))
            {
                return false;
            }
            if (boost::get<2>(p1) < boost::get<2>(p2))
            {
                return true;
            }
            if (boost::get<2>(p1) > boost::get<2>(p2))
            {
                return false;
            }

            if (boost::get<3>(p1) < boost::get<3>(p2))
            {
                return true;
            }
            if (boost::get<3>(p1) > boost::get<3>(p2))
            {
                return false;
            }
            
            return false;
        }
        
        std::ostream &operator<<(std::ostream &os, OperatorKey const &p)
        {
            os <<                       boost::get<0>(p)  << " "
               << OperatorTypeMap      [boost::get<1>(p)] << " "
               << ImplementationTypeMap[boost::get<2>(p)] << " "
               << ImplementationTypeMap[boost::get<3>(p)];
            return os;
        }

        Operator::~Operator(void)
        {
        }

        
        OperatorFactory& GetOperatorFactory()
        {
            typedef Loki::SingletonHolder<OperatorFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy > Type;
            return Type::Instance();
        }
        
        /*
         * ----------------------------------------------------------
         * BwdTrans operators
         * ----------------------------------------------------------
         */
        class BwdTrans_StdMat : public Operator
        {
        public:
            BwdTrans_StdMat(StdRegions::StdExpansionSharedPtr pExp,
                            vector<SpatialDomains::GeometrySharedPtr> pGeom,
                            CoalescedGeomDataSharedPtr GeomData)
                : Operator(pExp, pGeom, GeomData)
            {
                StdRegions::StdMatrixKey  key(StdRegions::eBwdTrans, pExp->DetShapeType(), *pExp);
                m_mat = m_stdExp->GetStdMatrix(key);
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                Blas::Dgemm('N', 'N', m_mat->GetRows(), m_numElmt,
                            m_mat->GetColumns(), 1.0, m_mat->GetRawPtr(),
                            m_mat->GetRows(), input.get(), m_stdExp->GetNcoeffs(),
                            0.0, output.get(), m_stdExp->GetTotPoints());
            }
            
            OPERATOR_CREATE(BwdTrans_StdMat)
            
            DNekMatSharedPtr m_mat;
        };
        
        OperatorKey BwdTrans_StdMat::m_typeArr[] =
            {
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eSegment, eBwdTrans, eStdMat,false),
                      BwdTrans_StdMat::create, "BwdTrans_StdMat_Seg"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, eBwdTrans, eStdMat,false),
                      BwdTrans_StdMat::create, "BwdTrans_StdMat_Tri"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, eBwdTrans, eStdMat,true),
                      BwdTrans_StdMat::create, "BwdTrans_StdMat_NodalTri"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, eBwdTrans, eSumFac,true),
                      BwdTrans_StdMat::create, "BwdTrans_SumFac_NodalTri"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eQuadrilateral, eBwdTrans, eStdMat,false),
                      BwdTrans_StdMat::create, "BwdTrans_StdMat_Quad"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTetrahedron, eBwdTrans, eStdMat,false),
                      BwdTrans_StdMat::create, "BwdTrans_StdMat_Tet"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTetrahedron, eBwdTrans, eStdMat,true),
                      BwdTrans_StdMat::create, "BwdTrans_StdMat_NodalTet"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTetrahedron, eBwdTrans, eSumFac,true),
                      BwdTrans_StdMat::create, "BwdTrans_SumFac_NodalTet"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePyramid, eBwdTrans, eStdMat,false),
                      BwdTrans_StdMat::create, "BwdTrans_StdMat_Pyr"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePyramid, eBwdTrans, eSumFac,false),
                      BwdTrans_StdMat::create, "BwdTrans_SumFac_Pyr"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePrism, eBwdTrans, eStdMat,false),
                      BwdTrans_StdMat::create, "BwdTrans_StdMat_Prism"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePrism, eBwdTrans, eStdMat,true),
                      BwdTrans_StdMat::create, "BwdTrans_StdMat_NodalPrism"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePrism, eBwdTrans, eSumFac,true),
                      BwdTrans_StdMat::create, "BwdTrans_SumFac_NodalPrism"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eHexahedron, eBwdTrans, eStdMat,false),
                      BwdTrans_StdMat::create, "BwdTrans_StdMat_Hex"),
            };
        
        class BwdTrans_IterPerExp : public Operator
        {
        public:
            BwdTrans_IterPerExp(StdRegions::StdExpansionSharedPtr pExp,
                                vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                CoalescedGeomDataSharedPtr GeomData)
                : Operator(pExp, pGeom,GeomData)
            {
            }
            
            virtual void operator()(
                                    const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                const int nCoeffs = m_stdExp->GetNcoeffs();
                const int nPhys   = m_stdExp->GetTotPoints();
                Array<OneD, NekDouble> tmp;
                
                for (int i = 0; i < m_numElmt; ++i)
                {
                    m_stdExp->BwdTrans(input + i*nCoeffs, tmp = output + i*nPhys);
                }
            }
            
            OPERATOR_CREATE(BwdTrans_IterPerExp)
        };
        
        OperatorKey BwdTrans_IterPerExp::m_typeArr[] =
            {
                GetOperatorFactory().RegisterCreatorFunction(
                    OperatorKey(LibUtilities::eSegment, eBwdTrans, eIterPerExp,false),
                    BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Seg"),
                GetOperatorFactory().RegisterCreatorFunction(
                    OperatorKey(LibUtilities::eTriangle, eBwdTrans, eIterPerExp,false),
                    BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Tri"),
                GetOperatorFactory().RegisterCreatorFunction(
                    OperatorKey(LibUtilities::eTriangle, eBwdTrans, eIterPerExp,true),
                    BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_NodalTri"),
                GetOperatorFactory().RegisterCreatorFunction(
                    OperatorKey(LibUtilities::eTriangle, eBwdTrans, eIterPerExp,true),
                    BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_NodalTri"),
                GetOperatorFactory().RegisterCreatorFunction(
                    OperatorKey(LibUtilities::eQuadrilateral, eBwdTrans, eIterPerExp,false),
                    BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Quad"),
                GetOperatorFactory().RegisterCreatorFunction(
                    OperatorKey(LibUtilities::eTetrahedron, eBwdTrans, eIterPerExp,false),
                    BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Tet"),
                GetOperatorFactory().RegisterCreatorFunction(
                    OperatorKey(LibUtilities::eTetrahedron, eBwdTrans, eIterPerExp,true),
                    BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_NodalTet"),
                GetOperatorFactory().RegisterCreatorFunction(
                    OperatorKey(LibUtilities::ePyramid, eBwdTrans, eIterPerExp,false),
                    BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Pyr"),
                GetOperatorFactory().RegisterCreatorFunction(
                    OperatorKey(LibUtilities::ePrism, eBwdTrans, eIterPerExp,false),
                    BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Prism"),
                GetOperatorFactory().RegisterCreatorFunction(
                    OperatorKey(LibUtilities::ePrism, eBwdTrans, eIterPerExp,true),
                    BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_NodalPrism"),
                GetOperatorFactory().RegisterCreatorFunction(
                    OperatorKey(LibUtilities::eHexahedron, eBwdTrans, eIterPerExp,false),
                    BwdTrans_IterPerExp::create, "BwdTrans_IterPerExp_Hex"),
            };
        
        /*
         * ----------------------------------------------------------
         * IProductWRTBase operators
         * ----------------------------------------------------------
         */
        
        class IProductWRTBase_StdMat : public Operator
        {
        public:
            IProductWRTBase_StdMat(StdRegions::StdExpansionSharedPtr pExp,
                                   vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                   CoalescedGeomDataSharedPtr GeomData)
                : Operator(pExp, pGeom,GeomData)
            {
                m_jac = GeomData->GetJac(pExp,pGeom);
                StdRegions::StdMatrixKey key(StdRegions::eIProductWRTBase, pExp->DetShapeType(), *pExp);
                m_mat = m_stdExp->GetStdMatrix(key);
                m_wspSize = pExp->GetTotPoints()*m_numElmt;
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                ASSERTL1(wsp.num_elements() == m_wspSize,
                         "Incorrect workspace size");
                
                Vmath::Vmul(m_jac.num_elements(),m_jac,1,input,1,wsp,1);
                
                Blas::Dgemm('N', 'N', m_mat->GetRows(), m_numElmt,
                            m_mat->GetColumns(), 1.0, m_mat->GetRawPtr(),
                            m_mat->GetRows(), wsp.get(), m_stdExp->GetTotPoints(),
                            0.0, output.get(), m_stdExp->GetNcoeffs());
            }
            
            OPERATOR_CREATE(IProductWRTBase_StdMat)
            
            DNekMatSharedPtr m_mat;
            Array<OneD, const NekDouble> m_jac;
            
        };
        
        OperatorKey IProductWRTBase_StdMat::m_typeArr[] =
            {
                GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eSegment, eIProductWRTBase, eStdMat,false),
            IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Seg"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTriangle, eIProductWRTBase, eStdMat,false),
            IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Tri"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTriangle, eIProductWRTBase, eStdMat,true),
            IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_NodalTri"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTriangle, eIProductWRTBase, eSumFac,true),
            IProductWRTBase_StdMat::create, "IProductWRTBase_SumFac_NodalTri"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eQuadrilateral, eIProductWRTBase, eStdMat,false),
            IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Quad"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTetrahedron, eIProductWRTBase, eStdMat,false),
            IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Tet"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTetrahedron, eIProductWRTBase, eStdMat,true),
            IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_NodalTet"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTetrahedron, eIProductWRTBase, eSumFac,true),
            IProductWRTBase_StdMat::create, "IProductWRTBase_SumFac_NodalTet"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::ePyramid, eIProductWRTBase, eStdMat,false),
            IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Pyr"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::ePyramid, eIProductWRTBase, eSumFac,false),
            IProductWRTBase_StdMat::create, "IProductWRTBase_SumFac_Pyr"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::ePrism, eIProductWRTBase, eStdMat,false),
            IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Prism"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::ePrism, eIProductWRTBase, eStdMat,true),
            IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_NodalPrism"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::ePrism, eIProductWRTBase, eSumFac,true),
            IProductWRTBase_StdMat::create, "IProductWRTBase_SumFac_NodalPrism"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eHexahedron, eIProductWRTBase, eStdMat,false),
            IProductWRTBase_StdMat::create, "IProductWRTBase_StdMat_Hex"),
    };


    class IProductWRTBase_IterPerExp : public Operator
    {
    public:
        IProductWRTBase_IterPerExp(StdRegions::StdExpansionSharedPtr pExp,
                                   vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                   CoalescedGeomDataSharedPtr GeomData)
            : Operator(pExp, pGeom, GeomData)
        {
            int nqtot = 1;
            LibUtilities::PointsKeyVector PtsKey = pExp->GetPointsKeys();
            for(int i = 0; i < PtsKey.size(); ++i)
            {
                nqtot *= PtsKey[i].GetNumPoints();
            }
            
            m_jac = GeomData->GetJacWithStdWeights(pExp,pGeom);
            
            m_wspSize = nqtot*m_numElmt;
        }

        virtual void operator()(const Array<OneD, const NekDouble> &input,
                                Array<OneD,       NekDouble> &output,
                                Array<OneD,       NekDouble> &output1,
                                Array<OneD,       NekDouble> &output2,
                                Array<OneD,       NekDouble> &wsp)
        {
            ASSERTL1(wsp.num_elements() == m_wspSize,
                     "Incorrect workspace size");

            const int nCoeffs = m_stdExp->GetNcoeffs();
            const int nPhys   = m_stdExp->GetTotPoints();
            Array<OneD, NekDouble> tmp;

            Vmath::Vmul(m_jac.num_elements(),m_jac,1,input,1,wsp,1);

            for (int i = 0; i < m_numElmt; ++i)
            {
                m_stdExp->IProductWRTBase_SumFac(wsp + i*nPhys, tmp = output + i*nCoeffs,false);
            }
        }

        OPERATOR_CREATE(IProductWRTBase_IterPerExp)
        
        Array<OneD, NekDouble> m_jac;

    };

    OperatorKey IProductWRTBase_IterPerExp::m_typeArr[] =
    {
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eSegment, eIProductWRTBase, eIterPerExp,false),
            IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Seg"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTriangle, eIProductWRTBase, eIterPerExp,false),
            IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Tri"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTriangle, eIProductWRTBase, eIterPerExp,true),
            IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_NodalTri"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eQuadrilateral, eIProductWRTBase, eIterPerExp,false),
            IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Quad"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTetrahedron, eIProductWRTBase, eIterPerExp,false),
            IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Tet"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eTetrahedron, eIProductWRTBase, eIterPerExp,true),
            IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_NodalTet"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::ePyramid, eIProductWRTBase, eIterPerExp,false),
            IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Pyr"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::ePrism, eIProductWRTBase, eIterPerExp,false),
            IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Prism"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::ePrism, eIProductWRTBase, eIterPerExp,true),
            IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_NodalPrism"),
        GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eHexahedron, eIProductWRTBase, eIterPerExp,false),
            IProductWRTBase_IterPerExp::create, "IProductWRTBase_IterPerExp_Hex"),
    };



        /*
         * ----------------------------------------------------------
         * PhysDeriv operators
         * ----------------------------------------------------------
         */
        
        class PhysDeriv_IterPerExp : public Operator
        {
        public:
            PhysDeriv_IterPerExp(StdRegions::StdExpansionSharedPtr pExp,
                                 vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                 CoalescedGeomDataSharedPtr GeomData)
                : Operator(pExp, pGeom, GeomData)
            {
                int nqtot = 1;
                LibUtilities::PointsKeyVector PtsKey = pExp->GetPointsKeys();
                m_dim = PtsKey.size();
                m_coordim = m_stdExp->GetCoordim();

                for(int i = 0; i < m_dim; ++i)
                {
                    nqtot *= PtsKey[i].GetNumPoints();
                }
                m_derivFac = GeomData->GetDerivFactors(pExp,pGeom);
                m_wspSize = 3*nqtot*m_numElmt;
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output0,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
#if 1
                int nPhys = m_stdExp->GetTotPoints();
                int ntot = m_numElmt*nPhys;
                Array<OneD, NekDouble> tmp0,tmp1,tmp2;
                Array<OneD, Array<OneD, NekDouble> > Diff(3);
                Array<OneD, Array<OneD, NekDouble> > out(3);
                out[0] = output0;  out[1] = output1;  out[2] = output2;

                for(int i = 0; i < m_dim; ++i)
                {
                    Diff[i] = wsp + i*ntot;
                }

                // calculate local derivatives
                for (int i = 0; i < m_numElmt; ++i)
                {
                    m_stdExp->PhysDeriv(input + i*nPhys, tmp0 = Diff[0] + i*nPhys,
                                        tmp1 = Diff[1] + i*nPhys, tmp2 = Diff[2] + i*nPhys);
                }
#if 1
                // calculate full derivative 
                for(int i = 0; i < m_coordim; ++i)
                {
                    Vmath::Vmul(ntot,m_derivFac[i*m_dim],1,Diff[0],1,out[i],1);
                    for(int j = 1; j < m_dim; ++j)
                    {
                        Vmath::Vvtvp (ntot,m_derivFac[i*m_dim+j],1,Diff[j],1, out[i], 1, out[i],1);
                    }
                }
#endif
#endif          
            }
            
            OPERATOR_CREATE(PhysDeriv_IterPerExp)

            Array<TwoD, const NekDouble> m_derivFac;
            int m_dim;
            int m_coordim;
        };

        OperatorKey PhysDeriv_IterPerExp::m_typeArr[] =
            {
                GetOperatorFactory().RegisterCreatorFunction(
            OperatorKey(LibUtilities::eSegment, ePhysDeriv, eIterPerExp,false),
            PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Seg"),
                GetOperatorFactory().RegisterCreatorFunction(
                               OperatorKey(LibUtilities::eTriangle, ePhysDeriv, eIterPerExp,false),
                               PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Tri"),
                GetOperatorFactory().RegisterCreatorFunction(
                               OperatorKey(LibUtilities::eTriangle, ePhysDeriv, eIterPerExp,true),
                               PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_NodalTri"),
                GetOperatorFactory().RegisterCreatorFunction(
                               OperatorKey(LibUtilities::eQuadrilateral, ePhysDeriv, eIterPerExp,false),
                               PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Quad"),
                GetOperatorFactory().RegisterCreatorFunction(
                               OperatorKey(LibUtilities::eTetrahedron, ePhysDeriv, eIterPerExp,false),
                               PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Tet"),
                GetOperatorFactory().RegisterCreatorFunction(
                               OperatorKey(LibUtilities::eTetrahedron, ePhysDeriv, eIterPerExp,true),
                               PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_NodalTet"),
                GetOperatorFactory().RegisterCreatorFunction(
                               OperatorKey(LibUtilities::ePyramid, ePhysDeriv, eIterPerExp,false),
                               PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Pyr"),
                GetOperatorFactory().RegisterCreatorFunction(
                               OperatorKey(LibUtilities::ePrism, ePhysDeriv, eIterPerExp,false),
                               PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Prism"),
                GetOperatorFactory().RegisterCreatorFunction(
                               OperatorKey(LibUtilities::ePrism, ePhysDeriv, eIterPerExp,true),
                               PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_NodalPrism"),
                GetOperatorFactory().RegisterCreatorFunction(
                               OperatorKey(LibUtilities::eHexahedron, ePhysDeriv, eIterPerExp,false),
                               PhysDeriv_IterPerExp::create, "PhysDeriv_IterPerExp_Hex")
            };

        class PhysDeriv_StdMat : public Operator
        {
        public:
            PhysDeriv_StdMat(StdRegions::StdExpansionSharedPtr pExp,
                                 vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                 CoalescedGeomDataSharedPtr GeomData)
                : Operator(pExp, pGeom, GeomData)
            {
                int nqtot = 1;
                LibUtilities::PointsKeyVector PtsKey = pExp->GetPointsKeys();
                m_dim = PtsKey.size();
                m_coordim = m_stdExp->GetCoordim();

                for(int i = 0; i < m_dim; ++i)
                {
                    nqtot *= PtsKey[i].GetNumPoints();
                }
                // set up a PhysDeriv StdMat. 
                m_derivMat = Array<OneD, DNekMatSharedPtr>(m_dim);
                for(int i = 0; i < m_dim; ++i)
                {
                    Array<OneD, NekDouble> tmp(nqtot),tmp1(nqtot);
                    m_derivMat[i] = MemoryManager<DNekMat>::AllocateSharedPtr(nqtot,nqtot);
                    for(int j = 0; j < nqtot; ++j)
                    {
                        Vmath::Zero(nqtot,tmp,1);
                        tmp[j] = 1.0;
                        m_stdExp->PhysDeriv(i,tmp,tmp1);
                        Vmath::Vcopy(nqtot,&tmp1[0],1,&(m_derivMat[i]->GetPtr())[0]+j*nqtot,1);
                    }
                }
                m_derivFac = GeomData->GetDerivFactors(pExp,pGeom);
                m_wspSize = 3*nqtot*m_numElmt;
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &input,
                                    Array<OneD,       NekDouble> &output0,
                                    Array<OneD,       NekDouble> &output1,
                                    Array<OneD,       NekDouble> &output2,
                                    Array<OneD,       NekDouble> &wsp)
            {
                int nPhys = m_stdExp->GetTotPoints();
                int ntot = m_numElmt*nPhys;
                Array<OneD, NekDouble> tmp0,tmp1,tmp2;
                Array<OneD, Array<OneD, NekDouble> > Diff(3);
                Array<OneD, Array<OneD, NekDouble> > out(3);
                out[0] = output0;  out[1] = output1;    out[2] = output2;

                for(int i = 0; i < m_dim; ++i)
                {
                    Diff[i] = wsp + i*ntot;
                }

                // calculate local derivatives
                for(int i = 0; i < m_dim; ++i)
                {
                    Blas::Dgemm('N', 'N', m_derivMat[i]->GetRows(), m_numElmt,
                                m_derivMat[i]->GetColumns(), 1.0, m_derivMat[i]->GetRawPtr(),
                                m_derivMat[i]->GetRows(), input.get(), nPhys,
                                0.0, &Diff[i][0],nPhys);
                }

                // calculate full derivative 
                for(int i = 0; i < m_coordim; ++i)
                {
                    Vmath::Zero(ntot,out[i],1);
                    for(int j = 0; j < m_dim; ++j)
                    {
                        Vmath::Vvtvp (ntot,m_derivFac[i*m_dim+j],1,Diff[j],1, out[i], 1, out[i],1);
                    }
                }
            }
            
            OPERATOR_CREATE(PhysDeriv_StdMat)

            Array<OneD, DNekMatSharedPtr> m_derivMat;
            Array<TwoD, const NekDouble> m_derivFac;
            int m_dim;
            int m_coordim;
        };

        OperatorKey PhysDeriv_StdMat::m_typeArr[] =
            {
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eSegment, ePhysDeriv, eStdMat,false),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Seg"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, ePhysDeriv, eStdMat,false),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Tri"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, ePhysDeriv, eStdMat,true),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_NodalTri"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, ePhysDeriv, eSumFac,true),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_NodalTri"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eQuadrilateral, ePhysDeriv, eStdMat,false),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Quad"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTetrahedron, ePhysDeriv, eStdMat,false),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Tet"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTetrahedron, ePhysDeriv, eStdMat,false),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Tet"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTetrahedron, ePhysDeriv, eStdMat,true),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_NodalTet"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTetrahedron, ePhysDeriv, eSumFac,true),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_NodalTet"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePyramid, ePhysDeriv, eStdMat,false),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Pyr"),
                GetOperatorFactory().RegisterCreatorFunction(
                                                             OperatorKey(LibUtilities::ePyramid, ePhysDeriv, eSumFac,false),
                      PhysDeriv_StdMat::create, "PhysDeriv_SumFac_Pyr"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePrism, ePhysDeriv, eStdMat,false),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Prism"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePrism, ePhysDeriv, eStdMat,true),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_NodalPrism"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePrism, ePhysDeriv, eSumFac,true),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_NodalPrism"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eHexahedron, ePhysDeriv, eStdMat,false),
                      PhysDeriv_StdMat::create, "PhysDeriv_StdMat_Hex")
            };

        
        
        class IProductWRTDerivBase_IterPerExp : public Operator
        {
        public:
            IProductWRTDerivBase_IterPerExp(StdRegions::StdExpansionSharedPtr pExp,
                                 vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                 CoalescedGeomDataSharedPtr GeomData)
                : Operator(pExp, pGeom, GeomData)
            {
                LibUtilities::PointsKeyVector PtsKey = pExp->GetPointsKeys();
                m_dim = PtsKey.size();
                m_coordim = m_stdExp->GetCoordim();
                
                int nqtot  = pExp->GetTotPoints(); 

                m_derivFac = GeomData->GetDerivFactors(pExp,pGeom);
                m_jac = GeomData->GetJac(pExp,pGeom);
                m_wspSize = m_dim*nqtot*m_numElmt;
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &entry0,
                                    Array<OneD, NekDouble> &entry1,
                                    Array<OneD, NekDouble> &entry2,
                                    Array<OneD, NekDouble> &entry3,
                                    Array<OneD, NekDouble> &wsp)
            {
                unsigned int nPhys  = m_stdExp->GetTotPoints();
                unsigned int ntot   = m_numElmt*nPhys;
                unsigned int nmodes = m_stdExp->GetNcoeffs();
                unsigned int nmax   = max(ntot,m_numElmt*nmodes);
                Array<OneD, Array<OneD, const NekDouble> > in(3);
                Array<OneD, NekDouble> output, tmp1;
                Array<OneD, Array<OneD, NekDouble> > tmp(3);

                in[0] = entry0; in[1] = entry1; in[2] = entry2; 

                output = (m_coordim == 3)? entry3: (m_coordim == 2)?
                    entry2: entry1;

                for(int i = 0; i < m_dim; ++i)
                {
                    tmp[i] = wsp + i*nmax; 
                }
                
                // calculate dx/dxi in[0] + dy/dxi in[2] + dz/dxi in[3]
                for(int i = 0; i < m_dim; ++i)
                {
                    Vmath::Vmul (ntot,m_derivFac[i],1, in[0],1, 
                                 tmp[i],1);
                    for(int j = 1; j < m_coordim; ++j)
                    {
                        Vmath::Vvtvp (ntot,m_derivFac[i +j*m_dim],1,
                                      in[j],1, tmp[i], 1, tmp[i],1);
                    }
                }

                // calculate Iproduct WRT Std Deriv
                // first component
                Vmath::Vmul(ntot,m_jac,1,tmp[0],1,tmp[0],1);
                for(int n = 0; n < m_numElmt; ++n)
                {
                    m_stdExp->IProductWRTDerivBase(0,tmp[0]+n*nPhys,
                                                   tmp1 = output + n*nmodes);
                }

                // other components
                for(int i = 1; i < m_dim; ++i)
                {
                    // multiply by Jacobian
                    Vmath::Vmul(ntot,m_jac,1,tmp[i],1,tmp[i],1);
                    for(int n = 0; n < m_numElmt; ++n)
                    {
                        m_stdExp->IProductWRTDerivBase(i,tmp[i]+n*nPhys,tmp[0]);
                        Vmath::Vadd(nmodes,tmp[0],1,output+n*nmodes,1,
                                    tmp1 = output+n*nmodes,1);
                    }
                }
            }
            
            OPERATOR_CREATE(IProductWRTDerivBase_IterPerExp)

            Array<TwoD, const NekDouble>  m_derivFac;
            Array<OneD, const NekDouble> m_jac;
            int m_dim;
            int m_coordim;
        };

        OperatorKey IProductWRTDerivBase_IterPerExp::m_typeArr[] =
            {
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eSegment, 
                                  eIProductWRTDerivBase, eIterPerExp,false),
                      IProductWRTDerivBase_IterPerExp::create, 
                      "IProductWRTDerivBase_IterPerExp_Seg"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, 
                                  eIProductWRTDerivBase, eIterPerExp,false),
                      IProductWRTDerivBase_IterPerExp::create, 
                      "IProductWRTDerivBase_IterPerExp_Tri"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, 
                                  eIProductWRTDerivBase, eIterPerExp,true),
                      IProductWRTDerivBase_IterPerExp::create, 
                      "IProductWRTDerivBase_IterPerExp_NodalTri"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eQuadrilateral, 
                                  eIProductWRTDerivBase, eIterPerExp,false),
                      IProductWRTDerivBase_IterPerExp::create, 
                      "IProductWRTDerivBase_IterPerExp_Quad"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTetrahedron, 
                                  eIProductWRTDerivBase, eIterPerExp,false),
                      IProductWRTDerivBase_IterPerExp::create, 
                      "IProductWRTDerivBase_IterPerExp_Tet"), 
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTetrahedron, 
                                  eIProductWRTDerivBase, eIterPerExp,true),
                      IProductWRTDerivBase_IterPerExp::create, 
                      "IProductWRTDerivBase_IterPerExp_NodalTet"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePyramid, 
                                  eIProductWRTDerivBase, eIterPerExp,false),
                      IProductWRTDerivBase_IterPerExp::create, 
                      "IProductWRTDerivBase_IterPerExp_Pyr"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePrism, 
                                  eIProductWRTDerivBase, eIterPerExp,false),
                      IProductWRTDerivBase_IterPerExp::create, 
                      "IProductWRTDerivBase_IterPerExp_Prism"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePrism, 
                                  eIProductWRTDerivBase, eIterPerExp,true),
                      IProductWRTDerivBase_IterPerExp::create, 
                      "IProductWRTDerivBase_IterPerExp_NodalPrism"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eHexahedron, 
                                  eIProductWRTDerivBase, eIterPerExp,false),
                      IProductWRTDerivBase_IterPerExp::create, 
                      "IProductWRTDerivBase_IterPerExp_Hex")
            };
        
        class IProductWRTDerivBase_StdMat : public Operator
        {
        public:
            IProductWRTDerivBase_StdMat(StdRegions::StdExpansionSharedPtr pExp,
                                 vector<SpatialDomains::GeometrySharedPtr> pGeom,
                                 CoalescedGeomDataSharedPtr GeomData)
                : Operator(pExp, pGeom, GeomData)
            {
                LibUtilities::PointsKeyVector PtsKey = pExp->GetPointsKeys();
                m_dim = PtsKey.size();
                m_coordim = m_stdExp->GetCoordim();
                
                int nqtot  = pExp->GetTotPoints(); 
                int nmodes = pExp->GetNcoeffs(); 

                // set up a IProductWRTDerivBase StdMat. 
                m_iProdWRTStdDBase = Array<OneD, DNekMatSharedPtr>(m_dim);
                for(int i = 0; i < m_dim; ++i)
                {
                    Array<OneD, NekDouble> tmp(nqtot),tmp1(nmodes);
                    m_iProdWRTStdDBase[i] = MemoryManager<DNekMat>::AllocateSharedPtr(nmodes,nqtot);
                    for(int j = 0; j < nqtot; ++j)
                    {
                        Vmath::Zero(nqtot,tmp,1);
                        tmp[j] = 1.0;
                        m_stdExp->IProductWRTDerivBase(i,tmp,tmp1);
                        Vmath::Vcopy(nmodes,&tmp1[0],1,
                                     &(m_iProdWRTStdDBase[i]->GetPtr())[0]+j*nmodes,1);
                    }
                }
                m_derivFac = GeomData->GetDerivFactors(pExp,pGeom);
                m_jac      = GeomData->GetJac(pExp,pGeom);
                m_wspSize = m_dim*nqtot*m_numElmt;
            }
            
            virtual void operator()(const Array<OneD, const NekDouble> &entry0,
                                    Array<OneD,       NekDouble> &entry1,
                                    Array<OneD,       NekDouble> &entry2,
                                    Array<OneD,       NekDouble> &entry3,
                                    Array<OneD,       NekDouble> &wsp)
            {
                int nPhys = m_stdExp->GetTotPoints();
                int ntot = m_numElmt*nPhys;
                int nmodes = m_stdExp->GetNcoeffs();
                Array<OneD, Array<OneD, const NekDouble> > in(3);
                Array<OneD, NekDouble> output;
                Array<OneD, Array<OneD, NekDouble> > tmp(3);

                in[0] = entry0; in[1] = entry1; 
                in[2] = entry2; 
                
                output = (m_coordim == 3)? entry3: (m_coordim == 2)?
                    entry2: entry1;

                for(int i = 0; i < m_dim; ++i)
                {
                    tmp[i] = wsp + i*ntot; 
                }
                
                // calculate dx/dxi in[0] + dy/dxi in[1] + dz/dxi in[2]
                for(int i = 0; i < m_dim; ++i)
                {
                    Vmath::Vmul (ntot,m_derivFac[i],1, in[0],1, 
                                 tmp[i],1);
                    for(int j = 1; j < m_coordim; ++j)
                    {
                        Vmath::Vvtvp (ntot,m_derivFac[i +j*m_dim],1,
                                      in[j],1, tmp[i], 1, tmp[i],1);
                    }
                }

                // calculate Iproduct WRT Std Deriv

                // First component
                Vmath::Vmul(ntot,m_jac,1,tmp[0],1,tmp[0],1);
                Blas::Dgemm('N', 'N', m_iProdWRTStdDBase[0]->GetRows(), 
                            m_numElmt,m_iProdWRTStdDBase[0]->GetColumns(), 
                            1.0, m_iProdWRTStdDBase[0]->GetRawPtr(),
                            m_iProdWRTStdDBase[0]->GetRows(), 
                            tmp[0].get(), nPhys, 0.0,
                            output.get(), nmodes);

                // Other components
                for(int i = 1; i < m_dim; ++i)
                {
                    Vmath::Vmul(ntot,m_jac,1,tmp[i],1,tmp[i],1);
                    Blas::Dgemm('N', 'N', m_iProdWRTStdDBase[i]->GetRows(), 
                                m_numElmt,m_iProdWRTStdDBase[i]->GetColumns(), 
                                1.0, m_iProdWRTStdDBase[i]->GetRawPtr(),
                                m_iProdWRTStdDBase[i]->GetRows(), 
                                tmp[i].get(), nPhys, 1.0,
                                output.get(), nmodes);
                }
            }
            
            OPERATOR_CREATE(IProductWRTDerivBase_StdMat)

            Array<OneD, DNekMatSharedPtr> m_iProdWRTStdDBase;
            Array<TwoD, const NekDouble>  m_derivFac;
            Array<OneD, const NekDouble>  m_jac;
            int m_dim;
            int m_coordim;
        };

        OperatorKey IProductWRTDerivBase_StdMat::m_typeArr[] =
            {
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eSegment, 
                                  eIProductWRTDerivBase, eStdMat,false),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_StdMat_Seg"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, 
                                  eIProductWRTDerivBase, eStdMat,false),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_StdMat_Tri"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, 
                                  eIProductWRTDerivBase, eStdMat,true),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_StdMat_NodalTri"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTriangle, 
                                  eIProductWRTDerivBase, eSumFac,true),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_SumFac_NodalTri"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eQuadrilateral, 
                                  eIProductWRTDerivBase, eStdMat,false),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_StdMat_Quad"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTetrahedron, 
                                  eIProductWRTDerivBase, eStdMat,false),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_StdMat_Tet"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTetrahedron, 
                                  eIProductWRTDerivBase, eStdMat,true),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_StdMat_NodalTet"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eTetrahedron, 
                                  eIProductWRTDerivBase, eSumFac,true),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_SumFac_NodalTet"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePyramid, 
                                  eIProductWRTDerivBase, eStdMat,false),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_StdMat_Pyr"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePyramid, 
                                  eIProductWRTDerivBase, eSumFac,false),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_SumFac_Pyr"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePrism, 
                                  eIProductWRTDerivBase, eStdMat,false),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_StdMat_Prism"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePrism, 
                                  eIProductWRTDerivBase, eStdMat,true),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_StdMat_NodalPrism"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::ePrism, 
                                  eIProductWRTDerivBase, eSumFac,true),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_SumFac_NodalPrism"),
                GetOperatorFactory().RegisterCreatorFunction(
                      OperatorKey(LibUtilities::eHexahedron, 
                                  eIProductWRTDerivBase, eStdMat,false),
                      IProductWRTDerivBase_StdMat::create, 
                      "IProductWRTDerivBase_StdMat_Hex")
            };
        
        // simple operator Map evaluation
        OperatorImpMap SetFixedImpType(ImplementationType defaultType)
        {
            OperatorImpMap opMap;
            
            for(int i = 0; i < SIZE_OperatorType; ++i)
            {
                opMap[(OperatorType)i] = defaultType;
            }
            
            return opMap;
        }
       
    }
}
