///////////////////////////////////////////////////////////////////////////////
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#include <LocalRegions/SegExp.h>
#include <SpatialDomains/MeshGraph.h>
#include <Collections/Collection.h>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace SegCollectionTests
    {
        SpatialDomains::SegGeomSharedPtr CreateSegGeom(unsigned int id, 
                                                       SpatialDomains::PointGeomSharedPtr v0,
                                                       SpatialDomains::PointGeomSharedPtr v1)
        {
            SpatialDomains::PointGeomSharedPtr vertices[] = {v0, v1};
            SpatialDomains::SegGeomSharedPtr result(new SpatialDomains::SegGeom(id, 2, vertices));
            return result;
        }
        
        
        BOOST_AUTO_TEST_CASE(TestSegBwdTrans_StdMat_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            
            SpatialDomains::SegGeomSharedPtr segGeom = CreateSegGeom(0, v0, v1);
            
            Nektar::LibUtilities::PointsType segPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType  basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numSegPoints = 6;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(numSegPoints, segPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);
            
            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1, segGeom);

            Nektar::StdRegions::StdSegExpSharedPtr stdExp = 
                MemoryManager<Nektar::StdRegions::StdSegExp>::AllocateSharedPtr(basisKeyDir1);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            GeomVec.push_back(segGeom);
            
            Collections::Collection c(stdExp, GeomVec,Collections::eStdMat);

            Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(Exp->GetTotPoints());
            

            Exp->BwdTrans(coeffs, phys1);
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }

        
        BOOST_AUTO_TEST_CASE(TestSegBwdTrans_StdMat_UniformP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            
            SpatialDomains::SegGeomSharedPtr segGeom = CreateSegGeom(0, v0, v1);
            
            Nektar::LibUtilities::PointsType segPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numSegPoints = 6;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(numSegPoints, segPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);
            
            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1, segGeom);

            Nektar::StdRegions::StdSegExpSharedPtr stdExp = 
                MemoryManager<Nektar::StdRegions::StdSegExp>::AllocateSharedPtr(basisKeyDir1);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;

            int nelmts = 10;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(segGeom);
            }
            
            Collections::Collection c(stdExp, GeomVec,Collections::eStdMat);


            Array<OneD, NekDouble> coeffs(nelmts*Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(nelmts*Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(nelmts*Exp->GetTotPoints());
            
            for(int i = 0; i < nelmts; ++i)
            {
                Exp->BwdTrans(coeffs + i*Exp->GetNcoeffs(), tmp = phys1+i*Exp->GetTotPoints());
            }
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }
        

        BOOST_AUTO_TEST_CASE(TestSegBwdTrans_IterPerExp_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            
            SpatialDomains::SegGeomSharedPtr segGeom = CreateSegGeom(0, v0, v1);
            
            Nektar::LibUtilities::PointsType segPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numSegPoints = 6;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(numSegPoints, segPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);
            
            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1, segGeom);

            Nektar::StdRegions::StdSegExpSharedPtr stdExp = 
                MemoryManager<Nektar::StdRegions::StdSegExp>::AllocateSharedPtr(basisKeyDir1);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            GeomVec.push_back(segGeom);
            
            Collections::Collection c(stdExp, GeomVec,Collections::eIterPerExp);

            Array<OneD, NekDouble> coeffs(Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(Exp->GetTotPoints());
            

            Exp->BwdTrans(coeffs, phys1);
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }
        
        
        BOOST_AUTO_TEST_CASE(TestSegBwdTrans_SumFac_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            
            SpatialDomains::SegGeomSharedPtr segGeom = CreateSegGeom(0, v0, v1);
            
            Nektar::LibUtilities::PointsType segPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numSegPoints = 6;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(numSegPoints, segPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);
            
            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1, segGeom);

            Nektar::StdRegions::StdSegExpSharedPtr stdExp = 
                MemoryManager<Nektar::StdRegions::StdSegExp>::AllocateSharedPtr(basisKeyDir1);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;

            int nelmts = 1;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(segGeom);
            }
            
            Collections::Collection c(stdExp, GeomVec, Collections::eSumFac);

            Array<OneD, NekDouble> coeffs(nelmts*Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(nelmts*Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(nelmts*Exp->GetTotPoints());
            
            for(int i = 0; i < nelmts; ++i)
            {
                Exp->BwdTrans(coeffs + i*Exp->GetNcoeffs(), tmp = phys1+i*Exp->GetTotPoints());
            }
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestSegBwdTrans_SumFac_UniformP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.0, -1.0, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            
            SpatialDomains::SegGeomSharedPtr segGeom = CreateSegGeom(0, v0, v1);
            
            Nektar::LibUtilities::PointsType segPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numSegPoints = 6;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(numSegPoints, segPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);
            
            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1, segGeom);

            Nektar::StdRegions::StdSegExpSharedPtr stdExp = 
                MemoryManager<Nektar::StdRegions::StdSegExp>::AllocateSharedPtr(basisKeyDir1);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            
            int nelmts = 10;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(segGeom);
            }
            
            Collections::Collection c(stdExp, GeomVec, Collections::eSumFac);

            Array<OneD, NekDouble> coeffs(nelmts*Exp->GetNcoeffs(), 1.0), tmp;
            Array<OneD, NekDouble> phys1(nelmts*Exp->GetTotPoints());
            Array<OneD, NekDouble> phys2(nelmts*Exp->GetTotPoints());
            
            for(int i = 0; i < nelmts; ++i)
            {
                Exp->BwdTrans(coeffs + i*Exp->GetNcoeffs(), tmp = phys1+i*Exp->GetTotPoints());
            }
            c.ApplyOperator(Collections::eBwdTrans, coeffs, phys2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < phys1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(phys1[i],phys2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestSegIProductWRTBase_IterPerExp_UniformP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            
            SpatialDomains::SegGeomSharedPtr segGeom = CreateSegGeom(0, v0, v1);
            
            Nektar::LibUtilities::PointsType segPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(5, segPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);

            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1,
                                                                               segGeom);

            Nektar::StdRegions::StdSegExpSharedPtr stdExp = 
                MemoryManager<Nektar::StdRegions::StdSegExp>::AllocateSharedPtr(basisKeyDir1);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;

            int nelmts = 10;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(segGeom);
            }
            
            Collections::Collection c(stdExp, GeomVec,Collections::eIterPerExp);

            const int nq = Exp->GetTotPoints();
            const int nc = Exp->GetNcoeffs();
            Array<OneD, NekDouble> xc(nq),yc(nq);
            Array<OneD, NekDouble> phys(nelmts*nq),tmp;
            Array<OneD, NekDouble> coeffs1(nelmts*nc);
            Array<OneD, NekDouble> coeffs2(nelmts*nc);
            
            Exp->GetCoords(xc,yc);
        
            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i]);
            }

            for(int i = 0; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
                Exp->IProductWRTBase(phys + i*nq ,tmp = coeffs1 + i*nc);
            }
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestSegIProductWRTBase_StdMat_UniformP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            
            SpatialDomains::SegGeomSharedPtr segGeom = CreateSegGeom(0, v0, v1);
            
            Nektar::LibUtilities::PointsType segPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(5, segPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);

            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1,
                                                                               segGeom);

            Nektar::StdRegions::StdSegExpSharedPtr stdExp = 
                MemoryManager<Nektar::StdRegions::StdSegExp>::AllocateSharedPtr(basisKeyDir1);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;

            int nelmts = 10;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(segGeom);
            }
            
            Collections::Collection c(stdExp, GeomVec,Collections::eStdMat);

            const int nq = Exp->GetTotPoints();
            const int nc = Exp->GetNcoeffs();
            Array<OneD, NekDouble> xc(nq),yc(nq);
            Array<OneD, NekDouble> phys(nelmts*nq),tmp;
            Array<OneD, NekDouble> coeffs1(nelmts*nc);
            Array<OneD, NekDouble> coeffs2(nelmts*nc);
            
            Exp->GetCoords(xc,yc);
        
            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i]);
            }

            for(int i = 0; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
                Exp->IProductWRTBase(phys + i*nq ,tmp = coeffs1 + i*nc);
            }
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }

        BOOST_AUTO_TEST_CASE(TestSegIProductWRTBase_SumFac_UniformP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            
            SpatialDomains::SegGeomSharedPtr segGeom = CreateSegGeom(0, v0, v1);
            
            Nektar::LibUtilities::PointsType segPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(5, segPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);

            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1,
                                                                               segGeom);

            Nektar::StdRegions::StdSegExpSharedPtr stdExp = 
                MemoryManager<Nektar::StdRegions::StdSegExp>::AllocateSharedPtr(basisKeyDir1);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;

            int nelmts = 10;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(segGeom);
            }
            
            Collections::Collection c(stdExp, GeomVec,Collections::eSumFac);

            const int nq = Exp->GetTotPoints();
            const int nc = Exp->GetNcoeffs();
            Array<OneD, NekDouble> xc(nq),yc(nq);
            Array<OneD, NekDouble> phys(nelmts*nq),tmp;
            Array<OneD, NekDouble> coeffs1(nelmts*nc);
            Array<OneD, NekDouble> coeffs2(nelmts*nc);
            
            Exp->GetCoords(xc,yc);
        
            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i]);
            }

            for(int i = 0; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
                Exp->IProductWRTBase(phys + i*nq ,tmp = coeffs1 + i*nc);
            }
            c.ApplyOperator(Collections::eIProductWRTBase, phys, coeffs2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < coeffs1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(coeffs1[i],coeffs2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestSegPhysDeriv_IterPerExp_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            
            SpatialDomains::SegGeomSharedPtr segGeom = CreateSegGeom(0, v0, v1);
            
            Nektar::LibUtilities::PointsType segPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(5, segPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);

            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1,
                                                                               segGeom);

            Nektar::StdRegions::StdSegExpSharedPtr stdExp = 
                MemoryManager<Nektar::StdRegions::StdSegExp>::AllocateSharedPtr(basisKeyDir1);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            GeomVec.push_back(segGeom);
            
            Collections::Collection c(stdExp, GeomVec,Collections::eIterPerExp);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq),yc(nq);
            Array<OneD, NekDouble> phys(nq),tmp;
            Array<OneD, NekDouble> diff1(nq);
            Array<OneD, NekDouble> diff2(nq);
            
            Exp->GetCoords(xc,yc);
        
            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i]);
            }

            Exp->PhysDeriv(phys,diff1);
            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestSegPhysDeriv_IterPerExp_UniformP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            
            SpatialDomains::SegGeomSharedPtr segGeom = CreateSegGeom(0, v0, v1);
            
            Nektar::LibUtilities::PointsType segPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(5, segPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);

            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1,
                                                                               segGeom);

            Nektar::StdRegions::StdSegExpSharedPtr stdExp = 
                MemoryManager<Nektar::StdRegions::StdSegExp>::AllocateSharedPtr(basisKeyDir1);

            int nelmts = 10;
            
            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(segGeom);
            }

            Collections::Collection c(stdExp, GeomVec,Collections::eIterPerExp);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq);
            Array<OneD, NekDouble> phys(nelmts*nq),tmp;
            Array<OneD, NekDouble> diff1(nelmts*nq);
            Array<OneD, NekDouble> diff2(nelmts*nq);
            
            Exp->GetCoords(xc, yc);
        
            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }
            for(int i = 0; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
                Exp->PhysDeriv(phys, tmp = diff1+i*nq);
            }

            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.num_elements(); ++i)
            {
                diff1[i] = (fabs(diff1[i]) < 1e-14)? 0.0: diff1[i];
                diff2[i] = (fabs(diff2[i]) < 1e-14)? 0.0: diff2[i];
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestSegPhysDeriv_StdMat_UniformP)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            
            SpatialDomains::SegGeomSharedPtr segGeom = CreateSegGeom(0, v0, v1);
            
            Nektar::LibUtilities::PointsType segPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(5, segPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);

            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1,
                                                                               segGeom);

            Nektar::StdRegions::StdSegExpSharedPtr stdExp = 
                MemoryManager<Nektar::StdRegions::StdSegExp>::AllocateSharedPtr(basisKeyDir1);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            GeomVec.push_back(segGeom);
            
            Collections::Collection c(stdExp, GeomVec,Collections::eStdMat);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq),yc(nq);
            Array<OneD, NekDouble> phys(nq),tmp;
            Array<OneD, NekDouble> diff1(nq);
            Array<OneD, NekDouble> diff2(nq);
            
            Exp->GetCoords(xc,yc);
        
            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i]);
            }

            Exp->PhysDeriv(phys,diff1);
            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.num_elements(); ++i)
            {
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }


        BOOST_AUTO_TEST_CASE(TestSegPhysDeriv_StdMat_UniformP_MultiElmt)
        {
            SpatialDomains::PointGeomSharedPtr v0(new SpatialDomains::PointGeom(2u, 0u, -1.5, -1.5, 0.0));
            SpatialDomains::PointGeomSharedPtr v1(new SpatialDomains::PointGeom(2u, 1u,  1.0, -1.0, 0.0));
            
            SpatialDomains::SegGeomSharedPtr segGeom = CreateSegGeom(0, v0, v1);
            
            Nektar::LibUtilities::PointsType segPointsTypeDir1 = Nektar::LibUtilities::eGaussLobattoLegendre;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(5, segPointsTypeDir1);
            Nektar::LibUtilities::BasisType       basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            const Nektar::LibUtilities::BasisKey  basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);

            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1,
                                                                               segGeom);

            Nektar::StdRegions::StdSegExpSharedPtr stdExp = 
                MemoryManager<Nektar::StdRegions::StdSegExp>::AllocateSharedPtr(basisKeyDir1);

            int nelmts = 10;
            
            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(segGeom);
            }

            Collections::Collection c(stdExp, GeomVec,Collections::eStdMat);

            const int nq = Exp->GetTotPoints();
            Array<OneD, NekDouble> xc(nq), yc(nq);
            Array<OneD, NekDouble> phys(nelmts*nq),tmp;
            Array<OneD, NekDouble> diff1(nelmts*nq);
            Array<OneD, NekDouble> diff2(nelmts*nq);
            
            Exp->GetCoords(xc, yc);
        
            for (int i = 0; i < nq; ++i)
            {
                phys[i] = sin(xc[i])*cos(yc[i]);
            }
            for(int i = 0; i < nelmts; ++i)
            {
                Vmath::Vcopy(nq,phys,1,tmp = phys+i*nq,1);
                Exp->PhysDeriv(phys, tmp = diff1+i*nq);
            }

            c.ApplyOperator(Collections::ePhysDeriv, phys, diff2);

            double epsilon = 1.0e-8;
            for(int i = 0; i < diff1.num_elements(); ++i)
            {
                diff1[i] = (fabs(diff1[i]) < 1e-14)? 0.0: diff1[i];
                diff2[i] = (fabs(diff2[i]) < 1e-14)? 0.0: diff2[i];
                BOOST_CHECK_CLOSE(diff1[i],diff2[i], epsilon);
            }
        }

    }
}
