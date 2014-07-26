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
#include <Collections/collection.h>
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
            Nektar::LibUtilities::BasisType basisTypeDir1 = Nektar::LibUtilities::eModified_A;
            unsigned int numSegPoints = 6;
            const Nektar::LibUtilities::PointsKey segPointsKeyDir1(numSegPoints, segPointsTypeDir1);
            const Nektar::LibUtilities::BasisKey basisKeyDir1(basisTypeDir1,4,segPointsKeyDir1);
            
            Nektar::LocalRegions::SegExpSharedPtr Exp = 
                MemoryManager<Nektar::LocalRegions::SegExp>::AllocateSharedPtr(basisKeyDir1, segGeom);

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            GeomVec.push_back(segGeom);
            
            Collections::Collection c(Exp, GeomVec,Collections::eStdMat);

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

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            GeomVec.push_back(segGeom);
            
            Collections::Collection c(Exp, GeomVec,Collections::eIterPerExp);

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

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;

            int nelmts = 1;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(segGeom);
            }
            
            Collections::Collection c(Exp, GeomVec, Collections::eSumFac);

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

            std::vector<SpatialDomains::GeometrySharedPtr> GeomVec;
            
            int nelmts = 10;
            for(int i = 0; i < nelmts; ++i)
            {
                GeomVec.push_back(segGeom);
            }
            
            Collections::Collection c(Exp, GeomVec, Collections::eSumFac);

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
    }
}
