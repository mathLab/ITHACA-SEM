///////////////////////////////////////////////////////////////////////////////
//
// File: TestConsistentObjectAccess.cpp
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
// Description: 
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/BasicUtils/ConsistentObjectAccess.hpp>
#include <UnitTests/util.h>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/test/auto_unit_test.hpp>

namespace Nektar
{
    namespace ConsistentObjectAccessUnitTests
    {
        BOOST_AUTO_TEST_CASE(TestValueTypes)
        {
            double d1 = 1.0;
            const double d2 = 2.0;
            
            BOOST_CHECK_EQUAL(d1, ConsistentObjectAccess<double>::reference(d1));
            BOOST_CHECK_EQUAL(&d1, &ConsistentObjectAccess<double>::reference(d1));
            BOOST_CHECK_EQUAL(d1, ConsistentObjectAccess<double>::const_reference(d1));
            BOOST_CHECK_EQUAL(&d1, &ConsistentObjectAccess<double>::const_reference(d1));
            BOOST_CHECK_EQUAL(&d1, ConsistentObjectAccess<double>::pointer(d1));
            BOOST_CHECK_EQUAL(&d1, ConsistentObjectAccess<double>::const_pointer(d1));
            BOOST_CHECK(ConsistentObjectAccess<double>::ReferencesObject(d1));
            
            BOOST_CHECK_EQUAL(d2, ConsistentObjectAccess<const double>::reference(d2));
            BOOST_CHECK_EQUAL(&d2, &ConsistentObjectAccess<const double>::reference(d2));
            BOOST_CHECK_EQUAL(d2, ConsistentObjectAccess<const double>::const_reference(d2));
            BOOST_CHECK_EQUAL(&d2, &ConsistentObjectAccess<const double>::const_reference(d2));
            BOOST_CHECK_EQUAL(&d2, ConsistentObjectAccess<const double>::pointer(d2));
            BOOST_CHECK_EQUAL(&d2, ConsistentObjectAccess<const double>::const_pointer(d2));
            BOOST_CHECK(ConsistentObjectAccess<const double>::ReferencesObject(d2));
        }
        
        BOOST_AUTO_TEST_CASE(TestPointerTypes)
        {
            UnitTests::RedirectCerrIfNeeded();
            double* d1 = new double(1.0);
            const double* d2 = new double(2.0);
            
            BOOST_CHECK_EQUAL(*d1, ConsistentObjectAccess<double*>::reference(d1));
            BOOST_CHECK_EQUAL(d1, &ConsistentObjectAccess<double*>::reference(d1));
            BOOST_CHECK_EQUAL(*d1, ConsistentObjectAccess<double*>::const_reference(d1));
            BOOST_CHECK_EQUAL(d1, &ConsistentObjectAccess<double*>::const_reference(d1));
            BOOST_CHECK_EQUAL(d1, ConsistentObjectAccess<double*>::pointer(d1));
            BOOST_CHECK_EQUAL(d1, ConsistentObjectAccess<double*>::const_pointer(d1));
            BOOST_CHECK(ConsistentObjectAccess<double*>::ReferencesObject(d1));
            BOOST_CHECK(!ConsistentObjectAccess<double*>::ReferencesObject(static_cast<double*>(0)));
            
            BOOST_CHECK_EQUAL(*d2, ConsistentObjectAccess<const double*>::reference(d2));
            BOOST_CHECK_EQUAL(d2, &ConsistentObjectAccess<const double*>::reference(d2));
            BOOST_CHECK_EQUAL(*d2, ConsistentObjectAccess<const double*>::const_reference(d2));
            BOOST_CHECK_EQUAL(d2, &ConsistentObjectAccess<const double*>::const_reference(d2));
            BOOST_CHECK_EQUAL(d2, ConsistentObjectAccess<const double*>::pointer(d2));
            BOOST_CHECK_EQUAL(d2, ConsistentObjectAccess<const double*>::const_pointer(d2));
            BOOST_CHECK(ConsistentObjectAccess<const double*>::ReferencesObject(d2));
            BOOST_CHECK(!ConsistentObjectAccess<const double*>::ReferencesObject(static_cast<const double*>(0)));
            
            #if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
                BOOST_CHECK_THROW(ConsistentObjectAccess<double*>::const_reference(static_cast<double*>(0)), ErrorUtil::NekError);
                BOOST_CHECK_THROW(ConsistentObjectAccess<const double*>::const_reference(static_cast<const double*>(0)), ErrorUtil::NekError);
                BOOST_CHECK_THROW(ConsistentObjectAccess<double*>::reference(static_cast<double*>(0)), ErrorUtil::NekError);
                BOOST_CHECK_THROW(ConsistentObjectAccess<const double*>::reference(static_cast<const double*>(0)), ErrorUtil::NekError);
            #endif
            
        }
        
        BOOST_AUTO_TEST_CASE(TestSharedPointerTypes)
        {
            double* d1 = new double(1.0);
            const double* d2 = new double(2.0);
            
            std::shared_ptr<double> sd1(d1);
            std::shared_ptr<const double> sd2(d2);
            
            BOOST_CHECK_EQUAL(*d1, ConsistentObjectAccess<std::shared_ptr<double> >::reference(sd1));
            BOOST_CHECK_EQUAL(d1, &ConsistentObjectAccess<std::shared_ptr<double> >::reference(sd1));
            BOOST_CHECK_EQUAL(*d1, ConsistentObjectAccess<std::shared_ptr<double> >::const_reference(sd1));
            BOOST_CHECK_EQUAL(d1, &ConsistentObjectAccess<std::shared_ptr<double> >::const_reference(sd1));
            BOOST_CHECK_EQUAL(d1, ConsistentObjectAccess<std::shared_ptr<double> >::pointer(sd1));
            BOOST_CHECK_EQUAL(d1, ConsistentObjectAccess<std::shared_ptr<double> >::const_pointer(sd1));
            BOOST_CHECK(ConsistentObjectAccess<std::shared_ptr<double> >::ReferencesObject(sd1));
            BOOST_CHECK(!ConsistentObjectAccess<std::shared_ptr<double> >::ReferencesObject(std::shared_ptr<double>()));
            
            BOOST_CHECK_EQUAL(*d2, ConsistentObjectAccess<std::shared_ptr<const double> >::reference(sd2));
            BOOST_CHECK_EQUAL(d2, &ConsistentObjectAccess<std::shared_ptr<const double> >::reference(sd2));
            BOOST_CHECK_EQUAL(*d2, ConsistentObjectAccess<std::shared_ptr<const double> >::const_reference(sd2));
            BOOST_CHECK_EQUAL(d2, &ConsistentObjectAccess<std::shared_ptr<const double> >::const_reference(sd2));
            BOOST_CHECK_EQUAL(d2, ConsistentObjectAccess<std::shared_ptr<const double> >::pointer(sd2));
            BOOST_CHECK_EQUAL(d2, ConsistentObjectAccess<std::shared_ptr<const double> >::const_pointer(sd2));
            BOOST_CHECK(ConsistentObjectAccess<std::shared_ptr<const double> >::ReferencesObject(sd2));
            BOOST_CHECK(!ConsistentObjectAccess<std::shared_ptr<const double> >::ReferencesObject(std::shared_ptr<const double>()));
            
            #if defined(NEKTAR_DEBUG) || defined(NEKTAR_FULLDEBUG)
                BOOST_CHECK_THROW(ConsistentObjectAccess<std::shared_ptr<double> >::const_reference(std::shared_ptr<double>()), ErrorUtil::NekError);
                BOOST_CHECK_THROW(ConsistentObjectAccess<std::shared_ptr<const double> >::const_reference(std::shared_ptr<const double>()), ErrorUtil::NekError);
                BOOST_CHECK_THROW(ConsistentObjectAccess<std::shared_ptr<double> >::reference(std::shared_ptr<double>()), ErrorUtil::NekError);
                BOOST_CHECK_THROW(ConsistentObjectAccess<std::shared_ptr<const double> >::reference(std::shared_ptr<const double>()), ErrorUtil::NekError);
            #endif
            
        }

        
    }
}
