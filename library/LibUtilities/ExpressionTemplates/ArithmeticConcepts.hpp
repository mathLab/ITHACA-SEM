//////////////////////////////////////////////////////////////////////////////////
//// 
//// ArithmeticConcepts.hpp
//// Blake Nelson
////
//// Defines several arithmetic concepts classes for use with the boost 
//// concept checking library.
////
//////////////////////////////////////////////////////////////////////////////////
//
//#ifndef NEKTAR_LIB_UTILITIES_ARITHMETIC_CONCEPT_HPP
//#define NEKTAR_LIB_UTILITIES_ARITHMETIC_CONCEPT_HPP
//
//namespace Nektar
//{
//    template<typename DataType>
//    class AddableConcept
//    {
//        public:
//            void constraints()
//            {
//                // Use pointers so we can check objects without 
//                // default constructors.
//                DataType* lhs;
//                DataType* rhs;
//
//                // Don't assign the result to anything in case operator+
//                // returns void or something weird.
//                *lhs + *rhs;
//            }
//    };
//
//    template<typename DataType>
//    class PlusEqualConcept
//    {
//        public:
//            void constraints()
//            {
//                DataType* lhs;
//                DataType* rhs;
//
//                *lhs += *rhs;
//            }
//    };
//}
//
//
//#endif
//
///**
//    $Log: ArithmeticConcepts.hpp,v $
//    Revision 1.1  2006/06/01 09:20:55  kirby
//    *** empty log message ***
//
//    Revision 1.1  2006/05/04 18:57:40  kirby
//    *** empty log message ***
//
//    Revision 1.1  2006/01/31 13:51:12  bnelson
//    Updated for new configure.
//
//**/
//
//
