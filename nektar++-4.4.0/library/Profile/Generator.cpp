///////////////////////////////////////////////////////////////////////////////
//
// File: Generator.cpp
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
// Description: 
//
// Generates profile tests for expression templates.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <boost/bind.hpp>
#include <string>

enum TestType
{
    Matrix,
    Vector,
    IntWrapper
};

void InitializeMatrix(std::ostream& outFile, unsigned int i)
{
    outFile << "\t Nektar::NekMatrix<double> m" << i << "(matrixSize, matrixSize, 1.0);\n";
}
    
void InitializeIntWrapper(std::ostream& outFile, unsigned int i)
{
    outFile << "\tIntWrapper m" << i << ";\n";
}

void InitializeVariableCostObject(std::ostream& outFile, unsigned int i)
{
    outFile << "\tVariableCostObject m" << i << ";\n";
}

////    boost::timer t;
////    double elapsedTime = 0.0;
////    
////    Nektar::NekMatrix<double> strVals[] = { Nektar::NekMatrix<double>(stringSize, stringSize), 
////        Nektar::NekMatrix<double>(stringSize, stringSize),
////        Nektar::NekMatrix<double>(stringSize, stringSize)};
////        
////    Nektar::NekMatrix<double> result(stringSize, stringSize);
////    
////    cout << "3 Matrix Addition Tests" << endl;
////    cout << "---------------------" << endl;
////    t.restart();
////    for(unsigned int i = 0; i < numTests; ++i)
////    {
////       AddMatrices(result, strVals[0], strVals[1], strVals[2]);    
////    }
////    elapsedTime = t.elapsed();
////    cout << "Straight Addition - Total : " << elapsedTime << endl;
////    cout << "Straight Addition - PerOp : " << elapsedTime/(double)numTests << endl;

template<typename Initializer>
void GenerateMatrixTest(std::ostream& outFile, const std::string& testName, unsigned int opCount, const std::string& opType,
    const std::string& dataTypeName,
    const Initializer& f)
{
    // Create matrix addition and matrix multiplication.
    outFile << "void ";
    outFile << testName;
    outFile << opCount;
    outFile << "(unsigned int numIterations, unsigned int matrixSize)\n";
    outFile << "{\n";
    for(unsigned int i = 0; i < opCount; ++i)
    {
        f(outFile, i);
    }
    
    outFile << "\tboost::timer t;\n";
    outFile << "\n";
    outFile << "\tfor(unsigned int i = 0; i < numIterations; ++i)\n";
    outFile << "\t{\n";
    outFile << "\t\t" << dataTypeName << " result = ";
    for(unsigned int i = 0; i < opCount; ++i)
    {
        outFile << "m" << i;
        if( i != opCount-1 )
        {
            outFile << " + ";
        }
    }
    
    outFile << ";\n";
    outFile << "\t}\n";
    //outFile << "\tcout << \"" << testName << opCount << " Total Time: \" << t.elapsed() << \"\\n\";\n";
    //outFile << "\tcout << \"" << testName << opCount << " Per Op: \" << t.elapsed()/(double)numIterations << \"\\n\";\n";
    outFile << "\tcout << t.elapsed()/(double)numIterations << \"\\t\";\n";
    outFile << "}\n\n";
}

int main(int argc, char** argv)
{
    std::ofstream outFile("g:/Nektar++/library/Profile/main.cpp");
    
    outFile << "#include <Profile/IntWrapper.h>\n";
    outFile << "#include <boost/lexical_cast.hpp>\n";
    outFile << "#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>\n";
    outFile << "#include <Profile/VariableCostObject.h>\n";
//    outFile << "#include <Profile/StringConcat.h>
//    outFile << "#include <Profile/StringConcatExprTemp.h>
    outFile << "#include <boost/timer.hpp>\n";
    outFile << "#include <boost/progress.hpp>\n";
    outFile << "#include <iostream>\n";
    outFile << "using std::cout;\n";
    outFile << "using std::endl;\n";
    outFile << "\n\n\n";
    
    unsigned int maxOpCount = 15;
    for(unsigned int opCount = 2; opCount < maxOpCount; ++opCount)
    {
        GenerateMatrixTest(outFile, "MatrixAddition", opCount, "+", "Nektar::NekMatrix<double>", boost::bind(InitializeMatrix, _1, _2));
        //GenerateMatrixTest(outFile, "MatrixMultiplication", opCount, "*", "Nektar::NekMatrix<double>", boost::bind(InitializeMatrix, _1, _2));
        //GenerateMatrixTest(outFile, "IntWrapperAddition", opCount, "+", "IntWrapper", boost::bind(InitializeIntWrapper, _1, _2));
        //GenerateMatrixTest(outFile, "VariableCostOperationAddition", opCount, "+", "VariableCostObject", boost::bind(InitializeVariableCostObject, _1, _2));
    }
    
    
    outFile << "int main(int argc, char** argv)\n";
    outFile << "{\n";
    outFile << "\tif( argc != 3 ) { cout << \"Usage: Profile <numTests> <problemSize>\\n\"; return 1; }\n\n";
    outFile << "\tunsigned int numTests = boost::lexical_cast<unsigned int>(argv[1]);\n";
    outFile << "\tunsigned int problemSize = boost::lexical_cast<unsigned int>(argv[2]);\n";
    outFile << "\tVariableCostObject::size = problemSize;\n";
    for(unsigned int opCount = 2; opCount < maxOpCount; ++opCount)
    {
        outFile << "MatrixAddition" << opCount << "(numTests, problemSize);\n";
        //outFile << "MatrixMultiplication" << opCount << "(numTests, problemSize);\n";
        //outFile << "IntWrapperAddition" << opCount << "(numTests, problemSize);\n";
        //outFile << "VariableCostOperationAddition" << opCount << "(numTests, problemSize);\n";
    }
    outFile << "}\n";
    
    outFile.close();
    
    return 0;
}


