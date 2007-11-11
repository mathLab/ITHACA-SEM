#include <boost/lexical_cast.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>


int main(int argc, char** argv)
{
    if( argc != 3 )
    {
        std::cerr << "Usage: perf <ProblemSize> <NumTests>" << std::endl;
        return 1;
    }

    unsigned int n = boost::lexical_cast<unsigned int>(argv[1]);
    unsigned int numTests = boost::lexical_cast<unsigned int>(argv[2]);

    unsigned int* buf1 = new unsigned int[n*n];
    unsigned int* buf2 = new unsigned int[n*n];
    
    Nektar::NekMatrix<unsigned int> m1(n, n, buf1);
    Nektar::NekMatrix<unsigned int> m2(n, n, buf2);
    Nektar::NekMatrix<unsigned int> result(n, n);
    
    for(unsigned int i = 0; i < numTests; ++i)
    {
        result = m1*m2;
    }
    
    return 0;
}
