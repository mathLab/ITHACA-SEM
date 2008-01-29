#include <Profile/IntWrapper.h>
#include <boost/lexical_cast.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/timer.hpp>
#include <boost/progress.hpp>
#include <iostream>
using std::cout;
using std::endl;



void MatrixAddition2(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1;
	}
	cout << "MatrixAddition2 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition2 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication2(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1;
	}
	cout << "MatrixMultiplication2 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication2 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition2(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1;
	}
	cout << "IntWrapperAddition2 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition2 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixAddition3(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2;
	}
	cout << "MatrixAddition3 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition3 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication3(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2;
	}
	cout << "MatrixMultiplication3 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication3 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition3(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	IntWrapper m2;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1 + m2;
	}
	cout << "IntWrapperAddition3 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition3 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixAddition4(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3;
	}
	cout << "MatrixAddition4 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition4 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication4(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3;
	}
	cout << "MatrixMultiplication4 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication4 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition4(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	IntWrapper m2;
	IntWrapper m3;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1 + m2 + m3;
	}
	cout << "IntWrapperAddition4 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition4 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixAddition5(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4;
	}
	cout << "MatrixAddition5 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition5 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication5(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4;
	}
	cout << "MatrixMultiplication5 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication5 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition5(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	IntWrapper m2;
	IntWrapper m3;
	IntWrapper m4;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1 + m2 + m3 + m4;
	}
	cout << "IntWrapperAddition5 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition5 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixAddition6(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5;
	}
	cout << "MatrixAddition6 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition6 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication6(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5;
	}
	cout << "MatrixMultiplication6 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication6 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition6(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	IntWrapper m2;
	IntWrapper m3;
	IntWrapper m4;
	IntWrapper m5;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1 + m2 + m3 + m4 + m5;
	}
	cout << "IntWrapperAddition6 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition6 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixAddition7(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6;
	}
	cout << "MatrixAddition7 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition7 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication7(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6;
	}
	cout << "MatrixMultiplication7 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication7 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition7(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	IntWrapper m2;
	IntWrapper m3;
	IntWrapper m4;
	IntWrapper m5;
	IntWrapper m6;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1 + m2 + m3 + m4 + m5 + m6;
	}
	cout << "IntWrapperAddition7 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition7 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixAddition8(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7;
	}
	cout << "MatrixAddition8 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition8 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication8(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7;
	}
	cout << "MatrixMultiplication8 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication8 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition8(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	IntWrapper m2;
	IntWrapper m3;
	IntWrapper m4;
	IntWrapper m5;
	IntWrapper m6;
	IntWrapper m7;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7;
	}
	cout << "IntWrapperAddition8 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition8 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixAddition9(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m8(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8;
	}
	cout << "MatrixAddition9 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition9 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication9(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m8(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8;
	}
	cout << "MatrixMultiplication9 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication9 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition9(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	IntWrapper m2;
	IntWrapper m3;
	IntWrapper m4;
	IntWrapper m5;
	IntWrapper m6;
	IntWrapper m7;
	IntWrapper m8;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8;
	}
	cout << "IntWrapperAddition9 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition9 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixAddition10(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m8(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m9(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9;
	}
	cout << "MatrixAddition10 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition10 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication10(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m8(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m9(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9;
	}
	cout << "MatrixMultiplication10 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication10 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition10(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	IntWrapper m2;
	IntWrapper m3;
	IntWrapper m4;
	IntWrapper m5;
	IntWrapper m6;
	IntWrapper m7;
	IntWrapper m8;
	IntWrapper m9;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9;
	}
	cout << "IntWrapperAddition10 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition10 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixAddition11(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m8(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m9(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m10(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10;
	}
	cout << "MatrixAddition11 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition11 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication11(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m8(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m9(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m10(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10;
	}
	cout << "MatrixMultiplication11 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication11 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition11(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	IntWrapper m2;
	IntWrapper m3;
	IntWrapper m4;
	IntWrapper m5;
	IntWrapper m6;
	IntWrapper m7;
	IntWrapper m8;
	IntWrapper m9;
	IntWrapper m10;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10;
	}
	cout << "IntWrapperAddition11 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition11 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixAddition12(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m8(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m9(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m10(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m11(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11;
	}
	cout << "MatrixAddition12 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition12 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication12(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m8(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m9(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m10(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m11(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11;
	}
	cout << "MatrixMultiplication12 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication12 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition12(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	IntWrapper m2;
	IntWrapper m3;
	IntWrapper m4;
	IntWrapper m5;
	IntWrapper m6;
	IntWrapper m7;
	IntWrapper m8;
	IntWrapper m9;
	IntWrapper m10;
	IntWrapper m11;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11;
	}
	cout << "IntWrapperAddition12 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition12 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixAddition13(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m8(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m9(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m10(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m11(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m12(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11 + m12;
	}
	cout << "MatrixAddition13 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition13 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication13(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m8(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m9(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m10(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m11(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m12(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11 + m12;
	}
	cout << "MatrixMultiplication13 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication13 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition13(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	IntWrapper m2;
	IntWrapper m3;
	IntWrapper m4;
	IntWrapper m5;
	IntWrapper m6;
	IntWrapper m7;
	IntWrapper m8;
	IntWrapper m9;
	IntWrapper m10;
	IntWrapper m11;
	IntWrapper m12;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11 + m12;
	}
	cout << "IntWrapperAddition13 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition13 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixAddition14(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m8(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m9(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m10(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m11(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m12(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m13(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11 + m12 + m13;
	}
	cout << "MatrixAddition14 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixAddition14 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void MatrixMultiplication14(unsigned int numIterations, unsigned int matrixSize)
{
	 Nektar::NekMatrix<double> m0(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m1(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m2(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m3(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m4(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m5(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m6(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m7(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m8(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m9(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m10(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m11(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m12(matrixSize, matrixSize, 1.0);
	 Nektar::NekMatrix<double> m13(matrixSize, matrixSize, 1.0);
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		Nektar::NekMatrix<double> result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11 + m12 + m13;
	}
	cout << "MatrixMultiplication14 Total Time: " << t.elapsed() << "\n";
	cout << "MatrixMultiplication14 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void IntWrapperAddition14(unsigned int numIterations, unsigned int matrixSize)
{
	IntWrapper m0;
	IntWrapper m1;
	IntWrapper m2;
	IntWrapper m3;
	IntWrapper m4;
	IntWrapper m5;
	IntWrapper m6;
	IntWrapper m7;
	IntWrapper m8;
	IntWrapper m9;
	IntWrapper m10;
	IntWrapper m11;
	IntWrapper m12;
	IntWrapper m13;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		IntWrapper result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11 + m12 + m13;
	}
	cout << "IntWrapperAddition14 Total Time: " << t.elapsed() << "\n";
	cout << "IntWrapperAddition14 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

int main(int argc, char** argv)
{
	if( argc != 3 ) { cout << "Usage: Profile <numTests> <problemSize>\n"; return 1; }

	unsigned int numTests = boost::lexical_cast<unsigned int>(argv[1]);
	unsigned int problemSize = boost::lexical_cast<unsigned int>(argv[2]);
MatrixAddition2(numTests, problemSize);
MatrixAddition3(numTests, problemSize);
MatrixAddition4(numTests, problemSize);
MatrixAddition5(numTests, problemSize);
MatrixAddition6(numTests, problemSize);
MatrixAddition7(numTests, problemSize);
MatrixAddition8(numTests, problemSize);
MatrixAddition9(numTests, problemSize);
MatrixAddition10(numTests, problemSize);
MatrixAddition11(numTests, problemSize);
MatrixAddition12(numTests, problemSize);
MatrixAddition13(numTests, problemSize);
MatrixAddition14(numTests, problemSize);
}
