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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
    cout << t.elapsed() << endl;
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
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
	cout << t.elapsed()/(double)numIterations << "\t";
}

int main(int argc, char** argv)
{
	if( argc != 3 ) { cout << "Usage: Profile <numTests> <problemSize>\n"; return 1; }
    cout.precision(20);
	unsigned int numTests = boost::lexical_cast<unsigned int>(argv[1]);
	unsigned int problemSize = boost::lexical_cast<unsigned int>(argv[2]);
// MatrixAddition2(numTests, problemSize);
// MatrixAddition3(numTests, problemSize);
// MatrixAddition4(numTests, problemSize);
// MatrixAddition5(numTests, problemSize);
// MatrixAddition6(numTests, problemSize);
// MatrixAddition7(numTests, problemSize);
// MatrixAddition8(numTests, problemSize);
// MatrixAddition9(numTests, problemSize);
// MatrixAddition10(numTests, problemSize);
// MatrixAddition11(numTests, problemSize);
// MatrixAddition12(numTests, problemSize);
// MatrixAddition13(numTests, problemSize);
MatrixAddition14(numTests, problemSize);
    cout << "\n";
}
