#include <Profile/IntWrapper.h>
#include <boost/lexical_cast.hpp>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <LibUtilities/LinearAlgebra/Blas.hpp>
#include <Profile/VariableCostObject.h>
#include <boost/timer.hpp>
#include <boost/progress.hpp>
#include <iostream>
#include <sys/time.h>
#include <tr1/memory>
#include <loki/SmartPtr.h>

using std::cout;
using std::endl;



void VariableCostOperationAddition2(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1;
	}
	cout << "VariableCostOperationAddition2 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition2 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition3(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2;
	}
	cout << "VariableCostOperationAddition3 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition3 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition4(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3;
	}
	cout << "VariableCostOperationAddition4 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition4 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition5(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4;
	}
	cout << "VariableCostOperationAddition5 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition5 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition6(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5;
	}
	cout << "VariableCostOperationAddition6 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition6 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition7(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6;
	}
	cout << "VariableCostOperationAddition7 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition7 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition8(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7;
	}
	cout << "VariableCostOperationAddition8 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition8 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition9(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	VariableCostObject m8;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8;
	}
	cout << "VariableCostOperationAddition9 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition9 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition10(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	VariableCostObject m8;
	VariableCostObject m9;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9;
	}
	cout << "VariableCostOperationAddition10 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition10 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition11(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	VariableCostObject m8;
	VariableCostObject m9;
	VariableCostObject m10;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10;
	}
	cout << "VariableCostOperationAddition11 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition11 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition12(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	VariableCostObject m8;
	VariableCostObject m9;
	VariableCostObject m10;
	VariableCostObject m11;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11;
	}
	cout << "VariableCostOperationAddition12 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition12 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition13(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	VariableCostObject m8;
	VariableCostObject m9;
	VariableCostObject m10;
	VariableCostObject m11;
	VariableCostObject m12;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11 + m12;
	}
	cout << "VariableCostOperationAddition13 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition13 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void VariableCostOperationAddition14(unsigned int numIterations, unsigned int matrixSize)
{
	VariableCostObject m0;
	VariableCostObject m1;
	VariableCostObject m2;
	VariableCostObject m3;
	VariableCostObject m4;
	VariableCostObject m5;
	VariableCostObject m6;
	VariableCostObject m7;
	VariableCostObject m8;
	VariableCostObject m9;
	VariableCostObject m10;
	VariableCostObject m11;
	VariableCostObject m12;
	VariableCostObject m13;
	boost::timer t;

	for(unsigned int i = 0; i < numIterations; ++i)
	{
		VariableCostObject result = m0 + m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9 + m10 + m11 + m12 + m13;
	}
	cout << "VariableCostOperationAddition14 Total Time: " << t.elapsed() << "\n";
	cout << "VariableCostOperationAddition14 Per Op: " << t.elapsed()/(double)numIterations << "\n";
}

void TensProdBwdTrans(const Nektar::ConstArray<Nektar::OneD, Nektar::NekDouble>& inarray, 
            Nektar::Array<Nektar::OneD, Nektar::NekDouble> &outarray,
            boost::shared_ptr<Nektar::NekMatrix<double> >& m)
{
    int nq = m->GetRows();
    int m_ncoeffs = m->GetColumns();
    
    Blas::Dgemv('N',nq,m_ncoeffs,1.0,m->GetPtr().get(),
               nq, &inarray[0], 1.0, 0.0, &outarray[0], 1.0);
    
//    Nektar::NekVector<double> v_in_0(m_ncoeffs, inarray);
//    Nektar::NekVector<double> v_out_0(m_ncoeffs, outarray, Nektar::eWrapper);
//    v_out_0 = (*m)*v_in_0;
//    
//    Nektar::NekVector<const double> v_in_1(m_ncoeffs, inarray, Nektar::eWrapper);
//    Nektar::NekVector<double> v_out_1(m_ncoeffs, outarray, Nektar::eWrapper);
//    v_out_1 = (*m)*v_in_1;
}

template<typename T>
class NekPtr
{
    public:
        NekPtr(T* value) : m_data(value), m_count(new ref_count()) 
        {
        }
        
        NekPtr(const NekPtr<T>& rhs) :
            m_data(rhs.m_data),
            m_count(rhs.m_count)
        {
            ++(*m_count).count;
        }
    
        NekPtr()
        {
            --(*m_count).count;
            if( m_count->count <= 0 )
            {
                delete m_count;
                m_count = 0;
                delete m_data;
            }
        }

        T& operator*() { return *m_data; }

    private:
        NekPtr<T>& operator=(const NekPtr<T>& rhs)
        {
            m_data = rhs.m_data;
            m_count = rhs.m_count;
            ++(*m_count).count;
        }

        struct ref_count
        {
            ref_count() : count(1) {}
            ref_count(const ref_count& r) : count(r.count) {}
            ref_count& operator=(const ref_count& r) { count = r.count; return *this; }

            int count;
        };
        T* m_data;
        ref_count* m_count;
};

class Foo
{
    public:
        Foo() :
            m_rawPointer(new int),
            m_boost(new int),
            m_tr1(new int),
            m_loki(new int),
            m_nekptr(new int)
        {
            *m_rawPointer = 1;
            *m_boost = 1;
            *m_tr1 = 1;
            *m_loki = 1;
            *m_nekptr = 1;
        }
        
        int* GetRawPointer() { return m_rawPointer; }
        boost::shared_ptr<int> GetBoost() { return m_boost; }
        std::tr1::shared_ptr<int> GetTR1() { return m_tr1; }
        Loki::SmartPtr<int> GetLoki() { return m_loki; }
        NekPtr<int> GetNekPatr() { return m_nekptr; }
        
    private:
        int* m_rawPointer;
        boost::shared_ptr<int> m_boost;
        std::tr1::shared_ptr<int> m_tr1;
        Loki::SmartPtr<int> m_loki;
        NekPtr<int> m_nekptr;
        
};

int main(int argc, char** argv)
{
//	if( argc != 3 ) { cout << "Usage: Profile <numTests> <problemSize>\n"; return 1; }
//    cout.precision(20);
//	unsigned int numTests = boost::lexical_cast<unsigned int>(argv[1]);
//	unsigned int problemSize = boost::lexical_cast<unsigned int>(argv[2]);
//	VariableCostObject::size = problemSize;
//
//    boost::shared_ptr<Nektar::NekMatrix<double> > m(new Nektar::NekMatrix<double>(problemSize, problemSize));
//    Nektar::ConstArray<Nektar::OneD, Nektar::NekDouble> inarray(problemSize, 1.0);
//    Nektar::Array<Nektar::OneD, Nektar::NekDouble> outarray(problemSize, 1.0);
//    
//    unsigned int numIterations = numTests;
//    boost::timer t;
//    for(unsigned int i = 0; i < numIterations; ++i)
//    {
//        TensProdBwdTrans(inarray, outarray, m);
//    }
//    cout << "VariableCostOperationAddition14 Total Time: " << t.elapsed() << "\n";
//	cout << "VariableCostOperationAddition14 Per Op: " << t.elapsed()/(double)numIterations << "\n";
//
//
//    
////VariableCostOperationAddition2(numTests, problemSize);
////VariableCostOperationAddition3(numTests, problemSize);
////VariableCostOperationAddition4(numTests, problemSize);
////VariableCostOperationAddition5(numTests, problemSize);
////VariableCostOperationAddition6(numTests, problemSize);
////VariableCostOperationAddition7(numTests, problemSize);
////VariableCostOperationAddition8(numTests, problemSize);
////VariableCostOperationAddition9(numTests, problemSize);
////VariableCostOperationAddition10(numTests, problemSize);
////VariableCostOperationAddition11(numTests, problemSize);
////VariableCostOperationAddition12(numTests, problemSize);
////VariableCostOperationAddition13(numTests, problemSize);
////VariableCostOperationAddition14(numTests, problemSize);

    unsigned int numTrials = boost::lexical_cast<unsigned int>(argv[1]);
    
    timeval timer1, timer2;
    double time1, time2;
    double exeTime;
    double rawTime = 0.0;
    Foo obj;
    
    gettimeofday(&timer1, NULL);
    int result = 0;
    for(unsigned int i = 0; i < numTrials; ++i)
    {
        result += *(obj.GetRawPointer());
    }
    
    gettimeofday(&timer2, NULL);
    time1 = timer1.tv_sec*1000000.0+(timer1.tv_usec);
    time2 = timer2.tv_sec*1000000.0+(timer2.tv_usec);
    exeTime = (time2-time1);
    rawTime = exeTime;
    std::cout << "Raw Pointer Time: " << exeTime << std::endl;
    std::cout << result << std::endl;
    
    gettimeofday(&timer1, NULL);
    result = 0;
    for(unsigned int i = 0; i < numTrials; ++i)
    {
        result += *(obj.GetBoost());
    }
    
    gettimeofday(&timer2, NULL);
    time1 = timer1.tv_sec*1000000.0+(timer1.tv_usec);
    time2 = timer2.tv_sec*1000000.0+(timer2.tv_usec);
    exeTime = (time2-time1);
    std::cout << "Boost Pointer Time: " << exeTime << std::endl;
    std::cout << exeTime/rawTime << std::endl;
    std::cout << result << std::endl;
    
    
    
    gettimeofday(&timer1, NULL);
    result = 0;
    for(unsigned int i = 0; i < numTrials; ++i)
    {
        result += *(obj.GetTR1());
    }
    
    gettimeofday(&timer2, NULL);
    time1 = timer1.tv_sec*1000000.0+(timer1.tv_usec);
    time2 = timer2.tv_sec*1000000.0+(timer2.tv_usec);
    exeTime = (time2-time1);
    std::cout << "TR1 Pointer Time: " << exeTime << std::endl;
    std::cout << exeTime/rawTime << std::endl;
    std::cout << result << std::endl;
    
    
    gettimeofday(&timer1, NULL);
    result = 0;
    for(unsigned int i = 0; i < numTrials; ++i)
    {
        result += *(obj.GetLoki());
    }
    
    gettimeofday(&timer2, NULL);
    time1 = timer1.tv_sec*1000000.0+(timer1.tv_usec);
    time2 = timer2.tv_sec*1000000.0+(timer2.tv_usec);
    exeTime = (time2-time1);
    std::cout << "Loki Pointer Time: " << exeTime << std::endl;
    std::cout << exeTime/rawTime << std::endl;
    std::cout << result << std::endl;


    gettimeofday(&timer1, NULL);
    result = 0;
    for(unsigned int i = 0; i < numTrials; ++i)
    {
        result += *(obj.GetNekPatr());
    }
    
    gettimeofday(&timer2, NULL);
    time1 = timer1.tv_sec*1000000.0+(timer1.tv_usec);
    time2 = timer2.tv_sec*1000000.0+(timer2.tv_usec);
    exeTime = (time2-time1);
    std::cout << "NekPtr Pointer Time: " << exeTime << std::endl;
    std::cout << exeTime/rawTime << std::endl;
    std::cout << result << std::endl;

      
}
