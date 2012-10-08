
#include <LibUtilities/Memory/ThreadSpecificPool.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>

#include <boost/shared_ptr.hpp>

#include <boost/mpl/size_t.hpp>



#include <iostream>

using namespace std;
using namespace Nektar;

class Disabled
{
    public:
        Disabled()
        {
            cout << "Creating a disabled." << endl;
        }

        ~Disabled()
        {
            cout << "Destroying a disabled." << endl;
        }
};

class Enabled
{
    public:
        Enabled()
        {
            cout << "Creating a Enabled." << endl;
        }

        ~Enabled()
        {
            cout << "Destroying a Enabled." << endl;
        }
};



int main()
{
//     // Object allocation.
//     Disabled* t1 = MemoryManager::Allocate<Disabled>();
//     MemoryManager::Deallocate(t1);
//     assert(t1 == NULL);
// 
//     const Disabled* t2 = MemoryManager::Allocate<Disabled>();
//     MemoryManager::Deallocate(t2);
//     assert(t2 == NULL);
// 
//     Enabled* t3 = MemoryManager::Allocate<Enabled>();
//     MemoryManager::Deallocate(t3);
//     assert(t3 == NULL);
// 
//     boost::shared_ptr<Disabled> t4 = MemoryManager::AllocateSharedPtr<Disabled>();
//     boost::shared_ptr<Enabled> t5 = MemoryManager::AllocateSharedPtr<Enabled>();
// 
//     // This doesn't compile, it really doesn't make sense to allow it to.
//     //double* d = MemoryManager::Allocate<double>();
// 
//     cout << "\nDouble Array." << endl;
//     double* d = MemoryManager::AllocateArray<10, double>();
//     MemoryManager::DeallocateArray<10>(d);
// 
//     cout << "\nDisabled Array" << endl;
//     Disabled* a1 = MemoryManager::AllocateArray<2, Disabled>();
//     MemoryManager::DeallocateArray<2>(a1);
// 
//     cout << "\nEnabled Array." << endl;
//     Enabled* a2 = MemoryManager::AllocateArray<3, Enabled>();
//     MemoryManager::DeallocateArray<3>(a2);
// 
//     boost::shared_array<Disabled> a3 = MemoryManager::AllocateSharedArray<3, Disabled>();
//     boost::shared_array<Enabled> a4 = MemoryManager::AllocateSharedArray<4, Enabled>();
//     boost::shared_array<double> a5 = MemoryManager::AllocateSharedArray<5, double>();
// 
//     boost::shared_array<double> a6 = MemoryManager::AllocateSharedArray<double>(12);
// 
//     cout << "\nDestroy all shared pointers." << endl;


//     boost::shared_array<Class2> a = MemoryManager::AllocateArray<Class2>();
//     boost::shared_array<double> d = MemoryManager::AllocateArray<double>();
    // Allocate a raw
//     boost::thread_group g;
//     for(unsigned int i = 0; i < 2; ++i)
//     {
//         //g.create_thread(testThread);
//         g.create_thread(testNakedThread);
//     }
//
//     g.join_all();



    return 0;
}

