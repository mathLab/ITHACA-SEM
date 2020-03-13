#include <boost/multi_array.hpp>
#include <cassert>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <iostream>

typedef boost::multi_array<double, 3> ThreeDArray;

using boost::shared_ptr;

void PassLocalCopy(shared_ptr<const ThreeDArray> &ar)
{
    static shared_ptr<ThreeDArray> localArray = Nektar::MemoryManager::AllocateSharedPtr<ThreeDArray>(boost::extents[4][5][6]);

    (*localArray)[3][4][5] = 10.0;

    ar = localArray;
}

shared_ptr<const ThreeDArray> ReturnLocalCopy(void)
{
    static shared_ptr<ThreeDArray> localArray = Nektar::MemoryManager::AllocateSharedPtr<ThreeDArray>(boost::extents[4][5][6]);

    (*localArray)[3][4][5] = 10.0;

    return localArray;
}

int main(void)
{
    // The three different types of access to multiarrays.
    typedef boost::multi_array<double, 1> OneDArrayT;
    typedef boost::multi_array_ref<double, 1> OneDArrayRefT;
    typedef boost::const_multi_array_ref<double, 1> OneDConstArrayRefT;


    // Allocate 50 element array with indices ranging from 0 to 49 (default).
    OneDArrayT my1dArray(boost::extents[50]);
    my1dArray[20] = 5;

    // Share the data from the above array with another array of the smaller extent.
    // object.data() retrieves the multi_array's underlying data.  Wrap the last
    // 30 elements in the reference.
    OneDArrayRefT my1dArrayRef(my1dArray.data()+20, boost::extents[30]);

    my1dArrayRef[0] = 7;

    // Create a view of the matrix from the range in the original matrix of 20 to 39.
    // The view will still use 0 to 19 as indices.  However, you can't use *.data() to
    // retrieve the data of the view.

    // Allows a range to be specified.  This is used to provide a view of my1dArray,
    // below.  A view is a way to look at a part of the array, or to remap the array
    // into different shape.
    typedef OneDArrayT::index_range range;

    OneDArrayT::array_view<1>::type my1Dview = my1dArray[boost::indices[range(20,40)]];

    // This refers to the same element used above when it was assigned to 5 then changed
    // to 7 in the reference.
    std::cout << my1Dview[0] << std::endl;

    // Shape just tells me the extent of the array in each dimension.  If we did not
    // use a base of 0 then we would have to take that into account.  We can use
    // iterators here that simplifies this.
    const OneDArrayT::size_type *sizes = my1Dview.shape();
    OneDArrayT::index index;
    for (index = 0; index<sizes[0]; ++index)
    {
        my1Dview[index] = index;
    }

    ThreeDArray my3dArray(boost::extents[3][4][5]);

    // If an array is multidimensional, then a subarray can be generated that reduces
    // the dimension to anything less than the original array dimension and greater
    // than or equal to one.  In this example we are reducing the three-dimensional
    // array down to two dimensions and we are specifically looking at the subarray
    // existing at index 1 of the first array's most significant (not sure that is
    // the correct term) dimension.  In our case it is the one that has extent 3.
    ThreeDArray::subarray<2>::type my2DSubarray = my3dArray[1];

    // To pass such an array we need to wrap it in a shared_ptr, then all
    // the usual syntax applies including constifying.
    boost::shared_ptr<const ThreeDArray> ar;
    PassLocalCopy(ar);

    std::cout << (*ar)[3][4][5] << std::endl;

    // Return a shared_ptr to the one allocated internally.  Pass back
    // as the return value.
    boost::shared_ptr<const ThreeDArray> ar2;
    ar2 = ReturnLocalCopy();

    std::cout << (*ar)[3][4][5] << std::endl;

    return 0;
}
