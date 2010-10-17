namespace Nektar
{
namespace LibUtilities
{
/**
 * @namespace Nektar::LibUtilities
 * @brief The namespace associated with the the LibUtilities library
 * (@ref pageLibUtilities "LibUtilities introduction")
 *
 * @page pageLibUtilities LibUtilities sublibrary
 * This contains the underlying building blocks for constructing a Spectral
 * Element formulation including linear algebra, polynomial routines and memory
 * management.
 *
 * @section sectionLibUtilitiesContents Contents
 * - Basic Constants
 * - Basic Utilities (@subpage pageNekArrays)
 * - Expression Templates
 * - Foundations
 * - Interpreter
 * - Kernel
 * - Linear Algebra
 * - Memory Management
 * - Nodal Data
 * - Polynomial Subroutines
 * - Time Integration
 */

/**
 * @page pageNekArrays Nektar++ Arrays
 * An Array is a thin wrapper around native arrays.  Arrays provide
 * all the functionality of native arrays, with the additional benefits of
 * automatic use of the Nektar++ memory pool, automatic memory allocation and
 * deallocation, bounds checking in debug mode, and easier to use multi-dimensional
 * arrays.
 *
 * Arrays are templated to allow compile-time customization of its dimensionality
 * and data type.
 *
 * @param Dim Must be a type with a static unsigned integer called Value that specifies the
 *        array's dimensionality.  For example
 *        @code
 *        struct TenD
 *        {
 *            static unsigned int Value = 10;
 *        };
 *        @endcode
 *
 * @param DataType The type of data to store in the array.
 *
 * It is often useful to create a class member Array that is shared with
 * users of the object without letting the users modify the array.  To allow this behavior,
 * Array<Dim, DataType> inherits from Array<Dim, const DataType>.  The following example
 * shows what is possible using this approach:
 *
 * @code
 * class Sample
 * {
 *     public:
 *         Array<OneD, const double>& getData() const { return m_data; }
 *         void getData(Array<OneD, const double>& out) const { out = m_data; }
 *
 *     private:
 *         Array<OneD, double> m_data;
 * };
 * @endcode
 *
 * In this example, each instance of Sample contains an array.  The getData method
 * gives the user access to the array values, but does not allow modification of those
 * values.
 *
 * @section efficiency Efficiency Considerations
 * Tracking memory so it is deallocated only when no more Arrays reference it does
 * introduce overhead when copying and assigning Arrays.  In most cases this loss
 * of efficiency is not noticeable.  There are some cases, however, where it can cause
 * a significant performance penalty (such as in tight inner loops).  If needed, Arrays
 * allow access to the C-style array through the Array::data member function.
 *
 * Because of this, it is sometimes advisable to return Arrays by reference:
 *
 * @code
 * // Instead of this:
 * Array<OneD, const NekDouble> lessEfficient(points->GetZ());
 *
 * // Do this
 * const Array<OneD, const NekDouble>& moreEfficient = points->GetZ();
 * @endcode
 */
}
}
