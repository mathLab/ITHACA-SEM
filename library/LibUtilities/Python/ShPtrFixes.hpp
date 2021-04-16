///////////////////////////////////////////////////////////////////////////////
//
// File: ShPtrFixes.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Original code from Boost (c) David Abrahams 2002, Stefan Seefeld 2016.
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
// Description: Patches for boost.python to support std::shared_ptr for boost
//              versions < 1.63.
//
///////////////////////////////////////////////////////////////////////////////

/*
 * This file is an unpleasant series of includes from boost.python which is
 * designed to support std::shared_ptr on versions of boost lower than 1.63
 * (where support was first introduced). This file overrides or augments
 * various headers from the release to enable this.
 *
 * Although one can mostly handle std::shared_ptr with the use of the
 * py::implicitly_convertible, in some cases the allocation/deallocation process
 * seems to lead to memory issues. This file therefore takes the patch made in
 * boost python changeset 97e4b34a159 which enables this. This is a pretty ugly
 * approach, but various other strategies did not seem to work!
 */

//
// From <boost/python/detail/is_shared_ptr.hpp>
//

#define IS_SHARED_PTR_DWA2003224_HPP

#include <boost/python/detail/is_xxx.hpp>
#include <boost/shared_ptr.hpp>

namespace boost { namespace python { namespace detail {

BOOST_PYTHON_IS_XXX_DEF(shared_ptr, shared_ptr, 1)
template <typename T>
struct is_shared_ptr<std::shared_ptr<T> > : std::true_type {};

}}} // namespace boost::python::detail

//
// From <boost/python/converter/registered.hpp>
//

#include <boost/python/type_id.hpp>
#include <boost/python/converter/registry.hpp>
#include <boost/python/converter/registrations.hpp>
#include <boost/type_traits/transform_traits.hpp>
#include <boost/type_traits/cv_traits.hpp>
#include <boost/type_traits/is_void.hpp>
#include <boost/detail/workaround.hpp>
#include <boost/type.hpp>
#include <memory>

namespace boost { namespace python { namespace converter { namespace detail {

template <class T>
inline void
register_shared_ptr0(std::shared_ptr<T>*)
{
    registry::lookup_shared_ptr(type_id<std::shared_ptr<T> >());
}

}}}}

// From <boost/python/detail/value_is_shared_ptr.hpp>

#define VALUE_IS_SHARED_PTR_DWA2003224_HPP
#include <boost/python/detail/value_is_xxx.hpp>
#include <boost/python/detail/is_shared_ptr.hpp>

namespace boost { namespace python { namespace detail {

template <class X_>
struct value_is_shared_ptr
{
  static bool const value = is_shared_ptr<typename remove_cv<
                                            typename remove_reference<X_>
                                              ::type>
                                            ::type>
    ::value;
  typedef mpl::bool_<value> type;
};

}}} // namespace boost::python::detail

//
// From <boost/python/converter/shared_ptr_from_python.hpp>
//

#define SHARED_PTR_FROM_PYTHON_DWA20021130_HPP

#include <boost/python/handle.hpp>
#include <boost/python/converter/shared_ptr_deleter.hpp>
#include <boost/python/converter/from_python.hpp>
#include <boost/python/converter/rvalue_from_python_data.hpp>
#include <boost/python/converter/registered.hpp>
#ifndef BOOST_PYTHON_NO_PY_SIGNATURES
# include <boost/python/converter/pytype_function.hpp>
#endif
#include <boost/shared_ptr.hpp>

namespace boost { namespace python { namespace converter {

template <class T, template <typename> class SP>
struct shared_ptr_from_python_new
{
  shared_ptr_from_python_new()
  {
    converter::registry::insert(&convertible, &construct, type_id<SP<T> >()
#ifndef BOOST_PYTHON_NO_PY_SIGNATURES
                                , &converter::expected_from_python_type_direct<T>::get_pytype
#endif
                                );
  }

 private:
  static void* convertible(PyObject* p)
  {
    if (p == Py_None)
      return p;

    return converter::get_lvalue_from_python(p, registered<T>::converters);
  }

  static void construct(PyObject* source, rvalue_from_python_stage1_data* data)
  {
    void* const storage = ((converter::rvalue_from_python_storage<SP<T> >*)data)->storage.bytes;
    // Deal with the "None" case.
    if (data->convertible == source)
      new (storage) SP<T>();
    else
    {
      SP<void> hold_convertible_ref_count(
         (void*)0, shared_ptr_deleter(handle<>(borrowed(source))) );
      // use aliasing constructor
      new (storage) SP<T>(hold_convertible_ref_count,
                          static_cast<T*>(data->convertible));
    }

    data->convertible = storage;
  }
};

// DM: Define original shared_ptr_from_python struct so that we can
// avoid changes to class_metadata.hpp.
template <class T>
struct shared_ptr_from_python
{
  shared_ptr_from_python()
  {
    shared_ptr_from_python_new<T, boost::shared_ptr>();
    shared_ptr_from_python_new<T, std::shared_ptr>();
  }
};



}}} // namespace boost::python::converter


//
// From <boost/python/converter/shared_ptr_to_python.hpp>
//

#include <boost/python/refcount.hpp>
#include <boost/python/converter/shared_ptr_deleter.hpp>
#include <boost/python/detail/none.hpp>
#include <boost/get_pointer.hpp>

namespace boost { namespace python { namespace converter {

template <class T>
PyObject* shared_ptr_to_python(std::shared_ptr<T> const& x)
{
  if (!x)
    return python::detail::none();
  else if (shared_ptr_deleter* d = std::get_deleter<shared_ptr_deleter>(x))
    return incref(get_pointer(d->owner));
  else
    return converter::registered<std::shared_ptr<T> const&>::converters.to_python(&x);
}

}}} // namespace boost::python::converter

//
// From <boost/python/to_python_value.hpp>: this is a pretty big hack, but
// needs to be fully included due to the definition of new private member
// variables.
//

#define TO_PYTHON_VALUE_DWA200221_HPP
#include <boost/python/detail/prefix.hpp>

#include <boost/python/refcount.hpp>
#include <boost/python/tag.hpp>
#include <boost/python/handle.hpp>

#include <boost/python/converter/registry.hpp>
#include <boost/python/converter/registered.hpp>
#include <boost/python/converter/builtin_converters.hpp>
#include <boost/python/converter/object_manager.hpp>
#include <boost/python/converter/shared_ptr_to_python.hpp>

#include <boost/python/detail/value_is_shared_ptr.hpp>
#include <boost/python/detail/value_arg.hpp>

#include <boost/type_traits/transform_traits.hpp>

#include <boost/mpl/if.hpp>
#include <boost/mpl/or.hpp>
#include <boost/type_traits/is_const.hpp>

namespace boost { namespace python { 

namespace detail
{
#ifndef BOOST_PYTHON_NO_PY_SIGNATURES

template <bool is_const_ref>
struct object_manager_get_pytype
{
   template <class U>
   static PyTypeObject const* get( U& (*)() =0)
   {
      return converter::object_manager_traits<U>::get_pytype();
   }
};

template <>
struct object_manager_get_pytype<true>
{
   template <class U>
   static PyTypeObject const* get( U const& (*)() =0)
   {
      return converter::object_manager_traits<U>::get_pytype();
   }
};

#endif

  template <class T>
  struct object_manager_to_python_value
  {
      typedef typename value_arg<T>::type argument_type;
    
      PyObject* operator()(argument_type) const;
#ifndef BOOST_PYTHON_NO_PY_SIGNATURES
      typedef boost::mpl::bool_<is_handle<T>::value> is_t_handle;
      typedef boost::detail::indirect_traits::is_reference_to_const<T> is_t_const;
      PyTypeObject const* get_pytype() const {
          return get_pytype_aux((is_t_handle*)0);
      }

      inline static PyTypeObject const* get_pytype_aux(mpl::true_*) {return converter::object_manager_traits<T>::get_pytype();}
      
      inline static PyTypeObject const* get_pytype_aux(mpl::false_* ) 
      {
          return object_manager_get_pytype<is_t_const::value>::get((T(*)())0);
      }
      
#endif 

      // This information helps make_getter() decide whether to try to
      // return an internal reference or not. I don't like it much,
      // but it will have to serve for now.
      BOOST_STATIC_CONSTANT(bool, uses_registry = false);
  };

  
  template <class T>
  struct registry_to_python_value
  {
      typedef typename value_arg<T>::type argument_type;
    
      PyObject* operator()(argument_type) const;
#ifndef BOOST_PYTHON_NO_PY_SIGNATURES
      PyTypeObject const* get_pytype() const {return converter::registered<T>::converters.to_python_target_type();}
#endif

      // This information helps make_getter() decide whether to try to
      // return an internal reference or not. I don't like it much,
      // but it will have to serve for now.
      BOOST_STATIC_CONSTANT(bool, uses_registry = true);
  };

  template <class T>
  struct shared_ptr_to_python_value
  {
      typedef typename value_arg<T>::type argument_type;
    
      PyObject* operator()(argument_type) const;
#ifndef BOOST_PYTHON_NO_PY_SIGNATURES
      PyTypeObject const* get_pytype() const {return get_pytype((boost::type<argument_type>*)0);}
#endif 
      // This information helps make_getter() decide whether to try to
      // return an internal reference or not. I don't like it much,
      // but it will have to serve for now.
      BOOST_STATIC_CONSTANT(bool, uses_registry = false);
  private:
#ifndef BOOST_PYTHON_NO_PY_SIGNATURES
    template <class U>
    PyTypeObject const* get_pytype(boost::type<shared_ptr<U> &> *) const {return converter::registered<U>::converters.to_python_target_type();}
    template <class U>
    PyTypeObject const* get_pytype(boost::type<const shared_ptr<U> &> *) const {return converter::registered<U>::converters.to_python_target_type();}
# if __cplusplus >= 201103L
    template <class U>
    PyTypeObject const* get_pytype(boost::type<std::shared_ptr<U> &> *) const {return converter::registered<U>::converters.to_python_target_type();}
    template <class U>
    PyTypeObject const* get_pytype(boost::type<const std::shared_ptr<U> &> *) const {return converter::registered<U>::converters.to_python_target_type();}
# endif
#endif
  };
}

template <class T>
struct to_python_value
    : mpl::if_<
          detail::value_is_shared_ptr<T>
        , detail::shared_ptr_to_python_value<T>
        , typename mpl::if_<
              mpl::or_<
                  converter::is_object_manager<T>
                , converter::is_reference_to_object_manager<T>
              >
            , detail::object_manager_to_python_value<T>
            , detail::registry_to_python_value<T>
          >::type
      >::type
{
};

//
// implementation 
//
namespace detail
{
  template <class T>
  inline PyObject* registry_to_python_value<T>::operator()(argument_type x) const
  {
      return converter::registered<argument_type>::converters.to_python(&x);
  }

  template <class T>
  inline PyObject* object_manager_to_python_value<T>::operator()(argument_type x) const
  {
      return python::upcast<PyObject>(
          python::xincref(
              get_managed_object(x, tag))
          );
  }

  template <class T>
  inline PyObject* shared_ptr_to_python_value<T>::operator()(argument_type x) const
  {
      return converter::shared_ptr_to_python(x);
  }
}

}} // namespace boost::python


