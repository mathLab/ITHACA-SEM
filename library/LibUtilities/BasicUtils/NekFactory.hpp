///////////////////////////////////////////////////////////////////////////////
//
// File: NekFactory.hpp
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
// Description: Factory pattern class for Nektar
//
///////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_FACTORY_HPP
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_NEK_FACTORY_HPP

#include <boost/shared_ptr.hpp>

#include <iostream>
#include <map>

namespace Nektar
{
    namespace LibUtilities
    {

        /**
         * @class NekFactory
         * \brief Provides a generic Factory class.
         *
         * Implements a generic class factory. To allow a class to be
         * instantiated by the factory, the following are required in each
         * class definition:
         * \code
         *   static [baseclass]* create([paramtype1] &P) {
         *      return new [derivedclass](P);
         *   }
         *   static std::string className;
         * \endcode
         * and outside the class definition in the implementation:
         * \code
         *   std::string [derivedclass]::className
         *      = Factory<std::string,[baseclass],[paramtype1]>::
         *          RegisterCreatorFunction("[derivedclass]",
         *                              [derivedclass]::create,"Description");
         * \endcode
         * The assignment of the static variable className is done through the
         * call to RegisterCreatorFunction, which registers the class with the
         * factory prior to the start of the main() routine.
         *
         * To create an instance of a derived class, for instance:
         * \code
         *   [baseclass]* var_name =
         *      Factory<std::string,[baseclass],[paramtype1]>
         *              ::CreateInstance("[derivedclass]",Param1);
         * \endcode
         */

        template < typename tKey,	       	// reference tag (e.g. string, int)
                   typename tBase,	       	// base class
                   typename tParam1,        // First parameter type
                   typename tDescription = std::string,
                   typename tPredicator = std::less<tKey> >
        class NekFactory
        {
        public:
            /// Shared pointer to an object of baseclass type.
            typedef boost::shared_ptr<tBase> tBaseSharedPtr;
            /// CreatorFunction type which takes parameter and returns base
            /// class shared pointer.
            typedef tBaseSharedPtr (*CreatorFunction) (tParam1);

            /// Define a struct to hold the information about a module.
            struct ModuleEntry
            {
                ModuleEntry(CreatorFunction pFunc, const tDescription pDesc)
                    : m_func(pFunc),
                      m_desc(pDesc)
                {
                }

                /// Function used to create instance of class.
                CreatorFunction m_func;
                /// Description of class for use in listing available classes.
                tDescription m_desc;
            };

            /// Factory map between key and module data.
            typedef std::map<tKey, ModuleEntry, tPredicator> TMapFactory;
            /// Iterator for factory map
            typedef typename TMapFactory::iterator TMapFactoryIterator;


            /**
             * @brief Create an instance of the class referred to by \c idKey.
             *
             * Searches the factory's map for the given key and returns a shared
             * base class pointer to a new instance of the associated class.
             * @param   idKey           Key of class to create.
             * @param   x               Parameter to pass to class constructor.
             * @returns                 Base class pointer to new instance.
             */
            static tBaseSharedPtr CreateInstance(tKey idKey, tParam1 x)
            {
                // Check we've not been given a blank key. If so, display an
                // error and list the available modules.
                if (idKey == "")
                {
                    std::cout << "ERROR no module name given." << std::endl;
                    PrintAvailableClasses();
                    abort();
                }

                // Now try and find the key in the map.
                TMapFactoryIterator it = getMapFactory()->find(idKey);

                // If successful, check the CreatorFunction is defined and
                // create a new instance of the class.
                if (it != getMapFactory()->end())
                {
                    if (it->second.m_func)
                    {
                        try
                        {
                            return it->second.m_func(x);
                        }
                        catch (const std::string& s)
                        {
                            std::cout << "ERROR Creating module: " << s << std::endl;
                            abort();
                        }
                    }
                }

                // If we get this far, the key does not exist, so display an
                // error and print the list of available modules.
                std::cout << "No such module: " << idKey << std::endl;
                PrintAvailableClasses();
                abort();
            }


            /**
             * @brief Register a class with the factory.
             *
             * This function is called by each class in a static context (prior
             * to the execution of main()) and creates an entry for the class
             * in the factory's map.
             * @param   idKey           Key used to reference the class.
             * @param   classCreator    Function to call to create an instance
             *                          of this class.
             * @param   pDesc           Optional description of class.
             * @returns                 The given key \c idKey.
             */
            static tKey RegisterCreatorFunction(tKey idKey,
                                                CreatorFunction classCreator,
                                                tDescription pDesc = "") {
                ModuleEntry e(classCreator, pDesc);
                getMapFactory()->insert(std::pair<tKey,ModuleEntry>(idKey, e));
                return idKey;
            }


            /**
             * @brief Prints the available classes to stdout.
             */
            static void PrintAvailableClasses()
            {
                std::cout << std::endl << "Available classes: " << std::endl;
                TMapFactoryIterator it;
                for (it = getMapFactory()->begin(); it != getMapFactory()->end(); ++it)
                {
                    std::cout << "  " << it->first;
                    if (it->second.m_desc != "")
                    {
                        std::cout << ":" << std::endl << "    "
                                  << it->second.m_desc << std::endl;
                    }
                    else
                    {
                        std::cout << std::endl;
                    }
                }
            }

        protected:
            /**
             * @brief Ensure the factory's map is created.
             * @returns                 The factory's map.
             */
            static TMapFactory * getMapFactory() {
                static TMapFactory mMapFactory;
                return &mMapFactory;
            }

        private:
            /// Constructor - we never explicitly instantiate this class
            NekFactory() {}
            /// Destructor
            ~NekFactory() {}

        };

    }
}

#endif
