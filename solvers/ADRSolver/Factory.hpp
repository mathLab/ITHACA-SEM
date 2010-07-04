/*
 * $Id: Factory.h 128 2008-07-06 18:30:47Z cc $
 *
 * Copyright (c) 2005 - $Date: 2008-07-06 19:30:47 +0100 (Sun, 06 Jul 2008) $ Chris Cantwell and Paul Clifford
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef CLASS_FACTORY
#define CLASS_FACTORY

#include <boost/shared_ptr.hpp>

#include <iostream>
#include <map>

/**
 * @class Factory
 * \brief Provides a generic Factory class.
 *
 * Implements a generic class factory. Classes should register themselves with
 * the class factory immediately, so they are available for use in the main()
 * routine. To allow a class to be instantiated by the factory, the following
 * are required in each class definition:
 * \code
 *   static [baseclass]* create([paramtype1] &P, [paramtype2] &Q) {
 *      return new [derivedclass](P,Q);
 *   }
 *   static std::string className;
 * \endcode
 * and outside the class in the implementation file:
 * \code
 *   std::string [derivedclass]::className
 *      = Factory<[paramtype1],[paramtype2],std::string,[baseclass]>::
 *          RegisterCreatorFunction("[derivedclass]",[derivedclass]::create);
 * \endcode
 * So, we have a static create function to actually create an instance of this
 * class, and the static init dummy variable is used to force the class to
 * register immediately on execution.
 *
 * To create an instance of a derived class, for instance:
 * \code
 *   [baseclass]* var_name =
 *      Factory<[paramtype1],[paramtype2],std::string,[baseclass]>
 *              ::CreateInstance("[derivedclass]",Param1,Param2)
 * \endcode
 */
template < typename tKey,				// reference tag (e.g. string, int)
		   typename tBase,				// base class
		   typename tParam1,            // First parameter type
		   typename tParam2,
		   typename tPredicator = std::less<tKey> >
class Factory {
	public:
		/// Constructor
		Factory() {}
		/// Destructor
		~Factory() {}

		typedef boost::shared_ptr<tBase> tBaseSharedPtr;
		/// typedef for the creator function, return of type base class.
    	typedef tBaseSharedPtr (*CreatorFunction) (tParam1, tParam2);
		/// Shorthand for our map
	    typedef std::map<tKey, CreatorFunction, tPredicator> TMapFactory;
		/// typename for our maps iterator
		typedef typename TMapFactory::iterator TMapFactoryIterator;

		/// Create an instance of the class referred to by \c idKey
	    static tBaseSharedPtr CreateInstance(tKey idKey, tParam1 x, tParam2 y) {
			TMapFactoryIterator it = getMapFactory()->find(idKey);
	        if (it != getMapFactory()->end()) {
	            if (it->second) {
	                try
	                {
	                    return it->second(x, y);
	                }
	                catch (const std::string& s)
	                {
	                    std::cout << "ERROR Creating module: " << s << std::endl;
	                    abort();
	                }
    	        }
        	}
        	std::cout << "No such class: " << idKey << std::endl;
			throw -1;
    	}

		/// Register a class with the factory
		static tKey RegisterCreatorFunction(tKey idKey,
											CreatorFunction classCreator) {
		    getMapFactory()->insert(std::pair<tKey,CreatorFunction>
													(idKey, classCreator));
        	return idKey;
		}

	protected:
		/// Static function to return the map, ensures it is created first.
	    static TMapFactory * getMapFactory() {
	        static TMapFactory mMapFactory;
	        return &mMapFactory;
	    }
};

/**
 * /def CREATE_OBJ(U,V)
 * Shortcut to avoid the REALLY long line that's normally required!
 * \li \a U base class
 * \li \a V class name
 */
#define CREATE_OBJ(U,V) Factory< std::string, U >::CreateInstance(V)
#endif
