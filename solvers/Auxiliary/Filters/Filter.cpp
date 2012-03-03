/*
 * Filter.cpp
 *
 *  Created on: 24 Feb 2012
 *      Author: cc
 */

#include <Auxiliary/Filters/Filter.h>

namespace Nektar
{
    FilterFactory& GetFilterFactory()
    {
        typedef Loki::SingletonHolder<FilterFactory,
            Loki::CreateUsingNew,
            Loki::NoDestroy > Type;
        return Type::Instance();
    }


    Filter::Filter()
    {

    }

    Filter::~Filter()
    {

    }

}
