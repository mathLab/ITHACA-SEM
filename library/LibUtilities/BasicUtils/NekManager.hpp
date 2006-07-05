#include <utility>
#include <vector>
#include <algorithm>

#include <functional> // for bind1st
#include <iostream>

#include <boost/mem_fn.hpp>
#include <boost/bind.hpp>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

using namespace std;

template <typename keyT, typename valueT,
    typename valueContainer = vector< pair<keyT, valueT> >,
    typename creatorContainer = vector< pair<keyT, valueT (*)(const keyT &key)> > >
class NekManager
{
public:
    typedef std::pair<keyT, valueT> KeyValuePair;
    typedef valueT (*CreateFunc)(const keyT &key);
    typedef std::pair<keyT, CreateFunc> KeyCreatorPair;

    NekManager() {};
    ~NekManager()
	{
		m_keyValues.clear();
		m_keyCreators.clear();
	};

    void RegisterCreator(keyT key, CreateFunc createFunc)
    {
        KeyCreatorPair createPair;

        createPair.first = key;
        createPair.second = createFunc;
        m_keyCreators.push_back(createPair);
    }

    void AddValue(keyT key, valueT value)
    {
		KeyValuePair keyValuePair(key, value);

        m_keyValues.push_back(keyValuePair);
    }

    valueT GetValue(const keyT key)
    {
        valueT returnval;

        // Search the values first, if not there
        // find a creator.  If no creator is found,
        // it is a fatal error.  If a creator is found
        // use the creator to create a value, store it,
        // then return the value.

        valueContainer::iterator loc;

        loc = find_if(m_keyValues.begin(), m_keyValues.end(),
            boost::bind<bool>(&NekManager<keyT, valueT>::IsFound<valueT>, _1, key));

        if (loc != m_keyValues.end())
        {
            returnval = (*loc).second;
        }
        else
        {
			// Not found.  Try to create one using the registered creator.
            creatorContainer::iterator creatorIter;

            creatorIter = find_if(m_keyCreators.begin(), m_keyCreators.end(),
                boost::bind<bool>(&NekManager<keyT, valueT>::IsFound<CreateFunc>, _1, key));

            if (creatorIter != m_keyCreators.end())
            {
				// Create the instance using the registered creator.
                returnval = (*((*creatorIter).second))(key);

                KeyValuePair pairT;
                pairT.first = key;
                pairT.second = returnval;
                m_keyValues.push_back(pairT);
            }
            else
            {
				// No creator registered for this key.
				NEKERROR(ErrorUtil::efatal, "Unable to find a creator registered for key.");
            }
        }

        return returnval;
    }

private:
	template <typename itemT>
    static bool IsFound(const pair<keyT, itemT> itemPair, keyT key) 
    {
        return (itemPair.first == key);
    }

    valueContainer m_keyValues;
    creatorContainer m_keyCreators;
};