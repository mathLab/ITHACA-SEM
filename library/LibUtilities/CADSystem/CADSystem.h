





#ifndef NEKTAR_LIB_UTILITIES_CADSYSTEM_CADSYSTEM_H
#define NEKTAR_LIB_UTILITIES_CADSYSTEM_CADSYSTEM_H

#include <string>

#include <boost/shared_ptr.hpp>

#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

namespace Nektar {
    namespace LibUtilities {

	class CADSystem
	{
	public:
	    friend class MemoryManager<CADSystem>;

	    LIB_UTILITIES_EXPORT CADSystem(std::string name) : m_name(name)
	    {
	    }
	    LIB_UTILITIES_EXPORT std::string GetName();

	private:
	    std::string m_name;
	};

	typedef boost::shared_ptr<CADSystem> CADSystemSharedPtr;
    }
}

#endif
