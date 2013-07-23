/*
 * Core.cpp
 *
 *  Created on: 17 Jul 2013
 *      Author: cc
 */

#include <utility>
#include <vector>
using namespace std;

#include <boost/lexical_cast.hpp>

#include <SolverUtils/Core/Misc.h>

namespace Nektar {
namespace SolverUtils {
    /**
     * Adds an item to a SummaryList
     */
    void AddSummaryItem(
            SummaryList& l,
            const std::string& name,
            const std::string& value)
    {
        l.push_back(std::make_pair<std::string, std::string>(name, value));
    }

    /// Adds a summary item to the summary info list
    void AddSummaryItem(
            SummaryList& l,
            const std::string& name,
            const int& value)
    {
        l.push_back(std::make_pair<std::string, std::string>(name, boost::lexical_cast<std::string>(value)));
    }

    /// Adds a summary item to the summary info list
    void AddSummaryItem(
            SummaryList& l,
            const std::string& name,
            const NekDouble& value)
    {
        l.push_back(std::make_pair<std::string, std::string>(name, boost::lexical_cast<std::string>(value)));
    }

}
}
