/*
 * Core.h
 *
 *  Created on: 17 Jul 2013
 *      Author: cc
 */

#ifndef MISC_H_
#define MISC_H_

#include <string>

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar {
namespace SolverUtils {

    typedef std::vector<std::pair<std::string, std::string> > SummaryList;

    /// Adds a summary item to the summary info list
    SOLVER_UTILS_EXPORT void AddSummaryItem(
            SummaryList& l,
            const std::string& name,
            const std::string& value);

    /// Adds a summary item to the summary info list
    SOLVER_UTILS_EXPORT void AddSummaryItem(
            SummaryList& l,
            const std::string& name,
            const int& value);

    /// Adds a summary item to the summary info list
    SOLVER_UTILS_EXPORT void AddSummaryItem(
            SummaryList& l,
            const std::string& name,
            const NekDouble& value);

}
}

#endif /* CORE_H_ */
