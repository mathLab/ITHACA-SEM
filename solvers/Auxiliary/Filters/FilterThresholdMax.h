/*
 * FilterThresholdMax.h
 *
 *  Created on: 24 Feb 2012
 *      Author: cc
 */

#ifndef FILTERTHRESHOLDMAX_H_
#define FILTERTHRESHOLDMAX_H_

#include <Auxiliary/Filters/Filter.h>

namespace Nektar
{

    class FilterThresholdMax : public Filter
    {
    public:
        friend class MemoryManager<FilterThresholdMax>;

        /// Creates an instance of this class
        static FilterSharedPtr create(const std::map<std::string, std::string> &pParams) {
            FilterSharedPtr p = MemoryManager<FilterThresholdMax>::AllocateSharedPtr(pParams);
            //p->InitObject();
            return p;
        }

        ///Name of the class
        static std::string className;

        FilterThresholdMax(const std::map<std::string, std::string> &pParams);
        ~FilterThresholdMax();

    protected:
        virtual void v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
        virtual void v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
        virtual void v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
        virtual bool v_IsTimeDependent();

    private:
        Array<OneD, NekDouble> m_threshold;
        NekDouble m_thresholdValue;
        NekDouble m_initialValue;
        std::string m_outputFile;
    };

}

#endif /* FILTERTHRESHOLDMAX_H_ */
