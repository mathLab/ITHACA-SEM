/*
 * Filter.h
 *
 *  Created on: 24 Feb 2012
 *      Author: cc
 */

#ifndef FILTER_H_
#define FILTER_H_

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
    class Filter;

    /// A shared pointer to a Driver object
    typedef boost::shared_ptr<Filter> FilterSharedPtr;

    /// Datatype of the NekFactory used to instantiate classes derived from
    /// the Driver class.
    typedef LibUtilities::NekFactory<
                std::string, Filter,
                const std::map<std::string, std::string>&
            > FilterFactory;
    FilterFactory& GetFilterFactory();

    class Filter
    {
    public:
        Filter();
        ~Filter();

        inline void Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
        inline void Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
        inline void Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time);
        inline bool IsTimeDependent();

    protected:
        virtual void v_Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time) = 0;
        virtual void v_Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time) = 0;
        virtual void v_Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time) = 0;
        virtual bool v_IsTimeDependent() = 0;

    };

    inline void Filter::Initialise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
    {
        v_Initialise(pFields, time);
    }

    inline void Filter::Update(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
    {
        v_Update(pFields, time);
    }

    inline void Filter::Finalise(const Array<OneD, const MultiRegions::ExpListSharedPtr> &pFields, const NekDouble &time)
    {
        v_Finalise(pFields, time);
    }

    inline bool Filter::IsTimeDependent()
    {
        return v_IsTimeDependent();
    }

}
#endif /* FILTER_H_ */
