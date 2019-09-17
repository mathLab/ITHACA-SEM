///////////////////////////////////////////////////////////////////////////////
//
// File Interpolator.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2016 Kilian Lackhove
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
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
// Description: Interpolator
//
///////////////////////////////////////////////////////////////////////////////

#ifndef FIELDUTILS_INTERPOLATOR_H
#define FIELDUTILS_INTERPOLATOR_H

#include <LibUtilities/BasicUtils/Interpolator.h>
#include <MultiRegions/ExpList.h>
#include <FieldUtils/FieldUtilsDeclspec.h>

namespace Nektar
{
namespace FieldUtils
{

/// A class that contains algorithms for interpolation between pts fields,
/// expansions and different meshes
class Interpolator : public LibUtilities::Interpolator
{
public:
    /**
    * @brief Constructor of the Interpolator class
    *
    * @param method    interpolation method, defaults to a sensible value if not
    * set
    * @param coordId   coordinate id along which the interpolation should be
    * performed
    * @param filtWidth filter width, required by some algorithms such as eGauss
    * @param maxPts    limit number of considered points
    *
    * if method is not specified, the best algorithm is chosen autpomatically.
    *
    * If coordId is not specified, a full 1D/2D/3D interpolation is performed
    * without
    * collapsing any coordinate.
    *
    * filtWidth must be specified for the eGauss algorithm only.
    */
    Interpolator(LibUtilities::InterpMethod method    = LibUtilities::eNoMethod,
                 short int                  coordId   = -1,
                 NekDouble                  filtWidth = 0.0,
                 int                        maxPts    = 1000)
        : LibUtilities::Interpolator(method, coordId, filtWidth, maxPts)
    {
    }

    /// Interpolate from an expansion to an expansion
    FIELD_UTILS_EXPORT void Interpolate(
        const std::vector<MultiRegions::ExpListSharedPtr> expInField,
        std::vector<MultiRegions::ExpListSharedPtr> &expOutField,
        NekDouble def_value = 0.0);

    /// Interpolate from an expansion to a pts field
    FIELD_UTILS_EXPORT void Interpolate(
        const std::vector<MultiRegions::ExpListSharedPtr> expInField,
        LibUtilities::PtsFieldSharedPtr &ptsOutField,
        NekDouble def_value = 0.0);

    /// Interpolate from a pts field to an expansion
    FIELD_UTILS_EXPORT void Interpolate(
        const LibUtilities::PtsFieldSharedPtr ptsInField,
        std::vector<MultiRegions::ExpListSharedPtr> &expOutField);

    /// Interpolate from a pts field to a pts field
    FIELD_UTILS_EXPORT void Interpolate(
        const LibUtilities::PtsFieldSharedPtr ptsInField,
        LibUtilities::PtsFieldSharedPtr &ptsOutField);

protected:
    /// input field
    std::vector<MultiRegions::ExpListSharedPtr> m_expInField;
    /// output field
    std::vector<MultiRegions::ExpListSharedPtr> m_expOutField;
};

typedef std::shared_ptr<Interpolator> InterpolatorSharedPtr;

}
}

#endif
