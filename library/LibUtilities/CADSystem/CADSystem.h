////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSystem.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_LIB_UTILITIES_CADSYSTEM_CADSYSTEM_H
#define NEKTAR_LIB_UTILITIES_CADSYSTEM_CADSYSTEM_H

#include <boost/shared_ptr.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <Standard_Macro.hxx>

#include <LibUtilities/CADSystem/OpenCascade.h>
#include <LibUtilities/CADSystem/CADCurve.h>
#include <LibUtilities/CADSystem/CADSurf.h>

namespace Nektar {
namespace LibUtilities {

/**
 * @brief Base class for CAD interface system.
 *
 * A class which can load and interact with cad for Nektar++ using opencascade
 * This class contains maps to subclasses surface and curves
 */

class CADSystem
{
    public:
        friend class MemoryManager<CADSystem>;

        /**
         * @brief Defualt constructor.
         */

        LIB_UTILITIES_EXPORT CADSystem(const std::string &name) : m_name(name)
        {
        };

        LIB_UTILITIES_EXPORT std::string GetName();

        LIB_UTILITIES_EXPORT bool LoadCAD();

        LIB_UTILITIES_EXPORT void Report();

        LIB_UTILITIES_EXPORT Array<OneD, NekDouble> GetBoundingBox();

        LIB_UTILITIES_EXPORT int GetNumSurf()
        {
            return m_surfs.size();
        }

        LIB_UTILITIES_EXPORT int GetNumCurve()
        {
            return m_curves.size();
        }

        /**
         * @brief Gets curve type from map.
         */
        LIB_UTILITIES_EXPORT const CADCurveSharedPtr GetCurve(int i)
        {
            std::map<int,CADCurveSharedPtr>::iterator
                                    search = m_curves.find(i);
            ASSERTL0(search != m_curves.end(), "curve does not exist");

            return search->second;
        }
        /**
         * @brief Gets suface from map.
         */
        LIB_UTILITIES_EXPORT CADSurfSharedPtr GetSurf(int i)
        {
            std::map<int,CADSurfSharedPtr>::iterator
                            search = m_surfs.find(i);
            ASSERTL0(search != m_surfs.end(), "surface does not exist");

            return search->second;
        }

        LIB_UTILITIES_EXPORT int GetEPC()
        {
            return m_epc;
        }

    private:

        /// private function to add curve to map
        void AddCurve(int i, TopoDS_Shape in);
        /// private function to add surf to map
        void AddSurf(int i, TopoDS_Shape in,
                     std::vector<std::vector<std::pair<int,int> > > ein);
        /// name of cad file to be opened including ext
        std::string m_name;
        /// euler poincare number of the cad
        int m_epc;
        /// map of curves
        std::map<int,CADCurveSharedPtr> m_curves;
        /// map of surfaces
        std::map<int,CADSurfSharedPtr> m_surfs;
};

typedef boost::shared_ptr<CADSystem> CADSystemSharedPtr;

}
}

#endif
