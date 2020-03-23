////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSystemCFI.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
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

#ifndef NekMeshUtils_CADSYSTEM_CFI_CADSYSTEMCFI
#define NekMeshUtils_CADSYSTEM_CFI_CADSYSTEMCFI

#include "../CADSystem.h"

#ifndef NEK_CADFIXAPI_HXX
#define NEK_CADFIXAPI_HXX
#include "cadfixapi.hxx"
#endif

namespace Nektar
{
namespace NekMeshUtils
{

class CADSystemCFI : public CADSystem
{
public:
    static CADSystemSharedPtr create(std::string name)
    {
        return MemoryManager<CADSystemCFI>::AllocateSharedPtr(name);
    }

    static std::string key;

    /**
     * @brief Default constructor.
     */
    CADSystemCFI(std::string name) : CADSystem(name, "CFI")
    {
    }

    ~CADSystemCFI()
    {
        if (m_model != nullptr)
        {
            delete m_model;
        }
    }

    bool LoadCAD();

    Array<OneD, NekDouble> GetBoundingBox();

    cfi::Model *GetCFIModel()
    {
        return m_model;
    }
    std::map<std::string, int> GetCFICurveId()
    {
        return m_nameToCurveId;
    }
    std::map<std::string, int> GetCFIFaceId()
    {
        return m_nameToFaceId;
    }
    std::map<std::string, int> GetCFIVertId()
    {
        return m_nameToVertId;
    }
    NekDouble GetScaling()
    {
        return m_scal;
    }

private:
    void AddVert(int i, cfi::Point *in);
    void AddCurve(int i, cfi::Line *in);
    void AddSurf(int i, cfi::Face *in);

    cfi::Cfi m_cfiHandle;
    cfi::Model *m_model = nullptr;
    std::vector<cfi::Body* > m_bodies;
    std::map<std::string, int> m_nameToVertId;
    std::map<std::string, int> m_nameToCurveId;
    std::map<std::string, int> m_nameToFaceId;
    std::map<std::string, std::vector<std::string> > m_mapVertToListEdge;

    NekDouble m_scal;
    bool m_useCFIMesh = false;
};

typedef std::shared_ptr<CADSystemCFI> CADSystemCFISharedPtr;
}
}

#endif
