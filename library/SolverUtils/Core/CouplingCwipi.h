///////////////////////////////////////////////////////////////////////////////
//
// File CouplingCwipi.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2017 Kilian Lackhove
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
// Description: CWIPI Coupling class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_COUPLINGCWIPI
#define NEKTAR_COUPLINGCWIPI

#include <SolverUtils/Core/Coupling.h>

#include <FieldUtils/Interpolator.h>

#include <cwipi.h>

#include <functional>

namespace Nektar
{
namespace SolverUtils
{

class CouplingCwipi;

class CouplingCwipi : public Coupling
{

public:
    static std::string className;

    /// Creates an instance of this class
    static CouplingSharedPtr create(MultiRegions::ExpListSharedPtr field)
    {
        CouplingSharedPtr p =
            MemoryManager<CouplingCwipi>::AllocateSharedPtr(field);
        p->Init();
        return p;
    }

    SOLVER_UTILS_EXPORT CouplingCwipi(MultiRegions::ExpListSharedPtr field);

    SOLVER_UTILS_EXPORT virtual ~CouplingCwipi();

    SOLVER_UTILS_EXPORT void SendCallback(
        Array<OneD, Array<OneD, NekDouble> > &interpField,
        Array<OneD, Array<OneD, NekDouble> > &distCoords);

    SOLVER_UTILS_EXPORT static void InterpCallback(
        const int entities_dim,
        const int n_local_vertex,
        const int n_local_element,
        const int n_local_polhyedra,
        const int n_distant_point,
        const double local_coordinates[],
        const int local_connectivity_index[],
        const int local_connectivity[],
        const int local_polyhedra_face_index[],
        const int local_polyhedra_cell_to_face_connectivity[],
        const int local_polyhedra_face_connectivity_index[],
        const int local_polyhedra_face_connectivity[],
        const double distant_points_coordinates[],
        const int distant_points_location[],
        const float distant_points_distance[],
        const int distant_points_barycentric_coordinates_index[],
        const double distant_points_barycentric_coordinates[],
        const int stride,
        const cwipi_solver_type_t solver_type,
        const void *local_field,
        void *distant_field);

protected:
    NekDouble m_filtWidth;

    Array<OneD, Array<OneD, NekDouble> > m_sendField;

    MultiRegions::ExpListSharedPtr m_recvField;

    Array<OneD, Array<OneD, NekDouble> > m_oldFields;
    Array<OneD, Array<OneD, NekDouble> > m_newFields;

    int m_sendHandle;
    int m_recvHandle;

    int m_sendTag;
    int m_recvTag;

    int m_lastSend;
    int m_lastReceive;

    int m_spacedim;
    int m_nPoints;
    double *m_points;
    double *m_coords;
    int *m_connecIdx;
    int *m_connec;
    double *m_rValsInterl;
    double *m_sValsInterl;

    std::map<int, int> m_vertMap;

    FieldUtils::InterpolatorSharedPtr m_sendInterpolator;

    FieldUtils::InterpolatorSharedPtr m_extrapInterpolator;

    SOLVER_UTILS_EXPORT virtual void v_Init();

    SOLVER_UTILS_EXPORT virtual void v_Send(
        const int step,
        const NekDouble time,
        const Array<OneD, const Array<OneD, NekDouble> > &field,
        std::vector<std::string> &varNames);

    SOLVER_UTILS_EXPORT virtual void v_Receive(
        const int step,
        const NekDouble time,
        Array<OneD, Array<OneD, NekDouble> > &field,
        std::vector<std::string> &varNames);

    SOLVER_UTILS_EXPORT virtual void v_Finalize();

    SOLVER_UTILS_EXPORT NekDouble GetSendField(const int i,
                                               const int j) const
    {
        return m_sendField[i][j];
    }

private:
    void ReadConfig(LibUtilities::SessionReaderSharedPtr session);

    void SetupReceive();

    void SetupSend();

    void SendComplete();

    void ReceiveStart();

    void ReceiveCwipi(const int step,
                      const NekDouble time,
                      Array<OneD, Array<OneD, NekDouble> > &field);

    void EvaluateFields(Array<OneD, Array<OneD, NekDouble> > interpField,
                        Array<OneD, Array<OneD, NekDouble> > distCoords);

    void SetupSendInterpolation();

    void AnnounceMesh();

    void DumpRawFields(const NekDouble time,
                       Array<OneD, Array<OneD, NekDouble> > rVals);

    void ExtrapolateFields(Array<OneD, Array<OneD, NekDouble> > &rVals,
                           Array<OneD, int> &notLoc);

    template <typename T>
    void AddElementsToMesh(T geom,
                           int &coordsPos,
                           int &connecPos,
                           int &conidxPos);
};

typedef std::function<void(Array<OneD, Array<OneD, NekDouble> > &interpField,
                           Array<OneD, Array<OneD, NekDouble> > &distCoords)>
    SendCallbackType;

static std::map<std::string, SendCallbackType> SendCallbackMap;
}
}

#endif
