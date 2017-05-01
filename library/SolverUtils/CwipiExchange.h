///////////////////////////////////////////////////////////////////////////////
//
// File CwipiExchange.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2015 Kilian Lackhove
//
// License for the specific language governing rights and limitations under
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
// Description: CWIPI Exchange class
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_CWIPIEXCHANGE
#define NEKTAR_CWIPIEXCHANGE

#include <FieldUtils/Interpolator.h>

#include <cwipi.h>

namespace Nektar
{
namespace SolverUtils
{

class CwipiCoupling;

typedef boost::shared_ptr<CwipiCoupling> CwipiCouplingSharedPointer;

/// Declaration of the Coupling factory
typedef LibUtilities::NekFactory<std::string,
                                 CwipiCoupling,
                                 MultiRegions::ExpListSharedPtr>
    CouplingFactory;

/// Declaration of the Coupling factory singleton
CouplingFactory &GetCouplingFactory();

class CwipiCoupling
{

public:
    typedef std::map<std::string, std::string> CouplingConfigMap;

    static std::string className;

    /// Creates an instance of this class
    static CwipiCouplingSharedPointer create(
        MultiRegions::ExpListSharedPtr field)
    {
        CwipiCouplingSharedPointer p =
            MemoryManager<CwipiCoupling>::AllocateSharedPtr(field);
        return p;
    }

    CwipiCoupling(MultiRegions::ExpListSharedPtr field);

    ~CwipiCoupling();

    const std::map<std::string, std::string> GetConfig()
    {
        return m_config;
    }

    std::vector<std::string> GetSendFieldNames()
    {
        return m_sendFieldNames;
    }

    std::vector<std::string> GetRecvFieldNames()
    {
        return m_recvFieldNames;
    }

    inline void FinalizeCoupling()
    {
        v_FinalizeCoupling();
    }

    void Send(const int step,
              const NekDouble time,
              const Array<OneD, const Array<OneD, NekDouble> > &field,
              LibUtilities::FieldMetaDataMap &fieldMetaDataMap);

    void ReceiveInterp(const int step,
                       const NekDouble time,
                       Array<OneD, Array<OneD, NekDouble> > &field,
                       LibUtilities::FieldMetaDataMap &fieldMetaDataMap);

    void SendCallback(Array<OneD, Array<OneD, NekDouble> > &interpField,
                      Array<OneD, Array<OneD, NekDouble> > &distCoords);

    void PrintProgressbar(const int position, const int goal) const;


    static void InterpCallback(
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
    std::string m_couplingName;

    CouplingConfigMap m_config;
    NekDouble m_filtWidth;

    Array<OneD, Array<OneD, NekDouble> > m_sendField;

    MultiRegions::ExpListSharedPtr m_evalField;
    MultiRegions::ExpListSharedPtr m_recvField;

    Array<OneD, Array<OneD, NekDouble> > m_oldFields;
    Array<OneD, Array<OneD, NekDouble> > m_newFields;

    int m_nSendVars;
    std::vector<std::string> m_sendFieldNames;
    int m_sendSteps;

    int m_nRecvVars;
    std::vector<std::string> m_recvFieldNames;
    int m_recvSteps;

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

    virtual void v_FinalizeCoupling(void);

    const NekDouble GetSendField(const int i, const int j) const
    {
        return m_sendField[i][j];
    }

private:
    void ReadConfig(LibUtilities::SessionReaderSharedPtr session);

    std::vector<int> GenerateVariableMapping(std::vector<std::string> &vars, std::vector<std::string> &transVars);

    void SetupReceive();

    void SetupSend();

    void SendComplete();

    void ReceiveStart();

    void Receive(const int step,
                 const NekDouble time,
                 Array<OneD, Array<OneD, NekDouble> > &field);

    void EvaluateFields(Array<OneD, Array<OneD, NekDouble> > interpField,
                        Array<OneD, Array<OneD, NekDouble> > distCoords);

    void SetupSendInterpolation();

    void AnnounceMesh();

    void DumpRawFields(const NekDouble time,
                       Array<OneD, Array<OneD, NekDouble> > rVals);

    void OverrrideFields(Array<OneD, Array<OneD, NekDouble> > &rVals);

    template <typename T>
    void AddElementsToMesh(T geom,
                           int &coordsPos,
                           int &connecPos,
                           int &conidxPos);
};

typedef boost::function<void(Array<OneD, Array<OneD, NekDouble> > &interpField,
                             Array<OneD, Array<OneD, NekDouble> > &distCoords)>
    SendCallbackType;

static std::map<std::string, SendCallbackType> SendCallbackMap;

}
}

#endif
