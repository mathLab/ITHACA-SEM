#ifndef FIELDUTILS_FIELDCONVERTCOMM
#define FIELDUTILS_FIELDCONVERTCOMM

#include <LibUtilities/Communication/CommSerial.h>

namespace Nektar
{
namespace FieldUtils
{

class FieldConvertComm : public LibUtilities::CommSerial
{
public:
    FieldConvertComm(int argc, char *argv[], int size, int rank)
        : CommSerial(argc, argv)
    {
        m_size = size;
        m_rank = rank;
		m_type = "FieldConvert parallel";
    }
    FieldConvertComm(int size, int rank) : CommSerial(0, NULL)
    {
        m_size = size;
        m_rank = rank;
		m_type = "FieldConvert parallel";
    }
    virtual ~FieldConvertComm()
    {
    }
    void v_SplitComm(int pRows, int pColumns)
    {
        // Compute row and column in grid.
        m_commRow = std::shared_ptr<FieldConvertComm>(
            new FieldConvertComm(pColumns, m_rank));
        m_commColumn =
            std::shared_ptr<FieldConvertComm>(new FieldConvertComm(pRows, 0));
    }

protected:
    int v_GetRank(void)
    {
        return m_rank;
    }

    bool v_TreatAsRankZero(void)
    {
        return true;
    }

    bool v_IsSerial(void)
    {
        return true;
    }

    bool v_RemoveExistingFiles(void)
    {
        return false;
    }

private:
    int m_rank;
};
}
}
#endif
