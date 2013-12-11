#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/MeshPartition.h>
#include <LibUtilities/Communication/CommSerial.h>

#include <iostream>

using namespace std;
using namespace Nektar;
using namespace Nektar::LibUtilities;

class BrokenComm : public CommSerial
{
public:
    LIB_UTILITIES_EXPORT BrokenComm(int argc, char* argv[], int size) : CommSerial(argc, argv)
    {
        m_size = size;
        m_type = "Broken parallel";
    }
    LIB_UTILITIES_EXPORT BrokenComm(int size) : CommSerial(0, NULL)
    {
        m_size = size;
        m_type = "Broken parallel";
    }
    LIB_UTILITIES_EXPORT virtual ~BrokenComm() {}
    void v_SplitComm(int pRows, int pColumns)
    {
        m_commRow    = boost::shared_ptr<BrokenComm>(new BrokenComm(pColumns));
        m_commColumn = boost::shared_ptr<BrokenComm>(new BrokenComm(pRows));
    }
};

int main(int argc, char *argv[])
{
    vector<string> filenames(1);
    filenames[0] = argv[argc-1];

    CommSharedPtr vComm = boost::shared_ptr<BrokenComm>(new BrokenComm(argc, argv, 256));

    LibUtilities::SessionReaderSharedPtr vSession
        = LibUtilities::SessionReader::CreateInstance(argc, argv, filenames, vComm);

    return 0;
}

