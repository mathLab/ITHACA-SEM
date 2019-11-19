#include <Collections/Collection.h>

using namespace Nektar;

void HexFlops(Collections::ImplementationType  impType,
              int                              order,
              int                              nElmt,
              NekDouble                       &gflop,
              NekDouble                       &matSize)
{
    const int nM = order + 1;
    const int nQ = order + 2;

    if (impType == Collections::eIterPerExp ||
        impType == Collections::eSumFac)
    {
        gflop = 2.0 * nElmt * (nM * nM * nM * nQ + nQ * nQ * nM * nM +
                               nQ * nQ * nQ * nM);
        matSize =
            ((nQ * nM) + (nM * nM * nM)) +// + (nQ * nM * nM)) + // m_funcs[0]
            ((nQ * nM) + (nM * nQ)) * nM +// + (nQ * nQ)) * nM      + // m_funcs[1]
            ((nQ * nQ * nM) + (nM * nQ));// + (nQ * nQ * nQ));  // m_funcs[2]
        matSize *= nElmt;
    }
    else if (impType == Collections::eStdMat)
    {
        gflop = 2.0 * nElmt * nM * nM * nM * nQ * nQ * nQ;
        matSize = nQ*nQ*nQ * nM*nM*nM * nElmt;
    }

    gflop *= 1e-9;
    matSize *= sizeof(NekDouble);
    matSize /= 1024.0 * 1024.0 * 1024.0;
}

void TetFlops(Collections::ImplementationType  impType,
              int                              order,
              int                              nElmt,
              NekDouble                       &gflop,
              NekDouble                       &matSize)
{
    const int nM = order + 1;
    const int nQ = order + 2;

    if (impType == Collections::eIterPerExp ||
        impType == Collections::eSumFac)
    {
        gflop = 2.0 * nM * nQ * nQ * nQ;

        for (int i = 0; i < nM; ++i)
        {
            for (int j = 0; j < nM - i; ++j)
            {
                gflop += 2.0 * nQ * (nM - i - j);
            }

            gflop += 2.0 * (nM - i) * nQ * nQ;
        }

        gflop *= nElmt;
    }
    else if (impType == Collections::eStdMat)
    {
        gflop = 2.0 * nElmt * nQ * nQ * nQ *
            LibUtilities::StdTetData::getNumberOfCoefficients(nM, nM, nM);
    }

    gflop *= 1e-9;
    matSize = 0;
}

void TetFlopsIProduct(Collections::ImplementationType  impType,
                      int                              order,
                      int                              nElmt,
                      NekDouble                       &gflop,
                      NekDouble                       &matSize)
{
    const int nM = order + 1;
    const int nQ = order + 2;
    if (impType == Collections::eIterPerExp ||
        impType == Collections::eSumFac)
    {
        gflop = 2.0 * nQ*nQ * nM * nQ;

        for (int i = 0; i < nM; ++i)
        {
            gflop += 2.0 * nQ * (nM-i) * nQ;
            gflop += 2.0 * nQ * nQ;

            for (int j = 0; j < nM - i; ++j)
            {
                gflop += 2.0 * nQ * (nM-i-j);
            }
        }

        gflop += 2.0 * nQ;

        // Account for multiplication by weights & jacobian
        gflop += 4.0 * nQ * nQ * nQ;
        gflop *= nElmt;
    }
    else
    {
        // haven't othered yet
        gflop = 0.0;
    }

    gflop *= 1e-9;
    matSize = 0.0;
}
