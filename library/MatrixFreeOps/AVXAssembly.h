#ifndef AVXASSEMBLY_H
#define AVXASSEMBLY_H

#include <MultiRegions/ContField.h>
#include <Collections/Collection.h>
#include <iostream>
#include "VecData.hpp"
#include "Operator.hpp"
#include "AVXUtil.hpp"

using namespace Nektar;

template<int VW>
struct AVXAssembly
{
public:
    AVXAssembly(MultiRegions::AssemblyMapSharedPtr asmMap,
                MultiRegions::ExpListSharedPtr expList,
                const int nElmt)
        : m_asmMap    (asmMap),
          m_signChange(asmMap->GetSignChange()),
          m_nLocal    (asmMap->GetNumLocalCoeffs())
    {
        int nLocalPad = nElmt * expList->GetExp(0)->GetNcoeffs();

        // Copy & pad local to global array
        m_l2g = Array<OneD, int>(nLocalPad, asmMap->GetNumGlobalCoeffs());
        Vmath::Vcopy(asmMap->GetNumLocalCoeffs(),
                     asmMap->GetLocalToGlobalMap(), 1, m_l2g, 1);
        TransposeData(nElmt, VW, m_l2g);

        // Coyp & pad local to global sign array
        if(m_signChange)
        {
            m_l2gSign = Array<OneD, double>(nLocalPad, 0.0);
            Vmath::Vcopy(asmMap->GetNumLocalCoeffs(),
                         asmMap->GetLocalToGlobalSign(), 1, m_l2gSign, 1);
            TransposeData(nElmt, VW, m_l2gSign);
        }
    }

    template<int NM>
    inline int AssembleBlock(const double* local, double* global, int map_index)
    {
        if (m_signChange)
        {
#if 1
            for (int i = 0; i < NM * VW; ++i)
            {
                global[m_l2g[map_index + i]] += m_l2gSign[map_index + i]
                    * local[i];
            }
            map_index += NM * VW;

#else
            for(int i = 0; i < NM; i++, map_index += VW)
            {
                VecData<int, VW> index(&m_l2g[map_index]);
                VecData<double, VW> sign(&m_l2gSign[map_index]);
                VecData<double, VW> local_val(&local[i*VW]);
                VecData<double, VW> out = gather(global, index);
                out.fma(sign,local_val);
                out.scatter(global, index);
            }
#endif
        }
        else
        {
#if 1
            for (int i = 0; i < NM * VW; ++i)
            {
                global[m_l2g[map_index + i]] += local[i];
            }
            map_index += NM * VW;
#else
            for(int i = 0; i < NM; i++, map_index += VW)
            {
                VecData<int, VW> index(&m_l2g[map_index]);
                VecData<double, VW> local_val(&local[i*VW]);
                VecData<double, VW> out = gather(global, index) + local_val;
                out.scatter(global, index);
            }
#endif
        }

        return map_index;

    }

    void Assemble(const Array<OneD, const NekDouble> &local,
                  Array<OneD, NekDouble> &global)
    {
        if(m_signChange) //Consider templating based bool
        {
            for(int i = 0; i < m_nLocal; i += VW)
            {
                VecData<int, VW> index(&m_l2g[i]);
                VecData<double, VW> sign(&m_l2gSign[i]);
                VecData<double, VW> local_val(&local[i]);
                VecData<double, VW> out = gather(&global[0], index);
                out.fma(sign,local_val);
                out.scatter(&global[0], index);
            }
        }
        else
        {
            for(int i = 0; i < m_nLocal; i += VW)
            {
                VecData<int, VW> index(&m_l2g[i]);
                VecData<double, VW> local_val(&local[i]);
                VecData<double, VW> out = gather(&global[0], index) + local_val;
                out.scatter(&global[0], index);
            }
        }
        m_asmMap->UniversalAssemble(global);
    }

    template<int NM>
    inline void ScatterBlock(const double* global,
                      double *local,
                      int map_index)
    {
        if(m_signChange)
        {
#if 1
            for (int i = 0; i < NM * VW; ++i)
            {
                local[i] = m_l2gSign[map_index + i] *
                    global[m_l2g[map_index + i]];
            }
            map_index += NM * VW;
#else
            for(int i = 0; i < NM; i++, map_index += VW)
            {
                VecData<int,VW> index(&m_l2g[map_index]);
                VecData<double,VW> sign(&m_l2gSign[map_index]);
                VecData<double,VW> out = sign * gather(global, index);
                out.store(&local[i*VW]);
            }
#endif
        }
        else
        {
#if 1
            for (int i = 0; i < NM * VW; ++i)
            {
                local[i] = global[m_l2g[map_index + i]];
            }
            map_index += NM * VW;
#else
            for(int i = 0; i < NM; i++, map_index += VW)
            {
                VecData<int,VW> index(&m_l2g[map_index]);
                VecData<double,VW> out = gather(global, index);
                out.store(&local[i*VW]);
            }
#endif
        }
    }

    void Scatter(const Array<OneD, const NekDouble> &global,
                 Array<OneD, NekDouble> &local)
    {
        if (m_signChange)
        {
            for (int i = 0; i < m_nLocal; i += VW)
            {
                VecData<int,VW> index(&m_l2g[i]);
                VecData<double,VW> sign(&m_l2gSign[i]);
                VecData<double,VW> out = sign * gather(&global[0], index);
                out.store(&local[i]);
            }
        }
        else
        {
            for (int i = 0; i < m_nLocal; i += VW)
            {
                VecData<int,VW> index(&m_l2g[i]);
                VecData<double,VW> out = gather(&global[0], index);
                out.store(&local[i]);
            }
        }
    }

    MultiRegions::AssemblyMapSharedPtr GetAssemblyMap()
    {
        return m_asmMap;
    }

    int GetNGlobal()
    {
        return m_asmMap->GetNumGlobalCoeffs();
    }

    MultiRegions::AssemblyMapSharedPtr m_asmMap;
    Array<OneD, int> m_l2g;
    Array<OneD, double> m_l2gSign;
    bool m_signChange;
    int m_nLocal;
};
template<int VW, int NM>
struct AVXAssemblyOld
{
public:
    AVXAssemblyOld(MultiRegions::AssemblyMapSharedPtr asmMap,
                const int nElmt)
        : m_asmMap(asmMap), m_l2g(asmMap->GetLocalToGlobalMap()),
          m_signChange(asmMap->GetSignChange()),
          m_nLocal(asmMap->GetNumLocalCoeffs()),
          m_vecWidth(VW)
    {
        // Grab original local to global map & transpose
        TransposeData(nElmt, VW, m_l2g);

        if (m_signChange)
        {
            m_l2gSign = asmMap->GetLocalToGlobalSign();
            TransposeData(nElmt, VW, m_l2gSign);
        }

    }

    void CheckIndices()
    {
        std::cout << "it is " << m_nLocal << std::endl;
        for(int i = 0; i < m_nLocal; i += m_vecWidth)
        {
            int idx0 = m_l2g[i];
            int idx1 = m_l2g[i + 1];
            int idx2 = m_l2g[i + 2];
            int idx3 = m_l2g[i + 3];

            if( idx0 == idx1 || idx0 == idx2 || idx0 == idx3 || idx1 == idx2 || idx1 == idx3 || idx2 == idx3 )
            {
                std::cout << i << std::endl;
            }
        }
    }

    int AssembleBlock(
        const double* local,
        double* global,
        int map_index)
    {
        if(m_signChange)
        {
            for(int i = 0; i < NM*VW; i++, map_index++)
            {
                int index = m_l2g[map_index];
                double sign = m_l2gSign[map_index];
                double local_val = local[i];
                global[index] += local_val * sign;
            }
        }
        else
        {
            for(int i = 0; i < NM*VW; i++, map_index++)
            {
                int index = m_l2g[map_index];
                double local_val = local[i];
                global[index] += local_val;
            }
        }

        return map_index;

    }

    void Assemble(const Array<OneD, const NekDouble> &local,
                  Array<OneD, NekDouble> &global)
    {
        if (m_signChange)
        {
            Vmath::Assmb(m_nLocal, &m_l2gSign[0], &local[0], &m_l2g[0],
                         &global[0]);
        }
        else
        {
            Vmath::Assmb(m_nLocal, &local[0], &m_l2g[0], &global[0]);
        }
        m_asmMap->UniversalAssemble(global);
    }

    void ScatterBlock(const double* global,
                    double *local,
                    int map_index)
    {
        if(m_signChange)
        {
            for(int i = 0; i < NM*VW; i++, map_index++)
            {
                int index = m_l2g[map_index];
                double sign = m_l2gSign[map_index];
                local[i] = sign * global[index];
            }
        }
        else
        {
            for(int i = 0; i < NM*VW; i++, map_index++)
            {
                int index = m_l2g[map_index];
                local[i] = global[index];
            }
        }
    }

    void Scatter(const Array<OneD, const NekDouble> &global,
                 Array<OneD, NekDouble> &local)
    {
        if (m_signChange)
        {
            Vmath::Gathr(m_nLocal, &m_l2gSign[0], &global[0],
                         &m_l2g[0], &local[0]);
        }
        else
        {
            for (int i = 0; i < m_nLocal; ++i)
            {
                local[i] = global[m_l2g[i]];
            }
        }
    }

    MultiRegions::AssemblyMapSharedPtr GetAssemblyMap()
    {
        return m_asmMap;
    }

private:
    MultiRegions::AssemblyMapSharedPtr m_asmMap;
    Array<OneD, int> m_l2g;
    Array<OneD, double> m_l2gSign;
    bool m_signChange;
    int m_nLocal;
    int m_vecWidth;
};

#endif
