#include"Svtvp-explicit-inst.hpp"

namespace SvtvpExplInst
{

    template<class T> void Svtvp(int n, const T alpha, const T *x, const T *y, T *z)
    {
        ++n;
            while( --n )
            {
                *z = alpha * (*x) + (*y);
                ++x;
                ++y;
                ++z;
            }
    }

    // explicit instantiation
    template void Svtvp(int n, const double alpha,
                               const double *x,
                               const double *y,
                                     double *z);

}