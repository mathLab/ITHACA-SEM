#include"Svtvp-explicit-inst-inc.hpp"

namespace SvtvpExplInstInc
{

    template<class T> void Svtvp(int n, const T alpha, const T *x,
                 const int incx, const T *y, const int incy,
                 T *z, const int incz)
    {
        ++n;
        if (incx == 1 && incy == 1 && incz == 1)
        {
            while( --n )
            {
                *z = alpha * (*x) + (*y);
                ++x;
                ++y;
                ++z;
            }
        }
        else
        {
            while( --n )
            {
                *z = alpha * (*x) + (*y);
                x += incx;
                y += incy;
                z += incz;
            }
        }
    }

    // explicit instantiation
    template void Svtvp(int n, const double alpha, const double *x,
                 const int incx, const double *y, const int incy,
                 double *z, const int incz);

}