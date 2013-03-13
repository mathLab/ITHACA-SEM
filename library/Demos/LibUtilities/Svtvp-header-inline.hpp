namespace SvtvpHeaderInline
{
    /// \brief  svtvp (scalar times vector plus vector): z = alpha*x + y
    template<class T> inline void Svtvp(int n, const T alpha, const T *x,
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
}