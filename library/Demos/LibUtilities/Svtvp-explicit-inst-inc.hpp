namespace SvtvpExplInstInc
{
    /// \brief  svtvp (scalar times vector plus vector): z = alpha*x + y
    template<class T> void Svtvp(int n, const T alpha, const T *x,
                 const int incx, const T *y, const int incy,
                 T *z, const int incz);
}