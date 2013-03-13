namespace SvtvpExplInst
{
    /// \brief  svtvp (scalar times vector plus vector): z = alpha*x + y
    template<class T> void Svtvp(int n, const T alpha, const T *x, const T *y, T *z);
}