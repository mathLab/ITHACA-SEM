#include"Svtvp-header.hpp"
#include<LibUtilities/BasicUtils/Vmath.hpp>
#include<vector>
#include<iostream>
#include<cstdlib>

using namespace std;

double experiment(int size)
{
    vector<double> a(size, 1.0/rand());
    vector<double> b(size, 2.0/rand());
    vector<double> c(size, 3.0/rand());

    SvtvpHeader::Svtvp<double>(size, 5.0, &a[0], 1, &b[0], 1, &c[0], 1);
    Vmath::Svtvp(size, 5.0, &a[0], 1, &b[0], 1, &c[0], 1);
    return c[size-1];
}

double experimentVmath(int size)
{
    vector<double> a(size, 1.0/rand());
    vector<double> b(size, 2.0/rand());
    vector<double> c(size, 3.0/rand());

    Vmath::Svtvp(size, 5.0, &a[0], 1, &b[0], 1, &c[0], 1);
    return c[size-1];
}


int main()
{
    int n = 1000000;

    for (int i = 0; i < n; i++)
    {
        experiment(10000);
    }

    
}

