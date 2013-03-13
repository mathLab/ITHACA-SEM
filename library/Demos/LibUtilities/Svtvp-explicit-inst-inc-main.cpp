#include"Svtvp-explicit-inst-inc.hpp"
#include<vector>
#include<iostream>
#include<cstdlib>

using namespace std;

double experiment(int size)
{
    vector<double> a(size, 1.0/rand());
    vector<double> b(size, 2.0/rand());
    vector<double> c(size, 3.0/rand());

    SvtvpExplInstInc::Svtvp<double>(size, 5.0, &a[0], 1, &b[0], 1, &c[0], 1);
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

