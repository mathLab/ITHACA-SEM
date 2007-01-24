
namespace Nektar
{
    namespace LibUtilities 
    {
        static const unsigned int perm3A_2d[3][3] = {{0,1,2},{2,0,1},{0,2,1}};
        static const unsigned int perm3B_2d[3][3] = {{0,1,2},{1,0,2},{1,2,0}};
        static const unsigned int perm6_2d [6][3] = {{0,1,2},{1,0,2},{2,0,1},
                                                     {2,1,0},{0,2,1},{1,2,0}};
        const unsigned int NodalTetElecAvailable = 16;
        static const unsigned int NodalTetElecNPTS[NodalTetElecAvailable] = {1,2,3,4,5,7,8,10,12,14,16,19,21,24,27,30};
        static const double NodalTetElecData[][6] = {};

    }
}