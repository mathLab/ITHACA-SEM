#ifndef H_LOC2GLOMAP

#include <vector>

namespace MultiRegions{
  
  class Loc2GloMap{
  private:
    int *_nbndry;
    std::vector <int *> _locid;
    std::vector <int *> _gmap;
  protected:

  public:
    Loc2GloMap();
    ~Loc2GloMap();

  };

}

#define H_LOC2GLOMAP
#endif
